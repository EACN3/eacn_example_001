"""PoC v3: Cell-Neighborhood Expression Mismatch (CNEM) framework.

Core idea: After integration, measure how well each cell's own expression
matches its new neighborhood's expression. Rare cells absorbed into large
clusters will have HIGH mismatch (their expression differs from neighbors).

This is fundamentally different from v1 (neighbor identity) and v2 (neighbor
coherence): we measure the cell-TO-neighborhood fit, not within-neighborhood
consistency.

Three metrics:
1. CNEM: cosine distance between cell expression and mean neighbor expression
2. LNTEC: ratio of post/pre neighborhood expression variance
3. DCR: displacement / neighborhood quality ratio
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics import roc_auc_score
import warnings, os, time
warnings.filterwarnings('ignore')

DATA_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/poc_v1/data'
OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/poc_v3/results'
os.makedirs(OUT_DIR, exist_ok=True)

K = 30

# ============================================================
# 1. Load and preprocess
# ============================================================
print("=" * 60)
print("Step 1: Load and preprocess")
print("=" * 60)
adata = ad.read_h5ad(f'{DATA_DIR}/pancreas.h5ad')

adata.obs['celltype_true'] = adata.obs['celltype'].copy()
rare_types = ['epsilon', 'schwann', 't_cell', 'mast']
adata.obs['is_rare'] = adata.obs['celltype'].isin(rare_types).astype(int)

sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='tech')

# Normalized expression for mismatch computation (NOT scaled)
X_expr = adata[:, adata.var['highly_variable']].X.copy()
if hasattr(X_expr, 'toarray'):
    X_expr = X_expr.toarray()
X_expr = np.array(X_expr, dtype=np.float32)

# Scaled + PCA for integration
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

print(f"Cells: {adata.n_obs}, HVG: {X_expr.shape[1]}")
print(f"Rare types: {adata.obs[adata.obs['is_rare']==1]['celltype'].value_counts().to_dict()}")

# ============================================================
# 2. Pre-integration baselines (per-batch)
# ============================================================
print("\n" + "=" * 60)
print("Step 2: Pre-integration baselines")
print("=" * 60)

batches = adata.obs['tech'].unique()
pre_cnem = np.zeros(adata.n_obs)  # cell-neighborhood mismatch
pre_var = np.zeros(adata.n_obs)   # neighborhood expression variance

for batch in batches:
    mask = adata.obs['tech'] == batch
    idx = np.where(mask)[0]
    X_pca_b = adata.obsm['X_pca'][idx]
    X_expr_b = X_expr[idx]

    k_use = min(K, len(idx) - 1)
    nn = NearestNeighbors(n_neighbors=k_use + 1, metric='euclidean')
    nn.fit(X_pca_b)
    _, indices = nn.kneighbors(X_pca_b)

    for i_local, i_global in enumerate(idx):
        nbr_local = indices[i_local, 1:]
        cell_e = X_expr_b[i_local]
        nbr_e = X_expr_b[nbr_local]
        # CNEM: 1 - cosine_sim(cell, mean_neighbor)
        mean_nbr = nbr_e.mean(axis=0, keepdims=True)
        pre_cnem[i_global] = 1.0 - cosine_similarity(cell_e.reshape(1, -1), mean_nbr)[0, 0]
        # Neighborhood variance
        dists = 1.0 - cosine_similarity(cell_e.reshape(1, -1), nbr_e)[0]
        pre_var[i_global] = np.var(dists)

print(f"Pre-CNEM: mean={pre_cnem.mean():.4f}, std={pre_cnem.std():.4f}")

# ============================================================
# 3. Integration + post-integration metrics
# ============================================================
print("\n" + "=" * 60)
print("Step 3: Integration methods + CNEM/LNTEC/DCR")
print("=" * 60)

def compute_metrics(embedding, X_expr, pre_cnem, pre_var, X_pca_pre, k=K):
    """Compute CNEM, LNTEC, DCR for all cells."""
    n = embedding.shape[0]
    nn = NearestNeighbors(n_neighbors=k + 1, metric='euclidean')
    nn.fit(embedding)
    _, indices = nn.kneighbors(embedding)

    post_cnem = np.zeros(n)
    post_var = np.zeros(n)
    displacement = np.zeros(n)

    for i in range(n):
        nbr_idx = indices[i, 1:]
        cell_e = X_expr[i]
        nbr_e = X_expr[nbr_idx]

        # CNEM: mismatch between cell and its post-integration neighborhood
        mean_nbr = nbr_e.mean(axis=0, keepdims=True)
        post_cnem[i] = 1.0 - cosine_similarity(cell_e.reshape(1, -1), mean_nbr)[0, 0]

        # Neighborhood variance (for LNTEC)
        dists = 1.0 - cosine_similarity(cell_e.reshape(1, -1), nbr_e)[0]
        post_var[i] = np.var(dists)

    # CNEM change: positive = mismatch increased after integration
    cnem_delta = post_cnem - pre_cnem

    # LNTEC: ratio of post/pre variance
    lntec = post_var / np.maximum(pre_var, 1e-8)

    # DCR: we approximate displacement in PCA space
    # Since embeddings have different dims, use the post_cnem as quality proxy
    # DCR = cnem_delta * (1 / neighborhood_quality)
    # neighborhood_quality = 1 - post_cnem (higher = better match)
    nbr_quality = 1.0 - post_cnem
    dcr = cnem_delta / np.maximum(nbr_quality, 0.01)

    return pd.DataFrame({
        'pre_cnem': pre_cnem,
        'post_cnem': post_cnem,
        'cnem_delta': cnem_delta,
        'lntec': lntec,
        'dcr': dcr,
    })

def compute_metrics_from_graph(conn_matrix, X_expr, pre_cnem, pre_var, k=K):
    """Compute metrics from connectivity graph (BBKNN)."""
    n = conn_matrix.shape[0]
    post_cnem = np.zeros(n)
    post_var = np.zeros(n)

    for i in range(n):
        row = conn_matrix[i].toarray().flatten()
        nbr_idx = np.argsort(row)[-k:]
        nbr_idx = nbr_idx[row[nbr_idx] > 0]
        if len(nbr_idx) == 0:
            post_cnem[i] = pre_cnem[i]
            post_var[i] = pre_var[i]
            continue
        cell_e = X_expr[i]
        nbr_e = X_expr[nbr_idx]
        mean_nbr = nbr_e.mean(axis=0, keepdims=True)
        post_cnem[i] = 1.0 - cosine_similarity(cell_e.reshape(1, -1), mean_nbr)[0, 0]
        dists = 1.0 - cosine_similarity(cell_e.reshape(1, -1), nbr_e)[0]
        post_var[i] = np.var(dists)

    cnem_delta = post_cnem - pre_cnem
    lntec = post_var / np.maximum(pre_var, 1e-8)
    nbr_quality = 1.0 - post_cnem
    dcr = cnem_delta / np.maximum(nbr_quality, 0.01)

    return pd.DataFrame({
        'pre_cnem': pre_cnem,
        'post_cnem': post_cnem,
        'cnem_delta': cnem_delta,
        'lntec': lntec,
        'dcr': dcr,
    })

results = {}

# --- Harmony ---
print("\n--- Harmony ---")
t0 = time.time()
import harmonypy as hm
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'tech')
Z = np.array(ho.Z_corr)
if Z.shape[0] != adata.n_obs:
    Z = Z.T
metrics = compute_metrics(Z, X_expr, pre_cnem, pre_var, adata.obsm['X_pca'])
results['Harmony'] = metrics
print(f"  Done in {time.time()-t0:.1f}s")

# --- Scanorama ---
print("\n--- Scanorama ---")
t0 = time.time()
import scanorama
adata_s = adata.copy()
adata_s = adata_s[:, adata_s.var['highly_variable']].copy()
sc.pp.scale(adata_s, max_value=10)
sc.pp.pca(adata_s, n_comps=50)
batch_cats = list(adata_s.obs['tech'].cat.categories) if hasattr(adata_s.obs['tech'], 'cat') else list(adata_s.obs['tech'].unique())
adatas_list = [adata_s[adata_s.obs['tech'] == b].copy() for b in batch_cats]
scanorama.integrate_scanpy(adatas_list)
adata_sc = ad.concat(adatas_list)
adata_sc = adata_sc[adata.obs_names].copy()
metrics = compute_metrics(adata_sc.obsm['X_scanorama'], X_expr, pre_cnem, pre_var, adata.obsm['X_pca'])
results['Scanorama'] = metrics
print(f"  Done in {time.time()-t0:.1f}s")

# --- BBKNN ---
print("\n--- BBKNN ---")
t0 = time.time()
import bbknn
adata_b = adata.copy()
adata_b = adata_b[:, adata_b.var['highly_variable']].copy()
sc.pp.scale(adata_b, max_value=10)
sc.pp.pca(adata_b, n_comps=50)
bbknn.bbknn(adata_b, batch_key='tech', n_pcs=50)
metrics = compute_metrics_from_graph(adata_b.obsp['connectivities'], X_expr, pre_cnem, pre_var)
results['BBKNN'] = metrics
print(f"  Done in {time.time()-t0:.1f}s")

# ============================================================
# 4. Evaluation
# ============================================================
print("\n" + "=" * 60)
print("Step 4: Evaluation")
print("=" * 60)

celltypes = adata.obs['celltype_true'].values
is_rare = adata.obs['is_rare'].values

for method_name, metrics_df in results.items():
    print(f"\n{'='*40}")
    print(f"  {method_name}")
    print(f"{'='*40}")

    metrics_df['celltype'] = celltypes
    metrics_df['is_rare'] = is_rare

    for metric in ['post_cnem', 'cnem_delta', 'lntec', 'dcr']:
        print(f"\n--- {metric} ---")

        # Summary by celltype
        summary = metrics_df.groupby('celltype')[metric].agg(['mean', 'median', 'count'])
        summary = summary.sort_values('mean', ascending=False)
        print(summary.round(4))

        # Wilcoxon test
        rare_vals = metrics_df[metrics_df['is_rare'] == 1][metric].dropna()
        common_vals = metrics_df[metrics_df['is_rare'] == 0][metric].dropna()
        if len(rare_vals) > 0 and len(common_vals) > 0:
            stat, pval = stats.mannwhitneyu(rare_vals, common_vals, alternative='greater')
            print(f"Wilcoxon rare>common: p={pval:.2e}")

        # ROC-AUC for rare detection
        valid = ~np.isnan(metrics_df[metric]) & ~np.isinf(metrics_df[metric])
        if valid.sum() > 100:
            try:
                auc = roc_auc_score(metrics_df.loc[valid, 'is_rare'], metrics_df.loc[valid, metric])
                print(f"ROC-AUC (rare detection): {auc:.3f}")
            except:
                pass

        # Top 5% enrichment
        threshold = metrics_df[metric].quantile(0.95)
        top5 = metrics_df[metrics_df[metric] >= threshold]
        if len(top5) > 0:
            rare_in_top5 = top5['is_rare'].sum()
            rare_frac = is_rare.sum() / len(is_rare)
            enrichment = (rare_in_top5 / len(top5)) / rare_frac if rare_frac > 0 else 0
            print(f"Top 5%: {rare_in_top5}/{len(top5)} rare, OR={enrichment:.2f}")

            # Epsilon-specific
            eps_in_top5 = (top5['celltype'] == 'epsilon').sum()
            eps_total = (celltypes == 'epsilon').sum()
            print(f"  Epsilon: {eps_in_top5}/{eps_total} in top 5%")

        # Top 10%
        threshold_10 = metrics_df[metric].quantile(0.90)
        top10 = metrics_df[metrics_df[metric] >= threshold_10]
        if len(top10) > 0:
            eps_in_top10 = (top10['celltype'] == 'epsilon').sum()
            rare_in_top10 = top10['is_rare'].sum()
            rare_frac = is_rare.sum() / len(is_rare)
            enrichment_10 = (rare_in_top10 / len(top10)) / rare_frac if rare_frac > 0 else 0
            print(f"Top 10%: epsilon {eps_in_top10}/{eps_total}, OR={enrichment_10:.2f}")

# ============================================================
# 5. Cross-method summary
# ============================================================
print("\n" + "=" * 60)
print("Step 5: Cross-method CNEM summary")
print("=" * 60)

all_dfs = []
for method_name, metrics_df in results.items():
    df = metrics_df[['celltype', 'post_cnem', 'cnem_delta']].copy()
    df['method'] = method_name
    all_dfs.append(df)
combined = pd.concat(all_dfs)

print("\n=== post_cnem (absolute mismatch after integration) ===")
pivot = combined.groupby(['celltype', 'method'])['post_cnem'].mean().unstack()
pivot['is_rare'] = pivot.index.isin(rare_types)
pivot = pivot.sort_values('is_rare', ascending=False)
print(pivot.round(4).to_string())

print("\n=== cnem_delta (mismatch change: positive = worse after integration) ===")
pivot2 = combined.groupby(['celltype', 'method'])['cnem_delta'].mean().unstack()
pivot2['is_rare'] = pivot2.index.isin(rare_types)
pivot2 = pivot2.sort_values('is_rare', ascending=False)
print(pivot2.round(4).to_string())

# Save
for method_name, metrics_df in results.items():
    metrics_df.to_csv(f'{OUT_DIR}/cnem_scores_{method_name.lower()}.csv', index=False)
combined.to_csv(f'{OUT_DIR}/cnem_summary_all.csv', index=False)

print("\nPoC v3 complete!")
