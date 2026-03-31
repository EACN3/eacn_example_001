"""PoC v2: Improved disruption detection via expression coherence.

Key insight: Instead of comparing neighbor identity (NPS), compare the
EXPRESSION COHERENCE of post-integration neighborhoods.

Normal integration (convergence): cells gain same-type cross-batch neighbors
-> neighborhood remains expression-coherent.

Rare subgroup destruction (dispersion): rare cells get pushed into
neighborhoods of different cell types -> expression-INCOHERENT neighborhood.

Metric: coherence_drop = pre_coherence - post_coherence
High coherence_drop = potential disruption.
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics.pairwise import cosine_similarity
import warnings, os, time
warnings.filterwarnings('ignore')

DATA_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/poc_v1/data'
OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/poc_v2/results'
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
adata.obs['is_rare'] = adata.obs['celltype'].isin(['epsilon', 'schwann', 't_cell', 'mast']).astype(int)

sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='tech')
adata.raw = adata

# Keep normalized expression for coherence computation
X_expr = adata[:, adata.var['highly_variable']].X.copy()
if hasattr(X_expr, 'toarray'):
    X_expr = X_expr.toarray()
X_expr = np.array(X_expr, dtype=np.float32)

adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

print(f"Loaded: {adata.shape}, HVG expression matrix: {X_expr.shape}")

# ============================================================
# 2. Pre-integration per-batch neighborhood coherence
# ============================================================
print("\n" + "=" * 60)
print("Step 2: Pre-integration per-batch neighborhood coherence")
print("=" * 60)

batches = adata.obs['tech'].unique()
pre_coherence = np.zeros(adata.n_obs)

for batch in batches:
    mask = adata.obs['tech'] == batch
    idx = np.where(mask)[0]
    X_pca_batch = adata.obsm['X_pca'][idx]
    X_expr_batch = X_expr[idx]

    k_use = min(K, len(idx) - 1)
    nn = NearestNeighbors(n_neighbors=k_use + 1, metric='euclidean')
    nn.fit(X_pca_batch)
    _, indices = nn.kneighbors(X_pca_batch)

    for i_local, i_global in enumerate(idx):
        nbr_local = indices[i_local, 1:]
        cell_expr = X_expr_batch[i_local:i_local+1]
        nbr_expr = X_expr_batch[nbr_local]
        sims = cosine_similarity(cell_expr, nbr_expr)[0]
        pre_coherence[i_global] = np.mean(sims)

print(f"Pre-integration coherence: mean={pre_coherence.mean():.4f}, std={pre_coherence.std():.4f}")

# ============================================================
# 3. Run integration methods
# ============================================================
print("\n" + "=" * 60)
print("Step 3: Run integration methods")
print("=" * 60)

def compute_post_coherence(embedding, X_expr, k=K):
    nn = NearestNeighbors(n_neighbors=k + 1, metric='euclidean')
    nn.fit(embedding)
    _, indices = nn.kneighbors(embedding)
    post_coh = np.zeros(embedding.shape[0])
    for i in range(embedding.shape[0]):
        nbr_idx = indices[i, 1:]
        cell_expr = X_expr[i:i+1]
        nbr_expr = X_expr[nbr_idx]
        sims = cosine_similarity(cell_expr, nbr_expr)[0]
        post_coh[i] = np.mean(sims)
    return post_coh

def compute_post_coherence_from_graph(conn_matrix, X_expr, k=K):
    n = conn_matrix.shape[0]
    post_coh = np.zeros(n)
    for i in range(n):
        row = conn_matrix[i].toarray().flatten()
        nbr_idx = np.argsort(row)[-k:]
        nbr_idx = nbr_idx[row[nbr_idx] > 0]
        if len(nbr_idx) == 0:
            continue
        cell_expr = X_expr[i:i+1]
        nbr_expr = X_expr[nbr_idx]
        sims = cosine_similarity(cell_expr, nbr_expr)[0]
        post_coh[i] = np.mean(sims)
    return post_coh

results = {}

# --- Harmony ---
print("\n--- Harmony ---")
t0 = time.time()
import harmonypy as hm
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'tech')
Z = np.array(ho.Z_corr)
if Z.shape[0] != adata.n_obs:
    Z = Z.T
post_coh = compute_post_coherence(Z, X_expr)
results['Harmony'] = {'post_coherence': post_coh, 'embedding': Z}
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
post_coh = compute_post_coherence(adata_sc.obsm['X_scanorama'], X_expr)
results['Scanorama'] = {'post_coherence': post_coh, 'embedding': adata_sc.obsm['X_scanorama']}
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
post_coh = compute_post_coherence_from_graph(adata_b.obsp['connectivities'], X_expr)
results['BBKNN'] = {'post_coherence': post_coh}
print(f"  Done in {time.time()-t0:.1f}s")

# ============================================================
# 4. Compute coherence drop scores
# ============================================================
print("\n" + "=" * 60)
print("Step 4: Compute coherence drop scores")
print("=" * 60)

for method_name, res in results.items():
    print(f"\n--- {method_name} ---")
    post_coh = res['post_coherence']
    coherence_drop = pre_coherence - post_coh
    coherence_ratio = post_coh / np.maximum(pre_coherence, 1e-8)

    scores = pd.DataFrame({
        'pre_coherence': pre_coherence,
        'post_coherence': post_coh,
        'coherence_drop': coherence_drop,
        'coherence_ratio': coherence_ratio,
        'celltype': adata.obs['celltype_true'].values,
        'is_rare': adata.obs['is_rare'].values,
        'batch': adata.obs['tech'].values,
    })
    res['scores'] = scores

    summary = scores.groupby('celltype')[['coherence_drop', 'post_coherence']].agg(['mean', 'median'])
    summary.columns = ['_'.join(c) for c in summary.columns]
    summary = summary.sort_values('coherence_drop_mean', ascending=False)
    print(summary.round(4))

    rare_scores = scores[scores['is_rare'] == 1]['coherence_drop']
    common_scores = scores[scores['is_rare'] == 0]['coherence_drop']
    stat, pval = stats.mannwhitneyu(rare_scores, common_scores, alternative='greater')
    print(f"\nWilcoxon rare vs common: U={stat:.0f}, p={pval:.2e}")

    threshold = scores['coherence_drop'].quantile(0.95)
    top5 = scores[scores['coherence_drop'] >= threshold]
    rare_in_top5 = top5['is_rare'].sum()
    rare_total = scores['is_rare'].sum()
    enrichment = (rare_in_top5 / len(top5)) / (rare_total / len(scores)) if len(top5) > 0 else 0
    print(f"Top 5% enrichment OR: {enrichment:.2f}")
    print(f"  Rare in top 5%: {rare_in_top5}/{len(top5)}")
    print(f"  Types in top 5%: {top5['celltype'].value_counts().head(10).to_dict()}")

    epsilon_mask = scores['celltype'] == 'epsilon'
    epsilon_in_top5 = (scores['coherence_drop'] >= threshold)[epsilon_mask].sum()
    print(f"\n  Epsilon: {epsilon_in_top5}/{epsilon_mask.sum()} in top 5%")

    threshold_10 = scores['coherence_drop'].quantile(0.90)
    epsilon_in_top10 = (scores['coherence_drop'] >= threshold_10)[epsilon_mask].sum()
    top10 = scores[scores['coherence_drop'] >= threshold_10]
    rare_in_top10 = top10['is_rare'].sum()
    enrichment_10 = (rare_in_top10 / len(top10)) / (rare_total / len(scores)) if len(top10) > 0 else 0
    print(f"  Top 10%: epsilon {epsilon_in_top10}/{epsilon_mask.sum()}, OR={enrichment_10:.2f}")

# ============================================================
# 5. Save
# ============================================================
print("\n" + "=" * 60)
print("Step 5: Save results")
print("=" * 60)

for method_name, res in results.items():
    res['scores'].to_csv(f'{OUT_DIR}/coherence_scores_{method_name.lower()}.csv', index=False)

all_scores = []
for method_name, res in results.items():
    df = res['scores'].copy()
    df['method'] = method_name
    all_scores.append(df)
combined = pd.concat(all_scores)
combined.to_csv(f'{OUT_DIR}/coherence_scores_all.csv', index=False)

print("\n=== CROSS-METHOD SUMMARY (mean coherence_drop) ===")
pivot = combined.groupby(['celltype', 'method'])['coherence_drop'].mean().unstack()
pivot['is_rare'] = pivot.index.isin(['epsilon', 'schwann', 't_cell', 'mast'])
pivot = pivot.sort_values('is_rare', ascending=False)
print(pivot.round(4).to_string())

print("\nPoC v2 complete!")
