"""PoC Experiment: Detect rare subpopulation disruption in batch integration
using neighborhood structure comparison (convergence vs dispersion)."""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.neighbors import NearestNeighbors
import warnings, os, time
warnings.filterwarnings('ignore')

DATA_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/poc_v1/data'
OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/poc_v1/results'
os.makedirs(OUT_DIR, exist_ok=True)

K = 30  # k for kNN

# ============================================================
# 1. Load and preprocess
# ============================================================
print("=" * 60)
print("Step 1: Load and preprocess")
print("=" * 60)
adata = ad.read_h5ad(f'{DATA_DIR}/pancreas.h5ad')
print(f"Loaded: {adata.shape}")

# Store raw cell type labels
adata.obs['celltype_true'] = adata.obs['celltype'].copy()
adata.obs['is_rare'] = adata.obs['celltype'].isin(['epsilon', 'schwann', 't_cell', 'mast']).astype(int)

# Standard preprocessing
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='tech')
adata.raw = adata
adata = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata, n_comps=50)
print(f"After preprocessing: {adata.shape}")

# ============================================================
# 2. Compute pre-integration per-batch kNN
# ============================================================
print("\n" + "=" * 60)
print("Step 2: Pre-integration per-batch kNN")
print("=" * 60)

batches = adata.obs['tech'].unique()
pre_neighbors = {}  # cell_idx -> set of neighbor indices (global)

for batch in batches:
    mask = adata.obs['tech'] == batch
    idx = np.where(mask)[0]
    X_batch = adata.obsm['X_pca'][idx]
    
    k_use = min(K, len(idx) - 1)
    nn = NearestNeighbors(n_neighbors=k_use + 1, metric='euclidean')
    nn.fit(X_batch)
    distances, indices = nn.kneighbors(X_batch)
    
    for i, global_i in enumerate(idx):
        # indices are local, convert to global
        neighbors_global = set(idx[indices[i, 1:]])  # exclude self
        pre_neighbors[global_i] = neighbors_global
    
    print(f"  Batch {batch}: {len(idx)} cells, k={k_use}")

print(f"Pre-integration neighbors computed for {len(pre_neighbors)} cells")

# ============================================================
# 3. Run integration methods and compute post-integration kNN
# ============================================================
print("\n" + "=" * 60)
print("Step 3: Run integration methods")
print("=" * 60)

def compute_post_neighbors(adata_int, embedding_key, k=K):
    """Compute post-integration kNN from embedding."""
    X = adata_int.obsm[embedding_key]
    nn = NearestNeighbors(n_neighbors=k + 1, metric='euclidean')
    nn.fit(X)
    distances, indices = nn.kneighbors(X)
    post_nbrs = {}
    for i in range(X.shape[0]):
        post_nbrs[i] = set(indices[i, 1:])
    return post_nbrs, distances[:, 1:]

results = {}

# --- Harmony ---
print("\n--- Harmony ---")
t0 = time.time()
import harmonypy as hm
adata_h = adata.copy()
ho = hm.run_harmony(adata_h.obsm['X_pca'], adata_h.obs, 'tech')
Z = np.array(ho.Z_corr)
if Z.shape[0] != adata_h.n_obs:
    Z = Z.T
adata_h.obsm['X_harmony'] = Z
post_nbrs_h, post_dists_h = compute_post_neighbors(adata_h, 'X_harmony')
results['Harmony'] = {'post_nbrs': post_nbrs_h, 'post_dists': post_dists_h, 'adata': adata_h}
print(f"  Done in {time.time()-t0:.1f}s")

# --- Scanorama ---
print("\n--- Scanorama ---")
t0 = time.time()
import scanorama
adata_s = adata.copy()
batch_cats = adata_s.obs['tech'].cat.categories if hasattr(adata_s.obs['tech'], 'cat') else adata_s.obs['tech'].unique()
adatas_list = [adata_s[adata_s.obs['tech'] == b].copy() for b in batch_cats]
scanorama.integrate_scanpy(adatas_list)
adata_s = ad.concat(adatas_list)
# Reorder to match original
adata_s = adata_s[adata.obs_names].copy()
post_nbrs_s, post_dists_s = compute_post_neighbors(adata_s, 'X_scanorama')
results['Scanorama'] = {'post_nbrs': post_nbrs_s, 'post_dists': post_dists_s, 'adata': adata_s}
print(f"  Done in {time.time()-t0:.1f}s")

# --- BBKNN ---
print("\n--- BBKNN ---")
t0 = time.time()
import bbknn
adata_b = adata.copy()
bbknn.bbknn(adata_b, batch_key='tech', n_pcs=50)
# BBKNN modifies connectivities directly; extract neighbors from the graph
conn = adata_b.obsp['connectivities']
post_nbrs_b = {}
post_dists_b_list = []
for i in range(conn.shape[0]):
    row = conn[i].toarray().flatten()
    # Get top-K neighbors by connectivity weight
    top_k = np.argsort(row)[-K:]
    top_k = top_k[row[top_k] > 0]
    post_nbrs_b[i] = set(top_k)
    # Use 1-connectivity as pseudo-distance
    if len(top_k) > 0:
        dists = 1.0 - row[top_k]
    else:
        dists = np.array([])
    post_dists_b_list.append(dists)
results['BBKNN'] = {'post_nbrs': post_nbrs_b, 'post_dists': post_dists_b_list, 'adata': adata_b}
print(f"  Done in {time.time()-t0:.1f}s")

# ============================================================
# 4. Compute disruption scores
# ============================================================
print("\n" + "=" * 60)
print("Step 4: Compute disruption scores")
print("=" * 60)

batch_labels = adata.obs['tech'].values

def compute_disruption_scores(pre_nbrs, post_nbrs, post_dists, batch_labels, k=K):
    """Compute NPS and Disruption score for each cell."""
    n_cells = len(pre_nbrs)
    nps = np.zeros(n_cells)
    neighbor_dispersion = np.zeros(n_cells)
    cbr = np.zeros(n_cells)
    
    for i in range(n_cells):
        pre_set = pre_nbrs.get(i, set())
        post_set = post_nbrs.get(i, set())
        
        if len(pre_set) == 0 or len(post_set) == 0:
            continue
        
        # NPS: overlap of same-batch neighbors
        my_batch = batch_labels[i]
        post_same_batch = {j for j in post_set if batch_labels[j] == my_batch}
        
        if len(post_same_batch) > 0:
            nps[i] = len(pre_set & post_same_batch) / min(k, len(post_same_batch))
        
        # Cross-batch ratio
        post_other = {j for j in post_set if batch_labels[j] != my_batch}
        cbr[i] = len(post_other) / max(len(post_set), 1)
        
        # Neighbor dispersion: variance of distances to post-integration neighbors
        if isinstance(post_dists, np.ndarray) and i < len(post_dists):
            dists_i = post_dists[i]
            if len(dists_i) > 1:
                neighbor_dispersion[i] = np.std(dists_i)
        elif isinstance(post_dists, list) and i < len(post_dists):
            dists_i = post_dists[i]
            if len(dists_i) > 1:
                neighbor_dispersion[i] = np.std(dists_i)
    
    # Normalize dispersion to [0,1]
    if neighbor_dispersion.max() > 0:
        neighbor_dispersion = neighbor_dispersion / neighbor_dispersion.max()
    
    # Disruption = (1 - NPS) * NeighborDispersion
    disruption = (1 - nps) * neighbor_dispersion
    
    return pd.DataFrame({
        'NPS': nps,
        'CBR': cbr,
        'NeighborDispersion': neighbor_dispersion,
        'Disruption': disruption
    })

for method_name, res in results.items():
    print(f"\n--- {method_name} ---")
    scores = compute_disruption_scores(
        pre_neighbors, res['post_nbrs'], res['post_dists'], 
        batch_labels, k=K
    )
    scores['celltype'] = adata.obs['celltype_true'].values
    scores['is_rare'] = adata.obs['is_rare'].values
    scores['batch'] = batch_labels
    
    results[method_name]['scores'] = scores
    
    # Summary by cell type
    summary = scores.groupby('celltype')['Disruption'].agg(['mean', 'median', 'std', 'count'])
    summary = summary.sort_values('mean', ascending=False)
    print(summary)
    
    # Statistical test: rare vs common
    rare_scores = scores[scores['is_rare'] == 1]['Disruption']
    common_scores = scores[scores['is_rare'] == 0]['Disruption']
    stat, pval = stats.mannwhitneyu(rare_scores, common_scores, alternative='greater')
    print(f"\nWilcoxon rare vs common: U={stat:.0f}, p={pval:.2e}")
    
    # Enrichment in top 5%
    threshold = scores['Disruption'].quantile(0.95)
    top5 = scores[scores['Disruption'] >= threshold]
    rare_in_top5 = top5['is_rare'].sum()
    rare_total = scores['is_rare'].sum()
    enrichment = (rare_in_top5 / len(top5)) / (rare_total / len(scores))
    print(f"Top 5% enrichment OR for rare types: {enrichment:.2f}")
    print(f"  Rare in top 5%: {rare_in_top5}/{len(top5)}")
    print(f"  Top 5% cell types: {top5['celltype'].value_counts().to_dict()}")

# ============================================================
# 5. Blind test: hide epsilon labels
# ============================================================
print("\n" + "=" * 60)
print("Step 5: Blind test (epsilon hidden)")
print("=" * 60)

for method_name, res in results.items():
    scores = res['scores']
    threshold = scores['Disruption'].quantile(0.95)
    top5_mask = scores['Disruption'] >= threshold
    
    # How many epsilon cells are in top 5%?
    epsilon_mask = scores['celltype'] == 'epsilon'
    epsilon_in_top5 = (top5_mask & epsilon_mask).sum()
    epsilon_total = epsilon_mask.sum()
    recall = epsilon_in_top5 / epsilon_total if epsilon_total > 0 else 0
    
    print(f"\n{method_name}: epsilon recall in top 5% = {epsilon_in_top5}/{epsilon_total} = {recall:.1%}")

# ============================================================
# 6. Save results
# ============================================================
print("\n" + "=" * 60)
print("Step 6: Save results")
print("=" * 60)

for method_name, res in results.items():
    res['scores'].to_csv(f'{OUT_DIR}/disruption_scores_{method_name.lower()}.csv', index=False)
    print(f"Saved {method_name} scores")

# Combined summary
all_scores = []
for method_name, res in results.items():
    df = res['scores'].copy()
    df['method'] = method_name
    all_scores.append(df)
combined = pd.concat(all_scores)
combined.to_csv(f'{OUT_DIR}/disruption_scores_all.csv', index=False)

# Summary table
print("\n=== FINAL SUMMARY ===")
summary = combined.groupby(['method', 'celltype'])['Disruption'].mean().unstack(level=0)
summary['is_rare'] = summary.index.isin(['epsilon', 'schwann', 't_cell', 'mast'])
summary = summary.sort_values('is_rare', ascending=False)
print(summary.to_string())

print("\nPoC experiment complete!")
