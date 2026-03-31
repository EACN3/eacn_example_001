"""Stratified detection: per-celltype subclustering + R(S) + survival analysis.
Matches data science methodology: subcluster within each celltype, then
check which subclusters are disrupted by integration.

Also searches for alternative rare targets: cDC3/AS-DC, Tpex.
GPU-accelerated.
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import adjusted_rand_score
import torch
import faiss
import warnings, os, time, gc
warnings.filterwarnings('ignore')

OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/stratified/results'
os.makedirs(OUT_DIR, exist_ok=True)

K = 30
GPU_ID = 0

def faiss_knn_gpu(X, k, gpu_id=GPU_ID):
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    index = faiss.index_cpu_to_gpu(res, gpu_id, faiss.IndexFlatL2(X.shape[1]))
    index.add(X)
    _, indices = index.search(X, k + 1)
    return indices[:, 1:]

def compute_RS(X_pre, X_post, pre_indices):
    n = X_post.shape[0]
    DS = np.zeros(n, dtype=np.float32)
    for i in range(n):
        nbr = pre_indices[i]
        nbr_post = X_post[nbr]
        mean_p = nbr_post.mean(axis=0)
        DS[i] = np.mean(np.sum((nbr_post - mean_p) ** 2, axis=1))
    D_ref = np.median(DS)
    return DS / max(D_ref, 1e-8)

# ============================================================
# 1. Load the Phase 1 subset (same as data science used)
# ============================================================
print("=" * 60)
print("Step 1: Load subset")
print("=" * 60)
t_total = time.time()

cache = '/ssd/data/agent/bio/eacn_example_001/experiments/phase1/data/immune_subset_phase1.h5ad'
adata = ad.read_h5ad(cache)
print(f"Loaded: {adata.shape}")

# Preprocess
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ============================================================
# 2. Global HVG + PCA + Harmony (for post-integration embedding)
# ============================================================
print("\n" + "=" * 60)
print("Step 2: Global integration (Harmony)")
print("=" * 60)

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)

import harmonypy as hm
ho = hm.run_harmony(adata_hvg.obsm['X_pca'], adata.obs, 'Dataset')
Z = np.array(ho.Z_corr)
if Z.shape[0] != adata.n_obs:
    Z = Z.T
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.obsm['X_harmony'] = Z
print(f"Harmony done: {Z.shape}")

# ============================================================
# 3. Per-celltype stratified analysis
# ============================================================
print("\n" + "=" * 60)
print("Step 3: Stratified per-celltype subclustering + survival")
print("=" * 60)

celltypes = adata.obs['Celltype'].unique()
all_results = []

for ct in celltypes:
    ct_mask = adata.obs['Celltype'] == ct
    ct_idx = np.where(ct_mask)[0]
    n_ct = len(ct_idx)
    if n_ct < 50:
        continue

    print(f"\n--- {ct} ({n_ct} cells) ---")

    # Extract celltype-specific PCA and Harmony embeddings
    X_pca_ct = adata.obsm['X_pca'][ct_idx]
    X_harm_ct = adata.obsm['X_harmony'][ct_idx]

    # Pre-integration subclustering (within celltype)
    adata_ct = ad.AnnData(obs=adata.obs.iloc[ct_idx].copy())
    adata_ct.obsm['X_pca'] = X_pca_ct

    k_use = min(K, n_ct - 1)
    sc.pp.neighbors(adata_ct, use_rep='X_pca', n_neighbors=k_use)
    sc.tl.leiden(adata_ct, resolution=1.0, key_added='pre_sub')
    n_pre = adata_ct.obs['pre_sub'].nunique()

    # Post-integration subclustering (within celltype, using Harmony)
    adata_ct.obsm['X_harmony'] = X_harm_ct
    sc.pp.neighbors(adata_ct, use_rep='X_harmony', n_neighbors=k_use)
    sc.tl.leiden(adata_ct, resolution=1.0, key_added='post_sub')
    n_post = adata_ct.obs['post_sub'].nunique()

    # ARI
    ari = adjusted_rand_score(adata_ct.obs['pre_sub'], adata_ct.obs['post_sub'])

    # Per-subcluster survival
    pre_labels = adata_ct.obs['pre_sub'].values
    post_labels = adata_ct.obs['post_sub'].values
    pre_counts = pd.Series(pre_labels).value_counts()

    n_dispersed = 0
    for cl_id, cl_size in pre_counts.items():
        cl_local_idx = np.where(pre_labels == cl_id)[0]
        post_assign = post_labels[cl_local_idx]
        post_counts = pd.Series(post_assign).value_counts()
        survival = post_counts.iloc[0] / cl_size
        n_frags = len(post_counts)

        if survival < 0.5:
            n_dispersed += 1

        all_results.append({
            'celltype': ct,
            'pre_cluster': cl_id,
            'size': cl_size,
            'survival': survival,
            'n_fragments': n_frags,
        })

    dispersed_pct = n_dispersed / n_pre * 100 if n_pre > 0 else 0
    print(f"  Pre-clusters: {n_pre}, Post-clusters: {n_post}")
    print(f"  ARI: {ari:.3f}")
    print(f"  Dispersed (<50% survival): {n_dispersed}/{n_pre} ({dispersed_pct:.0f}%)")

# ============================================================
# 4. R(S) per celltype (stratified)
# ============================================================
print("\n" + "=" * 60)
print("Step 4: Stratified R(S)")
print("=" * 60)

rs_all = np.zeros(adata.n_obs, dtype=np.float32)

for ct in celltypes:
    ct_mask = adata.obs['Celltype'] == ct
    ct_idx = np.where(ct_mask)[0]
    n_ct = len(ct_idx)
    if n_ct < 50:
        continue

    X_pca_ct = adata.obsm['X_pca'][ct_idx]
    X_harm_ct = adata.obsm['X_harmony'][ct_idx]

    k_use = min(K, n_ct - 1)
    pre_idx_ct = faiss_knn_gpu(X_pca_ct, k_use)
    rs_ct = compute_RS(X_pca_ct, X_harm_ct, pre_idx_ct)
    rs_all[ct_idx] = rs_ct

    print(f"  {ct}: R(S) mean={rs_ct.mean():.3f}, median={np.median(rs_ct):.3f}, >3: {(rs_ct>3).sum()}")

# ============================================================
# 5. Search for alternative rare targets
# ============================================================
print("\n" + "=" * 60)
print("Step 5: Alternative rare target search")
print("=" * 60)

# Check marker genes
alt_markers = {
    'cDC3_ASDC': ['AXL', 'SIGLEC6'],
    'Tpex': ['TCF7', 'PDCD1', 'CD8A', 'CXCR5'],
}

for target, markers in alt_markers.items():
    print(f"\n--- {target} ---")
    found = [m for m in markers if m in adata.var_names]
    missing = [m for m in markers if m not in adata.var_names]
    print(f"  Found markers: {found}, Missing: {missing}")

    if len(found) >= 2:
        # Check co-expression
        expr_data = {}
        for m in found:
            e = adata[:, m].X
            if hasattr(e, 'toarray'):
                e = e.toarray()
            expr_data[m] = e.flatten()

        if target == 'cDC3_ASDC' and 'AXL' in found and 'SIGLEC6' in found:
            mask = (expr_data['AXL'] > 0) & (expr_data['SIGLEC6'] > 0)
            dc_mask = adata.obs['Celltype'] == 'Dendritic cell'
            target_cells = mask & dc_mask.values
            print(f"  AXL+SIGLEC6+ in DC: {target_cells.sum()} cells")
            if target_cells.sum() > 0:
                print(f"  Datasets: {adata.obs.loc[target_cells, 'Dataset'].value_counts().to_dict()}")
                # R(S) for these cells
                target_rs = rs_all[target_cells]
                print(f"  R(S): mean={target_rs.mean():.3f}")

        if target == 'Tpex' and 'TCF7' in found and 'PDCD1' in found:
            t_mask = adata.obs['Celltype'] == 'T cell'
            cd8_expr = expr_data.get('CD8A', np.zeros(adata.n_obs))
            mask = (expr_data['TCF7'] > 0) & (expr_data['PDCD1'] > 0) & (cd8_expr > 0) & t_mask.values
            print(f"  TCF7+PD1+CD8A+ in T cell: {mask.sum()} cells")
            if mask.sum() > 0:
                print(f"  Datasets: {adata.obs.loc[mask, 'Dataset'].value_counts().head(5).to_dict()}")
                target_rs = rs_all[mask]
                print(f"  R(S): mean={target_rs.mean():.3f}")

# ============================================================
# 6. Save
# ============================================================
print("\n" + "=" * 60)
print("Step 6: Save")
print("=" * 60)

results_df = pd.DataFrame(all_results)
results_df.to_csv(f'{OUT_DIR}/stratified_survival.csv', index=False)

# Summary
summary = results_df.groupby('celltype').agg(
    n_clusters=('pre_cluster', 'count'),
    n_dispersed=('survival', lambda x: (x < 0.5).sum()),
    mean_survival=('survival', 'mean'),
    min_survival=('survival', 'min'),
).reset_index()
summary['dispersed_pct'] = summary['n_dispersed'] / summary['n_clusters'] * 100
summary = summary.sort_values('dispersed_pct', ascending=False)
print("\n=== Stratified survival summary ===")
print(summary.to_string(index=False))

# Most disrupted subclusters
print("\n=== Top 10 most disrupted subclusters ===")
print(results_df.nsmallest(10, 'survival').to_string(index=False))

rs_df = pd.DataFrame({'RS': rs_all, 'celltype': adata.obs['Celltype'].values})
rs_df.to_csv(f'{OUT_DIR}/stratified_rs.csv', index=False)

print(f"\nTotal time: {time.time()-t_total:.0f}s")
print("Stratified detection complete!")
