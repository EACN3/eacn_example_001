"""Phase 1b: Comprehensive detection with NHL batches.
Implements R(S) dispersion score + stratified CNEM + small cluster survival scan.
All GPU-accelerated with chunked computation to avoid OOM.
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
import torch
import faiss
import warnings, os, time, gc
warnings.filterwarnings('ignore')

OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/phase1b/results'
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

def compute_cnem_chunked(X_expr, indices, chunk=5000):
    """CPU-chunked CNEM to avoid GPU OOM."""
    n = X_expr.shape[0]
    cnem = np.zeros(n, dtype=np.float32)
    for s in range(0, n, chunk):
        e = min(s + chunk, n)
        X_chunk = torch.from_numpy(X_expr[s:e]).cuda(GPU_ID)
        idx_chunk = indices[s:e]
        # Gather neighbor expressions on CPU, compute mean, then GPU cosine
        nbr_mean = np.zeros((e - s, X_expr.shape[1]), dtype=np.float32)
        for i in range(e - s):
            nbr_mean[i] = X_expr[idx_chunk[i]].mean(axis=0)
        nbr_mean_t = torch.from_numpy(nbr_mean).cuda(GPU_ID)
        cos = torch.nn.functional.cosine_similarity(X_chunk, nbr_mean_t, dim=1)
        cnem[s:e] = (1.0 - cos).cpu().numpy()
        del X_chunk, nbr_mean_t
    return cnem

def compute_RS(X_embed_pre, X_embed_post, pre_indices):
    """R(S) dispersion score: measures how much pre-integration neighborhoods
    spread out after integration. Works in embedding space only."""
    n = X_embed_post.shape[0]
    DS = np.zeros(n, dtype=np.float32)
    for i in range(n):
        nbr_idx = pre_indices[i]
        nbr_post = X_embed_post[nbr_idx]
        mean_post = nbr_post.mean(axis=0)
        DS[i] = np.mean(np.sum((nbr_post - mean_post) ** 2, axis=1))
    D_ref = np.median(DS)
    RS = DS / max(D_ref, 1e-8)
    return RS, DS

# ============================================================
# 1. Extract subset with NHL batches
# ============================================================
print("=" * 60)
print("Step 1: Extract subset with NHL + non-NHL batches")
print("=" * 60)
t_total = time.time()

# Selected: 3 NHL (D61,D63,D65) + 5 non-NHL with pDC (D5,D100,D38,D88,D83)
SELECTED = ['D61', 'D63', 'D65', 'D5', 'D100', 'D38', 'D88', 'D83']

cache_path = '/ssd/data/agent/bio/eacn_example_001/experiments/phase1b/data/subset_nhl.h5ad'
os.makedirs(os.path.dirname(cache_path), exist_ok=True)

if os.path.exists(cache_path):
    print("Loading cached subset...")
    adata = ad.read_h5ad(cache_path)
else:
    print("Extracting from full atlas...")
    full = ad.read_h5ad('/ssd/data/agent/bio/atlas_merged_immune.h5ad')
    adata = full[full.obs['Dataset'].isin(SELECTED)].copy()
    adata.write(cache_path)
    del full
    gc.collect()

print(f"Subset: {adata.shape}")
print(f"Celltypes: {adata.obs['Celltype'].value_counts().to_dict()}")
print(f"Datasets: {adata.obs['Dataset'].value_counts().to_dict()}")

# ============================================================
# 2. Preprocess
# ============================================================
print("\n" + "=" * 60)
print("Step 2: Preprocess")
print("=" * 60)

sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')

X_expr = adata[:, adata.var['highly_variable']].X.copy()
if hasattr(X_expr, 'toarray'):
    X_expr = X_expr.toarray()
X_expr = np.array(X_expr, dtype=np.float32)

adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)
X_pca = adata_hvg.obsm['X_pca'].copy()
adata.obsm['X_pca'] = X_pca

print(f"Cells: {adata.n_obs}, HVG: {X_expr.shape[1]}")

# ============================================================
# 3. Pre-integration: high-res clustering + per-batch kNN
# ============================================================
print("\n" + "=" * 60)
print("Step 3: Pre-integration clustering and kNN")
print("=" * 60)

# High-res clustering (pre-integration)
sc.pp.neighbors(adata_hvg, n_pcs=50)
sc.tl.leiden(adata_hvg, resolution=2.0, key_added='leiden_pre_hr')
n_clusters_pre = adata_hvg.obs['leiden_pre_hr'].nunique()
print(f"Pre-integration high-res clusters: {n_clusters_pre}")

# Record small clusters (<200 cells)
cluster_sizes = adata_hvg.obs['leiden_pre_hr'].value_counts()
small_clusters = cluster_sizes[cluster_sizes < 200]
print(f"Small clusters (<200 cells): {len(small_clusters)}")
print(f"  Sizes: {small_clusters.values[:20]}")

# Pre-integration global kNN (for R(S))
print("Computing pre-integration kNN (FAISS GPU)...")
pre_indices = faiss_knn_gpu(X_pca, K)
print("Done.")

# ============================================================
# 4. Harmony integration
# ============================================================
print("\n" + "=" * 60)
print("Step 4: Harmony integration")
print("=" * 60)

import harmonypy as hm
t0 = time.time()
ho = hm.run_harmony(X_pca, adata.obs, 'Dataset')
Z = np.array(ho.Z_corr)
if Z.shape[0] != adata.n_obs:
    Z = Z.T
adata.obsm['X_harmony'] = Z
print(f"Harmony done in {time.time()-t0:.0f}s")

# Post-integration clustering
adata_hvg.obsm['X_harmony'] = Z
sc.pp.neighbors(adata_hvg, use_rep='X_harmony')
sc.tl.leiden(adata_hvg, resolution=2.0, key_added='leiden_post_hr')
n_clusters_post = adata_hvg.obs['leiden_post_hr'].nunique()
print(f"Post-integration high-res clusters: {n_clusters_post}")

# Post-integration kNN
print("Computing post-integration kNN (FAISS GPU)...")
post_indices = faiss_knn_gpu(Z, K)

# ============================================================
# 5. Compute R(S) dispersion score
# ============================================================
print("\n" + "=" * 60)
print("Step 5: R(S) dispersion score")
print("=" * 60)

t0 = time.time()
RS, DS = compute_RS(X_pca, Z, pre_indices)
print(f"R(S) computed in {time.time()-t0:.0f}s")
print(f"R(S): mean={RS.mean():.3f}, median={np.median(RS):.3f}, >3: {(RS>3).sum()}")

# R(S) by celltype
rs_df = pd.DataFrame({'RS': RS, 'celltype': adata.obs['Celltype'].values})
print("\nR(S) by celltype:")
print(rs_df.groupby('celltype')['RS'].agg(['mean', 'median', 'count']).sort_values('mean', ascending=False).round(3))

# ============================================================
# 6. Compute CNEM (global)
# ============================================================
print("\n" + "=" * 60)
print("Step 6: Global CNEM")
print("=" * 60)

t0 = time.time()
post_cnem = compute_cnem_chunked(X_expr, post_indices)
print(f"CNEM computed in {time.time()-t0:.0f}s")

cnem_df = pd.DataFrame({'CNEM': post_cnem, 'celltype': adata.obs['Celltype'].values})
print("\nCNEM by celltype:")
print(cnem_df.groupby('celltype')['CNEM'].agg(['mean', 'median', 'count']).sort_values('mean', ascending=False).round(4))

# ============================================================
# 7. Small cluster survival analysis
# ============================================================
print("\n" + "=" * 60)
print("Step 7: Small cluster survival analysis")
print("=" * 60)

pre_labels = adata_hvg.obs['leiden_pre_hr'].values
post_labels = adata_hvg.obs['leiden_post_hr'].values

survival_results = []
for cl_id, cl_size in small_clusters.items():
    cl_mask = pre_labels == cl_id
    cl_idx = np.where(cl_mask)[0]

    # Where do these cells end up after integration?
    post_assignments = post_labels[cl_idx]
    post_counts = pd.Series(post_assignments).value_counts()

    # Survival: what fraction ends up in the same post-integration cluster?
    max_post_frac = post_counts.iloc[0] / len(cl_idx)
    n_post_clusters = len(post_counts)

    # Purity of the best post-cluster (from this pre-cluster's perspective)
    best_post_cl = post_counts.index[0]
    best_post_total = (post_labels == best_post_cl).sum()
    purity = post_counts.iloc[0] / best_post_total

    # Mean R(S) and CNEM for this cluster
    mean_rs = RS[cl_idx].mean()
    mean_cnem = post_cnem[cl_idx].mean()

    # Dominant celltype
    celltypes_in_cl = adata.obs['Celltype'].values[cl_idx]
    dominant_ct = pd.Series(celltypes_in_cl).value_counts().index[0]
    ct_frac = pd.Series(celltypes_in_cl).value_counts().iloc[0] / len(cl_idx)

    survival_results.append({
        'pre_cluster': cl_id,
        'size': cl_size,
        'survival_frac': max_post_frac,
        'n_post_clusters': n_post_clusters,
        'purity_in_best_post': purity,
        'mean_RS': mean_rs,
        'mean_CNEM': mean_cnem,
        'dominant_celltype': dominant_ct,
        'ct_purity': ct_frac,
    })

survival_df = pd.DataFrame(survival_results).sort_values('survival_frac')
print("\n=== Small clusters most disrupted by integration ===")
print(survival_df.head(20).to_string(index=False))

# Destroyed clusters: survival < 50%
destroyed = survival_df[survival_df['survival_frac'] < 0.5]
survived = survival_df[survival_df['survival_frac'] >= 0.5]
print(f"\nDestroyed (<50% survival): {len(destroyed)} clusters")
print(f"Survived (>=50% survival): {len(survived)} clusters")

if len(destroyed) > 0 and len(survived) > 0:
    # Can R(S) or CNEM distinguish destroyed vs survived?
    dest_rs = destroyed['mean_RS'].values
    surv_rs = survived['mean_RS'].values
    stat, pval = stats.mannwhitneyu(dest_rs, surv_rs, alternative='greater')
    print(f"\nR(S) destroyed vs survived: Wilcoxon p={pval:.2e}")
    print(f"  Destroyed mean R(S)={dest_rs.mean():.3f}, Survived={surv_rs.mean():.3f}")

    dest_cnem = destroyed['mean_CNEM'].values
    surv_cnem = survived['mean_CNEM'].values
    stat2, pval2 = stats.mannwhitneyu(dest_cnem, surv_cnem, alternative='greater')
    print(f"CNEM destroyed vs survived: Wilcoxon p={pval2:.2e}")
    print(f"  Destroyed mean CNEM={dest_cnem.mean():.4f}, Survived={surv_cnem.mean():.4f}")

# ============================================================
# 8. pDC-specific analysis
# ============================================================
print("\n" + "=" * 60)
print("Step 8: pDC-specific analysis")
print("=" * 60)

pdc_mask = adata.obs['Celltype'] == 'pDC'
pdc_idx = np.where(pdc_mask)[0]
print(f"Total pDC: {len(pdc_idx)}")

# pDC pre-integration sub-clustering
pdc_pre_clusters = pre_labels[pdc_idx]
pdc_cluster_sizes = pd.Series(pdc_pre_clusters).value_counts()
print(f"pDC pre-integration clusters: {len(pdc_cluster_sizes)}")
print(pdc_cluster_sizes.head(10))

# Check survival of each pDC sub-cluster
print("\npDC sub-cluster survival:")
for cl_id, cl_size in pdc_cluster_sizes.items():
    if cl_size < 5:
        continue
    cl_mask_local = pdc_pre_clusters == cl_id
    cl_global_idx = pdc_idx[cl_mask_local]
    post_assign = post_labels[cl_global_idx]
    post_counts = pd.Series(post_assign).value_counts()
    surv = post_counts.iloc[0] / len(cl_global_idx)
    mean_rs_cl = RS[cl_global_idx].mean()
    mean_cnem_cl = post_cnem[cl_global_idx].mean()
    print(f"  Cluster {cl_id} (n={cl_size}): survival={surv:.1%}, R(S)={mean_rs_cl:.2f}, CNEM={mean_cnem_cl:.4f}, scattered to {len(post_counts)} post-clusters")

# ============================================================
# 9. Save results
# ============================================================
print("\n" + "=" * 60)
print("Step 9: Save results")
print("=" * 60)

results_df = pd.DataFrame({
    'RS': RS,
    'CNEM': post_cnem,
    'celltype': adata.obs['Celltype'].values,
    'dataset': adata.obs['Dataset'].values,
    'leiden_pre': pre_labels,
    'leiden_post': post_labels,
})
results_df.to_csv(f'{OUT_DIR}/phase1b_scores.csv', index=False)
survival_df.to_csv(f'{OUT_DIR}/phase1b_survival.csv', index=False)

print(f"\nTotal time: {time.time()-t_total:.0f}s")
print("Phase 1b complete!")
