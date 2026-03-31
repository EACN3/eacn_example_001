"""Artificial ablation experiment: Test if integration WOULD destroy pDC
if they were fewer and less distinct.

Strategy:
1. Subsample pDC to ~50-100 cells (simulating extreme rarity)
2. Merge them into a single batch (simulating batch-specific rare subgroup)
3. Run integration and check if they survive as independent cluster
4. Compute CNEM on the subsampled pDC to see if the metric detects disruption

Also: Check if pDC cluster 24 (93 cells, NHL-only) gets disrupted.
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
import torch
import faiss
import warnings, os, time
warnings.filterwarnings('ignore')

DATA_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/phase1/data'
OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/phase1/results'
os.makedirs(OUT_DIR, exist_ok=True)

K = 30
GPU_ID = 0

def faiss_knn_gpu(X, k, gpu_id=GPU_ID):
    X = np.ascontiguousarray(X, dtype=np.float32)
    d = X.shape[1]
    res = faiss.StandardGpuResources()
    index = faiss.index_cpu_to_gpu(res, gpu_id, faiss.IndexFlatL2(d))
    index.add(X)
    _, indices = index.search(X, k + 1)
    return indices[:, 1:]

def compute_cnem_gpu(X_expr_torch, indices_np):
    indices_t = torch.from_numpy(indices_np.astype(np.int64)).cuda(GPU_ID)
    nbr_expr = X_expr_torch[indices_t]
    mean_nbr = nbr_expr.mean(dim=1)
    cos_sim = torch.nn.functional.cosine_similarity(X_expr_torch, mean_nbr, dim=1)
    return (1.0 - cos_sim).cpu().numpy()

# ============================================================
# 1. Load and preprocess
# ============================================================
print("=" * 60)
print("Step 1: Load and preprocess")
print("=" * 60)
t_total = time.time()

adata = ad.read_h5ad(f'{DATA_DIR}/immune_subset_phase1.h5ad')

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
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

print(f"Cells: {adata.n_obs}")

# ============================================================
# 2. Artificial ablation: subsample pDC to extreme rarity
# ============================================================
print("\n" + "=" * 60)
print("Step 2: Artificial ablation experiments")
print("=" * 60)

import harmonypy as hm

pdc_idx = np.where(adata.obs['Celltype'] == 'pDC')[0]
non_pdc_idx = np.where(adata.obs['Celltype'] != 'pDC')[0]
print(f"Total pDC: {len(pdc_idx)}")

for n_keep in [50, 100, 200, 500]:
    print(f"\n--- Ablation: keep {n_keep} pDC ---")
    t0 = time.time()

    # Randomly subsample pDC
    np.random.seed(42)
    pdc_keep = np.random.choice(pdc_idx, size=min(n_keep, len(pdc_idx)), replace=False)

    # Combine with all non-pDC
    keep_idx = np.sort(np.concatenate([non_pdc_idx, pdc_keep]))

    # Subset data
    X_pca_sub = adata.obsm['X_pca'][keep_idx]
    X_expr_sub = X_expr[keep_idx]
    obs_sub = adata.obs.iloc[keep_idx].copy()

    # Run Harmony
    ho = hm.run_harmony(X_pca_sub, obs_sub, 'Dataset', verbose=False)
    Z = np.array(ho.Z_corr)
    if Z.shape[0] != len(keep_idx):
        Z = Z.T

    # Post-integration clustering
    adata_sub = ad.AnnData(obs=obs_sub)
    adata_sub.obsm['X_harmony'] = Z
    sc.pp.neighbors(adata_sub, use_rep='X_harmony')
    sc.tl.leiden(adata_sub, resolution=1.0)

    # Check pDC cluster survival
    pdc_mask_sub = obs_sub['Celltype'] == 'pDC'
    pdc_clusters = adata_sub.obs.loc[pdc_mask_sub, 'leiden'].value_counts()
    best_cl = pdc_clusters.index[0]
    best_count = pdc_clusters.iloc[0]
    total_in_best = (adata_sub.obs['leiden'] == best_cl).sum()
    purity = best_count / total_in_best
    recall = best_count / pdc_mask_sub.sum()

    # Compute CNEM on GPU
    X_expr_torch = torch.from_numpy(X_expr_sub).cuda(GPU_ID)
    indices = faiss_knn_gpu(Z, K)
    post_cnem = compute_cnem_gpu(X_expr_torch, indices)

    pdc_cnem = post_cnem[pdc_mask_sub.values]
    non_pdc_cnem = post_cnem[~pdc_mask_sub.values]

    stat, pval = stats.mannwhitneyu(pdc_cnem, non_pdc_cnem, alternative='greater')
    auc = roc_auc_score(pdc_mask_sub.astype(int).values, post_cnem)

    # Top 5% enrichment
    threshold = np.percentile(post_cnem, 95)
    top5_mask = post_cnem >= threshold
    pdc_in_top5 = top5_mask[pdc_mask_sub.values].sum()
    n_top5 = top5_mask.sum()
    rare_frac = pdc_mask_sub.sum() / len(pdc_mask_sub)
    enrichment = (pdc_in_top5 / n_top5) / rare_frac if n_top5 > 0 else 0

    print(f"  Cluster: purity={purity:.1%}, recall={recall:.1%}")
    print(f"  pDC scattered to {len(pdc_clusters)} clusters (top: {pdc_clusters.head(3).to_dict()})")
    print(f"  CNEM: pDC mean={pdc_cnem.mean():.4f}, others mean={non_pdc_cnem.mean():.4f}")
    print(f"  Wilcoxon pDC>others: p={pval:.2e}")
    print(f"  ROC-AUC: {auc:.3f}")
    print(f"  Top 5%: pDC {pdc_in_top5}/{pdc_mask_sub.sum()}, OR={enrichment:.2f}")
    print(f"  Time: {time.time()-t0:.1f}s")

# ============================================================
# 3. Extreme ablation: put all pDC in single batch only
# ============================================================
print("\n" + "=" * 60)
print("Step 3: Single-batch pDC (batch-specific rare subgroup)")
print("=" * 60)

# Keep only 50 pDC, all from one batch (D5)
pdc_d5 = np.where((adata.obs['Celltype'] == 'pDC') & (adata.obs['Dataset'] == 'D5'))[0]
print(f"pDC in D5: {len(pdc_d5)}")

if len(pdc_d5) > 50:
    np.random.seed(42)
    pdc_d5_keep = np.random.choice(pdc_d5, size=50, replace=False)
else:
    pdc_d5_keep = pdc_d5

# Remove ALL pDC except selected D5 ones
non_pdc_all = np.where(adata.obs['Celltype'] != 'pDC')[0]
keep_idx = np.sort(np.concatenate([non_pdc_all, pdc_d5_keep]))

X_pca_sub = adata.obsm['X_pca'][keep_idx]
X_expr_sub = X_expr[keep_idx]
obs_sub = adata.obs.iloc[keep_idx].copy()

print(f"Subset: {len(keep_idx)} cells, pDC={len(pdc_d5_keep)} (all in D5)")

# Run Harmony
t0 = time.time()
ho = hm.run_harmony(X_pca_sub, obs_sub, 'Dataset', verbose=False)
Z = np.array(ho.Z_corr)
if Z.shape[0] != len(keep_idx):
    Z = Z.T

# Clustering
adata_sub = ad.AnnData(obs=obs_sub)
adata_sub.obsm['X_harmony'] = Z
sc.pp.neighbors(adata_sub, use_rep='X_harmony')
sc.tl.leiden(adata_sub, resolution=1.0)

pdc_mask_sub = obs_sub['Celltype'] == 'pDC'
pdc_clusters = adata_sub.obs.loc[pdc_mask_sub, 'leiden'].value_counts()
best_cl = pdc_clusters.index[0]
best_count = pdc_clusters.iloc[0]
total_in_best = (adata_sub.obs['leiden'] == best_cl).sum()
purity = best_count / total_in_best
recall = best_count / pdc_mask_sub.sum()

# CNEM
X_expr_torch = torch.from_numpy(X_expr_sub).cuda(GPU_ID)
indices = faiss_knn_gpu(Z, K)
post_cnem = compute_cnem_gpu(X_expr_torch, indices)

pdc_cnem = post_cnem[pdc_mask_sub.values]
non_pdc_cnem = post_cnem[~pdc_mask_sub.values]
stat, pval = stats.mannwhitneyu(pdc_cnem, non_pdc_cnem, alternative='greater')
auc = roc_auc_score(pdc_mask_sub.astype(int).values, post_cnem)
threshold = np.percentile(post_cnem, 95)
top5_mask = post_cnem >= threshold
pdc_in_top5 = top5_mask[pdc_mask_sub.values].sum()
n_top5 = top5_mask.sum()
rare_frac = pdc_mask_sub.sum() / len(pdc_mask_sub)
enrichment = (pdc_in_top5 / n_top5) / rare_frac if n_top5 > 0 else 0

print(f"  Cluster: purity={purity:.1%}, recall={recall:.1%}")
print(f"  pDC in clusters: {pdc_clusters.head(5).to_dict()}")
print(f"  CNEM: pDC mean={pdc_cnem.mean():.4f}, others={non_pdc_cnem.mean():.4f}")
print(f"  Wilcoxon: p={pval:.2e}, ROC-AUC: {auc:.3f}")
print(f"  Top 5%: pDC {pdc_in_top5}/{pdc_mask_sub.sum()}, OR={enrichment:.2f}")
print(f"  Time: {time.time()-t0:.1f}s")

print(f"\nTotal time: {time.time()-t_total:.0f}s")
print("Ablation experiment complete!")
