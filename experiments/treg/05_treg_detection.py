"""Detect TIGIT+CCR8- Treg fate in batch integration.
GPU-accelerated. Works on T cell subset to manage memory.

Strategy:
1. Extract T cells, identify Treg (FOXP3+) subset
2. Within Treg, classify TIGIT+CCR8- vs others
3. Run Harmony integration on T cell subset
4. Compute R(S) and sub-cluster survival for TIGIT+CCR8- Treg
5. Check if this subgroup is destroyed/absorbed
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
import torch
import faiss
import warnings, os, time, gc
warnings.filterwarnings('ignore')

OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/treg/results'
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
    """R(S): neighborhood dispersion score."""
    n = X_post.shape[0]
    DS = np.zeros(n, dtype=np.float32)
    for i in range(n):
        nbr = pre_indices[i]
        nbr_post = X_post[nbr]
        mean_p = nbr_post.mean(axis=0)
        DS[i] = np.mean(np.sum((nbr_post - mean_p) ** 2, axis=1))
    D_ref = np.median(DS)
    RS = DS / max(D_ref, 1e-8)
    return RS

# ============================================================
# 1. Extract T cells from a manageable subset
# ============================================================
print("=" * 60)
print("Step 1: Extract T cell subset")
print("=" * 60)
t_total = time.time()

# Use 10 batches with diverse cancer types for good coverage
# Include batches known to have Treg diversity
SELECTED = ['D5', 'D38', 'D40', 'D47', 'D61', 'D83', 'D88', 'D92']

cache_path = '/ssd/data/agent/bio/eacn_example_001/experiments/treg/data/tcell_subset.h5ad'
os.makedirs(os.path.dirname(cache_path), exist_ok=True)

if os.path.exists(cache_path):
    print("Loading cached T cell subset...")
    adata = ad.read_h5ad(cache_path)
else:
    print("Extracting from full atlas...")
    full = ad.read_h5ad('/ssd/data/agent/bio/atlas_merged_immune.h5ad')
    # Get selected batches, T cells only
    mask = full.obs['Dataset'].isin(SELECTED) & (full.obs['Celltype'] == 'T cell')
    adata = full[mask].copy()
    adata.write(cache_path)
    del full
    gc.collect()

print(f"T cell subset: {adata.shape}")
print(f"Datasets: {adata.obs['Dataset'].value_counts().to_dict()}")

# ============================================================
# 2. Identify Treg and TIGIT+CCR8- subgroup
# ============================================================
print("\n" + "=" * 60)
print("Step 2: Identify Treg subgroups")
print("=" * 60)

# Normalize for gene expression analysis
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Extract marker expressions
markers = ['FOXP3', 'TIGIT', 'CCR8', 'CTLA4', 'IL2RA', 'CD4']
for g in markers:
    if g in adata.var_names:
        expr = adata[:, g].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray()
        adata.obs[f'{g}_expr'] = expr.flatten()
        pct_pos = (expr.flatten() > 0).mean() * 100
        print(f"  {g}: {pct_pos:.1f}% positive, mean={expr.mean():.3f}")

# Define Treg: FOXP3 > 0 (log-normalized)
treg_mask = adata.obs['FOXP3_expr'] > 0
print(f"\nTreg (FOXP3+): {treg_mask.sum()} cells ({treg_mask.mean()*100:.1f}%)")

# Within Treg: TIGIT+CCR8- vs others
treg_tigit_pos = treg_mask & (adata.obs['TIGIT_expr'] > 0)
treg_ccr8_pos = treg_mask & (adata.obs['CCR8_expr'] > 0)
treg_target = treg_mask & (adata.obs['TIGIT_expr'] > 0) & (adata.obs['CCR8_expr'] == 0)
treg_ccr8_main = treg_mask & (adata.obs['CCR8_expr'] > 0)

print(f"Treg TIGIT+: {treg_tigit_pos.sum()}")
print(f"Treg CCR8+: {treg_ccr8_pos.sum()}")
print(f"Treg TIGIT+CCR8-: {treg_target.sum()} ({treg_target.mean()*100:.2f}%)")
print(f"Treg CCR8+ (mainstream): {treg_ccr8_main.sum()}")

# Label subgroups
adata.obs['treg_subgroup'] = 'non-Treg'
adata.obs.loc[treg_mask, 'treg_subgroup'] = 'Treg_other'
adata.obs.loc[treg_target, 'treg_subgroup'] = 'TIGIT+CCR8-'
adata.obs.loc[treg_ccr8_main, 'treg_subgroup'] = 'CCR8+'

print(f"\nSubgroup distribution:")
print(adata.obs['treg_subgroup'].value_counts())

# Cancer type distribution of TIGIT+CCR8- Treg
if treg_target.sum() > 0:
    print(f"\nTIGIT+CCR8- Treg by dataset:")
    print(adata.obs.loc[treg_target, 'Dataset'].value_counts())

# ============================================================
# 3. Preprocess and integrate
# ============================================================
print("\n" + "=" * 60)
print("Step 3: Preprocess and Harmony integration")
print("=" * 60)

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)
X_pca = adata_hvg.obsm['X_pca'].copy()

# Pre-integration kNN and clustering
print("Pre-integration kNN (FAISS GPU)...")
pre_indices = faiss_knn_gpu(X_pca, K)

sc.pp.neighbors(adata_hvg, n_pcs=50)
sc.tl.leiden(adata_hvg, resolution=2.0, key_added='leiden_pre')
print(f"Pre-integration clusters: {adata_hvg.obs['leiden_pre'].nunique()}")

# Harmony
print("Running Harmony...")
import harmonypy as hm
t0 = time.time()
ho = hm.run_harmony(X_pca, adata.obs, 'Dataset')
Z = np.array(ho.Z_corr)
if Z.shape[0] != adata.n_obs:
    Z = Z.T
print(f"Harmony done in {time.time()-t0:.0f}s")

# Post-integration clustering
adata_hvg.obsm['X_harmony'] = Z
sc.pp.neighbors(adata_hvg, use_rep='X_harmony')
sc.tl.leiden(adata_hvg, resolution=2.0, key_added='leiden_post')
print(f"Post-integration clusters: {adata_hvg.obs['leiden_post'].nunique()}")

# ============================================================
# 4. Compute R(S) for all cells
# ============================================================
print("\n" + "=" * 60)
print("Step 4: R(S) dispersion score")
print("=" * 60)

t0 = time.time()
RS = compute_RS(X_pca, Z, pre_indices)
print(f"R(S) computed in {time.time()-t0:.0f}s")

# R(S) by subgroup
rs_df = pd.DataFrame({
    'RS': RS,
    'subgroup': adata.obs['treg_subgroup'].values,
    'dataset': adata.obs['Dataset'].values,
})
print("\nR(S) by Treg subgroup:")
print(rs_df.groupby('subgroup')['RS'].agg(['mean', 'median', 'std', 'count']).round(3))

# Statistical tests
if treg_target.sum() > 0:
    target_rs = RS[treg_target.values]
    nontreg_rs = RS[(adata.obs['treg_subgroup'] == 'non-Treg').values]
    treg_other_rs = RS[(adata.obs['treg_subgroup'] == 'Treg_other').values]
    ccr8_rs = RS[(adata.obs['treg_subgroup'] == 'CCR8+').values]

    stat1, p1 = stats.mannwhitneyu(target_rs, nontreg_rs, alternative='greater')
    print(f"\nTIGIT+CCR8- vs non-Treg: p={p1:.2e}, mean RS {target_rs.mean():.3f} vs {nontreg_rs.mean():.3f}")

    if len(ccr8_rs) > 0:
        stat2, p2 = stats.mannwhitneyu(target_rs, ccr8_rs, alternative='greater')
        print(f"TIGIT+CCR8- vs CCR8+ Treg: p={p2:.2e}, mean RS {target_rs.mean():.3f} vs {ccr8_rs.mean():.3f}")

    # Fraction with R(S) > 3
    frac_dispersed = (target_rs > 3).mean()
    print(f"\nTIGIT+CCR8- with R(S)>3: {(target_rs>3).sum()}/{len(target_rs)} ({frac_dispersed:.1%})")

# ============================================================
# 5. Sub-cluster survival analysis for TIGIT+CCR8- Treg
# ============================================================
print("\n" + "=" * 60)
print("Step 5: TIGIT+CCR8- Treg cluster survival")
print("=" * 60)

pre_labels = adata_hvg.obs['leiden_pre'].values
post_labels = adata_hvg.obs['leiden_post'].values

# Which pre-clusters contain TIGIT+CCR8- Treg?
target_idx = np.where(treg_target.values)[0]
target_pre_clusters = pd.Series(pre_labels[target_idx]).value_counts()
print(f"TIGIT+CCR8- Treg found in {len(target_pre_clusters)} pre-clusters")
print(target_pre_clusters.head(10))

# For each cluster containing target cells, check survival
print("\nCluster-level survival:")
for cl_id, n_target in target_pre_clusters.items():
    cl_mask = pre_labels == cl_id
    cl_idx = np.where(cl_mask)[0]
    cl_size = len(cl_idx)

    post_assign = post_labels[cl_idx]
    post_counts = pd.Series(post_assign).value_counts()
    survival = post_counts.iloc[0] / cl_size
    n_scattered = len(post_counts)

    # Mean R(S)
    mean_rs = RS[cl_idx].mean()

    # What fraction of this cluster is TIGIT+CCR8-?
    target_frac = n_target / cl_size

    print(f"  Cluster {cl_id}: size={cl_size}, target={n_target} ({target_frac:.0%}), "
          f"survival={survival:.1%}, scattered_to={n_scattered}, R(S)={mean_rs:.2f}")

# Overall: are TIGIT+CCR8- cells scattered after integration?
target_post = post_labels[target_idx]
target_post_counts = pd.Series(target_post).value_counts()
max_frac = target_post_counts.iloc[0] / len(target_idx)
n_post_clusters = len(target_post_counts)
print(f"\nTIGIT+CCR8- Treg after integration:")
print(f"  Scattered to {n_post_clusters} clusters")
print(f"  Largest cluster contains {target_post_counts.iloc[0]}/{len(target_idx)} ({max_frac:.1%})")
print(f"  Top clusters: {target_post_counts.head(5).to_dict()}")

# ARI between pre and post for Treg cells
from sklearn.metrics import adjusted_rand_score
treg_idx = np.where(treg_mask.values)[0]
if len(treg_idx) > 0:
    ari = adjusted_rand_score(pre_labels[treg_idx], post_labels[treg_idx])
    print(f"\nTreg sub-cluster ARI (pre vs post): {ari:.3f}")

# ============================================================
# 6. Save
# ============================================================
print("\n" + "=" * 60)
print("Step 6: Save results")
print("=" * 60)

results = pd.DataFrame({
    'RS': RS,
    'subgroup': adata.obs['treg_subgroup'].values,
    'dataset': adata.obs['Dataset'].values,
    'leiden_pre': pre_labels,
    'leiden_post': post_labels,
})
results.to_csv(f'{OUT_DIR}/treg_rs_scores.csv', index=False)

print(f"\nTotal time: {time.time()-t_total:.0f}s")
print("Treg detection experiment complete!")
