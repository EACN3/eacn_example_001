"""NP-Guard: Adaptive protection via batch label shuffling for high-risk cells.
Based on ML agent's design. GPU-accelerated.

Two-phase approach (no scVI source modification needed):
Phase A: Standard scVI training -> get latent -> compute NP risk
Phase B: Shuffle batch labels for high-risk cells -> retrain -> compare survival
"""
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import adjusted_rand_score
import faiss
import torch
import warnings, os, time, gc
warnings.filterwarnings('ignore')

OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/np_guard/results'
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

def compute_np_fast(pre_indices, post_indices):
    n = pre_indices.shape[0]
    k = pre_indices.shape[1]
    np_scores = np.zeros(n, dtype=np.float32)
    for i in range(n):
        np_scores[i] = len(set(pre_indices[i]) & set(post_indices[i])) / k
    return np_scores

def compute_survival(pre_labels, post_labels):
    results = []
    for cl_id in pd.Series(pre_labels).unique():
        cl_mask = pre_labels == cl_id
        cl_idx = np.where(cl_mask)[0]
        if len(cl_idx) < 10:
            continue
        post_assign = post_labels[cl_idx]
        post_counts = pd.Series(post_assign).value_counts()
        survival = post_counts.iloc[0] / len(cl_idx)
        results.append({'cluster': cl_id, 'size': len(cl_idx), 'survival': survival})
    return pd.DataFrame(results)

# ============================================================
# 1. Load and preprocess
# ============================================================
print("=" * 60)
print("Step 1: Load and preprocess")
t_total = time.time()

CACHE = '/ssd/data/agent/bio/eacn_example_001/experiments/np_guard/data/subset.h5ad'
os.makedirs(os.path.dirname(CACHE), exist_ok=True)
if os.path.exists(CACHE):
    adata = ad.read_h5ad(CACHE)
else:
    SELECTED = ['D5', 'D38', 'D80', 'D83', 'D88', 'D94', 'D100', 'D102']
    full = ad.read_h5ad('/ssd/data/agent/bio/atlas_merged_immune.h5ad')
    adata = full[full.obs['Dataset'].isin(SELECTED)].copy()
    adata.write(CACHE)
    del full
    gc.collect()
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')

adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)
X_pca = adata_hvg.obsm['X_pca'].astype(np.float32)

print(f"Cells: {adata.n_obs}")

# Pre-integration kNN
pre_indices = faiss_knn_gpu(X_pca, K)

# Pre-integration subclustering (per celltype, for survival evaluation)
pre_labels = np.full(adata.n_obs, '', dtype=object)
for ct in adata.obs['Celltype'].unique():
    ct_mask = adata.obs['Celltype'] == ct
    ct_idx = np.where(ct_mask)[0]
    if len(ct_idx) < 50:
        continue
    adata_ct = ad.AnnData(obs=adata.obs.iloc[ct_idx].copy())
    adata_ct.obsm['X_pca'] = X_pca[ct_idx]
    sc.pp.neighbors(adata_ct, use_rep='X_pca', n_neighbors=min(15, len(ct_idx)-1))
    sc.tl.leiden(adata_ct, resolution=1.0)
    for j, gi in enumerate(ct_idx):
        pre_labels[gi] = f"{ct}_{adata_ct.obs['leiden'].iloc[j]}"

# ============================================================
# 2. Phase A: Standard scVI
# ============================================================
print("\nStep 2: Phase A - Standard scVI")
t0 = time.time()

try:
    from scvi.model import SCVI
    SCVI.setup_anndata(adata, batch_key='Dataset')
    model_std = SCVI(adata, n_latent=30, n_layers=2)
    model_std.train(max_epochs=100, accelerator='gpu', devices=[GPU_ID])
    Z_std = model_std.get_latent_representation().astype(np.float32)
    scvi_available = True
    print(f"Standard scVI done in {time.time()-t0:.0f}s")
except Exception as e:
    print(f"scVI failed: {e}")
    print("Falling back to Harmony for NP-Guard demonstration")
    import harmonypy as hm
    ho = hm.run_harmony(X_pca, adata.obs, 'Dataset')
    Z_std = np.array(ho.Z_corr, dtype=np.float32)
    if Z_std.shape[0] != adata.n_obs:
        Z_std = Z_std.T
    scvi_available = False
    print(f"Harmony done in {time.time()-t0:.0f}s")

# Post kNN and NP
post_indices_std = faiss_knn_gpu(Z_std, K)
np_std = compute_np_fast(pre_indices, post_indices_std)
print(f"Standard NP: mean={np_std.mean():.4f}")

# ============================================================
# 3. Compute risk scores
# ============================================================
print("\nStep 3: Compute risk scores")

# Risk = cells with low NP (their neighborhoods changed the most)
# Use percentile-based threshold
risk = 1.0 - np_std  # high risk = low NP
# Normalize to [0, 1]
risk = (risk - risk.min()) / (risk.max() - risk.min() + 1e-8)

high_risk_mask = risk > 0.5
print(f"High-risk cells (risk>0.5): {high_risk_mask.sum()} ({high_risk_mask.mean()*100:.1f}%)")
print(f"High-risk by celltype:")
for ct in adata.obs['Celltype'].unique():
    ct_mask = adata.obs['Celltype'] == ct
    hr = (high_risk_mask & ct_mask.values).sum()
    total = ct_mask.sum()
    print(f"  {ct}: {hr}/{total} ({hr/total*100:.1f}%)")

# ============================================================
# 4. Phase B: NP-Guard protected integration
# ============================================================
print("\nStep 4: Phase B - NP-Guard protected integration")
t0 = time.time()

# Shuffle batch labels for high-risk cells
adata_protected = adata.copy()
np.random.seed(42)
batches = adata.obs['Dataset'].unique()

n_shuffled = 0
for i in np.where(high_risk_mask)[0]:
    if np.random.random() < risk[i]:
        other_batches = [b for b in batches if b != adata.obs['Dataset'].iloc[i]]
        adata_protected.obs['Dataset'].iloc[i] = np.random.choice(other_batches)
        n_shuffled += 1

print(f"Shuffled {n_shuffled} batch labels ({n_shuffled/adata.n_obs*100:.1f}%)")

if scvi_available:
    SCVI.setup_anndata(adata_protected, batch_key='Dataset')
    model_prot = SCVI(adata_protected, n_latent=30, n_layers=2)
    model_prot.train(max_epochs=100, accelerator='gpu', devices=[GPU_ID])
    Z_prot = model_prot.get_latent_representation().astype(np.float32)
else:
    ho2 = hm.run_harmony(
        X_pca, adata_protected.obs, 'Dataset')
    Z_prot = np.array(ho2.Z_corr, dtype=np.float32)
    if Z_prot.shape[0] != adata.n_obs:
        Z_prot = Z_prot.T

print(f"Protected integration done in {time.time()-t0:.0f}s")

# Protected NP
post_indices_prot = faiss_knn_gpu(Z_prot, K)
np_prot = compute_np_fast(pre_indices, post_indices_prot)
print(f"Protected NP: mean={np_prot.mean():.4f} (was {np_std.mean():.4f})")

# ============================================================
# 5. Compare survival: standard vs protected
# ============================================================
print("\n" + "=" * 60)
print("Step 5: Compare survival")

# Post-integration clustering for both
for name, Z in [('standard', Z_std), ('protected', Z_prot)]:
    adata_tmp = ad.AnnData(obs=adata.obs.copy())
    adata_tmp.obsm['X_int'] = Z
    # Per-celltype post-clustering
    post_labels = np.full(adata.n_obs, '', dtype=object)
    for ct in adata.obs['Celltype'].unique():
        ct_mask = adata.obs['Celltype'] == ct
        ct_idx = np.where(ct_mask)[0]
        if len(ct_idx) < 50:
            continue
        adata_ct = ad.AnnData(obs=adata.obs.iloc[ct_idx].copy())
        adata_ct.obsm['X_int'] = Z[ct_idx]
        sc.pp.neighbors(adata_ct, use_rep='X_int', n_neighbors=min(15, len(ct_idx)-1))
        sc.tl.leiden(adata_ct, resolution=1.0)
        for j, gi in enumerate(ct_idx):
            post_labels[gi] = f"{ct}_{adata_ct.obs['leiden'].iloc[j]}"

    surv = compute_survival(pre_labels, post_labels)
    surv['method'] = name
    if name == 'standard':
        surv_std = surv
    else:
        surv_prot = surv

# Merge and compare
merged = surv_std[['cluster', 'size', 'survival']].merge(
    surv_prot[['cluster', 'survival']], on='cluster', suffixes=('_std', '_prot'))

merged['improvement'] = merged['survival_prot'] - merged['survival_std']
merged['dispersed_std'] = merged['survival_std'] < 0.5

print("\n=== Dispersed clusters: standard vs NP-Guard ===")
dispersed = merged[merged['dispersed_std']].sort_values('survival_std')
print(dispersed[['cluster', 'size', 'survival_std', 'survival_prot', 'improvement']].to_string(index=False))

if len(dispersed) > 0:
    mean_imp = dispersed['improvement'].mean()
    n_improved = (dispersed['improvement'] > 0).sum()
    print(f"\nDispersed clusters: {len(dispersed)}")
    print(f"Improved: {n_improved}/{len(dispersed)} ({n_improved/len(dispersed)*100:.0f}%)")
    print(f"Mean survival improvement: {mean_imp:+.3f}")

# Overall stats
print("\n=== Overall ===")
print(f"Standard: mean survival={merged['survival_std'].mean():.3f}")
print(f"Protected: mean survival={merged['survival_prot'].mean():.3f}")

# Save
merged.to_csv(f'{OUT_DIR}/np_guard_comparison.csv', index=False)

print(f"\nTotal time: {time.time()-t_total:.0f}s")
print("NP-Guard validation complete!")
