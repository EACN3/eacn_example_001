"""BD-Selective Integration: weighted blending based on Batch Diversity.

z_final(i) = BD(i) * z_harmony(i) + (1-BD(i)) * z_pre(i)
BD(i) = |unique batches in kNN(i)| / min(k, B)

Smooth version of selective integration - no hard threshold.
"""
import sys; sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)
import scanpy as sc, anndata as ad, numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import silhouette_score
import faiss, time, os, warnings
warnings.filterwarnings('ignore')

OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/bd_integration/results'
SHARED = '/ssd/data/agent/bio/shared'
os.makedirs(OUT_DIR, exist_ok=True); os.makedirs(SHARED, exist_ok=True)
K = 30; GPU_ID = 0

def faiss_knn_gpu(X, k, gpu_id=GPU_ID):
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    idx = faiss.index_cpu_to_gpu(res, gpu_id, faiss.IndexFlatL2(X.shape[1]))
    idx.add(X); _, indices = idx.search(X, k + 1)
    return indices[:, 1:]

def compute_np(pre_idx, post_idx, k=K):
    return np.array([len(set(pre_idx[i]) & set(post_idx[i])) / k for i in range(pre_idx.shape[0])], dtype=np.float32)

t_total = time.time()

# Load
print("Load...")
adata = ad.read_h5ad('/ssd/data/agent/bio/eacn_example_001/experiments/np_guard/data/subset.h5ad')
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)
X_pca = adata_hvg.obsm['X_pca'].astype(np.float32)
batches = adata.obs['Dataset'].values
B = len(np.unique(batches))
ct = adata.obs['Celltype'].values

# 1. Compute BD
print("Step 1: BD scores...")
t0 = time.time()
knn = faiss_knn_gpu(X_pca, K)
BD = np.zeros(adata.n_obs, dtype=np.float32)
for i in range(adata.n_obs):
    BD[i] = len(set(batches[knn[i]])) / min(K, B)
print(f"BD: mean={BD.mean():.3f}, time={time.time()-t0:.0f}s")
print(f"BD by celltype:")
for c in pd.Series(ct).unique():
    print(f"  {c}: BD={BD[ct==c].mean():.3f}")

# 2. Standard Harmony
print("\nStep 2: Harmony...")
t0 = time.time()
import harmonypy as hm
ho = hm.run_harmony(X_pca, adata.obs, 'Dataset')
Z_harm = np.array(ho.Z_corr, dtype=np.float32)
if Z_harm.shape[0] != adata.n_obs: Z_harm = Z_harm.T
t_harmony = time.time() - t0
print(f"Harmony: {t_harmony:.0f}s")

# 3. BD-weighted blending
print("\nStep 3: BD blending...")
Z_bd = np.zeros_like(Z_harm)
for i in range(adata.n_obs):
    Z_bd[i] = BD[i] * Z_harm[i] + (1 - BD[i]) * X_pca[i]
t_bd_total = time.time() - t_total
print(f"BD-Harmony total: {t_bd_total:.0f}s (Harmony overhead: ~0s)")

# 4. Compare
print("\nStep 4: Compare")
pre_indices = faiss_knn_gpu(X_pca, K)

for name, Z in [('Standard', Z_harm), ('BD-Selective', Z_bd)]:
    print(f"\n--- {name} ---")
    post_idx = faiss_knn_gpu(Z, K)
    np_scores = compute_np(pre_indices, post_idx)
    print(f"  NP: {np_scores.mean():.4f}")

    # Subclustering survival
    pre_lab = np.full(adata.n_obs, '', dtype=object)
    post_lab = np.full(adata.n_obs, '', dtype=object)
    for c in pd.Series(ct).unique():
        ci = np.where(ct == c)[0]
        if len(ci) < 50: continue
        a = ad.AnnData(obs=adata.obs.iloc[ci].copy())
        a.obsm['p'] = X_pca[ci]; a.obsm['z'] = Z[ci]
        sc.pp.neighbors(a, use_rep='p', n_neighbors=min(15, len(ci)-1))
        sc.tl.leiden(a, resolution=1.0, key_added='pre')
        sc.pp.neighbors(a, use_rep='z', n_neighbors=min(15, len(ci)-1))
        sc.tl.leiden(a, resolution=1.0, key_added='post')
        for j, gi in enumerate(ci):
            pre_lab[gi] = f"{c}_{a.obs['pre'].iloc[j]}"
            post_lab[gi] = f"{c}_{a.obs['post'].iloc[j]}"

    surv = []
    for cl in pd.Series(pre_lab).unique():
        ci = np.where(pre_lab == cl)[0]
        if len(ci) < 10: continue
        pc = pd.Series(post_lab[ci]).value_counts()
        surv.append({'cluster': cl, 'size': len(ci), 'survival': pc.iloc[0]/len(ci)})
    sdf = pd.DataFrame(surv)
    sdf['dispersed'] = (sdf['survival'] < 0.5).astype(int)
    print(f"  Dispersed: {sdf['dispersed'].sum()}/{len(sdf)}, mean survival: {sdf['survival'].mean():.3f}")

    if name == 'Standard': surv_std = sdf
    else: surv_bd = sdf

# Merge
merged = surv_std[['cluster','size','survival']].merge(
    surv_bd[['cluster','survival']], on='cluster', suffixes=('_std','_bd'))
merged['improvement'] = merged['survival_bd'] - merged['survival_std']
disp = merged[merged['survival_std'] < 0.5]

print(f"\n=== BD-Selective vs Standard ===")
if len(disp) > 0:
    n_imp = (disp['improvement'] > 0).sum()
    print(f"Dispersed improved: {n_imp}/{len(disp)} ({n_imp/len(disp)*100:.0f}%)")
    print(f"Mean improvement: {disp['improvement'].mean():+.3f}")
    print(disp.nlargest(10, 'improvement')[['cluster','size','survival_std','survival_bd','improvement']].to_string(index=False))

# Batch mixing
np.random.seed(42)
sample = np.random.choice(adata.n_obs, min(5000, adata.n_obs), replace=False)
asw_std = silhouette_score(Z_harm[sample], batches[sample])
asw_bd = silhouette_score(Z_bd[sample], batches[sample])
print(f"\nBatch mixing ASW: Standard={asw_std:.3f}, BD={asw_bd:.3f}")
print(f"Time: Standard Harmony={t_harmony:.0f}s, BD total={t_bd_total:.0f}s (ratio: {t_bd_total/t_harmony:.1f}x)")

merged.to_csv(f'{OUT_DIR}/bd_vs_standard.csv', index=False)
merged.to_csv(f'{SHARED}/bd_vs_standard.csv', index=False)
print(f"\nTotal: {time.time()-t_total:.0f}s")
print("BD-Selective Integration complete!")
