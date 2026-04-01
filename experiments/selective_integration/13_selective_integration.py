"""Selective Integration: only align shared regions, preserve unique regions.

Algorithm:
1. For each cell, compute cross-batch kNN diversity → shared vs unique classification
2. Only integrate "shared" cells with Harmony
3. Unique cells keep original embedding
4. Unified embedding = corrected shared + original unique

This is a new paradigm: NOT "integrate everything then fix", but "classify first, integrate selectively".
"""
import sys; sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)
import scanpy as sc, anndata as ad, numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import adjusted_rand_score
import faiss, time, os, gc, warnings
warnings.filterwarnings('ignore')

OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/selective_integration/results'
SHARED_DIR = '/ssd/data/agent/bio/shared'
os.makedirs(OUT_DIR, exist_ok=True); os.makedirs(SHARED_DIR, exist_ok=True)

K = 30; GPU_ID = 0

def faiss_knn_gpu(X, k, gpu_id=GPU_ID):
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    idx = faiss.index_cpu_to_gpu(res, gpu_id, faiss.IndexFlatL2(X.shape[1]))
    idx.add(X); _, indices = idx.search(X, k + 1)
    return indices[:, 1:]

def compute_np(pre_idx, post_idx, k=K):
    n = pre_idx.shape[0]
    return np.array([len(set(pre_idx[i]) & set(post_idx[i])) / k for i in range(n)], dtype=np.float32)

t_total = time.time()

# 1. Load
print("=" * 60)
print("Step 1: Load")
adata = ad.read_h5ad('/ssd/data/agent/bio/eacn_example_001/experiments/np_guard/data/subset.h5ad')
print(f"Loaded: {adata.shape}")

sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)
X_pca = adata_hvg.obsm['X_pca'].astype(np.float32)
batches = adata.obs['Dataset'].values

print(f"Cells: {adata.n_obs}, Batches: {len(np.unique(batches))}")

# ============================================================
# 2. Step 1 — Region Classification: shared vs unique
# ============================================================
print("\n" + "=" * 60)
print("Step 2: Classify shared vs unique regions")
t0 = time.time()

# For each cell, compute cross-batch kNN diversity
# Use global kNN in PCA space
global_knn = faiss_knn_gpu(X_pca, K)

# For each cell: what fraction of its k neighbors come from OTHER batches?
batch_diversity = np.zeros(adata.n_obs, dtype=np.float32)
n_batches_in_knn = np.zeros(adata.n_obs, dtype=np.int32)

for i in range(adata.n_obs):
    nbr_batches = batches[global_knn[i]]
    my_batch = batches[i]
    other_batch_frac = (nbr_batches != my_batch).mean()
    batch_diversity[i] = other_batch_frac
    n_batches_in_knn[i] = len(set(nbr_batches))

# Classification: shared = high cross-batch diversity, unique = low
# Threshold: cells with <20% cross-batch neighbors are "unique"
UNIQUE_THRESHOLD = 0.2
is_unique = batch_diversity < UNIQUE_THRESHOLD
is_shared = ~is_unique

n_unique = is_unique.sum()
n_shared = is_shared.sum()
print(f"Shared: {n_shared} ({n_shared/adata.n_obs*100:.1f}%)")
print(f"Unique: {n_unique} ({n_unique/adata.n_obs*100:.1f}%)")
print(f"Classification time: {time.time()-t0:.0f}s")

# Which celltypes are unique?
print("\nUnique cells by celltype:")
ct = adata.obs['Celltype'].values
for c in pd.Series(ct).unique():
    c_mask = ct == c
    n_uniq = (is_unique & c_mask).sum()
    print(f"  {c}: {n_uniq}/{c_mask.sum()} unique ({n_uniq/c_mask.sum()*100:.0f}%)")

# ============================================================
# 3. Step 2 — Selective Integration: only integrate shared cells
# ============================================================
print("\n" + "=" * 60)
print("Step 3: Selective integration (Harmony on shared only)")
t0 = time.time()

import harmonypy as hm

# Harmony on shared cells only
shared_idx = np.where(is_shared)[0]
X_pca_shared = X_pca[shared_idx]
obs_shared = adata.obs.iloc[shared_idx]

ho = hm.run_harmony(X_pca_shared, obs_shared, 'Dataset')
Z_shared = np.array(ho.Z_corr, dtype=np.float32)
if Z_shared.shape[0] != len(shared_idx):
    Z_shared = Z_shared.T

print(f"Harmony on {len(shared_idx)} shared cells: {time.time()-t0:.0f}s")

# ============================================================
# 4. Step 3 — Unified Embedding
# ============================================================
print("\n" + "=" * 60)
print("Step 4: Unified embedding")

# Shared cells: use Harmony-corrected embedding
# Unique cells: use original PCA (projected to same space)
Z_selective = np.zeros((adata.n_obs, 50), dtype=np.float32)
Z_selective[shared_idx] = Z_shared
Z_selective[np.where(is_unique)[0]] = X_pca[np.where(is_unique)[0]]

# Also compute standard Harmony (all cells) for comparison
print("Standard Harmony (all cells)...")
t0 = time.time()
ho_all = hm.run_harmony(X_pca, adata.obs, 'Dataset')
Z_standard = np.array(ho_all.Z_corr, dtype=np.float32)
if Z_standard.shape[0] != adata.n_obs:
    Z_standard = Z_standard.T
print(f"Standard Harmony: {time.time()-t0:.0f}s")

# ============================================================
# 5. Compare: Standard vs Selective
# ============================================================
print("\n" + "=" * 60)
print("Step 5: Compare standard vs selective integration")

# Pre-integration reference kNN
pre_indices = faiss_knn_gpu(X_pca, K)

methods = {'Standard Harmony': Z_standard, 'Selective Integration': Z_selective}
results_all = {}

for name, Z in methods.items():
    print(f"\n--- {name} ---")
    post_idx = faiss_knn_gpu(Z, K)
    np_scores = compute_np(pre_indices, post_idx)
    print(f"  Global NP: {np_scores.mean():.4f}")

    # NP by celltype
    for c in pd.Series(ct).unique():
        c_mask = ct == c
        print(f"    {c}: NP={np_scores[c_mask].mean():.4f}")

    # Per-celltype subclustering + survival
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
        surv.append({'cluster': cl, 'size': len(ci), 'survival': pc.iloc[0]/len(ci), 'mean_NP': np_scores[ci].mean()})
    sdf = pd.DataFrame(surv)
    sdf['dispersed'] = (sdf['survival'] < 0.5).astype(int)

    n_disp = sdf['dispersed'].sum()
    mean_surv = sdf['survival'].mean()
    print(f"  Subclusters: {len(sdf)}, dispersed: {n_disp}, mean survival: {mean_surv:.3f}")
    results_all[name] = sdf

# Direct comparison
std_surv = results_all['Standard Harmony']
sel_surv = results_all['Selective Integration']
merged = std_surv[['cluster','size','survival']].merge(
    sel_surv[['cluster','survival']], on='cluster', suffixes=('_std','_sel'))
merged['improvement'] = merged['survival_sel'] - merged['survival_std']

disp = merged[merged['survival_std'] < 0.5]
if len(disp) > 0:
    n_imp = (disp['improvement'] > 0).sum()
    print(f"\n=== Dispersed clusters: Standard → Selective ===")
    print(f"Improved: {n_imp}/{len(disp)} ({n_imp/len(disp)*100:.0f}%)")
    print(f"Mean improvement: {disp['improvement'].mean():+.3f}")
    print(disp.nlargest(10, 'improvement')[['cluster','size','survival_std','survival_sel','improvement']].to_string(index=False))

# Batch mixing quality (for shared cells only)
# Check: does selective integration still mix batches well for shared cells?
from sklearn.metrics import silhouette_score
np.random.seed(42)
sample = np.random.choice(shared_idx, size=min(5000, len(shared_idx)), replace=False)
asw_std = silhouette_score(Z_standard[sample], batches[sample])
asw_sel = silhouette_score(Z_selective[sample], batches[sample])
print(f"\nBatch mixing (shared cells, ASW): Standard={asw_std:.3f}, Selective={asw_sel:.3f}")

# Save
merged.to_csv(f'{OUT_DIR}/selective_vs_standard.csv', index=False)
merged.to_csv(f'{SHARED_DIR}/selective_vs_standard.csv', index=False)

print(f"\nTotal: {time.time()-t_total:.0f}s")
print("Selective Integration experiment complete!")
