"""Full-scale 4.9M NP on complete pan-cancer atlas (102 datasets).
Strategy: Load all h5ad files, concat, HVG+PCA+Harmony+NP.
GPU-accelerated with FAISS.
"""
import sys
sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)

import scanpy as sc, anndata as ad, numpy as np, pandas as pd
import faiss, os, time, gc, glob
import warnings; warnings.filterwarnings('ignore')

OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/fullscale_4.9m/results'
os.makedirs(OUT_DIR, exist_ok=True)
SHARED = '/ssd/data/agent/bio/shared'
os.makedirs(SHARED, exist_ok=True)

K = 30; GPU_ID = 0
DATA_DIR = '/ssd/data/agent/bio/atlas_full/atlas_dataset'

def faiss_knn_gpu_chunked(X, k, gpu_id=GPU_ID, chunk_size=200000):
    X = np.ascontiguousarray(X, dtype=np.float32)
    n, d = X.shape
    res = faiss.StandardGpuResources()
    index = faiss.index_cpu_to_gpu(res, gpu_id, faiss.IndexFlatL2(d))
    index.add(X)
    all_idx = np.zeros((n, k), dtype=np.int64)
    for s in range(0, n, chunk_size):
        e = min(s + chunk_size, n)
        _, idx = index.search(X[s:e], k + 1)
        all_idx[s:e] = idx[:, 1:]
        if s % 500000 == 0 and s > 0:
            print(f"    kNN: {s}/{n}")
    return all_idx

# ============================================================
# 1. Load and concatenate all datasets
# ============================================================
print("=" * 60)
print("Step 1: Load all datasets")
t_total = time.time()

files = sorted(glob.glob(os.path.join(DATA_DIR, '*.h5ad')))
print(f"Found {len(files)} h5ad files")

adatas = []
total_cells = 0
for i, f in enumerate(files):
    a = ad.read_h5ad(f)
    total_cells += a.n_obs
    adatas.append(a)
    if (i + 1) % 20 == 0:
        print(f"  Loaded {i+1}/{len(files)} files, {total_cells:,} cells so far")

print(f"Loaded all {len(files)} files: {total_cells:,} cells total")
print(f"Load time: {time.time()-t_total:.0f}s")

# Concatenate
print("Concatenating...")
t0 = time.time()
adata = ad.concat(adatas, join='inner')
del adatas; gc.collect()
print(f"Concatenated: {adata.shape}, took {time.time()-t0:.0f}s")

# ============================================================
# 2. Preprocess
# ============================================================
print("\nStep 2: Preprocess")
t0 = time.time()

sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')
adata_hvg = adata[:, adata.var['highly_variable']].copy()

# Free adata to save memory (keep obs)
obs_full = adata.obs.copy()
del adata; gc.collect()

sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)
X_pca = adata_hvg.obsm['X_pca'].astype(np.float32)
print(f"Preprocessed: {X_pca.shape}, took {time.time()-t0:.0f}s")

del adata_hvg; gc.collect()

# ============================================================
# 3. Pre-integration kNN
# ============================================================
print("\nStep 3: Pre-integration kNN (FAISS GPU)")
t0 = time.time()
pre_indices = faiss_knn_gpu_chunked(X_pca, K)
print(f"Pre-kNN done in {time.time()-t0:.0f}s")

# ============================================================
# 4. Harmony
# ============================================================
print("\nStep 4: Harmony integration")
t0 = time.time()
import harmonypy as hm
ho = hm.run_harmony(X_pca, obs_full, 'Dataset')
Z = np.array(ho.Z_corr, dtype=np.float32)
if Z.shape[0] != len(obs_full):
    Z = Z.T
print(f"Harmony done in {time.time()-t0:.0f}s")

# ============================================================
# 5. Post-integration kNN
# ============================================================
print("\nStep 5: Post-integration kNN (FAISS GPU)")
t0 = time.time()
post_indices = faiss_knn_gpu_chunked(Z, K)
print(f"Post-kNN done in {time.time()-t0:.0f}s")

# ============================================================
# 6. Compute NP
# ============================================================
print("\nStep 6: Compute NP")
t0 = time.time()
n = len(obs_full)
np_scores = np.zeros(n, dtype=np.float32)
for i in range(n):
    np_scores[i] = len(set(pre_indices[i]) & set(post_indices[i])) / K
    if i % 500000 == 0 and i > 0:
        print(f"    NP: {i}/{n}")
print(f"NP done in {time.time()-t0:.0f}s")
print(f"NP: mean={np_scores.mean():.4f}, median={np.median(np_scores):.4f}")

# ============================================================
# 7. Results
# ============================================================
print("\nStep 7: Results")

celltypes = obs_full['Celltype'].values
np_df = pd.DataFrame({'NP': np_scores, 'celltype': celltypes, 'dataset': obs_full['Dataset'].values})

print("\nNP by celltype:")
ct_stats = np_df.groupby('celltype')['NP'].agg(['mean', 'median', 'count']).sort_values('mean')
print(ct_stats.round(4))

print(f"\nGlobal stats:")
print(f"  Total cells: {n:,}")
print(f"  Mean NP: {np_scores.mean():.4f}")
print(f"  <0.3: {(np_scores<0.3).sum():,} ({(np_scores<0.3).mean()*100:.1f}%)")
print(f"  <0.1: {(np_scores<0.1).sum():,} ({(np_scores<0.1).mean()*100:.1f}%)")

# Save
np_df.to_csv(f'{OUT_DIR}/fullscale_4.9m_np.csv', index=False, float_format='%.4f')
# Also to shared
np_df.to_csv(f'{SHARED}/fullscale_4.9m_np.csv', index=False, float_format='%.4f')

print(f"\nTotal time: {time.time()-t_total:.0f}s")
print(f"Cells: {n:,}")
print("4.9M full-scale NP complete!")
