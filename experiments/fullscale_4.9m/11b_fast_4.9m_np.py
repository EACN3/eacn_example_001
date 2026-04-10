"""Fast 4.9M NP: avoid dense matrix, use sparse PCA + chunked processing.
Key optimizations:
1. Use TruncatedSVD (sparse) instead of sc.pp.pca (dense)
2. Skip scale (use log-normalized sparse directly)
3. All GPU-accelerated
"""
import sys; sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)
import scanpy as sc, anndata as ad, numpy as np, pandas as pd
import faiss, os, time, gc, glob
from scipy.sparse import issparse
from sklearn.decomposition import TruncatedSVD
import warnings; warnings.filterwarnings('ignore')

OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/fullscale_4.9m/results'
SHARED = '/ssd/data/agent/bio/shared'
os.makedirs(OUT_DIR, exist_ok=True); os.makedirs(SHARED, exist_ok=True)

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
        if s % 1000000 == 0 and s > 0:
            print(f"    kNN: {s:,}/{n:,}")
    return all_idx

t_total = time.time()

# 1. Load all
print("Step 1: Load")
files = sorted(glob.glob(os.path.join(DATA_DIR, '*.h5ad')))
adatas = [ad.read_h5ad(f) for f in files]
print(f"Loaded {len(files)} files")

adata = ad.concat(adatas, join='inner')
del adatas; gc.collect()
n = adata.n_obs
print(f"Total: {n:,} cells, {adata.shape[1]} genes")

# 2. Fast preprocess: normalize + log + HVG + sparse TruncatedSVD (NO scale/dense)
print("\nStep 2: Fast preprocess (sparse)")
t0 = time.time()
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')

X_hvg = adata[:, adata.var['highly_variable']].X
if not issparse(X_hvg):
    from scipy.sparse import csr_matrix
    X_hvg = csr_matrix(X_hvg)

# TruncatedSVD works on sparse directly — no need for scale/dense
svd = TruncatedSVD(n_components=50, random_state=42)
X_pca = svd.fit_transform(X_hvg).astype(np.float32)
print(f"Sparse SVD done: {X_pca.shape}, explained var: {svd.explained_variance_ratio_.sum():.3f}")

obs_full = adata.obs.copy()
del adata, X_hvg; gc.collect()
print(f"Preprocess: {time.time()-t0:.0f}s")

# 3. Pre-kNN
print("\nStep 3: Pre-kNN (FAISS GPU)")
t0 = time.time()
pre_indices = faiss_knn_gpu_chunked(X_pca, K)
print(f"Pre-kNN: {time.time()-t0:.0f}s")

# 4. Harmony
print("\nStep 4: Harmony")
t0 = time.time()
import harmonypy as hm
ho = hm.run_harmony(X_pca, obs_full, 'Dataset')
Z = np.array(ho.Z_corr, dtype=np.float32)
if Z.shape[0] != n: Z = Z.T
print(f"Harmony: {time.time()-t0:.0f}s")

# 5. Post-kNN
print("\nStep 5: Post-kNN (FAISS GPU)")
t0 = time.time()
post_indices = faiss_knn_gpu_chunked(Z, K)
print(f"Post-kNN: {time.time()-t0:.0f}s")

# 6. NP
print("\nStep 6: NP")
t0 = time.time()
np_scores = np.zeros(n, dtype=np.float32)
for i in range(n):
    np_scores[i] = len(set(pre_indices[i]) & set(post_indices[i])) / K
    if i % 1000000 == 0 and i > 0:
        print(f"    NP: {i:,}/{n:,}")
print(f"NP: {time.time()-t0:.0f}s, mean={np_scores.mean():.4f}")

# 7. Results
print("\nStep 7: Results")
ct = obs_full['Celltype'].values
np_df = pd.DataFrame({'NP': np_scores, 'celltype': ct, 'dataset': obs_full['Dataset'].values})
print("\nNP by celltype:")
print(np_df.groupby('celltype')['NP'].agg(['mean', 'median', 'count']).sort_values('mean').round(4))
print(f"\n<0.3: {(np_scores<0.3).sum():,} ({(np_scores<0.3).mean()*100:.1f}%)")
print(f"<0.1: {(np_scores<0.1).sum():,} ({(np_scores<0.1).mean()*100:.1f}%)")

np_df.to_csv(f'{OUT_DIR}/fullscale_4.9m_np.csv', index=False, float_format='%.4f')
np_df.to_csv(f'{SHARED}/fullscale_4.9m_np.csv', index=False, float_format='%.4f')

print(f"\nTotal: {time.time()-t_total:.0f}s ({(time.time()-t_total)/60:.0f} min)")
print(f"Cells: {n:,}")
print("DONE!")
