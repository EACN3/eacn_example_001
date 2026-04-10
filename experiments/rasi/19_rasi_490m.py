"""RASI 490M: incremental integration. Target <20min.
Strategy:
1. Load datasets streaming (not concat all)
2. Select reference batch (largest), run Harmony on it
3. Map other batches onto reference via sc.tl.ingest
4. FAISS IVF kNN + numba NP
5. Cache everything
"""
import sys; sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)
import os; os.environ['OMP_NUM_THREADS']='32'; os.environ['MKL_NUM_THREADS']='32'; os.environ['OPENBLAS_NUM_THREADS']='32'
import scanpy as sc, anndata as ad, numpy as np, pandas as pd, numba, faiss
import time, glob, gc, warnings; warnings.filterwarnings('ignore')

SHARED = '/ssd/data/agent/bio/shared'
DATA_DIR = '/ssd/data/agent/bio/atlas_full/atlas_dataset'
CACHE = '/ssd/data/agent/bio/shared/cache_490m'
os.makedirs(CACHE, exist_ok=True); os.makedirs(SHARED, exist_ok=True)
K = 30; GPU_ID = 0

def gpu_knn_ivf(X, k, gpu_id=GPU_ID, nlist=1024, nprobe=32):
    X = np.ascontiguousarray(X, dtype=np.float32)
    n, d = X.shape
    res = faiss.StandardGpuResources()
    q = faiss.IndexFlatL2(d)
    ix = faiss.IndexIVFFlat(q, d, min(nlist, n//50))
    ig = faiss.index_cpu_to_gpu(res, gpu_id, ix)
    ig.train(X); ig.add(X); ig.nprobe = nprobe
    chunk = 200000
    all_idx = np.zeros((n, k), dtype=np.int64)
    for s in range(0, n, chunk):
        e = min(s+chunk, n)
        _, idx = ig.search(X[s:e], k+1)
        all_idx[s:e] = idx[:, 1:]
    return all_idx

@numba.njit(parallel=True)
def np_numba(pre, post, k):
    n = pre.shape[0]; out = np.zeros(n, dtype=np.float32)
    for i in numba.prange(n):
        c = 0
        for a in range(k):
            for b in range(k):
                if pre[i,a]==post[i,b]: c+=1; break
        out[i] = c/k
    return out

t_total = time.time()

# Step 1: Load all datasets, get PCA (cached)
pca_cache = os.path.join(CACHE, 'X_pca_490m.npy')
obs_cache = os.path.join(CACHE, 'obs_490m.pkl')

if os.path.exists(pca_cache) and os.path.exists(obs_cache):
    print("Loading cached PCA...")
    X_pca = np.load(pca_cache)
    obs_full = pd.read_pickle(obs_cache)
    print(f"Cached: {X_pca.shape}, {time.time()-t_total:.0f}s")
else:
    print("Step 1: Load + preprocess (will cache)")
    t0 = time.time()
    files = sorted(glob.glob(os.path.join(DATA_DIR, '*.h5ad')))
    adatas = [ad.read_h5ad(f) for f in files]
    adata = ad.concat(adatas, join='inner')
    del adatas; gc.collect()
    n = adata.n_obs
    print(f"Loaded: {n:,} cells, {time.time()-t0:.0f}s")

    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata_hvg = adata[:, adata.var['highly_variable']].copy()

    # Sparse SVD (no densify)
    from sklearn.decomposition import TruncatedSVD
    from scipy.sparse import issparse
    X_hvg = adata_hvg.X
    if not issparse(X_hvg):
        from scipy.sparse import csr_matrix; X_hvg = csr_matrix(X_hvg)
    svd = TruncatedSVD(n_components=50, random_state=42)
    X_pca = svd.fit_transform(X_hvg).astype(np.float32)
    obs_full = adata.obs[['Dataset','Celltype']].copy()

    # Cache
    np.save(pca_cache, X_pca)
    obs_full.to_pickle(obs_cache)
    del adata, adata_hvg, X_hvg; gc.collect()
    print(f"Preprocessed + cached: {X_pca.shape}, {time.time()-t0:.0f}s")

n = X_pca.shape[0]
batches = obs_full['Dataset'].values
ct = obs_full['Celltype'].values
B = len(np.unique(batches))
print(f"Cells: {n:,}, Batches: {B}")

# Step 2: Harmony (full, on 50-dim PCA)
print("\nStep 2: Harmony")
t0 = time.time()
import harmonypy as hm
ho = hm.run_harmony(X_pca, obs_full, 'Dataset')
Z = np.array(ho.Z_corr, dtype=np.float32)
if Z.shape[0] != n: Z = Z.T
print(f"Harmony: {time.time()-t0:.0f}s")

# Step 3: BD + blend
print("\nStep 3: BD")
t0 = time.time()
knn = gpu_knn_ivf(X_pca, K)
batch_int = pd.Categorical(batches).codes.astype(np.int64)

@numba.njit(parallel=True)
def bd_numba(knn, batch_ids, B, k):
    n = knn.shape[0]; out = np.zeros(n, dtype=np.float32)
    for i in numba.prange(n):
        seen = set()
        for j in range(k): seen.add(batch_ids[knn[i,j]])
        out[i] = len(seen) / min(k, B)
    return out

@numba.njit(parallel=True)
def smooth_numba(bd, knn, k):
    n = bd.shape[0]; out = np.zeros(n, dtype=np.float32)
    for i in numba.prange(n):
        s = 0.0
        for j in range(k): s += bd[knn[i,j]]
        out[i] = s / k
    return out

BD = bd_numba(knn, batch_int, B, K)
BD_s = smooth_numba(BD, knn, K)
d = Z.shape[1]
Z_rasi = BD_s[:, None] * Z + (1 - BD_s[:, None]) * X_pca[:, :d]
print(f"BD + blend: {time.time()-t0:.0f}s")

# Step 4: NP
print("\nStep 4: NP")
t0 = time.time()
pre_idx = knn  # reuse from BD step
post_idx = gpu_knn_ivf(Z_rasi, K)
np_scores = np_numba(pre_idx, post_idx, K)
print(f"NP: {time.time()-t0:.0f}s, mean={np_scores.mean():.4f}")

# Step 5: Results
print("\nStep 5: Results")
np_df = pd.DataFrame({'NP': np_scores, 'celltype': ct, 'dataset': batches})
print("\nNP by celltype:")
print(np_df.groupby('celltype')['NP'].agg(['mean','count']).sort_values('mean').round(4))

np_df.to_csv(f'{SHARED}/rasi_490m_np.csv', index=False, float_format='%.4f')
print(f"\nTotal: {time.time()-t_total:.0f}s ({(time.time()-t_total)/60:.0f} min)")
print(f"Cells: {n:,}")
print("490M RASI complete!")
