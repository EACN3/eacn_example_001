"""RASI 490M FAST: Random Projection + IVF kNN + numba. Target <20min."""
import sys; sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)
import os; os.environ['OMP_NUM_THREADS']='32'; os.environ['MKL_NUM_THREADS']='32'; os.environ['OPENBLAS_NUM_THREADS']='32'
import scanpy as sc, anndata as ad, numpy as np, pandas as pd, numba, faiss
from sklearn.random_projection import SparseRandomProjection
from scipy.sparse import issparse
import time, glob, gc, warnings; warnings.filterwarnings('ignore')

SHARED = '/ssd/data/agent/bio/shared'
CACHE = '/ssd/data/agent/bio/shared/cache_490m'
DATA_DIR = '/ssd/data/agent/bio/atlas_full/atlas_dataset'
os.makedirs(CACHE, exist_ok=True)
K = 30; GPU_ID = 0

def gpu_knn_ivf(X, k, gpu_id=GPU_ID):
    X = np.ascontiguousarray(X, dtype=np.float32)
    n, d = X.shape
    res = faiss.StandardGpuResources()
    q = faiss.IndexFlatL2(d); ix = faiss.IndexIVFFlat(q, d, min(1024, n//50))
    ig = faiss.index_cpu_to_gpu(res, gpu_id, ix)
    ig.train(X); ig.add(X); ig.nprobe = 32
    chunk = 200000; all_idx = np.zeros((n, k), dtype=np.int64)
    for s in range(0, n, chunk):
        e = min(s+chunk, n); _, idx = ig.search(X[s:e], k+1); all_idx[s:e] = idx[:,1:]
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

t_total = time.time()

# Step 1: Load + Random Projection (cached)
rp_cache = os.path.join(CACHE, 'X_rp50.npy')
obs_cache = os.path.join(CACHE, 'obs_490m.pkl')

if os.path.exists(rp_cache) and os.path.exists(obs_cache):
    print("Loading cached RP...")
    X_rp = np.load(rp_cache)
    obs_full = pd.read_pickle(obs_cache)
    print(f"Cached: {X_rp.shape}, {time.time()-t_total:.0f}s")
else:
    print("Step 1: Load + preprocess")
    t0 = time.time()
    files = sorted(glob.glob(os.path.join(DATA_DIR, '*.h5ad')))
    adatas = [ad.read_h5ad(f) for f in files]
    adata = ad.concat(adatas, join='inner')
    del adatas; gc.collect()
    print(f"Loaded: {adata.n_obs:,} cells, {time.time()-t0:.0f}s")

    t0 = time.time()
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    X_hvg = adata[:, adata.var['highly_variable']].X
    if not issparse(X_hvg):
        from scipy.sparse import csr_matrix; X_hvg = csr_matrix(X_hvg)

    # Random Projection: ~30s for 4.8M (vs 16min SVD)
    rp = SparseRandomProjection(n_components=50, random_state=42)
    X_rp = rp.fit_transform(X_hvg).astype(np.float32)
    obs_full = adata.obs[['Dataset','Celltype']].copy()
    np.save(rp_cache, X_rp); obs_full.to_pickle(obs_cache)
    del adata, X_hvg; gc.collect()
    print(f"RP done: {X_rp.shape}, {time.time()-t0:.0f}s")

n = X_rp.shape[0]
batches = obs_full['Dataset'].values
batch_int = pd.Categorical(batches).codes.astype(np.int64)
ct = obs_full['Celltype'].values
B = len(np.unique(batches))
print(f"Cells: {n:,}, Batches: {B}")

# Step 2: Pre-kNN + BD
print("\nStep 2: kNN + BD")
t0 = time.time()
knn = gpu_knn_ivf(X_rp, K)
BD = bd_numba(knn, batch_int, B, K)
BD_s = smooth_numba(BD, knn, K)
print(f"kNN+BD: {time.time()-t0:.0f}s")

# Step 3: Harmony
print("\nStep 3: Harmony")
t0 = time.time()
import harmonypy as hm
Z_harm = np.array(hm.run_harmony(X_rp, obs_full, 'Dataset').Z_corr, dtype=np.float32)
if Z_harm.shape[0] != n: Z_harm = Z_harm.T
print(f"Harmony: {time.time()-t0:.0f}s")

# Step 4: BD blend
d = Z_harm.shape[1]
Z_rasi = BD_s[:,None]*Z_harm + (1-BD_s[:,None])*X_rp[:,:d]

# Step 5: Post-kNN + NP
print("\nStep 5: NP")
t0 = time.time()
post_knn = gpu_knn_ivf(Z_rasi, K)
np_scores = np_numba(knn, post_knn, K)
print(f"NP: {time.time()-t0:.0f}s, mean={np_scores.mean():.4f}")

# Results
print("\nResults:")
np_df = pd.DataFrame({'NP': np_scores, 'celltype': ct, 'dataset': batches})
print(np_df.groupby('celltype')['NP'].agg(['mean','count']).sort_values('mean').round(4))
np_df.to_csv(f'{SHARED}/rasi_490m_np.csv', index=False, float_format='%.4f')

total = time.time() - t_total
print(f"\nTotal: {total:.0f}s ({total/60:.1f} min)")
print(f"Cells: {n:,}")
print("490M RASI FAST complete!")
