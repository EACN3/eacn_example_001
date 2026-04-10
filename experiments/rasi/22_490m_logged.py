"""490M RASI with FULL logging on every step."""
import os, sys
# Force unbuffered
sys.stdout.reconfigure(line_buffering=True)
os.environ['OMP_NUM_THREADS']='32'; os.environ['MKL_NUM_THREADS']='32'

import numpy as np, pandas as pd, numba, faiss, anndata as ad
from scipy import sparse
from sklearn.random_projection import SparseRandomProjection
import time, gc, warnings; warnings.filterwarnings('ignore')

def log(msg): print(f"[{time.time()-T0:.0f}s] {msg}", flush=True)
T0 = time.time()

SHARED = '/ssd/data/agent/bio/shared'
CACHE = '/ssd/data/agent/bio/shared/cache_490m'
K = 30; GPU_DEV = 0

rp_cache = os.path.join(CACHE, 'X_rp50_v4.npy')
obs_cache = os.path.join(CACHE, 'obs_490m_v4.pkl')

if os.path.exists(rp_cache) and os.path.exists(obs_cache):
    log("Loading cached RP...")
    X_rp = np.load(rp_cache); obs_full = pd.read_pickle(obs_cache)
    log(f"Cached loaded: {X_rp.shape}")
else:
    log("Reading cached concat h5ad...")
    adata = ad.read_h5ad(os.path.join(CACHE, 'atlas_concat.h5ad'))
    log(f"Read done: {adata.shape}")

    log("Filter genes...")
    gc_arr = np.array((adata.X > 0).sum(axis=0)).flatten()
    X_sp = adata.X[:, gc_arr >= 10].tocsr()
    obs_full = adata.obs[['Dataset','Celltype']].copy()
    del adata; gc.collect()
    log(f"Filter done: {X_sp.shape}")

    log("Normalize (sparse diag)...")
    row_sums = np.array(X_sp.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1
    D = sparse.diags(1e4 / row_sums)
    X_sp = D @ X_sp
    log("Normalize done")

    log("Log1p...")
    X_sp.data = np.log1p(X_sp.data)
    log("Log1p done")

    log("HVG (variance)...")
    mean = np.array(X_sp.mean(axis=0)).flatten()
    sq_mean = np.array(X_sp.power(2).mean(axis=0)).flatten()
    var = sq_mean - mean**2
    hvg_idx = np.argsort(var)[-2000:]
    X_hvg = X_sp[:, hvg_idx]
    del X_sp; gc.collect()
    log(f"HVG done: {X_hvg.shape}")

    log("Random Projection...")
    rp = SparseRandomProjection(n_components=50, random_state=42)
    X_rp = rp.fit_transform(X_hvg).astype(np.float32)
    del X_hvg; gc.collect()
    log(f"RP done: {X_rp.shape}")

    log("Caching...")
    np.save(rp_cache, X_rp); obs_full.to_pickle(obs_cache)
    log("Cached!")

n = X_rp.shape[0]; batches = obs_full['Dataset'].values
batch_int = pd.Categorical(batches).codes.astype(np.int64)
ct = obs_full['Celltype'].values; B = len(np.unique(batches))
log(f"Ready: {n:,} cells, {B} batches")

log("Pre-kNN (FAISS IVF GPU)...")
def gpu_knn_ivf(X, k):
    X = np.ascontiguousarray(X, dtype=np.float32); n, d = X.shape
    res = faiss.StandardGpuResources()
    q = faiss.IndexFlatL2(d); ix = faiss.IndexIVFFlat(q, d, min(1024, n//50))
    ig = faiss.index_cpu_to_gpu(res, GPU_DEV, ix); ig.train(X); ig.add(X); ig.nprobe = 32
    all_idx = np.zeros((n, k), dtype=np.int64)
    for s in range(0, n, 200000):
        e = min(s+200000, n); _, idx = ig.search(X[s:e], k+1); all_idx[s:e] = idx[:,1:]
    return all_idx
knn = gpu_knn_ivf(X_rp, K)
log("Pre-kNN done")

log("BD...")
@numba.njit(parallel=True)
def bd_n(knn, bi, B, k):
    n=knn.shape[0]; o=np.zeros(n,dtype=np.float32)
    for i in numba.prange(n):
        s=set()
        for j in range(k): s.add(bi[knn[i,j]])
        o[i]=len(s)/min(k,B)
    return o
@numba.njit(parallel=True)
def sm(bd, knn, k):
    n=bd.shape[0]; o=np.zeros(n,dtype=np.float32)
    for i in numba.prange(n):
        s=0.
        for j in range(k): s+=bd[knn[i,j]]
        o[i]=s/k
    return o
BD=bd_n(knn,batch_int,B,K); BD_s=sm(BD,knn,K)
log("BD done")

log("Harmony...")
import harmonypy as hm
Z = np.array(hm.run_harmony(X_rp, obs_full, 'Dataset').Z_corr, dtype=np.float32)
if Z.shape[0] != n: Z = Z.T
log("Harmony done")

Z_r = BD_s[:,None]*Z + (1-BD_s[:,None])*X_rp[:,:Z.shape[1]]

log("Post-kNN...")
post_knn = gpu_knn_ivf(Z_r, K)
log("Post-kNN done")

log("NP...")
@numba.njit(parallel=True)
def np_n(pre,post,k):
    n=pre.shape[0]; o=np.zeros(n,dtype=np.float32)
    for i in numba.prange(n):
        c=0
        for a in range(k):
            for b in range(k):
                if pre[i,a]==post[i,b]: c+=1; break
        o[i]=c/k
    return o
np_s = np_n(knn, post_knn, K)
log(f"NP done: mean={np_s.mean():.4f}")

df = pd.DataFrame({'NP': np_s, 'celltype': ct, 'dataset': batches})
log("Results:")
print(df.groupby('celltype')['NP'].agg(['mean','count']).sort_values('mean').round(4), flush=True)
df.to_csv(f'{SHARED}/rasi_490m_np.csv', index=False, float_format='%.4f')
log(f"DONE! Total: {time.time()-T0:.0f}s ({(time.time()-T0)/60:.1f}min), {n:,} cells")
