"""490M: fast CSR normalize (0.1s/file vs 17s). Full logging."""
import os, sys, time
sys.stdout.reconfigure(line_buffering=True)
os.environ['OMP_NUM_THREADS']='32'
import numpy as np, pandas as pd, numba, faiss, anndata as ad
from scipy import sparse
from sklearn.random_projection import SparseRandomProjection
import glob, gc, warnings; warnings.filterwarnings('ignore')

def log(msg): print(f"[{time.time()-T0:.0f}s] {msg}", flush=True)
T0 = time.time()
SHARED = '/ssd/data/agent/bio/shared'; CACHE = '/ssd/data/agent/bio/shared/cache_490m'
os.makedirs(CACHE, exist_ok=True)
DATA_DIR = '/ssd/data/agent/bio/atlas_full/atlas_dataset'
K = 30; GPU_DEV = 0

rp_cache = os.path.join(CACHE, 'X_rp50_v6.npy')
obs_cache = os.path.join(CACHE, 'obs_490m_v6.pkl')

def fast_normalize_csr(X):
    """Normalize CSR rows to sum=1e4, then log1p. O(nnz), ~0.1s/25k cells."""
    X = X.tocsr().copy()
    rs = np.array(X.sum(axis=1)).flatten(); rs[rs==0]=1
    for i in range(X.shape[0]):
        X.data[X.indptr[i]:X.indptr[i+1]] *= 1e4/rs[i]
    X.data = np.log1p(X.data)
    return X

if os.path.exists(rp_cache) and os.path.exists(obs_cache):
    log("Loading cached RP..."); X_rp = np.load(rp_cache); obs_full = pd.read_pickle(obs_cache)
    log(f"Loaded: {X_rp.shape}")
else:
    files = sorted(glob.glob(os.path.join(DATA_DIR, '*.h5ad')))
    log(f"Pass 1: gene stats from {len(files)} files (fast CSR normalize)")
    gene_sum = None; gene_sqsum = None; total = 0; gnames = None; all_obs = []
    for i, f in enumerate(files):
        a = ad.read_h5ad(f); X = fast_normalize_csr(a.X)
        if gnames is None: gnames = a.var_names; gene_sum = np.zeros(len(gnames)); gene_sqsum = np.zeros(len(gnames))
        gene_sum += np.array(X.sum(axis=0)).flatten()
        gene_sqsum += np.array(X.power(2).sum(axis=0)).flatten()
        total += a.n_obs; all_obs.append(a.obs[['Dataset','Celltype']])
        del a, X
        if (i+1) % 10 == 0: log(f"  {i+1}/{len(files)}, {total:,} cells")
    gene_mean = gene_sum/total; gene_var = gene_sqsum/total - gene_mean**2
    hvg_idx = np.argsort(gene_var)[-2000:]
    obs_full = pd.concat(all_obs, ignore_index=True)
    log(f"Pass 1 done: {total:,} cells, HVG selected")

    log("Pass 2: extract HVG + RP")
    chunks = []
    for i, f in enumerate(files):
        a = ad.read_h5ad(f); X = fast_normalize_csr(a.X)[:, hvg_idx]
        chunks.append(X); del a
        if (i+1) % 10 == 0: log(f"  {i+1}/{len(files)}")
    X_hvg = sparse.vstack(chunks); del chunks; gc.collect()
    log(f"HVG matrix: {X_hvg.shape}")
    rp = SparseRandomProjection(n_components=50, random_state=42)
    X_rp = rp.fit_transform(X_hvg).astype(np.float32); del X_hvg; gc.collect()
    log(f"RP: {X_rp.shape}")
    np.save(rp_cache, X_rp); obs_full.to_pickle(obs_cache)
    log("Cached!")

n = X_rp.shape[0]; batches = obs_full['Dataset'].values
batch_int = pd.Categorical(batches).codes.astype(np.int64)
ct = obs_full['Celltype'].values; B = len(np.unique(batches))
log(f"Ready: {n:,} cells, {B} batches")

def gpu_knn_ivf(X, k):
    X = np.ascontiguousarray(X, dtype=np.float32); n, d = X.shape
    res = faiss.StandardGpuResources(); q = faiss.IndexFlatL2(d)
    ix = faiss.IndexIVFFlat(q, d, min(1024, n//50))
    ig = faiss.index_cpu_to_gpu(res, GPU_DEV, ix); ig.train(X); ig.add(X); ig.nprobe = 32
    r = np.zeros((n, k), dtype=np.int64)
    for s in range(0, n, 200000):
        e = min(s+200000, n); _, idx = ig.search(X[s:e], k+1); r[s:e] = idx[:,1:]
    return r

log("Pre-kNN..."); knn = gpu_knn_ivf(X_rp, K); log("Pre-kNN done")

@numba.njit(parallel=True)
def bd_n(knn,bi,B,k):
    n=knn.shape[0];o=np.zeros(n,dtype=np.float32)
    for i in numba.prange(n):
        s=set()
        for j in range(k):s.add(bi[knn[i,j]])
        o[i]=len(s)/min(k,B)
    return o
@numba.njit(parallel=True)
def sm(bd,knn,k):
    n=bd.shape[0];o=np.zeros(n,dtype=np.float32)
    for i in numba.prange(n):
        s=0.
        for j in range(k):s+=bd[knn[i,j]]
        o[i]=s/k
    return o
log("BD..."); BD=bd_n(knn,batch_int,B,K); BD_s=sm(BD,knn,K); log("BD done")

log("Harmony...")
import harmonypy as hm
Z = np.array(hm.run_harmony(X_rp, obs_full, 'Dataset').Z_corr, dtype=np.float32)
if Z.shape[0] != n: Z = Z.T
log("Harmony done")

Z_r = BD_s[:,None]*Z + (1-BD_s[:,None])*X_rp[:,:Z.shape[1]]
log("Post-kNN..."); post_knn = gpu_knn_ivf(Z_r, K); log("Post-kNN done")

@numba.njit(parallel=True)
def np_n(pre,post,k):
    n=pre.shape[0];o=np.zeros(n,dtype=np.float32)
    for i in numba.prange(n):
        c=0
        for a in range(k):
            for b in range(k):
                if pre[i,a]==post[i,b]:c+=1;break
        o[i]=c/k
    return o
log("NP..."); np_s = np_n(knn, post_knn, K); log(f"NP done: mean={np_s.mean():.4f}")

df = pd.DataFrame({'NP': np_s, 'celltype': ct, 'dataset': batches})
print(df.groupby('celltype')['NP'].agg(['mean','count']).sort_values('mean').round(4), flush=True)
df.to_csv(f'{SHARED}/rasi_490m_np.csv', index=False, float_format='%.4f')
log(f"DONE! {n:,} cells, {(time.time()-T0)/60:.1f}min")
