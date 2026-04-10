"""490M: parallel file reading (8 workers) + single pass. Full logging."""
import os, sys, time
sys.stdout.reconfigure(line_buffering=True)
os.environ['OMP_NUM_THREADS']='4'  # per-worker threads

import numpy as np, pandas as pd, numba, faiss
from scipy import sparse
from sklearn.random_projection import SparseRandomProjection
import anndata as ad, glob, gc, warnings
from multiprocessing import Pool
warnings.filterwarnings('ignore')

def log(msg): print(f"[{time.time()-T0:.0f}s] {msg}", flush=True)
T0 = time.time()
SHARED = '/ssd/data/agent/bio/shared'; CACHE = '/ssd/data/agent/bio/shared/cache_490m'
os.makedirs(CACHE, exist_ok=True)
DATA_DIR = '/ssd/data/agent/bio/atlas_full/atlas_dataset'
K = 30; GPU_DEV = 0

rp_cache = os.path.join(CACHE, 'X_rp50_v7.npy')
obs_cache = os.path.join(CACHE, 'obs_490m_v7.pkl')

def process_file(f):
    """Single-pass: read, normalize, compute stats, extract HVG. Returns (stats, obs, hvg_X)."""
    a = ad.read_h5ad(f)
    X = a.X.tocsr().copy()
    rs = np.array(X.sum(axis=1)).flatten(); rs[rs==0]=1
    # Fast CSR row scale
    for i in range(X.shape[0]):
        X.data[X.indptr[i]:X.indptr[i+1]] *= 1e4/rs[i]
    X.data = np.log1p(X.data)
    gene_sum = np.array(X.sum(axis=0)).flatten()
    gene_sqsum = np.array(X.power(2).sum(axis=0)).flatten()
    obs = a.obs[['Dataset','Celltype']].copy()
    return {'gene_sum': gene_sum, 'gene_sqsum': gene_sqsum, 'n': a.n_obs,
            'obs': obs, 'X': X, 'var_names': list(a.var_names)}

if os.path.exists(rp_cache) and os.path.exists(obs_cache):
    log("Loading cached RP..."); X_rp = np.load(rp_cache); obs_full = pd.read_pickle(obs_cache)
    log(f"Loaded: {X_rp.shape}")
else:
    files = sorted(glob.glob(os.path.join(DATA_DIR, '*.h5ad')))
    log(f"Single-pass parallel processing {len(files)} files (8 workers)")

    # Parallel read + process
    with Pool(8) as pool:
        results = []
        for i, r in enumerate(pool.imap_unordered(process_file, files)):
            results.append(r)
            if (i+1) % 10 == 0:
                n_so_far = sum(r2['n'] for r2 in results)
                log(f"  {i+1}/{len(files)}, {n_so_far:,} cells")

    # Aggregate gene stats
    log("Aggregating stats + HVG selection...")
    total = sum(r['n'] for r in results)
    gene_sum = sum(r['gene_sum'] for r in results)
    gene_sqsum = sum(r['gene_sqsum'] for r in results)
    gene_mean = gene_sum / total
    gene_var = gene_sqsum / total - gene_mean**2
    hvg_idx = np.argsort(gene_var)[-2000:]

    # Extract HVG columns and concat
    log("Extracting HVG + concat...")
    hvg_chunks = [r['X'][:, hvg_idx] for r in results]
    X_hvg = sparse.vstack(hvg_chunks)
    obs_full = pd.concat([r['obs'] for r in results], ignore_index=True)
    del results, hvg_chunks; gc.collect()
    log(f"HVG matrix: {X_hvg.shape}")

    # Random Projection
    log("Random Projection...")
    rp = SparseRandomProjection(n_components=50, random_state=42)
    X_rp = rp.fit_transform(X_hvg).astype(np.float32)
    del X_hvg; gc.collect()
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
log("NP..."); np_s=np_n(knn,post_knn,K); log(f"NP done: mean={np_s.mean():.4f}")

df = pd.DataFrame({'NP': np_s, 'celltype': ct, 'dataset': batches})
print(df.groupby('celltype')['NP'].agg(['mean','count']).sort_values('mean').round(4), flush=True)
df.to_csv(f'{SHARED}/rasi_490m_np.csv', index=False, float_format='%.4f')
log(f"DONE! {n:,} cells, {(time.time()-T0)/60:.1f}min")
