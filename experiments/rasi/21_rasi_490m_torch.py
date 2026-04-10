"""RASI 490M: torch GPU preprocessing + FAISS IVF + numba. Target <20min."""
import sys; sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)
import os; os.environ['OMP_NUM_THREADS']='32'; os.environ['MKL_NUM_THREADS']='32'
import numpy as np, pandas as pd, numba, faiss, torch
import anndata as ad, scipy.sparse as sp
import time, glob, gc, warnings; warnings.filterwarnings('ignore')

SHARED = '/ssd/data/agent/bio/shared'
CACHE = '/ssd/data/agent/bio/shared/cache_490m'
DATA_DIR = '/ssd/data/agent/bio/atlas_full/atlas_dataset'
os.makedirs(CACHE, exist_ok=True)
K = 30; GPU_DEV = 0  # CUDA_VISIBLE_DEVICES maps physical GPU

def gpu_knn_ivf(X, k, nlist=1024, nprobe=32):
    X = np.ascontiguousarray(X, dtype=np.float32)
    n, d = X.shape; res = faiss.StandardGpuResources()
    q = faiss.IndexFlatL2(d); ix = faiss.IndexIVFFlat(q, d, min(nlist, n//50))
    ig = faiss.index_cpu_to_gpu(res, GPU_DEV, ix)
    ig.train(X); ig.add(X); ig.nprobe = nprobe
    all_idx = np.zeros((n, k), dtype=np.int64)
    for s in range(0, n, 200000):
        e = min(s+200000, n); _, idx = ig.search(X[s:e], k+1); all_idx[s:e] = idx[:,1:]
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

# ============ Step 1: Load + GPU preprocess (cached) ============
rp_cache = os.path.join(CACHE, 'X_rp50_v2.npy')
obs_cache = os.path.join(CACHE, 'obs_490m_v2.pkl')

if os.path.exists(rp_cache) and os.path.exists(obs_cache):
    print("Loading cached...")
    X_rp = np.load(rp_cache); obs_full = pd.read_pickle(obs_cache)
    print(f"Cached: {X_rp.shape}, {time.time()-t_total:.0f}s")
else:
    print("Step 1: Load")
    t0 = time.time()
    files = sorted(glob.glob(os.path.join(DATA_DIR, '*.h5ad')))
    adatas = [ad.read_h5ad(f) for f in files]
    adata = ad.concat(adatas, join='inner')
    del adatas; gc.collect()
    n = adata.n_obs; genes = adata.var_names
    obs_full = adata.obs[['Dataset','Celltype']].copy()

    # Filter genes (CPU, fast)
    gene_counts = np.array((adata.X > 0).sum(axis=0)).flatten()
    keep = gene_counts >= 10
    X_sp = adata.X[:, keep]
    gene_names = genes[keep]
    del adata; gc.collect()
    print(f"Loaded: {n:,} cells, {X_sp.shape[1]} genes, {time.time()-t0:.0f}s")

    # GPU normalize + log1p in chunks
    print("GPU normalize + log1p...")
    t0 = time.time()
    if not sp.issparse(X_sp): X_sp = sp.csr_matrix(X_sp)

    # Compute per-gene mean/var for HVG (accumulate on GPU)
    chunk = 100000  # 100k x 35k = 13GB, fits in 15GB free
    gene_sum = np.zeros(X_sp.shape[1], dtype=np.float64)
    gene_sqsum = np.zeros(X_sp.shape[1], dtype=np.float64)

    # Also build normalized log1p sparse for RP
    norm_chunks = []
    for s in range(0, n, chunk):
        e = min(s+chunk, n)
        X_dense = torch.from_numpy(X_sp[s:e].toarray().astype(np.float32)).cuda()
        # normalize_total
        row_sums = X_dense.sum(dim=1, keepdim=True).clamp(min=1)
        X_dense = X_dense / row_sums * 1e4
        # log1p
        X_dense = torch.log1p(X_dense)
        # accumulate stats for HVG
        X_cpu = X_dense.cpu().numpy()
        gene_sum += X_cpu.sum(axis=0)
        gene_sqsum += (X_cpu**2).sum(axis=0)
        norm_chunks.append(sp.csr_matrix(X_cpu))
        del X_dense; torch.cuda.empty_cache()

    X_norm = sp.vstack(norm_chunks); del norm_chunks; gc.collect()
    print(f"GPU normalize+log1p: {time.time()-t0:.0f}s")

    # HVG: top 2000 by variance
    t0 = time.time()
    gene_mean = gene_sum / n
    gene_var = gene_sqsum / n - gene_mean**2
    hvg_idx = np.argsort(gene_var)[-2000:]
    X_hvg = X_norm[:, hvg_idx]
    del X_norm; gc.collect()
    print(f"HVG: {time.time()-t0:.0f}s")

    # Random Projection
    t0 = time.time()
    from sklearn.random_projection import SparseRandomProjection
    rp = SparseRandomProjection(n_components=50, random_state=42)
    X_rp = rp.fit_transform(X_hvg).astype(np.float32)
    del X_hvg; gc.collect()
    print(f"RP: {time.time()-t0:.0f}s, shape={X_rp.shape}")

    # Cache
    np.save(rp_cache, X_rp); obs_full.to_pickle(obs_cache)
    print(f"Preprocess total: {time.time()-t_total:.0f}s")

n = X_rp.shape[0]
batches = obs_full['Dataset'].values
batch_int = pd.Categorical(batches).codes.astype(np.int64)
ct = obs_full['Celltype'].values
B = len(np.unique(batches))
print(f"\nCells: {n:,}, Batches: {B}")

# ============ Step 2: kNN + BD ============
print("\nStep 2: kNN + BD")
t0 = time.time()
knn = gpu_knn_ivf(X_rp, K)
BD = bd_numba(knn, batch_int, B, K)
BD_s = smooth_numba(BD, knn, K)
print(f"kNN+BD: {time.time()-t0:.0f}s")

# ============ Step 3: Harmony ============
print("\nStep 3: Harmony")
t0 = time.time()
import harmonypy as hm
Z = np.array(hm.run_harmony(X_rp, obs_full, 'Dataset').Z_corr, dtype=np.float32)
if Z.shape[0] != n: Z = Z.T
print(f"Harmony: {time.time()-t0:.0f}s")

# ============ Step 4: BD blend ============
Z_rasi = BD_s[:,None]*Z + (1-BD_s[:,None])*X_rp[:,:Z.shape[1]]

# ============ Step 5: NP ============
print("\nStep 5: NP")
t0 = time.time()
post_knn = gpu_knn_ivf(Z_rasi, K)
np_scores = np_numba(knn, post_knn, K)
print(f"NP: {time.time()-t0:.0f}s, mean={np_scores.mean():.4f}")

# ============ Results ============
print("\nResults:")
np_df = pd.DataFrame({'NP': np_scores, 'celltype': ct, 'dataset': batches})
print(np_df.groupby('celltype')['NP'].agg(['mean','count']).sort_values('mean').round(4))
np_df.to_csv(f'{SHARED}/rasi_490m_np.csv', index=False, float_format='%.4f')

total = time.time() - t_total
print(f"\nTotal: {total:.0f}s ({total/60:.1f} min)")
print(f"Cells: {n:,}")
print("DONE!")
