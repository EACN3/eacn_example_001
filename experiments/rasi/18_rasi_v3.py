"""RASI v3: 50-dim Harmony + PCA+UMAP concat for BD/NP. Target <1min."""
import sys; sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)
import os; os.environ['OMP_NUM_THREADS']='32'; os.environ['MKL_NUM_THREADS']='32'; os.environ['OPENBLAS_NUM_THREADS']='32'
import scanpy as sc, anndata as ad, numpy as np, pandas as pd, numba, faiss
from sklearn.metrics import silhouette_score
import time, warnings; warnings.filterwarnings('ignore')

SHARED = '/ssd/data/agent/bio/shared'
K = 30; GPU_ID = 0

def gpu_knn(X, k, gpu_id=GPU_ID):
    X = np.ascontiguousarray(X, dtype=np.float32)
    n, d = X.shape
    res = faiss.StandardGpuResources()
    if n > 500000:
        q = faiss.IndexFlatL2(d)
        ix = faiss.IndexIVFFlat(q, d, min(1024, n//100))
        ig = faiss.index_cpu_to_gpu(res, gpu_id, ix)
        ig.train(X); ig.add(X); ig.nprobe = 32
    else:
        ig = faiss.index_cpu_to_gpu(res, gpu_id, faiss.IndexFlatL2(d))
        ig.add(X)
    _, idx = ig.search(X, k + 1)
    return idx[:, 1:]

@numba.njit(parallel=True)
def np_numba(pre, post, k):
    n = pre.shape[0]
    out = np.zeros(n, dtype=np.float32)
    for i in numba.prange(n):
        count = 0
        for a in range(k):
            for b in range(k):
                if pre[i, a] == post[i, b]:
                    count += 1; break
        out[i] = count / k
    return out

@numba.njit(parallel=True)
def bd_numba(knn, batch_ids, B, k):
    n = knn.shape[0]
    out = np.zeros(n, dtype=np.float32)
    for i in numba.prange(n):
        seen = set()
        for j in range(k):
            seen.add(batch_ids[knn[i, j]])
        out[i] = len(seen) / min(k, B)
    return out

@numba.njit(parallel=True)
def smooth_bd(bd, knn, k):
    n = bd.shape[0]
    out = np.zeros(n, dtype=np.float32)
    for i in numba.prange(n):
        s = 0.0
        for j in range(k): s += bd[knn[i, j]]
        out[i] = s / k
    return out

t_total = time.time()

# Load
print("Load...")
t0 = time.time()
adata = ad.read_h5ad('experiments/np_guard/data/subset.h5ad')
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
batches = adata.obs['Dataset'].values
batch_int = pd.Categorical(batches).codes.astype(np.int64)
ct = adata.obs['Celltype'].values
B = len(np.unique(batches)); n = adata.n_obs
markers = ['CLEC4C','IL3RA','IRF7','CCR8','CFTR','ASCL3','GHRL','SOX10','AXL','SIGLEC6']
print(f"Load: {time.time()-t0:.0f}s")

# === A: Standard ===
print("\n=== A: Standard ===")
t0 = time.time()
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')
hvg_std = set(adata.var_names[adata.var['highly_variable']])
adata_a = adata[:, list(hvg_std)].copy()
sc.pp.scale(adata_a, max_value=10)
sc.pp.pca(adata_a, n_comps=50, random_state=42)
X_a = adata_a.obsm['X_pca'].astype(np.float32)
import harmonypy as hm
Z_a = np.array(hm.run_harmony(X_a, adata.obs, 'Dataset').Z_corr, dtype=np.float32)
if Z_a.shape[0] != n: Z_a = Z_a.T
pre_a = gpu_knn(X_a, K); post_a = gpu_knn(Z_a, K)
np_a = np_numba(pre_a, post_a, K)
cov_a = sum(1 for m in markers if m in hvg_std and m in adata.var_names)
s = np.random.RandomState(42).choice(n, min(5000,n), replace=False)
asw_a = silhouette_score(Z_a[s], batches[s])
t_a = time.time() - t0
print(f"Standard: NP={np_a.mean():.4f}, ASW={asw_a:.3f}, markers={cov_a}/10, time={t_a:.0f}s")

# === B: RASI v3 ===
print("\n=== B: RASI v3 ===")
t0 = time.time()

# Step 1: Rare-HVG (t-test, fast)
sc.pp.neighbors(adata_a, n_pcs=50, random_state=42)
sc.tl.leiden(adata_a, resolution=0.5, key_added='coarse', random_state=42)
adata.obs['coarse'] = adata_a.obs['coarse']
sc.tl.rank_genes_groups(adata, 'coarse', method='t-test', key_added='rg')
rare_genes = set()
for cl in adata.obs['coarse'].cat.categories:
    try: rare_genes.update(sc.get.rank_genes_groups_df(adata, group=cl, key='rg').head(50)['names'])
    except: pass
hvg_rasi = hvg_std | (rare_genes & set(adata.var_names))
cov_b = sum(1 for m in markers if m in hvg_rasi and m in adata.var_names)
rescued = [m for m in markers if m in hvg_rasi and m not in hvg_std and m in adata.var_names]
print(f"Rare-HVG: {len(hvg_rasi)}, markers={cov_b}/10, rescued={rescued}")

# Step 2: PCA 50 (for Harmony) + UMAP 10 (for BD/NP)
adata_b = adata[:, list(hvg_rasi)].copy()
sc.pp.scale(adata_b, max_value=10)
sc.pp.pca(adata_b, n_comps=50, random_state=42)
sc.pp.neighbors(adata_b, n_pcs=50, random_state=42)
sc.tl.umap(adata_b, random_state=42, n_components=10)
X_pca50 = adata_b.obsm['X_pca'].astype(np.float32)
X_hybrid = np.hstack([X_pca50, adata_b.obsm['X_umap'].astype(np.float32)])  # 60 dims
print(f"Hybrid embedding: {X_hybrid.shape}")

# Step 3: Harmony on 50-dim PCA (fast), BD/NP on 60-dim hybrid
Z_harm = np.array(hm.run_harmony(X_pca50, adata.obs, 'Dataset').Z_corr, dtype=np.float32)
if Z_harm.shape[0] != n: Z_harm = Z_harm.T

knn_h = gpu_knn(X_hybrid, K)
BD = bd_numba(knn_h, batch_int, B, K)
BD_s = smooth_bd(BD, knn_h, K)

d = Z_harm.shape[1]
Z_b = BD_s[:, None] * Z_harm + (1 - BD_s[:, None]) * X_pca50[:, :d]

# Step 4: NP on hybrid space
pre_b = gpu_knn(X_hybrid, K)
post_b = gpu_knn(Z_b, K)
np_b = np_numba(pre_b, post_b, K)
asw_b = silhouette_score(Z_b[s], batches[s])
t_b = time.time() - t0
print(f"RASI v3: NP={np_b.mean():.4f}, ASW={asw_b:.3f}, markers={cov_b}/10, time={t_b:.0f}s")

# Step 5: Low-NP discovery
thr = np.percentile(np_b, 5)
low = np_b < thr
adata.obs['np_grp'] = pd.Categorical(np.where(low, 'low_NP', 'normal'))
print(f"Low-NP: {low.sum()} cells, types: {dict(pd.Series(ct[low]).value_counts().head(5))}")

print(f"\n{'='*50}")
print(f"| Metric  | Standard | RASI v3  |")
print(f"|---------|----------|----------|")
print(f"| NP      | {np_a.mean():.4f}   | {np_b.mean():.4f}   |")
print(f"| ASW     | {asw_a:.3f}    | {asw_b:.3f}    |")
print(f"| Markers | {cov_a}/10     | {cov_b}/10     |")
print(f"| Time    | {t_a:.0f}s      | {t_b:.0f}s      |")
print(f"Total: {time.time()-t_total:.0f}s")
