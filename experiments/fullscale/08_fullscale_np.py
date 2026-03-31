"""Standard 3: Full-scale NP on 2.25M immune cells.
Uses FAISS GPU for kNN, Harmony on GPU.
Memory-efficient: chunked processing.
"""
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import faiss
import warnings, os, time, gc
warnings.filterwarnings('ignore')

OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/fullscale/results'
os.makedirs(OUT_DIR, exist_ok=True)

K = 30
GPU_ID = 0

def faiss_knn_gpu_chunked(X, k, gpu_id=GPU_ID, chunk_size=100000):
    """Memory-efficient GPU kNN: build index once, query in chunks."""
    X = np.ascontiguousarray(X, dtype=np.float32)
    n, d = X.shape
    res = faiss.StandardGpuResources()
    index = faiss.index_cpu_to_gpu(res, gpu_id, faiss.IndexFlatL2(d))
    index.add(X)
    all_indices = np.zeros((n, k), dtype=np.int64)
    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        _, idx = index.search(X[start:end], k + 1)
        all_indices[start:end] = idx[:, 1:]
        if start % 500000 == 0:
            print(f"    kNN: {start}/{n}")
    return all_indices

# ============================================================
# 1. Load full atlas
# ============================================================
print("=" * 60)
print("Step 1: Load full atlas")
t_total = time.time()

adata = ad.read_h5ad('/ssd/data/agent/bio/atlas_merged_immune.h5ad')
print(f"Loaded: {adata.shape}")
print(f"Memory: ~{adata.X.data.nbytes / 1e9:.1f} GB" if hasattr(adata.X, 'data') else "dense")

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
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)
X_pca = adata_hvg.obsm['X_pca'].astype(np.float32)

print(f"Preprocessed in {time.time()-t0:.0f}s: {X_pca.shape}")

# Free memory
del adata_hvg
gc.collect()

# ============================================================
# 3. Pre-integration kNN (FAISS GPU)
# ============================================================
print("\nStep 3: Pre-integration kNN")
t0 = time.time()
pre_indices = faiss_knn_gpu_chunked(X_pca, K)
print(f"Pre-kNN done in {time.time()-t0:.0f}s")

# ============================================================
# 4. Harmony integration
# ============================================================
print("\nStep 4: Harmony integration")
t0 = time.time()
import harmonypy as hm
ho = hm.run_harmony(X_pca, adata.obs, 'Dataset')
Z = np.array(ho.Z_corr, dtype=np.float32)
if Z.shape[0] != adata.n_obs:
    Z = Z.T
print(f"Harmony done in {time.time()-t0:.0f}s")

# ============================================================
# 5. Post-integration kNN (FAISS GPU)
# ============================================================
print("\nStep 5: Post-integration kNN")
t0 = time.time()
post_indices = faiss_knn_gpu_chunked(Z, K)
print(f"Post-kNN done in {time.time()-t0:.0f}s")

# ============================================================
# 6. Compute global NP
# ============================================================
print("\nStep 6: Compute global NP")
t0 = time.time()

np_scores = np.zeros(adata.n_obs, dtype=np.float32)
chunk = 50000
for start in range(0, adata.n_obs, chunk):
    end = min(start + chunk, adata.n_obs)
    for i in range(start, end):
        pre_set = set(pre_indices[i])
        post_set = set(post_indices[i])
        np_scores[i] = len(pre_set & post_set) / K
    if start % 500000 == 0:
        print(f"    NP: {start}/{adata.n_obs}")

print(f"NP computed in {time.time()-t0:.0f}s")
print(f"NP: mean={np_scores.mean():.4f}, median={np.median(np_scores):.4f}")

# ============================================================
# 7. Results by celltype
# ============================================================
print("\nStep 7: Results")

np_df = pd.DataFrame({
    'NP': np_scores,
    'celltype': adata.obs['Celltype'].values,
    'dataset': adata.obs['Dataset'].values,
})

print("\nNP by celltype:")
ct_stats = np_df.groupby('celltype')['NP'].agg(['mean', 'median', 'std', 'count'])
ct_stats = ct_stats.sort_values('mean')
print(ct_stats.round(4))

# NP distribution stats
print(f"\nGlobal NP stats:")
print(f"  Mean: {np_scores.mean():.4f}")
print(f"  <0.3 (high risk): {(np_scores < 0.3).sum()} ({(np_scores < 0.3).mean()*100:.1f}%)")
print(f"  <0.2 (very high risk): {(np_scores < 0.2).sum()} ({(np_scores < 0.2).mean()*100:.1f}%)")
print(f"  <0.1 (extreme risk): {(np_scores < 0.1).sum()} ({(np_scores < 0.1).mean()*100:.1f}%)")

# NP by dataset
print("\nNP by dataset (bottom 10):")
ds_stats = np_df.groupby('dataset')['NP'].mean().sort_values()
print(ds_stats.head(10).round(4))

# Save (compact)
np_df.to_csv(f'{OUT_DIR}/fullscale_np_summary.csv', index=False,
             float_format='%.4f')

print(f"\nTotal time: {time.time()-t_total:.0f}s")
print(f"Cells: {adata.n_obs:,}")
print("Full-scale NP complete!")
