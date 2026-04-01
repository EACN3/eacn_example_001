"""NP vs batch count: Phase transition curve.
For each batch count (5,10,20,30,50,70,103), subsample batches,
run Harmony + NP, compute per-celltype mean NP.
GPU-accelerated.
"""
import sys
sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)  # line-buffered
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import faiss
import harmonypy as hm
import warnings, os, time, gc
warnings.filterwarnings('ignore')

OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/np_vs_batches/results'
os.makedirs(OUT_DIR, exist_ok=True)

K = 30
GPU_ID = 0
SEED = 42

def faiss_knn_gpu_chunked(X, k, gpu_id=GPU_ID, chunk_size=200000):
    X = np.ascontiguousarray(X, dtype=np.float32)
    n, d = X.shape
    res = faiss.StandardGpuResources()
    index = faiss.index_cpu_to_gpu(res, gpu_id, faiss.IndexFlatL2(d))
    index.add(X)
    all_indices = np.zeros((n, k), dtype=np.int64)
    for s in range(0, n, chunk_size):
        e = min(s + chunk_size, n)
        _, idx = index.search(X[s:e], k + 1)
        all_indices[s:e] = idx[:, 1:]
    return all_indices

def compute_np_fast(pre_indices, post_indices, k=K):
    n = pre_indices.shape[0]
    np_scores = np.zeros(n, dtype=np.float32)
    for i in range(n):
        np_scores[i] = len(set(pre_indices[i]) & set(post_indices[i])) / k
    return np_scores

# ============================================================
# 1. Load full atlas (preprocessed)
# ============================================================
print("=" * 60)
print("Loading and preprocessing full atlas...")
t_total = time.time()

adata = ad.read_h5ad('/ssd/data/agent/bio/atlas_merged_immune.h5ad')
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)
X_pca_full = adata_hvg.obsm['X_pca'].astype(np.float32)
celltypes_full = adata.obs['Celltype'].values
datasets_full = adata.obs['Dataset'].values

all_batches = sorted(adata.obs['Dataset'].unique())
print(f"Full atlas: {adata.n_obs} cells, {len(all_batches)} batches")
print(f"Preprocessing done in {time.time()-t_total:.0f}s")

del adata_hvg
gc.collect()

# ============================================================
# 2. Run NP for each batch count
# ============================================================
BATCH_COUNTS = [5, 10, 20, 30, 50, 70, len(all_batches)]
all_results = []

np.random.seed(SEED)

for n_batches in BATCH_COUNTS:
    print(f"\n{'='*60}")
    print(f"Batch count: {n_batches}")
    t0 = time.time()

    if n_batches >= len(all_batches):
        # Reuse existing fullscale result
        fullscale_csv = '/ssd/data/agent/bio/eacn_example_001/experiments/fullscale/results/fullscale_np_summary.csv'
        if os.path.exists(fullscale_csv):
            print("  Reusing fullscale result for 103 batches")
            fs = pd.read_csv(fullscale_csv)
            ct_means = fs.groupby('celltype')['NP'].agg(['mean', 'median', 'count'])
            for celltype, row in ct_means.iterrows():
                all_results.append({
                    'n_batches': n_batches, 'n_cells': len(fs),
                    'celltype': celltype, 'np_mean': row['mean'],
                    'np_median': row['median'], 'n_cells_ct': int(row['count']),
                })
            pd.DataFrame(all_results).to_csv(f'{OUT_DIR}/np_vs_batches_partial.csv', index=False)
            print(f"  Done (reused)")
            continue
        selected = all_batches
    else:
        selected = list(np.random.choice(all_batches, size=n_batches, replace=False))

    # Subset
    mask = np.isin(datasets_full, selected)
    idx = np.where(mask)[0]
    X_pca = X_pca_full[idx]
    ct = celltypes_full[idx]
    ds = datasets_full[idx]
    # Ensure Dataset is plain string (not categorical with 103 levels)
    obs_df = pd.DataFrame({'Dataset': pd.Categorical(ds).remove_unused_categories()}, index=range(len(idx)))

    print(f"  Cells: {len(idx)}")

    # Pre-integration kNN
    pre_indices = faiss_knn_gpu_chunked(X_pca, K)

    # Harmony
    ho = hm.run_harmony(X_pca, obs_df, 'Dataset', verbose=False)
    Z = np.array(ho.Z_corr, dtype=np.float32)
    if Z.shape[0] != len(idx):
        Z = Z.T

    # Post-integration kNN
    post_indices = faiss_knn_gpu_chunked(Z, K)

    # NP
    np_scores = compute_np_fast(pre_indices, post_indices)

    # Per-celltype stats
    np_df = pd.DataFrame({'NP': np_scores, 'celltype': ct})
    ct_means = np_df.groupby('celltype')['NP'].agg(['mean', 'median', 'count'])

    for celltype, row in ct_means.iterrows():
        all_results.append({
            'n_batches': n_batches,
            'n_cells': len(idx),
            'celltype': celltype,
            'np_mean': row['mean'],
            'np_median': row['median'],
            'n_cells_ct': int(row['count']),
        })

    elapsed = time.time() - t0
    print(f"  NP mean: {np_scores.mean():.4f}")
    print(f"  pDC NP: {np_df[np_df['celltype']=='pDC']['NP'].mean():.4f}" if 'pDC' in ct else "  No pDC")
    print(f"  Time: {elapsed:.0f}s")

    # Incremental save after each point
    pd.DataFrame(all_results).to_csv(f'{OUT_DIR}/np_vs_batches_partial.csv', index=False)
    print(f"  Saved partial results ({len(all_results)} rows)")

    # Free memory
    del X_pca, Z, pre_indices, post_indices, np_scores
    gc.collect()

# ============================================================
# 3. Save and summarize
# ============================================================
results_df = pd.DataFrame(all_results)
results_df.to_csv(f'{OUT_DIR}/np_vs_batches.csv', index=False)

print("\n" + "=" * 60)
print("=== NP vs Batch Count Summary ===")
pivot = results_df.pivot_table(index='celltype', columns='n_batches', values='np_mean')
print(pivot.round(4).to_string())

# Highlight pDC and Mast
print("\n=== pDC NP trajectory ===")
pdc = results_df[results_df['celltype'] == 'pDC'].sort_values('n_batches')
print(pdc[['n_batches', 'n_cells_ct', 'np_mean']].to_string(index=False))

print("\n=== Mast NP trajectory ===")
mast = results_df[results_df['celltype'] == 'Mast'].sort_values('n_batches')
print(mast[['n_batches', 'n_cells_ct', 'np_mean']].to_string(index=False))

print(f"\nTotal time: {time.time()-t_total:.0f}s")
print("Phase transition analysis complete!")
