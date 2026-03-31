"""Phase 1: CNEM validation on immune atlas subset (~105k cells).
GPU-accelerated: FAISS GPU for kNN, PyTorch CUDA for cosine similarity.
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
import torch
import faiss
import warnings, os, time
warnings.filterwarnings('ignore')

DATA_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/phase1/data'
OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/phase1/results'
os.makedirs(OUT_DIR, exist_ok=True)

K = 30
GPU_ID = 0

def faiss_knn_gpu(X, k, gpu_id=GPU_ID):
    """GPU-accelerated kNN using FAISS."""
    X = np.ascontiguousarray(X, dtype=np.float32)
    d = X.shape[1]
    index = faiss.IndexFlatL2(d)
    gpu_res = faiss.StandardGpuResources()
    index_gpu = faiss.index_cpu_to_gpu(gpu_res, gpu_id, index)
    index_gpu.add(X)
    distances, indices = index_gpu.search(X, k + 1)
    return indices[:, 1:]  # exclude self

def compute_cnem_gpu(X_expr_torch, indices):
    """Compute CNEM using PyTorch on GPU. Fully vectorized."""
    n, k = indices.shape
    # Gather neighbor expressions: (n, k, d)
    nbr_expr = X_expr_torch[indices]
    # Mean neighbor expression: (n, d)
    mean_nbr = nbr_expr.mean(dim=1)
    # Cosine similarity between each cell and its mean neighbor
    cos_sim = torch.nn.functional.cosine_similarity(X_expr_torch, mean_nbr, dim=1)
    cnem = 1.0 - cos_sim
    return cnem.cpu().numpy()

# ============================================================
# 1. Load and preprocess
# ============================================================
print("=" * 60)
print("Step 1: Load and preprocess")
print("=" * 60)
t_total = time.time()

adata = ad.read_h5ad(f'{DATA_DIR}/immune_subset_phase1.h5ad')
print(f"Loaded: {adata.shape}")

adata.obs['celltype_true'] = adata.obs['Celltype'].copy()
rare_types = ['pDC', 'Mast']
adata.obs['is_rare'] = adata.obs['Celltype'].isin(rare_types).astype(int)

sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')

X_expr = adata[:, adata.var['highly_variable']].X.copy()
if hasattr(X_expr, 'toarray'):
    X_expr = X_expr.toarray()
X_expr = np.array(X_expr, dtype=np.float32)

adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

# Move expression to GPU
X_expr_torch = torch.from_numpy(X_expr).cuda(GPU_ID)

print(f"Cells: {adata.n_obs}, HVG: {X_expr.shape[1]}")
print(f"GPU: {torch.cuda.get_device_name(GPU_ID)}")

# ============================================================
# 2. Pre-integration per-batch CNEM (GPU)
# ============================================================
print("\n" + "=" * 60)
print("Step 2: Pre-integration per-batch CNEM (GPU)")
print("=" * 60)

batches = adata.obs['Dataset'].unique()
pre_cnem = np.zeros(adata.n_obs, dtype=np.float32)

for batch in batches:
    t0 = time.time()
    mask = adata.obs['Dataset'] == batch
    idx = np.where(mask)[0]
    X_pca_b = adata.obsm['X_pca'][idx]
    X_expr_b = X_expr_torch[idx]

    k_use = min(K, len(idx) - 1)
    indices_local = faiss_knn_gpu(X_pca_b, k_use)
    indices_local_torch = torch.from_numpy(indices_local.astype(np.int64)).cuda(GPU_ID)

    cnem_b = compute_cnem_gpu(X_expr_b, indices_local_torch)
    pre_cnem[idx] = cnem_b
    print(f"  {batch}: {len(idx)} cells, {time.time()-t0:.1f}s")

print(f"Pre-CNEM: mean={pre_cnem.mean():.4f}")

# ============================================================
# 3. Integration + post-integration CNEM (GPU)
# ============================================================
print("\n" + "=" * 60)
print("Step 3: Integration + post-CNEM (GPU)")
print("=" * 60)

def run_post_cnem(embedding, label):
    t0 = time.time()
    indices = faiss_knn_gpu(embedding, K)
    indices_torch = torch.from_numpy(indices.astype(np.int64)).cuda(GPU_ID)
    post = compute_cnem_gpu(X_expr_torch, indices_torch)
    print(f"  {label}: kNN + CNEM in {time.time()-t0:.1f}s")
    return post

results = {}

# --- Harmony ---
print("\n--- Harmony ---")
t0 = time.time()
import harmonypy as hm
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'Dataset')
Z = np.array(ho.Z_corr)
if Z.shape[0] != adata.n_obs:
    Z = Z.T
results['Harmony'] = run_post_cnem(Z, 'Harmony')
print(f"  Total: {time.time()-t0:.1f}s")

# --- Scanorama ---
print("\n--- Scanorama ---")
t0 = time.time()
import scanorama
adata_s = adata.copy()
adata_s = adata_s[:, adata_s.var['highly_variable']].copy()
sc.pp.scale(adata_s, max_value=10)
sc.pp.pca(adata_s, n_comps=50)
batch_cats = list(adata_s.obs['Dataset'].cat.categories) if hasattr(adata_s.obs['Dataset'], 'cat') else sorted(adata_s.obs['Dataset'].unique())
adatas_list = [adata_s[adata_s.obs['Dataset'] == b].copy() for b in batch_cats]
scanorama.integrate_scanpy(adatas_list)
adata_sc = ad.concat(adatas_list)
adata_sc = adata_sc[adata.obs_names].copy()
results['Scanorama'] = run_post_cnem(adata_sc.obsm['X_scanorama'], 'Scanorama')
print(f"  Total: {time.time()-t0:.1f}s")

# --- BBKNN (graph-based, extract kNN then GPU CNEM) ---
print("\n--- BBKNN ---")
t0 = time.time()
import bbknn
adata_b = adata.copy()
adata_b = adata_b[:, adata_b.var['highly_variable']].copy()
sc.pp.scale(adata_b, max_value=10)
sc.pp.pca(adata_b, n_comps=50)
bbknn.bbknn(adata_b, batch_key='Dataset', n_pcs=50)
# Extract top-K neighbors from graph
conn = adata_b.obsp['connectivities']
bbknn_indices = np.zeros((adata.n_obs, K), dtype=np.int64)
for i in range(adata.n_obs):
    row = conn[i].toarray().flatten()
    top_k = np.argsort(row)[-K:]
    bbknn_indices[i] = top_k
bbknn_indices_torch = torch.from_numpy(bbknn_indices).cuda(GPU_ID)
results['BBKNN'] = compute_cnem_gpu(X_expr_torch, bbknn_indices_torch)
print(f"  Total: {time.time()-t0:.1f}s")

# ============================================================
# 4. Evaluation
# ============================================================
print("\n" + "=" * 60)
print("Step 4: Evaluation")
print("=" * 60)

celltypes = adata.obs['celltype_true'].values
is_rare = adata.obs['is_rare'].values

for method_name, post_cnem_vals in results.items():
    print(f"\n{'='*50}")
    print(f"  {method_name} — post_cnem")
    print(f"{'='*50}")

    scores = pd.DataFrame({
        'post_cnem': post_cnem_vals,
        'pre_cnem': pre_cnem,
        'celltype': celltypes,
        'is_rare': is_rare,
        'batch': adata.obs['Dataset'].values,
    })

    summary = scores.groupby('celltype')['post_cnem'].agg(['mean', 'median', 'std', 'count'])
    summary = summary.sort_values('mean', ascending=False)
    print(summary.round(4))

    # Wilcoxon
    rare_vals = scores[scores['is_rare'] == 1]['post_cnem']
    common_vals = scores[scores['is_rare'] == 0]['post_cnem']
    stat, pval = stats.mannwhitneyu(rare_vals, common_vals, alternative='greater')
    print(f"\nWilcoxon pDC+Mast > others: p={pval:.2e}")

    pdc_vals = scores[scores['celltype'] == 'pDC']['post_cnem']
    stat2, pval2 = stats.mannwhitneyu(pdc_vals, common_vals, alternative='greater')
    print(f"Wilcoxon pDC > others: p={pval2:.2e}")

    # ROC-AUC
    auc_rare = roc_auc_score(is_rare, post_cnem_vals)
    auc_pdc = roc_auc_score((celltypes == 'pDC').astype(int), post_cnem_vals)
    print(f"ROC-AUC (rare=pDC+Mast): {auc_rare:.3f}")
    print(f"ROC-AUC (pDC only): {auc_pdc:.3f}")

    # Top 5%
    threshold = scores['post_cnem'].quantile(0.95)
    top5 = scores[scores['post_cnem'] >= threshold]
    rare_in_top5 = top5['is_rare'].sum()
    rare_frac = is_rare.sum() / len(is_rare)
    enrichment = (rare_in_top5 / len(top5)) / rare_frac
    pdc_in_top5 = (top5['celltype'] == 'pDC').sum()
    mast_in_top5 = (top5['celltype'] == 'Mast').sum()
    print(f"\nTop 5%: {len(top5)} cells, rare={rare_in_top5} (OR={enrichment:.2f})")
    print(f"  pDC: {pdc_in_top5}/{(celltypes=='pDC').sum()}")
    print(f"  Mast: {mast_in_top5}/{(celltypes=='Mast').sum()}")
    print(f"  Types: {top5['celltype'].value_counts().to_dict()}")

    # Top 10%
    threshold_10 = scores['post_cnem'].quantile(0.90)
    top10 = scores[scores['post_cnem'] >= threshold_10]
    pdc_in_top10 = (top10['celltype'] == 'pDC').sum()
    rare_in_top10 = top10['is_rare'].sum()
    enrichment_10 = (rare_in_top10 / len(top10)) / rare_frac
    print(f"Top 10%: pDC={pdc_in_top10}/{(celltypes=='pDC').sum()}, OR={enrichment_10:.2f}")

    scores.to_csv(f'{OUT_DIR}/phase1_cnem_{method_name.lower()}.csv', index=False)

# Cross-method summary
print("\n" + "=" * 60)
print("Cross-method post_cnem summary")
print("=" * 60)
all_dfs = []
for method_name, vals in results.items():
    df = pd.DataFrame({'celltype': celltypes, 'post_cnem': vals, 'method': method_name})
    all_dfs.append(df)
combined = pd.concat(all_dfs)
pivot = combined.groupby(['celltype', 'method'])['post_cnem'].mean().unstack()
pivot['is_rare'] = pivot.index.isin(rare_types)
pivot = pivot.sort_values('is_rare', ascending=False)
print(pivot.round(4).to_string())

print(f"\nTotal time: {time.time()-t_total:.0f}s")
print("Phase 1 complete!")
