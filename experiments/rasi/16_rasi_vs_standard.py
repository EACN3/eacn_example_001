"""RASI vs Standard pipeline comparison."""
import sys; sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)
import scanpy as sc, anndata as ad, numpy as np, pandas as pd
from sklearn.metrics import silhouette_score
import faiss, time, os, warnings
warnings.filterwarnings('ignore')

SHARED = '/ssd/data/agent/bio/shared'
OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/rasi/results'
os.makedirs(OUT_DIR, exist_ok=True); os.makedirs(SHARED, exist_ok=True)
K = 30; GPU_ID = 0

def faiss_knn_gpu(X, k, gpu_id=GPU_ID):
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    idx = faiss.index_cpu_to_gpu(res, gpu_id, faiss.IndexFlatL2(X.shape[1]))
    idx.add(X); _, indices = idx.search(X, k + 1)
    return indices[:, 1:]

def compute_np(pre_idx, post_idx, k=K):
    return np.array([len(set(pre_idx[i]) & set(post_idx[i])) / k for i in range(pre_idx.shape[0])], dtype=np.float32)

def eval_pipeline(Z, X_pca, pre_indices, batches, ct, name, k=K):
    post_idx = faiss_knn_gpu(Z, k)
    np_s = compute_np(pre_indices, post_idx)
    pre_lab = np.full(len(ct), '', dtype=object)
    post_lab = np.full(len(ct), '', dtype=object)
    for c in pd.Series(ct).unique():
        ci = np.where(ct == c)[0]
        if len(ci) < 50: continue
        a = ad.AnnData(obs=pd.DataFrame({'ct': ct[ci]}))
        a.obsm['p'] = X_pca[ci]; a.obsm['z'] = Z[ci]
        sc.pp.neighbors(a, use_rep='p', n_neighbors=min(15, len(ci)-1))
        sc.tl.leiden(a, resolution=1.0, key_added='pre')
        sc.pp.neighbors(a, use_rep='z', n_neighbors=min(15, len(ci)-1))
        sc.tl.leiden(a, resolution=1.0, key_added='post')
        for j, gi in enumerate(ci):
            pre_lab[gi] = f"{c}_{a.obs['pre'].iloc[j]}"; post_lab[gi] = f"{c}_{a.obs['post'].iloc[j]}"
    surv = []
    for cl in pd.Series(pre_lab).unique():
        ci = np.where(pre_lab == cl)[0]
        if len(ci) < 10: continue
        pc = pd.Series(post_lab[ci]).value_counts()
        surv.append({'cluster': cl, 'survival': pc.iloc[0]/len(ci), 'size': len(ci)})
    sdf = pd.DataFrame(surv)
    np.random.seed(42)
    sample = np.random.choice(len(ct), min(5000, len(ct)), replace=False)
    asw = silhouette_score(Z[sample], batches[sample])
    n_disp = (sdf['survival'] < 0.5).sum()
    print(f"{name}: NP={np_s.mean():.4f}, dispersed={n_disp}/{len(sdf)}, surv={sdf['survival'].mean():.3f}, ASW={asw:.3f}")
    return sdf, np_s

t_total = time.time()

# Load
print("Load...")
adata = ad.read_h5ad('/ssd/data/agent/bio/eacn_example_001/experiments/np_guard/data/subset.h5ad')
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
batches = adata.obs['Dataset'].values
ct = adata.obs['Celltype'].values
B = len(np.unique(batches))

# ============================================================
# Pipeline A: Standard (HVG 2000 + PCA 50 + Harmony)
# ============================================================
print("\n=== Pipeline A: Standard ===")
t0 = time.time()
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')
hvg_std = set(adata.var_names[adata.var['highly_variable']])

adata_a = adata[:, list(hvg_std)].copy()
sc.pp.scale(adata_a, max_value=10)
sc.pp.pca(adata_a, n_comps=50)
X_pca_a = adata_a.obsm['X_pca'].astype(np.float32)

import harmonypy as hm
ho = hm.run_harmony(X_pca_a, adata.obs, 'Dataset')
Z_a = np.array(ho.Z_corr, dtype=np.float32)
if Z_a.shape[0] != adata.n_obs: Z_a = Z_a.T
t_a = time.time() - t0
print(f"Time: {t_a:.0f}s")

# Check marker coverage
markers_rare = ['CLEC4C','IL3RA','IRF7','CCR8','CFTR','ASCL3','GHRL','SOX10','AXL','SIGLEC6']
cov_a = sum(1 for m in markers_rare if m in hvg_std and m in adata.var_names)
print(f"Rare marker coverage: {cov_a}/{len(markers_rare)}")

pre_idx_a = faiss_knn_gpu(X_pca_a, K)
surv_a, np_a = eval_pipeline(Z_a, X_pca_a, pre_idx_a, batches, ct, "Standard")

# ============================================================
# Pipeline B: RASI (Rare-HVG + PCA 200 + BD-Harmony)
# ============================================================
print("\n=== Pipeline B: RASI ===")
t0 = time.time()

# Step 1: Rare-HVG
# Start with standard HVG, then add cluster-specific genes
sc.pp.neighbors(adata_a, n_pcs=50)
sc.tl.leiden(adata_a, resolution=0.5, key_added='coarse')
adata.obs['coarse'] = adata_a.obs['coarse']
sc.tl.rank_genes_groups(adata, 'coarse', groups=adata.obs['coarse'].cat.categories.tolist(),
                        reference='rest', method='wilcoxon', key_added='rare_hvg')
# Collect top 50 genes per cluster
rare_genes = set()
for cl in adata_a.obs['coarse'].unique():
    try:
        degs = sc.get.rank_genes_groups_df(adata, group=cl, key='rare_hvg')
        rare_genes.update(degs.head(50)['names'].tolist())
    except:
        pass

# Merge with standard HVG
hvg_rasi = hvg_std | (rare_genes & set(adata.var_names))
print(f"Rare-HVG: {len(hvg_std)} standard + {len(rare_genes & set(adata.var_names))} rare-specific = {len(hvg_rasi)} total")

# Check marker coverage
cov_b = sum(1 for m in markers_rare if m in hvg_rasi and m in adata.var_names)
print(f"Rare marker coverage: {cov_b}/{len(markers_rare)} (was {cov_a})")
for m in markers_rare:
    if m in adata.var_names:
        in_std = m in hvg_std
        in_rasi = m in hvg_rasi
        if in_rasi and not in_std:
            print(f"  RESCUED: {m}")

# Step 2: PCA 200
adata_b = adata[:, list(hvg_rasi)].copy()
sc.pp.scale(adata_b, max_value=10)
sc.pp.pca(adata_b, n_comps=min(200, len(hvg_rasi)-1))
X_pca_b = adata_b.obsm['X_pca'].astype(np.float32)

# Step 3: BD-Harmony
knn_b = faiss_knn_gpu(X_pca_b, K)
BD = np.array([len(set(batches[knn_b[i]])) / min(K, B) for i in range(adata.n_obs)], dtype=np.float32)
BD_smooth = np.array([BD[knn_b[i]].mean() for i in range(adata.n_obs)], dtype=np.float32)

ho_b = hm.run_harmony(X_pca_b, adata.obs, 'Dataset')
Z_harm_b = np.array(ho_b.Z_corr, dtype=np.float32)
if Z_harm_b.shape[0] != adata.n_obs: Z_harm_b = Z_harm_b.T
Z_b = np.array([BD_smooth[i]*Z_harm_b[i] + (1-BD_smooth[i])*X_pca_b[i,:Z_harm_b.shape[1]] for i in range(adata.n_obs)], dtype=np.float32)

t_b = time.time() - t0
print(f"Time: {t_b:.0f}s (ratio: {t_b/t_a:.1f}x)")

pre_idx_b = faiss_knn_gpu(X_pca_b, K)
surv_b, np_b = eval_pipeline(Z_b, X_pca_b, pre_idx_b, batches, ct, "RASI")

# Step 5: Low-NP discovery
print("\n=== Low-NP Discovery ===")
np_threshold = np.percentile(np_b, 5)
low_np_mask = np_b < np_threshold
adata.obs['np_group'] = 'normal'
adata.obs.loc[adata.obs.index[low_np_mask], 'np_group'] = 'low_NP'
print(f"Low NP cells: {low_np_mask.sum()}")
print(f"Celltypes: {adata.obs.loc[low_np_mask, 'Celltype'].value_counts().to_dict()}")

try:
    sc.tl.rank_genes_groups(adata, 'np_group', groups=['low_NP'], reference='normal', method='wilcoxon')
    deg = sc.get.rank_genes_groups_df(adata, group='low_NP')
    deg = deg[deg['pvals_adj'] < 0.05].head(15)
    print(f"\nTop DEGs in low-NP region:")
    print(deg[['names','logfoldchanges','pvals_adj']].to_string(index=False))
except Exception as e:
    print(f"DEG failed: {e}")

# Summary
print(f"\n{'='*60}")
print(f"=== RASI vs Standard Summary ===")
print(f"{'='*60}")
print(f"| Metric | Standard | RASI |")
print(f"|--------|----------|------|")
print(f"| HVG count | {len(hvg_std)} | {len(hvg_rasi)} |")
print(f"| Rare marker coverage | {cov_a}/10 | {cov_b}/10 |")
print(f"| Time | {t_a:.0f}s | {t_b:.0f}s ({t_b/t_a:.1f}x) |")

# Save
pd.DataFrame({'metric': ['hvg_count','rare_marker_cov','time_s'],
              'standard': [len(hvg_std), cov_a, t_a],
              'rasi': [len(hvg_rasi), cov_b, t_b]}).to_csv(f'{SHARED}/rasi_vs_standard.csv', index=False)

print(f"\nTotal: {time.time()-t_total:.0f}s")
print("RASI validation complete!")
