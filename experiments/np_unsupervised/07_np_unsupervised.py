"""Fully unsupervised NP: use Leiden clusters instead of celltype labels.
Also: cross-method validation (Harmony + Scanorama) and GITR+ Treg verification.
GPU-accelerated.
"""
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, adjusted_rand_score
import faiss
import warnings, os, time
warnings.filterwarnings('ignore')

OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/np_unsupervised/results'
os.makedirs(OUT_DIR, exist_ok=True)

K = 30
GPU_ID = 0

def faiss_knn_gpu(X, k, gpu_id=GPU_ID):
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    index = faiss.index_cpu_to_gpu(res, gpu_id, faiss.IndexFlatL2(X.shape[1]))
    index.add(X)
    _, indices = index.search(X, k + 1)
    return indices[:, 1:]

def compute_stratified_np(X_pre, X_post, strat_labels, k=K):
    """Compute NP stratified by group labels (celltype or Leiden)."""
    n = len(strat_labels)
    np_arr = np.zeros(n, dtype=np.float32)
    groups = pd.Series(strat_labels).unique()
    for g in groups:
        g_mask = strat_labels == g
        g_idx = np.where(g_mask)[0]
        if len(g_idx) < k + 1:
            continue
        k_use = min(k, len(g_idx) - 1)
        pre_idx = faiss_knn_gpu(X_pre[g_idx], k_use)
        post_idx = faiss_knn_gpu(X_post[g_idx], k_use)
        for i in range(len(g_idx)):
            pre_set = set(pre_idx[i])
            post_set = set(post_idx[i])
            np_arr[g_idx[i]] = len(pre_set & post_set) / k_use
    return np_arr

def compute_subcluster_survival(pre_labels, post_labels, strat_labels):
    """Compute per-subcluster survival within each stratum."""
    results = []
    for g in pd.Series(strat_labels).unique():
        g_mask = strat_labels == g
        g_idx = np.where(g_mask)[0]
        pre_sub = pre_labels[g_idx]
        post_sub = post_labels[g_idx]
        for cl_id in pd.Series(pre_sub).unique():
            cl_mask = pre_sub == cl_id
            cl_local_idx = np.where(cl_mask)[0]
            if len(cl_local_idx) < 10:
                continue
            post_assign = post_sub[cl_local_idx]
            post_counts = pd.Series(post_assign).value_counts()
            survival = post_counts.iloc[0] / len(cl_local_idx)
            results.append({
                'stratum': g, 'pre_cluster': cl_id,
                'size': len(cl_local_idx), 'survival': survival,
                'mean_NP': 0.0,  # filled later
            })
    return pd.DataFrame(results)

# ============================================================
# 1. Load and preprocess
# ============================================================
print("=" * 60)
print("Step 1: Load and preprocess")
t_total = time.time()

adata = ad.read_h5ad('experiments/phase1/data/immune_subset_phase1.h5ad')
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)
X_pca = adata_hvg.obsm['X_pca'].copy()

# ============================================================
# 2. Unsupervised stratification: coarse Leiden
# ============================================================
print("\nStep 2: Unsupervised stratification (coarse Leiden)")
sc.pp.neighbors(adata_hvg, n_pcs=50)
sc.tl.leiden(adata_hvg, resolution=0.3, key_added='coarse_leiden')
n_coarse = adata_hvg.obs['coarse_leiden'].nunique()
print(f"Coarse clusters (res=0.3): {n_coarse}")

# Check how coarse Leiden aligns with celltype
ct_vs_leiden = pd.crosstab(adata.obs['Celltype'], adata_hvg.obs['coarse_leiden'])
print("\nCelltype vs coarse Leiden (top per cluster):")
for col in ct_vs_leiden.columns:
    top_ct = ct_vs_leiden[col].idxmax()
    purity = ct_vs_leiden[col].max() / ct_vs_leiden[col].sum()
    print(f"  Leiden {col}: {top_ct} ({purity:.0%}), n={ct_vs_leiden[col].sum()}")

# ============================================================
# 3. Integration methods
# ============================================================
print("\nStep 3: Integration")

import harmonypy as hm
ho = hm.run_harmony(X_pca, adata.obs, 'Dataset')
Z_harmony = np.array(ho.Z_corr)
if Z_harmony.shape[0] != adata.n_obs:
    Z_harmony = Z_harmony.T

import scanorama
adata_s = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_s, max_value=10)
sc.pp.pca(adata_s, n_comps=50)
batch_cats = sorted(adata_s.obs['Dataset'].unique())
adatas_list = [adata_s[adata_s.obs['Dataset'] == b].copy() for b in batch_cats]
scanorama.integrate_scanpy(adatas_list)
adata_sc = ad.concat(adatas_list)[adata.obs_names].copy()
Z_scanorama = adata_sc.obsm['X_scanorama']

methods = {'Harmony': Z_harmony, 'Scanorama': Z_scanorama}

# ============================================================
# 4. Compute NP: supervised (celltype) vs unsupervised (Leiden)
# ============================================================
print("\nStep 4: NP computation")

coarse_labels = adata_hvg.obs['coarse_leiden'].values
ct_labels = adata.obs['Celltype'].values

for method_name, Z in methods.items():
    print(f"\n--- {method_name} ---")

    # Supervised NP (celltype)
    np_sup = compute_stratified_np(X_pca, Z, ct_labels)
    # Unsupervised NP (Leiden)
    np_unsup = compute_stratified_np(X_pca, Z, coarse_labels)
    # Global NP (no stratification)
    pre_global = faiss_knn_gpu(X_pca, K)
    post_global = faiss_knn_gpu(Z, K)
    np_global = np.zeros(adata.n_obs, dtype=np.float32)
    for i in range(adata.n_obs):
        np_global[i] = len(set(pre_global[i]) & set(post_global[i])) / K

    # Subcluster survival (using celltype stratification as ground truth)
    # Pre/post subclustering within each celltype
    pre_sub = np.full(adata.n_obs, '', dtype=object)
    post_sub = np.full(adata.n_obs, '', dtype=object)
    for ct in pd.Series(ct_labels).unique():
        ct_mask = ct_labels == ct
        ct_idx = np.where(ct_mask)[0]
        if len(ct_idx) < 50:
            continue
        adata_ct = ad.AnnData(obs=adata.obs.iloc[ct_idx].copy())
        adata_ct.obsm['X_pca'] = X_pca[ct_idx]
        adata_ct.obsm['X_int'] = Z[ct_idx]
        k_use = min(15, len(ct_idx) - 1)
        sc.pp.neighbors(adata_ct, use_rep='X_pca', n_neighbors=k_use)
        sc.tl.leiden(adata_ct, resolution=1.0, key_added='pre')
        sc.pp.neighbors(adata_ct, use_rep='X_int', n_neighbors=k_use)
        sc.tl.leiden(adata_ct, resolution=1.0, key_added='post')
        for j, gi in enumerate(ct_idx):
            pre_sub[gi] = f"{ct}_{adata_ct.obs['pre'].iloc[j]}"
            post_sub[gi] = f"{ct}_{adata_ct.obs['post'].iloc[j]}"

    # Build survival DataFrame
    surv_list = []
    for cl_id in pd.Series(pre_sub).unique():
        if cl_id == '':
            continue
        cl_mask = pre_sub == cl_id
        cl_idx = np.where(cl_mask)[0]
        if len(cl_idx) < 10:
            continue
        post_assign = post_sub[cl_idx]
        post_counts = pd.Series(post_assign).value_counts()
        survival = post_counts.iloc[0] / len(cl_idx)
        surv_list.append({
            'cluster': cl_id, 'size': len(cl_idx), 'survival': survival,
            'np_sup': np_sup[cl_idx].mean(),
            'np_unsup': np_unsup[cl_idx].mean(),
            'np_global': np_global[cl_idx].mean(),
        })

    surv_df = pd.DataFrame(surv_list)
    surv_df['dispersed'] = (surv_df['survival'] < 0.5).astype(int)

    # ROC-AUC for each NP variant
    for np_col in ['np_sup', 'np_unsup', 'np_global']:
        if surv_df['dispersed'].sum() > 0 and surv_df['dispersed'].sum() < len(surv_df):
            auc = roc_auc_score(surv_df['dispersed'], -surv_df[np_col])
            corr, p = stats.spearmanr(surv_df[np_col], surv_df['survival'])
            print(f"  {np_col}: AUC={auc:.3f}, Spearman r={corr:.3f} p={p:.2e}")

    surv_df.to_csv(f'{OUT_DIR}/np_comparison_{method_name.lower()}.csv', index=False)

# ============================================================
# 5. GITR+ Treg verification
# ============================================================
print("\n" + "=" * 60)
print("Step 5: GITR+ Treg (TNFRSF18) verification")

if 'TNFRSF18' in adata.var_names:
    gitr_expr = adata[:, 'TNFRSF18'].X
    if hasattr(gitr_expr, 'toarray'):
        gitr_expr = gitr_expr.toarray()
    gitr_expr = gitr_expr.flatten()

    foxp3_expr = adata[:, 'FOXP3'].X
    if hasattr(foxp3_expr, 'toarray'):
        foxp3_expr = foxp3_expr.toarray()
    foxp3_expr = foxp3_expr.flatten()

    # GITR+ Treg = FOXP3+ & TNFRSF18+
    gitr_treg = (foxp3_expr > 0) & (gitr_expr > 0) & (ct_labels == 'T cell')
    print(f"GITR+ Treg (FOXP3+TNFRSF18+): {gitr_treg.sum()} cells")

    if gitr_treg.sum() > 0:
        # NP for GITR+ Treg vs other T cells
        t_mask = ct_labels == 'T cell'
        gitr_np = np_sup[gitr_treg]  # using Harmony supervised NP
        other_t_np = np_sup[t_mask & ~gitr_treg]
        stat, p = stats.mannwhitneyu(gitr_np, other_t_np, alternative='less')
        print(f"GITR+ Treg NP: mean={gitr_np.mean():.4f}")
        print(f"Other T cell NP: mean={other_t_np.mean():.4f}")
        print(f"Wilcoxon GITR+ < others: p={p:.2e}")

        # Dataset distribution
        print(f"Datasets: {adata.obs.loc[gitr_treg, 'Dataset'].value_counts().to_dict()}")
else:
    print("TNFRSF18 not found in var_names")

print(f"\nTotal time: {time.time()-t_total:.0f}s")
print("Unsupervised NP validation complete!")
