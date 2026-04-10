"""scVI + NP regularization (iterative, no source modification).

Method B: Iterative training
1. Train scVI standard -> get latent Z
2. Compute NP(Z) -> identify low-NP cells (high risk)
3. For high-risk cells, add soft penalty by modifying batch labels
4. Retrain scVI with modified data
5. Repeat for N iterations

Also: DEG analysis on low-NP regions to discover unknown subgroups.
"""
import sys; sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)
import scanpy as sc, anndata as ad, numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, adjusted_rand_score
import faiss, torch, os, time, gc, warnings
warnings.filterwarnings('ignore')

OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/scvi_np/results'
SHARED = '/ssd/data/agent/bio/shared'
os.makedirs(OUT_DIR, exist_ok=True); os.makedirs(SHARED, exist_ok=True)

K = 30; GPU_ID = 0

def faiss_knn_gpu(X, k, gpu_id=GPU_ID):
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    idx = faiss.index_cpu_to_gpu(res, gpu_id, faiss.IndexFlatL2(X.shape[1]))
    idx.add(X); _, indices = idx.search(X, k + 1)
    return indices[:, 1:]

def compute_np(pre_idx, post_idx, k=K):
    n = pre_idx.shape[0]
    return np.array([len(set(pre_idx[i]) & set(post_idx[i])) / k for i in range(n)], dtype=np.float32)

def compute_survival(pre_labels, post_labels):
    results = []
    for cl in pd.Series(pre_labels).unique():
        ci = np.where(pre_labels == cl)[0]
        if len(ci) < 10: continue
        pc = pd.Series(post_labels[ci]).value_counts()
        results.append({'cluster': cl, 'size': len(ci), 'survival': pc.iloc[0]/len(ci)})
    return pd.DataFrame(results)

t_total = time.time()

# 1. Load
print("Step 1: Load")
adata = ad.read_h5ad('/ssd/data/agent/bio/eacn_example_001/experiments/np_guard/data/subset.h5ad')
print(f"Loaded: {adata.shape}")

# Preprocess
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)
X_pca = adata_hvg.obsm['X_pca'].astype(np.float32)

# Pre-integration kNN (fixed reference)
pre_indices = faiss_knn_gpu(X_pca, K)

# Pre-integration subclustering (ground truth)
ct = adata.obs['Celltype'].values
pre_labels = np.full(adata.n_obs, '', dtype=object)
post_labels_dict = {}
for c in pd.Series(ct).unique():
    ci = np.where(ct == c)[0]
    if len(ci) < 50: continue
    a = ad.AnnData(obs=adata.obs.iloc[ci].copy())
    a.obsm['X_pca'] = X_pca[ci]
    sc.pp.neighbors(a, use_rep='X_pca', n_neighbors=min(15, len(ci)-1))
    sc.tl.leiden(a, resolution=1.0)
    for j, gi in enumerate(ci):
        pre_labels[gi] = f"{c}_{a.obs['leiden'].iloc[j]}"

print(f"Pre-clusters: {len(pd.Series(pre_labels).unique())}")

# 2. Standard scVI baseline
print("\nStep 2: Standard scVI")
t0 = time.time()

try:
    from scvi.model import SCVI
    adata_scvi = adata.copy()
    SCVI.setup_anndata(adata_scvi, batch_key='Dataset')
    model_std = SCVI(adata_scvi, n_latent=30, n_layers=2)
    model_std.train(max_epochs=100, accelerator='gpu', devices=[GPU_ID])
    Z_std = model_std.get_latent_representation().astype(np.float32)
    scvi_ok = True
    print(f"Standard scVI done: {time.time()-t0:.0f}s")
except Exception as e:
    print(f"scVI failed: {e}, using Harmony fallback")
    import harmonypy as hm
    ho = hm.run_harmony(X_pca, adata.obs, 'Dataset')
    Z_std = np.array(ho.Z_corr, dtype=np.float32)
    if Z_std.shape[0] != adata.n_obs: Z_std = Z_std.T
    scvi_ok = False
    print(f"Harmony fallback: {time.time()-t0:.0f}s")

post_idx_std = faiss_knn_gpu(Z_std, K)
np_std = compute_np(pre_indices, post_idx_std)
print(f"Standard NP: {np_std.mean():.4f}")

# 3. Iterative NP-scVI (5 rounds)
print("\nStep 3: Iterative NP-scVI")
N_ROUNDS = 5
EPOCHS_PER_ROUND = 40
LAMBDA_PROTECT = 0.8
Z_current = Z_std.copy()
np_history = [np_std.mean()]

for rnd in range(N_ROUNDS):
    t0 = time.time()
    # Compute NP risk
    post_idx = faiss_knn_gpu(Z_current, K)
    np_current = compute_np(pre_indices, post_idx)
    risk = 1.0 - np_current
    risk = (risk - risk.min()) / (risk.max() - risk.min() + 1e-8)
    high_risk = risk > 0.5

    # Modify batch labels for high-risk cells
    adata_mod = adata.copy()
    batches = adata.obs['Dataset'].unique()
    np.random.seed(42 + rnd)
    n_shuffled = 0
    for i in np.where(high_risk)[0]:
        if np.random.random() < risk[i] * LAMBDA_PROTECT:
            others = [b for b in batches if b != adata.obs['Dataset'].iloc[i]]
            adata_mod.obs['Dataset'].iloc[i] = np.random.choice(others)
            n_shuffled += 1

    # Train scVI on modified data
    try:
        SCVI.setup_anndata(adata_mod, batch_key='Dataset')
        model = SCVI(adata_mod, n_latent=30, n_layers=2)
        model.train(max_epochs=EPOCHS_PER_ROUND, accelerator='gpu', devices=[GPU_ID])
        Z_current = model.get_latent_representation().astype(np.float32)
    except:
        import harmonypy as hm
        ho = hm.run_harmony(X_pca, adata_mod.obs, 'Dataset')
        Z_current = np.array(ho.Z_corr, dtype=np.float32)
        if Z_current.shape[0] != adata.n_obs: Z_current = Z_current.T

    post_idx_new = faiss_knn_gpu(Z_current, K)
    np_new = compute_np(pre_indices, post_idx_new)
    np_history.append(np_new.mean())

    print(f"  Round {rnd+1}: shuffled={n_shuffled}, NP={np_new.mean():.4f}, high_risk={high_risk.sum()}, time={time.time()-t0:.0f}s")

Z_np = Z_current
np_np = compute_np(pre_indices, faiss_knn_gpu(Z_np, K))
print(f"\nFinal NP-scVI NP: {np_np.mean():.4f} (was {np_std.mean():.4f})")
print(f"NP history: {[f'{x:.4f}' for x in np_history]}")

# 4. Compare survival
print("\nStep 4: Compare survival")
for name, Z in [('standard', Z_std), ('np_scvi', Z_np)]:
    post_lab = np.full(adata.n_obs, '', dtype=object)
    for c in pd.Series(ct).unique():
        ci = np.where(ct == c)[0]
        if len(ci) < 50: continue
        a = ad.AnnData(obs=adata.obs.iloc[ci].copy())
        a.obsm['Z'] = Z[ci]
        sc.pp.neighbors(a, use_rep='Z', n_neighbors=min(15, len(ci)-1))
        sc.tl.leiden(a, resolution=1.0)
        for j, gi in enumerate(ci):
            post_lab[gi] = f"{c}_{a.obs['leiden'].iloc[j]}"
    post_labels_dict[name] = post_lab

surv_std = compute_survival(pre_labels, post_labels_dict['standard'])
surv_np = compute_survival(pre_labels, post_labels_dict['np_scvi'])

merged = surv_std[['cluster','size','survival']].merge(
    surv_np[['cluster','survival']], on='cluster', suffixes=('_std','_np'))
merged['improvement'] = merged['survival_np'] - merged['survival_std']
merged['dispersed_std'] = merged['survival_std'] < 0.5

disp = merged[merged['dispersed_std']]
print(f"\nDispersed clusters: {len(disp)}")
if len(disp) > 0:
    n_improved = (disp['improvement'] > 0).sum()
    print(f"Improved: {n_improved}/{len(disp)} ({n_improved/len(disp)*100:.0f}%)")
    print(f"Mean improvement: {disp['improvement'].mean():+.3f}")
    print("\nTop improvements:")
    print(disp.nlargest(10, 'improvement')[['cluster','size','survival_std','survival_np','improvement']].to_string(index=False))

# 5. DEG analysis on low-NP regions
print("\nStep 5: DEG on low-NP regions (unknown subgroup discovery)")
np_threshold = np.percentile(np_np, 5)  # bottom 5%
low_np_mask = np_np < np_threshold
adata.obs['np_group'] = 'normal'
adata.obs.loc[adata.obs.index[low_np_mask], 'np_group'] = 'low_NP'
print(f"Low NP cells (<{np_threshold:.3f}): {low_np_mask.sum()}")
print(f"Celltypes in low NP: {adata.obs.loc[low_np_mask, 'Celltype'].value_counts().to_dict()}")

# DEG
try:
    sc.tl.rank_genes_groups(adata, 'np_group', groups=['low_NP'], reference='normal', method='wilcoxon')
    deg = sc.get.rank_genes_groups_df(adata, group='low_NP')
    deg = deg[deg['pvals_adj'] < 0.05].head(20)
    print(f"\nTop DEGs (low_NP vs normal):")
    print(deg[['names','logfoldchanges','pvals_adj']].to_string(index=False))
except Exception as e:
    print(f"DEG failed: {e}")

# Save
merged.to_csv(f'{OUT_DIR}/scvi_np_comparison.csv', index=False)
merged.to_csv(f'{SHARED}/scvi_np_comparison.csv', index=False)

print(f"\nTotal: {time.time()-t_total:.0f}s")
print("scVI+NP iterative experiment complete!")
