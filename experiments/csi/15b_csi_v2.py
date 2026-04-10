"""CSI: Contrastive Selective Integration.
MNN-based contrastive learning — only cells with cross-batch MNN pairs get aligned.
Unique cells stay in place. GPU-accelerated.
"""
import sys; sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)
import scanpy as sc, anndata as ad, numpy as np, pandas as pd
import torch, torch.nn as nn, torch.optim as optim
from sklearn.metrics import silhouette_score
import faiss, time, os, warnings
warnings.filterwarnings('ignore')

SHARED = '/ssd/data/agent/bio/shared'
OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/csi/results'
os.makedirs(OUT_DIR, exist_ok=True); os.makedirs(SHARED, exist_ok=True)
K = 30; GPU_ID = 0; D_LATENT = 30

def faiss_knn_gpu(X, k, gpu_id=GPU_ID):
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    idx = faiss.index_cpu_to_gpu(res, gpu_id, faiss.IndexFlatL2(X.shape[1]))
    idx.add(X); _, indices = idx.search(X, k + 1)
    return indices[:, 1:]

def compute_mnn_pairs(X, batches, k=10):
    """Compute MNN pairs across all batch pairs."""
    unique_b = np.unique(batches)
    pairs = []
    for i, b1 in enumerate(unique_b):
        for b2 in unique_b[i+1:]:
            idx1 = np.where(batches == b1)[0]
            idx2 = np.where(batches == b2)[0]
            if len(idx1) < k or len(idx2) < k: continue
            # kNN from b1 to b2
            X1 = np.ascontiguousarray(X[idx1], dtype=np.float32)
            X2 = np.ascontiguousarray(X[idx2], dtype=np.float32)
            res = faiss.StandardGpuResources()
            idx_gpu = faiss.index_cpu_to_gpu(res, GPU_ID, faiss.IndexFlatL2(X2.shape[1]))
            idx_gpu.add(X2); _, nn12 = idx_gpu.search(X1, k)
            idx_gpu2 = faiss.index_cpu_to_gpu(res, GPU_ID, faiss.IndexFlatL2(X1.shape[1]))
            idx_gpu2.add(X1); _, nn21 = idx_gpu2.search(X2, k)
            # MNN: mutual nearest neighbors
            for a in range(len(idx1)):
                for nb in nn12[a]:
                    if a in nn21[nb]:
                        pairs.append((idx1[a], idx2[nb]))
    return pairs

class CSIEncoder(nn.Module):
    def __init__(self, in_dim, latent_dim=D_LATENT):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, 512), nn.GELU(), nn.BatchNorm1d(512),
            nn.Linear(512, 256), nn.GELU(), nn.BatchNorm1d(256),
            nn.Linear(256, latent_dim),
        )
    def forward(self, x):
        z = self.net(x)
        return nn.functional.normalize(z, dim=1)

def compute_np(pre_idx, post_idx, k=K):
    return np.array([len(set(pre_idx[i]) & set(post_idx[i])) / k for i in range(pre_idx.shape[0])], dtype=np.float32)

t_total = time.time()

# 1. Load
print("Step 1: Load")
adata = ad.read_h5ad('/ssd/data/agent/bio/eacn_example_001/experiments/np_guard/data/subset.h5ad')
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='Dataset')
adata_hvg = adata[:, adata.var['highly_variable']].copy()

# HVG expression (normalized, not scaled) as encoder input
X_hvg = adata_hvg.X
if hasattr(X_hvg, 'toarray'): X_hvg = X_hvg.toarray()
X_hvg = np.array(X_hvg, dtype=np.float32)

sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)
X_pca = adata_hvg.obsm['X_pca'].astype(np.float32)
batches = adata.obs['Dataset'].values
ct = adata.obs['Celltype'].values
n = adata.n_obs
print(f"Cells: {n}, HVG dim: {X_hvg.shape[1]}")

# 2. Pre-integration kNN (reference)
print("\nStep 2: Pre-integration kNN + MNN pairs")
t0 = time.time()
pre_indices = faiss_knn_gpu(X_pca, K)

# Compute MNN pairs in PCA space
mnn_pairs = compute_mnn_pairs(X_pca, batches, k=10)
print(f"MNN pairs: {len(mnn_pairs)}, time: {time.time()-t0:.0f}s")

# Intra-batch kNN for structure loss
intra_knn = {}
for b in np.unique(batches):
    bi = np.where(batches == b)[0]
    if len(bi) < K: continue
    knn_b = faiss_knn_gpu(X_pca[bi], min(K, len(bi)-1))
    for j, gi in enumerate(bi):
        intra_knn[gi] = bi[knn_b[j]]

# 3. Train CSI
print("\nStep 3: Train CSI")
t0 = time.time()

X_tensor = torch.from_numpy(X_hvg).cuda(GPU_ID)
model = CSIEncoder(X_hvg.shape[1], D_LATENT).cuda(GPU_ID)
optimizer = optim.Adam(model.parameters(), lr=1e-3)
tau = 0.1; lambda_struct = 1.0
mnn_arr = np.array(mnn_pairs)

# Initial embedding for structure reference
with torch.no_grad():
    Z_init = model(X_tensor).cpu().numpy()
    init_dists = {}
    for gi, nbrs in intra_knn.items():
        init_dists[gi] = np.linalg.norm(Z_init[gi] - Z_init[nbrs], axis=1)

N_EPOCHS = 200; BATCH_SIZE = 4096

for epoch in range(N_EPOCHS):
    model.train()
    # Sample MNN pairs for contrastive loss
    perm = np.random.permutation(len(mnn_arr))[:BATCH_SIZE]
    anchors = mnn_arr[perm, 0]
    positives = mnn_arr[perm, 1]

    z_all = model(X_tensor)
    z_a = z_all[anchors]
    z_p = z_all[positives]

    # InfoNCE: positive pair similarity vs negatives from same batch
    pos_sim = (z_a * z_p).sum(dim=1) / tau
    # Negatives: random cells from positive's batch
    neg_idx = np.random.choice(n, size=(len(anchors), 64))
    z_neg = z_all[neg_idx]
    neg_sim = torch.bmm(z_a.unsqueeze(1), z_neg.transpose(1, 2)).squeeze(1) / tau
    logits = torch.cat([pos_sim.unsqueeze(1), neg_sim], dim=1)
    labels = torch.zeros(len(anchors), dtype=torch.long, device=z_all.device)
    L_contrastive = nn.functional.cross_entropy(logits, labels)

    # Structure loss: sample intra-batch pairs
    struct_idx = np.random.choice(list(intra_knn.keys()), size=min(BATCH_SIZE, len(intra_knn)))
    L_struct = torch.tensor(0.0, device=z_all.device)
    for gi in struct_idx[:256]:  # subsample for speed
        nbrs = intra_knn[gi][:5]
        z_i = z_all[gi]
        z_n = z_all[nbrs]
        curr_dist = torch.norm(z_i - z_n, dim=1)
        ref_dist = torch.from_numpy(init_dists[gi][:5]).cuda(GPU_ID)
        margin = 0.01
        L_struct += torch.relu(curr_dist - ref_dist - margin).mean()
    L_struct /= 256

    loss = L_contrastive + lambda_struct * L_struct
    optimizer.zero_grad(); loss.backward(); optimizer.step()

    if (epoch + 1) % 10 == 0:
        print(f"  Epoch {epoch+1}: L_contra={L_contrastive.item():.3f}, L_struct={L_struct.item():.4f}")

# Get final embedding
model.eval()
with torch.no_grad():
    Z_csi = model(X_tensor).cpu().numpy().astype(np.float32)
print(f"CSI training: {time.time()-t0:.0f}s")

# 4. Standard Harmony baseline
print("\nStep 4: Harmony baseline")
t0 = time.time()
import harmonypy as hm
ho = hm.run_harmony(X_pca, adata.obs, 'Dataset')
Z_harm = np.array(ho.Z_corr, dtype=np.float32)
if Z_harm.shape[0] != n: Z_harm = Z_harm.T
print(f"Harmony: {time.time()-t0:.0f}s")

# 5. Compare
print("\nStep 5: Compare")
methods = {'Harmony': Z_harm, 'CSI': Z_csi}
for name, Z in methods.items():
    post_idx = faiss_knn_gpu(Z, K)
    np_s = compute_np(pre_indices, post_idx)

    pre_lab = np.full(n, '', dtype=object)
    post_lab = np.full(n, '', dtype=object)
    for c in pd.Series(ct).unique():
        ci = np.where(ct == c)[0]
        if len(ci) < 50: continue
        a = ad.AnnData(obs=adata.obs.iloc[ci].copy())
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
        surv.append({'cluster': cl, 'size': len(ci), 'survival': pc.iloc[0]/len(ci)})
    sdf = pd.DataFrame(surv)
    n_disp = (sdf['survival'] < 0.5).sum()
    np.random.seed(42)
    sample = np.random.choice(n, min(5000, n), replace=False)
    asw = silhouette_score(Z[sample], batches[sample])
    print(f"{name}: NP={np_s.mean():.4f}, dispersed={n_disp}/{len(sdf)}, surv={sdf['survival'].mean():.3f}, ASW={asw:.3f}")

print(f"\nTotal: {time.time()-t_total:.0f}s")
print("CSI experiment complete!")
