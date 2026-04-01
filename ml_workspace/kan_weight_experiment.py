"""
KAN-Weight实验：可解释的自适应整合权重学习
在105k免疫子集上验证KAN vs 线性BD加权

用法：python kan_weight_experiment.py
依赖：pip install pykan scanpy harmonypy faiss-gpu matplotlib
"""

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from pathlib import Path
import os
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# 配置
# ============================================================
DATA_PATH = '/ssd/data/agent/bio/eacn_example_001/experiments/phase1/data/immune_subset_phase1.h5ad'
OUTPUT_DIR = Path('/ssd/data/agent/bio/shared/kan_results/')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
K = 30  # kNN k值
N_PCS = 50  # PCA维度

# ============================================================
# Step 1: 数据加载和预处理
# ============================================================
print("Loading data...")
adata = sc.read_h5ad(DATA_PATH)
print(f"Data: {adata.shape[0]} cells, {adata.shape[1]} genes")

# 预处理（如果尚未完成）
if 'X_pca' not in adata.obsm:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch')
    sc.pp.pca(adata, n_comps=N_PCS)

# ============================================================
# Step 2: 计算Harmony整合和BD
# ============================================================
print("Running Harmony...")
import harmonypy
ho = harmonypy.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
z_harmony = ho.Z_corr.T
z_pre = adata.obsm['X_pca']

# kNN和BD计算
print("Computing kNN and BD...")
nn = NearestNeighbors(n_neighbors=K, metric='euclidean')
nn.fit(z_pre)
_, pre_indices = nn.kneighbors()

batches = adata.obs['batch'].values
n_batches = len(set(batches))

bd = np.zeros(len(adata))
for i in range(len(adata)):
    neighbor_batches = set(batches[pre_indices[i]])
    bd[i] = len(neighbor_batches) / min(K, n_batches)

# ============================================================
# Step 3: 计算5维KAN输入特征
# ============================================================
print("Computing KAN input features...")

# 1. BD（已计算）
# 2. log密度
from sklearn.neighbors import KernelDensity
# 用kNN距离的倒数作为密度代理
nn_dist = NearestNeighbors(n_neighbors=K)
nn_dist.fit(z_pre)
distances, _ = nn_dist.kneighbors()
log_density = -np.log(distances[:, -1] + 1e-8)  # 第k个邻居距离的负log

# 3. 类型纯度（粗聚类后邻域中最大类型占比）
sc.pp.neighbors(adata, use_rep='X_pca')
sc.tl.leiden(adata, resolution=0.5, key_added='coarse_cluster')
clusters = adata.obs['coarse_cluster'].values
type_purity = np.zeros(len(adata))
for i in range(len(adata)):
    neighbor_clusters = clusters[pre_indices[i]]
    _, counts = np.unique(neighbor_clusters, return_counts=True)
    type_purity[i] = counts.max() / K

# 4. MNN对数量（简化：跨批次邻居数）
n_mnn = np.zeros(len(adata))
for i in range(len(adata)):
    my_batch = batches[i]
    n_mnn[i] = sum(1 for j in pre_indices[i] if batches[j] != my_batch) / K

# 5. 表达特异性（离全局均值的距离）
global_mean = z_pre.mean(axis=0)
expr_specificity = np.linalg.norm(z_pre - global_mean, axis=1)
expr_specificity = expr_specificity / expr_specificity.max()  # 归一化到[0,1]

# 组装特征矩阵
features = np.column_stack([bd, log_density, type_purity, n_mnn, expr_specificity])
feature_names = ['BD', 'log_density', 'type_purity', 'cross_batch_ratio', 'expr_specificity']

# 归一化到[0,1]
for j in range(features.shape[1]):
    fmin, fmax = features[:, j].min(), features[:, j].max()
    if fmax > fmin:
        features[:, j] = (features[:, j] - fmin) / (fmax - fmin)

print(f"Features shape: {features.shape}")

# ============================================================
# Step 4: 计算NP作为训练目标
# ============================================================
print("Computing NP for training target...")

# 线性BD加权的嵌入
z_linear = bd[:, None] * z_harmony + (1 - bd[:, None]) * z_pre

# 计算各方案的NP
def compute_np(z_post, pre_indices, k=K):
    nn_post = NearestNeighbors(n_neighbors=k)
    nn_post.fit(z_post)
    _, post_indices = nn_post.kneighbors()
    np_scores = np.zeros(len(z_post))
    for i in range(len(z_post)):
        pre_set = set(pre_indices[i])
        post_set = set(post_indices[i])
        np_scores[i] = len(pre_set & post_set) / k
    return np_scores

np_harmony = compute_np(z_harmony, pre_indices)
np_linear = compute_np(z_linear, pre_indices)

print(f"NP Harmony: {np_harmony.mean():.4f}")
print(f"NP Linear BD: {np_linear.mean():.4f}")

# ============================================================
# Step 5: KAN训练
# ============================================================
print("Training KAN...")

try:
    from kan import KAN
    import torch

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")

    # 准备数据
    X_train = torch.tensor(features, dtype=torch.float32).to(device)

    # 训练目标：每个细胞的"最优权重"
    # 通过grid search找：对每个细胞，什么权重w使得加权后的NP最高
    # 简化：用NP_harmony和NP_pre的比值作为代理
    np_pre_scores = compute_np(z_pre, pre_indices)

    # 最优权重近似：如果harmony NP高→权重应高，如果pre NP高→权重应低
    # w_optimal ≈ NP_harmony / (NP_harmony + NP_pre + eps)
    w_optimal = np_harmony / (np_harmony + np_pre_scores + 1e-8)
    w_optimal = np.clip(w_optimal, 0, 1)
    y_train = torch.tensor(w_optimal, dtype=torch.float32).unsqueeze(1).to(device)

    # 创建KAN [5→3→1]
    model = KAN(width=[5, 3, 1], grid=15, k=3, device=device)

    # 训练
    dataset = {'train_input': X_train, 'train_label': y_train,
               'test_input': X_train[:1000], 'test_label': y_train[:1000]}

    model.fit(dataset, opt='LBFGS', steps=100, lamb=0.01)

    # 获取KAN预测的权重
    with torch.no_grad():
        w_kan = model(X_train).cpu().numpy().flatten()
    w_kan = np.clip(w_kan, 0, 1)

    # KAN加权嵌入
    z_kan = w_kan[:, None] * z_harmony + (1 - w_kan[:, None]) * z_pre
    np_kan = compute_np(z_kan, pre_indices)

    print(f"NP KAN: {np_kan.mean():.4f}")

    kan_success = True

except Exception as e:
    print(f"KAN training failed: {e}")
    print("Falling back to simple MLP approximation...")

    # Fallback: 用简单的numpy多项式拟合
    from numpy.polynomial import polynomial as P

    # 对BD维度做多项式拟合
    w_optimal = np_harmony / (np_harmony + np_pre_scores + 1e-8)
    coeffs = np.polyfit(bd, w_optimal, deg=5)
    w_poly = np.clip(np.polyval(coeffs, bd), 0, 1)

    z_poly = w_poly[:, None] * z_harmony + (1 - w_poly[:, None]) * z_pre
    np_poly = compute_np(z_poly, pre_indices)

    print(f"NP Polynomial: {np_poly.mean():.4f}")
    kan_success = False

# ============================================================
# Step 6: 结果对比
# ============================================================
print("\n=== Results Comparison ===")
results = {
    'Method': ['Standard Harmony', 'Linear BD', 'KAN/Poly BD'],
    'Mean NP': [np_harmony.mean(), np_linear.mean(),
                np_kan.mean() if kan_success else np_poly.mean()],
}

# Batch ASW
from sklearn.metrics import silhouette_score
for name, z in [('harmony', z_harmony), ('linear', z_linear),
                ('kan', z_kan if kan_success else z_poly)]:
    try:
        asw = silhouette_score(z[:5000], adata.obs['batch'].values[:5000],
                               sample_size=5000)
        results.setdefault('Batch ASW', []).append(f"{asw:.4f}")
    except:
        results.setdefault('Batch ASW', []).append('N/A')

import pandas as pd
df_results = pd.DataFrame(results)
print(df_results.to_string(index=False))
df_results.to_csv(OUTPUT_DIR / 'kan_vs_linear_comparison.csv', index=False)

# ============================================================
# Step 7: 样条曲线可视化（Nature Figure级别）
# ============================================================
print("\nGenerating spline curve visualizations...")

fig, axes = plt.subplots(1, 5, figsize=(25, 5))
fig.suptitle('KAN-Learned Integration Weight Functions', fontsize=16, y=1.02)

for dim, (ax, name) in enumerate(zip(axes, feature_names)):
    # 扫描该维度，固定其他维度在中位数
    x_sweep = np.tile(np.median(features, axis=0), (200, 1))
    x_range = np.linspace(0, 1, 200)
    x_sweep[:, dim] = x_range

    if kan_success:
        with torch.no_grad():
            x_tensor = torch.tensor(x_sweep, dtype=torch.float32).to(device)
            w_sweep = model(x_tensor).cpu().numpy().flatten()
    else:
        # Fallback: 画BD vs optimal weight的散点+拟合
        if dim == 0:  # BD
            w_sweep = np.polyval(coeffs, x_range)
        else:
            # 对其他维度做简单拟合
            c = np.polyfit(features[:, dim], w_optimal, deg=3)
            w_sweep = np.polyval(c, x_range)

    w_sweep = np.clip(w_sweep, 0, 1)

    # 画散点（采样）
    idx = np.random.choice(len(features), 2000, replace=False)
    ax.scatter(features[idx, dim], w_optimal[idx], alpha=0.1, s=1, c='gray', label='data')
    # 画KAN学到的曲线
    ax.plot(x_range, w_sweep, 'r-', linewidth=2.5, label='KAN-learned')
    # 画线性基线
    if dim == 0:
        ax.plot(x_range, x_range, 'b--', linewidth=1.5, alpha=0.7, label='linear')

    ax.set_xlabel(name, fontsize=12)
    ax.set_ylabel('Integration weight' if dim == 0 else '', fontsize=12)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'kan_spline_curves.png', dpi=300, bbox_inches='tight')
plt.savefig(OUTPUT_DIR / 'kan_spline_curves.pdf', bbox_inches='tight')
print(f"Saved to {OUTPUT_DIR / 'kan_spline_curves.png'}")

# ============================================================
# Step 8: 符号回归尝试
# ============================================================
if kan_success:
    print("\nAttempting symbolic regression...")
    try:
        model.auto_symbolic()
        formula = model.symbolic_formula()
        print(f"Discovered formula: {formula}")
        with open(OUTPUT_DIR / 'kan_symbolic_formula.txt', 'w') as f:
            f.write(str(formula))
    except Exception as e:
        print(f"Symbolic regression failed: {e}")

print("\n=== KAN Experiment Complete ===")
print(f"Results saved to {OUTPUT_DIR}")
