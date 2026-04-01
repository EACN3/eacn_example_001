"""
NP-Guard: scVI 上的自适应保护策略实现伪代码
作者：机器学习智能体 (agent-mnez8qvx)
执行：计算生物学智能体 (agent-mneys6aw)

核心思想：在 scVI 训练中，用 NP（邻域保留率）作为 per-cell 信号
调节批次对齐损失权重。NP 低的区域降低整合力度。
"""

import numpy as np
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from scvi.model import SCVI
import torch

# ============================================================
# Step 1: 计算整合前的 kNN 图（只算一次）
# ============================================================

def compute_pre_knn(adata, n_neighbors=30, n_pcs=50):
    """在整合前的表达空间中计算 kNN 图"""
    # 使用 PCA 降维后的空间
    sc.pp.pca(adata, n_comps=n_pcs)
    nn = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean')
    nn.fit(adata.obsm['X_pca'])
    distances, indices = nn.kneighbors()
    return indices  # shape: (n_cells, n_neighbors)

# ============================================================
# Step 2: 计算 NP（邻域保留率）
# ============================================================

def compute_np(pre_knn_indices, post_embedding, n_neighbors=30):
    """
    计算每个细胞的邻域保留率
    pre_knn_indices: 整合前的 kNN 索引
    post_embedding: 整合后的嵌入（如 scVI latent space）
    """
    nn_post = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean')
    nn_post.fit(post_embedding)
    _, post_indices = nn_post.kneighbors()

    n_cells = pre_knn_indices.shape[0]
    np_scores = np.zeros(n_cells)

    for i in range(n_cells):
        pre_set = set(pre_knn_indices[i])
        post_set = set(post_indices[i])
        np_scores[i] = len(pre_set & post_set) / n_neighbors

    return np_scores

# ============================================================
# Step 3: 计算风险分数（相对 NP 变化）
# ============================================================

def compute_risk(np_pre, np_post, pre_knn_indices):
    """
    风险分数 = 1 - NP_post/NP_pre（局部平均版）
    采纳肿瘤生物学建议：用相对变化而非绝对阈值
    """
    n_cells = len(np_pre)
    risk = np.zeros(n_cells)

    for i in range(n_cells):
        # 局部平均 NP（用整合前邻域聚合降噪）
        neighbors = pre_knn_indices[i]
        np_pre_local = np.mean(np_pre[neighbors])
        np_post_local = np.mean(np_post[neighbors])

        if np_pre_local > 0:
            risk[i] = max(0, 1 - np_post_local / np_pre_local)
        else:
            risk[i] = 0  # 本来就没有同类邻居，不保护

    # 归一化到 [0, 1]
    if risk.max() > 0:
        risk = risk / risk.max()

    return risk

# ============================================================
# Step 4: NP-Guard 训练流程
# ============================================================

def train_scvi_with_np_guard(
    adata,
    pre_knn_indices,
    n_epochs=100,
    monitor_every=50,      # 每多少步监控 NP
    lambda_protect=0.8,    # 最大保护强度
    n_latent=30,
    n_neighbors=30,
):
    """
    NP-Guard 增强的 scVI 训练

    方案一（推荐）：两阶段训练
    - Phase A: 标准 scVI 训练（获取初始嵌入）
    - Phase B: 计算 NP 风险 → 带权重重新训练

    方案二（高级）：在线监控
    - 需要修改 scVI 源码，在训练循环中插入 NP 计算
    - 更复杂但理论上更优
    """

    # === 方案一：两阶段训练（推荐，不需改scVI源码）===

    # Phase A: 标准训练
    SCVI.setup_anndata(adata, batch_key="batch")
    model_standard = SCVI(adata, n_latent=n_latent)
    model_standard.train(max_epochs=n_epochs)

    # 获取标准整合嵌入
    latent_standard = model_standard.get_latent_representation()

    # 计算 NP
    # NP_pre: 整合前每个细胞的"同类邻居比例"（用表达空间 kNN）
    np_pre = compute_np(pre_knn_indices, adata.obsm['X_pca'], n_neighbors)
    # NP_post: 整合后的邻域保留率
    np_post = compute_np(pre_knn_indices, latent_standard, n_neighbors)

    # 计算风险分数
    risk = compute_risk(np_pre, np_post, pre_knn_indices)

    # Phase B: 带权重重新训练
    # 方法：用 risk 作为 per-cell 权重，降低高风险细胞的批次校正
    # scVI 支持通过 training_plan 自定义损失权重

    # 计算 per-cell 权重
    cell_weights = 1.0 - lambda_protect * risk  # [0.2, 1.0]
    adata.obs['np_guard_weight'] = cell_weights

    # 用加权损失重新训练
    # 注意：标准 scVI 不直接支持 per-cell batch loss 权重
    # 实现方式A：修改 scVI TrainingPlan（见下方）
    # 实现方式B（更简单）：对高风险细胞，在训练数据中人为减少其batch标签的信息量
    #   - 例如：对 risk > 0.5 的细胞，将其 batch 标签随机替换为其他批次
    #   - 这样模型无法从这些细胞学到batch信号，相当于降低了整合力度

    # 实现方式B（推荐，不需改源码）：
    adata_protected = adata.copy()
    high_risk_mask = risk > 0.5
    n_high_risk = high_risk_mask.sum()
    print(f"NP-Guard: {n_high_risk} cells ({n_high_risk/len(risk)*100:.1f}%) marked as high-risk")

    # 对高风险细胞，按 risk 概率打乱 batch 标签
    for i in np.where(high_risk_mask)[0]:
        if np.random.random() < risk[i]:
            # 随机分配一个其他批次的标签
            other_batches = [b for b in adata.obs['batch'].unique()
                          if b != adata.obs['batch'].iloc[i]]
            adata_protected.obs['batch'].iloc[i] = np.random.choice(other_batches)

    SCVI.setup_anndata(adata_protected, batch_key="batch")
    model_protected = SCVI(adata_protected, n_latent=n_latent)
    model_protected.train(max_epochs=n_epochs)

    latent_protected = model_protected.get_latent_representation()

    return {
        'latent_standard': latent_standard,
        'latent_protected': latent_protected,
        'risk': risk,
        'np_pre': np_pre,
        'np_post_standard': np_post,
        'np_post_protected': compute_np(pre_knn_indices, latent_protected, n_neighbors),
        'cell_weights': cell_weights,
    }

# ============================================================
# Step 5: 验证保护效果
# ============================================================

def validate_protection(adata, results, pre_cluster_key='pre_cluster'):
    """
    对比标准整合 vs NP-Guard 整合：
    1. 被打散亚群的 survival rate 是否提升
    2. 主流类型的批次校正是否保持
    """
    from sklearn.metrics import adjusted_rand_score

    # 在两种嵌入上做 Leiden 聚类
    for name, latent in [('standard', results['latent_standard']),
                          ('protected', results['latent_protected'])]:
        adata.obsm[f'X_{name}'] = latent
        sc.pp.neighbors(adata, use_rep=f'X_{name}', key_added=f'{name}_nn')
        sc.tl.leiden(adata, neighbors_key=f'{name}_nn', key_added=f'{name}_leiden')

    # 对每个 pre-cluster 计算 survival rate
    pre_clusters = adata.obs[pre_cluster_key].unique()
    results_table = []

    for pc in pre_clusters:
        mask = adata.obs[pre_cluster_key] == pc
        n_cells = mask.sum()

        for method in ['standard', 'protected']:
            post_labels = adata.obs[f'{method}_leiden'][mask]
            # survival = 最大同聚类比例
            survival = post_labels.value_counts().iloc[0] / n_cells
            results_table.append({
                'pre_cluster': pc,
                'n_cells': n_cells,
                'method': method,
                'survival_rate': survival,
            })

    import pandas as pd
    df = pd.DataFrame(results_table)

    # 关键对比：被打散亚群（survival < 0.5）是否在 NP-Guard 下改善
    low_survival = df[df['method'] == 'standard'].query('survival_rate < 0.5')['pre_cluster']
    for pc in low_survival:
        standard = df[(df['pre_cluster']==pc) & (df['method']=='standard')]['survival_rate'].values[0]
        protected = df[(df['pre_cluster']==pc) & (df['method']=='protected')]['survival_rate'].values[0]
        print(f"Cluster {pc}: standard={standard:.3f} → protected={protected:.3f} "
              f"({'↑ IMPROVED' if protected > standard else '↓ worse'})")

    return df

# ============================================================
# 主流程
# ============================================================

if __name__ == '__main__':
    # 1. 加载数据
    adata = sc.read_h5ad('/ssd/data/agent/bio/eacn_example_001/experiments/phase1/data/immune_subset_phase1.h5ad')

    # 2. 预处理（如果尚未完成）
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch')

    # 3. 生物学安全过滤（采纳肿瘤生物学建议）
    # - Doublet 过滤（假设已完成）
    # - 排除 HBA/HBB/MALAT1 主导群
    # - 标记应激群（FOS/JUN/HSP）

    # 4. 计算整合前 kNN
    pre_knn = compute_pre_knn(adata, n_neighbors=30)

    # 5. 运行 NP-Guard
    results = train_scvi_with_np_guard(
        adata,
        pre_knn,
        n_epochs=100,
        lambda_protect=0.8,
    )

    # 6. 验证
    df = validate_protection(adata, results, pre_cluster_key='pre_leiden')

    # 7. 保存
    df.to_csv('np_guard_validation_results.csv', index=False)
    print("Done! Results saved to np_guard_validation_results.csv")
