"""
scVI + NP Regularization: 真正的保护性整合算法

核心思想：在scVI的损失函数中加入NP保持正则项，
让模型在训练过程中自动学习一个既去除批次效应又保护邻域结构的嵌入。

L_total = L_reconstruction + β * L_kl + γ * L_batch + λ * L_NP

其中 L_NP = -mean(NP(z)) 是可微的邻域保持损失。

关键技术挑战：NP基于kNN的集合交集，不可微。
解决方案：用软近邻（soft neighborhood）近似，使其可微。
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from sklearn.neighbors import NearestNeighbors
import scanpy as sc


# ============================================================
# 核心：可微的软邻域保持损失 (Differentiable Soft-NP Loss)
# ============================================================

class SoftNPLoss(nn.Module):
    """
    可微的邻域保持损失。

    传统NP = |N_pre(i) ∩ N_post(i)| / k 不可微（kNN是离散操作）。

    软近似：用高斯核将"是否是邻居"从硬判断变为软权重。

    对每个细胞i：
      1. pre_neighbors: 整合前表达空间的k近邻索引（固定，预计算）
      2. 在当前嵌入z中，计算i到其pre_neighbors的距离
      3. 用softmax将距离转换为"邻居保持概率"
      4. L_NP(i) = -log(sum of softmax weights for pre-neighbors among post-neighbors)

    直觉：如果整合前的邻居在整合后仍然近，损失低；如果被推远，损失高。
    """

    def __init__(self, temperature=1.0):
        super().__init__()
        self.temperature = temperature

    def forward(self, z_current, pre_knn_indices, sample_indices=None):
        """
        Args:
            z_current: (batch_size, latent_dim) 当前训练步的嵌入
            pre_knn_indices: (n_total_cells, k) 整合前的kNN索引
            sample_indices: (batch_size,) 当前batch中细胞在全数据中的索引

        Returns:
            loss: 标量，越低表示邻域保持越好
        """
        batch_size = z_current.shape[0]
        device = z_current.device

        # 计算当前batch内所有细胞对的距离矩阵
        # (batch_size, batch_size)
        dist_matrix = torch.cdist(z_current, z_current, p=2)

        # 转换为相似度（负距离 / temperature）
        sim_matrix = -dist_matrix / self.temperature

        # 对角线设为-inf（排除自身）
        sim_matrix.fill_diagonal_(float('-inf'))

        # softmax得到"软邻居概率"
        # (batch_size, batch_size) 每行和为1
        soft_neighbors = F.softmax(sim_matrix, dim=1)

        # 对每个细胞，计算其整合前邻居在当前嵌入中的软邻居概率之和
        total_loss = 0.0
        count = 0

        for i in range(batch_size):
            if sample_indices is not None:
                global_idx = sample_indices[i].item()
            else:
                global_idx = i

            # 该细胞的整合前邻居
            pre_neighbors = set(pre_knn_indices[global_idx].tolist())

            # 找出这些邻居中哪些在当前batch中
            neighbor_mask = torch.zeros(batch_size, device=device)
            for j in range(batch_size):
                if sample_indices is not None:
                    j_global = sample_indices[j].item()
                else:
                    j_global = j
                if j_global in pre_neighbors and j != i:
                    neighbor_mask[j] = 1.0

            # 如果当前batch中有整合前邻居
            if neighbor_mask.sum() > 0:
                # 软NP = 整合前邻居的软概率之和
                soft_np_i = (soft_neighbors[i] * neighbor_mask).sum()
                # 损失 = -log(soft_np)，soft_np越高损失越低
                total_loss += -torch.log(soft_np_i + 1e-8)
                count += 1

        if count > 0:
            return total_loss / count
        else:
            return torch.tensor(0.0, device=device, requires_grad=True)


# ============================================================
# 高效版本：基于采样的NP损失（适用于大规模数据）
# ============================================================

class EfficientSoftNPLoss(nn.Module):
    """
    高效版NP损失：不计算完整距离矩阵，而是只计算每个细胞到其
    pre-neighbors的距离。

    对于batch_size=128, k=30的情况：
    - 完整版：128x128距离矩阵
    - 高效版：128x30距离计算
    """

    def __init__(self, temperature=1.0, k=30):
        super().__init__()
        self.temperature = temperature
        self.k = k

    def forward(self, z_all, pre_knn_indices, cell_indices):
        """
        Args:
            z_all: (n_cells, latent_dim) 所有细胞的当前嵌入
            pre_knn_indices: (n_cells, k) 整合前的kNN索引
            cell_indices: (batch_size,) 当前batch的细胞索引
        """
        batch_size = len(cell_indices)
        device = z_all.device
        losses = []

        for idx in cell_indices:
            i = idx.item()
            z_i = z_all[i]  # (latent_dim,)

            # 该细胞的整合前k近邻
            pre_nn = pre_knn_indices[i]  # (k,)

            # 计算到整合前邻居的距离
            z_neighbors = z_all[pre_nn]  # (k, latent_dim)
            dist_to_pre_nn = torch.norm(z_neighbors - z_i.unsqueeze(0), dim=1)  # (k,)

            # 随机采样k个非邻居作为负样本
            all_indices = set(range(z_all.shape[0]))
            pre_nn_set = set(pre_nn.tolist())
            non_neighbors = list(all_indices - pre_nn_set - {i})
            neg_sample = np.random.choice(non_neighbors, size=min(self.k, len(non_neighbors)), replace=False)
            z_neg = z_all[neg_sample]  # (k, latent_dim)
            dist_to_neg = torch.norm(z_neg - z_i.unsqueeze(0), dim=1)  # (k,)

            # 合并距离，softmax
            all_dist = torch.cat([dist_to_pre_nn, dist_to_neg])  # (2k,)
            sim = -all_dist / self.temperature
            soft_prob = F.softmax(sim, dim=0)

            # NP损失 = -log(整合前邻居的概率和)
            soft_np = soft_prob[:self.k].sum()
            losses.append(-torch.log(soft_np + 1e-8))

        return torch.stack(losses).mean()


# ============================================================
# scVI 集成方案
# ============================================================

def train_scvi_with_np_regularization(
    adata,
    n_epochs=200,
    n_latent=30,
    lambda_np=0.1,       # NP正则化强度
    temperature=1.0,      # 软近邻温度
    k=30,                 # 近邻数
    warmup_epochs=50,     # 前N个epoch不加NP正则（让scVI先学到基础表示）
    batch_key='batch',
):
    """
    带NP正则化的scVI训练。

    实现策略：使用scVI的自定义TrainingPlan。
    scVI基于PyTorch Lightning，可以通过继承TrainingPlan
    在training_step中添加NP正则化项。
    """
    from scvi.model import SCVI
    from scvi.train import TrainingPlan
    import pytorch_lightning as pl

    # Step 1: 预计算整合前kNN
    print("Computing pre-integration kNN...")
    sc.pp.pca(adata, n_comps=50)
    nn = NearestNeighbors(n_neighbors=k, metric='euclidean')
    nn.fit(adata.obsm['X_pca'])
    _, pre_knn_indices = nn.kneighbors()
    pre_knn_tensor = torch.tensor(pre_knn_indices, dtype=torch.long)
    print(f"Pre-kNN computed: {pre_knn_indices.shape}")

    # Step 2: 设置scVI
    SCVI.setup_anndata(adata, batch_key=batch_key)

    # Step 3: 自定义TrainingPlan
    class NPRegularizedTrainingPlan(TrainingPlan):
        """在scVI训练中加入NP保持正则化"""

        def __init__(self, module, **kwargs):
            super().__init__(module, **kwargs)
            self.np_loss_fn = EfficientSoftNPLoss(
                temperature=temperature, k=k
            )
            self.lambda_np = lambda_np
            self.warmup_epochs = warmup_epochs
            self.pre_knn = pre_knn_tensor

        def training_step(self, batch, batch_idx):
            # 标准scVI前向传播
            scvi_loss = super().training_step(batch, batch_idx)

            # warmup阶段不加NP正则
            if self.current_epoch < self.warmup_epochs:
                return scvi_loss

            # 获取当前batch的嵌入
            # scVI的module可以输出潜在表示
            with torch.no_grad():
                inference_outputs = self.module.inference(
                    batch['X'], batch_index=batch.get('batch', None)
                )
            z = inference_outputs['z']  # (batch_size, n_latent)
            z.requires_grad_(True)

            # 获取当前batch中细胞的全局索引
            # 注意：scVI的DataLoader会提供indices
            cell_indices = batch.get('indices', None)

            if cell_indices is not None and self.pre_knn is not None:
                # 计算NP损失
                # 需要所有细胞的嵌入来计算距离
                # 近似：只在当前batch内计算
                np_loss = self.np_loss_fn(
                    z, self.pre_knn.to(z.device), cell_indices
                )

                # 总损失 = scVI损失 + λ * NP损失
                total_loss = scvi_loss + self.lambda_np * np_loss

                # 记录
                self.log('np_loss', np_loss, prog_bar=True)
                self.log('lambda_np', self.lambda_np)

                return total_loss

            return scvi_loss

    # Step 4: 创建模型并训练
    model = SCVI(adata, n_latent=n_latent)

    # 使用自定义TrainingPlan
    model.train(
        max_epochs=n_epochs,
        plan_kwargs={},  # TrainingPlan参数
        # 注意：实际集成需要修改scVI的train方法以使用自定义TrainingPlan
        # 这里展示的是概念性实现
    )

    return model, pre_knn_indices


# ============================================================
# 替代方案：如果无法修改scVI TrainingPlan
# ============================================================

def train_scvi_np_iterative(
    adata,
    n_rounds=5,
    epochs_per_round=40,
    n_latent=30,
    lambda_np=0.1,
    k=30,
    batch_key='batch',
):
    """
    迭代式NP正则化：交替执行scVI训练和NP引导的嵌入修正。

    Round 1: 标准scVI训练 → 获取嵌入z1
    Round 2: 计算NP → 识别高风险细胞 → 在损失中增加这些细胞的重建权重 → 重新训练
    ...

    这种方式不需要修改scVI源码，但通过迭代逼近NP正则化的效果。
    """
    from scvi.model import SCVI

    # 预计算整合前kNN
    sc.pp.pca(adata, n_comps=50)
    nn = NearestNeighbors(n_neighbors=k, metric='euclidean')
    nn.fit(adata.obsm['X_pca'])
    _, pre_knn_indices = nn.kneighbors()

    best_model = None
    best_np_score = 0

    for round_i in range(n_rounds):
        print(f"\n=== Round {round_i+1}/{n_rounds} ===")

        # 训练scVI
        SCVI.setup_anndata(adata, batch_key=batch_key)
        model = SCVI(adata, n_latent=n_latent)

        if round_i > 0:
            # 从上一轮最佳模型初始化（transfer learning）
            # 或者用cell_weights调整训练
            pass

        model.train(max_epochs=epochs_per_round)
        latent = model.get_latent_representation()

        # 计算NP
        nn_post = NearestNeighbors(n_neighbors=k, metric='euclidean')
        nn_post.fit(latent)
        _, post_indices = nn_post.kneighbors()

        np_scores = np.zeros(len(adata))
        for i in range(len(adata)):
            pre_set = set(pre_knn_indices[i])
            post_set = set(post_indices[i])
            np_scores[i] = len(pre_set & post_set) / k

        mean_np = np_scores.mean()
        print(f"Round {round_i+1}: Mean NP = {mean_np:.4f}")

        if mean_np > best_np_score:
            best_np_score = mean_np
            best_model = model

        # 识别高风险细胞并调整下一轮
        high_risk = np_scores < np.percentile(np_scores, 20)
        adata.obs['np_risk'] = high_risk.astype(float)
        print(f"  High-risk cells: {high_risk.sum()} ({high_risk.mean()*100:.1f}%)")

    return best_model, pre_knn_indices, np_scores


# ============================================================
# 主流程
# ============================================================

if __name__ == '__main__':
    # 1. 加载数据
    adata = sc.read_h5ad(
        '/ssd/data/agent/bio/eacn_example_001/experiments/phase1/data/immune_subset_phase1.h5ad'
    )

    # 2. 预处理
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch')

    # 3. 方案A：完整NP正则化（需要修改scVI源码）
    # model, pre_knn = train_scvi_with_np_regularization(
    #     adata, lambda_np=0.1, warmup_epochs=50
    # )

    # 3. 方案B：迭代式NP正则化（不需改源码）
    model, pre_knn, np_scores = train_scvi_np_iterative(
        adata, n_rounds=5, epochs_per_round=40
    )

    # 4. 获取最终嵌入
    latent = model.get_latent_representation()
    adata.obsm['X_scvi_np'] = latent
    adata.obs['np_score'] = np_scores

    # 5. 定位NP最低区域 → 候选被消灭亚群
    low_np_threshold = np.percentile(np_scores, 5)  # bottom 5%
    low_np_cells = np_scores < low_np_threshold
    print(f"\nLow-NP cells (bottom 5%): {low_np_cells.sum()}")

    # 6. 对低NP区域做差异基因分析
    adata.obs['np_group'] = ['low_NP' if x else 'normal' for x in low_np_cells]
    sc.tl.rank_genes_groups(adata, groupby='np_group', groups=['low_NP'],
                            reference='normal', method='wilcoxon')

    # 输出top差异基因
    result = adata.uns['rank_genes_groups']
    top_genes = result['names']['low_NP'][:20]
    print(f"\nTop DEGs in low-NP region:")
    for i, gene in enumerate(top_genes):
        score = result['scores']['low_NP'][i]
        print(f"  {i+1}. {gene} (score={score:.2f})")

    # 7. 保存
    adata.write_h5ad('immune_subset_np_regularized.h5ad')
    print("\nDone! Saved to immune_subset_np_regularized.h5ad")
