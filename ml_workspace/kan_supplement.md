# KAN补充方案：可解释的自适应整合权重学习

> 作者：机器学习智能体 (agent-mnez8qvx)
> 定位：RASI的补充/升级方案，不替代核心流程

---

## 动机

RASI当前用线性BD加权：`z = BD * z_harmony + (1-BD) * z_pre`。这假设BD与最优整合强度是线性关系。但实际可能：
- 存在一个BD临界点（如BD<0.2时应完全不整合，BD>0.5时可全力整合）
- 最优权重还依赖局部密度、细胞类型纯度等其他特征
- 非线性映射可能进一步提升亚群保护效果

KAN（Kolmogorov-Arnold Networks）的独特优势：学到的映射函数可以可视化为样条曲线，直接揭示"什么样的细胞应该被整合、整合多少"的最优策略。

## 架构设计

### KAN-Weight: 可解释的整合权重生成器

```
输入特征（per-cell, 5维）：
  x_i = [BD(i),                    # batch diversity [0,1]
         log_density(i),           # 局部密度的log值
         type_purity(i),           # kNN中最大类型占比 [0,1]
         n_mnn(i),                 # 跨批次MNN对数量
         expr_specificity(i)]      # 表达特异性（离平均profile的距离）

KAN网络：
  [5] → [3] → [1]    # 两层KAN，极简

  输入层→隐藏层：5×3 = 15个可学习样条函数
  隐藏层→输出层：3×1 = 3个可学习样条函数
  总参数：~180个（每个样条用10个控制点）

输出：
  w(i) = sigmoid(KAN(x_i)) ∈ [0,1]  # 整合权重

最终嵌入：
  z_RASI_KAN(i) = w(i) * z_harmony(i) + (1-w(i)) * z_pre_aligned(i)
```

### 为什么用KAN而不是MLP

| 特性 | MLP [5→3→1] | KAN [5→3→1] |
|------|-------------|-------------|
| 参数量 | ~25 | ~180（样条控制点） |
| 可解释性 | 黑盒 | **每条边是一个可视化的样条曲线** |
| 函数发现 | 不能 | **可以自动发现BD→权重的最优函数形式** |
| 精度（小数据） | 中 | **高（KAN在低维小数据上优于MLP）** |
| 推理速度 | 快 | 稍慢（但输入仅5维，可忽略） |

关键优势：训练后可以**提取每个输入特征的单变量影响曲线**——直接画出"BD vs 整合权重"的最优函数形状，这是Nature级的可解释性结果。

## 训练方案

### 训练信号：NP作为自监督目标

不需要标签。用NP作为训练目标——学习使整合后NP最大化的权重函数：

```python
# 训练循环
for epoch in range(100):
    # 前向：用KAN生成权重
    features = compute_cell_features(adata)  # (n_cells, 5)
    weights = kan_model(features)             # (n_cells, 1)

    # 加权整合
    z_kan = weights * z_harmony + (1 - weights) * z_pre_aligned

    # 计算NP（可微近似版）
    # 用soft kNN：对每个细胞，计算到pre-neighbors的平均距离
    np_loss = 0
    for i in sample(n_cells, batch_size=1024):
        pre_neighbors = pre_knn_indices[i]
        dist_to_pre_nn = torch.norm(z_kan[pre_neighbors] - z_kan[i], dim=1)
        np_loss += dist_to_pre_nn.mean()  # 越小=邻居越近=NP越高

    # 同时鼓励batch mixing（防止完全不整合）
    batch_loss = compute_batch_entropy(z_kan, batch_labels)

    # 总损失
    loss = np_loss - lambda_batch * batch_loss
    loss.backward()
    optimizer.step()
```

### 训练效率

- 输入：5维 × n_cells
- KAN：180参数，前向传播O(n)
- NP近似：采样1024个细胞，O(1024 × k)
- 总训练：~100 epoch × ~10秒/epoch = **~15分钟**
- 推理：所有细胞一次前向 = **<1秒**

### 可选：两阶段训练

1. **Stage 1**：先用RASI线性BD加权作为warm start
2. **Stage 2**：用KAN在RASI基础上微调权重

这样KAN只需学习"线性BD加权的残差"，收敛更快。

## 可解释性分析

训练后提取每个输入维度的样条曲线：

```python
# 提取KAN学到的每个输入特征的边际效应
for dim, name in enumerate(['BD', 'log_density', 'type_purity', 'n_mnn', 'expr_specificity']):
    # 固定其他维度在均值，扫描该维度
    x_sweep = torch.zeros(100, 5)
    x_sweep[:, :] = feature_means
    x_sweep[:, dim] = torch.linspace(feature_min[dim], feature_max[dim], 100)

    w_sweep = kan_model(x_sweep).detach()

    plt.plot(x_sweep[:, dim], w_sweep)
    plt.xlabel(name)
    plt.ylabel('Integration weight')
    plt.title(f'KAN-learned: {name} → weight')
    plt.savefig(f'kan_edge_{name}.png')
```

**预期发现**（可作为论文Figure）：
- BD曲线：可能是S形（存在临界BD阈值）而非线性
- density曲线：低密度→低权重（稀有细胞少整合）
- type_purity曲线：高纯度→高权重（同质邻域可安全整合）
- 这些曲线直接揭示"什么样的细胞应该被保护"的生物学规律

## 与RASI核心流程的关系

```
RASI标准版（当前）：
  z = BD * z_harmony + (1-BD) * z_pre
  → 简单、高效、可解释（线性）
  → 论文Main Method

RASI-KAN升级版（本方案）：
  z = KAN(BD, density, purity, mnn, specificity) * z_harmony + (1-w) * z_pre
  → 更精确、自适应、可解释（非线性样条）
  → 论文Extended Method / Supplementary
  → 额外Figure：KAN学到的最优权重函数

两者对比实验：
  → 如果KAN > 线性BD：说明最优整合策略是非线性的（新发现）
  → 如果KAN ≈ 线性BD：说明线性BD已足够好（简洁性的胜利）
  → 无论哪种结果都有论文价值
```

## 符号公式发现（KAN的独特能力）

KAN训练后可以尝试符号回归：

```python
# pykan的符号公式拟合
kan_model.auto_symbolic()
# 可能发现：w = sigmoid(3.2 * BD - 1.5) * (1 - exp(-density/0.3))
# 这给出一个人类可读的最优整合权重公式
```

如果成功，这个公式本身就是Nature级的贡献——它用一个简洁的数学表达式回答了"什么样的细胞应该被整合"。

---

## KAN应用点2：NP的可微近似（KAN-NP）

### 问题

NP = |N_pre ∩ N_post| / k 基于硬kNN和集合交集，不可微。之前尝试用softmax近似但效果不佳（CSI失败的一个原因）。

### KAN方案

用KAN学习一个NP的可微代理函数：给定细胞的嵌入特征，直接预测其NP值。

```
训练阶段：
  1. 对训练集计算真实NP（硬kNN）
  2. 提取每个细胞的特征向量：
     f_i = [z_i,                           # 当前嵌入 (30维)
            mean(z_j for j in N_pre(i)),   # pre-邻域中心 (30维)
            std(z_j for j in N_pre(i)),    # pre-邻域spread (30维)
            local_density(i)]               # 局部密度 (1维)
     → 总91维

  3. KAN网络：[91] → [20] → [5] → [1]
     训练目标：MSE(KAN(f_i), NP_true(i))

推理阶段：
  KAN(f_i) ≈ NP(i)，完全可微
  可直接作为整合损失的一部分反向传播
```

### 为什么KAN比MLP更适合

- NP是一个**分段平滑函数**（邻居重叠比例的连续近似）——KAN的样条基函数天然适合拟合分段平滑
- 91维输入中大部分信息冗余，KAN的自动特征选择（通过样条系数稀疏化）能自动找到NP最依赖的几个关键维度
- 训练后可以可视化：NP最依赖嵌入的哪些维度？→ 揭示"邻域破坏"的几何特征

### 应用场景

1. **作为scVI的可微正则化项**：L_total = L_scVI + λ * (1 - KAN_NP(z))
   - 之前的SoftNPLoss用softmax近似不够准确
   - KAN_NP直接从真实NP学习，近似精度更高
   - 训练分两阶段：先训练KAN_NP代理→再用KAN_NP作正则化训练scVI

2. **实时监控**：在整合训练过程中，用KAN_NP实时预测每个细胞的NP，无需重新计算kNN
   - 标准NP计算需要重建kNN图（O(n log n)），每步都做太慢
   - KAN_NP推理只需一次前向传播（O(n)），可每个batch都监控

3. **作为NP-Guard的核心引擎**：用KAN_NP替代硬NP计算，NP-Guard从"每50步监控"变为"每步监控"

### 精度预期

NP本质是一个关于嵌入位置的平滑函数（邻居在嵌入空间中近→NP高），KAN应该能达到R² > 0.9的近似精度。可在训练集上验证。

---

## KAN应用点3：稀有性感知的HVG选择（KAN-HVG）

### 问题

标准HVG选择基于全局均值-方差关系，稀有亚群的marker基因（如pDC的CLEC4C、epsilon的GHRL）因全局方差低被排除。当前RASI的方案（标准HVG + cluster-specific top genes）是手工规则，需要先做粗聚类。

### KAN方案

用KAN学习"一个基因对稀有亚群有多重要"的评分函数，替代方差-based选择。

```
核心思想：
  重要基因 = 在少数细胞中高表达但在多数细胞中低表达的基因
  这恰好是稀有亚群marker的定义

特征（per-gene, 6维）：
  g_j = [mean_expr(j),                # 全局平均表达
         var_expr(j),                  # 全局方差（标准HVG用的）
         zero_fraction(j),             # 零表达比例
         max_expr(j),                  # 最大表达值
         skewness(j),                  # 偏度（右偏=少数高表达）
         bimodality_coef(j)]           # 双峰系数

KAN网络：[6] → [4] → [1]
  输出：rare_importance(j) ∈ [0,1]

训练信号（自监督）：
  不需要标签！用以下代理目标：
  - 正样本：已知在粗聚类小群（<5%细胞）中差异表达的基因
  - 负样本：在所有粗聚类中均匀表达的housekeeping基因
  - 或更简单：skewness > 3 且 zero_fraction > 0.95 的基因作为正样本代理
```

### KAN vs 手工规则的优势

| 方案 | 做法 | 缺点 |
|------|------|------|
| 标准HVG | top 2000 by variance | 遗漏稀有marker |
| RASI当前 | HVG + cluster-specific top50 | 需要先聚类，依赖聚类质量 |
| **KAN-HVG** | 学习的评分函数 | 自动、连续、可解释 |

### 可解释性产出

训练后可视化6个输入特征的样条曲线：
- 预期发现：`skewness`和`bimodality_coef`对rare_importance贡献最大
- 可能发现一个"稀有基因特征签名"——高偏度+高零表达比+低全局方差的组合
- 这个签名本身就是对"什么基因定义稀有亚群"的新理解

### 实用输出

```python
# 训练后
gene_scores = kan_hvg(gene_features)  # (n_genes,)
rare_hvg = genes[gene_scores > threshold]  # 稀有性感知HVG
final_hvg = standard_hvg.union(rare_hvg)   # 合并
```

预期增加200-500个稀有相关基因，总HVG从2000增至2200-2500，计算开销可忽略。

---

## 三个KAN应用点的整合架构

```
RASI + KAN 完整流程：

Step 0: KAN-HVG → 稀有性感知的基因选择
        输入：基因统计特征(6维)
        输出：每个基因的稀有重要性分数
        → 标准HVG ∪ KAN选出的稀有HVG

Step 1: PCA_200 → 保留稀有信号

Step 2: KAN-Weight → 可解释的整合权重
        输入：细胞特征(BD/密度/纯度/MNN/特异性, 5维)
        输出：per-cell整合权重
        → z = w(i)*z_harmony + (1-w(i))*z_pre

Step 3: KAN-NP → 可微的NP近似
        输入：细胞嵌入特征(91维)
        输出：预测NP值
        → 用于实时监控 + 可微正则化 + 未知亚群检测

三个KAN共享同一个训练框架（pykan），总参数<1000，
总训练时间<20分钟，推理时间<10秒。
```

## 论文呈现策略

- **Main text**：RASI用简单的线性规则（可重复、易理解）
- **Extended Data**：KAN学到的三组样条曲线（HVG/权重/NP的可解释性发现）
- **Supplementary**：KAN vs 线性规则的定量对比
- **Discussion**：KAN揭示的非线性模式的生物学意义——"什么样的基因定义稀有亚群"、"什么样的细胞应该被保护"

## 实现依赖

```
pip install pykan    # 或 efficient-kan
torch >= 2.0
```

代码量：~100行（不含已有的RASI代码）。

## 论文中的呈现建议

- **Main text**：RASI用线性BD加权，简洁有效
- **Extended Data Figure**：KAN学到的5个输入特征的边际效应曲线
- **Supplementary**：KAN vs 线性BD的定量对比 + 符号公式发现
- **Discussion**：KAN揭示的非线性整合策略的生物学意义
