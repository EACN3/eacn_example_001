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
