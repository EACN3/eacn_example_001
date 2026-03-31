# PoC 实验设计：基于邻域结构变化检测批次整合中的稀有亚群破坏

## 核心假设

批次整合对常见细胞类型产生**聚拢效应**（同类型跨批次细胞邻域趋于一致），对稀有亚群产生**弥散效应**（邻域结构被打散）。这两种模式在邻域保持度指标上可被区分，且无需细胞标签即可检测。

---

## 1. 数据集

**人胰腺单细胞数据集**，通过 `scvi-tools` 或 `scanpy` 获取：

| 数据集 | 来源 | 预期细胞数 |
|--------|------|-----------|
| Baron | `scvi.data.baron()` | ~8,000 |
| Muraro | `scvi.data.muraro()` | ~2,000 |
| Segerstolpe | `scvi.data.segerstolpe()` | ~2,000 |

**已知稀有类型**（作为验证 ground truth）：
- Epsilon 细胞（<1%）
- Schwann 细胞（极少）
- Activated stellate 细胞（少量）

**常见类型**（作为对照）：Alpha、Beta、Ductal、Acinar 等。

---

## 2. 预处理流程

```
合并三个数据集 → 标准 QC（基因/细胞过滤）→ 高变基因选择（2000 HVG）
→ 标准化 + log1p → PCA（50 components）→ 记录 batch 标签和 cell type 标签
```

保留原始 cell type 标签用于最终验证，但框架运行过程中**不使用**标签信息。

---

## 3. 批次整合方法

对同一预处理数据分别运行四种方法：

| 方法 | 输出空间 | Python 包 |
|------|---------|-----------|
| Harmony | 校正后 PCA embedding | `harmonypy` |
| scVI | 潜在空间 embedding | `scvi-tools` |
| BBKNN | 校正后 kNN 图 | `bbknn` |
| Scanorama | 校正后 embedding | `scanorama` |

每种方法使用默认参数。BBKNN 直接输出 kNN 图，其余方法在 embedding 上构建 kNN 图。

---

## 4. 邻域结构比较框架

### 4.1 整合前邻域（per-batch kNN）

对每个批次独立计算：
- 在该批次内部的 PCA 空间上构建 **k=30 的 kNN 图**
- 每个细胞得到一组"批次内邻居"集合 $N_{\text{pre}}(i)$

### 4.2 整合后邻域（global kNN）

在整合后的 embedding 上：
- 对全部细胞构建 **k=30 的 kNN 图**
- 每个细胞得到一组"整合后邻居"集合 $N_{\text{post}}(i)$
- 对于 BBKNN，直接使用其输出的 kNN 图

### 4.3 逐细胞邻域保持度评分（Neighborhood Preservation Score, NPS）

对每个细胞 $i$，计算：

$$
\text{NPS}(i) = \frac{|N_{\text{pre}}(i) \cap N_{\text{post}}(i)|}{k}
$$

即整合前后共享邻居的比例。

**注意**：$N_{\text{pre}}$ 仅包含同批次细胞，而 $N_{\text{post}}$ 包含所有批次细胞。因此需要归一化处理——仅在 $N_{\text{post}}$ 中取与细胞 $i$ 同批次的邻居子集进行比较：

$$
\text{NPS}(i) = \frac{|N_{\text{pre}}(i) \cap N_{\text{post,same\_batch}}(i)|}{\min(k, |N_{\text{post,same\_batch}}(i)|)}
$$

### 4.4 聚拢-弥散分类

在 NPS 基础上，进一步计算每个细胞的**跨批次邻居比例**：

$$
\text{CBR}(i) = \frac{|N_{\text{post,other\_batch}}(i)|}{k}
$$

构建二维指标空间 (NPS, CBR)：
- **聚拢模式**：NPS 较低 + CBR 高 → 邻域变化大，但获得了大量跨批次同类型邻居（正常整合效果）
- **弥散模式**：NPS 较低 + CBR 高，但新邻居类型异质性高 → 邻域被打散到不同类型中

为量化"新邻居的异质性"，对每个细胞计算其整合后邻域的**局部密度一致性**：
- 整合后 k 邻居的 embedding 距离方差（弥散 → 方差大）
- 整合后 k 邻居之间的互为邻居比例（弥散 → 比例低）

### 4.5 综合异常评分

$$
\text{Disruption}(i) = (1 - \text{NPS}(i)) \times \text{NeighborDispersion}(i)
$$

其中 NeighborDispersion 用邻居间距离方差的标准化值表示。高 Disruption 且低 NPS 的细胞被标记为"潜在受损亚群"。

---

## 5. 验证实验

### 5.1 盲测验证

1. 从标签中移除 epsilon 细胞的类型信息（标记为 "unknown"）
2. 运行上述框架，获得每个细胞的 Disruption 评分
3. 对 Disruption 评分进行聚类或阈值检测（如 top 5% 高分细胞）
4. **检验**：被检测出的高 Disruption 细胞中，epsilon 细胞的富集程度

### 5.2 对照比较

| 组别 | 细胞 | 预期 Disruption |
|------|------|----------------|
| 稀有型（epsilon, schwann, activated stellate） | ~50-100 个 | 高 |
| 常见型（alpha, beta, ductal） | ~5000+ 个 | 低 |
| 中等型（delta, gamma/PP） | ~200-500 个 | 中等 |

用 Wilcoxon 秩和检验比较稀有型 vs 常见型的 Disruption 分布。

### 5.3 跨方法一致性

检查四种整合方法中，epsilon 等稀有细胞是否一致地获得高 Disruption 评分。一致性高则说明该指标捕捉的是数据结构特征而非方法特异性偏差。

---

## 6. 预期输出

1. **每种整合方法的 UMAP 图**：着色为 Disruption 评分，叠加 cell type 标签
2. **Disruption 评分箱线图**：按 cell type 分组，展示稀有 vs 常见的差异
3. **ROC 曲线**：以 Disruption 评分预测"是否为稀有类型"的分类性能
4. **富集分析表**：top N% 高 Disruption 细胞中各类型的占比
5. **跨方法热力图**：细胞类型 x 整合方法，值为平均 Disruption 评分

---

## 7. 成功标准

| 标准 | 指标 | 阈值 |
|------|------|------|
| 稀有类型可检测 | top 5% Disruption 细胞中稀有类型富集 OR | > 3.0 |
| 稀有 vs 常见可区分 | Wilcoxon 检验 p-value | < 0.01 |
| 无标签检测有效 | 盲测 epsilon 召回率 | > 50% |
| 跨方法稳健 | 至少 3/4 方法中稀有类型 Disruption 显著高于常见类型 | -- |

---

## 8. 实现要点

- **计算框架**：Python，基于 `scanpy`、`anndata`、`sklearn.neighbors`
- **kNN 计算**：使用 `scanpy.pp.neighbors` 或 `sklearn.neighbors.NearestNeighbors`
- **统计检验**：`scipy.stats.mannwhitneyu`
- **可视化**：`scanpy.pl.umap` + `matplotlib`
- **预计运行时间**：< 30 分钟（含四种整合方法）
- **内存需求**：< 16 GB

---

## 9. 风险与备选方案

| 风险 | 应对 |
|------|------|
| 稀有细胞数量过少导致统计功效不足 | 可通过下采样常见类型平衡比例后重新验证 |
| NPS 对 k 值敏感 | 扫描 k=10,20,30,50 进行敏感性分析 |
| 某些整合方法不破坏稀有亚群 | 这本身是有价值的发现，记录方法间差异 |
| 弥散模式定义不够精确 | 可引入更多拓扑指标（如局部连通性、持续同调） |
