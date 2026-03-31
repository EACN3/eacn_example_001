# 免疫图谱「聚拢 vs 弥散」检测框架验证实验方案

## 实验目标

在 atlas_merged_immune.h5ad（225万免疫细胞，103批次）上验证：批次整合是否系统性地破坏已知稀有亚群（pDC、Mast），并量化这种破坏的结构特征。

---

## 前期实验经验总结（来自胰腺数据 PoC）

在 16,382 细胞的胰腺数据（9批次，14种细胞类型，含 epsilon/schwann/mast 等稀有类型）上，我们已完成两轮 PoC 实验，获得以下关键发现：

1. **NPS（邻域保持度）指标失败**：跨批次整合后所有类型的同批次邻居都减少，NPS 对所有类型均低，无法区分稀有与常见类型。
2. **表达一致性（coherence_drop）指标部分失败**：稀有类型（mast、macrophage、t_cell）整合后 coherence_drop 为负 — 被吸收进大群后邻域反而更"一致"了，但这恰好说明它们的独特性被抹平。
3. **核心洞察**：检测信号不在"邻域是否被打散"，而在**"细胞自身表达与整合后邻域是否匹配"** — 即 cell-neighborhood mismatch。

---

## 实验设计

### Phase 1：子集概念验证（~5万细胞）

#### 1.1 子集选取策略

从103个批次中选取 **8-10个批次**，满足：
- 至少包含 3 个含有 pDC 的批次（确保 pDC 跨批次存在）
- 至少包含 2 个不含 pDC 的批次（作为批次特异性对照）
- 包含 Mast 细胞的批次
- 批次大小中等（500-5000 细胞/批次），避免极端不平衡
- 目标总细胞数：3-5万

```python
# 选取逻辑（伪代码）
batch_stats = adata.obs.groupby('batch')['celltype'].value_counts().unstack(fill_value=0)
pdc_batches = batch_stats[batch_stats['pDC'] > 10].index  # 至少10个pDC
mast_batches = batch_stats[batch_stats['Mast'] > 10].index
# 选取 5 个 pDC 批次 + 3 个无 pDC 批次 + 2 个含 Mast 批次
```

#### 1.2 预处理流程

```
子集提取 → QC（min_genes=200, min_cells=3, mito<20%）
→ Normalize（target_sum=1e4）→ log1p
→ HVG 选择（2000 genes, batch_key=batch_col）
→ Scale（max_value=10）→ PCA（50 comps）
→ 保存归一化表达矩阵 X_expr（用于后续 mismatch 计算）
```

#### 1.3 整合方法（3种）

| 方法 | 输出 | 参数 | 选择理由 |
|------|------|------|---------|
| **Harmony** | 校正后 PCA embedding | 默认（max_iter=10, theta=2） | 最常用，速度快，可扩展到百万级 |
| **scVI** | 潜在空间 embedding（10-30 dims） | n_latent=30, n_layers=2, max_epochs=200 | 深度学习代表，scIB benchmark 表现好 |
| **Scanorama** | 校正后 embedding | 默认 | 基于 MNN 的经典方法，与 Harmony 机制不同 |

#### 1.4 核心指标：三维指标体系

基于 PoC 经验，设计三个互补指标：

**指标 1：Cell-Neighborhood Expression Mismatch (CNEM)**

核心思想：整合后，细胞自身的表达与其新邻域的平均表达之间的距离。

```
对每个细胞 i：
1. 在整合后 embedding 上找 k=30 近邻 N_post(i)
2. 计算细胞 i 的表达向量与 N_post(i) 中所有邻居的平均表达向量的余弦距离
   CNEM(i) = 1 - cosine_sim(x_i, mean(x_j for j in N_post(i)))
```

- 正常整合（聚拢）：细胞移向同类型跨批次细胞 → CNEM 低
- 稀有亚群被消灭（弥散）：细胞被推入异类型大群 → CNEM 高

**指标 2：Local Neighborhood Type Entropy Change (LNTEC)**

核心思想：无监督版本 — 不用标签，用整合前后邻域的表达谱多样性变化。

```
对每个细胞 i：
1. 整合前：在同批次 PCA 空间中找 k 近邻，计算邻域内表达距离的方差 var_pre
2. 整合后：在整合 embedding 中找 k 近邻，计算邻域内表达距离的方差 var_post
3. LNTEC(i) = var_post / var_pre

LNTEC > 1：整合后邻域表达更异质（被推入混合区域）
LNTEC < 1：整合后邻域表达更均质（正常聚拢）
LNTEC ≈ 1：邻域变化不大
```

**指标 3：Displacement-Coherence Ratio (DCR)**

核心思想：综合考虑细胞在 embedding 中的位移量和位移后的邻域质量。

```
对每个细胞 i：
1. 位移量：D(i) = ||embed_post(i) - embed_pre(i)||（需对齐坐标系，用 Procrustes）
2. 位移后邻域质量：Q(i) = mean cosine_sim(x_i, x_j) for j in N_post(i)
3. DCR(i) = D(i) / Q(i)

高位移 + 低质量 = 高 DCR → 被强行移到不匹配的区域（弥散）
高位移 + 高质量 = 低 DCR → 移到了同类型细胞附近（正常聚拢）
低位移 = 低 DCR → 位置没怎么变（无影响）
```

#### 1.5 验证流程

**Step A：有标签验证（ground truth）**

1. 用已知标签，比较 pDC vs 常见类型（T cell、NK）的三个指标分布
2. 统计检验：Wilcoxon 秩和检验，pDC 的 CNEM/LNTEC/DCR 是否显著高于常见类型
3. ROC-AUC：以指标值预测"是否为稀有类型"的分类性能
4. 富集分析：top 5%/10% 异常细胞中 pDC 的富集 OR
5. 跨方法一致性：3种方法中 pDC 是否一致获得高分

**Step B：盲测验证**

1. 移除 pDC 的标签（标记为 "unknown"）
2. 运行三个指标，获取每个细胞的得分
3. 检查 top N% 高分细胞中 pDC 的召回率
4. 对高分细胞进行无监督聚类，看是否能形成独立 cluster
5. 对该 cluster 做差异表达分析，验证 marker 基因（如 CLEC4C, IL3RA, IRF7）

**Step C：阴性对照**

1. 对常见类型（如 CD4 T cell）施加同样的分析流程
2. 确认常见类型不会被错误标记为"受损"
3. 计算假阳性率

### Phase 2：全数据集扩展（225万细胞）

Phase 1 验证成功后：

1. 使用 Harmony（最快）在全量数据上运行
2. 仅计算 CNEM（最稳健的指标），跳过 LNTEC/DCR 以节约计算
3. 使用批量化计算：分 chunk 处理 kNN（每次 10 万细胞）
4. 输出：全量细胞的 CNEM 分数 + 高分细胞的聚类和 marker 分析

---

## 参数选择

| 参数 | 值 | 理由 |
|------|-----|------|
| k (kNN) | 30 | scIB benchmark 标准值 |
| n_hvg | 2000 | 平衡信息量与计算效率 |
| n_pcs | 50 | 标准 PCA 维度 |
| top N% 阈值 | 5%, 10% | 两个阈值交叉验证 |
| cosine_sim 基础 | HVG 归一化表达 | 原始表达空间，不受整合影响 |

敏感性分析：k = [10, 20, 30, 50]，n_hvg = [1000, 2000, 3000]

---

## 预期结果

### 理想结果

1. **CNEM**: pDC 的 CNEM 显著高于 T cell/NK（p < 0.001），top 5% 富集 OR > 5
2. **盲测**: pDC 召回率 > 50%（top 10% 中）
3. **跨方法一致性**: 至少 2/3 方法中 pDC 指标显著
4. **阴性对照**: CD4 T cell 假阳性率 < 5%

### 可接受结果

1. 至少 1 个指标在至少 2 种方法中对 pDC 有显著检测力
2. 盲测 pDC 召回率 > 30%

### 需要调整方向的结果

1. 所有指标均无法区分 pDC 与常见类型 → 需要重新审视"聚拢 vs 弥散"假设
2. 仅在特定方法下有效 → 指标可能对方法特异性偏差敏感，需引入方法无关的指标

---

## 计算资源估算

| 阶段 | 内存 | GPU | 时间 |
|------|------|-----|------|
| Phase 1 子集（5万细胞） | ~16 GB | Harmony/Scanorama: CPU; scVI: 1×A800 | ~1小时 |
| Phase 2 全量（225万细胞） | ~128 GB | Harmony: CPU; scVI: 4×A800 | ~6-12小时 |

---

## 风险与备选

| 风险 | 应对 |
|------|------|
| pDC 在该数据中不够"稀有"（12898个=0.57%） | 用正常组织的 pDC 子集（2528个=0.11%）做更严格的验证 |
| 225万细胞的 kNN 计算内存爆炸 | 使用 FAISS GPU 加速，或分 chunk 近似计算 |
| scVI 训练不收敛 | 增大 max_epochs，调整 lr，或使用 scANVI 替代 |
| CNEM 对 HVG 选择敏感 | 多组 HVG 做敏感性分析 |
| 批次间细胞类型分布极度不平衡 | 子集选取时控制平衡性 |
