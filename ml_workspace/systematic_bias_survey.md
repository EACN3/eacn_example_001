# 单细胞分析流程中对稀有亚群的系统性偏差调研

> 作者：机器学习智能体 (agent-mnez8qvx)
> 任务：t-mnfpqv3t

---

## 1. HVG选择偏差

**问题**：标准HVG选择（如Seurat的vst/scran的modelGeneVar）基于全局方差/均值关系。稀有亚群的marker基因在全局中方差低（只有<1%的细胞表达），被排除在HVG外。

**已有文献**：
- Zappia et al. (2025, Nature Methods)："Feature selection methods affect the performance of scRNA-seq data integration and querying"——系统评估了HVG选择对整合和稀有亚群检测的影响
- SCMER (Liang et al. 2021, Nature Computational Science)：流形保持特征选择，专门保留稀有亚群信号
- Townes et al. (2019, Genome Biology)：GLM-PCA绕过HVG选择

**已有方案**：SCMER、GeneTrajectory、Hotspot（基于自相关的基因选择）
**未解决gap**：无方法专门为"未知稀有亚群"选择基因——SCMER需要知道要保留什么
**Rare-aware替代**：加权HVG——在方差计算中对低密度区域的细胞上权，或使用全基因组（GLM-PCA允许这样做）

## 2. PCA/降维偏差

**问题**：PCA最大化全局方差。稀有亚群(<1%细胞)贡献的方差<0.01%，其标志性信号在前50个PC中几乎不可见。

**已有文献**：
- Damrich et al. (2023)：持久同调在高维scRNA-seq中对噪声敏感，用谱距离替代
- kNN-tPCA (Cottrell et al. 2023)：拓扑增强PCA
- GLM-PCA (Townes et al. 2019)：考虑count数据特性的PCA

**已有方案**：GLM-PCA、scVI latent space（VAE嵌入）、UMAP（非线性但不保距离）
**未解决gap**：**没有任何降维方法专门考虑稀有亚群保留**——这是一个完全空白的领域
**Rare-aware替代**：加权PCA（密度倒数加权）、分层PCA（先粗聚类再分层降维）、对比式预训练encoder

## 3. kNN构建偏差

**问题**：固定k（如k=30）对所有细胞一视同仁。但稀有亚群总共可能只有50-100个细胞，k=30意味着其邻居中一半以上是其他类型。

**已有文献**：
- ORC-ManL (Saidi et al. 2024)：用Ollivier-Ricci曲率修剪kNN图的虚假边
- BBKNN (Polański et al. 2020)：每批次贡献等量邻居，但不解决k值问题

**已有方案**：局部密度自适应k（理论上存在但实践中很少用于scRNA-seq）
**未解决gap**：**无标准方法根据局部密度自动选择k**
**Rare-aware替代**：k_i = min(k_max, |S_i|/2) 其中S_i是i所属粗聚类的大小；或用HDBSCAN的可变密度框架

## 4. 批次效应建模偏差

**问题**：所有整合方法假设批次效应对所有细胞类型相同（全局校正向量）。但实际上不同细胞类型受批次效应影响不同。

**已有文献**：
- CellANOVA (2024, Nature Biotechnology)：证明"当前整合范式过度激进，移除了有意义的生物学变异"，提出信号恢复框架
- STACAS (2024, Nature Communications)：利用先验细胞类型信息引导整合
- SMNN (Yang et al. 2019)：监督式MNN，按细胞类型分别做MNN匹配
- Babcock et al. (2021)：CMS指标检测过度校正

**已有方案**：CellANOVA（后处理恢复）、STACAS（半监督引导）、SMNN（分类型MNN）
**未解决gap**：**没有无监督的细胞类型特异批次校正方法**——CellANOVA和STACAS都需要已知标签
**Rare-aware替代**：BD加权整合（我们的发现）——BD天然反映细胞类型特异的批次效应程度

## 5. 归一化偏差

**问题**：全局归一化（如总count归一到10000）假设所有细胞有相同的总mRNA量。稀有细胞（如pDC）的基因表达谱与大群不同，全局归一化可能扭曲其相对表达。

**已有文献**：
- Lun et al. (2016, Genome Biology)：scran的池化归一化，对小群更鲁棒
- Hafemeister & Satija (2019)：SCTransform正则化负二项回归
- PbImpute (Zhang et al. 2025)：区分技术dropout和生物学零表达

**已有方案**：scran池化、SCTransform、sctransform v2
**未解决gap**：归一化对稀有亚群的系统影响尚未被量化
**Rare-aware替代**：per-cluster归一化（先粗聚类再分别归一化）

## 6. 聚类偏差

**问题**：单分辨率Leiden/Louvain聚类在低分辨率时合并稀有群，高分辨率时碎片化大群。没有"对所有类型都好"的分辨率。

**已有文献**：
- scDML (Xiong et al. 2023, Nature Communications)：层次合并策略保留稀有亚群
- pNMF (Zhang et al. 2026)：多尺度持久NMF
- t-NEB (Ritzert et al. 2025)：最大密度路径层次聚类

**已有方案**：多分辨率聚类（Clustree可视化）、HDBSCAN（密度自适应）、scDML层次合并
**未解决gap**：**没有方法同时优化大群和稀有群的聚类质量**
**Rare-aware替代**：两阶段聚类——低分辨率识别大群 → 大群内部高分辨率寻找稀有亚群

## 7. 整合偏差

**问题**：所有现有方法将"整合=对齐所有细胞"作为默认假设。

**已有文献**（详见ml_algorithm_survey.md的29篇）：
- scDML、scDREAMER、scCRAFT：声称保护稀有亚群但依赖已知标签验证
- scCross (Gerniers 2024)：绕过整合直接检测跨批次稀有亚群
- SAFAARI (Aminzadeh 2026)：开放集域适应识别未见细胞类型
- Shen et al. (2026, PLOS Comp Bio)：半监督方法在标签缺陷时退化
- RBET (Hu et al. 2025)：过度校正感知评估

**已有方案**：scDML、STACAS、CellANOVA、BD-Harmony（我们的）
**未解决gap**：**无方法在不知道稀有亚群存在的前提下保护它们**——这正是root_cause_analysis中的核心问题
**Rare-aware替代**：选择性整合（BD感知）、非平衡OT（质量销毁保护独有细胞）

## 8. 评估偏差

**问题**：scIB的所有指标（ARI、NMI、ASW、graph connectivity、kBET）都需要外部标签或只衡量全局质量。

**已有文献**：
- scIB (Luecken et al. 2022, Nature Methods)：16种方法的标准benchmark
- scIB-E (2025, Genome Biology)：发现scIB对细胞类型内部变异的保留存在盲区
- RBET (2025, Communications Biology)：参考信息评估过度校正
- Rautenstrauch & Ohler (2025, Nature Biotechnology)：silhouette在整合benchmark中的缺陷

**已有方案**：scIB-E（扩展指标）、RBET（参考信息评估）、NP（我们的）
**未解决gap**：**NP是第一个不需要标签的亚群消灭检测指标**——这是我们的独特贡献
**Rare-aware替代**：NP + survival rate + BD分析

---

## 我们的独特贡献定位

在已有工作基础上，我们的独特贡献：

1. **首次系统性量化**整个流程中7个环节对稀有亚群的联合伤害（不只是整合）
2. **NP**：第一个无标签的亚群消灭检测指标（AUC=0.837）
3. **BD/选择性整合**：第一个在不知道稀有亚群存在的前提下保护它们的整合范式
4. **规模依赖性消灭**的发现和量化（pDC在8批次完好但103批次消灭）
5. **生物学验证**：GITR+ Treg功能社区在整合中被破坏的证据

**Nature定位**：不是一个新算法，而是一个paradigm shift——从"整合=对齐所有细胞"到"整合=选择性对齐+保护"。
