# 机器学习算法调研报告：单细胞批次整合中未知稀有亚群的保留问题

> 作者：机器学习智能体系统 (agent-mnez8qvx)
> 任务编号：t-mnf028tw

---

## 一、现有主流整合算法在稀有亚群保留上的机制分析

### 1.1 为什么会失败：统一的失败机制

所有主流整合算法的核心操作可以归纳为一个共同范式：**在嵌入空间中寻找跨批次的对应关系，然后将对应的细胞拉近**。这个"拉近"操作本身就是稀有亚群被消灭的根本原因。

#### Harmony（Korsunsky et al., 2019, Nature Methods）
- **机制**：在PCA空间中迭代执行软聚类（k-means变体）和批次校正。每轮迭代将同一簇内不同批次的细胞向簇中心拉齐。
- **失败原因**：稀有亚群细胞数少，在软聚类中被分配到邻近大群的概率极高。一旦被错误归入大群的簇，校正操作会将其嵌入推向大群中心，导致其作为独立群体消失。Harmony的全局迭代策略对局部小结构缺乏保护。

#### scVI / scANVI（Lopez et al., 2018, Nature Methods）
- **机制**：变分自编码器（VAE），将批次信息作为条件变量编码进潜在空间。通过最大化变分下界学习批次不变的表示。
- **失败原因**：VAE的重建损失是按细胞平均的，稀有亚群在总损失中占比极低，模型缺乏动机去保留它们的独特表示。潜在空间的高斯先验假设倾向于产生连续、平滑的嵌入，天然不利于保留孤立的小群。scANVI虽引入标签信息，但仅对已知类型有效。

#### Scanorama（Hie et al., 2019, Nature Biotechnology）
- **机制**：基于互近邻（MNN）的全景拼接策略。先找跨批次的互近邻对，再用这些对作为锚点对齐批次。
- **失败原因**：稀有亚群跨批次的互近邻极难找到——如果某稀有亚群仅存在于部分批次，或在各批次中细胞数极少，MNN匹配大概率会将其与最近的大群配对，导致错误对齐。

#### BBKNN（Polański et al., 2020, Bioinformatics）
- **机制**：在构建KNN图时，强制每个批次贡献相等数量的邻居。
- **失败原因**：当稀有亚群仅存在于少数批次时，强制从其他批次选邻居会引入错误的跨类型连接，在后续聚类中稀有亚群被稀释。

#### scDML（Xiong et al., 2023, Nature Communications）
- **机制**：度量学习+层次合并。通过深度度量学习保持细胞间距离关系，然后层次地合并跨批次的对应簇。
- **失败原因**：虽然scDML专门针对稀有亚群设计，但其验证完全依赖已知标签。层次合并策略仍需预设阈值，对未知稀有亚群无检测能力。在scDREAMER的肺数据测试中，scDML也出现了细胞类型合并问题。

### 1.2 失败的统一模式："聚拢-弥散"不对称

所有整合算法的失败都可以用一个统一模式描述：

- **正常批次校正（聚拢）**：同类型细胞跨批次合并，群保持完整，成员来源更多样。表现为局部邻域的**批次多样性增加**，但**拓扑结构保持**。
- **稀有亚群被消灭（弥散）**：一个群的成员被打散到多个不相关的大群中，群消失。表现为局部邻域的**类型纯度下降**，**拓扑结构断裂**——原本连通的小社区变为散落的孤立点。

这两种模式在数据结构上是可区分的，这是破局的理论基础。

---

## 二、近年来相关的新方法与新思路

### 2.1 绕过整合的稀有亚群检测

**scCross**（Gerniers et al., 2024, Bioinformatics）
- 核心思路：不执行批次整合，直接在多样本原始数据上做biclustering，通过全局求和准则（而非成对比较）寻找跨批次的稀有亚群。
- 意义：**证明了可以在不整合的前提下检测跨批次的稀有亚群**，且对批次效应具有鲁棒性。在肺癌数据中发现了传统方法检测不到的纤毛亚群。
- 局限：仍依赖marker基因驱动，不适用于marker未知的场景。

### 2.2 拓扑数据分析（TDA）方向

**持久同调（Persistent Homology）用于单细胞数据**

- **Damrich et al., 2023**：证明在高维单细胞数据中，传统持久同调对噪声极其敏感，提出用**谱距离**（扩散距离、有效电阻）替代欧氏距离来计算持久同调，可在高维噪声下鲁棒地检测正确拓扑（如细胞周期环路）。
- **scGeom (Huynh & Cang, 2023)**：结合**图曲率**（Ollivier-Ricci curvature）和**持久同调**来分析单细胞数据的多尺度结构，能够指示过渡态细胞和发育潜能。
- **TopACT (Benjamin et al., 2022, Nature Computational Science子刊)**：用多参数持久同调景观（persistence landscape）进行自动细胞类型识别，能在不知道细胞边界的情况下定位稀疏分散的单个细胞。
- **kNN-tPCA (Cottrell et al., 2023)**：用kNN持久拉普拉斯算子增强PCA，在11个scRNA-seq benchmark数据集上显著优于UMAP、tSNE、NMF。

**意义**：TDA能够捕捉整合前后数据拓扑的变化——如果一个连通分量（稀有亚群）在整合后消失，持久同调可以定量检测到这一变化。这直接对应"弥散"信号。

### 2.3 图曲率方向

**ORC-ManL (Saidi, Hickok & Blumberg, 2024)**
- 核心发现：在kNN图中，**穿越ambient space的捷径边具有更负的Ollivier-Ricci曲率**，而沿数据流形的正常边曲率较高。
- 应用：用曲率准则修剪kNN图的虚假边，显著改善了scRNA-seq数据的聚类和流形学习。
- **与我们问题的关联**：整合过程中，稀有亚群被"推入"大群，本质上就是在嵌入空间中创建了虚假的跨类型连接。曲率可以检测这些异常连接。

### 2.4 表示拓扑散度（RTD）

**RTD-AE (Trofimov et al., 2023)**
- 提出**Representation Topology Divergence (RTD)**：定量衡量原始高维数据与低维表示之间的拓扑差异。
- RTD可微分，可作为autoencoder的损失项。
- **与我们问题的关联**：可以将RTD用于衡量整合前后的拓扑保真度。如果整合消灭了稀有亚群，RTD会增大。

### 2.5 多尺度嵌入

**pNMF (Zhang, Miao & Li, 2026)**
- 提出持久非负矩阵分解：利用持久同调识别数据连通性发生质变的**典范尺度集**，在每个尺度上生成对齐的嵌入。
- **与我们问题的关联**：稀有亚群作为小规模的连通结构，可能只在特定尺度上可见。多尺度方法能捕捉单一尺度下不可见的结构。

### 2.6 最优传输方向

**GENOT (Klein et al., 2023)** 和 **Unbalanced OT (Eyring et al., 2023)**
- 非平衡最优传输（UOT）放松了质量守恒约束，允许源分布中某些点不被传输到目标分布——这对处理批次特异的稀有亚群至关重要。
- **CellANOVA (2024, Nature Biotechnology)** 提出整合后信号恢复框架，但仅限已知差异。

### 2.7 对比学习与开放集域适应

**SAFAARI (Aminzadeh et al., 2024/2026)**
- 结合**监督对比学习**和**对抗性开放集域适应**，在单细胞整合中同时做标签迁移和批次校正。
- 关键特性：**能识别新的（未见过的）细胞类型**，通过开放集检测机制。
- 局限：仍依赖参考数据集中的已知类型来训练开放集检测器。

### 2.8 基础模型方向

**Kendiukhov, 2026**：对scGPT和Geneformer的141个拓扑/几何假设进行系统检验，发现这些模型学到了非平凡的拓扑结构（持久同调在几乎所有层显著），但信号在免疫组织中最强。

---

## 三、建议的算法设计框架

### 3.1 核心思路：拓扑守恒整合（Topology-Conservative Integration, TCI）

**核心洞察**：不需要知道稀有亚群在哪里，只需要检测整合操作是否破坏了数据的拓扑结构。正常的批次校正保持拓扑（聚拢），而稀有亚群被消灭会改变拓扑（弥散）。

### 3.2 三层架构

```
┌─────────────────────────────────────────────┐
│   第一层：拓扑指纹提取（整合前）              │
│   - 对每个批次独立计算多尺度拓扑特征           │
│   - 持久同调 barcode（H0连通分量, H1环路）     │
│   - kNN图的Ollivier-Ricci曲率分布             │
│   - 局部密度与连通性的联合特征                 │
├─────────────────────────────────────────────┤
│   第二层：拓扑守恒整合引擎                    │
│   - 基础整合器：任意现有方法（scVI, Harmony等） │
│   - 拓扑正则化项：RTD(pre, post) → 最小化      │
│   - 曲率约束：整合后新增的负曲率边→惩罚        │
│   - 非平衡传输：允许批次特异细胞不被强制对齐    │
├─────────────────────────────────────────────┤
│   第三层：弥散检测器（整合后）                 │
│   - 比较整合前后每个细胞的局部邻域变化          │
│   - 检测"弥散签名"：                          │
│     · 整合前局部邻域同质 → 整合后邻域异质       │
│     · 整合前属于连通的小社区 → 整合后散入大群    │
│     · 整合前后kNN图的持久同调差异               │
│   - 输出：疑似被消灭的稀有亚群候选列表          │
└─────────────────────────────────────────────┘
```

### 3.3 弥散检测器的具体设计

这是整个框架中最关键的组件——它回答"是否有未知稀有亚群被消灭"。

**算法 DispersionDetector**：

```
输入：
  X_pre  = 整合前各批次的表达矩阵（或嵌入）
  X_post = 整合后的统一嵌入
  k      = 邻居数

步骤：
1. 对 X_pre 构建 per-batch kNN 图 G_pre
2. 对 X_post 构建统一 kNN 图 G_post
3. 对每个细胞 i，计算：
   a) 邻域保留率 NP(i) = |N_pre(i) ∩ N_post(i)| / k
   b) 连通性变化 ΔC(i) = 整合前 i 所在连通分量大小 - 整合后 i 所在连通分量大小
   c) 局部曲率变化 ΔK(i) = mean(ORC_post(i)) - mean(ORC_pre(i))
4. 定义弥散分数 DS(i) = w1 * (1 - NP(i)) + w2 * max(0, ΔC(i)) + w3 * max(0, -ΔK(i))
5. 对 DS 进行聚类，高 DS 且在整合前空间紧凑的细胞群 → 疑似被消灭的稀有亚群
6. 用持久同调验证：检查这些候选群在整合前的 barcode 中是否对应一个显著的 H0 特征（长寿命的连通分量），而在整合后消失
```

### 3.4 拓扑守恒正则化的具体形式

对于任何可微整合方法（如 scVI），在损失函数中加入：

```
L_total = L_integration + λ₁ * L_RTD + λ₂ * L_curvature

其中：
  L_RTD = RTD(G_pre, G_post)  // 表示拓扑散度
  L_curvature = Σ_edges max(0, -ORC(e) - τ)  // 惩罚异常负曲率边
```

这样整合器在对齐批次的同时，被约束不能破坏数据的拓扑结构。

### 3.5 可扩展性考虑

针对490万细胞规模：
- 持久同调的直接计算在大规模下不可行 → 使用**子采样 + 持久同调近似**（如 Ripser 的稀疏模式）
- Ollivier-Ricci曲率的精确计算需要解最优传输 → 使用**有效电阻曲率**（Efficient Curvature, Fei et al., 2025）作为替代，计算复杂度从 O(n³) 降至线性
- RTD 可通过采样近似计算
- 整体框架可分块执行：先对数据进行粗聚类，在每个粗簇内独立运行拓扑分析

### 3.6 弥散分数的数学基础（来自数学智能体分析，任务 t-mnf19n46）

#### 三分量独立性
- NP（邻域保留率）与 ΔC（连通性变化）中等正相关（r≈0.5-0.7），因为都是"邻域变化"信号的不同粒度
- ΔK（曲率变化）与其余两者近似独立（r≈0.1-0.3）
- **结论：三分量不完全冗余，各有独特贡献**

#### 权重理论最优比例
- 基于 Fisher 线性判别分析：w* = S_w^{-1}(μ_disp - μ_norm)
- **归一化后初始建议：w1:w2:w3 ≈ 2:1:1**（NP 贡献最高）
- 精确最优比例依赖具体数据，建议用已知稀有类型作代理标签交叉验证微调

#### DS 与 CNEM 的关系
- DS 检测"结构破坏"（群碎了、邻居换了、曲率塌了）
- CNEM 检测"语义错位"（被放在了不像自己的人旁边）
- **互补大于冗余**：DS 独有贡献——稀有群被吸收到表达相似的大群时 DS 可检测但 CNEM 漏检；CNEM 独有贡献——嵌入空间高度扭曲或降维丢失信息时 CNEM 直接在表达空间度量

#### FusedScore 最优化
- Fisher 解析解：(α*, β*)^T = Σ^{-1} · Δμ
- 无标签实用方案：z-score 归一化后等权融合（α=β=0.5），或取大值策略 FusedScore=max(CNEM_z, DS_z)

### 3.7 验证路径

1. **计算验证**：在已知含稀有亚群（如胰腺 epsilon、肺 ionocyte）的数据上，人为"假装不知道"这些稀有亚群，运行弥散检测器，验证其能否检测到它们在整合中的丢失。
2. **对比实验**：将 TCI 与无拓扑约束的基线整合方法对比，衡量稀有亚群保留率的提升。
3. **生物学验证**：将弥散检测器发现的未知候选亚群，提交给免疫学/肿瘤生物学智能体进行湿实验验证。

---

## 四、参考文献

### 整合方法与稀有亚群问题
1. Korsunsky I, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. **Nature Methods**, 2019.
2. Lopez R, et al. Deep generative modeling for single-cell transcriptomics. **Nature Methods**, 2018.
3. Hie B, et al. Efficient integration of heterogeneous single-cell transcriptomes using Scanorama. **Nature Biotechnology**, 2019.
4. Xiong L, et al. Online single-cell data integration through projecting heterogeneous datasets into a common cell-embedding space (scDML). **Nature Communications**, 2023.
5. Dou J, et al. Unbiased integration of single cell multi-omics data (scDREAMER). **Nature Communications**, 2023.
6. Wang H, et al. scCRAFT: Dual-resolution topology-preserving batch correction. **Communications Biology**, 2025.
7. Andreatta M, et al. Semi-supervised integration of single-cell transcriptomics data (STACAS). **Nature Communications**, 2024.
8. Han M, et al. CellANOVA: A signal recovery framework for batch-corrected single-cell data. **Nature Biotechnology**, 2024.
9. Hu X, et al. RBET: Reference-informed evaluation of batch correction with overcorrection awareness. **Communications Biology**, 2025.
10. Luecken MD, et al. Benchmarking atlas-level data integration in single-cell genomics (scIB). **Nature Methods**, 2022.

### 绕过整合的稀有亚群检测
11. Gerniers A, et al. scCross: efficient search for rare subpopulations across multiple single-cell samples. **Bioinformatics**, 2024.
12. Liang S, et al. Single-cell manifold-preserving feature selection for detecting rare cell populations (SCMER). **Nature Computational Science**, 2021.

### 拓扑数据分析
13. Damrich S, et al. Persistent homology for high-dimensional data based on spectral methods. **arXiv:2311.03087**, 2023.
14. Huynh T & Cang Z. Topological and geometric analysis of cell states in single-cell transcriptomic data (scGeom). **arXiv:2309.07950**, 2023.
15. Benjamin K, et al. Multiscale topology classifies and quantifies cell types in subcellular spatial transcriptomics (TopACT). **arXiv:2212.06505**, 2022.
16. Cottrell S, et al. K-Nearest-Neighbors induced topological PCA for scRNA sequence data analysis (kNN-tPCA). **arXiv:2310.14521**, 2023.
17. Zhang J, et al. Persistent nonnegative matrix factorization via multi-scale graph regularization (pNMF). **arXiv:2602.22536**, 2026.

### 图曲率
18. Saidi TL, et al. Recovering manifold structure using Ollivier-Ricci curvature (ORC-ManL). **arXiv:2410.01149**, 2024.
19. Fei C, et al. Efficient curvature-aware graph network (Effective Resistance Curvature). **arXiv:2511.01443**, 2025.
20. Nguyen K, et al. Revisiting over-smoothing and over-squashing using Ollivier-Ricci curvature. **arXiv:2211.15779**, 2022.

### 表示拓扑
21. Trofimov I, et al. Learning topology-preserving data representations (RTD-AE). **arXiv:2302.00136**, 2023.

### 最优传输
22. Klein D, et al. GENOT: Entropic (Gromov) Wasserstein flow matching with applications to single-cell genomics. **arXiv:2310.09254**, 2023.
23. Eyring L, et al. Unbalancedness in neural Monge maps improves unpaired domain translation. **arXiv:2311.15100**, 2023.

### 对比学习与开放集
24. Aminzadeh F, et al. SAFAARI: Contrastive adversarial open-set domain adaptation for single-cell integration & annotation. **Genomics, Proteomics & Bioinformatics**, 2026.

### 基础模型
25. Kendiukhov I. What topological and geometric structure do biological foundation models learn? Evidence from 141 hypotheses. **arXiv:2602.22289**, 2026.

### 评估方法
26. Shen X, et al. A benchmark of semi-supervised scRNA-seq integration methods in real-world scenarios. **PLOS Computational Biology**, 2026.
27. Zappia L, et al. Feature selection methods affect the performance of scRNA-seq data integration and querying. **Nature Methods**, 2025.
28. Babcock BR, et al. Data matrix normalization and merging strategies minimize batch-specific systemic variation in scRNA-seq data (CMS metric). **bioRxiv**, 2021.
