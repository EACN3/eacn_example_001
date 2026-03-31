# CNEM 失效的理论诊断与修正方案

## 1. CNEM 失效的充要条件

### 1.1 形式化设定

设细胞类型集合 $\mathcal{T} = \{T_1, \ldots, T_M\}$，其中 $T_r$ 为稀有亚群。每个类型 $T_m$ 的表达分布为 $P_m$，均值 $\mu_m \in \mathbb{R}^G$，类内协方差 $\Sigma_m$。

CNEM 的判别力取决于稀有群细胞与其整合后邻居之间的信号：
$$\text{CNEM}(i) = 1 - \cos(x_i, \bar{x}_{N(i)})$$

### 1.2 失效的充要条件

**定理（CNEM 失效的充要条件）**：CNEM 无法区分稀有亚群 $T_r$ 与常见类型 $T_c$（即 $\mathbb{E}[\text{CNEM}|T_r] \leq \mathbb{E}[\text{CNEM}|T_c]$），当且仅当以下条件成立：

$$\frac{\mathbb{E}\left[\left\|\frac{x_r}{\|x_r\|} - \frac{\bar{x}_{N_r}}{\|\bar{x}_{N_r}\|}\right\|^2\right]}{\mathbb{E}\left[\left\|\frac{x_c}{\|x_c\|} - \frac{\bar{x}_{N_c}}{\|\bar{x}_{N_c}\|}\right\|^2\right]} \leq 1$$

其中 $x_r \in T_r$，$x_c \in T_c$，$N_r, N_c$ 为各自的整合后邻域。

展开后，这等价于三个条件的组合：

**(C1) 表达锥角不足**：稀有群均值与邻域均值的夹角不大于常见类型的对应夹角：
$$\theta(T_r, \text{neighbors of } T_r) \leq \theta(T_c, \text{neighbors of } T_c)$$

**(C2) 高维集中效应**：当有效维度 $G_{\text{eff}}$ 大时，所有成对余弦相似度集中在 $1 - O(1/\sqrt{G_{\text{eff}}})$ 附近，类间差异被压缩。

**(C3) 邻域纯度悖论**：稀有群整合后的邻居可能仍以同类为主（整合成功保留了该群），导致 CNEM 低——此时低 CNEM 是**正确的**，不是失效。

### 1.3 免疫图谱的具体失效机制

免疫数据 CNEM 失效的最可能原因是 **(C1)+(C2)** 的叠加：

- 免疫细胞（pDC, T cell, B cell）在 HVG 空间中共享大量免疫相关基因的高表达（如 MHC、细胞因子受体等）→ 方向高度相似（ρ 大）
- HVG 数量大（G=2000+），集中效应使余弦差异进一步收窄
- 结果：pDC 与 T/B cell 的 CNEM 差异被淹没在噪声中

但 **(C3)** 也不能排除：如果 pDC 整合后邻居仍以 pDC 为主，低 CNEM 恰恰说明整合没有消灭 pDC。**建议计算生物学检查 pDC 的邻域类型纯度（同类占比）作为诊断。**

### 1.4 对比：胰腺数据为什么成功

胰腺数据中，epsilon/schwann/mast 等稀有类型与主流类型（acinar, ductal）的表达模式差异更大（不同胚层来源、不同功能系统），在 HVG 空间中 ρ 显著小于 1 → CNEM 信号足够强。

---

## 2. 修正方案：度量空间变换

### 2.1 方案 A：判别性子空间投影（Contrastive PCA）

**核心思想**：找到一个投影矩阵 $W \in \mathbb{R}^{G \times d'}$，使得在投影空间 $W^T x$ 中，类间差异最大化而邻域内共享变异最小化。

**数学定义**：对每个细胞 $i$，定义邻域内协方差 $\Sigma_{\text{local}}(i) = \text{Cov}(\{x_j : j \in N(i)\})$ 和全局协方差 $\Sigma_{\text{global}} = \text{Cov}(\{x_j : j = 1, \ldots, N\})$。

求解对比 PCA：
$$W^* = \arg\max_W \text{tr}\left(W^T(\Sigma_{\text{global}} - \alpha \bar{\Sigma}_{\text{local}})W\right) \quad \text{s.t. } W^T W = I$$

其中 $\bar{\Sigma}_{\text{local}} = \frac{1}{N}\sum_i \Sigma_{\text{local}}(i)$，$\alpha > 0$ 为对比强度参数。

在投影空间中计算 CNEM：
$$\text{cPCA-CNEM}(i) = 1 - \cos(W^T x_i, W^T \bar{x}_{N(i)})$$

**理论保证**：在投影空间中，$\rho' = \cos(W^T\mu_r, W^T\mu_c)$ 被显式最小化（最大化类间差异），从而 CNEM 的信噪比提升。

**计算复杂度**：$O(NG^2)$ 建协方差 + $O(G^2 d')$ 特征分解，可通过随机化 SVD 加速。

### 2.2 方案 B：基因加权 CNEM

$$\text{gwCNEM}(i) = 1 - \cos(g \odot x_i, \, g \odot \bar{x}_{N(i)})$$

基因权重 $g \in \mathbb{R}^G$ 的三种选取方案：

**(B1) 方差比加权**：$g_j = \sigma_{\text{between},j} / \sigma_{\text{within},j}$，其中 between 为邻域间方差，within 为邻域内方差。放大区分邻域的基因。

**(B2) 稀疏性加权**：$g_j = 1 - p_j$，$p_j$ 为基因 $j$ 的非零表达细胞比例。稀有亚群的 marker 基因通常在大多数细胞中不表达，此加权放大这些基因。

**(B3) Laplacian 信息加权**：构建整合后 kNN 图的图 Laplacian $L$，计算前 $k$ 个非平凡特征向量 $\{v_1, \ldots, v_k\}$。定义 $g_j = \sum_{l=1}^k |v_l^T e_j|$（基因 $j$ 在图结构特征向量上的总投影），优先保留携带图结构信息的基因。

### 2.3 方案 C：RBF 核映射

$$\text{kCNEM}(i) = 1 - \frac{1}{k}\sum_{j \in N(i)} \exp\left(-\frac{\|x_i - x_j\|^2}{2\sigma^2}\right)$$

**理论保证**：RBF 核在 RKHS 中的映射 $\phi: \mathbb{R}^G \to \mathcal{H}$ 将原空间中接近的点展开到高维空间中更可分的位置。当 $\sigma$ 适中时，$\text{kCNEM}$ 对局部表达差异的灵敏度显著高于余弦。

**带宽选取**：$\sigma = \text{median}(\{\|x_i - x_j\| : (i,j) \in \text{random pairs}\})$（中位数启发式）。

### 2.4 方案 D：局部自适应马氏距离

$$\text{mCNEM}(i) = \sqrt{(x_i - \bar{x}_{N(i)})^T \Sigma_{N(i)}^{-1} (x_i - \bar{x}_{N(i)})}$$

其中 $\Sigma_{N(i)}$ 为邻域的局部协方差矩阵（需正则化以保证可逆）。

**直觉**：马氏距离自适应地放大邻域内低方差方向上的差异——如果邻域内某些基因高度一致，而细胞 $i$ 在这些基因上不同，马氏距离会放大这个信号。这正是我们需要的：在邻域的"共识基因"上检测偏差。

---

## 3. R(S) 作为 CNEM 的替代/补充

### 3.1 R(S) 的独立性分析

回顾：$R(S) = D(S)/D_{\text{ref}}$，其中 $D(S) = \frac{1}{|S|}\sum_{i \in S}\|y_i - \bar{y}_S\|^2$ 为整合后位移的弥散度。

**R(S) 与 CNEM 工作在完全不同的空间**：
- CNEM：在**表达空间**度量不匹配（语义层面）
- R(S)：在**嵌入空间**度量位移弥散（几何层面）

### 3.2 CNEM 失效时 R(S) 是否仍有效？

**是的**，在大多数情况下。原因：

CNEM 的失效条件（ρ≈1，表达空间中类间差异小）**不影响** R(S) 的有效性。R(S) 依赖的是整合前后嵌入空间中的位移模式，与表达空间中的余弦相似度无关。

具体而言：即使 pDC 与 T cell 的表达向量方向接近（ρ≈1），只要整合算法在嵌入空间中将 pDC 从原来的聚类位置移动到了 T/B cell 的区域（弥散），R(S) 就能检测到这种异常位移。

**R(S) 的失效条件**与 CNEM 不同：
- R(S) 失效当且仅当整合在嵌入空间中对稀有群做了"整体平移"（所有成员向同一方向移动相同距离）→ D(S)≈0 → 无弥散信号。
- 这对应于整合算法正确地将稀有群作为整体重新定位，不是病理性的。

### 3.3 组合策略

**推荐策略**：

$$\text{CombinedScore}(i) = \max\left(\text{rank}(\text{CNEM}_z(i)), \, \text{rank}(R_z(S_i))\right)$$

取排名的最大值（取并集），确保两种信号源中任何一个检测到的异常都不遗漏。

或者更精细：
- 在胰腺类数据（表达差异大）上：CNEM 为主，R(S) 为辅
- 在免疫类数据（表达差异小）上：R(S) 为主，改进版 CNEM（kCNEM 或 cPCA-CNEM）为辅

---

## 4. 自适应指标设计

### 4.1 核心思想

设计一个元指标，自动检测数据特性并选择/组合最优指标。

### 4.2 数据特性探针

在运行主指标之前，先计算三个"探针统计量"：

**(P1) 全局余弦集中度**：$\gamma = \text{std}(\{\cos(x_i, x_j) : (i,j) \in \text{random pairs}\})$
- $\gamma$ 小 → 余弦空间集中 → CNEM 可能失效 → 切换到 kCNEM 或 R(S)
- $\gamma$ 大 → 余弦空间有足够变异 → CNEM 有效

**(P2) 邻域类型纯度估计**：$\pi = \frac{1}{N}\sum_i \frac{|\{j \in N(i) : \text{CNEM}_j < \tau\}|}{k}$（邻域中表达匹配的比例）
- $\pi$ 高 → 大多数细胞与邻居匹配 → 正常；关注少数低 $\pi$ 的细胞
- $\pi$ 低 → 全局整合不充分 → 需要先修复整合

**(P3) 嵌入-表达空间一致性**：$\kappa = \text{Spearman}(\|z_i - z_j\|, \|x_i - x_j\|)$ for neighbor pairs
- $\kappa$ 高 → 两个空间一致 → CNEM 和 R(S) 冗余，用 CNEM 即可
- $\kappa$ 低 → 两个空间不一致 → 两个指标互补，应组合使用

### 4.3 自适应调度算法

```
输入：整合前后的细胞嵌入和表达数据
1. 计算探针 (P1), (P2), (P3)
2. if γ > γ_threshold:
     primary = CNEM
   else:
     primary = kCNEM(σ = median distance)
3. secondary = R(S)
4. if κ < κ_threshold:
     score = max(rank(primary_z), rank(secondary_z))  # 组合
   else:
     score = primary  # 冗余，用主指标即可
5. 输出 score 排名前 5% 的细胞为疑似弥散
```

### 4.4 理论保证

**自适应指标的检测力下界**：

$$\text{AUC}(\text{adaptive}) \geq \max(\text{AUC}(\text{CNEM}), \text{AUC}(R(S)), \text{AUC}(\text{kCNEM})) - \epsilon$$

其中 $\epsilon$ 为探针估计误差导致的损失。当探针样本量 $n_{\text{probe}} \geq O(\log M)$（$M$ 为类型数）时，$\epsilon \to 0$。

即：自适应指标的性能不会比任何单个子指标差太多，且在不同数据特性下自动切换到最优子指标。
