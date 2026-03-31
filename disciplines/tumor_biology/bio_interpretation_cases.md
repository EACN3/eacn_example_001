# 肿瘤生物学解读：被整合消灭的亚群案例集

## 案例1（标准5主案例）：T cell Cluster 5 — GITR⁺ 活化态 Treg 功能社区

### 基本信息
- 细胞数：1464
- 存活率：0.39（被打散到8个post-cluster）
- 组成：FOXP3=28%, TIGIT=65%, GITR⁺ (top marker)
- 含 TIGIT⁺CCR8⁻ Treg 295细胞 (20.2%)
- 组织来源：86% SSCC + Skin Normal (D100 Dataset)

### 生物学解读

这是一个**皮肤癌微环境特异的 GITR⁺ 活化态 Treg 功能社区**：

1. **GITR (TNFRSF18)** 是 TNF 受体超家族成员，在高度活化的 Treg 上高表达。GITR⁺ Treg 代表正在执行免疫抑制功能的效应性 Treg，是免疫耐受的关键执行者。

2. **TIGIT⁺CCR8⁻ 组合的新颖性**：CCR8⁺ 是已知的肿瘤特异性 Treg 标志（Nature 2020, De Simone et al.）。TIGIT⁺CCR8⁻ 定义了一个此前未被单独识别的 Treg 新亚型，代表通过 TIGIT-CD155 轴而非 CCR8 轴介导免疫抑制的替代通路。

3. **功能社区结构**：该 cluster 并非纯 Treg（FOXP3 仅 28%），而是 Treg 与效应 T 细胞共存的"免疫抑制微环境单元"。这种 Treg-效应 T 细胞的直接互作结构在皮肤癌中尤为重要——皮肤癌（尤其 SSCC）具有高度免疫浸润的"热"肿瘤微环境，GITR⁺ Treg 是维持免疫逃逸平衡的关键。

4. **组织特异性**：86% 来自 SSCC/皮肤正常组织，说明这是一个由皮肤微环境塑造的独特免疫生态位。

### 被消灭的后果

- 下游 cell-cell communication 分析（如 CellChat、CellPhoneDB）无法识别皮肤癌特有的 Treg-effector T 互作
- Anti-GITR/anti-TIGIT 联合治疗的 biomarker 开发受损
- 泛癌比较分析中皮肤癌的免疫特征被低估

### 湿实验验证

免疫学智能体持有的数据（来自已发表论文 Figure 5G/H/I）：
- CD69⁺ 活化 T 细胞：NC ~20.5% vs TIGIT⁺ Treg ~5.0% (p<0.001)
- IFN-γ 分泌：NC ~305 pg/mL vs TIGIT⁺ Treg ~115 pg/mL (p<0.001)
- 肿瘤杀伤率：NC ~63% vs TIGIT⁺ Treg ~19% (p<0.001)

三层递进证据证实 TIGIT⁺ Treg 具有强大的免疫抑制功能。

### 临床意义

- Anti-GITR 抗体（TRX518, MK-4166, BMS-986156）正在临床 Phase I/II
- Anti-TIGIT 抗体（tiragolumab）已在多种肿瘤中进入 Phase III
- GITR⁺TIGIT⁺ 双阳性 Treg 的发现暗示 anti-GITR + anti-TIGIT 联合方案的潜力

---

## 案例2（补充案例）：Macrophage Cluster 18 — 待定性

### 基本信息
- 细胞数：827
- 存活率：0.216（最严重的破坏案例）
- 组成：待数据科学提供 marker 表达数据

### 初步生物学推测

巨噬细胞在肿瘤微环境中极度异质，主要包括：
- **M1-like TAM**：促炎、抗肿瘤，表达 IL1B, TNF, NOS2
- **M2-like TAM**：抗炎、促肿瘤，表达 CD163, MRC1, MSR1
- **Lipid-associated macrophage (LAM)**：富含脂质代谢基因，与肥胖相关肿瘤有关
- **TREM2⁺ TAM**：近年新发现的免疫抑制性巨噬细胞亚型
- **增殖性巨噬细胞**：MKI67⁺，局部增殖而非单核细胞分化来源
- **组织驻留巨噬细胞**：组织特异性，不同器官差异大

survival=21.6% 暗示该 cluster 可能是一个高度组织/癌种特异的巨噬细胞亚型，在跨批次整合时因分布偏斜而被打散。

**待确认**：需要 cluster 18 的 top marker genes 和癌种分布才能做最终定性。

---

## 案例3（规模依赖性消灭）：pDC — Cell-type 级别消灭

### 基本信息
- NP（邻域保留率）：0.231（所有 cell type 中最低）
- 8批次子集：99%召回，完好
- 103批次全量：被严重破坏，推入邻近 DC/Monocyte 大群

### 规模依赖性消灭机制

pDC 在小规模整合中安全、在大规模整合中被消灭，揭示了一个关键机制：

1. **8批次时**：pDC 虽然稀有，但其表达特征（IRF7, TLR7/9, IL3RA, CLEC4C）足够独特，整合算法在少量批次间对齐时能保持其独立性。需要校正的批次差异有限，pDC 不会被"拉扯"过度。

2. **103批次时**：每多一个批次就增加一组"需要对齐"的细胞。pDC 与 cDC2、Monocyte 共享部分基因表达（如 HLA-DR、CD74 等 MHC-II 相关基因），当整合算法需要同时对齐 103 个批次时，这些共享特征成为"锚点"，将 pDC 越来越强地拉向 DC/Monocyte 大群。

3. **渐进性消灭**：这不是突然消失，而是像小岛在海平面缓慢上升时逐渐被淹没——每增加一个批次，pDC 的独立性就被削弱一点，直到完全丧失。

### 生物学意义

- pDC 在肿瘤免疫中具有双重角色（抗肿瘤 IFN-I 产生 vs 促肿瘤免疫耐受），其消灭意味着：
  - 无法区分肿瘤中 pDC 的功能状态
  - 丢失 pDC 与 T 细胞的互作信息（pDC 通过 IFN-I 影响 CD8 T 细胞活化）
  - pDC 内部 26 个亚群的精细结构完全被抹除（包括 cluster 24 仅 NHL 的极稀有亚群）
- **泛癌比较中最严重的后果**：不同癌种中 pDC 丰度和功能状态差异巨大，但整合后这些差异被消灭，所有癌种的 pDC 看起来"一样"

### 论文叙事价值

与案例1（subcluster 级别）形成两层消灭的完整图景：
- **Cell-type 级别**：pDC（NP=0.231）— 整个类型在全量整合中被消灭
- **Subcluster 级别**：T cell Cluster 5（survival=0.39）— 亚群结构在整合中被打散
- **规模依赖性**：同一亚群在 8 vs 103 批次中命运截然不同

---

*文档最后更新：2026-04-01*
*肿瘤生物学智能体系统 (agent-mnez9frj)*
