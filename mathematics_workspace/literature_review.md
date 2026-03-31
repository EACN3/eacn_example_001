# 数学智能体 — 文献检索结果

## 一、最优传输 (OT) 在单细胞中的应用

### 已有工作（OT 用于轨迹推断与扰动预测，非批次整合检测）

| 论文 | 年份 | 核心思路 | 与我们问题的关系 |
|------|------|----------|------------------|
| Waddington-OT (Schiebinger et al.) | 2019 | 用 OT 在时间切片间匹配细胞，推断分化轨迹 | 提供了 OT 在单细胞间建立对应关系的范式 |
| CellOT (W₂ neural OT) | 2024 | 用 Wasserstein-2 OT 预测单细胞扰动响应 | 展示了 OT 映射可以捕捉细胞分布变化 |
| W₁ Neural OT (2411.00614) | 2024 | 用 W₁ 对偶简化为单函数最大化，25-45x 加速 | 可扩展到大规模数据集 |
| TrajectoryNet (2002.04461) | 2020 | 连续正则化流 + 动态 OT 建模细胞动力学 | 连续路径建模，可追踪细胞位移 |
| Geodesic Sinkhorn (2211.00805) | 2022 | 流形上的 Sinkhorn，O(n log n) 复杂度 | 高效计算流形测地距离上的 OT |
| MIOFlow (2206.14928) | 2022 | 流形插值 OT 流，处理分叉与合并 | 可建模群的分裂/合并过程 |
| LOT for single-cell (2510.22033) | 2025 | 线性 OT 嵌入，将点云映射到固定维向量空间 | 提供了可解释的 OT 表示 |
| Unbalanced Neural OT (2506.11969) | 2025 | 非平衡 OT + Fréchet 回归，允许质量变化 | **关键**：非平衡 OT 自然处理细胞数不守恒 |
| WFR-MFM (2601.20606) | 2026 | Wasserstein-Fisher-Rao 几何下的非平衡 OT | **关键**：WFR 度量同时量化传输和质量增减 |

### 关键发现

**没有任何现有工作将 OT 用于检测批次整合是否消灭了稀有亚群。** OT 目前被用于：
- 轨迹推断（时间点间匹配）
- 扰动响应预测
- 空间转录组重建
- 分布插值

**我们提出的用途是全新的：用 OT transport plan 的结构来区分"聚拢"与"弥散"。**

## 二、持久同调 (TDA/PH) 在单细胞中的应用

| 论文 | 年份 | 核心思路 | 与我们问题的关系 |
|------|------|----------|------------------|
| Hernández-Lemus, "TDA in single cell biology" | 2025 | 综述：PH 和 Mapper 在单细胞中检测稀有群、过渡态、分支轨迹 | **直接相关**：明确提到 TDA 可检测 rare cell populations |
| Daneshmand et al., "Impact of integration on PH clustering" | 2025 | 用 PH 评估批次整合对拓扑特征的影响，比较 Betti 曲线、Euler 特征、persistence landscapes | **高度相关**：首次将 PH 用于评估整合质量 |
| Bokor Bleile et al., "Persistence diagrams as morphological signatures" | 2023 | 用 PH 量化单细胞形态，2-Wasserstein 距离比较 | PH + Wasserstein 距离的组合方法 |

### 关键发现

**Daneshmand et al. (2025) 是最接近我们问题的工作**——他们用 PH 来评估整合是否破坏了生物学信号。但他们的工作：
- 聚焦于整体拓扑结构变化，未专门针对稀有亚群
- 使用已知标签来验证，仍在"已知"的框架内
- 未提出针对"未知稀有亚群消灭"的检测指标

**Hernández-Lemus (2025)** 明确指出 TDA 能检测 rare cell populations，但这是在非整合场景下。

## 三、关键差距（我们的切入点）

现有文献中存在一个明确的空白：

1. **OT 领域**：大量工作将 OT 用于单细胞分布间的映射，但**没有人用 transport plan 的结构模式来诊断整合病理**（聚拢 vs 弥散）。
2. **TDA 领域**：已有人用 PH 评估整合质量，但**没有人专门设计检测未知稀有亚群消灭的拓扑指标**。
3. **OT + TDA 组合**：**完全未被探索**。

这正是数学智能体应该贡献的方向。
