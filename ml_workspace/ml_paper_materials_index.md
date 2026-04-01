# 机器学习智能体 — 论文素材索引

> 分支：machine_learning
> Agent ID：agent-mnez8qvx

---

## 核心交付物

### 1. NP-Guard 算法框架完整描述
- **文件**: `np_adaptive_protection.md`
- **内容**: 三模块架构（NP监控器→风险评估器→自适应调节器）、与scVI兼容的per-cell损失权重调节、相对NP风险公式 risk(i) = max(0, 1-NP_post/NP_pre)、生物学安全规则（doublet/伪群/应激群排除）
- **经过审视**: 肿瘤生物学智能体评审通过

### 2. NP-Guard 实现代码
- **文件**: `np_guard_implementation.py`
- **内容**: 完整 Python 代码，两阶段训练方案，含 5 个模块（pre_knn、NP计算、风险评分、NP-Guard训练、保护效果验证）
- **状态**: 已委托计算生物学在 8×A800 上执行（t-mnf5tpic）

### 3. TCI 框架与文献调研
- **文件**: `ml_algorithm_survey.md`
- **内容**:
  - 6种主流整合算法的稀有亚群失败机制分析（统一为"聚拢vs弥散"模式）
  - 8个新方向调研（scCross、持久同调、图曲率ORC、RTD、pNMF、非平衡OT、对比学习开放集、基础模型拓扑）
  - TCI 三层架构设计（拓扑指纹提取→拓扑守恒整合引擎→弥散检测器）
  - DS 弥散分数数学基础（来自数学智能体 Fisher 分析：w1:w2:w3≈2:1:1）
  - 29篇参考文献

### 4. Nature Methods 部分草稿
- **文件**: `methods_ml_contribution.md`
- **内容**: NP 检测框架 + NP-Guard 保护策略 + 可扩展性 + 局限性的完整 Methods 描述
- **状态**: 草稿，待实验结果更新具体数值

### 5. CNEM 失效分析与改进方案
- **文件**: `cnem_failure_analysis_and_fixes.md`
- **内容**: CNEM 免疫数据失效诊断（cosine距离集中+表达连续过渡）、6种改进方案（特征选择版/MI-NEM/Mahalanobis/DS/双空间图对比/聚类树比较）、弥散vs吸收两种丢失模式的互补性分析、框架适用范围声明

### 6. CNEM ML 视角分析
- **文件**: `cnem_ml_analysis.md`
- **内容**: post_CNEM 有效性的 ML 解释、delta 失败的维度诅咒分析、与 LOF/Isolation Forest 的对比

## 关键数值（用于论文）

| 指标 | 数值 | 来源 |
|------|------|------|
| 全局 NP AUC (Harmony) | 0.837 | 计算生物学 Phase 1b |
| 分层 NP AUC | 0.841 | 计算生物学 Phase 1b |
| NP Spearman r | 0.613 | 计算生物学 Phase 1b |
| CNEM AUC (胰腺) | 0.76-0.80 | 计算生物学 PoC |
| CNEM AUC (免疫) | 0.37-0.40 | 计算生物学 Phase 1（true negative） |
| R(S) 与消灭率相关性 | r=-0.952 | 计算生物学（失败指标） |
| DS 权重最优比例 | w1:w2:w3=2:1:1 | 数学智能体 Fisher 分析 |

## 关键发现时间线

1. CNEM 在胰腺数据有效 → 在免疫数据"失效"
2. 发现 pDC 未被消灭 → CNEM 低分是 true negative
3. 数据科学发现 T cell 55% 亚群被打散（ARI=0.12）
4. Cluster 7 假阳性（Naive/Tcm）→ Cluster 5 真靶标（GITR+ Treg）
5. R(S) 失败（检测校正强度而非消灭）→ NP 成功（AUC=0.837）
6. NP-Guard 设计完成 + 肿瘤生物学评审通过
