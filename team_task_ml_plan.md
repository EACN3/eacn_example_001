# 团队核心任务 — 机器学习智能体工作计划

> 任务ID：t-mnf2lkh2
> 团队：team-mnezv3cg

---

## ML 智能体的职责范围

从 ML 文献和算法思想层面贡献：
- 标准1（能检测）：设计不依赖标签的检测框架
- 标准2（能保护）：设计拓扑守恒整合策略
- 标准3（能扩展）：确保框架在490万细胞规模可行

## 已完成工作

1. **TCI 框架设计**（ml_algorithm_survey.md）— 三层架构：拓扑指纹提取→拓扑正则化整合→弥散检测器
2. **DS 弥散分数**：DS(i) = 2*(1-NP) + max(0,ΔC) + max(0,-ΔK)，权重经数学智能体 Fisher 分析验证
3. **CNEM 失效诊断**（cnem_failure_analysis_and_fixes.md）— 5种改进方案
4. **关键发现**：pDC 未被消灭，CNEM 低分是 true negative，需要找真正被消灭的目标

## 当前进行中

- **子任务 t-mnf2tdww**：委托计算生物学执行人工消灭实验 + 系统性扫描真正被消灭的亚群
  - 实验A：cluster 24（93细胞）的 DS/CNEM 验证
  - 实验B：系统性扫描整合中被打散的小簇

## 关键进展（来自数据科学+肿瘤生物学交叉分析）

- T cell cluster 7（4095细胞，survival=0.29）是真正被消灭的亚群，完美阳性对照
- T cell cluster 4（survival=0.93）是完美阴性对照
- TIGIT+CCR8- Treg（1922细胞）未被消灭（R(S)=0.73）→ 细胞数足够+表达独特=安全
- 肿瘤生物学6个高危候选中，cDC3/AS-DC（<0.5%，边界模糊）和 pDC cluster 24 最可能被消灭

## 验证策略（两层）

- **计算验证**：T cell cluster 7 vs cluster 4 → R(S)+CNEM 的 ROC-AUC
- **生物学验证**：cDC3/AS-DC 或 pDC cluster 24 → 湿实验靶标

## 重要转折

- **Cluster 7 是假阳性**：Naive/Tcm，75.4%实际仍聚一起，survival=0.29因混入其他细胞
- **Cluster 5 是真靶标**：GITR+ TIGIT+CCR8- 活化态 Treg（295细胞），被打散到8个post-cluster
- **R(S) 失败**：与亚群消灭率强负相关(r=-0.952)，检测的是嵌入校正强度而非消灭
- **NP 成功！**：分层NP（per-celltype内邻域保留率）AUC=0.841, Spearman r=0.613
- NP 是 DS 弥散分数中权重最高的分量（2:1:1），也是最简单最有效的单指标

## 下一步计划

1. 等待计算生物学测试 cohesion_retention 在所有 pre-cluster 上的表现
2. 如果 cohesion_retention 能区分 cluster 5（被打散）vs cluster 4（完好）→ 标准1达成
3. Cluster 5 的 TIGIT+CCR8- Treg 与湿实验数据直接对应 → 标准5可达成
4. 完善 TCI 框架：用 cohesion_retention 替代 R(S) 作为核心检测指标
5. 最终输出：Nature 论文 Methods 部分
