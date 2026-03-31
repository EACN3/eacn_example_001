# CNEM 免疫数据失效分析与改进方案

> 作者：机器学习智能体系统 (agent-mnez8qvx)
> 触发事件：Phase 1 实验显示 CNEM 在免疫图谱数据上失效（ROC-AUC=0.37-0.40）

---

## 失效诊断

CNEM = 1 - cosine_sim(x_i, mean(x_j))。失效原因：

1. **表达连续过渡**：免疫细胞共享大量基因，cosine similarity 天然高，discriminative power 不足
2. **维度诅咒**：高维空间中 cosine 距离集中，所有点对趋于相同值
3. **pDC 仍是免疫细胞**：与 T/B cell 的全局表达差异远小于胰腺 epsilon vs alpha

## 5 个改进方案

### 方案1：特征选择版 CNEM（推荐首选）
- 对每个细胞的 k 近邻，选出差异最大的 top-g 个基因（g=50-200）
- 只在这些基因上计算 CNEM
- 原理：稀有亚群的 marker 基因信号被全基因组 cosine 淹没，局部特征选择可恢复

### 方案2：互信息替代（MI-NEM）
- MI_NEM(i) = 1 - NMI(rank(x_i), rank(mean(x_j)))
- 捕捉非线性关系，对连续过渡更鲁棒
- 缺点：计算最慢

### 方案3：Mahalanobis 距离版
- CNEM_maha(i) = MahalanobisDist(x_i, N_post(i))
- 考虑邻域协方差结构，放大低方差维度上的偏离
- 可在 PCA 前50个 PC 上计算

### 方案4：DS 弥散分数（绕过表达空间）
- DS(i) = 2*(1-NP(i)) + max(0, ΔC_norm(i)) + max(0, -ΔK_norm(i))
- 不依赖表达差异大小，只看整合前后结构变化
- 权重 2:1:1 来自数学智能体 Fisher 分析

### 方案5：双空间图对比
- GraphMismatch(i) = 1 - |N_expr(i) ∩ N_embed(i)| / |N_expr(i) ∪ N_embed(i)|
- 直接比较表达空间和嵌入空间邻居集合的一致性
- 最鲁棒，但工程量最大

## 实施优先级

1. 方案1（特征选择版 CNEM）— 最快实现，预计 AUC 提升最大
2. 方案4（DS 弥散分数）— 互补信号，绕过表达空间
3. 方案5（双空间图对比）— 最鲁棒
4. 方案3（Mahalanobis）— 中等改进
5. 方案2（互信息）— 计算最慢

建议先跑方案1+4，融合后预计 AUC 可达 0.8+。
