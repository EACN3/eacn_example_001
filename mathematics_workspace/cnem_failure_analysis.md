# CNEM 免疫图谱失效分析与改进方案

## 实验事实

CNEM 在免疫图谱数据上失效：ROC-AUC=0.37-0.40，pDC 的 post_cnem 反而低于 T cell/B cell。

## 根因

失效条件 F1 被验证：免疫细胞在 HVG 空间中余弦相似度过高（ρ≈1），加上高维集中现象（concentration of measure），CNEM 的区分力被消灭。

## 五个改进方向

1. **判别性子空间投影**：contrastive PCA，去掉邻域内共享变异，保留类间差异
2. **基因加权 CNEM**：变异系数加权 / 稀疏基因加权 / Laplacian 特征向量筛选
3. **非线性核映射**：RBF 核展开高维空间，带宽 σ = median pairwise distance
4. **欧氏距离替换余弦**：保留模长信息
5. **局部 MMD**：分布级比较，有统计检验保证

## 关键诊断问题

pDC 的 post_cnem 低可能说明整合成功保留了 pDC（邻居仍是同类）。需检查 pDC 整合后邻居的类型组成。
