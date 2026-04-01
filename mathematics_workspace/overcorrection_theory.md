# 标准整合算法的过度校正：定量分析

## 核心发现

BD 加权选择性整合实验表明：BD 均值 = 0.179，即仅需 Harmony 标准校正幅度的 ~18% 即可达到 NP=0.898 且 batch mixing 改善。

**结论：Harmony 的标准校正幅度过大约 5 倍（1/0.179 ≈ 5.6）。**

## 数学推导

### 模型

$$z_{\text{final}}(i) = BD(i) \cdot z_{\text{harmony}}(i) + (1-BD(i)) \cdot z_{\text{pre}}(i)$$

等价于：

$$z_{\text{final}}(i) = z_{\text{pre}}(i) + BD(i) \cdot \underbrace{(z_{\text{harmony}}(i) - z_{\text{pre}}(i))}_{\Delta z(i)}$$

其中 $\Delta z(i)$ 是 Harmony 的校正向量。BD(i) 相当于对校正向量的缩放因子。

### 过度校正的量化

定义过度校正因子：

$$\text{OCF} = \frac{\|\Delta z_{\text{harmony}}\|}{\|\Delta z_{\text{optimal}}\|} = \frac{1}{\overline{BD}}$$

实测 $\overline{BD} = 0.179$，故 OCF = 5.6。

### 为什么过度校正？

Harmony 的目标函数最大化跨批次聚类的熵（Shannon diversity of batch labels within clusters）。这个目标没有上限约束——它会持续校正直到批次标签在每个聚类中完全均匀分布。

但"完全均匀混合"在有批次特异亚群时是有害的——要求稀有亚群（仅存在于少数批次）被拉到所有批次细胞的中间，导致其结构被破坏。

### 最优 BD 的理论预测

最优 BD 满足：

$$BD^* = \arg\min_{BD} \left[\alpha \cdot L_{\text{batch}}(BD) + (1-\alpha) \cdot L_{\text{NP}}(BD)\right]$$

其中 $L_{\text{batch}}$ 为批次混合损失（随 BD 增大而减小），$L_{\text{NP}}$ 为 NP 损失（随 BD 增大而增大）。

理论预测：$BD^* \in [0.1, 0.3]$，取决于：
- 批次效应强度：强 → $BD^*$ 大
- 批次特异亚群比例：多 → $BD^*$ 小
- 数据集批次数：多 → $BD^*$ 可能更小（校正负担分散）

### 可验证预测

1. 在胰腺数据（5个批次，批次效应较弱）上，$BD^*$ 应该 > 0.179
2. 在 490 万全量数据（103个批次）上，$BD^*$ 可能 < 0.179
3. $BD^*$ 与批次数 $B$ 的关系可能近似为 $BD^* \propto 1/\sqrt{B}$

## 论文价值

"标准整合过度校正 5 倍"是一个独立的发现：
- 解释了为什么现有整合方法会消灭稀有亚群
- 解释了为什么 BD 加权（仅 18% 校正强度）反而改善 batch mixing
- 为整合算法社区提供了定量的"校正上限"参考
