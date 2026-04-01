# 选择性整合的数学理论框架

## 核心公式

### BD(i)：批次多样性分数

$$BD(i) = \frac{|\{\text{unique batches in } N_{\text{pre}}(i)\}|}{\min(k, B)}$$

- BD(i) ≈ 1：共享细胞（邻居来自多个批次）
- BD(i) ≈ 1/B：独有细胞（邻居几乎全来自同一批次）
- 连续值，避免硬阈值的敏感性

### 选择性整合公式

$$z_{\text{integrated}}(i) = BD(i) \cdot f_{\text{align}}(z_{\text{pre}}(i)) + (1 - BD(i)) \cdot z_{\text{pre}}(i)$$

### 与具体算法的兼容

**Harmony**：对校正向量加权
$$\Delta z(i) = BD(i) \cdot \Delta z_{\text{harmony}}(i)$$

**scVI**：对批次损失加权
$$L_{\text{total}} = \sum_i \left[L_{\text{recon}}(i) + \beta \cdot BD(i) \cdot L_{\text{batch}}(i)\right]$$

## 理论性质

### 连续性保证

若 f_align 连续（所有主流整合器满足），则 z_integrated(i) 关于 BD(i) 连续。不会出现嵌入空间断裂。

### 保护机制

独有细胞（BD≈0）的校正幅度趋近 0 → 局部邻域结构完全保留 → 稀有亚群不被消灭。

### 效率

额外开销仅为一次 kNN 查询：O(N·k·log N)。对 Harmony（O(N·k·T)）额外 <10%；对 scVI（训练时间主导）可忽略。

## 与超网络的关系

BD 加权是超网络 h(x; φ) → θ_i 的 1 维特例（输出为标量权重）。超网络可输出高维参数向量提供更精细控制，但需要端到端训练，实现复杂度高。推荐 BD 作为基线，超网络作为进阶方向。

## 与 CellANOVA 的区别

- CellANOVA：先全局整合 → 再后处理修复（后验方法，依赖整合结果）
- 我们：先分类(BD) → 再选择性整合（先验方法，仅依赖原始数据）
- 数学优势：先验方法不会被整合错误污染
