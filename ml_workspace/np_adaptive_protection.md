# NP驱动的自适应保护策略设计（标准2）

> 作者：机器学习智能体系统 (agent-mnez8qvx)
> 任务编号：t-mnf582gc
> 基础：TCI框架 + NP检测成功（分层AUC=0.841，全局AUC=0.837）
>
> **实现说明**：设计文档描述了两种实现模式——(A) per-cell L_batch权重调节（理论最优，需改scVI源码）和 (B) batch标签打乱（实际采用，不需改源码）。代码 np_guard_implementation.py 和实验均使用方案B。

---

## 1. 核心思想

**预防原则的算法化**：在整合过程中实时监控每个细胞的邻域保留率（NP），当某区域 NP 骤降时，自动降低该区域的整合力度。优先保护结构，其次追求对齐。

## 2. 算法框架：NP-Guard

### 2.1 概览

```
输入：多批次单细胞数据 {X_b}，基础整合器 M（如scVI）
输出：保护性整合结果

NP-Guard 在基础整合器的训练循环中插入三个模块：
  [A] NP 监控器（每 N 步计算一次局部 NP）
  [B] 风险评估器（识别 NP 骤降的区域）
  [C] 自适应调节器（对高风险区域降低整合强度）
```

### 2.2 模块 A：NP 监控器

```
每 T_monitor 个训练步（如 T_monitor=50）：
1. 从当前模型获取所有细胞的嵌入 Z_current
2. 构建 kNN 图 G_current（k=30）
3. 对每个细胞 i，计算：
   NP(i) = |N_pre(i) ∩ N_current(i)| / k
   其中 N_pre(i) 是整合前表达空间的 k 近邻（只计算一次，缓存）
4. 输出 NP 向量：np = [NP(1), NP(2), ..., NP(n)]
```

**计算效率**：N_pre 只算一次（O(nk log n)），之后每次监控只需重建 G_current 和计算交集。对于 scVI 每 epoch 数十步的训练循环，每 50 步监控一次开销可接受。

### 2.3 模块 B：风险评估器

```
对 NP 向量进行局部聚合：
1. 使用整合前的 kNN 图 G_pre 定义局部邻域
2. 对每个细胞 i，计算邻域平均 NP：
   NP_local(i) = mean(NP(j) for j in N_pre(i))
3. 定义风险分数（相对变化版，采纳肿瘤生物学建议）：
   risk(i) = max(0, 1 - NP_local_post(i) / NP_local_pre(i))
   即检测NP的相对下降幅度，而非绝对值。
   这样本来就分散的稀有类型（如cDC3，NP_pre本身低）不会被永久保护。
   NP_pre(i)→NP_post(i) 变化小则risk≈0，变化大则risk→1。
4. 归一化：risk_norm(i) = risk(i) / max(risk) ∈ [0,1]
```

**为什么用局部聚合**：单个细胞的 NP 噪声大，局部平均更稳健。如果一个区域（而非单个点）的 NP 同时下降，说明是系统性结构破坏，不是随机波动。

### 2.4 模块 C：自适应调节器

与 scVI 的兼容方式——**per-cell 损失权重调节**：

```
scVI 的标准损失函数：
  L = Σ_i L_recon(i) + β * L_kl(i) + γ * L_batch(i)

NP-Guard 修改为：
  L = Σ_i [L_recon(i) + β * L_kl(i) + γ * α(i) * L_batch(i)]

其中 α(i) = 1 - λ * risk_norm(i)

  - risk_norm(i) = 0（安全区域）→ α(i) = 1 → 正常整合
  - risk_norm(i) = 1（高风险区域）→ α(i) = 1-λ → 整合力度降低
  - λ ∈ [0,1] 控制最大降低幅度（默认 λ=0.8，即高风险区域整合力度降至 20%）
```

**直觉**：L_batch 是驱动跨批次对齐的损失项。对高风险区域降低 L_batch 权重，等于告诉模型"这个区域不需要那么用力对齐"——保护了局部结构。

### 2.5 与其他可微整合器的兼容

| 整合器 | 批次对齐损失项 | NP-Guard 调节方式 |
|--------|---------------|-------------------|
| scVI | L_batch（条件VAE的批次条件项） | 降低条件强度 |
| scANVI | L_classify + L_batch | 只调节 L_batch，保留 L_classify |
| Harmony | 迭代软聚类中的批次校正步 | 对高风险细胞减小校正向量的幅度 |
| scDREAMER | 对抗损失 L_adv | 降低对抗损失权重 |

### 2.6 非可微整合器的适配（Scanorama, BBKNN）

对于非迭代的方法，NP-Guard 改为**后处理保护模式**：

```
1. 执行标准整合
2. 计算 NP
3. 对 NP < NP_threshold 的细胞，将其嵌入拉回整合前位置：
   Z_protected(i) = (1 - risk_norm(i)) * Z_post(i) + risk_norm(i) * Z_pre_aligned(i)
   其中 Z_pre_aligned 是将整合前嵌入通过全局刚体变换对齐到整合后空间
4. 重新计算 NP，验证保护效果
```

## 3. 超参数

| 参数 | 默认值 | 含义 |
|------|--------|------|
| k | 30 | kNN 的 k 值 |
| T_monitor | 50 | 每多少训练步监控一次 |
| NP_threshold | 0.3 | NP 低于此值视为高风险 |
| λ | 0.8 | 最大整合力度降低幅度 |

**NP_threshold 的设置依据**：
- 在无整合时，同一类型细胞跨批次的 NP 通常在 0.5-0.8（因为批次效应导致部分邻居变化）
- 整合后，正常校正的细胞 NP 应提高到 0.6-0.9
- NP < 0.3 意味着 70% 以上的邻居都变了，高概率是结构破坏
- 可用计算验证数据（T cell cluster 5 vs cluster 4 的 NP 分布）校准

## 4. 理论保证

### 4.1 不牺牲主流细胞类型的批次校正

NP-Guard 只对高风险区域降低整合力度。主流细胞类型在整合时 NP 保持高（邻居不变），risk_norm ≈ 0，α(i) ≈ 1——整合行为与未加保护时完全相同。

**数学表述**：设主流类型集合为 C_major，稀有类型集合为 C_rare。如果 ∀i ∈ C_major: NP(i) > NP_threshold，则 risk_norm(i) = 0，α(i) = 1，L_batch(i) 不被调节。主流类型的整合质量不受影响。

### 4.2 稀有亚群保护的上界

最好情况：NP-Guard 将高风险细胞的整合力度降至 (1-λ)，这些细胞的嵌入保持接近整合前位置，亚群结构完全保留。

最差情况：如果稀有亚群的破坏在第一次监控（T_monitor 步）之前就发生了，NP-Guard 来不及保护。解决方案：T_monitor 设小（如 10），或在训练早期整合力度全局降低（warmup）。

## 5. 可扩展性（标准3）

对于 490 万细胞：
- N_pre 的 kNN 图：用 FAISS 或 Annoy 近似最近邻，O(n log n)，490 万细胞约 5 分钟
- 每次 NP 监控：重建 G_current 的 kNN（O(n log n)） + 集合交集（O(nk)）
- 总额外开销：约为基础整合器训练时间的 10-20%

## 6. 生物学安全规则（采纳肿瘤生物学智能体建议）

NP-Guard 的预处理模块，防止保护伪群：

1. **Doublet 过滤**：在 NP-Guard 之前用 Scrublet/DoubletFinder 清除 doublet。Doublet 形成"假稀有群"，NP 低但不应保护。
2. **技术伪群排除**：如果某群的 top markers 被 HBA/HBB/MALAT1 等 ambient RNA 污染基因主导，排除保护。
3. **应激群降级**：如果某群的 top markers 被解离诱导的应激基因（FOS/JUN/HSP90 等 immediate-early genes）主导，降低保护优先级。

## 7. 验证计划

1. 在 T cell 数据上：比较 scVI vs scVI+NP-Guard 在 cluster 5 上的 survival_rate
2. 预期：NP-Guard 应显著提高 cluster 5 的 survival_rate（从 ~0.4 到 >0.7）
3. 同时验证主流类型（T cell cluster 4, survival=0.93）的整合质量不下降
4. 如果成功，标准2达成
