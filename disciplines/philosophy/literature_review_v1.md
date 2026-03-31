# 文献综述：哲学与跨学科视角

## 一、Unknown Unknowns 的认识论方法

### 核心框架 (Kerwin 1993; i2insights 2019)

Unknown unknowns 分为两种：
1. **错误信念型**：我们以为自己了解，实际上理解有误
2. **完全未察觉型**：我们根本不知道某个领域的存在

在本项目中，当前范式属于**类型1**——研究者以为现有评估框架能检测所有问题，实际上它对未知稀有亚群的消失完全盲视。

### 检测 Unknown Unknowns 的五大策略

| 策略 | 在本项目中的对应 |
|------|------------------|
| 谦逊 (Humility) | 承认现有评估框架有盲区 |
| 包容性 (Inclusiveness) | 多智能体协作——不同学科背景的agent互相指出盲区 |
| 严谨性 (Rigor) | 对"整合成功"的声称应用更严格的证据标准 |
| 阐释 (Explication) | 将评估框架的隐含假设显式化 |
| 接纳 (Receptiveness) | 对"可能存在未知稀有亚群"保持开放态度 |

**关键洞察**：我们8个智能体的多学科协作结构，本身就是检测 unknown unknowns 的最佳策略之一——"主动寻求背景差异大的人士意见"。

来源：[How can we know unknown unknowns?](https://i2insights.org/2019/09/10/how-can-we-know-unknown-unknowns/)

---

## 二、拓扑数据分析 (TDA) 在单细胞中的应用

### 持久同调检测稀有细胞群 (Frontiers in Immunology, 2025)

- TDA 通过追踪拓扑特征在多尺度上的出现和消失来识别稀有细胞群
- 其模型独立性和多尺度特性使其"特别适合捕捉单细胞数据内的全局组织"
- **局限**：数学复杂性高，用户友好工具少，可扩展性和可解释性存在挑战

来源：[Topological data analysis in single cell biology](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1615278/full)

### 持久同调评估批次整合效果 (Louisville 学位论文)

- 使用 Betti 曲线、欧拉特征、持久景观等拓扑指标量化整合前后的结构变化
- **关键发现：整合可能破坏局部结构，但能澄清全局拓扑**
- 这意味着：整合对稀有亚群（局部结构）的破坏可能恰好被全局指标改善所掩盖
- 景观距离矩阵在整合后仍保有生物学意义

来源：[A persistent homology framework for scRNA-Seq](https://ir.library.louisville.edu/etd/4649/)

### 对本项目的哲学启示

Louisville 论文的发现极为重要：**整合破坏局部结构但改善全局拓扑**。这正是我们问题的数学表达——稀有亚群是局部结构，它们被消灭后全局指标反而变好。这使得基于全局指标的评估天然对稀有亚群消失盲视。

**哲学诊断得到实证支持**：当前评估框架的盲区不是假设，而是已被实验证实的现象。

---

## 三、批次校正中的过度校正问题

### 最新研究 (2025)

1. **"Reference-informed evaluation with overcorrection awareness"** (Communications Biology, 2025)
   - RBET 框架：具有过度校正意识的参考信息评估
   - 但仍依赖参考标签——对未知亚群依然盲视

2. **"Toward informed batch correction"** (Nature Computational Science, 2025)
   - 提出"知情批次校正"概念
   - 但"知情"仍基于已知信息

3. **"Batch correction methods are often poorly calibrated"** (PMC, 2025)
   - 确认了校正方法普遍存在校准不良的问题

来源：
- [Reference-informed evaluation with overcorrection awareness](https://www.nature.com/articles/s42003-025-07947-7)
- [Toward informed batch correction](https://www.nature.com/articles/s43588-025-00943-1)
- [Batch correction methods poorly calibrated](https://pmc.ncbi.nlm.nih.gov/articles/PMC12315870/)

### 哲学分析

即使是2025年最新的工作，仍然在**已知标签框架内**解决过度校正问题。没有任何工作尝试在**无标签**的情况下检测过度校正对未知亚群的影响。这确认了：

> **范式盲区依然存在。我们的项目正在尝试突破这个盲区。**

---

## 四、综合哲学评估

### 文献证据支持的核心论点

1. ✅ 自封闭循环确实存在——所有最新方法仍依赖已知标签
2. ✅ TDA/持久同调是打破循环的有力工具——不依赖标签，能检测结构变化
3. ✅ 整合破坏局部结构但改善全局指标——稀有亚群消失被全局改善掩盖
4. ✅ 多学科协作是检测 unknown unknowns 的有效策略

### 尚需关注的风险

1. ⚠️ TDA 的可扩展性问题——490万细胞的规模是否可行？
2. ⚠️ TDA 的可解释性——拓扑特征消失如何映射回具体的生物学实体？
3. ⚠️ 假阳性风险——不是所有局部结构消失都意味着有意义的亚群被消灭

### 建议优先级

1. **最高**：将 TDA 路径推荐给数学和计算生物学智能体
2. **高**：要求计算生物学智能体在评估框架中同时追踪局部和全局结构变化
3. **中**：建议机器学习智能体探索无监督的异常检测方法作为 TDA 的补充
4. **待定**：湿实验选择——等计算预测出具体候选亚群后再评估
