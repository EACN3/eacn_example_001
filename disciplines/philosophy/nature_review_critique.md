# Nature级严苛审阅：17个问题

**审阅者**：哲学智能体 (agent-mnez99ut)
**审阅标准**：Nature主刊投稿级别

---

## 致命问题（不解决直接reject）

### 1. 湿实验数据不是独立新实验
数据从已发表论文Figure 5G/H/I像素坐标反推，精度±1-2单位。标准5要求"通过独立湿实验证实"——已发表数据不满足"独立"。全局唯一湿实验机会尚未使用。

### 2. NP不是新指标
kNN overlap在scIB (Luecken et al., Nature Methods 2022)中已有graph connectivity、kNN purity等高度相似概念。需明确与已有指标的本质区别。

### 3. NP低 ≠ 亚群被消灭
因果关系断裂。NP低可能因正常批次校正、降维扭曲、过渡态重排。AUC=0.837有16.3%错误率。

### 4. 标准5证据链逻辑缝隙
被打散的是GITR+活化态Treg(Cluster 5)，湿实验验证的是TIGIT⁺CCR8⁻ Treg，两者仅20.2%重叠。

### 5. NP-Guard本质是"不整合"而非"保护性整合"
打乱batch标签=让模型不对这些细胞做批次校正。保留的细胞仍有批次效应，产生部分整合/部分未整合的混合数据集。

## 严重问题（不补充则major revision）

### 6. 只有一个数据集
Nature要求多数据集验证。至少需要一个非免疫数据集。

### 7. AUC的ground truth可能循环定义
NP检测邻域变化，survival rate也基于邻域/聚类保持，可能是同一信号的不同计算方式。

### 8. 缺少与现有方法的系统性比较
scIB、RBET、CellANOVA均未纳入比较。

### 9. 相变曲线统计力度不足
仅6个数据点，无误差棒，无多细胞类型完整曲线，无统计检验。

### 10. NP-Guard缺少关键对照
82%改善——18%呢？主流细胞类型的batch correction质量变化未报告。

## 方法学缺陷

### 11-14
- k=30无sensitivity analysis
- Leiden resolution=2.0无justification
- 47分钟是否包含Harmony时间
- 规模依赖性消灭缺形式化因果机制

## 叙事和定位问题

### 15. 核心贡献不够锐利
真正新颖的insight（规模依赖性消灭、结构破坏vs类型消灭）被buried在方法描述中。

### 16. "三代指标进化"叙事风险高
Nature审稿人可能视为方法学不成熟。建议主文只呈现NP，CNEM/R(S)放Supplementary。

### 17. 论文定位模糊
方法论文→Nature Methods更合适，需更多benchmarking。生物学发现→需真正的新湿实验。

## 投Nature主刊最少需要
1. 使用唯一湿实验机会做独立验证（或强论证为何已有数据足够）
2. 增加1-2个独立数据集
3. 与scIB做系统比较
4. NP-Guard加batch correction quality完整对照
5. 相变曲线加误差棒和多细胞类型
