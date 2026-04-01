# 机器学习智能体完整工作日志

> Agent ID: agent-mnez8qvx
> 分支: machine_learning
> 注册时间: 2026-04-01
> 总commit数: 22

---

## 一、启动与环境准备

1. 克隆仓库到工作区 `E:/eacn_example/CC_1/eacn_example_001`
2. 连接eacn3网络 `http://175.102.130.69:37892`
3. 注册为机器学习智能体（agent-mnez8qvx），domains: machine-learning, single-cell-integration, batch-correction等
4. 设置3分钟定时轮询next
5. 创建`machine_learning`分支
6. 完成团队握手（team-mnezv3cg，生物科学智能体发起）

## 二、第一阶段：文献调研与框架设计

### 2.1 ML算法调研（任务 t-mnf028tw）
- 分析6种主流整合算法（Harmony/scVI/Scanorama/BBKNN/scDML）的稀有亚群失败机制
- 统一归纳为"聚拢vs弥散"模式
- 调研8个新方向：scCross、持久同调、图曲率ORC、RTD、pNMF、非平衡OT、对比学习开放集、基础模型拓扑
- 提出TCI（拓扑守恒整合）三层架构
- 29篇参考文献
- 输出：`ml_algorithm_survey.md`

### 2.2 DS弥散分数权重优化（委托数学智能体，t-mnf19n46）
- 数学智能体完成Fisher判别分析：w1:w2:w3≈2:1:1
- DS与CNEM互补大于冗余
- FusedScore Fisher最优解

### 2.3 CNEM分析
- ML视角分析post_CNEM有效性（双空间不一致检测）
- delta失败的维度诅咒解释
- 与LOF/Isolation Forest对比
- 输出：`cnem_ml_analysis.md`

## 三、第二阶段：CNEM失效与NP发现

### 3.1 CNEM免疫数据失效诊断
- Phase 1显示CNEM在免疫数据AUC=0.37-0.40
- 诊断：cosine距离集中+表达连续过渡
- 提出6种改进方案（特征选择版/MI-NEM/Mahalanobis/DS/双空间图对比/聚类树比较）
- 发现pDC未被消灭→CNEM低分是true negative
- 输出：`cnem_failure_analysis_and_fixes.md`

### 3.2 弥散vs吸收两种丢失模式
- 采纳哲学智能体审视：DS专检弥散但漏检吸收
- CNEM可检测两种模式
- FusedScore改为max策略
- 聚类树比较（Robinson-Foulds距离）作为吸收检测备选

### 3.3 NP成功（AUC=0.837）
- 建议NP（邻域保留率）作为DS的简化版
- 全局NP AUC=0.837，分层NP AUC=0.841
- R(S)失败（r=-0.952，检测校正强度而非消灭）
- NP确立为核心检测指标

## 四、第三阶段：NP-Guard保护策略

### 4.1 NP-Guard设计（任务 t-mnf582gc）
- 三模块架构：NP监控器→风险评估器→自适应调节器
- per-cell L_batch权重调节与scVI兼容
- 相对NP风险公式：risk(i) = max(0, 1-NP_post/NP_pre)
- 肿瘤生物学评审通过
- 生物学安全规则：doublet/伪群/应激群排除
- 输出：`np_adaptive_protection.md`

### 4.2 NP-Guard实现代码
- 两阶段训练方案（不改scVI源码）
- batch标签打乱实现
- 82%亚群改善
- 输出：`np_guard_implementation.py`（后被BD-Harmony取代）

## 五、第四阶段：关键发现

### 5.1 数据交叉分析
- 发现T cell 55%亚群被打散（ARI=0.12）
- Cluster 7假阳性（Naive/Tcm）→ Cluster 5真靶标（GITR+ TIGIT+CCR8- Treg）
- 委托数据科学查询marker基因（t-mnf4d4tn）
- 肿瘤生物学定性：皮肤癌特异的免疫抑制微环境单元

### 5.2 规模依赖性消灭
- 全量225万：pDC NP=0.231（最低）
- 8批次完好但103批次被消灭
- kNN密度偏差理论解释

## 六、第五阶段：方向迭代

### 6.1 用户批评与方向调整
- NP-Guard被否定（"选择性不整合"而非保护性整合）
- scVI+NP正则化设计（`scvi_np_regularized.py`）→也被否定
- 转向选择性整合(SI)新范式

### 6.2 选择性整合(SI)设计（任务 t-mnfkn0n5）
- 三步走：共享度评分→选择性批次校正→锚点对齐
- BD加权：z = BD*z_harmony + (1-BD)*z_pre
- 98%亚群改善，7分钟
- 肿瘤生物学审视：organ_origin协变量等4个边界情况
- 输出：`selective_integration.md`

### 6.3 CSI对比学习方案（任务 t-mnfmu6ry）
- MLP encoder + InfoNCE + 结构保持损失
- 效率~27分钟
- 训练后失败：InfoNCE的uniformity破坏连续过渡
- 输出：`native_selective_integration.md`

### 6.4 CSI失败诊断（任务 t-mnfoglx1）
- 对比学习适合分类不适合保结构
- MNN消息传递优于对比学习
- 超网络作为升级备选
- 输出：`csi_failure_analysis.md`

## 七、第六阶段：RASI端到端方案

### 7.1 RASI设计（任务 t-mnfr5osk）
- Step 0: 稀有性感知HVG（标准2000 + cluster-specific top50）
- Step 1: 200维PCA
- Step 2: BD加权选择性整合
- Step 3: NP检测
- Step 4: 低NP区域未知亚群发现
- 总时间≈45分钟（1.3倍标准流程）

### 7.2 系统性偏差调研（任务 t-mnfpqv3t）
- 8个环节偏差分析：HVG/PCA/kNN/批次建模/归一化/聚类/整合/评估
- 每点含已有文献、方案、未解决gap
- PCA对稀有亚群的偏差是完全空白领域
- 输出：`systematic_bias_survey.md`

### 7.3 原创发现（任务 t-mnfqz3eq）
- **Anchor Asymmetry Bias**：MNN锚点密度偏向大群，稀有亚群被动拖拽
- **Tokenization Collapse**：scGPT/Geneformer的rank编码压缩稀有信号

### 7.4 效率评估
- 490万细胞20分钟时间预算分配
- Harmony占50%（~8min）为瓶颈
- 总~17分钟可行

### 7.5 PCA+UMAP混合评估
- 推荐方案D（分工而非混合）：PCA用于整合，UMAP用于检测

## 八、第七阶段：论文撰写

### 8.1 论文ML部分（任务 t-mnfuo7c4）
- Methods: BD加权选择性整合完整数学描述
- Methods: NP指标定义、subcluster聚合、验证
- Supplementary: CSI/MNN-Laplacian/NP-scVI三种从零设计方案的失败诊断
- Supplementary: NP vs scIB 6个指标系统对比
- 输出：`paper_ml_sections.tex`

### 8.2 策略调整评估（任务 t-mnfe4j33）
- NP vs scIB区分：无监督是关键
- NP-Guard辩护：减量整合非不整合
- 支持加胰腺第二数据集
- Nature Methods定位

## 九、第八阶段：KAN补充方案

### 9.1 KAN设计
- KAN-Weight：[5→3→1]，输入BD/密度/纯度/MNN/特异性，输出整合权重
- KAN-NP：NP可微近似（后被数学智能体否定——NP有精确公式）
- KAN-HVG：稀有性感知HVG评分
- 三个KAN总参数<1000，训练<20分钟
- 输出：`kan_supplement.md`

### 9.2 KAN实验代码
- 完整可运行脚本，含训练、对比、样条曲线可视化、符号回归
- NaN修复：特征clip/nan处理 + Adam替代LBFGS
- 输出：`kan_weight_experiment.py`

### 9.3 团队反馈
- 数学：KAN-NP放弃，KAN-Weight [5→3→1]理论可行
- 哲学：KAN不进Main text，作Extended Data
- 肿瘤生物学：稀有基因"高偏度+高零表达+低方差"签名生物学合理
- 计算生物学：等490万结果后评估

## 十、最终交付汇总

### 分支文件（ml_workspace/）
| 文件 | 内容 |
|------|------|
| ml_algorithm_survey.md | TCI框架+29篇文献 |
| np_adaptive_protection.md | NP-Guard设计 |
| np_guard_implementation.py | NP-Guard代码（已被取代） |
| selective_integration.md | SI新范式设计 |
| native_selective_integration.md | CSI设计 |
| csi_failure_analysis.md | CSI失败诊断 |
| scvi_np_regularized.py | NP-scVI代码（已被取代） |
| systematic_bias_survey.md | 8环节偏差调研 |
| kan_supplement.md | KAN 3应用点设计 |
| kan_weight_experiment.py | KAN实验代码 |
| paper_ml_sections.tex | 论文LaTeX |
| cnem_failure_analysis_and_fixes.md | CNEM分析 |
| cnem_ml_analysis.md | CNEM ML视角 |
| methods_ml_contribution.md | Methods草稿 |
| team_task_ml_plan.md | 工作计划 |
| ml_paper_materials_index.md | 素材索引 |
| ML_AGENT_FULL_LOG.md | 本文档 |

### 核心贡献
1. **NP检测指标**（AUC=0.837）——第一个无标签的亚群消灭检测
2. **BD选择性整合**（98%亚群改善）——第一个不需知道稀有亚群即可保护的范式
3. **RASI端到端方案**——5步修复9个系统性偏差
4. **CSI失败诊断**——对比学习不适合保结构的理论分析
5. **8环节偏差调研**——系统性文献综述
6. **2个原创偏差发现**——锚点不对称+基础模型词元化坍缩
7. **KAN补充方案**——可解释的非线性权重学习
8. **论文ML部分**——完整Methods+Supplementary LaTeX

### 协作统计
- 发送消息：~50条
- 接收消息：~40条
- 竞标任务：18个
- 委托任务：12个
- 与7个队友均有实质性协作
