# NP (Neighborhood Preservation) 方法完整描述

## 算法定义

对每个细胞 $i$，NP 定义为整合前后 k 近邻的交集比例：

$$NP(i) = \frac{|N_{pre}(i) \cap N_{post}(i)|}{k}$$

其中：
- $N_{pre}(i)$：整合前 PCA 空间中细胞 $i$ 的 $k$ 近邻集合
- $N_{post}(i)$：整合后 embedding 空间中细胞 $i$ 的 $k$ 近邻集合
- $k$：近邻数（默认 30）

## 算法伪代码

```
Input: X_pre (n×d pre-integration PCA), X_post (n×d' post-integration embedding), k
Output: NP (n×1 per-cell scores)

1. Build kNN index on X_pre using FAISS GPU
   N_pre = FAISS_kNN(X_pre, k)  # n×k index matrix

2. Build kNN index on X_post using FAISS GPU
   N_post = FAISS_kNN(X_post, k)  # n×k index matrix

3. For each cell i = 1..n:
   NP[i] = |set(N_pre[i]) ∩ set(N_post[i])| / k

Return NP
```

## 亚群级聚合

对整合前的 Leiden 聚类簇 $S$：
$$\overline{NP}(S) = \frac{1}{|S|} \sum_{i \in S} NP(i)$$

$\overline{NP}(S)$ 低 → 簇 $S$ 的成员在整合后失去了原有邻居 → 簇可能被消灭

## 参数选择

| 参数 | 默认值 | 选择理由 |
|------|--------|---------|
| k | 30 | scIB benchmark 标准值，平衡分辨率与稳定性 |
| PCA dims | 50 | 标准 scRNA-seq 降维 |
| 距离度量 | L2 (Euclidean) | FAISS GPU 默认，与 scanpy 一致 |

## 计算复杂度

- kNN 构建：$O(n \cdot d \cdot \log n)$ with FAISS GPU
- NP 计算：$O(n \cdot k)$（集合交集）
- 总体：$O(n \cdot d \cdot \log n)$，225万细胞 < 1小时

## 关键性质

1. **完全无监督**：不需要任何 celltype 标签
2. **隐式分层**：kNN 天然在同类型细胞内部，不需要显式分层
3. **跨方法一致**：Harmony AUC=0.837, Scanorama AUC=0.728
4. **检测模式**：对"吸收"（absorbed into large group）和"弥散"（scattered）两种消灭模式均敏感

## 实验结果摘要

### Phase 1 子集（105k cells, 8 datasets）

| 指标 | Harmony | Scanorama |
|------|---------|-----------|
| AUC (subcluster survival prediction) | 0.837 | 0.728 |
| Spearman r (NP vs survival) | 0.612 | 0.443 |

### 亚群消灭检测

| Celltype | 亚群数 | 被消灭(survival<50%) | 比例 | ARI |
|----------|--------|---------------------|------|-----|
| T cell | 24 | 10 | 42% | 0.167 |
| Mast | 13 | 5 | 38% | 0.359 |
| Macrophage | 32 | 9 | 28% | 0.345 |

最严重消灭案例：
- Macrophage_16: 835 cells, survival=20.0%, NP=0.314
- T cell_12: 1034 cells, survival=25.4%, NP=0.252

### 失败指标对比

| 指标 | AUC | 失败原因 |
|------|-----|---------|
| CNEM (表达空间) | 0.37-0.40 | 免疫细胞间余弦距离集中 |
| R(S) (嵌入弥散度) | 反向(-0.952) | 检测校正强度而非消灭 |
| **NP (邻域保留率)** | **0.837** | ✓ 有效 |
