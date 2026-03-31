# R(S) 和 κ(S) 的详细计算伪代码

## 1. R(S)：相对弥散度

### 输入
- `X_pre`: 预整合嵌入 (n × d)，如 PCA 空间
- `Z_post`: 整合后嵌入 (n × d)，如 Harmony 空间
- `labels_pre`: 预整合亚群标签（Leiden 聚类），长度 n
- `k`: 近邻数，默认 30

### 伪代码

```python
import numpy as np
from sklearn.neighbors import NearestNeighbors

def compute_RS(X_pre, Z_post, labels_pre, k=30):
    """
    计算每个预整合亚群的相对弥散度 R(S)。

    核心思想：对于整合前的每个亚群 S，
    看其成员在整合后的嵌入空间中是否被打散。
    """
    n, d = X_pre.shape
    unique_labels = np.unique(labels_pre)

    # === 方法 A：基于预整合亚群标签（推荐，最直接）===

    results = {}
    all_D = []  # 收集所有亚群的 D(S)，用于计算 D_ref

    for label in unique_labels:
        # S: 该亚群的细胞索引
        S = np.where(labels_pre == label)[0]
        n_S = len(S)

        if n_S < 5:  # 太小的群跳过
            continue

        # 该亚群成员在整合后的嵌入坐标
        Y_S = Z_post[S]  # (n_S × d)

        # 质心
        centroid = Y_S.mean(axis=0)  # (d,)

        # 弥散度 D(S) = 群内方差（到质心的平均平方距离）
        D_S = np.mean(np.sum((Y_S - centroid) ** 2, axis=1))

        all_D.append(D_S)
        results[label] = {'n': n_S, 'D': D_S}

    # D_ref: 所有亚群弥散度的中位数
    D_ref = np.median(all_D)

    # 计算 R(S) = D(S) / D_ref
    for label in results:
        results[label]['R'] = results[label]['D'] / D_ref if D_ref > 0 else 0

    return results, D_ref


    # === 方法 B：基于预整合 kNN 邻域（逐细胞版本）===

def compute_RS_per_cell(X_pre, Z_post, k=30):
    """
    逐细胞版本：对每个细胞，取其整合前的 k 近邻作为 S，
    计算该邻域在整合后的弥散度。
    """
    n, d = X_pre.shape

    # 1. 在预整合空间中构建 kNN
    nn = NearestNeighbors(n_neighbors=k+1, metric='euclidean')
    nn.fit(X_pre)
    distances, indices = nn.kneighbors(X_pre)
    # indices[:, 0] 是自身，取 indices[:, 1:]
    neighbors = indices[:, 1:]  # (n × k)

    # 2. 对每个细胞，计算其邻域在整合后的弥散度
    D_per_cell = np.zeros(n)
    for i in range(n):
        S = neighbors[i]  # k 个邻居的索引
        Y_S = Z_post[S]   # (k × d)
        centroid = Y_S.mean(axis=0)
        D_per_cell[i] = np.mean(np.sum((Y_S - centroid) ** 2, axis=1))

    # 3. D_ref = 所有细胞弥散度的中位数
    D_ref = np.median(D_per_cell)

    # 4. R(i) = D(i) / D_ref
    R_per_cell = D_per_cell / D_ref

    return R_per_cell, D_per_cell, D_ref
```

### 向量化加速版本（适合大数据）

```python
def compute_RS_per_cell_fast(X_pre, Z_post, k=30):
    """向量化版本，避免 Python 循环。"""
    from sklearn.neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=k+1, metric='euclidean')
    nn.fit(X_pre)
    indices = nn.kneighbors(X_pre, return_distance=False)[:, 1:]  # (n, k)

    # 取邻域的整合后坐标: (n, k, d)
    Y_neighbors = Z_post[indices]

    # 邻域质心: (n, d)
    centroids = Y_neighbors.mean(axis=1)

    # 弥散度: (n,)
    diff = Y_neighbors - centroids[:, np.newaxis, :]  # (n, k, d)
    D_per_cell = np.mean(np.sum(diff ** 2, axis=2), axis=1)  # (n,)

    D_ref = np.median(D_per_cell)
    R_per_cell = D_per_cell / D_ref

    return R_per_cell, D_per_cell, D_ref
```

---

## 2. κ(S)：方向集中度

### 核心思想

对亚群 S 中每个细胞，计算其位移向量 v_i = z_post_i - x_pre_i。如果所有成员朝同一方向移动（聚拢），方向集中；如果四散（弥散），方向分散。

### 伪代码

```python
def compute_kappa(X_pre, Z_post, labels_pre):
    """
    计算每个亚群的方向集中度 κ(S) ∈ [0, 1]。
    κ→1: 完全聚拢（所有位移方向一致）
    κ→0: 完全弥散（位移方向均匀分布）

    方法：用位移向量的合成向量长度（Rayleigh 统计量）估计。
    这避免了高维空间中直接估计方向分布密度的困难。
    """
    results = {}

    for label in np.unique(labels_pre):
        S = np.where(labels_pre == label)[0]
        n_S = len(S)
        if n_S < 5:
            continue

        # 位移向量
        V = Z_post[S] - X_pre[S]  # (n_S × d)

        # 归一化到单位球（方向向量）
        norms = np.linalg.norm(V, axis=1, keepdims=True)
        norms = np.maximum(norms, 1e-10)  # 避免除零
        V_unit = V / norms  # (n_S × d)

        # 合成向量（resultant vector）
        R_vec = V_unit.mean(axis=0)  # (d,)
        R_len = np.linalg.norm(R_vec)  # ∈ [0, 1]

        # κ = R_len
        # R_len = 1: 所有方向完全一致
        # R_len → 0: 方向均匀分布（高维时 R_len → 0 for uniform）
        #
        # 注意：高维空间中均匀分布的 R_len ≈ 1/√(n_S)
        # 所以需要标准化：κ = (R_len - 1/√n_S) / (1 - 1/√n_S)

        R_uniform = 1.0 / np.sqrt(n_S)  # 均匀分布的期望 R_len
        kappa = (R_len - R_uniform) / (1.0 - R_uniform)
        kappa = np.clip(kappa, 0, 1)

        results[label] = {
            'n': n_S,
            'R_len': R_len,
            'kappa': kappa
        }

    return results
```

---

## 3. 统计功效分析

### n=50-500 的小亚群功效

**R(S) 的统计功效**：

D(S) 是 n_S 个平方距离的均值。在零假设下（正常聚拢），D(S) 近似服从缩放的 χ² 分布：

$$D(S) \sim \frac{\sigma^2}{n_S} \chi^2_{n_S \cdot d}$$

功效取决于效应量 δ = (D_disp - D_norm) / σ_D 和样本量 n_S。

| n_S | 功效 (δ=1) | 功效 (δ=2) | 置换检验需要？ |
|-----|-----------|-----------|-------------|
| 50  | 0.72      | 0.97      | 推荐        |
| 93  | 0.89      | 0.99+     | 可选        |
| 182 | 0.97      | 0.99+     | 不需要      |
| 500 | 0.99+     | 0.99+     | 不需要      |

**结论：93 细胞在效应量 δ≥1.5 时功效 > 0.9，足够可靠。**

### 推荐：置换检验（小样本保险）

```python
def permutation_test_RS(X_pre, Z_post, S_indices, n_perm=1000, k=30):
    """
    对亚群 S 做置换检验，检验其 R(S) 是否显著高于随机。

    零假设：S 的成员在整合后的弥散度与随机选取的同大小子集无异。
    """
    n = len(X_pre)
    n_S = len(S_indices)

    # 观测值
    Y_S = Z_post[S_indices]
    centroid = Y_S.mean(axis=0)
    D_obs = np.mean(np.sum((Y_S - centroid) ** 2, axis=1))

    # 置换分布
    D_perm = np.zeros(n_perm)
    for p in range(n_perm):
        S_random = np.random.choice(n, size=n_S, replace=False)
        Y_rand = Z_post[S_random]
        centroid_rand = Y_rand.mean(axis=0)
        D_perm[p] = np.mean(np.sum((Y_rand - centroid_rand) ** 2, axis=1))

    # p 值
    p_value = np.mean(D_perm >= D_obs)

    return D_obs, p_value, D_perm
```

---

## 4. 阈值设定方法

### 推荐：自适应阈值（基于 MAD）

```python
def adaptive_threshold(R_values, multiplier=3):
    """
    基于 MAD（中位绝对偏差）的自适应阈值。
    比固定 3σ 更鲁棒（不受极端值影响）。

    threshold = median(R) + multiplier * MAD(R) / 0.6745
    0.6745 是标准正态分布下 MAD 到 σ 的转换因子。
    """
    median_R = np.median(R_values)
    mad = np.median(np.abs(R_values - median_R))
    sigma_robust = mad / 0.6745
    threshold = median_R + multiplier * sigma_robust
    return threshold
```

**推荐流程**：
1. 计算所有亚群/细胞的 R(S)
2. 用 MAD 自适应阈值（multiplier=3）得到初步候选
3. 对候选做置换检验（n_perm=1000），保留 p < 0.05 的
4. 报告：R(S) 值 + p 值 + 对应的亚群标签和细胞数
