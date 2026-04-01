# 选择性整合算法：Selective Integration (SI)

> 作者：机器学习智能体 (agent-mnez8qvx)
> 任务：t-mnfkn0n5
> 范式：先分类再整合，而非先整合再修复

---

## 核心思想

**打破"整合=对齐所有细胞"的默认假设。**

现有方法将所有细胞强制对齐到统一空间，稀有亚群因为没有跨批次对应物而被错误对齐到邻近大群。选择性整合反转这个逻辑：先判断每个细胞是否有跨批次对应物，只对有对应物的细胞做批次校正。

## 算法框架：三步走

### Step 1：共享度评分（Sharedness Score）

对每个细胞 i，度量它在多少个批次中有"同类"：

```
输入：多批次数据 {X_b}，合并后PCA嵌入 X_pca
参数：k=30

1. 在 X_pca 上构建 kNN 图（跨所有批次）
2. 对每个细胞 i，统计其 k 个邻居来自多少个不同批次：
   batch_diversity(i) = |{batch(j) : j ∈ N_k(i)}| / n_batches
3. 共享度：
   sharedness(i) = batch_diversity(i)

   - sharedness ≈ 1：该细胞的邻居来自几乎所有批次 → "共享"细胞
   - sharedness ≈ 1/n_batches：邻居只来自自己的批次 → "独有"细胞
```

**直觉**：如果一个T cell在所有批次都有对应的T cell，它的邻居会来自多个批次（高共享度）。如果一个pDC只在少数批次中出现，它的邻居大多来自同批次（低共享度）。

**关键区分**：低共享度有两个原因——(a) 真正的稀有独有亚群，(b) 批次效应太大导致同类细胞在PCA中不近。为了区分两者：

```
4. 对低共享度细胞，检查其邻居的细胞类型多样性：
   type_homogeneity(i) = max(celltype_freq in N_k(i))

   - 低共享度 + 高类型同质性 → 独有亚群（保护）
   - 低共享度 + 低类型同质性 → 批次效应导致（需要整合）
```

注意：这里的"细胞类型"不需要精确注释——粗聚类（如Leiden resolution=0.5）即可。

### Step 2：选择性批次校正

```
输入：sharedness(i) 和 type_homogeneity(i)
参数：
  τ_share = 0.3  # 共享度阈值
  τ_homo = 0.7   # 同质性阈值

分类规则：
  if sharedness(i) > τ_share:
      label(i) = "shared"      # 需要整合
  elif type_homogeneity(i) > τ_homo:
      label(i) = "unique"      # 独有亚群，不整合
  else:
      label(i) = "ambiguous"   # 批次效应导致，温和整合

整合操作：
  - "shared" 细胞：用标准Harmony/scVI/Scanorama全力整合
  - "unique" 细胞：保留原始PCA嵌入，不做任何批次校正
  - "ambiguous" 细胞：用降低强度的整合（如Harmony with theta=0.5）
```

### Step 3：统一嵌入空间

关键挑战：shared细胞在整合后的嵌入空间中，unique细胞在原始PCA空间中，两者不在同一个坐标系。

```
解决方案：锚点对齐

1. 对 shared 细胞执行标准整合 → 得到整合后嵌入 Z_integrated
2. 找到 shared 细胞中与 unique 细胞最近的"锚点"
3. 用这些锚点学一个线性变换 T，将原始PCA坐标映射到整合后坐标系：
   Z_unique = T(X_pca_unique)
   其中 T 通过最小化锚点的变换误差学习：
   T* = argmin_T Σ_{anchor} ||T(X_pca_anchor) - Z_integrated_anchor||²
4. 统一嵌入：
   Z_final = {Z_integrated  for shared cells,
              T(X_pca)      for unique cells,
              Z_integrated  for ambiguous cells (with reduced correction)}
```

**为什么线性变换**：整合主要是去除线性批次效应（均值偏移+缩放），线性变换足够且计算快。

## 效率分析

| 步骤 | 操作 | 时间（225万细胞） |
|------|------|-------------------|
| Step 1a | 全局kNN（FAISS GPU） | ~2分钟 |
| Step 1b | batch_diversity统计 | ~30秒（O(nk)） |
| Step 2 | 标准整合（仅shared细胞） | <30分钟（细胞数减少） |
| Step 3 | 线性变换拟合+映射 | ~10秒 |
| **总计** | | **~33分钟（vs 标准30分钟≈1.1倍）** |

关键加速：只整合shared细胞（通常80-90%），unique细胞跳过整合步骤。

## 与现有方法的区别

| 方法 | 策略 | 稀有亚群命运 |
|------|------|-------------|
| Harmony/scVI | 对齐所有细胞 | 被推入大群 |
| CellANOVA | 先整合，再后处理修复 | 先被消灭，再试图恢复（不完美） |
| STACAS | 用先验标签引导 | 需要已知标签 |
| **SI (本方法)** | **先分类，再选择性整合** | **直接保护——不对齐就不会被消灭** |

## NP 作为检测指标的角色

NP在SI框架中的角色从"正则化项"变为：
1. **验证指标**：SI整合后，unique细胞的NP应该显著高于标准整合
2. **调参指标**：用NP来选择τ_share和τ_homo的最优值
3. **检测指标**：整合后NP仍低的细胞 → 可能是ambiguous分类错误，需要人工检查

## 与标准5的衔接

选择性整合天然支持未知亚群发现：
1. 被分类为"unique"的细胞群本身就是候选的未知亚群
2. 对unique细胞做聚类 → 发现此前未注释的亚群
3. 这些亚群在标准整合中会被消灭，但SI保护了它们
4. 对这些候选亚群做DEG分析 → 生物学定性 → 湿实验验证

## 实现伪代码

```python
def selective_integration(adata, batch_key='batch', k=30,
                          tau_share=0.3, tau_homo=0.7,
                          integration_method='harmony'):
    """选择性整合主函数"""

    # Step 1: 共享度评分
    sc.pp.pca(adata, n_comps=50)

    # FAISS GPU kNN
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=k)
    nn.fit(adata.obsm['X_pca'])
    _, indices = nn.kneighbors()

    # batch diversity
    batches = adata.obs[batch_key].values
    n_batches = len(set(batches))
    sharedness = np.zeros(len(adata))
    for i in range(len(adata)):
        neighbor_batches = set(batches[indices[i]])
        sharedness[i] = len(neighbor_batches) / n_batches

    # 粗聚类用于type homogeneity
    sc.pp.neighbors(adata, use_rep='X_pca')
    sc.tl.leiden(adata, resolution=0.5, key_added='coarse_cluster')

    type_homo = np.zeros(len(adata))
    clusters = adata.obs['coarse_cluster'].values
    for i in range(len(adata)):
        neighbor_clusters = clusters[indices[i]]
        _, counts = np.unique(neighbor_clusters, return_counts=True)
        type_homo[i] = counts.max() / k

    # 分类
    cell_class = np.full(len(adata), 'shared')
    for i in range(len(adata)):
        if sharedness[i] <= tau_share and type_homo[i] > tau_homo:
            cell_class[i] = 'unique'
        elif sharedness[i] <= tau_share:
            cell_class[i] = 'ambiguous'

    adata.obs['si_class'] = cell_class
    n_unique = (cell_class == 'unique').sum()
    n_shared = (cell_class == 'shared').sum()
    n_ambig = (cell_class == 'ambiguous').sum()
    print(f"SI classification: {n_shared} shared, {n_unique} unique, {n_ambig} ambiguous")

    # Step 2: 只对shared+ambiguous细胞做整合
    mask_integrate = cell_class != 'unique'
    adata_to_integrate = adata[mask_integrate].copy()

    if integration_method == 'harmony':
        import harmonypy
        ho = harmonypy.run_harmony(
            adata_to_integrate.obsm['X_pca'],
            adata_to_integrate.obs, batch_key
        )
        Z_integrated = ho.Z_corr.T
    elif integration_method == 'scvi':
        from scvi.model import SCVI
        SCVI.setup_anndata(adata_to_integrate, batch_key=batch_key)
        model = SCVI(adata_to_integrate)
        model.train(max_epochs=100)
        Z_integrated = model.get_latent_representation()

    # Step 3: 统一嵌入
    # 找锚点：shared细胞中与unique细胞最近的
    shared_pca = adata_to_integrate.obsm['X_pca']
    unique_mask = cell_class == 'unique'
    unique_pca = adata[unique_mask].obsm['X_pca']

    # 学习线性变换：PCA → integrated space
    # 用所有integrated细胞的PCA和Z做最小二乘
    from numpy.linalg import lstsq
    T, _, _, _ = lstsq(shared_pca, Z_integrated, rcond=None)

    # 映射unique细胞
    Z_unique = unique_pca @ T

    # 组装最终嵌入
    Z_final = np.zeros((len(adata), Z_integrated.shape[1]))
    Z_final[mask_integrate] = Z_integrated
    Z_final[unique_mask] = Z_unique

    adata.obsm['X_si'] = Z_final
    adata.obs['sharedness'] = sharedness

    return adata
```
