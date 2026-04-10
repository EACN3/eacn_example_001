# Methods: Computational Experiments

## Datasets

| Dataset | Cells | Batches | Celltypes | Source |
|---------|-------|---------|-----------|--------|
| Pan-cancer immune atlas | 2,256,276 | 103 | 8 immune | Kang et al. 2024 Nat Commun |
| Pan-cancer full atlas | 4,827,716 | 101 | 11 (immune+stromal) | Kang et al. 2024, Zenodo DOI:10.5281/zenodo.10651059 |
| Human pancreas | 16,382 | 9 | 14 | Baron/Muraro/Segerstolpe, figshare 24539828 |
| Human lung airway | 32,472 | 16 | 17 | Vieira Braga et al. 2019 Nat Med |

## Preprocessing

All datasets processed with scanpy (v1.12):
- Gene filtering: min_cells=10
- Normalization: sc.pp.normalize_total(target_sum=1e4)
- Log-transformation: sc.pp.log1p()
- HVG selection: sc.pp.highly_variable_genes(n_top_genes=2000, batch_key='Dataset')
- Scaling: sc.pp.scale(max_value=10)
- PCA: sc.pp.pca(n_comps=50, random_state=42)

## RASI Pipeline

### Step 0: Rare-HVG Selection
1. Standard HVG (2000 genes, batch-aware)
2. Coarse Leiden clustering (resolution=0.5, random_state=42)
3. Per-cluster DEG (sc.tl.rank_genes_groups, method='t-test', top 50 per cluster)
4. Merge: HVG_RASI = HVG_standard ∪ cluster-specific DEGs
5. Typical result: 2000 + 400-1000 rare-specific = 2400-3000 total genes

### Step 1: BD-Selective Integration
1. Compute kNN in PCA space (k=30, FAISS GPU)
2. BD(i) = |unique batches in N_k(i)| / min(k, B)
3. BD_smooth(i) = mean(BD(j) for j in N_k(i))
4. Run Harmony on PCA (50 dims, max_iter=10)
5. Blend: z_final(i) = BD_smooth(i) * z_harmony(i) + (1-BD_smooth(i)) * z_pca(i)

### Step 2: NP Detection
NP(i) = |N_pre(i) ∩ N_post(i)| / k
- N_pre: k-nearest neighbors in pre-integration PCA space
- N_post: k-nearest neighbors in post-integration embedding
- k=30 for all experiments

### Step 3: Low-NP Discovery
- Bottom 5% NP cells → DEG analysis vs remaining cells
- High-resolution Leiden on low-NP region for subgroup identification

## Integration Methods Compared
- Harmony (v0.2.0, PyTorch CUDA backend)
- BD-Harmony (RASI Step 1)
- Scanorama
- BBKNN

## GPU Acceleration
- kNN: FAISS GPU (IndexFlatL2 for <500k, IndexIVFFlat for >500k, nlist=1024, nprobe=32)
- NP computation: numba.njit(parallel=True)
- BD computation: numba.njit(parallel=True)
- Harmony: PyTorch CUDA
- Hardware: 8×NVIDIA A800 80GB

## Evaluation Metrics
- NP (Neighborhood Preservation): per-cell, unsupervised
- Subcluster survival rate: per-celltype Leiden pre/post, survival = max overlap fraction
- Batch mixing: silhouette score (ASW)
- scIB comparison: Graph Connectivity, ASW bio/batch

## Key Parameters
| Parameter | Value | Rationale |
|-----------|-------|-----------|
| k (kNN) | 30 | scIB benchmark standard |
| n_hvg | 2000 | Standard, supplemented by Rare-HVG |
| n_pcs | 50 | Standard for Harmony input |
| Leiden resolution | 0.5 (coarse), 1.0 (fine) | |
| BD alpha | 1.0 | Pareto-optimal (6-point sweep) |
| NP threshold | 5th percentile | Bottom 5% for discovery |

## Reproducibility
All random seeds set to 42. All experiments available at [GitHub repo].
