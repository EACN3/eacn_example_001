# Systematic Biases Against Rare Cell Subpopulations in Single-Cell RNA-seq Analysis Pipelines

## Literature Review

---

## 1. HVG Selection Bias Against Rare Cell Types (Marker Genes Missed)

### Key Papers

- **Luecken et al. (2025)** - "A systematic evaluation of highly variable gene selection methods for single-cell RNA-sequencing" (Genome Biology). Comprehensive benchmark of HVG methods revealing that computationally annotated cell type labels used as ground truth introduce "double dipping" bias, potentially favoring the HVG method used to derive annotations.
- **Townes et al. (2019)** - "Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model" (Genome Biology). Demonstrates that standard HVG selection and log-CPM normalization produce false variability; proposes deviance-based feature selection using a multinomial null model and GLM-PCA as alternatives.
- **Jiang et al. (2016)** - "GiniClust: detecting rare cell types from single-cell gene expression data with Gini index" (Genome Biology). Shows that Fano factor-based HVG selection misses rare-cell-specific genes; proposes Gini index-based feature selection that specifically captures genes expressed in small subsets of cells.
- **Imoto et al. (2025)** - "Accurate highly variable gene selection using RECODE in scRNA-seq data analysis" (bioRxiv). Proposes RECODE-denoised variance for HVG selection, outperforming conventional approaches in capturing known marker genes with robustness to normalization parameters.
- **Su et al. (2026)** - "A New Framework for Explainable Rare Cell Identification in Single-Cell Transcriptomics Data" (arXiv 2601.01358). Identifies that the standard analytical framework hinders rare cell detection by relying on PCA which transforms gene expression into uninterpretable features; proposes eliminating PCA and using anomaly detection with gene-level explanations.

### Known Biases

- Standard HVG methods (Seurat, Cell Ranger flavors) select genes based on mean-variance relationships, systematically favoring genes that distinguish large populations. Genes marking rare populations (expressed in few cells) are often filtered out.
- The mean-variance trend fitting can be inflated by strongly upregulated cell-type-specific genes at high abundances, compromising detection of relevant genes in other abundance intervals.
- Genes expressed randomly across all cell types may show high variance and be selected as false positives, while rare-cell markers with sparse but biologically meaningful expression patterns are missed.
- HVG selection is confounded by technical noise, leading to misidentification of biologically relevant genes.

### Known Solutions

- **Gini index-based selection** (GiniClust, GiniClust3): Captures genes with highly unequal expression across cells, naturally selecting rare-cell markers.
- **Deviance-based selection** (Townes et al.): Uses multinomial null model rather than normal distribution, better suited to count data.
- **RECODE-denoised variance**: Models and removes technical noise before variance estimation, improving marker gene capture.
- **Hybrid strategies**: Combining standard HVG with specialized rare-cell feature selection (e.g., GiniClust3 + HVG) for both global and rare cell type levels.

### Remaining Gaps

- No consensus on how many HVGs to select; the optimal number varies by dataset and affects rare cell detection differently than abundant cell detection.
- Benchmark studies themselves are biased by the "double dipping" problem -- using HVG-derived annotations to evaluate HVG methods.
- Most methods still assume a single set of features is optimal for all analytical tasks, when in reality rare and abundant cell type detection may require fundamentally different feature sets.
- Integration of feature selection with downstream rare cell detection methods remains ad hoc.

---

## 2. PCA/Dimensionality Reduction Bias (Rare Signals Compressed)

### Key Papers

- **DeMeo & Berger (2023)** - "SCA: recovering single-cell heterogeneity through information-based dimensionality reduction" (Genome Biology). Demonstrates that PCA favors large populations defined by many genes at the expense of smaller populations; SCA requires 38% fewer marker genes than PCA to recover a rare population and can detect populations 29% smaller.
- **Sun et al. (2021)** - "A Comparison for Dimensionality Reduction Methods of Single-Cell RNA-seq Data" (Frontiers in Genetics). Shows PCA performance degrades faster than other methods as the number of cell types increases.
- **Townes et al. (2019)** - "Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model" (Genome Biology). Proposes GLM-PCA as an alternative that respects the count nature of data, avoiding artifacts from log-normalization + PCA.
- **Su et al. (2026)** - arXiv 2601.01358. Argues PCA is a primary source of deficiency in rare cell detection, transforming meaningful gene expression into abstract uninterpretable features.
- **Torabi Goudarzi & Baran Pouyan (2025)** - "Enhanced Single-Cell RNA-seq Embedding through Gene Expression and Data-Driven Gene-Gene Interaction Integration" (arXiv 2509.02639). Graph neural network approach integrating gene-gene interactions enhances rare cell population detection over standard dimensionality reduction.

### Known Biases

- PCA maximizes global variance, so principal components preferentially capture variation between large populations. Rare subpopulations contribute minimally to total variance and their signals are compressed into trailing components that are typically discarded.
- Log-normalization before PCA introduces additional distortion: the fraction of zeros becomes the largest source of variability, biasing downstream clustering.
- As dataset complexity increases (more cell types), PCA performance degrades faster than alternatives.
- The number of PCs retained is typically chosen by elbow plots or jackstraw methods that are insensitive to rare population signals.

### Known Solutions

- **Surprisal Component Analysis (SCA)**: Information-theoretic approach that focuses on statistical significance within local neighborhoods rather than global variance. Cleanly separates gamma-delta and MAIT T-cell subpopulations that are invisible to PCA, ICA, and scVI.
- **GLM-PCA**: Directly models count data without log-transformation, avoiding normalization-induced artifacts.
- **Graph-based embeddings**: Methods that incorporate gene-gene interaction structure (e.g., Cell-Leaf Graphs + GNNs) can preserve rare population signals.
- **Correspondence Analysis**: Alternative linear method that may better preserve rare cell signals.

### Remaining Gaps

- SCA and GLM-PCA are not yet standard in major pipelines (Scanpy, Seurat) and adoption remains limited.
- No systematic study has quantified exactly how many rare cells are needed for PCA to retain their signal as a function of dataset size and number of marker genes.
- The interaction between number of PCs retained and rare cell detection has not been thoroughly characterized.
- Nonlinear methods (UMAP, t-SNE) applied after PCA inherit its biases but add their own distortions that are poorly understood for rare populations.

---

## 3. kNN Graph Construction Bias (Fixed k Disadvantages Small Populations)

### Key Papers

- **Zhu et al. (2024)** - "aKNNO: single-cell and spatial transcriptomics clustering with an optimized adaptive k-nearest neighbor graph" (Genome Biology). Demonstrates that fixed k disadvantages rare cells; proposes adaptive k assignment based on local distance distributions.
- **Baran et al. (2019)** - "MetaCell: analysis of single-cell RNA-seq data using K-nn graph partitions" (Genome Biology). Uses kNN graph partitioning to create metacells, with implications for rare cell representation.
- **Wagner et al. (2018)** - "K-nearest neighbor smoothing for high-throughput single-cell RNA-Seq data" (bioRxiv). kNN smoothing method that can over-smooth rare populations when k is too large.
- **Dann et al. (2022)** - "Differential abundance testing on single-cell data using k-nearest neighbor graphs" (Nature Biotechnology). Milo method for differential abundance testing on kNN graphs, relevant to detecting compositional changes in rare populations.

### Known Biases

- Standard kNN graph construction uses a fixed k (typically 15-30) for all cells. For rare populations with fewer than k members, neighbors are drawn from other cell types, creating spurious inter-type edges that blur cluster boundaries.
- Fixed k creates asymmetric graph density: abundant populations form tight, well-connected subgraphs while rare populations are poorly connected or absorbed into neighboring clusters.
- Shared nearest neighbor (SNN) graphs amplify this bias because rare cells share few common neighbors.
- Downstream community detection algorithms (Louvain, Leiden) inherit these graph structure biases.

### Known Solutions

- **aKNNO**: Assigns small k for rare cells (removing spurious long-range connections) and large k for abundant cells (balancing local/global variance). Benchmarked on 38 simulated and 20 real datasets, outperforming both general and specialized methods.
- **Adaptive kNN strategies**: Several methods compute k based on local distance distributions, automatically adjusting graph connectivity to cell density.
- **Milo**: Uses kNN graph neighborhoods for differential abundance testing, providing statistical framework for rare population detection.

### Remaining Gaps

- Adaptive kNN methods are not yet integrated into standard Scanpy/Seurat workflows; users must implement them separately.
- The interaction between adaptive k and batch correction (which also modifies the kNN graph) is unexplored.
- No theoretical framework exists for choosing k as a function of expected rare population size.
- How adaptive kNN interacts with different community detection algorithms (Leiden vs. Louvain) for rare cell recovery is not well characterized.

---

## 4. Batch Effect Modeling Bias (Cell-Type-Specific Batch Effects)

### Key Papers

- **Zhang et al. (2024)** - "Recovery of biological signals lost in single-cell batch integration with CellANOVA" (Nature Biotechnology). Demonstrates that current integration paradigms are "unnecessarily aggressive" and proposes ANOVA-based decomposition to explicitly estimate cell- and gene-specific batch effects.
- **Hao et al. (2024)** - "Semi-supervised integration of single-cell transcriptomics data [STACAS]" (Nature Communications). Semi-supervised approach leveraging prior cell type knowledge to preserve biological variability during integration.
- **Xiong et al. (2023)** - "Batch alignment of single-cell transcriptomics data using deep metric learning [scDML]" (Nature Communications). Deep metric learning guided by initial clusters preserves subtle cell types and enables discovery of new subtypes.
- **Haghverdi et al. (2018)** - "Batch effects in single-cell RNA sequencing data are corrected by matching mutual nearest neighbours" (Nature Biotechnology). MNN-based approach that does not require equal population compositions across batches.
- **Tran et al. (2020)** - "A benchmark of batch-effect correction methods for single-cell RNA sequencing data" (Genome Biology). Shows that cell-type-specific batch effects limit correction performance even in balanced scenarios.

### Known Biases

- Most batch correction methods assume batch effects are uniform across cell types. In reality, batch effects are often cell-type-specific -- different cell types experience different magnitudes and directions of technical variation.
- When cell type compositions differ across batches, transcriptomic differences caused by compositional differences are mistakenly treated as technical biases, leading to overcorrection.
- Rare cell types present in only one batch are particularly vulnerable: they may be "corrected away" as batch-specific artifacts.
- Overcorrection is often undetectable from UMAP visualizations -- a biologically coherent-looking plot provides no guarantee the expression matrix reflects real biology.

### Known Solutions

- **CellANOVA**: Uses "pool-of-controls" design to explicitly separate batch variation from biological variation; estimates cell- and gene-specific batch effect terms.
- **STACAS**: Semi-supervised approach that uses prior cell type labels (even imprecise ones) to guide integration and preserve biological variability.
- **scDML**: Deep metric learning preserving subtle cell types through initial cluster guidance.
- **MNN (fastMNN)**: Does not assume equal compositions; only requires a shared subset between batches.

### Remaining Gaps

- No method fully solves the fundamental identifiability problem when batch and biology are perfectly confounded.
- Cell-type-specific batch effect estimation requires knowing cell types, creating a chicken-and-egg problem.
- CellANOVA's pool-of-controls design requires appropriate experimental design not always available in retrospective atlas-building.
- Overcorrection detection metrics remain inadequate -- methods can appear to work well on standard benchmarks while silently erasing rare population signals.
- The interaction between cell-type-specific batch effects and rare cell detection has received insufficient systematic study.

---

## 5. Normalization Bias Against Rare Cells

### Key Papers

- **Bacher et al. (2017)** - "SCnorm: robust normalization of single-cell RNA-seq data" (Nature Methods). Demonstrates that bulk RNA-seq normalization assumptions fail for single-cell data; proposes group-specific normalization that reduces bias for lowly expressed genes.
- **Hafemeister & Satija (2019)** - "Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression [sctransform]" (Genome Biology). Regularized negative binomial regression removing technical characteristics while preserving biological heterogeneity.
- **Lun et al. (2016)** - "Pooling across cells to normalize single-cell RNA sequencing data with many zero counts [scran]" (Bioinformatics). Deconvolution-based normalization using cell pools, with known over-normalization issues for lowly expressed genes.
- **Chen et al. (2025)** - "Transcriptome size matters for single-cell RNA-seq normalization and bulk deconvolution" (Nature Communications). Shows that variation in transcriptome size across cell types significantly impacts normalization but is often overlooked.
- **Lytal et al. (2020)** - "Normalization Methods on Single-Cell RNA-seq Data: An Empirical Survey" (Frontiers in Genetics). Comprehensive comparison finding performance varies substantially across data characteristics.

### Known Biases

- Global scaling normalization (e.g., CPM) assumes all cells have similar total mRNA content. Rare cell types with different transcriptome sizes are systematically distorted.
- Cell-wise normalization factors (including scran) introduce bias for lowly expressed genes because these genes are affected by library size to a lesser degree and become over-normalized. SCnorm comparison showed scran produced 684 false DE genes vs. zero for SCnorm between identical cell groups.
- Log-transformation after normalization makes the zero fraction the dominant source of variability, which can create artificial clusters of low-quality cells or merge rare populations with similar sparsity patterns.
- Rare cell types with distinctive transcriptome sizes (e.g., large secretory cells) are particularly affected by size-factor-based normalization.

### Known Solutions

- **SCnorm**: Groups genes by expression level and normalizes within groups, reducing over-normalization of lowly expressed genes.
- **sctransform**: Regularized negative binomial regression that models gene-level technical characteristics; analytic Pearson residuals are well suited for identifying rare cell types.
- **scran deconvolution**: Pool-based approach that handles zero counts better than global methods, though still has issues with low-expression genes.
- **scKWARN**: Kernel-weighted-average robust normalization addressing outlier sensitivity.

### Remaining Gaps

- No normalization method explicitly accounts for rare cell types or optimizes for their preservation.
- The impact of transcriptome size variation on rare cell detection has not been systematically quantified.
- Normalization and feature selection are typically treated as independent steps, but their interaction strongly affects rare cell detection.
- Current benchmarks for normalization methods do not specifically evaluate rare cell type preservation.
- The choice between sctransform and scran has large downstream effects but no clear guidelines exist for when each is preferable for rare cell recovery.

---

## 6. Clustering Bias (Resolution Affects Rare Subpopulation Detection)

### Key Papers

- **Traag et al. (2019)** - "From Louvain to Leiden: guaranteeing well-connected communities" (Scientific Reports). Leiden algorithm improving on Louvain's limitations for community detection; resolution parameter controls cluster granularity.
- **Wan et al. (2023)** - "Self-supervised deep clustering of single-cell RNA-seq data to hierarchically detect rare cell populations" (Briefings in Bioinformatics). Hierarchical approach using self-supervised learning to detect rare populations across multiple resolutions.
- **Wegmann et al. (2019)** - "CellSIUS provides sensitive and specific detection of rare cell populations from complex single-cell RNA-seq data" (Genome Biology). Post-clustering method identifying rare subpopulations within initial coarse clusters through bimodal gene expression analysis.
- **Crecitelli et al. (2025)** - "Artificial variables help to avoid over-clustering in single-cell RNA sequencing" (American Journal of Human Genetics). Proposes adding artificial null variables to detect and prevent over-clustering.
- **Zhu et al. (2024)** - "aKNNO" (Genome Biology). Optimization of hyperparameters including resolution for simultaneous detection of abundant and rare cell types.

### Known Biases

- The resolution parameter in Leiden/Louvain clustering presents a fundamental trade-off: low resolution merges rare populations into larger clusters; high resolution fragments abundant populations into artificial subclusters.
- No single resolution captures both rare and abundant populations optimally. This is sometimes called the "resolution limit" problem.
- Graph-based clustering inherits biases from kNN graph construction (see Topic 3).
- Louvain algorithm has known limitations: imprecise community division, density-dependent subpopulation identification, and badly connected communities (addressed partially by Leiden).
- Cluster stability metrics do not distinguish between biologically meaningful rare clusters and noise-driven splits.

### Known Solutions

- **Multi-resolution approaches**: Running clustering at multiple resolutions and reconciling results (e.g., Clustree visualization).
- **CellSIUS**: Post-hoc rare cell detection within coarse clusters by identifying genes with bimodal expression and cluster-specific patterns.
- **GapClust**: Lightweight algorithm for rare cell detection that operates independently of resolution parameter.
- **Hierarchical clustering**: Self-supervised deep clustering approaches that recursively identify rare populations.
- **DBSCAN-based methods** (used in GiniClust): Density-based clustering that naturally identifies outlier/rare populations without requiring a resolution parameter.
- **Artificial variable methods**: Adding null variables to provide a statistical reference for over-clustering detection.

### Remaining Gaps

- No automated method reliably selects resolution parameters that simultaneously preserve rare and abundant populations.
- The interaction between resolution parameter, kNN graph structure, and rare cell detection is poorly understood theoretically.
- Multi-resolution reconciliation methods lack formal statistical frameworks for deciding which clusters at which resolution are "real."
- Benchmark datasets for evaluating clustering of rare populations are limited and often use simulated data with unrealistic rare cell properties.
- Sub-clustering strategies (iteratively clustering within clusters) are widely used in practice but lack formal justification and stopping criteria.

---

## 7. Integration Bias (scDML, STACAS, CellANOVA, RBET Approaches to Rare Subpopulation Protection)

### Key Papers

- **Zhang et al. (2024)** - "Recovery of biological signals lost in single-cell batch integration with CellANOVA" (Nature Biotechnology). ANOVA-based framework using pool-of-controls design to recover erased biological signals after integration; explicitly quantifies batch variation at cell and gene level.
- **Hao et al. (2024)** - "Semi-supervised integration of single-cell transcriptomics data [STACAS]" (Nature Communications). Leverages prior cell type knowledge (even imprecise labels) for integration; outperforms unsupervised methods (Harmony, scVI) and supervised methods (scANVI, scGen) in preserving biological variability.
- **Xiong et al. (2023)** - "Batch alignment of single-cell transcriptomics data using deep metric learning [scDML]" (Nature Communications). Deep metric learning guided by initial clusters and nearest neighbor information preserves subtle cell types; enables discovery of new subtypes invisible when analyzing batches individually.
- **Liu et al. (2020)** - "Alignment of single-cell RNA-seq samples without overcorrection using kernel density matching" (Genome Biology). Kernel density matching approach specifically designed to avoid overcorrection.

### Known Biases in Standard Integration

- **Harmony**: Clustering step makes it "powerless to uncover rare cell types" -- in benchmarks, Harmony mistakenly mixes rare cell types with abundant populations.
- **scVI**: Deep generative model with enough capacity to learn complex batch structures can also learn and erase complex biological structures.
- **Scanorama**: Performs well for rare cell type preservation in some benchmarks but not consistently across all scenarios.
- **General pattern**: Only scDML-like, Seurat v3, fastMNN, Scanorama, and scVI (in some settings) can simultaneously remove batch effects and discern rare cell types. DESC, Liger, and Harmony fail at rare cell type preservation.

### How These Methods Protect Rare Subpopulations

| Method | Strategy | Key Mechanism |
|--------|----------|---------------|
| **CellANOVA** | Post-integration signal recovery | Uses control samples to estimate and remove only batch-related variation; recovers signals erased by other methods |
| **STACAS** | Semi-supervised anchoring | Prior cell type labels prevent rare types from being merged with abundant types during anchor finding |
| **scDML** | Deep metric learning | Initial cluster-guided learning preserves within-batch rare substructure; discovers new subtypes |
| **MNN/fastMNN** | Local matching | Does not require equal compositions; matches only shared populations, leaving batch-unique rare types less affected |

### Remaining Gaps

- No integration method has been specifically designed with rare cell preservation as a primary objective.
- CellANOVA requires a pool-of-controls experimental design that is not always available.
- STACAS requires prior cell type labels, which may not be available for truly novel rare populations.
- scDML depends on initial clustering quality, which itself may fail for rare populations (circular dependency).
- The field lacks a systematic comparison of integration methods specifically on rare cell type preservation across varying proportions (0.1%, 0.5%, 1%, 5%).
- RBET (Robust Brain Extraction Tool) appears in the literature as a neuroimaging tool rather than a single-cell integration method; no single-cell method by this name was found in the literature search.

---

## 8. Evaluation Bias (scIB Framework Blind Spots for Rare Subpopulations)

### Key Papers

- **Luecken et al. (2022)** - "Benchmarking atlas-level data integration in single-cell genomics" (Nature Methods). The foundational scIB benchmark using 14 metrics across 16 methods; introduced isolated label scores for rare cell evaluation but with significant limitations.
- **Bao et al. (2025)** - "Reference-informed evaluation of batch correction for single-cell omics data with overcorrection awareness" (Communications Biology). Proposes overcorrection-aware evaluation metrics addressing a key blind spot in scIB.
- **Zhang et al. (2023)** - "Signal recovery in single cell batch integration" (bioRxiv/Nature Biotechnology). Demonstrates that standard metrics fail to detect signal loss from overcorrection.
- **Sikkema et al. (2025)** - "Benchmarking deep learning methods for biologically conserved single-cell integration" (PMC). Enhanced benchmarking framework (scIB-E) with refined metrics for biological conservation.
- **Huang et al. (2020)** - "Accounting for cell type hierarchy in evaluating single cell RNA-seq clustering" (Genome Biology). Shows that ARI and NMI treat cell type groups as completely exchangeable, failing to account for hierarchical relationships.

### Known Blind Spots in scIB

1. **Isolated label metrics are insufficient**: The isolated label score (ASW-based and F1-based) only evaluates labels present in few batches, not truly rare cell types (those with few cells). A cell type present in all batches but comprising 0.1% of cells is not captured by this metric.

2. **ARI/NMI treat all errors equally**: Misclassifying 100 rare cells has the same impact as misclassifying 100 abundant cells, despite the biological significance being very different. These metrics are dominated by abundant populations.

3. **Overcorrection goes undetected**: Standard batch removal metrics (kBET, graph connectivity, batch ASW) reward aggressive mixing without penalizing loss of biological signal. A method that merges all cell types across batches scores well on batch removal.

4. **Intra-cell-type variation ignored**: scIB evaluates whether cell types are preserved as discrete clusters but not whether within-cell-type biological variation (e.g., activation states, developmental trajectories) survives integration.

5. **No cell-type-weighted metrics**: All scIB metrics weight cells equally, so performance is dominated by the most abundant cell types. A method could fail completely on a rare 0.5% population without meaningful impact on aggregate scores.

6. **Ground truth dependency**: Evaluation requires pre-defined cell type labels. Novel rare populations without labels cannot be evaluated.

### Known Solutions

- **Isolated label F1 score**: Specifically evaluates clustering of labels present in few batches, but only addresses batch-specific rarity, not numerical rarity.
- **scIB-E (enhanced)**: Refined integration framework with correlation-based loss functions for better biological conservation assessment.
- **Overcorrection-aware metrics**: Reference-informed approaches that explicitly test whether integration introduces distortion.
- **Cell-type-specific evaluation**: Breaking down metrics by cell type to identify failures on specific populations (not yet standard practice).
- **CellANOVA's signal recovery assessment**: Comparing pre- and post-recovery signals provides a complementary evaluation of integration quality.

### Remaining Gaps

- No standard metric specifically evaluates rare cell type preservation as a function of cell type proportion.
- The field needs weighted metrics where rare cell types contribute proportionally more to the overall score (inverse-frequency weighting).
- There is no benchmark dataset with carefully controlled rare populations at varying frequencies (0.01% to 5%) for systematic evaluation.
- Evaluation of whether integration preserves the ability to perform differential expression analysis on rare populations is absent from all current benchmarks.
- The interaction between multiple pipeline steps (HVG selection -> PCA -> kNN -> integration -> clustering) on rare cell detection has never been evaluated end-to-end.
- scIB metrics are computed on embeddings/graphs but not on corrected expression matrices, missing distortions that affect downstream gene-level analyses.

---

## Cross-Cutting Themes

### The Compounding Effect

Each step in the standard scRNA-seq pipeline introduces its own bias against rare populations, and these biases compound:

1. **HVG selection** removes rare-cell marker genes
2. **PCA** compresses remaining rare signals into discarded components
3. **kNN graph** connects rare cells to wrong neighbors
4. **Batch correction** may erase rare populations entirely
5. **Clustering** at any single resolution misses rare groups
6. **Evaluation** metrics are blind to these failures

No study has quantified this full compounding effect end-to-end.

### Recommendations from the Literature

1. Use **hybrid feature selection** combining standard HVG with Gini-based or deviance-based methods
2. Consider **SCA or GLM-PCA** instead of standard PCA for datasets where rare populations are expected
3. Use **adaptive kNN** (aKNNO) rather than fixed k
4. Apply **semi-supervised integration** (STACAS) or **post-integration signal recovery** (CellANOVA) when possible
5. Run clustering at **multiple resolutions** and use specialized rare cell detection tools (CellSIUS, GapClust)
6. Evaluate integration using **cell-type-stratified metrics**, not just aggregate scores
7. Design experiments with **appropriate controls** to enable CellANOVA-style batch effect decomposition

### Critical Gap: No Unified Rare-Cell-Aware Pipeline

Despite extensive work on individual pipeline components, no integrated pipeline exists that is specifically designed to preserve rare cell populations across all analytical steps. The field needs:

- A rare-cell-aware pipeline that chains optimized methods across all steps
- End-to-end benchmarks with controlled rare population frequencies
- Theoretical frameworks for minimum detectable population size as a function of dataset parameters
- Standard evaluation metrics weighted by cell type rarity

---

*Literature search conducted April 2026 using arXiv, Semantic Scholar, and web search across primary literature sources including Nature Methods, Nature Biotechnology, Genome Biology, and Bioinformatics.*
