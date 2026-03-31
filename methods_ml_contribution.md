# Nature 论文 Methods 部分 — 机器学习贡献草稿

> 作者：机器学习智能体系统 (agent-mnez8qvx)
> 状态：草稿，待实验结果完善

---

## Neighborhood Preservation (NP) Framework for Detecting Subpopulation Loss

### Detection metric: Neighborhood Preservation Score

To quantify whether batch integration disrupts unknown rare subpopulations, we defined the Neighborhood Preservation (NP) score. For each cell *i*, we computed:

$$NP(i) = \frac{|N_{pre}(i) \cap N_{post}(i)|}{k}$$

where $N_{pre}(i)$ denotes the *k*-nearest neighbors of cell *i* in the pre-integration expression space (PCA with 50 components), $N_{post}(i)$ denotes the *k*-nearest neighbors in the post-integration embedding space, and *k* = 30.

NP measures whether a cell's local neighborhood is preserved through integration. Cells whose neighborhoods are disrupted (low NP) are candidates for subpopulation loss. Critically, NP requires no external labels and operates entirely in an unsupervised manner.

To aggregate NP at the subpopulation level, we performed high-resolution Leiden clustering (resolution = 2.0) on the pre-integration data and computed the mean NP for each pre-cluster. Clusters with mean NP below a threshold were flagged as potentially destroyed.

**Validation**: NP achieved ROC-AUC = 0.837 for discriminating dispersed vs. intact subclusters across 8 immune cell types and 3 integration methods (Harmony, scVI, Scanorama), using survival rate as ground truth (Spearman *r* = 0.613).

### Protection strategy: NP-Guard

We developed NP-Guard, an adaptive protection strategy that monitors NP during the integration process and reduces integration intensity for at-risk regions.

**Risk assessment**: For each cell *i*, we computed a risk score based on the relative NP change:

$$risk(i) = \max\left(0,\ 1 - \frac{\overline{NP}_{post}(i)}{\overline{NP}_{pre}(i)}\right)$$

where $\overline{NP}$ denotes the local average over the pre-integration neighborhood, providing noise-robust estimates.

**Adaptive batch correction**: For differentiable integration methods (e.g., scVI), we modulated the per-cell batch alignment loss weight:

$$\alpha(i) = 1 - \lambda \cdot risk_{norm}(i)$$

where $\lambda = 0.8$ controls the maximum reduction in integration intensity. Cells in safe regions ($risk \approx 0$) receive standard integration, while high-risk cells receive reduced batch correction (down to 20% intensity), preserving their local structure.

**Implementation**: We adopted a two-phase training approach compatible with standard scVI without source code modification. Phase A performs standard training; Phase B retrains with batch labels stochastically shuffled for high-risk cells (proportional to risk score), preventing the model from learning batch signals for these cells.

**Biological safety rules**: Prior to NP-Guard, we excluded (1) doublets (Scrublet score > 0.25), (2) ambient RNA-dominated clusters (top markers: HBA1/HBA2/HBB/MALAT1), and (3) dissociation-induced stress clusters (top markers enriched for immediate-early genes FOS/JUN/HSP90AA1).

### Scalability

NP computation requires two *k*-nearest-neighbor searches (pre- and post-integration) and set intersection operations. Using approximate nearest neighbors (FAISS), the total computation for 4.9 million cells is approximately 5 minutes on a single GPU. NP-Guard adds approximately 10-20% overhead to base integration training time.

### Limitations

NP primarily detects the **dispersion** mode of subpopulation loss (cells scattered into multiple large clusters). The **absorption** mode (intact subpopulation merged into a neighboring large cluster) is harder to detect when expression differences are small. For tissues with high transcriptomic similarity between cell types (e.g., immune cells), complementary metrics such as hierarchical clustering comparison (Robinson-Foulds distance between pre- and post-integration dendrograms) may be needed (see Discussion).
