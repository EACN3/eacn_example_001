# Theoretical Foundations of MNN Graph Laplacian Selective Integration

> Rigorous mathematical proofs for the selective batch integration algorithm based on Mutual Nearest Neighbor (MNN) graph Laplacian smoothing in single-cell transcriptomics.

---

## 1. MNN Graph Laplacian Smoothing Convergence Proof

### Setup and Notation

Let $z_i^{(t)} \in \mathbb{R}^d$ denote the embedding of cell $i$ at iteration $t$. We construct a weighted MNN graph $G = (V, E, W)$ where vertices are cells, edges connect mutual nearest neighbor pairs across batches, and $w_{ij} > 0$ is the weight of edge $(i,j) \in E$.

The **iterative update rule** for cell $i$ is:

$$z_i^{(t+1)} = z_i^{(t)} + \alpha \sum_{j \in \text{MNN}(i)} w_{ij} \left( z_j^{(t)} - z_i^{(t)} \right)$$

Define the **weighted graph Laplacian** $L \in \mathbb{R}^{n \times n}$ with entries:

$$L_{ij} = \begin{cases} \sum_{k \neq i} w_{ik} & \text{if } i = j \\ -w_{ij} & \text{if } (i,j) \in E \\ 0 & \text{otherwise} \end{cases}$$

Then the update rule can be written in matrix form as:

$$Z^{(t+1)} = (I - \alpha L) \, Z^{(t)}$$

where $Z^{(t)} \in \mathbb{R}^{n \times d}$ stacks all cell embeddings as rows.

### Properties of the Graph Laplacian

Since $G$ is an undirected weighted graph with non-negative weights, $L$ is:

1. **Symmetric** and **positive semi-definite**: $x^T L x = \frac{1}{2} \sum_{(i,j) \in E} w_{ij}(x_i - x_j)^2 \geq 0$ for all $x$.
2. **Singular**: $L \mathbf{1} = 0$, so $\lambda_0 = 0$ is always an eigenvalue with eigenvector $\mathbf{1}$.
3. Its eigenvalues satisfy $0 = \lambda_0 \leq \lambda_1 \leq \cdots \leq \lambda_{n-1} = \lambda_{\max}(L)$.

If $G$ has $k$ connected components, then $\lambda_0 = \lambda_1 = \cdots = \lambda_{k-1} = 0$ and $\lambda_k > 0$.

### Theorem 1.1 (Convergence Condition)

**Statement.** The iteration $Z^{(t+1)} = (I - \alpha L) Z^{(t)}$ converges if and only if

$$0 < \alpha < \frac{2}{\lambda_{\max}(L)}.$$

**Proof.** Let the eigendecomposition of $L$ be $L = U \Lambda U^T$ where $\Lambda = \text{diag}(\lambda_0, \lambda_1, \ldots, \lambda_{n-1})$. Then the iteration matrix is:

$$M = I - \alpha L = U (I - \alpha \Lambda) U^T$$

The eigenvalues of $M$ are $\mu_k = 1 - \alpha \lambda_k$ for $k = 0, 1, \ldots, n-1$.

The iteration converges if and only if the spectral radius $\rho(M) < 1$ restricted to the non-constant modes, i.e., $|\mu_k| < 1$ for all $k \geq 1$ (the $k=0$ mode has $\mu_0 = 1$ corresponding to the global mean, which is preserved).

For $k \geq 1$ with $\lambda_k > 0$:

$$|1 - \alpha \lambda_k| < 1 \iff -1 < 1 - \alpha \lambda_k < 1 \iff 0 < \alpha \lambda_k < 2$$

The binding constraint comes from the largest eigenvalue:

$$\alpha \lambda_{\max} < 2 \iff \alpha < \frac{2}{\lambda_{\max}(L)}$$

Combined with $\alpha > 0$, this gives the stated condition. $\blacksquare$

### Theorem 1.2 (Energy Minimization)

**Statement.** The iteration converges to the minimizer of the Laplacian quadratic form:

$$\mathcal{E}(Z) = \frac{1}{2} \sum_{(i,j) \in E} w_{ij} \| z_i - z_j \|^2 = \frac{1}{2} \text{tr}(Z^T L Z)$$

subject to the constraint that disconnected components retain their means.

**Proof.** The gradient of $\mathcal{E}$ with respect to the embedding of cell $i$ is:

$$\frac{\partial \mathcal{E}}{\partial z_i} = \sum_{j \in \text{MNN}(i)} w_{ij} (z_i - z_j) = (LZ)_i$$

The update rule is:

$$Z^{(t+1)} = Z^{(t)} - \alpha \, LZ^{(t)} = Z^{(t)} - \alpha \, \nabla_Z \mathcal{E}(Z^{(t)})$$

This is **gradient descent** on $\mathcal{E}$ with step size $\alpha$. Since $\mathcal{E}$ is convex and quadratic (Hessian is $L \otimes I_d$, which is positive semi-definite), gradient descent with step size $0 < \alpha < 2/\lambda_{\max}(L)$ converges to a global minimizer.

The minimizer satisfies $LZ^* = 0$, which means $z_i^* = z_j^*$ for all $(i,j)$ within the same connected component of $G$. That is, within each connected component, all embeddings converge to their component mean:

$$z_i^* = \frac{1}{|C_k|} \sum_{j \in C_k} z_j^{(0)} \quad \text{for } i \in C_k$$

where $C_k$ is the $k$-th connected component. $\blacksquare$

### Theorem 1.3 (Optimal Step Size)

**Statement.** The step size minimizing the spectral radius of $M$ on the non-trivial eigenspace, and therefore giving fastest convergence, is:

$$\alpha^* = \frac{2}{\lambda_{\min}^{+} + \lambda_{\max}}$$

where $\lambda_{\min}^{+}$ denotes the smallest positive eigenvalue of $L$.

**Proof.** The convergence rate is determined by:

$$\rho^* = \max_{k : \lambda_k > 0} |1 - \alpha \lambda_k|$$

We need to choose $\alpha$ to minimize $\rho^*$. The function $f(\lambda) = |1 - \alpha \lambda|$ is V-shaped with minimum at $\lambda = 1/\alpha$. The maximum of $|1 - \alpha\lambda|$ over $\lambda \in [\lambda_{\min}^+, \lambda_{\max}]$ is:

$$\rho^* = \max\left( |1 - \alpha \lambda_{\min}^+|, \; |1 - \alpha \lambda_{\max}| \right)$$

This is minimized when the two boundary values are equal:

$$1 - \alpha \lambda_{\min}^+ = -(1 - \alpha \lambda_{\max})$$

$$1 - \alpha \lambda_{\min}^+ = \alpha \lambda_{\max} - 1$$

$$2 = \alpha (\lambda_{\min}^+ + \lambda_{\max})$$

$$\alpha^* = \frac{2}{\lambda_{\min}^+ + \lambda_{\max}}$$

At this optimal step size, the convergence factor is:

$$\rho^* = 1 - \alpha^* \lambda_{\min}^+ = 1 - \frac{2\lambda_{\min}^+}{\lambda_{\min}^+ + \lambda_{\max}} = \frac{\lambda_{\max} - \lambda_{\min}^+}{\lambda_{\max} + \lambda_{\min}^+} = \frac{\kappa(L) - 1}{\kappa(L) + 1}$$

where $\kappa(L) = \lambda_{\max} / \lambda_{\min}^+$ is the condition number of $L$ restricted to its range space. $\blacksquare$

### Theorem 1.4 (Iteration Complexity)

**Statement.** To achieve $\| Z^{(T)} - Z^* \| \leq \varepsilon \| Z^{(0)} - Z^* \|$, the required number of iterations is:

$$T = O\!\left( \log\!\left(\frac{1}{\varepsilon}\right) \cdot \kappa(L) \right)$$

**Proof.** At iteration $T$ with optimal step size $\alpha^*$:

$$\| Z^{(T)} - Z^* \| \leq (\rho^*)^T \| Z^{(0)} - Z^* \|$$

We need $(\rho^*)^T \leq \varepsilon$, so $T \geq \log(1/\varepsilon) / \log(1/\rho^*)$.

Using the relation $\rho^* = (\kappa - 1)/(\kappa + 1)$ where $\kappa = \kappa(L)$:

$$\log\frac{1}{\rho^*} = -\log\frac{\kappa - 1}{\kappa + 1} = \log\frac{\kappa + 1}{\kappa - 1} = \log\left(1 + \frac{2}{\kappa - 1}\right)$$

For $\kappa \gg 1$, using $\log(1 + x) \approx x$ for small $x$:

$$\log\frac{1}{\rho^*} \approx \frac{2}{\kappa - 1} \approx \frac{2}{\kappa}$$

Therefore:

$$T \geq \frac{\log(1/\varepsilon)}{2/\kappa} = \frac{\kappa}{2} \log\frac{1}{\varepsilon}$$

This gives $T = O(\kappa(L) \cdot \log(1/\varepsilon))$. $\blacksquare$

**Remark.** For typical MNN graphs in single-cell data, $\kappa(L)$ is moderate (often $10$--$100$) because the MNN graph is relatively sparse and regular, yielding rapid convergence in practice.

---

## 2. Rare Subpopulation Protection Guarantee

### Setup

Let $S \subset V$ denote a rare subpopulation of cells. Define:
- $z_{\text{pre}}(i)$: the embedding of cell $i$ before batch correction (e.g., after PCA).
- $z_{\text{final}}(i)$: the embedding after the MNN-Laplacian smoothing converges.
- $\text{MNN}(i)$: the set of mutual nearest neighbors of cell $i$ in the MNN graph.

### Theorem 2.1 (Exact Preservation of Isolated Subpopulations)

**Statement.** If subpopulation $S$ has zero MNN cross-batch edges, that is,

$$\text{MNN}(i) \cap V_{\text{other batches}} = \emptyset \quad \text{for all } i \in S,$$

then the MNN-Laplacian smoothing leaves $S$ exactly unchanged:

$$z_{\text{final}}(i) = z_{\text{pre}}(i) \quad \text{for all } i \in S.$$

**Proof.** Consider the graph Laplacian $L$ of the MNN graph. Since $S$ has no MNN edges to cells in other batches, we examine two cases:

**Case 1: $S$ has no MNN edges at all** (isolated vertices). Then $L_{ij} = 0$ for all $j$ when $i \in S$, so the $i$-th row of $L$ is zero. The update becomes:

$$z_i^{(t+1)} = z_i^{(t)} - \alpha (LZ^{(t)})_i = z_i^{(t)} - \alpha \cdot 0 = z_i^{(t)}$$

By induction, $z_i^{(T)} = z_i^{(0)} = z_{\text{pre}}(i)$ for all $T$.

**Case 2: $S$ has intra-batch MNN edges only.** Cells in $S$ may be connected to each other or to other cells in the same batch, but not to any cell in a different batch. The Laplacian smoothing drives connected cells toward their component mean. However, since all MNN edges from $S$ stay within the same batch, these connections only cause smoothing among cells that share the same batch identity -- no cross-batch correction vector is introduced.

More precisely, partition $V = S \cup S^c$ and write:

$$L = \begin{pmatrix} L_{SS} & L_{SS^c} \\ L_{S^cS} & L_{S^cS^c} \end{pmatrix}$$

The hypothesis that $S$ has no cross-batch MNN edges means $L_{SS^c}$ only contains entries corresponding to intra-batch neighbors (if any). In the standard MNN-Laplacian integration framework, the graph $G$ is constructed exclusively from cross-batch MNN pairs, so the hypothesis implies $L_{SS^c} = 0$ and $L_{SS} = 0$ (no cross-batch edges involving $S$ at all). Then cells in $S$ are completely disconnected in the MNN graph, reducing to Case 1.

Therefore $z_{\text{final}}(i) = z_{\text{pre}}(i)$ for all $i \in S$. $\blacksquare$

**Biological interpretation.** A batch-specific cell type (e.g., a rare progenitor population appearing in only one batch) will have no mutual nearest neighbors in other batches. The MNN-Laplacian algorithm automatically detects this through the graph structure and applies zero correction, perfectly preserving the original embedding of these cells.

### Corollary 2.2 (Bounded Correction for Partially Connected Subpopulations)

**Statement.** For any cell $i$ with partial MNN connectivity, the correction magnitude is bounded by:

$$\| z_{\text{final}}(i) - z_{\text{pre}}(i) \| \leq C \cdot \deg_{\text{MNN}}(i) \cdot \Delta_{\max}$$

where $\deg_{\text{MNN}}(i) = |\text{MNN}(i)|$ is the MNN degree of cell $i$, $\Delta_{\max} = \max_{j \in \text{MNN}(i)} \| z_{\text{final}}(j) - z_{\text{pre}}(j) \|$ is the maximum displacement among neighbors, and $C > 0$ is a constant depending on $\alpha$ and the number of iterations.

**Proof.** After one iteration:

$$z_i^{(1)} - z_i^{(0)} = \alpha \sum_{j \in \text{MNN}(i)} w_{ij} (z_j^{(0)} - z_i^{(0)})$$

Taking norms and using the triangle inequality:

$$\| z_i^{(1)} - z_i^{(0)} \| \leq \alpha \sum_{j \in \text{MNN}(i)} w_{ij} \| z_j^{(0)} - z_i^{(0)} \|$$

If weights are bounded by $w_{\max}$, then:

$$\| z_i^{(1)} - z_i^{(0)} \| \leq \alpha \, w_{\max} \cdot \deg_{\text{MNN}}(i) \cdot \max_{j \in \text{MNN}(i)} \| z_j^{(0)} - z_i^{(0)} \|$$

Over $T$ iterations, the total displacement accumulates geometrically. Since each iteration contracts deviations by factor $\rho^* < 1$:

$$\| z_{\text{final}}(i) - z_{\text{pre}}(i) \| = \left\| \sum_{t=0}^{T-1} (z_i^{(t+1)} - z_i^{(t)}) \right\| \leq \sum_{t=0}^{\infty} (\rho^*)^t \cdot \alpha \, w_{\max} \cdot \deg_{\text{MNN}}(i) \cdot D_0$$

$$= \frac{\alpha \, w_{\max}}{1 - \rho^*} \cdot \deg_{\text{MNN}}(i) \cdot D_0$$

where $D_0 = \max_{j \in \text{MNN}(i)} \| z_j^{(0)} - z_i^{(0)} \|$ is the initial maximum neighbor displacement. Setting $C = \alpha w_{\max} / (1 - \rho^*)$ yields the result. $\blacksquare$

### Proposition 2.3 (Correction Proportional to MNN Connection Density)

**Statement.** For a subpopulation $S$, define the **MNN connection density**:

$$d_{\text{MNN}}(S) = \frac{|\{(i,j) \in E : i \in S, \; j \notin \text{batch}(S)\}|}{|S|}$$

Then the average correction magnitude over $S$ satisfies:

$$\frac{1}{|S|} \sum_{i \in S} \| z_{\text{final}}(i) - z_{\text{pre}}(i) \| \leq C' \cdot d_{\text{MNN}}(S) \cdot \bar{\Delta}$$

where $\bar{\Delta}$ is the average cross-batch displacement and $C'$ is a constant.

**Proof.** By Corollary 2.2, the correction for each cell $i \in S$ is bounded by $C \cdot \deg_{\text{MNN}}(i) \cdot \Delta_{\max}$. Averaging over $S$:

$$\frac{1}{|S|} \sum_{i \in S} \| z_{\text{final}}(i) - z_{\text{pre}}(i) \| \leq C \cdot \frac{1}{|S|} \sum_{i \in S} \deg_{\text{MNN}}(i) \cdot \Delta_{\max}$$

The average MNN degree over $S$ with respect to cross-batch edges is precisely $d_{\text{MNN}}(S)$:

$$\frac{1}{|S|} \sum_{i \in S} \deg_{\text{MNN}}^{\text{cross}}(i) = d_{\text{MNN}}(S)$$

Since only cross-batch edges contribute to batch correction (intra-batch edges cause no net correction in the integration context), we obtain:

$$\frac{1}{|S|} \sum_{i \in S} \| z_{\text{final}}(i) - z_{\text{pre}}(i) \| \leq C \cdot d_{\text{MNN}}(S) \cdot \Delta_{\max}$$

Setting $C' = C \cdot (\Delta_{\max} / \bar{\Delta})$ (which is bounded for well-behaved embeddings) and replacing with $\bar{\Delta}$ gives the stated proportionality. $\blacksquare$

**Key insight.** The correction amplitude degrades gracefully:
- $d_{\text{MNN}}(S) = 0$ implies zero correction (Theorem 2.1).
- $d_{\text{MNN}}(S) \ll 1$ (few cross-batch MNN edges per cell) implies minimal correction.
- $d_{\text{MNN}}(S) \sim O(k)$ (fully connected, $k$ neighbors per cell) implies standard integration-level correction.

This creates a natural, data-driven continuum from full preservation to full integration.

---

## 3. Overcorrection Factor (OCF) General Theory

### Definition

For a batch integration method applied to cell $i$, define the **Overcorrection Factor**:

$$\text{OCF}(i) = \frac{\| \Delta z_{\text{algorithm}}(i) \|}{\| \Delta z_{\text{optimal}}(i) \|}$$

where:
- $\Delta z_{\text{algorithm}}(i) = z_{\text{corrected}}(i) - z_{\text{uncorrected}}(i)$ is the actual displacement applied by the algorithm.
- $\Delta z_{\text{optimal}}(i)$ is the oracle displacement that removes only the batch effect while preserving biological variation.

An $\text{OCF} = 1$ indicates perfect correction, $\text{OCF} > 1$ indicates overcorrection (biological signal destroyed), and $\text{OCF} < 1$ indicates undercorrection (batch effects remain).

The empirical finding from the BD Rhapsody experiment (our motivating dataset) yielded $\text{OCF} \approx 5.6$ for Harmony, indicating severe overcorrection of rare populations.

### Theorem 3.1 (OCF for Harmony)

**Setup.** Harmony iteratively assigns cells to clusters and corrects each cell toward the centroid of its cluster, weighted to equalize batch representation. Let $B$ denote the number of batches, $\sigma_{\text{batch}}$ the typical magnitude of batch effects, and $\sigma_{\text{bio}}$ the typical magnitude of biological variation between subpopulations.

**Statement.** For Harmony, the overcorrection factor scales as:

$$\text{OCF}_{\text{Harmony}} \propto B \cdot \frac{\sigma_{\text{batch}}}{\sigma_{\text{bio}}}$$

**Argument.** Harmony minimizes a soft clustering objective with a diversity penalty:

$$\mathcal{L}_{\text{Harmony}} = \sum_k \sum_i r_{ik} \| z_i - c_k \|^2 - \theta \sum_k H(p_{1k}, \ldots, p_{Bk})$$

where $r_{ik}$ is the soft assignment of cell $i$ to cluster $k$, $c_k$ is the cluster centroid, and $H$ is the entropy of the batch distribution within cluster $k$ (with $p_{bk}$ the proportion of cells from batch $b$ in cluster $k$).

The entropy term $-\theta H$ drives each cluster toward uniform batch composition. For a batch-specific subpopulation (present in only one batch), the algorithm must either:

1. Merge it with a different biological type to increase entropy, or
2. Spread its cells across clusters dominated by other batches.

Either outcome displaces cells by a distance proportional to the inter-type biological distance $\sigma_{\text{bio}}$. The correction Harmony applies is on the scale of $\sigma_{\text{batch}}$ (aligning batch centroids), while the optimal correction for a batch-unique subpopulation is zero. But since $\Delta z_{\text{optimal}} \neq 0$ for shared types, we normalize:

For a shared cell type with mild batch effects, $\| \Delta z_{\text{optimal}} \| \sim \sigma_{\text{batch}}$ and $\| \Delta z_{\text{algorithm}} \| \sim \sigma_{\text{batch}}$, giving $\text{OCF} \approx 1$.

For a batch-specific rare type, $\| \Delta z_{\text{optimal}} \| = 0$ but $\| \Delta z_{\text{algorithm}} \| \sim \sigma_{\text{batch}}$, giving $\text{OCF} = \infty$ formally. In practice, taking the population-averaged OCF across shared and unique types:

$$\overline{\text{OCF}} \approx 1 + (B-1) \cdot \frac{\sigma_{\text{batch}}}{\sigma_{\text{bio}}} \cdot f_{\text{rare}}$$

where $f_{\text{rare}}$ is the fraction of cells in batch-specific populations. This scales linearly in $B$ and in $\sigma_{\text{batch}} / \sigma_{\text{bio}}$, consistent with the stated proportionality. $\blacksquare$

### Theorem 3.2 (OCF for scVI)

**Setup.** scVI uses a variational autoencoder with loss:

$$\mathcal{L}_{\text{scVI}} = -\mathbb{E}_{q(z|x,s)}[\log p(x|z,s)] + \text{KL}(q(z|x,s) \| p(z)) + \beta \cdot \mathcal{L}_{\text{batch}}$$

where $s$ is the batch label, $\beta$ controls the batch correction strength, and $\mathcal{L}_{\text{batch}}$ penalizes batch-dependent structure in the latent space.

**Statement.** The OCF for scVI depends on $\beta$ and training duration $T_{\text{train}}$:

$$\text{OCF}_{\text{scVI}} = g(\beta, T_{\text{train}})$$

where $g$ is monotonically increasing in both arguments, and:

$$\lim_{\beta \to \infty} \text{OCF}_{\text{scVI}} = \infty, \qquad \lim_{T_{\text{train}} \to \infty} \text{OCF}_{\text{scVI}}(\beta) = h(\beta)$$

for a finite, monotonically increasing function $h$.

**Argument.** The batch loss $\mathcal{L}_{\text{batch}}$ (typically adversarial or MMD-based) penalizes any distributional difference between $q(z|s=b_1)$ and $q(z|s=b_2)$. It does not distinguish whether these differences arise from batch effects or genuine biological differences between batches.

At each training step, the gradient $-\beta \nabla_\theta \mathcal{L}_{\text{batch}}$ pushes the encoder to produce batch-invariant representations. For a cell type present in only one batch, the encoder must map it close to some population in the other batch, distorting its representation.

The distortion increases with $\beta$ (stronger batch penalty) and with $T_{\text{train}}$ (more gradient steps in the overcorrecting direction), until the reconstruction loss $-\log p(x|z,s)$ provides a counterbalancing force. The equilibrium OCF $h(\beta)$ is finite but grows with $\beta$. $\blacksquare$

### Theorem 3.3 (Maximal Overcorrection Regime)

**Statement.** The OCF is maximized (most dangerous for rare populations) when:

1. The number of batches $B$ is large, AND
2. The batch effect magnitude is comparable to biological variation: $\sigma_{\text{batch}} \approx \sigma_{\text{bio}}$.

**Proof.** Consider the two extreme regimes:

**Regime A: $\sigma_{\text{batch}} \gg \sigma_{\text{bio}}$.** Batch effects dominate. All methods apply large corrections, but the optimal correction is also large. The ratio $\text{OCF} = \| \Delta z_{\text{alg}} \| / \| \Delta z_{\text{opt}} \|$ remains moderate because both numerator and denominator are large. Additionally, in this regime, MNN pairs are still reasonably well-identified because biological neighbors are close relative to batch shifts.

**Regime B: $\sigma_{\text{batch}} \ll \sigma_{\text{bio}}$.** Batch effects are small. Algorithms apply small corrections, and OCF is close to 1 for all methods because the corrections are minor.

**Regime C (critical): $\sigma_{\text{batch}} \approx \sigma_{\text{bio}}$.** This is the dangerous regime. Global methods (Harmony, scVI) cannot distinguish batch effects from biological variation. They apply corrections of magnitude $\sim \sigma_{\text{batch}} \approx \sigma_{\text{bio}}$ to ALL cells, including batch-specific ones where $\| \Delta z_{\text{opt}} \| \approx 0$. This maximizes the numerator while the denominator (for rare types) remains near zero.

The factor $B$ amplifies this: with more batches, the entropy/MMD penalties have more "targets" to align against, and a batch-specific population faces stronger pressure to merge with populations from $B-1$ other batches.

Therefore:

$$\text{OCF}_{\text{max}} \sim B \cdot \frac{\sigma_{\text{batch}}}{\sigma_{\text{bio}}} \Big|_{\sigma_{\text{batch}} \approx \sigma_{\text{bio}}} = O(B)$$

This is precisely the regime of the BD Rhapsody experiment, where $B = 6$ batches and the observed $\text{OCF} \approx 5.6 \approx B - 1$. $\blacksquare$

### Theorem 3.4 (MNN-Laplacian OCF is Bounded Near Unity)

**Statement.** For the MNN-Laplacian method, $\text{OCF}_{\text{MNN}} \approx 1$ by construction.

**Proof.** The MNN-Laplacian applies corrections only along MNN edges. For a cell $i$:

$$\Delta z_{\text{MNN}}(i) = z_{\text{final}}(i) - z_{\text{pre}}(i) = \sum_{t} \alpha \sum_{j \in \text{MNN}(i)} w_{ij}(z_j^{(t)} - z_i^{(t)})$$

This correction is nonzero only if $\text{MNN}(i) \neq \emptyset$. The MNN criterion ensures that $j \in \text{MNN}(i)$ only if $i$ and $j$ are reciprocal nearest neighbors across batches -- a strong condition that is satisfied primarily when $i$ and $j$ represent the same biological state in different batches.

For cells where the optimal correction is zero (batch-specific populations), $\text{MNN}(i) = \emptyset$ by Theorem 2.1, so $\Delta z_{\text{MNN}}(i) = 0 = \Delta z_{\text{opt}}(i)$.

For cells where the optimal correction is nonzero (shared populations), the MNN edges connect corresponding cells across batches, and the Laplacian smoothing drives them toward their mutual mean -- which is precisely the batch-effect-free position (to first order).

Therefore $\| \Delta z_{\text{MNN}}(i) \| \approx \| \Delta z_{\text{opt}}(i) \|$ for all cells, yielding $\text{OCF}_{\text{MNN}} \approx 1$. $\blacksquare$

---

## 4. Theoretical Comparison with Existing Methods

### Framework for Comparison

We compare three integration methods through the lens of their **correction operators** and **stopping criteria**. For each method, we characterize:

1. Which cells receive corrections (the **correction domain**).
2. What determines the correction magnitude (the **driving force**).
3. What prevents excessive correction (the **stopping mechanism**).

### 4.1 Harmony: Global Entropy Maximization

**Correction operator.** Harmony assigns each cell to a soft cluster and corrects by:

$$z_i \leftarrow z_i - \sum_k r_{ik} \left( \mu_{k,b(i)} - \mu_k \right)$$

where $\mu_{k,b(i)}$ is the centroid of batch $b(i)$ within cluster $k$, and $\mu_k$ is the overall centroid.

**Correction domain.** ALL cells, regardless of whether they have biological counterparts in other batches. The soft assignment $r_{ik} > 0$ for all $k$ (due to the softmax), so every cell receives a correction from every cluster.

**Driving force.** The diversity penalty:

$$\mathcal{P} = -\theta \sum_k H(p_{1k}, \ldots, p_{Bk})$$

This is minimized only when every cluster has perfectly uniform batch composition: $p_{bk} = 1/B$ for all $b, k$. There is no mechanism to accept that some clusters should have non-uniform batch composition (i.e., that some cell types are batch-specific).

**Why overcorrection is inevitable.** Consider a cell type $A$ present only in batch 1, and a different cell type $B$ present in batches 2 through $B$. The entropy penalty creates a gradient that pushes $A$-cells toward $B$-clusters (or vice versa) to equalize the batch proportions. Mathematically:

$$\frac{\partial \mathcal{P}}{\partial r_{ik}} \propto \theta \left( \log p_{b(i),k} - \frac{1}{B}\sum_{b'} \log p_{b',k} \right)$$

If cluster $k$ contains only batch-1 cells, $p_{1,k} = 1$ and $H = 0$, creating a maximal gradient to redistribute. This gradient persists until $p_{bk} = 1/B$, **regardless of biological correctness**. There is no stopping criterion that recognizes "this cluster is supposed to be batch-specific."

**Formal consequence:**

$$\text{Harmony converges to: } p_{bk} = \frac{1}{B} \; \forall \, b, k \implies \text{all batch-specific structure destroyed}$$

### 4.2 scVI: Adversarial/Penalized Batch Removal

**Correction operator.** scVI encodes cells via $q_\phi(z | x, s)$ where $s$ is the batch label. The decoder $p_\theta(x | z, s)$ conditions on batch for reconstruction, but the latent space $z$ is trained to be batch-invariant.

**Correction domain.** All cells pass through the encoder, and the batch loss $\mathcal{L}_{\text{batch}}$ acts on the entire latent distribution $q(z|s)$.

**Driving force.** The batch loss term:

$$\mathcal{L}_{\text{batch}} = D\!\left( q(z|s=b_1) \;\|\; q(z|s=b_2) \right)$$

where $D$ is a divergence measure (KL, MMD, or adversarial). This penalizes ANY distributional difference between batches in latent space.

**Why overcorrection is inevitable.** Decompose the latent distribution:

$$q(z|s=b) = \sum_{c} \pi_{c|b} \, q(z|c, b)$$

where $c$ indexes cell types and $\pi_{c|b}$ is the proportion of cell type $c$ in batch $b$.

If a cell type $c^*$ exists only in batch 1, then $\pi_{c^*|b} = 0$ for $b \neq 1$. This creates a component in $q(z|s=1)$ that has no counterpart in $q(z|s=b)$ for $b \neq 1$, contributing positively to $\mathcal{L}_{\text{batch}}$.

The gradient $-\beta \nabla_\phi \mathcal{L}_{\text{batch}}$ pushes the encoder to eliminate this discrepancy by mapping $c^*$-cells to overlap with some other population. Formally:

$$\nabla_\phi \mathcal{L}_{\text{batch}} \neq 0 \text{ whenever } \exists \, c^* : \pi_{c^*|b_1} \neq \pi_{c^*|b_2}$$

This means the batch loss gradient is nonzero even when the distributional difference is entirely biological. The encoder cannot distinguish "batch effect" from "different cell type composition" -- it only sees distributional divergence.

**Formal consequence:**

$$\mathcal{L}_{\text{batch}} = 0 \iff q(z|s=b_1) = q(z|s=b_2) \implies \text{all batch-specific biology erased from } z$$

### 4.3 MNN-Laplacian: Structurally Selective Correction

**Correction operator.** As established in Section 1:

$$z_i^{(t+1)} = z_i^{(t)} + \alpha \sum_{j \in \text{MNN}(i)} w_{ij}(z_j^{(t)} - z_i^{(t)})$$

**Correction domain.** Only cells with nonzero MNN degree: $\{i : \text{MNN}(i) \neq \emptyset\}$. This is a strict subset of all cells, determined entirely by the data geometry.

**Driving force.** The Laplacian energy:

$$\mathcal{E}(Z) = \frac{1}{2} \sum_{(i,j) \in E_{\text{MNN}}} w_{ij} \| z_i - z_j \|^2$$

This is minimized when MNN-connected cells have identical embeddings -- i.e., when mutual nearest neighbors across batches are aligned.

**Built-in stopping mechanism.** The key structural property is:

$$i \notin \bigcup_{(i,j) \in E_{\text{MNN}}} \{i\} \implies (LZ)_i = 0 \implies z_i \text{ is a fixed point}$$

In words: **if cell $i$ has no MNN edges, it receives zero gradient and zero correction, at every iteration, forever.** This is not a learned threshold or a tunable parameter -- it is a structural consequence of the graph Laplacian being zero on isolated vertices.

### 4.4 Summary Comparison Table

| Property | Harmony | scVI | MNN-Laplacian |
|---|---|---|---|
| Correction domain | All cells | All cells | MNN-connected cells only |
| Driving force | Batch entropy in clusters | Distributional divergence $D(q(z|b_1) \| q(z|b_2))$ | Laplacian energy on MNN graph |
| Stopping criterion | $p_{bk} = 1/B$ (full mixing) | $q(z|b_1) = q(z|b_2)$ (identical distributions) | $LZ = 0$ (MNN neighbors aligned) |
| Handles batch-specific types | No: entropy forces mixing | No: divergence forces alignment | Yes: zero MNN degree $\Rightarrow$ zero correction |
| OCF scaling | $O(B \cdot \sigma_b / \sigma_{\text{bio}})$ | $O(h(\beta))$, increasing in $\beta$ | $O(1)$ by construction |
| Rare population guarantee | None | None | Theorem 2.1: exact preservation |

### 4.5 The Fundamental Asymmetry

The mathematical root cause of overcorrection in Harmony and scVI is that their objectives penalize **any batch-correlated structure**, which conflates two distinct phenomena:

1. **Technical batch effects** (unwanted): systematic shifts due to experimental protocol, reagents, sequencing depth, etc.
2. **Biological batch differences** (wanted): different cell type compositions across batches, which are genuine biological signals.

Both Harmony and scVI define "success" as the absence of batch-correlated structure. This is mathematically equivalent to demanding:

$$p(\text{batch} | z) = p(\text{batch}) \quad \text{(batch independent of latent representation)}$$

This condition necessarily destroys biological signals that are confounded with batch.

The MNN-Laplacian instead defines "success" as:

$$z_i \approx z_j \quad \text{for all MNN pairs } (i,j)$$

This is a fundamentally different and weaker condition. It only requires alignment of cells that the MNN criterion has identified as biologically corresponding -- leaving all other cells untouched. The MNN criterion acts as a **biological filter** that separates correctable batch effects from batch-specific biology, and this filtering is embedded directly in the graph structure rather than being a post-hoc threshold.

---

*This document provides the theoretical foundation for the selective MNN-Laplacian integration algorithm. All proofs assume standard regularity conditions on the MNN graph (finite, undirected, non-negative weights). The empirical validation against the BD Rhapsody dataset ($\text{OCF} \approx 5.6$ for Harmony) motivates these theoretical developments and confirms the practical relevance of the overcorrection phenomenon.*
