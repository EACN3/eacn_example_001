# Mathematics revision brief for the methods-only manuscript

## 1. Remove or replace from `paper_draft_v2.tex`

- Abstract, line 26: delete the sentence beginning `Wet-lab validation confirms ...`. Replace it with a computational-only closing sentence about label-free detection of rare-subpopulation collapse and improved retention under selective integration.
- Introduction, line 39: replace `and validate it across multiple datasets and scales with independent wet-lab experiments` with `and validate it computationally across multiple datasets and scales`.
- Results subsection `Biological validation: stepwise elimination of an immunosuppressive Treg community` (lines 143-163): remove in full.
  - If a case study is still needed, replace it with `Computational case study: fate of a rare T-cell candidate set during integration`.
  - Keep only transcriptomic/integration outcomes such as survival, neighborhood retention, leakage, and correction magnitude.
  - Remove the phenotype/function/pathway claims at lines 145, 149-157 (`immunosuppressive`, `healthy donor PBMCs`, `CD69`, `IFN-gamma`, `A549`, `TIGIT--CD155 axis`, `canonical CCR8 axis`).
- Figure 5 include/caption (lines 159-163) and Figure 5 legend (lines 287-293): remove or replace with a computational figure showing pre/post neighborhood preservation, leakage, and displacement for one candidate rare set.
- Discussion, `Clinical implications.` paragraph (line 173): remove in full. It depends on wet validation and therapeutic interpretation.
- Methods subsection `Wet-lab validation` (lines 251-257): remove in full.
- After those deletions, re-check unused citations and terms. `declerck2018` is likely removable if line 157 is removed.
- Keep `Enrichment strategy bias` (lines 70-72) only if framed as dataset-metadata heterogeneity. Otherwise rename it to `Sampling/composition heterogeneity` and avoid protocol-level wet-lab wording.

## 2. Ready-to-use LaTeX for a methods-only framing

Suggested preamble addition:

```latex
\usepackage{amsthm}
\newtheorem{definition}{Definition}
\newtheorem{proposition}{Proposition}
```

Suggested theory block:

```latex
Let $x_i \in \mathbb{R}^d$ denote the pre-integration embedding of cell $i$ and
let $z_i \in \mathbb{R}^d$ denote the post-integration embedding. For each
$i$, let $N_k^x(i)$ and $N_k^z(i)$ be the $k$-nearest-neighbor sets in the
pre-integration and post-integration spaces, respectively.

\begin{definition}[Rare-subpopulation retention]
For a candidate set $S \subseteq \{1,\dots,n\}$ with $|S|/n \le \rho$, define
\[
\operatorname{NP}_k(i) = \frac{|N_k^x(i) \cap N_k^z(i)|}{k},
\qquad
\operatorname{Ret}_k(S) = \frac{1}{|S|}\sum_{i \in S}\operatorname{NP}_k(i),
\]
\[
\operatorname{Leak}_k(S) =
\frac{1}{|S|k}\sum_{i \in S}|N_k^z(i)\setminus S|.
\]
We say that $S$ is retained at scale $k$ if
$\operatorname{Ret}_k(S) \ge \tau_{\mathrm{ret}}$ and
$\operatorname{Leak}_k(S) \le \tau_{\mathrm{leak}}$.
\end{definition}

\begin{definition}[Overcorrection score]
Let
\[
\Delta(i) = \|z_i - x_i\|_2,
\qquad
H_k^x(i) = -\sum_{b=1}^B p_b^x(i)\log p_b^x(i),
\qquad
H_k^z(i) = -\sum_{b=1}^B p_b^z(i)\log p_b^z(i),
\]
where $p_b^x(i)$ and $p_b^z(i)$ are the batch frequencies in $N_k^x(i)$ and
$N_k^z(i)$, respectively. Define the local mixing gain
\[
G_k(i) = H_k^z(i) - H_k^x(i)
\]
and the overcorrection score
\[
\operatorname{OC}_k(i) =
\frac{\Delta(i)}{\varepsilon + \max\{G_k(i),0\}}.
\]
Large $\operatorname{OC}_k(i)$ indicates large geometric displacement with
little local gain in cross-batch mixing.
\end{definition}

\begin{proposition}[Neighborhood stability under bounded displacement]
Let $r_\ell^x(i)$ be the distance from $x_i$ to its $\ell$-th nearest neighbor
in the pre-integration space. If $\|z_j - x_j\|_2 \le \eta$ for all cells $j$
and
\[
r_{k+1}^x(i) - r_k^x(i) > 4\eta,
\]
then $N_k^z(i) = N_k^x(i)$ and hence $\operatorname{NP}_k(i) = 1$.
\end{proposition}

\begin{proposition}[Selective integration bounds rare-cell drift]
If the final embedding is formed by
\[
z_{\lambda}(i) =
\lambda_i z_{\mathrm{int}}(i) + (1-\lambda_i)x_i,
\qquad 0 \le \lambda_i \le 1,
\]
then
\[
\|z_{\lambda}(i) - x_i\|_2 =
\lambda_i \|z_{\mathrm{int}}(i) - x_i\|_2.
\]
Therefore, for any candidate set $S$,
\[
\frac{1}{|S|}\sum_{i \in S}\|z_{\lambda}(i) - x_i\|_2
\le
\lambda_{\max}(S)
\frac{1}{|S|}\sum_{i \in S}\|z_{\mathrm{int}}(i) - x_i\|_2,
\]
where $\lambda_{\max}(S) = \max_{i \in S}\lambda_i$.
\end{proposition}

\begin{proposition}[Operational overcorrection flag]
For a candidate set $S$, if
\[
\operatorname{Ret}_k(S) < \tau_{\mathrm{ret}}
\quad\text{and}\quad
\frac{1}{|S|}\sum_{i \in S} G_k(i) \le \tau_{\mathrm{mix}},
\]
then the integration has rewritten the local geometry of $S$ without
substantial evidence of improved cross-batch alignment at scale $k$; $S$
should be flagged as overcorrected.
\end{proposition}
```

## 3. Computational-only validation implications

- The validation claim must change from `biologically functional rare population rescued` to `computationally coherent rare structure retained`.
- The minimum validation package should be:
  1. known rare positive controls from public data (for example pDC, epsilon, ionocyte);
  2. label-free candidate rare sets discovered in the pre-integration graph;
  3. scale tests across batch counts;
  4. ablations against the base integrator using `Ret`, `Leak`, `OC`, and batch-mixing metrics.
- The T-cell/Treg example can remain only as a transcriptomic case study. Do not claim function, therapeutic relevance, or a new biological subtype unless independent biology is restored elsewhere.
- If Figure 5 is replaced rather than deleted, make it a pure computation panel: candidate-set survival, NP distribution, leakage, displacement, and entropy gain before vs. after selective integration.
- The discussion should emphasize identifiability, geometric stability, and failure detection. It should not make clinical or mechanistic claims in the methods-only version.

## 4. Branch and commit recommendation

- Branch: `mathematics/methods-only-revision-brief`
- Commit: `Add mathematics revision brief for methods-only manuscript`
