# Machine-Learning Methods-Only Revision Brief

Scope: keep the paper computational and methodological. Remove all claims that require new wet-lab evidence, functional validation, or therapy-oriented biological interpretation.

## 1. Exact wet-experiment content to remove or replace from v2

| v2 location | Current content | Action for methods-only draft |
| --- | --- | --- |
| `biological_science_workspace/paper_draft_v2.tex:25-27` | Abstract sentence beginning `Wet-lab validation confirms...` with CD69, IFN-$\gamma$, and killing effect sizes. | Remove. Replace with a computational-only closing sentence about unsupervised detection, cross-dataset validation, and atlas-scale runtime. |
| `biological_science_workspace/paper_draft_v2.tex:39` | `...validate it across multiple datasets and scales with independent wet-lab experiments.` | Replace with `...validate it across multiple public datasets and scales using unsupervised computational diagnostics.` |
| `biological_science_workspace/paper_draft_v2.tex:143-164` | Entire subsection `Biological validation: stepwise elimination of an immunosuppressive Treg community`. | Remove as written. If the T cell case study is worth keeping, rewrite it as a purely computational fate-tracking example and keep only the integration-survival/dispersion observations from lines 145-147. |
| `biological_science_workspace/paper_draft_v2.tex:145-157` | Biological interpretation terms such as `functional community`, `previously uncharacterized Treg subtype`, and the TIGIT--CD155/CCR8 pathway claim. | Replace with neutral language such as `computationally defined Treg-associated subcluster` and avoid mechanistic or functional claims. |
| `biological_science_workspace/paper_draft_v2.tex:149-155` | Co-culture, CD69, IFN-$\gamma$, and A549 killing results. | Remove entirely. These are direct wet-lab claims. |
| `biological_science_workspace/paper_draft_v2.tex:159-163` and `:287-293` | Figure 5 and its legend. | Remove or replace with a computational-only figure showing pre/post-integration dispersion and RASI rescue for the retained T cell case study. No flow cytometry, ELISA, or killing panels should remain. |
| `biological_science_workspace/paper_draft_v2.tex:173` | Discussion paragraph `Clinical implications...` including anti-GITR / anti-TIGIT therapy claims. | Remove. Replace with a computational downstream-analysis implication paragraph. |
| `biological_science_workspace/paper_draft_v2.tex:251-257` | Entire `Wet-lab validation` subsection in Methods. | Remove entirely. |
| `biological_science_workspace/paper_draft_v2.tex:200` | `\bibitem{declerck2018}` | Drop if no longer cited after removing the wet-lab/mechanistic text. |

Keep `Enrichment strategy bias` (`paper_draft_v2.tex:70-72`): it discusses dataset acquisition heterogeneity as a computational confounder, not a new wet experiment performed for this paper.

## 2. Ready-to-use LaTeX paragraphs for the methods-only paper

```latex
\paragraph{RASI overview.}
RASI (Rare-Aware Single-cell Integration) is a rare-state-preserving meta-integration framework for atlas-scale scRNA-seq analysis. Its central premise is that cross-batch alignability is cell-specific rather than universal: cells embedded in broadly shared neighborhoods should be integrated strongly, whereas cells occupying sparse or batch-restricted neighborhoods should remain closer to the pre-integration manifold. This design directly targets the dominant failure mode of standard workflows, in which minority states are forced into majority neighborhoods during global correction. RASI therefore treats rare-state preservation as a first-class objective of integration rather than a downstream afterthought.

\paragraph{Algorithmic design.}
The pipeline combines four components. First, standard batch-aware highly variable genes are augmented with cluster-specific markers recovered from coarse unsupervised clustering, reducing the chance that rare-state identity genes are discarded before integration. Second, a hybrid pre-integration representation is constructed by concatenating PCA coordinates, which preserve global structure, with low-dimensional UMAP coordinates, which preserve local topology. Third, a per-cell batch-diversity score is computed from the pre-integration $k$-nearest-neighbor graph and used to interpolate between the unintegrated embedding and the base integrator output, so that strongly shared states receive more correction than weakly shared states. Fourth, Neighborhood Preservation (NP) quantifies how much of each cell's local neighborhood survives integration and serves as an unsupervised readout of subpopulation disruption.

\paragraph{Design objective.}
RASI is designed to optimize a tradeoff rather than a single metric. The goal is to improve batch removal for shared cell states while preserving local neighborhood structure for rare, continuous, or batch-restricted states that are most vulnerable to overcorrection. In practice, this means preferring solutions on the Pareto frontier of batch mixing and neighborhood stability, rather than maximizing mixing alone. Under this objective, an integration result is considered better only if it avoids artificial cross-type proximity while retaining or improving standard batch-correction quality.

\paragraph{Ablation strategy.}
Ablations should isolate each component of the framework against the same base integrator. The first comparison is standard HVG selection versus rare-aware HVG augmentation, reporting recovered marker coverage for rare states. The second is a PCA-only pre-integration space versus the hybrid PCA+UMAP representation, testing whether local topology improves NP. The third is full integration versus BD-weighted selective integration, quantifying whether selective correction improves NP without sacrificing batch mixing. The fourth is linear BD weighting versus nonlinear alternatives; the current results support retaining the linear form because it is both interpretable and empirically stronger than the tested polynomial alternative.
```

## 3. Computational-only validation and limitations

```latex
\paragraph{Computational validation.}
Validation should rely entirely on public datasets and unsupervised computational criteria. Across the immune atlas (105,682 cells, 8 batches), pancreas (16,382 cells, 9 batches), and full pan-cancer atlas (4,827,716 cells, 101 batches), RASI consistently improved Neighborhood Preservation relative to standard integration while maintaining or improving batch mixing. NP also functioned as an unsupervised failure detector, achieving ROC-AUC = 0.837 for identifying disrupted subclusters without requiring cell-type labels. At atlas scale, the method remained practical, completing the 4.83-million-cell run in 19.5 minutes with FAISS-based neighbor search and parallel NP computation. A methods-only manuscript can therefore support claims about detection, preservation, scalability, and cross-dataset robustness without invoking new biological experiments.

\paragraph{Limitations.}
The main limitation remains that RASI estimates batch diversity in a pre-integration space that is itself distorted by batch effects, so boundary cells can still receive imperfect correction. The hybrid embedding mitigates this bootstrap problem but does not eliminate it, and a minority of subpopulations can still worsen after selective integration. In addition, RASI is a meta-framework layered on top of an existing integrator rather than a fully generative model of batch-specific and batch-restricted variation. After removing wet-lab validation, the paper should avoid any claim that rescued subpopulations are functionally immunosuppressive, therapeutically actionable, or biologically confirmed; the defensible claim is narrower and computational: standard workflows can disperse rare states, and RASI reduces that failure mode.
```

Recommended computational-only replacement for the deleted biological-validation section: keep one short case-study paragraph on the T cell subcluster currently described in `paper_draft_v2.tex:145-147`, but present it only as an example of integration-induced dispersion and RASI-mediated rescue. Do not retain the healthy-donor PBMC assay, cytokine readout, tumor-killing assay, or mechanistic pathway interpretation.

## 4. Short branch / commit recommendation

- Branch: `ml/methods-only-revision-brief`
- Commit: `Add machine-learning brief for methods-only manuscript revision`

## 5. One cleanup item to resolve during redrafting

The manuscript currently uses inconsistent batch counts for the largest atlas: `101 batches` in the abstract and Table 1, but `103 batches` in the Results/Methods scale-analysis text. The methods-only draft should choose one definition and use it consistently.
