# Data Science Methods-Only Revision Brief

Goal: convert `paper_draft_v2.tex` into a methods-focused manuscript supported only by public datasets and simulations. Keep computational case studies; remove any claim that depends on new wet-lab evidence.

## 1. Exact wet-experiment content to remove or replace

| Location in `paper_draft_v2.tex` | Action | Data-science replacement |
| --- | --- | --- |
| `L25-L26` (Abstract) | Remove the clause beginning `Wet-lab validation confirms...` | Replace with one sentence on public-data and simulation-based validation. |
| `L39` (Introduction) | Replace `validate it across multiple datasets and scales with independent wet-lab experiments` | Use `validate it across multiple public datasets, scale sweeps, and simulation-based stress tests`. |
| `L143-L164` (Results subsection `Biological validation: stepwise elimination of an immunosuppressive Treg community`) | Replace the entire subsection. Remove the wet-lab framing, the functional claims, and the assay outcomes. | If this case study is kept, rewrite it as a public-data computational case study only. |
| `L145-L157` (within that subsection) | Remove or soften these exact claims: `functional community`, `previously uncharacterized Treg subtype`, `confirmed potent immunosuppressive function`, the three assay bullets, and `alternative immunosuppressive pathway...` | Keep only label-free fate-tracking claims that can be supported directly from the public atlas. |
| `L159-L163` (Figure 5 block) | Remove `Fig5_biological_validation` | Replace with a computational ablation/stress-test figure. |
| `L173-L174` (Discussion, `Clinical implications`) | Remove the therapy/combination-treatment paragraph tied to the wet-lab result | Replace with a paragraph about computational implications for atlas-scale benchmarking and downstream analyses. |
| `L251-L257` (Methods, `Wet-lab validation`) | Delete the entire subsection | No wet-lab methods section should remain in the methods-only draft. |
| `L287-L293` (Figure 5 legend) | Delete the full legend | Replace with a legend for a computational validation figure. |
| `L200` (`\bibitem{declerck2018}`) | Drop if unused after the cuts above | Keep only if cited elsewhere after rewrite. |

Optional phrasing cleanup: `L70-L72` can stay, but rename the subsection from `Enrichment strategy bias` to `Sample-processing heterogeneity as a metadata-defined confounder` so it reads as a public-metadata analysis rather than new bench work.

## 2. Ready-to-use LaTeX text

### One-line replacements

```latex
% Abstract replacement sentence
Across public atlases and simulation-based stress tests, RASI consistently improves rare-subpopulation preservation while maintaining or improving batch mixing, indicating that label-free protection of rare states is achievable without auxiliary wet-lab validation.

% Introduction replacement clause
...and validate it across multiple public datasets, scale sweeps, and simulation-based stress tests.
```

### Computational benchmarks, metrics, and ablations

```latex
\subsection*{Computational benchmarking across public datasets}
We benchmarked RASI exclusively on published datasets spanning distinct rare-cell scenarios: a pan-cancer immune atlas for atlas-scale integration, a pancreatic islet dataset containing epsilon cells, and a lung airway dataset containing ionocytes. For each dataset, we compared the standard preprocessing and integration workflow against RASI under matched preprocessing, identical random seeds, and the same batch definitions. This design isolates whether the rare-aware components improve preservation of low-frequency structure without relying on newly generated experimental measurements.

\paragraph{Evaluation metrics.}
We evaluated each method using complementary metrics that separate rare-state preservation from batch removal. Neighborhood Preservation (NP) quantified local structural stability before and after integration, batch average silhouette width (ASW) quantified mixing, graph connectivity provided a standard benchmark reference, rare-marker retention measured whether informative genes survived feature selection, and subcluster survival measured whether pre-integration rare subclusters remained recoverable after integration. For analyses involving known rare cell types, labels were withheld during method execution and restored only for post hoc auditing.

\paragraph{Ablation analysis.}
To quantify the contribution of each component, we performed incremental ablations of rare-aware HVG augmentation, hybrid PCA-UMAP embedding, and batch-diversity-weighted selective integration. Each ablation was evaluated with the same NP, ASW, marker-retention, and subcluster-survival metrics used in the main benchmark. We additionally tested prevalence sweeps and batch-count sweeps to determine whether each component remained effective when rare populations became smaller or more batch restricted.

\paragraph{Figure and table replacements.}
The wet-lab validation figure should be replaced by a computational stress-test summary showing: (i) rare-marker retention before and after HVG augmentation, (ii) NP-versus-ASW trade-offs across ablations and BD weights, (iii) recovery of known rare populations in pancreas, lung, and immune-atlas benchmarks, and (iv) runtime scaling from mid-sized cohorts to atlas scale. Table~1 should be expanded to include rare-marker retention, subcluster-survival rate, and runtime alongside NP and ASW so that preservation, integration quality, and scalability are reported in one place.
```

### Replacement captions for Figure 5 and Table 1

```latex
\caption{\textbf{Computational stress tests and ablations for rare-subpopulation preservation.} \textbf{a}, Rare-marker retention before and after rare-aware HVG augmentation across public datasets. \textbf{b}, NP versus batch ASW across ablations and BD weights. \textbf{c}, Recovery of known rare populations after integration in pancreas, lung airway, and immune-atlas benchmarks. \textbf{d}, Runtime scaling from matched subsamples to atlas-scale data.}
```

```latex
\caption{Benchmark performance across public datasets and simulation-based stress tests. Metrics should include cell count, batch count, known rare population used for audit, standard versus RASI NP, standard versus RASI batch ASW, rare-marker retention, subcluster-survival rate, and runtime.}
```

## 3. Validation suggestions using only published/public data and simulations

1. Use published rare-cell benchmarks with known audit targets: pancreas (`epsilon`), lung airway (`ionocyte`), and the pan-cancer immune atlas (`pDC`, `mast`, rare T-cell programs).
2. Hide rare-type labels during benchmarking, compute label-free NP-based detection, then restore labels only for post hoc ROC/AUC and survival-rate auditing.
3. Run prevalence sweeps by downsampling known rare populations to `0.1%`, `0.25%`, `0.5%`, `1%`, and `2%` of total cells.
4. Run batch-restriction simulations by confining a known rare population to a subset of batches, then measuring whether integration incorrectly absorbs it into abundant neighbors.
5. Run scale sweeps on the same public atlas (`5`, `10`, `20`, `50`, `103` batches) to preserve the current phase-transition argument without any new experiments.
6. Benchmark robustness across multiple base integrators (`Harmony`, `Scanorama`, `scVI`, `BBKNN`) so the methods-only paper does not read as Harmony-specific.
7. Add a pseudo-batch stress test by splitting public datasets by donor or technology and re-running the pipeline with identical labels to quantify false rare-state destruction induced purely by integration.

## 4. Branch and commit recommendation

- Branch: `data-science/methods-only-brief`
- Commit: `Add data science methods-only revision brief`
