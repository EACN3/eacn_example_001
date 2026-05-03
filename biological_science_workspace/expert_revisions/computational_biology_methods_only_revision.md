# Computational Biology Methods-Only Revision Brief

Scope: convert `biological_science_workspace/paper_draft_v2.tex` into a methods-focused manuscript that relies only on computational analyses of public scRNA-seq datasets. Do not carry forward wet-lab validation, ex vivo functional claims, or therapy-oriented interpretation that depends on those experiments.

## 1. Wet-experiment content to remove or replace from v2

| Location in `paper_draft_v2.tex` | Action | Notes |
| --- | --- | --- |
| line 26 (Abstract) | Replace | Remove the sentence beginning `Wet-lab validation confirms...`. Replace with a computational-validation sentence centered on public datasets, NP, ASW, and scale-up experiments. |
| line 39 (Introduction) | Replace | Change `validate it across multiple datasets and scales with independent wet-lab experiments` to a public-data and computational-validation formulation. |
| lines 143-164 (`\subsection*{Biological validation: stepwise elimination of an immunosuppressive Treg community}` and associated figure block) | Replace | Do not keep this as a biological-validation section. Either remove it entirely or retitle it as a computational case study using only public-data fate tracking. |
| lines 145-157 | Remove or rewrite tightly | The Treg case-study setup in lines 145-147 can be retained only if stripped of novelty/function claims. Remove the wet-lab sentence at line 149 and the functional/mechanistic claims in lines 151-157. |
| line 173 (`Clinical implications`) | Replace | This paragraph depends on the wet-lab-backed Treg interpretation. Replace with a computational implications paragraph about benchmark design, atlas-scale risk, and reproducible rare-state evaluation. |
| lines 251-257 (`\subsection*{Wet-lab validation}`) | Remove | Delete the full wet-lab Methods subsection, including PBMC isolation, sorting, co-culture, ELISA/LDH assays, and replicate statistics tied to those experiments. |
| lines 287-292 (Figure 5 legend) | Replace | Remove the biological-validation legend. If a new Figure 5 is needed, use a computational figure such as subcluster fate tracking, ablation, or sensitivity analysis. |
| line 200 (`\bibitem{declerck2018}`) | Drop if unused | After removing the TIGIT/CCR8 mechanism paragraph, this citation will likely become unused. |

Computational-biology note: line 72 can stay conceptually, but it should be reframed as a `sample-processing metadata confounder` derived from public dataset annotations rather than as a wet-lab narrative.

## 2. Ready-to-use LaTeX replacement paragraphs

Suggested replacement for the end of the Introduction:

```tex
Here we recast rare-subpopulation preservation as a computational reproducibility problem and present RASI (Rare-Aware Single-cell Integration), a workflow designed to preserve local neighborhood structure during batch correction. Using only publicly available scRNA-seq datasets spanning atlas-scale immune data, pancreas, and lung airway epithelium, we quantify where standard preprocessing and integration pipelines preferentially erase low-frequency states and test whether selective correction reduces this loss without compromising batch mixing.\cite{kang2024,segerstolpe2016,montoro2018,luecken2022}
```

Suggested replacement for a Methods or Results transition paragraph:

```tex
All analyses were performed on publicly available scRNA-seq datasets with explicit batch annotations.\cite{kang2024,segerstolpe2016,montoro2018} Count matrices were processed in Scanpy\cite{wolf2018} using gene filtering, library-size normalization, log-transformation, and batch-aware highly variable gene selection. We compared a standard integration workflow with the RASI workflow, which augments HVGs with cluster-enriched markers, combines PCA and UMAP coordinates into a hybrid representation, and applies batch-diversity-weighted selective correction on top of a base integrator. Neighbor graphs were computed with FAISS\cite{johnson2022}, and all random seeds were fixed for reproducibility.
```

Suggested replacement for the removed biological-validation subsection:

```tex
\subsection*{Computational validation across public datasets}

We evaluated RASI using only computational criteria that are reproducible from public data. Neighborhood Preservation (NP) quantified the fraction of pre-integration neighbors retained after correction, while batch average silhouette width measured batch mixing. Validation was performed across datasets containing established rare populations, including pancreatic epsilon cells and airway ionocytes, and across scale-up experiments in which batch number was increased systematically within the pan-cancer atlas. This design tests rare-state preservation, overcorrection detection, and computational scalability without requiring new wet-lab experiments.\cite{segerstolpe2016,montoro2018,kang2024,luecken2022}
```

Suggested replacement for the wet-lab-driven Discussion paragraph:

```tex
\textbf{Computational implications.} Rare-subpopulation preservation should be treated as a benchmarkable property of integration workflows rather than as an issue visible only after downstream biological follow-up. A methods-focused evaluation can therefore center on transparent public-data benchmarks, algorithmic ablations, and scale-dependent stress tests, while reserving mechanistic and therapeutic interpretation for future studies built on independent evidence.
```

## 3. Computational validation additions that do not require wet-lab work

1. Validate on public datasets with known rare populations: pancreatic epsilon cells, airway ionocytes, and atlas-scale immune rare states such as pDC and mast-cell subclusters.
2. Add an ablation table for RASI components: rare-aware HVG only, hybrid embedding only, BD weighting only, and full pipeline.
3. Run rare-cell spike-in or downsampling experiments from public data to quantify performance as rare-state frequency decreases.
4. Report robustness across `k`, PC count, random seeds, and batch number; include bootstrap confidence intervals for NP and ASW.
5. Compare NP against standard scIB-style metrics on the same datasets to make the evaluation-blindness argument purely computational.

## 4. Branch and commit recommendation

Branch: `computational-biology/methods-only-revision-brief`

Commit: `Draft computational biology methods-only revision brief`
