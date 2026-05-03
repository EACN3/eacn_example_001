# Tumor Biology Methods-Only Revision Brief

Scope: revise `biological_science_workspace/paper_draft_v2.tex` into a methods-focused manuscript that uses tumor-related examples only as public-atlas computational case studies. Remove all new wet-experiment content and all tumor-biology claims that depend on ex vivo functional validation.

## 1) Exact v2 content to remove or replace

- `paper_draft_v2.tex:26`  
  Remove the abstract clause beginning `Wet-lab validation confirms ...` through the CD69 / IFN-$\gamma$ / killing readouts. Replace it with a public-data summary sentence.

- `paper_draft_v2.tex:39`  
  Replace `validate it across multiple datasets and scales with independent wet-lab experiments` with language limited to public datasets, computational benchmarks, and atlas-scale stress tests.

- `paper_draft_v2.tex:143-164`  
  Replace the entire subsection `Biological validation: stepwise elimination of an immunosuppressive Treg community`. In a methods-only paper this should become a neutral computational case study, not a biological-validation section.

- `paper_draft_v2.tex:145-147`  
  Do not keep the phenotype-strengthening phrases:
  - `activated Treg functional community`
  - `previously uncharacterized Treg subtype`
  - any wording that treats transcript co-expression as proof of function or novelty

- `paper_draft_v2.tex:149-157`  
  Remove all wet-lab and mechanistic claims:
  - FACS-sorted donor PBMC validation
  - CD69 reduction
  - IFN-$\gamma$ reduction
  - A549 killing reduction
  - `alternative immunosuppressive pathway via the TIGIT--CD155 axis`

- `paper_draft_v2.tex:161-163` and `paper_draft_v2.tex:287-292`  
  Remove or fully replace Figure 5, its caption, and its legend. The current figure package is built around biological validation. If Figure 5 remains, it should be computational only (for example: pre/post integration neighborhood preservation, cluster dispersion, and cancer-type distribution of the same public-atlas subcluster).

- `paper_draft_v2.tex:173`  
  Replace the `Clinical implications` paragraph. Do not keep claims about therapeutic actionability, anti-GITR / anti-TIGIT rationale, or skin-cancer combination therapy.

- `paper_draft_v2.tex:251-257`  
  Delete the entire `Wet-lab validation` Methods subsection, including PBMC isolation, FACS sorting, co-culture, ELISA, LDH assay, and replicate statistics tied to those experiments.

- `paper_draft_v2.tex:200`  
  Drop `\bibitem{declerck2018}` if no longer cited after the mechanistic paragraph is removed.

- Keep but soften:
  - `paper_draft_v2.tex:72` can stay if framed strictly as source-study processing metadata (`CD45+`, `CD3+`, `Ficoll`, `no enrichment`) that may confound integration in public atlases.
  - `paper_draft_v2.tex:145-147` may survive only as transcript-defined fate tracking in the atlas, without functional, mechanistic, or clinical language.

## 2) Paste-ready LaTeX replacement paragraphs

Use the blocks below as direct replacements or starting text.

```tex
Here we evaluate biological relevance using only published scRNA-seq atlases and treat tumor-related examples as computational case studies. In this methods-focused framing, the purpose of the case studies is to test whether integration preserves or disperses low-frequency transcript-defined states in public data, not to establish cellular function, pathway activity, or therapeutic relevance by orthogonal experiment.
```

```tex
\subsection*{Pan-cancer atlas case study: a skin-enriched FOXP3/TIGIT/TNFRSF18-high T cell state is dispersed by standard integration}
As a tumor-focused public-data case study, we traced the fate of T cell Cluster~5 in the pan-cancer immune atlas. Before integration, this low-frequency state was enriched in squamous skin cancer and skin-normal samples and showed elevated \textit{FOXP3}, \textit{TIGIT}, and \textit{TNFRSF18} transcripts with relatively low \textit{CCR8}. Under the standard workflow, the cluster showed low neighborhood preservation and dispersed across multiple post-integration clusters, consistent with loss of local transcriptomic structure during alignment. Under RASI, the same cells retained higher neighborhood preservation and remained more topologically compact. We therefore use this example only to illustrate how atlas-scale integration can obscure tumor-associated immune states in public data; no claim of suppressive function, lineage novelty, or therapeutic actionability is made.
```

```tex
\subsection*{Public-atlas positive controls for rare-state preservation}
To complement the pan-cancer case study, we evaluated RASI on public datasets containing established rare populations with prior literature support, including pancreatic epsilon cells and airway ionocytes. These positive controls provide a label-aware benchmark for rare-state preservation without requiring new experiments. Across these datasets, standard integration reduced neighborhood preservation and marker continuity for the rare populations, whereas RASI preserved a larger fraction of local neighbors while maintaining acceptable batch mixing. Together with the pan-cancer atlas analyses, these results support the conclusion that rare-state loss is detectable in public data across both tumor-associated and non-tumor settings.
```

```tex
All tumor-biology interpretation in this manuscript is restricted to transcript-defined states observed in published atlases. Claims are therefore limited to sample distribution, marker enrichment, neighborhood preservation, cluster dispersion, and batch composition before and after integration. Any statement about suppressive activity, immune-evasion mechanism, ligand-receptor signaling, or treatment relevance requires independent validation and is outside the scope of the present methods paper.
```

## 3) Tumor-biology constraints on clinical or functional claims

- Use `Treg-like`, `transcript-defined state`, `FOXP3/TIGIT/TNFRSF18-high`, `CCR8-low`, or `skin-enriched subcluster`; do not use `immunosuppressive`, `functional community`, `previously uncharacterized subtype`, or `therapeutically actionable`.
- Do not infer suppressive function from marker co-expression alone. `FOXP3`, `TIGIT`, `TNFRSF18`, and `CCR8` support annotation hypotheses, not direct proof of regulatory potency or tumor-promoting activity.
- Do not claim TIGIT--CD155 signaling, an alternative CCR8-independent suppressive axis, immune-evasion mechanism, or tumor-killing inhibition without orthogonal evidence.
- Remove all claims based on donor PBMCs, FACS sorting, co-culture, cytokine assays, LDH assays, A549 targets, or any other ex vivo experiment.
- If the manuscript keeps the SSCC / skin-normal enrichment observation, present it only as source-atlas sample distribution, not as proof of skin-cancer specificity or clinical enrichment.
- Do not use transcript co-expression to motivate anti-GITR, anti-TIGIT, or combination-therapy hypotheses in the revised paper.
- Keep tumor-biology interpretation anchored to what the public atlases can support directly: cancer-type distribution, tumor-versus-normal composition, rare-state dispersion after integration, and reproducibility across datasets.

## 4) Branch / commit recommendation

- Branch: `tumor-biology/methods-only-revision-brief`
- Commit: `Add tumor biology methods-only revision brief`
