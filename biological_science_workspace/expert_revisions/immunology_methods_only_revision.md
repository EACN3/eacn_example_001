# Immunology Revision Brief for a Methods-Only Draft

Scope: remove all wet-experiment and functional-immunology claims from `biological_science_workspace/paper_draft_v2.tex`, while retaining immune examples only as public-data, in-silico case studies.

## 1) Exact v2 content to remove or replace

- `paper_draft_v2.tex:26`  
  Remove the abstract clause beginning `Wet-lab validation confirms ...` through the CD69 / IFN-$\gamma$ / killing claims. Replace with a purely computational summary sentence.

- `paper_draft_v2.tex:39`  
  Replace `validate it across multiple datasets and scales with independent wet-lab experiments` with language limited to public datasets and computational evaluation.

- `paper_draft_v2.tex:143-163`  
  Replace the entire subsection `Biological validation: stepwise elimination of an immunosuppressive Treg community`. Specific phrases that should not survive into the methods-only draft:
  - `immunosuppressive Treg community`
  - `activated Treg functional community`
  - `previously uncharacterized Treg subtype`
  - `confirmed potent immunosuppressive function`
  - `alternative immunosuppressive pathway`

- `paper_draft_v2.tex:149-155`  
  Remove all wet-assay outcome claims:
  - CD69 activation reduction
  - IFN-$\gamma$ secretion reduction
  - A549 killing reduction

- `paper_draft_v2.tex:157`  
  Remove the mechanistic pathway claim about `TIGIT--CD155` versus `CCR8` axis. Transcript patterns alone do not justify this functional interpretation here.

- `paper_draft_v2.tex:161-163` and `paper_draft_v2.tex:287-292`  
  Remove or fully replace Figure 5, its caption, and its legend. The current package is wet-lab framed and internally inconsistent (`Fig.~5d--f` in Results, but only panels `a--e` in the legend). If a Figure 5 remains, it should be a computational case-study figure only.

- `paper_draft_v2.tex:173`  
  Replace the `Clinical implications` paragraph. The current text makes unsupported therapeutic and combination-treatment claims from a transcript-defined cluster.

- `paper_draft_v2.tex:251-257`  
  Delete the entire `Wet-lab validation` methods subsection (`PBMC isolation and sorting`, `Co-culture`, wet-lab `Statistics`).

- Optional softening elsewhere:
  - `paper_draft_v2.tex:50`: keep `CCR8` as a marker example if needed, but avoid presenting a single marker as definitive proof of Treg identity or function.
  - `paper_draft_v2.tex:145-147`: retain only the computational fate-tracking logic; remove `functional community`, `metabolically active Tregs`, and similar phenotype-strengthening language unless directly supported by the public dataset.

## 2) Paste-ready LaTeX replacement paragraphs

Use these as drop-in replacements or starting text.

```tex
To demonstrate biological relevance without introducing new experiments, we use immune-cell examples from published single-cell atlases as in-silico case studies. In this framing, the role of the immune examples is to test whether integration preserves or disperses transcriptionally distinct low-frequency states in public data, not to establish cell function by orthogonal assay.
```

```tex
\subsection*{Public-data immune case study: integration disperses a rare T cell state}
As a public-data case study, we traced the fate of T cell Cluster~5 in the pan-cancer immune atlas, a small pre-integration subcluster enriched for \textit{TNFRSF18}, \textit{TIGIT}, and \textit{FOXP3} transcripts and relatively low in \textit{CCR8}. Under the standard workflow, this subcluster showed low neighborhood preservation and fragmented across multiple post-integration clusters, indicating loss of local transcriptomic structure during alignment. Under RASI, the same cells retained substantially higher neighborhood preservation and remained more topologically compact. We therefore use this cluster only as an in-silico example of how standard integration can obscure rare immune states in public atlases; no suppressive function, pathway activity, or therapeutic relevance is inferred from these transcript patterns alone.
```

```tex
\subsection*{Public-data immune case study analysis}
For the immune case study, we analyzed only publicly available single-cell transcriptomic data and quantified subcluster behavior before and after integration using neighborhood preservation, cluster dispersion, marker enrichment, and batch composition. Marker interpretation was limited to transcript abundance patterns observed in the source atlas. No donor-derived material, cell sorting, co-culture, cytokine measurement, cytotoxicity assay, or other wet-lab procedure was used in this manuscript.
```

```tex
In the methods-only version of the paper, immune examples should be interpreted as annotation-rich computational stress tests. The central claim is that integration can erase or disperse low-frequency transcriptomic states in immune datasets, thereby affecting downstream clustering and interpretation. Any statement about cellular function, suppressive capacity, signaling mechanism, or therapeutic action requires independent experimental follow-up and is outside the scope of the present manuscript.
```

## 3) Immunology constraints for the revised paper

- Use `Treg-like`, `FOXP3/TIGIT/TNFRSF18-enriched`, `CCR8-low`, or `transcript-defined subcluster`; do not use `immunosuppressive`, `activated`, `functional community`, or `previously uncharacterized subtype`.
- Do not infer function from marker expression alone. FOXP3, TIGIT, GITR/TNFRSF18, and CCR8 support annotation hypotheses, not direct proof of suppressive activity.
- Remove all claims tied to CD69, IFN-$\gamma$, killing, PBMC assays, FACS sorting, Ficoll isolation, ELISA, LDH release, A549 coculture, or donor experiments.
- Do not claim TIGIT--CD155 signaling, alternative suppressive axes, therapeutic actionability, or drug-combination rationale without independent validation.
- Keep immune evidence limited to public-data observables: transcript enrichment, neighborhood preservation, subcluster survival/dispersion, batch composition, and atlas-scale reproducibility.
- If sample-processing metadata such as `CD45+`, `CD3+`, `Ficoll`, or `no enrichment` remain in the manuscript, present them only as source-study metadata that may confound integration, not as experiments performed for this paper.

## 4) Branch / commit recommendation

- Branch: `immunology/methods-only-revision-brief`
- Commit: `Add immunology methods-only revision brief`
