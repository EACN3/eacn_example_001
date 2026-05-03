# Philosophy Revision Brief for a Methods-Only Manuscript

## 1. Exact v2 content to remove or replace

- `paper_draft_v2.tex:26`  
  Remove the abstract sentence beginning `Wet-lab validation confirms...`.  
  Reason: this is the strongest experimental-validation claim in the paper and is incompatible with a methods-only scope.

- `paper_draft_v2.tex:39`  
  Replace `...validate it across multiple datasets and scales with independent wet-lab experiments.`  
  With a computational-only claim such as: `...evaluate it across multiple datasets and atlas scales using label-aware and label-free computational benchmarks.`

- `paper_draft_v2.tex:143-163`  
  Remove the entire Results subsection `Biological validation: stepwise elimination of an immunosuppressive Treg community`, including the figure environment for `Fig5_biological_validation`.  
  If a case study is still wanted, replace the whole subsection with a short computational case study of a `candidate` or `Treg-like` cluster, with no functional or mechanistic claims.

- `paper_draft_v2.tex:145-157`  
  Do not reuse the following claims elsewhere without new experimental support: `previously uncharacterized Treg subtype`, `confirmed potent immunosuppressive function`, `alternative immunosuppressive pathway`, and `renders this mechanism invisible`.

- `paper_draft_v2.tex:173`  
  Replace the `Clinical implications.` paragraph.  
  It currently turns a computational pattern into therapeutic and disease-mechanism claims. A methods-only paper should discuss downstream analytical risk, not clinical actionability.

- `paper_draft_v2.tex:251-257`  
  Remove the entire Methods subsection `Wet-lab validation` (`PBMC isolation and sorting`, `Co-culture`, `Statistics` for those assays).

- `paper_draft_v2.tex:287-293`  
  Remove the entire Figure 5 legend. Drop Figure 5 from the figure plan.

- `paper_draft_v2.tex:200`  
  After the above cuts, check whether `\cite{declerck2018}` is still needed. If not, remove the unused reference.

## 2. Ready-to-use LaTeX paragraphs

### Intro framing paragraph

```latex
Current single-cell integration benchmarks evaluate success primarily against annotated cell identities and batch-mixing summaries\cite{luecken2022}. This leaves a blind spot for previously unknown rare subpopulations: if a population is absent from available labels, its disappearance cannot be directly registered by label-based metrics. We therefore frame unknown-rare-subpopulation preservation as a methods problem in failure detection under incomplete annotation, rather than as a claim that any specific biological population has already been experimentally established.
```

### Contribution and scope paragraph

```latex
In this work, we introduce RASI as a computational framework for diagnosing and reducing majority-biased distortion during single-cell integration. The contribution is methodological: a rare-aware integration strategy that attenuates correction for cells with weak cross-batch support, together with a label-free neighborhood-preservation statistic that quantifies local structural disruption before and after integration. Across public datasets and atlas scales, these tools improve preservation-oriented computational metrics while maintaining competitive batch mixing, supporting their use as methods for evaluation and mitigation.
```

### Explicit limitation paragraph

```latex
The evidence presented here is computational. It supports the claim that standard workflows can obscure candidate rare states and that selective integration can reduce this effect, but it does not by itself establish the biological reality, lineage status, function, or clinical relevance of any previously unknown population. Experimental confirmation of candidate states remains a necessary downstream step and is outside the scope of the present methods paper.
```

### Optional replacement for the removed Treg subsection

```latex
As a cell-state-level illustration, we tracked a candidate TIGIT$^+$CCR8$^-$ Treg-like cluster through the standard pipeline and through RASI. Under standard preprocessing and Harmony-based integration, the cluster showed substantial fragmentation and loss of neighborhood continuity, whereas RASI preserved a larger fraction of its pre-integration local structure. We present this example as a computational case study of integration behavior in a candidate rare state, not as experimental validation of a new subtype or its function.
```

## 3. Language constraints for a methods-only version

- Use `candidate`, `putative`, `Treg-like`, `computationally defined cluster`, `computational evidence`, `suggests`, and `is consistent with`.
- Do not use `validated`, `confirmed`, `demonstrated biological significance`, `functional community`, `potent immunosuppression`, `therapeutically actionable`, or `previously uncharacterized subtype`.
- State clearly that NP or related metrics detect `local structural disruption after integration`; do not say they detect `true biological loss` by themselves.
- Keep the epistemic boundary explicit: known rare labels can calibrate the method, but unknown rare populations remain hypotheses until experimentally tested.
- Recast downstream impact conservatively: standard integration may `obscure candidate states and alter downstream inference`; it should not be said to reveal a treatment strategy or clinical mechanism in this paper.
- If the manuscript keeps a Treg example, describe it as an `illustrative computational case study` and remove all assay-derived fold changes and mechanistic interpretation.
- A useful one-sentence scope marker for the Introduction or Discussion is: `This manuscript addresses the computational parts of the problem; experimental confirmation of previously unknown candidate populations is deferred to future work.`

## 4. Branch and commit recommendation

- Branch: `philosophy/methods-only-framing-brief`
- Commit: `Add philosophy brief for methods-only manuscript framing`
