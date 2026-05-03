# Methods-Only Revision Summary

## Files

- `biological_science_workspace/paper_draft_v3_methods_only.tex`

## Scope decisions

- Created a new `v3` manuscript instead of modifying `paper_draft_v2.tex`.
- Reframed the paper as a computational methods manuscript centered on RASI, NP, public-data benchmarks, scale sweeps, and ablations.
- Removed the biological-validation section and all assay-based claims from the new draft.
- Kept the immune and tumor examples only as transcript-defined public-atlas case studies with cautious wording.
- Replaced wet-lab Figure 5 with a computational Figure 5 built from existing selective-integration assets: `Fig5_BD_pareto` and `Fig5_BD_integration`.
- Kept only figures that exist in `biological_science_workspace/biological_science_workspace/figures/`.

## Validation performed

- Wet-lab-term audit on `paper_draft_v3_methods_only.tex` using `rg` for prohibited terms returned no matches.
- LaTeX compile succeeded with temporary output:
  - command run from `biological_science_workspace/`
  - output written to `/tmp/rasi_methods_only_GgJvfL/paper_draft_v3_methods_only.pdf`

## Notes

- The draft explicitly limits its claims to public-data and in-silico evidence.
- Batch-count wording was made coherent by using exact counts where stable (`101` for the full atlas benchmark) and broader phrasing for the larger immune scale sweep.
