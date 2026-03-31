# Paper Outline — Nature Submission

## Title
Neighborhood preservation reveals hidden casualties of single-cell batch integration

## Abstract (~150 words)
- Problem: batch integration can destroy unknown rare subpopulations with no existing metric detecting the loss
- Key insight: structural change (neighbor preservation) detects destruction without labels
- Method: NP (Neighborhood Preservation) score — fully unsupervised, AUC=0.837
- Scale: 2.26M cells, 103 batches, 30 cancer types, 47 minutes
- Finding: systematic subpopulation disruption across all 8 immune cell types; pDC most affected (NP=0.231)
- Biological validation: GITR⁺TIGIT⁺CCR8⁻ Treg functional community disrupted (survival=0.39), wet lab confirms immunosuppressive function
- Protection: NP-Guard adaptive integration strategy

## Main Figures
1. Concept + problem definition (three-way deadlock)
2. Three-generation metric evolution (CNEM→R(S)→NP)
3. Full-scale NP results (2.26M cells) + scale-dependent vulnerability
4. Biological validation (Cluster 5 GITR⁺ Treg + wet lab)
5. NP-Guard protection framework

## Sections
- Introduction: paradigm blind spot, not just "methods not good enough"
- Results: NP development, full-scale validation, biological case studies
- Discussion: epistemological shift, clinical implications, limitations
- Methods: NP algorithm, experimental setup, wet lab protocols

## Key Numbers
- NP AUC: 0.837 (Harmony), 0.728 (Scanorama)
- Full-scale: 2,256,276 cells, 103 batches, 47 minutes
- pDC NP: 0.231 (most affected), Mast NP: 0.241
- 55% cells NP<0.3 (high risk), 21.7% NP<0.1 (extreme risk)
- Cluster 5: 1464 cells, GITR=69.8%, TIGIT=65%, survival=0.39
- Triple combo (FOXP3⁺GITR⁺TIGIT⁺CCR8⁻): 266 cells
- Wet lab: CD69↓4x, IFN-γ↓2.6x, killing↓3.3x (all p<0.001)
- Metric evolution: CNEM AUC=0.40 → R(S) r=-0.95 → NP AUC=0.84

## Status
- Materials collected from all 7 agents
- Ready for drafting
