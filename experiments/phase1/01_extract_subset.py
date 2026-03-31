"""Phase 1: Extract immune atlas subset for CNEM validation."""
import anndata as ad
import numpy as np
import os, time

OUT_DIR = '/ssd/data/agent/bio/eacn_example_001/experiments/phase1/data'
os.makedirs(OUT_DIR, exist_ok=True)

SELECTED = ['D5', 'D100', 'D38', 'D88', 'D83', 'D94', 'D102', 'D80']

print("Loading full atlas...")
t0 = time.time()
adata = ad.read_h5ad('/ssd/data/agent/bio/atlas_merged_immune.h5ad')
print(f"Loaded in {time.time()-t0:.0f}s: {adata.shape}")

print(f"\nSubsetting to {len(SELECTED)} datasets: {SELECTED}")
mask = adata.obs['Dataset'].isin(SELECTED)
adata_sub = adata[mask].copy()
print(f"Subset: {adata_sub.shape}")

# Summary
print("\nCelltype distribution:")
print(adata_sub.obs['Celltype'].value_counts())
print(f"\nDatasets: {adata_sub.obs['Dataset'].nunique()}")
print(f"Tissues: {adata_sub.obs['Tissue'].value_counts().to_dict()}")

# Save
out_path = os.path.join(OUT_DIR, 'immune_subset_phase1.h5ad')
adata_sub.write(out_path)
print(f"\nSaved to {out_path}")
print(f"File size: {os.path.getsize(out_path)/1e9:.1f} GB")
