"""Step 1: Download and prepare human pancreas datasets for PoC experiment."""
import os, warnings
import numpy as np
warnings.filterwarnings('ignore')

os.environ['HTTP_PROXY'] = 'http://127.0.0.1:7890'
os.environ['HTTPS_PROXY'] = 'http://127.0.0.1:7890'

output_dir = '/ssd/data/agent/bio/eacn_example_001/experiments/poc_v1/data'
os.makedirs(output_dir, exist_ok=True)

import pooch, anndata as ad, scanpy as sc

# scIB pancreas benchmark dataset (Baron+Muraro+Segerstolpe+Wang)
url = "https://figshare.com/ndownloader/files/24539828"
print("Downloading pancreas dataset from figshare...")
fname = pooch.retrieve(url, known_hash=None, fname="pancreas.h5ad", path=output_dir)

adata = ad.read_h5ad(fname)
print(f"Shape: {adata.shape}")
print(f"Obs columns: {list(adata.obs.columns)}")
for col in adata.obs.columns:
    n = adata.obs[col].nunique()
    if n < 50:
        print(f"\n{col} ({n} unique):")
        print(adata.obs[col].value_counts())
print("\nDone!")
