#!/usr/bin/env python
"""
Convert Seurat object to AnnData format for Velorama analysis
"""

import subprocess
import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import tempfile
import json

def convert_seurat_to_anndata(seurat_path, output_path=None):
    """
    Convert a Seurat RDS file to AnnData h5ad format
    
    Parameters:
    -----------
    seurat_path : str
        Path to the Seurat RDS file
    output_path : str, optional
        Path where to save the h5ad file. If None, will save in same dir as input
    
    Returns:
    --------
    str : Path to the saved h5ad file
    """
    
    if output_path is None:
        output_path = seurat_path.replace('.rds', '.h5ad')
    
    print(f"Converting Seurat object from: {seurat_path}")
    print(f"Output will be saved to: {output_path}")
    
    # Create R script to convert Seurat to h5ad
    r_script = f"""
# Load required libraries
library(Seurat)
library(SeuratDisk)

# Read Seurat object
print("Loading Seurat object...")
seurat_obj <- readRDS("{seurat_path}")

# Convert to h5ad format
print("Converting to h5ad format...")
SaveH5Seurat(seurat_obj, filename="{output_path}", overwrite=TRUE)
print("Conversion complete!")
"""
    
    # Write R script to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
        f.write(r_script)
        r_script_path = f.name
    
    try:
        # Run R script
        print("\nRunning R conversion script...")
        result = subprocess.run(
            ['Rscript', r_script_path],
            capture_output=True,
            text=True,
            timeout=300
        )
        
        if result.returncode != 0:
            print(f"Error running R script:")
            print(result.stderr)
            return None
        
        print(result.stdout)
        
        if os.path.exists(output_path):
            print(f"\n✓ Successfully converted to: {output_path}")
            return output_path
        else:
            print(f"Error: Output file not created")
            return None
            
    finally:
        # Clean up temp file
        if os.path.exists(r_script_path):
            os.remove(r_script_path)

if __name__ == "__main__":
    seurat_file = "C:/Users/Floarea/Desktop/corebioinfo/cauzalitate/GSE202379_SeuratObject_AllCells.rds"
    
    if not os.path.exists(seurat_file):
        print(f"Error: Seurat file not found at {seurat_file}")
        sys.exit(1)
    
    # Convert
    h5ad_file = convert_seurat_to_anndata(seurat_file)
    
    if h5ad_file:
        print(f"\nNext step: Load this file in Python and prepare for Velorama:")
        print(f"  import scanpy as sc")
        print(f"  adata = sc.read('{h5ad_file}')")
