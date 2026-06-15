#!/usr/bin/env python
"""
Run Velorama analysis on a Seurat object
Converts Seurat RDS -> AnnData -> Velorama analysis
"""

import os
import sys
import numpy as np
import pandas as pd
import torch
import scanpy as sc
import subprocess
import tempfile

# Configuration
SEURAT_FILE = r"C:\Users\Floarea\Desktop\corebioinfo\cauzalitate\GSE202379_SeuratObject_AllCells.rds"
OUTPUT_DIR = r"C:\Users\Floarea\Desktop\corebioinfo\velorama_results"
H5AD_FILE = r"C:\Users\Floarea\Desktop\corebioinfo\cauzalitate\GSE202379_AnnData.h5ad"

def convert_seurat_to_h5ad(rds_path, h5ad_path):
    """Convert Seurat RDS to H5AD using R"""
    
    if os.path.exists(h5ad_path):
        print(f"✓ H5AD file already exists: {h5ad_path}")
        return h5ad_path
    
    print(f"Converting Seurat object to H5AD...")
    print(f"  Input: {rds_path}")
    print(f"  Output: {h5ad_path}")
    
    # Create R script
    r_script = f"""
library(Seurat)
suppressWarnings(library(SeuratDisk))

print("Loading Seurat object...")
seurat_obj <- readRDS('{rds_path}')

print(paste("Object dimensions:", nrow(seurat_obj), "x", ncol(seurat_obj)))
print(paste("Assays:", paste(names(seurat_obj@assays), collapse=", ")))

# Check for RNA velocity data
if ("velocity" %in% names(seurat_obj@assays)) {{
    print("✓ Found velocity data")
}}

print("Converting to H5AD format...")
tryCatch({{
    SaveH5Seurat(seurat_obj, filename='{h5ad_path}', overwrite=TRUE)
    print("✓ Conversion successful!")
}}, error = function(e) {{
    # Fallback: export manually
    print("Note: SeuratDisk not available, exporting manually...")
    
    # Get expression matrix
    X <- seurat_obj@assays[["RNA"]]@counts
    if (length(X@x) == 0) {{
        X <- seurat_obj@assays[["RNA"]]@data
    }}
    X <- as.matrix(X)
    
    # Get metadata
    metadata <- seurat_obj@meta.data
    
    # Write to CSV temporarily
    write.csv(t(X), file='{tempfile.gettempdir()}/expression.csv')
    write.csv(metadata, file='{tempfile.gettempdir()}/metadata.csv')
    print("Exported to CSV for Python conversion")
}})
"""
    
    # Write and execute R script
    with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
        f.write(r_script)
        script_path = f.name
    
    try:
        print("  Running R conversion...")
        result = subprocess.run(['Rscript', script_path], 
                              capture_output=True, text=True, timeout=600)
        print(result.stdout)
        if result.stderr:
            print("  Warnings:", result.stderr)
        
        if os.path.exists(h5ad_path):
            print(f"✓ Successfully created: {h5ad_path}")
            return h5ad_path
        else:
            print("✗ H5AD file not created")
            return None
    finally:
        os.remove(script_path)

def load_and_explore_data(h5ad_path):
    """Load and explore the AnnData object"""
    
    print(f"\nLoading AnnData object...")
    adata = sc.read(h5ad_path)
    
    print(f"\n{'='*60}")
    print(f"DATASET SUMMARY")
    print(f"{'='*60}")
    print(f"Shape: {adata.shape} (cells × genes)")
    print(f"\nMetadata columns: {list(adata.obs.columns)}")
    print(f"\nVar columns: {list(adata.var.columns)}")
    print(f"\nLayers: {list(adata.layers.keys())}")
    print(f"Obsm: {list(adata.obsm.keys())}")
    print(f"Varm: {list(adata.varm.keys())}")
    print(f"Obsp: {list(adata.obsp.keys())}")
    
    # Show sample metadata
    print(f"\nSample metadata (first 5 cells):")
    print(adata.obs.head())
    
    return adata

def prepare_for_velorama(adata, output_h5ad):
    """Prepare AnnData object for Velorama analysis"""
    
    print(f"\n{'='*60}")
    print(f"PREPARING DATA FOR VELORAMA")
    print(f"{'='*60}")
    
    # Normalize data
    print("Normalizing data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Select highly variable genes
    print("Finding highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable]
    
    print(f"Selected {adata.n_vars} highly variable genes")
    
    # PCA
    print("Computing PCA...")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50)
    
    # Neighbors
    print("Computing neighbors...")
    sc.pp.neighbors(adata, n_neighbors=30, use_rep='X_pca')
    
    # Get cell cycle phase if available
    if 'phase' in adata.obs.columns:
        print(f"✓ Cell cycle phase information available")
    
    # Check for pseudotime
    if 'pseudotime' in adata.obs.columns or 'Pseudotime' in adata.obs.columns:
        print(f"✓ Pseudotime information available")
        pt_col = 'pseudotime' if 'pseudotime' in adata.obs.columns else 'Pseudotime'
        print(f"  Using column: {pt_col}")
    else:
        print("⚠ No pseudotime column found")
        print("  Computing pseudotime using DPT...")
        sc.tl.dpt(adata)
        print(f"  DPT distances computed and saved to adata.obs['dpt_pseudotime']")
    
    # Check for velocity information
    if 'velocity' in adata.layers:
        print(f"✓ RNA velocity information available in layers")
    elif 'spliced' in adata.layers and 'unspliced' in adata.layers:
        print(f"✓ Spliced/unspliced counts available - can compute velocity")
    else:
        print("⚠ No velocity information available")
        print("  Velorama will use pseudotime mode instead")
    
    # Mark all genes as potential regulators and targets
    # (Can be refined later based on biological knowledge)
    print("\nMarking regulators and targets...")
    adata.var['is_reg'] = True
    adata.var['is_target'] = True
    print(f"  All {adata.n_vars} genes marked as potential regulators and targets")
    
    # Save prepared object
    print(f"\nSaving prepared AnnData object...")
    adata.write_h5ad(output_h5ad)
    print(f"✓ Saved to: {output_h5ad}")
    
    return adata

def print_velorama_commands(data_dir, dataset_name):
    """Print commands to run Velorama"""
    
    print(f"\n{'='*60}")
    print(f"NEXT STEPS: RUN VELORAMA")
    print(f"{'='*60}")
    
    print(f"\nTo run Velorama with pseudotime:")
    print(f"  python -m velorama \\")
    print(f"    -n GSE202379 \\")
    print(f"    -ds {dataset_name} \\")
    print(f"    -dyn pseudotime \\")
    print(f"    -dev cpu \\")
    print(f"    -cp 2 \\")
    print(f"    -l 5 \\")
    print(f"    -hd 32 \\")
    print(f"    -rd {data_dir} \\")
    print(f"    -sd ./results")
    
    print(f"\nOr with RNA velocity (if available):")
    print(f"  python -m velorama \\")
    print(f"    -n GSE202379 \\")
    print(f"    -ds {dataset_name} \\")
    print(f"    -dyn rna_velocity \\")
    print(f"    -dev cpu \\")
    print(f"    -cp 2 \\")
    print(f"    -l 5 \\")
    print(f"    -hd 32 \\")
    print(f"    -rd {data_dir} \\")
    print(f"    -sd ./results")

def main():
    print(f"{'='*60}")
    print(f"VELORAMA ANALYSIS PIPELINE")
    print(f"{'='*60}")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Step 1: Convert Seurat to H5AD
    print(f"\n[STEP 1/3] Converting Seurat to H5AD")
    h5ad_file = convert_seurat_to_h5ad(SEURAT_FILE, H5AD_FILE)
    
    if not h5ad_file or not os.path.exists(h5ad_file):
        print("✗ Failed to convert Seurat object")
        sys.exit(1)
    
    # Step 2: Load and explore
    print(f"\n[STEP 2/3] Exploring data")
    adata = load_and_explore_data(h5ad_file)
    
    # Step 3: Prepare for Velorama
    print(f"\n[STEP 3/3] Preparing for Velorama")
    prepared_h5ad = h5ad_file.replace('.h5ad', '_prepared.h5ad')
    adata = prepare_for_velorama(adata, prepared_h5ad)
    
    # Print next steps
    data_dir = os.path.dirname(prepared_h5ad)
    dataset_name = os.path.basename(prepared_h5ad).replace('_prepared.h5ad', '')
    print_velorama_commands(data_dir, dataset_name)
    
    print(f"\n{'='*60}")
    print(f"✓ PREPARATION COMPLETE!")
    print(f"{'='*60}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
