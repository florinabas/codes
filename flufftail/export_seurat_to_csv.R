# ── 1. Load Required Libraries ───────────────────────────────────────────────
library(Seurat)
library(org.Hs.eg.db)
library(scCustomize) # Used for robust, seamless .h5ad export

# ── 2. Load the Master Seurat Object ──────────────────────────────────────────
input_dir   <- "C:/Users/Floarea/Desktop/corebioinfo/codes/data_processing/input_data"
seurat_path <- file.path(input_dir, "GSE202379_SeuratObject_AllCells.rds")

seurat_object <- readRDS(seurat_path)
print(paste("Original Dataset:", nrow(seurat_object), "genes ×", ncol(seurat_object), "cells"))

# ── 3. Identify Protein-Coding Genes ──────────────────────────────────────────
all_genes <- rownames(seurat_object)

# Query the database for gene types
gene_mapping <- select(org.Hs.eg.db, 
                       keys = all_genes, 
                       columns = c("SYMBOL", "GENETYPE"), 
                       keytype = "SYMBOL")

# Extract only genes explicitly flagged as "protein-coding"
protein_coding_genes <- gene_mapping$SYMBOL[gene_mapping$GENETYPE == "protein-coding"]
protein_coding_genes <- protein_coding_genes[!is.na(protein_coding_genes)]

# ── 4. Filter the Seurat Object ───────────────────────────────────────────────
full_protein_coding <- subset(seurat_object, features = protein_coding_genes)
print(paste("Filtered Dataset:", nrow(full_protein_coding), "genes ×", ncol(full_protein_coding), "cells"))

# ── 5. OUTPUT 1: Save as R Data Structure (.rds) ──────────────────────────────
rds_output <- file.path(input_dir, "full_protein_coding.rds")
saveRDS(full_protein_coding, file = rds_output)
print(paste("📂 [1/2] R-ready file saved to:", rds_output))

# ── 6. OUTPUT 2: Save as Python Data Structure (.h5ad) ────────────────────────
# scCustomize handles the conversion of internal metadata (obs) and counts (X) 
# smoothly without breaking sparse matrix indexing.
print("Converting and writing .h5ad file for Python...")
as.anndata(
  x = full_protein_coding, 
  file_path = input_dir, 
  file_name = "full_protein_coding.h5ad",
  assay = "RNA" # Speeds up export by saving your standard raw/log count matrices
)
print(paste("📂 [2/2] Python-ready file saved to:", file.path(input_dir, "full_protein_coding.h5ad")))
print("🎉 Success! Both files are perfectly identical and saved.")