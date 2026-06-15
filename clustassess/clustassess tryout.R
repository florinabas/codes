rm(list = ls())

##### Step 0: Load libraries #####
library(Seurat)
library(SeuratDisk)    # convert h5ad/loom -> h5Seurat
library(hdf5r)
library(rhdf5)
library(dplyr)
library(Matrix)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(remotes)
library(data.table)  # for fwrite
#####


##### Step 1: Read H5AD file #####
path <- "C:/Users/Floarea/Desktop/app bioinfo/FGFR local/input data"
h5ad_file <- file.path(path, "MCA1.1_adata.h5ad")
adata <- H5File$new(h5ad_file, mode = "r")
h5ls(h5ad_file) # X: 34947 genes x 333778 cells

# Cell names (columns of X)
cell_names <- adata[["obs"]]$read()[["index"]]

# Find cell names that contain "lung" (case-insensitive)
lung_names <- cell_names[grepl("lung", cell_names, ignore.case = TRUE)]

# Convert to dataframe (helps with matching)
cell_names_df <- data.frame(cell_name = cell_names, stringsAsFactors = FALSE)

# Column indices of lung cells
lung_cols <- which(cell_names_df$cell_name %in% lung_names)

# Read gene names and make a dataframe
gene_names <- adata[["var"]]$read()[["index"]]
gene_names_df <- data.frame(gene_name = gene_names, stringsAsFactors = FALSE)

##### Step 2: Identify FGFR genes 
# Detect any occurrence of "fgfr" (case-insensitive) anywhere in the string
fgfr_genes_all <- gene_names[grepl("fgfr", gene_names, ignore.case = TRUE)]
print(fgfr_genes_all)
# output: "Fgfr1"      "Fgfr1op"    "Fgfr1op2"   "Fgfr2"      "Fgfr3"      "Fgfrl1"     "Fgfr4"      "Fgfr3-ps"   "Fgfr3-ps.1"

## focus on the classic gene names also used in the study trying to replicate
fgfr_genes <- c("Fgfr1", "Fgfr2", "Fgfr3", "Fgfr4")  # mouse names
fgfr_rows <- which(gene_names_df$gene_name %in% fgfr_genes)
# Fgfr genes (rows) 4802  4805  4806 23156
fgfr_genes_df <- gene_names_df[fgfr_rows, , drop = FALSE]
message("FGFR genes present in data: ", paste(fgfr_genes_df$gene_name, collapse = ", "))

##### Step 3: Subset X to FGFR genes and FGFR+ cells #####
# Full expression matrix
X <- adata[["X"]]

# Only FGFR rows
fgfr_expr <- X[fgfr_rows, lung_cols]


#### check how many tissues left 
# Identify cells where expression > 1 for any FGFR gene
cells_high_expr <- which(colSums(fgfr_expr > 1) >= 1)

# Extract corresponding cell types
cell_types_high_expr <- sub("_.*$", "", cell_names[cells_high_expr])


# Identify cells expressing at least one FGFR gene
fgfr_expr_per_cell <- colSums(fgfr_expr > 0)
fgfr_cells <- which(fgfr_expr_per_cell >= 1)
length(fgfr_cells)  # number of FGFR+ cells

# Keep all metadata for these cells
obs_df_fgfr <- adata[["obs"]]$read()[fgfr_cells, , drop = FALSE]
obs_df_fgfr$cell_type <- sub("_.*$", "", obs_df_fgfr$index)  # keep cell type column

# Subset expression matrix to FGFR genes x FGFR+ cells
fgfr_counts_small <- fgfr_expr[, fgfr_cells]
rownames(fgfr_counts_small) <- fgfr_genes_df$gene_name

##### Step 4: Create Seurat object #####
seurat_fgfr <- CreateSeuratObject(
  counts = fgfr_counts_small,
  meta.data = obs_df_fgfr
)

# Normalize and scale
seurat_fgfr <- NormalizeData(seurat_fgfr)
seurat_fgfr <- ScaleData(seurat_fgfr)

##### Step 5: PCA -----
seurat_fgfr <- RunPCA(seurat_fgfr, features = fgfr_genes, npcs = 4)

# Add tiny noise to PCA embeddings to avoid identical points
pca_data <- Embeddings(seurat_fgfr, "pca")
pca_data <- pca_data + matrix(rnorm(length(pca_data), sd = 1e-6), nrow = nrow(pca_data))
seurat_fgfr[["pca"]] <- CreateDimReducObject(embeddings = pca_data, key = "PC_", assay = "RNA")

# Quick check: Fgfr1 vs Fgfr2
FeatureScatter(seurat_fgfr, feature1 = "Fgfr1", feature2 = "Fgfr2", slot="count", raster = FALSE)

##### Step 6: UMAP -----
num_pcs <- ncol(Embeddings(seurat_fgfr, "pca"))
seurat_fgfr <- RunUMAP(seurat_fgfr, dims = 1:num_pcs)

# DimPlot colored by cell type (preserved from metadata)
DimPlot(seurat_fgfr, reduction = "umap", group.by = "cell_type", pt.size = 0.5)

# Extract embeddings for ggplot
umap_df <- as.data.frame(seurat_fgfr@reductions$umap@cell.embeddings)
umap_df$cell_type <- seurat_fgfr$cell_type

# Plot with ggplot
#if(!is.null(dev.list())) dev.off()
#dev.new()
ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cell_type)) +
  geom_point(size = 0.5, alpha = 0.6) +
  theme_minimal() +
  labs(title = "UMAP of FGFR+ cells colored by cell type")

##### Highlight Fgfr1 expression 
Fgfr1_expr <- GetAssayData(seurat_fgfr, layer = "data")["Fgfr1", ]
seurat_fgfr$Fgfr1_positive <- Fgfr1_expr > 0

DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "Fgfr1_positive",
  cols = c("grey80", "red"),
  pt.size = 0.5
) + ggtitle("Cells expressing Fgfr1")


##### Step 6A: Graph construction #####

seurat_fgfr <- FindNeighbors(
  seurat_fgfr,
  reduction = "pca",
  dims = 1:num_pcs,
  k.param = 20
)

##### Step 6B: Clustering #####

seurat_fgfr <- FindClusters(
  seurat_fgfr,
  resolution = 0.5
)

##### Visualize clusters #####
DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
)

#### Clustassess #####

#remotes::install_github("Core-Bioinformatics/ClustAssess")
library(ClustAssess)

##### Step 7: Prepare matrix for ClustAssess #####

expr_matrix <- as.matrix(
  GetAssayData(seurat_fgfr, layer = "data")
)

# ClustAssess expects:
# features = rows
# cells = columns

dim(expr_matrix)


##### Step 8: Run ClustAssess #####
options(scipen = 999)

library(ClustAssess)
library(RANN)

pca_embeddings <- Embeddings(seurat_fgfr, "pca")

nn_result <- RANN::nn2(
  pca_embeddings,
  k = 5
)

adj_matrix <- getNNmatrix(
  nn_result$nn.idx,
  k = 5,
  prune = 0
)$nn

rownames(adj_matrix) <- rownames(pca_embeddings)
colnames(adj_matrix) <- rownames(pca_embeddings)

clust_stability <- assess_clustering_stability(
  
  graph_adjacency_matrix = adj_matrix,
  
  resolution = c(0.2, 0.4, 0.6, 0.8),
  
  n_repetitions = 10
)

names(clust_stability)
str(clust_stability, max.level = 2)

plot_clustering_overall_stability(clust_stability)
plot_clustering_per_value_stability(
  clust_stability,
  value_type = "resolution"
)

best_clusters <- clust_stability$all$best_partition


seurat_fgfr$stable_cluster <- best_clusters

DimPlot(
  seurat_fgfr,
  group.by = "stable_cluster",
  label = TRUE
)


#####
############################################################
################### 0. LOAD LIBRARIES ######################
############################################################

library(Seurat)
library(hdf5r)
library(rhdf5)
library(Matrix)
library(ggplot2)
library(RANN)
library(ClustAssess)

############################################################
############### 1. LOAD H5AD SINGLE-CELL DATA ##############
############################################################

# Path to dataset
path <- "C:/Users/Floarea/Desktop/app bioinfo/FGFR local/input data"

# Full h5ad file
h5ad_file <- file.path(path, "MCA1.1_adata.h5ad")

# Open H5AD object
adata <- H5File$new(h5ad_file, mode = "r")

# Optional: inspect H5AD structure
h5ls(h5ad_file)

############################################################
############### 2. SELECT LUNG CELLS #######################
############################################################

gene_names <- adata[["var"]]$read()[["index"]]
gene_names_df <- data.frame(gene_name = gene_names, stringsAsFactors = FALSE)


# Extract cell names from metadata
cell_names <- adata[["obs"]]$read()[["index"]]

# Keep only cells containing "lung"
lung_cols <- grep(
  "lung",
  cell_names,
  ignore.case = TRUE
)

# Subsample 5000 lung cells for faster analysis
set.seed(123)
#fgfr_genes_df <- gene_names_df[fgfr_rows, , drop = FALSE]
fgfr_genes <- c("Fgfr1", "Fgfr2", "Fgfr3", "Fgfr4")  # mouse names
fgfr_rows <- which(gene_names_df$gene_name %in% fgfr_genes)
# Fgfr genes (rows) 4802  4805  4806 23156

#message("FGFR genes present in data: ", paste(fgfr_genes_df$gene_name, collapse = ", "))

##### Step 3: Subset X to FGFR genes and FGFR+ cells #####
# Full expression matrix
#X <- adata[["X"]]

# Only FGFR rows
#fgfr_expr <- X[fgfr_rows, lung_cols]
############################################################
############### 3. EXTRACT EXPRESSION MATRIX ###############
############################################################

# Extract expression matrix
# rows = genes
# columns = cells
expr <- adata[["X"]][fgfr_rows , lung_cols]

# Convert to standard matrix
counts <- as.matrix(expr)

############################################################
############### 4. CREATE SEURAT OBJECT ####################
############################################################

# Build Seurat object
seurat_object <- CreateSeuratObject(
  counts = counts,
  project = "lung_cols"
)

############################################################
############### 5. STANDARD PREPROCESSING ##################
############################################################

# Normalize counts
seurat_object <- NormalizeData(seurat_object)

# Detect highly variable genes
seurat_object <- FindVariableFeatures(
  seurat_object,
  nfeatures = 3000
)

# Scale expression values
seurat_object <- ScaleData(seurat_object)

############################################################
#################### 6. PCA REDUCTION ######################
############################################################

# Run PCA on variable genes
seurat_object <- RunPCA(
  seurat_object,
  npcs = 30
)

# Optional: inspect explained variance
ElbowPlot(seurat_object)

############################################################
################### 7. BUILD KNN GRAPH #####################
############################################################

# Build nearest-neighbor graph
# This graph is the basis for Louvain clustering
seurat_object <- FindNeighbors(
  seurat_object,
  dims = 1:20,
  k.param = 20
)

############################################################
################### 8. INITIAL CLUSTERING ##################
############################################################

# Initial Louvain clustering
seurat_object <- FindClusters(
  seurat_object,
  resolution = 0.5
)

############################################################
################### 9. COMPUTE UMAP ########################
############################################################

# Compute UMAP visualization
seurat_object <- RunUMAP(
  seurat_object,
  dims = 1:20
)

# Visualize initial clusters
DimPlot(
  seurat_object,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
)

############################################################
############ 10. EXTRACT PCA EMBEDDINGS ####################
############################################################

# PCA coordinates for each cell
pca_embeddings <- Embeddings(
  seurat_object,
  "pca"
)

############################################################
############### 11. BUILD ADJACENCY MATRIX #################
############################################################

# Compute nearest neighbors directly from PCA space
nn_result <- RANN::nn2(
  pca_embeddings,
  k = 10
)

# Convert neighbor indices into adjacency matrix
adj_matrix <- getNNmatrix(
  nn_result$nn.idx,
  k = 10,
  prune = 0
)$nn

# Add cell names to adjacency matrix
rownames(adj_matrix) <- rownames(pca_embeddings)
colnames(adj_matrix) <- rownames(pca_embeddings)

############################################################
############### 12. RUN CLUSTASSESS ########################
############################################################

# Prevent scientific notation issues
options(scipen = 999)

# Assess clustering stability
clust_stability <- assess_clustering_stability(
  
  # Graph used for clustering
  graph_adjacency_matrix = adj_matrix,
  
  # Louvain resolutions to test
  resolution = c(0.2, 0.4, 0.6, 0.8),
  
  # Number of repeated clustering runs
  n_repetitions = 20
)

############################################################
############### 13. INSPECT OUTPUT #########################
############################################################

# View structure of stability object
names(clust_stability)

str(clust_stability, max.level = 2)

############################################################
############### 14. VISUALIZE STABILITY ####################
############################################################

# Overall clustering stability
plot_clustering_overall_stability(
  clust_stability
)

# Stability per resolution value
plot_clustering_per_value_stability(
  clust_stability,
  value_type = "resolution"
)

############################################################
############### 15. EXTRACT BEST CLUSTERS ##################
############################################################

# Best partition chosen by ClustAssess
best_clusters <- clust_stability$all$best_partition

# Add stable clusters to metadata
seurat_object$stable_cluster <- best_clusters

############################################################
############### 16. VISUALIZE STABLE CLUSTERS ##############
############################################################

DimPlot(
  seurat_object,
  reduction = "umap",
  group.by = "stable_cluster",
  label = TRUE
)

############################################################
############### 17. IDENTIFY ROOT/END STATES ###############
############################################################

# Explore stable cluster sizes
table(seurat_object$stable_cluster)

############################################################
############### 18. BIOLOGICAL MARKER GENES ################
############################################################

# Visualize biologically meaningful genes
# Example markers:
# Alb = hepatocytes
# Krt19 = cholangiocytes

FeaturePlot(
  seurat_object,
  features = c("Alb", "Krt19")
)

############################################################
############### 19. DEFINE PSEUDOTIME ROOT #################
############################################################

# Example:
# Suppose cluster 1 expresses progenitor markers
root_cells <- Cells(seurat_object)[
  seurat_object$stable_cluster == 1
]

############################################################
############### 20. DEFINE TERMINAL STATE ##################
############################################################

# Example:
# Suppose cluster 5 expresses mature markers
terminal_cells <- Cells(seurat_object)[
  seurat_object$stable_cluster == 5
]

############################################################
############### 21. SAVE OBJECT ############################
############################################################

saveRDS(
  seurat_object,
  "clustassess_lung_analysis.rds"
)
