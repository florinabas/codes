

library(Seurat)
library(clustree)
library(ClustAssess)

set.seed(42)

#seurat_object <- readRDS("/Users/rafael/Desktop/RA/Flufftail-article/objects/gribben-et-al-end-stage.rds")
full <- readRDS("C:/Users/Floarea/Desktop/corebioinfo/cauzalitate/GSE202379_SeuratObject_AllCells.rds")

cells_keep <- colnames(full)[full$cell.annotation %in% c("Hepatocytes", "Cholangiocytes") &
                               full$Disease.status == "end stage"]

length(cells_keep) #25527

cells_keep <- sample(cells_keep, min(5000, length(cells_keep)))
seurat_object <- subset(full,cells = cells_keep)
#rm(full)

# Plot the UMAP reduction to visualize cell annotations
DimPlot(seurat_object, reduction="umap_harmony_t.0", group.by = "cell.annotation")


# Remove bottom island - all below -6 on UMAP_2
umap_df <- as.data.frame(seurat_object@reductions$umap_harmony_t.0@cell.embeddings)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
#cell_names_keep <- rownames(umap_df[umap_df$UMAP_2 > -6,]) #original
cell_names_keep <- rownames(umap_df[umap_df$UMAP_1 > 4 & umap_df$UMAP_2 > -10, ])
seurat_object <- subset(seurat_object, cells = cell_names_keep)
DimPlot(seurat_object, reduction="umap_harmony_t.0", group.by = "cell.annotation")

# Get sample names
sample_names <- colnames(seurat_object)



DefaultAssay(seurat_object) <- "RNA"

# Neighbors
seurat_object <- FindNeighbors(
  seurat_object,
  reduction = "harmony_t.0",
  dims = 1:30
)

# Multiple resolutions
resolutions <- seq(0.1, 1.2, by = 0.1)

for(res in resolutions){
  
  seurat_object <- FindClusters(
    seurat_object,
    graph.name = "SCT_snn",
    resolution = res,
    algorithm = 1
  )
  
}


# =========================================================
# 1. Extract clustering columns
# =========================================================

cluster_cols <- grep(
  "SCT_snn_res",
  colnames(seurat_object@meta.data),
  value = TRUE
)

cluster_cols

# =========================================================
# 2. Build dataframe for ClustAssess
# =========================================================

cluster_df <- seurat_object@meta.data[, cluster_cols]

head(cluster_df)

resolution_vals <- as.numeric(
  gsub("SCT_snn_res\\.", "", cluster_cols)
)

resolution_vals <- sort(resolution_vals)


# =========================================================
# 3. Run clustering stability assessment
# =========================================================



graph <- seurat_object@graphs$SCT_snn


ca <- assess_clustering_stability(
  graph_adjacency_matrix = graph,
  resolution = resolution_vals,
  n_repetitions = 10   # reduce for speed; use 100 later if needed
)


plot_clustering_overall_stability(ca)
plot_clustering_per_value_stability(ca)

Idents(seurat_object) <- "SCT_snn_res.0.5"

DimPlot(seurat_object, group.by = "SCT_snn_res.0.5")


chosen_res <- "SCT_snn_res.0.5"
Idents(seurat_object) <- chosen_res


markers <- FindAllMarkers(
  seurat_object,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


FeaturePlot(
  seurat_object,
  features = c(
    "ALB", "HNF4A", "APOA1",   # hepatocyte program
    "KRT19", "KRT7", "SOX9"    # cholangiocyte program
  )
)

cluster_summary <- markers %>%
  group_by(cluster) %>%
  summarise(top_genes = paste(head(gene, 10), collapse = ", "))
start_cluster <- "2"  # example only

start_genes <- markers %>%
  filter(cluster == start_cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  pull(gene) %>%
  head(20)


end_cluster <- "7"  # example only

end_genes <- markers %>%
  filter(cluster == end_cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  pull(gene) %>%
  head(20)


FeaturePlot(
  seurat_object,
  features = c(
    "ALB", "HNF4A",
    "KRT19", "SOX9"
  )
)


## data driven
avg_expr <- AverageExpression(seurat_object, group.by = chosen_res)$RNA


start_cluster <- colnames(avg_expr)[
  which.max(avg_expr["ALB", ] + avg_expr["HNF4A", ])
]

end_cluster <- colnames(avg_expr)[
  which.max(avg_expr["KRT19", ] + avg_expr["SOX9", ])
]


start_genes  # hepatocyte program
end_genes    # cholangiocyte program


## shiny
