## gribben_config.R — shared paths + gene sets for the R side of the pipeline.
## Source this at the top of every R script / Rmd:  source("gribben_config.R")
## Keep it in the same folder as the scripts and run R from that folder.

ROOT        <- "C:/Users/Floarea/Desktop/analyses/scRNA-Gribben"
INPUT_DIR   <- file.path(ROOT, "0_input_data")
INTERM_DIR  <- file.path(ROOT, "0_interm_data")
RESULTS_DIR <- file.path(ROOT, "99_results")
dir.create(INTERM_DIR, showWarnings = FALSE, recursive = TRUE)

## datasets / artifacts
DATASET_PREPARED <- "hep_chol_endstate_truncated_protein_coding"
RAW_RDS       <- file.path(INPUT_DIR, "GSE202379_SeuratObject_AllCells.rds")
PREPARED_RDS  <- file.path(INPUT_DIR, paste0(DATASET_PREPARED, ".rds"))
PREPARED_H5AD <- file.path(INPUT_DIR, paste0(DATASET_PREPARED, ".h5ad"))

## decision artifact: ClustAssess (step 01) -> Velorama build (step 04)
HVG_FILE <- file.path(INTERM_DIR, "clustassess_chosen_hvg.txt")

## biology (mirror of gribben_config.py)
HEPATO_TFs     <- c("HNF4A", "HNF1A", "FOXA2", "CEBPA", "NR5A2", "NR1H4")
CHOLANGIO_TFs  <- c("SOX9", "HNF1B", "FOXA1", "GATA4", "GATA6", "SP1")
FLUFFTAIL_HUBS <- c("KLF6", "SERPINE1", "CREB5", "FGF13")
