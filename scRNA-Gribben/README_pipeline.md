# scRNA-Gribben GRN pipeline — condensed

Parsimonious, freshly-commented rewrite of the pipeline. Same logic and same
decision flow as the full version, with exploratory/plotting cruft and old
comments stripped. Flufftail (`03`) is left exactly as-is.

Keep every file in one folder; run R from there (`source("gribben_config.R")`
resolves) and Jupyter/python from there (`from gribben_config import *` resolves).

## Run order

| step | file | format | does | reads → writes |
|------|------|--------|------|----------------|
| — | `gribben_config.py` / `.R` | module | single source of truth (paths, gene sets, `STABLE_TRAIN_CONFIG`, `CHOSEN_LAM`) | — |
| 00 | `00_prepare_data.Rmd` | R notebook | end-stage hep+chol, protein-coding | raw rds → prepared `.rds`/`.h5ad` |
| 01 | `01_clustassess_hvg.Rmd` | R notebook | pick #HVGs by feature stability | prepared rds → `clustassess_chosen_hvg.txt` |
| 02 | `02_noisyr_clean.Rmd` | R notebook | noise threshold → signal genes | prepared rds → `noisyr_signal_genes.txt` |
| 03 | `03_flufftail_gribben.Rmd` | R notebook | Flufftail tutorial (**untouched**) | prepared rds → `flufftail_*` hand-off |
| 04 | `04_velorama_grn.ipynb` | Jupyter | build curated h5ad + definitive causal GRN + null | prepared h5ad + HVG count → `liver_endstage.h5ad`, GRN CSVs |
| 05 | `05_velorama_stability.ipynb` | Jupyter | seed stability (curate→converge→lam sweep→verdict→ablation→noisyR) | curated h5ad + signal genes → consensus GRN, metrics |
| 06 | `06_pseudotime_reconciliation.ipynb` | Jupyter | DPT vs Monocle3 from a shared root (trajectory gate) | prepared h5ad → agreement summary |
| 07 | `07_flufftail_velorama_grn.ipynb` | Jupyter | stage-specific directional GRNs; do Flufftail hubs survive? | curated h5ad + Flufftail hand-off → stage GRNs, hub-influence |

R stages (`00`–`03`) are **.Rmd** — open in RStudio or run via `rmarkdown::render()`.
Python stages (`04`–`07`) are **.ipynb** — open in Jupyter/JupyterLab.

## Why this order

`04` trains in the stable regime by construction (config-driven). `05` certifies
that regime is seed-robust. `06` confirms Velorama's DPT axis and Flufftail's
Monocle3 axis are the same trajectory — the gate that licenses `07`'s comparison.
`07` then builds directional per-stage GRNs and tests whether Flufftail's
co-variation hubs re-emerge as causal, lag-resolved regulators.

## Format note

Stages 04–07 are `.ipynb` and 00–03 are `.Rmd`, matching the language each stage
is written in. Each `# === STEP ===` / `## === STEP ===` marker is its own
cell/chunk. Every file ends with a short `NEXT STEPS / GAPS` note on open
judgment calls (confirm `CHOSEN_LAM`, stricter/per-stage nulls, full cell-cycle
gene lists, the background-gene universe choice, optionally driving the DAG
with Monocle3).

## Visualizations

Every step that transforms the data now has a chart, plot, or table right after
it, so you can see the effect before moving on — e.g. in `00`: whole-object UMAP
→ after cell-type/disease-stage subsetting → the island-removal cutoff drawn on
the UMAP → after removal → gene-count bar chart for the protein-coding filter →
final composition table. Similar before/after visuals run through every stage:
HVG mean-variance plots and a decision UMAP in `01`; noise-threshold histograms
and signal/removed bars in `02`; QC histograms, cell-cycle scatter, HVG plot,
gene-role bar/pie charts, DPT-on-UMAP, GRN heatmaps, and null-vs-real histograms
in `04`; lam-sweep curves, stability verdict bars, and ablation/noisyR
comparison bars in `05`; the existing 3-panel DPT-vs-Monocle3 comparison plus a
summary bar in `06`; and stage-assignment histograms, GRN heatmaps, and
significant-edge-per-stage bars in `07`.
