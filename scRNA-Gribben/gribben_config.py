"""
gribben_config.py — single source of truth for the Python side of the
hepatocyte→cholangiocyte GRN pipeline (Gribben et al. end-stage liver).

Every notebook does `from gribben_config import *` so paths, gene sets, the lag,
the curated-regulator cap, the chosen lam and the *stable* training regime are
defined exactly once. Keep this file in the same folder as the notebooks and
launch Jupyter from that folder.

The numbers here encode pipeline DECISIONS made by earlier stages:
  - N_TOP_HVGS        ← ClustAssess        (step 01, via clustassess_chosen_hvg.txt)
  - STABLE_TRAIN_CONFIG, CAP_REGULATORS, CHOSEN_LAM ← seed-stability stage (step 05)
"""
import os

# ── project root & directories ───────────────────────────────────────────────
ROOT        = r"C:\Users\Floarea\Desktop\analyses\scRNA-Gribben"
INPUT_DIR   = os.path.join(ROOT, "0_input_data")
INTERM_DIR  = os.path.join(ROOT, "0_interm_data")
RESULTS_DIR = os.path.join(ROOT, "99_results")
VELO_DIR        = os.path.join(RESULTS_DIR, "velorama")              # steps 04, 05*
FLUFF_VELO_DIR  = os.path.join(RESULTS_DIR, "flufftail_velorama")    # step 07
PSEUDO_DIR      = os.path.join(RESULTS_DIR, "pseudotime_comparison") # step 06
for _d in (INTERM_DIR, VELO_DIR, FLUFF_VELO_DIR, PSEUDO_DIR):
    os.makedirs(_d, exist_ok=True)

# ── dataset names & the two AnnData files ────────────────────────────────────
# PREPARED: all protein-coding genes + Harmony embeddings (step 00). Used by the
#           Velorama BUILD half (step 04) and the pseudotime reconciliation (06).
# CURATED : reg/target/background subset, scaled, carries dpt_pseudotime / X_pca /
#           iroot and the velorama_*_names lists in .uns (written by step 04).
#           Used by the stability notebooks (05*) and the integration (07).
DATASET_PREPARED = "hep_chol_endstate_truncated_protein_coding"
DATASET_CURATED  = "liver_endstage"

PREPARED_RDS  = os.path.join(INPUT_DIR, f"{DATASET_PREPARED}.rds")
PREPARED_H5AD = os.path.join(INPUT_DIR, f"{DATASET_PREPARED}.h5ad")
CURATED_H5AD  = os.path.join(VELO_DIR,  f"{DATASET_CURATED}.h5ad")
TF_LIST_PATH  = os.path.join(INPUT_DIR, "allTFs_hg38.txt")

# ── decision artifacts handed between stages ─────────────────────────────────
HVG_FILE         = os.path.join(INTERM_DIR, "clustassess_chosen_hvg.txt")   # 01 → 04
SIGNAL_GENES_TXT = os.path.join(INTERM_DIR, "noisyr_signal_genes.txt")      # 02 → 05c
FLUFF_PSEUDOTIME = os.path.join(INTERM_DIR, "flufftail_pseudotime.csv")     # 03 → 06, 07
FLUFF_STATES     = os.path.join(INTERM_DIR, "flufftail_cell_states.csv")    # 03 → 07
FLUFF_DE_GENES   = os.path.join(INTERM_DIR, "flufftail_de_genes.txt")       # 03 → 07
FLUFF_HUBS       = os.path.join(INTERM_DIR, "flufftail_candidate_hubs.txt") # 03 → 07
FLUFF_DIR        = INTERM_DIR   # alias the integration notebook expects

# ── biology: lineage TFs + Flufftail-nominated hubs (defined once) ───────────
HEPATO_TFs    = {"HNF4A", "HNF1A", "FOXA2", "CEBPA", "NR5A2", "NR1H4"}
CHOLANGIO_TFs = {"SOX9", "HNF1B", "FOXA1", "GATA4", "GATA6", "SP1"}
LINEAGE_TFS   = HEPATO_TFs | CHOLANGIO_TFs
FLUFFTAIL_HUBS       = ["KLF6", "SERPINE1", "CREB5", "FGF13"]
FLUFFTAIL_HUBS_TF    = ["KLF6", "CREB5"]        # used as Velorama regulators
FLUFFTAIL_HUBS_NONTF = ["SERPINE1", "FGF13"]    # kept as candidate targets only

# ── Velorama: shared knobs + the regime proven stable in step 05 ─────────────
LAG                = 2
N_BACKGROUND_GENES = 500
CAP_REGULATORS     = 150     # curated regulator universe (variance-ranked + must-keep)
CHOSEN_LAM         = 0.05    # confirm/adjust via step 05a lam sweep; re-run 04 & 07 if changed

def n_top_hvgs(default=1500):
    """HVG count chosen by ClustAssess (step 01); falls back to `default`."""
    try:
        return int(open(HVG_FILE).read().strip())
    except OSError:
        return default

# Convergence-aware, full-cell, curated-regulator regime validated in step 05.
# Consumed by step 04 (definitive GRN) and step 07 (integration) so the headline
# result and its validation share one configuration. lam / seed / results_dir /
# dir_name are overridden per call.
STABLE_TRAIN_CONFIG = dict(
    reg_target=False, lr=0.005, lam=CHOSEN_LAM, lam_ridge=1e-4, penalty="H",
    lag=LAG, hidden=[16], max_iter=1000, device="cpu", lookback=20,
    check_every=25, verbose=False, dynamics="pseudotime",
)

# ── compute ──────────────────────────────────────────────────────────────────
N_CPUS  = max(1, (os.cpu_count() or 4) - 1)
WORKERS = max(1, min(4, N_CPUS))
