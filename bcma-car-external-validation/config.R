# -----------------------------------------------------------------------------
# Paths and analysis constants.
# DATA_DIR => Rade et al. Seurat object and the signature .rds files
# -----------------------------------------------------------------------------

# Root folder for the input data
DATA_DIR <- Sys.getenv("RADE_DATA_DIR", unset = "data")

# Rade et al. integrated T cell seurat obj (~1.8 GB)
EXT_RDS <- file.path(DATA_DIR, "06_seurat_harmony_t_all_new.Rds")

# Folder of gene-signature .rds files (from Tony)
SIG_DIR <- file.path(DATA_DIR, "Signatures")

# Intermediates written by 00_ and 01_ and read by later scripts
INTERIM_DIR <- "interim"

# Figures and stats tables
RESULTS_DIR <- "results"

dir.create(INTERIM_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# CAR+ barcode list written by 00_export_car_barcodes
CAR_BARCODES_RDS <- file.path(INTERIM_DIR, "car_pos_barcodes.rds")

# Per-cell module scores written by 01_signature_scoring.R
SCORED_META_RDS <- file.path(INTERIM_DIR, "rade_meta_scored_carpos.rds")

# Minimum CAR+ cells a patient must contribute to a given comparison before that patient's mean enters a pseudobulk test. 
# Patients with a handful of CAR+ cells give very noisy means.
# 10 is the threshold used throughout the paper; 
# 20 is reported alongside it as a sensitivity check in 03_regulon_pseudobulk.R
CAR_MIN <- 10
