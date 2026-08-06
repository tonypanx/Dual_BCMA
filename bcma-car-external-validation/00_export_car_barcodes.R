# -----------------------------------------------------------------------------
# Write out the CAR+ cell barcodes.
# Also prints CAR+ cell and patient counts per timepoint and product, which are
# the sample sizes behind the pseudobulk tests.
#
# Input:  Rade et al. integrated T-cell object
# Output: interim/car_pos_barcodes.rds
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(Seurat))

source("config.R")
source(file.path("R", "common.R"))
source(file.path("R", "car_gating.R"))

message("Loading object: ", EXT_RDS)
seu <- readRDS(EXT_RDS)

md <- seu@meta.data
if (CAR_COL %in% colnames(md)) {
  keep <- !is.na(md[[CAR_COL]]) & md[[CAR_COL]] == CAR_POS_VALUE
} else {
  cnts <- tryCatch(
    GetAssayData(seu, assay = "RNA", layer = "counts")[CAR_TRANSGENE, ],
    error = function(e) GetAssayData(seu, assay = "RNA", slot = "counts")[CAR_TRANSGENE, ])
  keep <- cnts > 0
}

car_pos_barcodes <- colnames(seu)[keep]
saveRDS(car_pos_barcodes, CAR_BARCODES_RDS)
message("Saved ", length(car_pos_barcodes), " CAR+ barcodes -> ", CAR_BARCODES_RDS)

message("\nCAR+ power table:")
car_report(md[keep, , drop = FALSE])
