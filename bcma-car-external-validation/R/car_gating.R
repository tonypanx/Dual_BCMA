# -----------------------------------------------------------------------------
# CAR+ gating...
# The Rade et al. atlas profiles all post-infusion T cells rather than sorted
# CAR T cells, so every comparison here is restricted to transgene-positive
# cells first. That keeps the external cohort on the same footing as the CAR+
# cells in our case study. CAR+ is the atlas authors' own call, CAR_BY_EXPRS == TRUE.
# -----------------------------------------------------------------------------

CAR_COL       <- "CAR_BY_EXPRS" # atlas CAR-by-expression call
CAR_POS_VALUE <- TRUE
CAR_TRANSGENE <- "CAR-BCMA" # backup: transgene feature in the RNA assay

# CAR+ cells, and CAR+ patients, per timepoint and product.
car_report <- function(md) {
  tp <- if ("TIMEPOINT" %in% colnames(md)) md[["TIMEPOINT"]] else rep(NA, nrow(md))
  prd <- if ("PRODUCT"   %in% colnames(md)) md[["PRODUCT"]]   else rep(NA, nrow(md))
  message("  CAR+ cells by TIMEPOINT x PRODUCT:")
  print(table(TIMEPOINT = tp, PRODUCT = prd, useNA = "ifany"))
  if ("PATIENT_ID" %in% colnames(md)) {
    pat <- unique(md[, c("PATIENT_ID", "TIMEPOINT", "PRODUCT")])
    message("  CAR+ patients (>=1 CAR+ cell) by TIMEPOINT x PRODUCT:")
    print(table(TIMEPOINT = pat$TIMEPOINT, PRODUCT = pat$PRODUCT, useNA = "ifany"))
  }
  invisible(NULL)
}

# drop patients contributing fewer than min_cells CAR+ cells to the subset
keep_min_carpos <- function(df, min_cells = CAR_MIN, patient_col = "PATIENT_ID", quiet = FALSE) {
  if (min_cells <= 1) return(df)
  if (!patient_col %in% colnames(df)) {
    if (!quiet) warning("keep_min_carpos: no '", patient_col, "' column; returning df unfiltered")
    return(df)
  }
  tab <- table(df[[patient_col]])
  keep_pt <- names(tab)[tab >= min_cells]
  before <- length(unique(df[[patient_col]]))
  out <- df[df[[patient_col]] %in% keep_pt, , drop = FALSE]
  if (!quiet)
    message(sprintf("  keep_min_carpos: %d/%d patients with >= %d CAR+ cells (%d cells kept)",
                    length(keep_pt), before, min_cells, nrow(out)))
  out
}

# subset
car_subset_obj <- function(seu) {
  md <- seu@meta.data
  if (CAR_COL %in% colnames(md)) {
    keep <- !is.na(md[[CAR_COL]]) & md[[CAR_COL]] == CAR_POS_VALUE
    src <- paste0("column ", CAR_COL)
  } else if (CAR_TRANSGENE %in% rownames(seu)) {
    cnts <- tryCatch(
      GetAssayData(seu, assay = "RNA", layer = "counts")[CAR_TRANSGENE, ],
      error = function(e) GetAssayData(seu, assay = "RNA", slot = "counts")[CAR_TRANSGENE, ])
    keep <- cnts > 0
    src <- paste0("feature ", CAR_TRANSGENE, " > 0")
  } else {
    stop("CAR+ gate: object has neither a '", CAR_COL, "' column nor a '", CAR_TRANSGENE, "' feature.")
  }
  message(sprintf("CAR+ gate [%s]: keeping %d / %d cells (%.1f%%)", src, sum(keep), length(keep), 100 * mean(keep)))
  seu <- subset(seu, cells = colnames(seu)[keep])
  car_report(seu@meta.data)
  seu
}

# same gate, applied to metadata data frame
car_filter_meta <- function(meta) {
  if (CAR_COL %in% colnames(meta)) {
    keep <- !is.na(meta[[CAR_COL]]) & meta[[CAR_COL]] == CAR_POS_VALUE
    src <- paste0("column ", CAR_COL)
  } else {
    if (!file.exists(CAR_BARCODES_RDS))
      stop("CAR+ gate: metadata has no '", CAR_COL, "' column and ",
           CAR_BARCODES_RDS, " does not exist. Run 00_export_car_barcodes.R first.")
    car_pos <- readRDS(CAR_BARCODES_RDS)
    idcol <- if ("cell_barcode" %in% colnames(meta)) meta[["cell_barcode"]] else rownames(meta)
    keep <- idcol %in% car_pos
    src <- paste0("barcode list (", basename(CAR_BARCODES_RDS), ")")
  }
  message(sprintf("CAR+ gate [%s]: keeping %d / %d cells (%.1f%%)", src, sum(keep), length(keep), 100 * mean(keep)))
  meta <- meta[keep, , drop = FALSE]
  car_report(meta)
  meta
}
