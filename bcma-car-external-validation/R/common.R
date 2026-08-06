# -----------------------------------------------------------------------------
# Metadata column names, cohort definitions and small shared helpers.
# The names below are the ones used in the Rade object; they are pulled
# out here so a change to the upstream schema only has to be made once
# -----------------------------------------------------------------------------

PATIENT_COL <- "PATIENT_ID"
PRODUCT_COL <- "PRODUCT"
TIMEPOINT_COL <- "TIMEPOINT"
RESPONSE_COL <- "BEST_RESPONSE"
T_LIN_COL <- "T_LIN"
SUBTYPE_COL <- "celltype"

CILTA_VALUE <- "cilta"
IDE_VALUE <- "ide"

# "Post-infusion" throughout the paper means the Late and Very Late timepoints
LATE_VALUES <- c("Late", "Very Late")

LINEAGES <- c("CD8", "CD4")

# Response grouping. Best response of CR or VGPR counts as a responder
RESPONDER_VALUES <- c("CR", "VGPR")
NONRESPONDER_VALUES <- c("PR", "SD", "PD")

add_response_group <- function(md) {
  md$response_group <- ifelse(
    md[[RESPONSE_COL]] %in% RESPONDER_VALUES, "Responder",
    ifelse(md[[RESPONSE_COL]] %in% NONRESPONDER_VALUES, "Non-Responder", NA_character_))
  md
}

# ggplot to pdf
save_pdf <- function(p, path, w = 12, h = 8) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  pdf(path, width = w, height = h)
  print(p)
  dev.off()
  message("Saved: ", path)
}

# significance stars, usual thresholds
sig_stars <- function(p) {
  if (is.na(p) || p >= 0.05) return("ns")
  if (p < 1e-4) return("****")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  "*"
}

format_pval <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "p<0.001", paste0("p=", signif(p, 2))))
}

use_arial <- function() {
  if (!requireNamespace("showtext", quietly = TRUE)) return(invisible(FALSE))
  candidates <- c(
    "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    "/usr/share/fonts/liberation/LiberationSans-Regular.ttf",
    "/Library/Fonts/Arial.ttf",
    "/System/Library/Fonts/Supplemental/Arial.ttf"
  )
  fp <- Filter(file.exists, candidates)
  if (!length(fp)) {
    message("Arial not found on this system; using the default sans font")
    return(invisible(FALSE))
  }
  showtext::font_add("Arial", regular = fp[1])
  showtext::showtext_auto()
  message("Font loaded: ", fp[1])
  invisible(TRUE)
}
