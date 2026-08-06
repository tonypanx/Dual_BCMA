# -----------------------------------------------------------------------------
# RUNX3 and KLF2 regulon scores per patient, cilta vs ide  (Fig 4g)
#
# Reads the per-cell module scores from 01_signature_scoring.R, averages each
# regulon score to one value per patient, and tests cilta against ide with a
# two-sided Wilcoxon. CD8 and CD4 are tested separately.
# The whole test is run twice, at a minimum of 10 and 20 CAR+ cells per patient.
# Ten is what the figure reports; twenty is the sensitivity check, since results
# at this sample size can be sensitive to where the threshold sits.
#
# Input:  interim/rade_meta_scored_carpos.rds
# Output: results/regulon/
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork)
})

source("config.R")
source(file.path("R", "common.R"))
source(file.path("R", "car_gating.R"))

use_arial()
theme_update(text = element_text(family = "Arial"))

FIG_DIR <- file.path(RESULTS_DIR, "regulon")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_LEVELS <- c(10, 20)

SIG_LBLS <- c(RUNX3_regulon_score = "RUNX3 Regulon",
              KLF2_regulon_score  = "KLF2 Regulon")

product_colors <- c(ide = "skyblue", cilta = "red")

make_plabel <- function(p) {
  if (is.na(p)) return("")
  s <- sig_stars(p)
  v <- if (p < 1e-3) "p<0.001" else paste0("p=", signif(p, 2))
  if (s == "ns") v else paste0(s, " ", v)
}

# --- Load --------------------------------------------------------------------
message("Loading: ", SCORED_META_RDS)
meta <- readRDS(SCORED_META_RDS)
meta <- car_filter_meta(meta)
meta <- meta[meta[[TIMEPOINT_COL]] %in% LATE_VALUES & !is.na(meta[[T_LIN_COL]]), ]

SIGS <- intersect(names(SIG_LBLS), colnames(meta))
if (length(SIGS) < length(SIG_LBLS))
  stop("Missing score columns: ",
       paste(setdiff(names(SIG_LBLS), SIGS), collapse = ", "),
       ". re-run 01_signature_scoring.R...")

# --- One test, One panel -----------------------------------------------------
pb_panel <- function(d, sig, min_cells, lineage) {
  d <- d[!is.na(d[[PRODUCT_COL]]) &
           d[[PRODUCT_COL]] %in% c(CILTA_VALUE, IDE_VALUE) & !is.na(d[[sig]]), ]
  d <- keep_min_carpos(d, min_cells = min_cells, patient_col = PATIENT_COL, quiet = TRUE)

  blank <- data.frame(min_carpos = min_cells, lineage = lineage,
                      signature = SIG_LBLS[[sig]], n_cilta = 0, n_ide = 0,
                      mean_cilta = NA_real_, mean_ide = NA_real_,
                      direction = NA_character_, p_value = NA_real_,
                      stringsAsFactors = FALSE)
  if (!nrow(d)) return(list(plot = NULL, row = blank))

  pb <- d %>%
    group_by(.data[[PATIENT_COL]], .data[[PRODUCT_COL]]) %>%
    summarize(m = mean(.data[[sig]], na.rm = TRUE), .groups = "drop")
  colnames(pb) <- c("patient", "grp", "m")
  pb$grp <- factor(pb$grp, levels = c(CILTA_VALUE, IDE_VALUE))

  v1 <- pb$m[pb$grp == CILTA_VALUE]
  v2 <- pb$m[pb$grp == IDE_VALUE]
  row <- blank
  row$n_cilta <- length(v1)
  row$n_ide <- length(v2)
  row$mean_cilta <- if (length(v1)) mean(v1) else NA_real_
  row$mean_ide <- if (length(v2)) mean(v2) else NA_real_
  if (length(v1) < 2 || length(v2) < 2) return(list(plot = NULL, row = row))

  p <- wilcox.test(v1, v2)$p.value
  row$p_value <- p
  row$direction <- if (mean(v1) > mean(v2)) "cilta higher" else "ide higher"

  pb$grp <- factor(as.character(pb$grp), levels = c(IDE_VALUE, CILTA_VALUE))
  ymax <- max(pb$m, na.rm = TRUE)
  yr <- max(ymax - min(pb$m, na.rm = TRUE), 1e-6)

  g <- ggplot(pb, aes(grp, m, fill = grp)) +
    geom_violin(trim = TRUE, scale = "width", alpha = .85, color = "black", linewidth = .4) +
    geom_boxplot(width = .12, fill = "white", outlier.shape = NA, linewidth = .4) +
    geom_jitter(width = .06, size = .8, alpha = .6, color = "grey20") +
    annotate("segment", x = 1, xend = 2, y = ymax + yr * .18, yend = ymax + yr * .18,
             linewidth = .4) +
    annotate("text", x = 1.5, y = ymax + yr * .26, label = make_plabel(p), size = 3) +
    scale_x_discrete(labels = c(ide = "Ide", cilta = "Cilta")) +
    scale_fill_manual(values = product_colors) +
    scale_y_continuous(expand = expansion(mult = c(.05, .30))) +
    labs(title = paste0(SIG_LBLS[[sig]], " - ", lineage, " (CAR+, >=", min_cells, ")"),
         subtitle = paste0("cilta n=", length(v1), " | ide n=", length(v2), " pts"),
         x = NULL, y = "Mean score per patient") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8.5, face = "bold"),
          plot.subtitle = element_text(size = 7, color = "grey45"))
  list(plot = g, row = row)
}

# --- Sweep -------------------------------------------------------------------
all_rows <- list()

for (mc in MIN_LEVELS) {
  message("\n===== minimum ", mc, " CAR+ cells per patient =====")
  panels  <- list()
  n_real  <- 0
  for (lin in LINEAGES) {
    L <- meta[meta[[T_LIN_COL]] == lin, ]
    for (sig in SIGS) {
      r <- pb_panel(L, sig, mc, lin)
      all_rows[[length(all_rows) + 1]] <- r$row
      if (is.null(r$plot)) {
        panels[[paste(lin, sig)]] <- patchwork::plot_spacer()
      } else {
        panels[[paste(lin, sig)]] <- r$plot
        n_real <- n_real + 1
      }
    }
  }

  if (n_real == 0) {
    message("\nno panel had at least 2 patients per arm, figure skipped")
    next
  }

  save_pdf(
    wrap_plots(panels, ncol = length(SIGS)) +
      plot_annotation(
        title = paste0("Regulon scores, cilta vs ide (>=", mc, " CAR+ cells/patient)"),
        subtitle = "CD8 top, CD4 bottom | per-patient means, two-sided Wilcoxon",
        theme = theme(plot.title = element_text(size = 12, face = "bold"),
                      plot.subtitle = element_text(size = 9, color = "grey45"))),
    file.path(FIG_DIR, paste0("regulon_pseudobulk_min", mc, ".pdf")),
    w = 2.7 * length(SIGS), h = 6.4)
}

stats <- do.call(rbind, all_rows)
write.csv(stats, file.path(FIG_DIR, "regulon_pseudobulk_stats.csv"), row.names = FALSE)
message("\nSaved stats: ", file.path(FIG_DIR, "regulon_pseudobulk_stats.csv"),
        " (", nrow(stats), " tests)")
print(subset(stats, min_carpos == 10))
