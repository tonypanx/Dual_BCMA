# -----------------------------------------------------------------------------
# Per-gene expression, cilta vs ide  (Fig S15b and S15d)
#
# Pseudobulk Wilcoxon on per-patient mean expression, run separately for CD8 and
# CD4 post-infusion CAR+ cells. One violin per gene with the p val ; the
# BH-adjusted values are written to the stats table alongside.
#
# The manuscript panels are cropped from these two multi-gene figures:
#   S15b  CD8: GNLY, NKG7, BAX, CD27, TNFRSF9, FYN
#   S15d  CD4: TNFRSF18, CX3CR1, BAX, FOXO1, ICOS
#
# Input:  Rade et al. integrated T-cell object
# Output: results/key_genes/
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr)
})

source("config.R")
source(file.path("R", "common.R"))
source(file.path("R", "car_gating.R"))

use_arial()

FIG_DIR <- file.path(RESULTS_DIR, "key_genes")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Memory and naive markers, cytotoxic effectors, exhaustion markers, TGF-beta
# receptors, activation and costimulatory receptors, apoptosis
KEY_GENES <- c("MKI67", "SELL", "LEF1", "FOXO1", "IL7R", "CCR7",
               "CXCR3", "FOXP3", "GNLY", "NKG7", "GZMA", "PRF1",
               "HAVCR2", "LAG3", "ENTPD1", "RUNX2", "CCL5", "RUNX3",
               "TGFBR1", "TGFBR2", "B2M", "CD40LG", "CD81", "BAX",
               "BAK1", "TNFRSF4", "TNFRSF18", "ICOS", "CD69", "CD27",
               "TNFRSF9", "FYN", "CX3CR1", "KLRC1", "KLRC2")

product_levels <- c(IDE_VALUE, CILTA_VALUE)
product_colors <- c(ide = "skyblue", cilta = "red")

# --- Pseudobulk --------------------------------------------------------------
make_pseudobulk <- function(seu_obj, genes) {
  gok <- intersect(genes, rownames(seu_obj))
  if (!length(gok)) {
    message("  none of the requested genes are in the object")
    return(NULL)
  }

  meta <- seu_obj@meta.data
  meta <- meta[!is.na(meta[[PRODUCT_COL]]) &
                 meta[[PRODUCT_COL]] %in% product_levels, , drop = FALSE]
  if (!nrow(meta)) return(NULL)
  meta <- keep_min_carpos(meta, patient_col = PATIENT_COL)
  if (!nrow(meta)) return(NULL)

  expr <- GetAssayData(seu_obj[gok, rownames(meta)], assay = "RNA", layer = "data")

  pb <- do.call(rbind, Filter(Negate(is.null), lapply(unique(meta[[PATIENT_COL]]), function(pt) {
    bc <- rownames(meta)[meta[[PATIENT_COL]] == pt]
    grp <- unique(meta[[PRODUCT_COL]][meta[[PATIENT_COL]] == pt])[1]
    if (!length(bc) || is.na(grp)) return(NULL)
    as.data.frame(c(list(patient = pt, group = grp),
                    as.list(rowMeans(as.matrix(expr[, bc, drop = FALSE])))), check.names = FALSE)
  })))
  pb$group <- as.character(pb$group)
  list(pb = pb, genes_ok = gok)
}

run_wilcoxon <- function(pb, genes_ok) {
  deg <- do.call(rbind, Filter(Negate(is.null), lapply(genes_ok, function(g) {
    if (!g %in% colnames(pb)) return(NULL)
    a <- as.numeric(pb[[g]][pb$group == CILTA_VALUE])
    b <- as.numeric(pb[[g]][pb$group == IDE_VALUE])
    wt <- if (length(a) >= 2 && length(b) >= 2)
      tryCatch(wilcox.test(a, b), error = function(e) NULL) else NULL
    data.frame(gene = g,
               log2FC = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
               mean_cilta = mean(a, na.rm = TRUE),
               mean_ide = mean(b, na.rm = TRUE),
               pval = if (!is.null(wt)) wt$p.value else NA_real_,
               n_cilta = length(a),
               n_ide = length(b))
  })))
  deg$padj <- p.adjust(deg$pval, method = "BH")
  deg
}

# --- Plotting -----------------------------------------------------------------
plot_pb_violin <- function(pb, genes_ok, stats_df, title, n_cilta, n_ide, ncols = 6) {
  pb_long <- do.call(rbind, lapply(genes_ok, function(g)
    data.frame(patient = pb$patient,
               group = pb$group,
               gene = g,
               expression = as.numeric(pb[[g]]),
               stringsAsFactors = FALSE)))
  pb_long$group <- factor(pb_long$group, levels = product_levels)
  pb_long$gene <- factor(pb_long$gene,  levels = genes_ok)

  ann_df <- do.call(rbind, lapply(genes_ok, function(g) {
    vals <- pb_long$expression[pb_long$gene == g & !is.na(pb_long$expression)]
    yrange <- max(max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE), 1e-6)
    pv <- stats_df$pval[stats_df$gene == g]
    lbl <- if (!length(pv) || is.na(pv)) "ns" else {
      stars <- sig_stars(pv)
      if (stars == "ns") format_pval(pv) else paste0(stars, " ", format_pval(pv))
    }
    data.frame(gene = g,
               brack_y = max(vals, na.rm = TRUE) + yrange * 0.20,
               tick_h = yrange * 0.06,
               label = lbl,
               stringsAsFactors = FALSE)
  }))
  ann_df$gene <- factor(ann_df$gene, levels = genes_ok)

  ggplot(pb_long, aes(x = group, y = expression, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.85, color = "black",
                linewidth = 0.4, scale = "width") +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA,
                 color = "black", linewidth = 0.4) +
    geom_segment(data = ann_df, aes(x = 1, xend = 2, y = brack_y, yend = brack_y),
                 inherit.aes = FALSE, linewidth = 0.4, color = "black") +
    geom_segment(data = ann_df, aes(x = 1, xend = 1, y = brack_y - tick_h, yend = brack_y),
                 inherit.aes = FALSE, linewidth = 0.4, color = "black") +
    geom_segment(data = ann_df, aes(x = 2, xend = 2, y = brack_y - tick_h, yend = brack_y),
                 inherit.aes = FALSE, linewidth = 0.4, color = "black") +
    geom_text(data = ann_df, aes(x = 1.5, y = brack_y + tick_h * 0.9, label = label),
              inherit.aes = FALSE, size = 3.2, color = "black") +
    scale_fill_manual(values = product_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
    facet_wrap(~gene, scales = "free_y", ncol = ncols) +
    labs(title = title,
         subtitle = paste0("Pseudobulk Wilcoxon | per-patient means | n=", n_cilta, " cilta / ", n_ide, " ide"),
         x = NULL, y = "Expression Level") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8, color = "grey50"),
          strip.text = element_text(face = "bold.italic"),
          axis.text.x = element_text(size = 9))
}

# --- Run ---------------------------------------------------------------------
message("\n-- Loading object --")
seu <- readRDS(EXT_RDS)
seu <- car_subset_obj(seu)
message("Loaded: ", ncol(seu), " cells | ", nrow(seu), " genes")

all_stats <- list()

for (lineage in LINEAGES) {
  message("\n===== ", lineage, " =====")

  bc <- rownames(seu@meta.data)[
    !is.na(seu@meta.data[[T_LIN_COL]]) &
      seu@meta.data[[T_LIN_COL]] == lineage &
      seu@meta.data[[TIMEPOINT_COL]] %in% LATE_VALUES]
  seu_lin <- seu[, bc]
  message("Post-infusion ", lineage, " cells: ", ncol(seu_lin))

  res <- make_pseudobulk(seu_lin, KEY_GENES)
  if (is.null(res)) {
    message("\nno pseudobulk result, skipping")
    next
  }

  n_cilta <- sum(res$pb$group == CILTA_VALUE)
  n_ide   <- sum(res$pb$group == IDE_VALUE)
  message("patients cilta: ", n_cilta, " | ide: ", n_ide)
  if (n_cilta < 2 || n_ide < 2) {
    message("too few patients, skipping")
    next
  }

  stats_df <- run_wilcoxon(res$pb, res$genes_ok)
  message("FDR<0.05: ", sum(!is.na(stats_df$padj) & stats_df$padj < 0.05), "/", nrow(stats_df))
  print(stats_df[, c("gene", "log2FC", "pval", "padj")])

  all_stats[[lineage]] <- cbind(lineage = lineage, stats_df)

  p <- plot_pb_violin(res$pb, res$genes_ok, stats_df, paste0("cilta vs ide - post-inf ", lineage, " (Rade)"), n_cilta, n_ide)
  save_pdf(p, file.path(FIG_DIR, paste0("key_genes_pb_cilta_vs_ide_", lineage, ".pdf")), w = 18, h = 14)
}

if (length(all_stats) > 0) {
  write.csv(do.call(rbind, all_stats), file.path(FIG_DIR, "key_genes_pb_stats.csv"), row.names = FALSE)
  message("\nSaved: key_genes_pb_stats.csv")
}
