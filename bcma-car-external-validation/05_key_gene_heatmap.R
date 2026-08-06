# -----------------------------------------------------------------------------
# Group-average expression heatmap (Fig 4h)
# Per-patient, per-lineage pseudobulk means for curated gene panel, z-scored
# per gene across all patient columns, then averaged within each of the four
# product-by-lineage groups. Z-scoring before averaging keeps the gradient that
# would otherwise be flattened by averaging first.
#
# Input:  Rade et al. integrated T-cell object
# Output: results/heatmap/key_genes_avg_4group.pdf
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(scales)
})

source("config.R")
source(file.path("R", "common.R"))
source(file.path("R", "car_gating.R"))

FIG_DIR <- file.path(RESULTS_DIR, "heatmap")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# raw order
KEY_GENES <- c(
  "FOXO1", "CD27", "LEF1", "TNFRSF9", "IL7R", "FYN",
  "RUNX3", "TGFBR1", "TGFBR2", "RUNX2", "SMAD2",
  "B2M", "CD81", "CD40LG",
  "IFNG", "GZMA", "GNLY", "GZMB", "NKG7",
  "BAX", "LGALS3", "FAS"
)

PATHWAYS <- list(
  "Memory/Costim" = c("FOXO1", "CD27", "LEF1", "TNFRSF9", "IL7R", "FYN"),
  "TGF-b Signaling" = c("RUNX3", "TGFBR1", "TGFBR2", "RUNX2", "SMAD2"),
  "T cell Activation" = c("B2M", "CD81", "CD40LG"),
  "T cell Cytotoxicity" = c("IFNG", "GZMA", "GNLY", "GZMB", "NKG7"),
  "Apoptosis" = c("BAX", "LGALS3", "FAS")
)

PATHWAY_COLORS <- c(
  "Memory/Costim" = "#E67E22",
  "TGF-b Signaling" = "#2C3E50",
  "T cell Activation" = "#1565C0",
  "T cell Cytotoxicity" = "#00897B",
  "Apoptosis" = "#7B1FA2"
)

PRODUCT_FILL <- c(cilta = "#F1948A", ide = "#AED6F1")
PRODUCT_TEXT <- c(cilta = "#922B21", ide = "#1A5276")
PRODUCT_LABEL <- c(cilta = "Cilta", ide = "Ide")

GROUP_LEVELS <- c("cilta_CD8", "cilta_CD4", "ide_CD8", "ide_CD4")

# --- Load --------------------------------------------------------------------
message("Loading object...")
seu <- readRDS(EXT_RDS)
seu <- car_subset_obj(seu)
message("Loaded: ", ncol(seu), " cells | ", nrow(seu), " genes")

genes_ok <- intersect(KEY_GENES, rownames(seu))
message("Genes present: ", length(genes_ok), "/", length(KEY_GENES))
if (length(genes_ok) < length(KEY_GENES))
  message(" missing: ", paste(setdiff(KEY_GENES, genes_ok), collapse = ", "))

bc_use <- rownames(seu@meta.data)[
  !is.na(seu@meta.data[[T_LIN_COL]]) &
    seu@meta.data[[T_LIN_COL]] %in% LINEAGES &
    seu@meta.data[[TIMEPOINT_COL]] %in% LATE_VALUES &
    !is.na(seu@meta.data[[PRODUCT_COL]])]
message("Post-infusion CD8 + CD4 cells: ", length(bc_use))

seu_sub <- seu[genes_ok, bc_use]
meta_sub <- seu_sub@meta.data
expr_mat <- GetAssayData(seu_sub, assay = "RNA", layer = "data")

# --- Pseudobulk per patient and lineage --------------------------------------
meta_sub$pt_lin <- paste0(meta_sub[[PATIENT_COL]], "_", meta_sub[[T_LIN_COL]])
meta_sub <- keep_min_carpos(meta_sub, patient_col = "pt_lin")
expr_mat <- expr_mat[, rownames(meta_sub), drop = FALSE]

pt_lin_ids <- unique(meta_sub$pt_lin)
pb_mat <- do.call(cbind, lapply(pt_lin_ids, function(pl) {
  bc <- rownames(meta_sub)[meta_sub$pt_lin == pl]
  rowMeans(as.matrix(expr_mat[, bc, drop = FALSE]))
}))
colnames(pb_mat) <- pt_lin_ids

ann_col <- do.call(rbind, lapply(pt_lin_ids, function(pl) {
  rows <- meta_sub[meta_sub$pt_lin == pl, ]
  data.frame(Product = unique(rows[[PRODUCT_COL]])[1],
             Lineage = unique(rows[[T_LIN_COL]])[1],
             row.names = pl, stringsAsFactors = FALSE)
}))
ann_col$gkey <- paste0(ann_col$Product, "_", ann_col$Lineage)

gene_order <- KEY_GENES[KEY_GENES %in% genes_ok]
pb_mat <- pb_mat[gene_order, ]
pb_z <- t(scale(t(pb_mat)))

pb_avg4 <- do.call(cbind, lapply(GROUP_LEVELS, function(g) {
  idx <- which(ann_col$gkey == g)
  if (!length(idx)) return(rep(NA_real_, length(gene_order)))
  rowMeans(pb_z[, idx, drop = FALSE], na.rm = TRUE)
}))
colnames(pb_avg4) <- GROUP_LEVELS
rownames(pb_avg4) <- gene_order

message("Patients per group: ", paste(GROUP_LEVELS, "=", as.integer(table(
  factor(ann_col$gkey, levels = GROUP_LEVELS))), collapse = "  "))

# --- Heatmap -----------------------------------------------------------------
make_heatmap <- function(mat, product_of_col, col_xlabels, title) {
  gene_df <- do.call(rbind, lapply(names(PATHWAYS), function(pw) {
    g <- intersect(PATHWAYS[[pw]], rownames(mat))
    if (!length(g)) return(NULL)
    data.frame(gene = g, pathway = pw, stringsAsFactors = FALSE)
  }))
  gene_df <- gene_df[order(match(gene_df$gene, rownames(mat))), , drop = FALSE]

  # walk down the rows, insert gap btwn pathways
  ROW_GAP <- 0.55
  y <- nrow(gene_df)
  prev_pw <- NULL
  y_vec <- numeric(nrow(gene_df))
  for (i in seq_len(nrow(gene_df))) {
    if (!is.null(prev_pw) && gene_df$pathway[i] != prev_pw) y <- y - ROW_GAP
    y_vec[i] <- y
    y <- y - 1
    prev_pw <- gene_df$pathway[i]
  }
  gene_df$y <- y_vec

  # walk down the rows, insert gap btwn products
  COL_GAP <- 0.55
  cols <- colnames(mat)
  x_vec <- numeric(length(cols))
  x <- 1
  for (j in seq_along(cols)) {
    if (j > 1 && product_of_col[[cols[j]]] != product_of_col[[cols[j - 1]]])
      x <- x + COL_GAP
    x_vec[j] <- x
    x <- x + 1
  }
  col_df <- data.frame(col = cols, x = x_vec, label = col_xlabels[cols],
                       product = product_of_col[cols], stringsAsFactors = FALSE)

  df <- do.call(rbind, lapply(gene_df$gene, function(g)
    do.call(rbind, lapply(seq_len(nrow(col_df)), function(j)
      data.frame(gene = g, x = col_df$x[j],
                 y = gene_df$y[gene_df$gene == g],
                 z = mat[g, col_df$col[j]], stringsAsFactors = FALSE)))))
  
  zlim <- min(max(abs(df$z), na.rm = TRUE), 2.5)
  inner_brks <- c(-1, 0, 1)
  cb_breaks <- sort(unique(c(-zlim, inner_brks[abs(inner_brks) < zlim - 0.05], zlim)))

  y_hdr_lo <- max(gene_df$y) + 0.58
  y_hdr_hi <- max(gene_df$y) + 1.08

  X_BAR_L <- max(col_df$x) + 0.85
  X_BAR_R <- X_BAR_L + 0.22
  X_TXT <- X_BAR_R + 0.20

  ann_df <- do.call(rbind, lapply(names(PATHWAYS), function(pw) {
    ys <- gene_df$y[gene_df$pathway == pw]
    if (!length(ys)) return(NULL)
    data.frame(pw = pw, ymin = min(ys) - 0.43, ymax = max(ys) + 0.43,
               ymid = mean(ys), stringsAsFactors = FALSE)
  }))

  p <- ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_tile(color = "black", linewidth = 0.28, width = 0.97, height = 0.97)

  for (pr in unique(col_df$product)) {
    xs <- col_df$x[col_df$product == pr]
    p <- p +
      annotate("rect", xmin = min(xs) - 0.44, xmax = max(xs) + 0.44,
               ymin = y_hdr_lo, ymax = y_hdr_hi,
               fill = PRODUCT_FILL[[pr]], color = "white", linewidth = 0.5) +
      annotate("text", x = mean(xs), y = (y_hdr_lo + y_hdr_hi) / 2,
               label = PRODUCT_LABEL[[pr]], fontface = "bold", size = 3.5,
               color = PRODUCT_TEXT[[pr]])
  }

  for (i in seq_len(nrow(ann_df))) {
    p <- p +
      annotate("rect", xmin = X_BAR_L, xmax = X_BAR_R,
               ymin = ann_df$ymin[i], ymax = ann_df$ymax[i],
               fill = PATHWAY_COLORS[[ann_df$pw[i]]], color = NA) +
      annotate("text", x = X_TXT, y = ann_df$ymid[i], label = ann_df$pw[i],
               hjust = 0, size = 2.9, color = "grey15")
  }

  p +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-zlim, zlim), oob = squish, name = "Z-score",
      breaks = cb_breaks, labels = sprintf("%.1f", cb_breaks),
      guide = guide_colorbar(barheight = unit(2.5, "cm"),
                             barwidth  = unit(0.32, "cm"),
                             title.position = "top", title.hjust = 0.5)) +
    scale_x_continuous(breaks = col_df$x, labels = col_df$label,
                       limits = c(min(col_df$x) - 0.55, X_TXT + 2.8), expand = c(0, 0)) +
    scale_y_continuous(breaks = gene_df$y, labels = gene_df$gene,
                       limits = c(min(gene_df$y) - 0.55, y_hdr_hi + 0.15), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_text(face = "italic", size = 9, hjust = 1, color = "black"),
          axis.text.x = element_text(size = 9, color = "black", vjust = 0.5),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 10),
          legend.position = "right",
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 7),
          plot.margin = margin(8, 5, 8, 5))
}

prod4 <- setNames(c("cilta", "cilta", "ide", "ide"), GROUP_LEVELS)
xlbl4 <- setNames(c("CD8", "CD4", "CD8", "CD4"), GROUP_LEVELS)

save_pdf(make_heatmap(pb_avg4, prod4, xlbl4, "Key gene expression: group average"),
         file.path(FIG_DIR, "key_genes_avg_4group.pdf"), w = 8, h = 11)
