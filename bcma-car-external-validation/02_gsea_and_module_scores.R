# -----------------------------------------------------------------------------
# GSEA and curated module scores on CAR+ post-infusion T cells.
#
# Produces:
#   Fig. 4e: enrichment curves for the Kaech naive-vs-day-8-effector set, CD8 and CD4, cilta vs ide
#   Fig. 4f: per-patient module scores, CD8 memory and CD8 cytotoxic effector
#   Fig. S15a: enrichment curves, CD8: RUNX3 regulon, KLF2 regulon, IFN-gamma
#   Fig. S15c: enrichment curves, CD4: RUNX3 regulon, KLF2 regulon, IFN-gamma,
#             inflammatory response
#
# Genes are ranked by sign(patient-level log2FC) * -log10(Wilcoxon p) on
# per-patient means, then fed to fgseaMultilevel. Ranking on per-patient means
# rather than per-cell keeps the test at the level the cohort is powered at
#
# Input:  Rade et al. integrated T-cell object, Signatures/
# Output: results/gsea/
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
  library(fgsea);  library(msigdbr)
})

source("config.R")
source(file.path("R", "common.R"))
source(file.path("R", "car_gating.R"))

FIG_DIR <- file.path(RESULTS_DIR, "gsea")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

GSEA_MIN_SIZE <- 10

source_colors <- c(Signature = "#C0392B", MSigDB = "#2980B9", Category = "#27AE60")
product_colors <- c(ide = "#2980B9", cilta = "#C0392B")

# --- Pathway sets ------------------------------------------------------------
message("\n-- Building pathway gene sets --")

pathway_sets <- list()
pathway_source <- character(0)
display_names <- character(0)

add_set <- function(name, genes, source, display = NULL) {
  genes <- unique(as.character(genes[!is.na(genes) & nzchar(genes)]))
  if (length(genes) < 3) return(invisible(NULL))
  pathway_sets[[name]] <<- genes
  pathway_source[[name]] <<- source
  display_names[[name]] <<- if (is.null(display)) gsub("_", " ", name) else display
}

extract_genes <- function(obj) {
  if (is.null(obj)) return(character(0))
  if (is.character(obj)) return(obj)
  if (is.data.frame(obj) && "gene" %in% colnames(obj)) return(obj$gene)
  if (is.list(obj) && length(obj) > 0) return(obj[[1]])
  character(0)
}

if (dir.exists(SIG_DIR)) {
  rds_files <- list.files(SIG_DIR, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
  raw_sigs <- setNames(lapply(rds_files, readRDS), tools::file_path_sans_ext(basename(rds_files)))
  g <- function(k) extract_genes(raw_sigs[[k]])
  gene_signatures <- list(
    CAR_TRM = head(g("CAR_TRM_gene_set"), 500),
    RUNX3_regulon = g("RUNX3_Regulon"),
    KLF2_regulon = g("KLF2_Regulon"),
    TGFb_signaling = g("TGFB"),
    Hallmark_IL2 = g("Hallmark_IL2"),
    IL2_up = g("IL2_up"),
    IL1_signaling = g("IL1A_up"),
    Naive_Tcell = g("Naive_CD8"),
    Activation_Effector = c("B2M","CD40LG","CD81","GZMA","GZMB","GNLY","PRF1","IFNG",
                            "TNF","LTA","KLRG1","CX3CR1","FCGR3A","BAX","BAK1","FAS",
                            "HAVCR2","LAG3","ENTPD1","TBX21","BHLHE40","RORC")
  )
  for (nm in names(gene_signatures)) add_set(nm, gene_signatures[[nm]], "Signature")
  message(" Signature sets: ", paste(names(gene_signatures), collapse = ", "))
} else {
  message(" ", SIG_DIR, " not found - skipping the signature sets.")
}

fetch_msigdb <- function(category, subcategory = NULL, gs_name) {
  args <- list(species = "Homo sapiens", category = category)
  if (!is.null(subcategory)) args$subcategory <- subcategory
  tbl <- tryCatch(do.call(msigdbr, args), error = function(e) NULL)
  if (is.null(tbl)) return(character(0))
  hits <- unique(tbl$gene_symbol[tbl$gs_name == gs_name])
  if (length(hits) == 0) {
    partial <- unique(tbl$gs_name[grepl(gs_name, tbl$gs_name, fixed = TRUE)])
    if (length(partial) > 0) hits <- unique(tbl$gene_symbol[tbl$gs_name == partial[1]])
  }
  as.character(hits)
}

add_set("HALLMARK_IL2_STAT5_SIGNALING",
        fetch_msigdb("H", NULL, "HALLMARK_IL2_STAT5_SIGNALING"),
        "MSigDB", "HALLMARK: IL2-STAT5 Signaling")
add_set("KAECH_NAIVE_VS_DAY8_EFF_CD8_UP",
        fetch_msigdb("C7", "IMMUNESIGDB", "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP"),
        "MSigDB", "Kaech: Naive vs D8 Eff CD8 (UP)")
add_set("WP_TGFBETA_RECEPTOR_SIGNALING",
        fetch_msigdb("C2", "CP:WIKIPATHWAYS", "WP_TGFBETA_RECEPTOR_SIGNALING"),
        "MSigDB", "WP: TGF-beta Receptor Signaling")
add_set("HALLMARK_TGF_BETA_SIGNALING",
        fetch_msigdb("H", NULL, "HALLMARK_TGF_BETA_SIGNALING"),
        "MSigDB", "HALLMARK: TGF-beta Signaling")
add_set("HALLMARK_APOPTOSIS",
        fetch_msigdb("H", NULL, "HALLMARK_APOPTOSIS"),
        "MSigDB", "HALLMARK: Apoptosis")
add_set("HALLMARK_INFLAMMATORY_RESPONSE",
        fetch_msigdb("H", NULL, "HALLMARK_INFLAMMATORY_RESPONSE"),
        "MSigDB", "HALLMARK: Inflammatory Response")
add_set("HALLMARK_INTERFERON_GAMMA_RESPONSE",
        fetch_msigdb("H", NULL, "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
        "MSigDB", "HALLMARK: IFN-gamma Response")

if (file.exists(CD8_SIG_RDS)) {
  cd8_df <- readRDS(CD8_SIG_RDS)
  add_set("IL1a_up", as.character(cd8_df$GENE[cd8_df$TERM == "IL1a_up"]),
          "MSigDB", "IL1-alpha (cytokine up)")
} else {
  message(" ", CD8_SIG_RDS, " not found - skipping IL1a_up.")
}

# curated sets, also the sets used for module scoring l8r
curated <- list(
  Cat_TRM_program = c("CD69","CXCR4","NR4A2","RUNX3","ZEB2","TGFBR1","TGFBR2", "TGFBR3","TGFB1","NFKBIA"),
  Cat_Memory = c("IL7R","CCR7","LEF1","SELL","FOXO1","BCL2","NCL","FYN", "ETS1","NELL2","TCF7"),
  Cat_Effector_cytotoxic = c("GZMA","GZMB","PRF1","GNLY","NKG7","CCL5","TNF","GZMK", "CX3CR1","KLRG1","FCGR3A","IFNG"),
  Cat_Exhaustion_inhib = c("HAVCR2","LAG3","ENTPD1","PDCD1","TIGIT","TOX","CTLA4"),
  Cat_Activation_markers = c("TNFRSF4","TNFRSF18","CD40LG","ICOS","CD27","TNFRSF9", "SLAMF7","CD81","B2M"),
  Cat_IL1_TGFb_signalling = c("SATB1","BATF","STAT3","NFKBIA","EP300","IKZF1","SMAD2","SMAD3"),
  Cat_IL2_family = c("IL2RA","IL7R","LARP1","CISH","PIM1","NCL","SOCS3"),
  Cat_Terminal_diff_TFs = c("TBX21","BHLHE40","RORC","MAF","ID2","PRDM1"),
  Cat_Memory_stemness_TFs = c("MYB","STAT1","LEF1","TCF7","ID3","FOXO1","BACH2"),
  Cat_Apoptosis = c("BAX","BAK1","BCL2","FAS","BCL2L11","LGALS3","CASP3","CASP8")
)
for (nm in names(curated))
  add_set(nm, curated[[nm]], "Category", gsub("^Cat_", "", gsub("_", " ", nm)))

message(" Total sets: ", length(pathway_sets),
        " (Signature=", sum(pathway_source == "Signature"),
        ", MSigDB=", sum(pathway_source == "MSigDB"),
        ", Category=", sum(pathway_source == "Category"), ")")

# --- Object ------------------------------------------------------------------
message("\n-- Loading object --")
seu <- readRDS(EXT_RDS)
seu <- car_subset_obj(seu)
message("Loaded: ", ncol(seu), " cells | ", nrow(seu), " genes")

seu@meta.data <- add_response_group(seu@meta.data)

# --- Comparisons -------------------------------------------------------------
post_cells <- function(lineage, product = NULL) {
  md <- seu@meta.data
  keep <- !is.na(md[[T_LIN_COL]]) & md[[T_LIN_COL]] == lineage &
          md[[TIMEPOINT_COL]] %in% LATE_VALUES
  if (!is.null(product))
    keep <- keep & !is.na(md[[PRODUCT_COL]]) & md[[PRODUCT_COL]] == product
  rownames(md)[keep]
}

comparisons <- list()
for (lin in LINEAGES) {
  comparisons[[paste0(lin, "_cilta_vs_ide")]] <- list(
    cells = post_cells(lin), comp_col = PRODUCT_COL, grp1 = CILTA_VALUE, grp2 = IDE_VALUE,
    label = paste0(lin, " post-inf - cilta vs ide"), lineage = lin)
  comparisons[[paste0(lin, "_cilta_R_vs_NR")]] <- list(
    cells = post_cells(lin, CILTA_VALUE), comp_col = "response_group",
    grp1 = "Responder", grp2 = "Non-Responder",
    label = paste0(lin, " cilta - R vs NR"), lineage = lin)
  comparisons[[paste0(lin, "_ide_R_vs_NR")]] <- list(
    cells = post_cells(lin, IDE_VALUE), comp_col = "response_group",
    grp1 = "Responder", grp2 = "Non-Responder",
    label = paste0(lin, " ide - R vs NR"), lineage = lin)
}

# --- Ranked list -------------------------------------------------------------
make_ranked_list <- function(cells, comp) {
  meta <- seu@meta.data[cells, ]
  meta <- meta[!is.na(meta[[comp$comp_col]]) & meta[[comp$comp_col]] %in% c(comp$grp1, comp$grp2), ]
  meta <- keep_min_carpos(meta, patient_col = PATIENT_COL)
  if (nrow(meta) < 4) return(NULL)

  expr <- GetAssayData(seu[, rownames(meta)], assay = "RNA", slot = "data")
  patients <- unique(meta[[PATIENT_COL]])
  pb <- do.call(rbind, lapply(patients, function(pt) {
    bc  <- rownames(meta)[meta[[PATIENT_COL]] == pt]
    grp <- unique(meta[[comp$comp_col]][meta[[PATIENT_COL]] == pt])[1]
    data.frame(patient = pt, group = grp,
               t(as.matrix(rowMeans(expr[, bc, drop = FALSE]))))
  }))
  v1 <- pb[pb$group == comp$grp1, -(1:2), drop = FALSE]
  v2 <- pb[pb$group == comp$grp2, -(1:2), drop = FALSE]
  message("  n patients - ", comp$grp1, ": ", nrow(v1), " | ", comp$grp2, ": ", nrow(v2))
  if (nrow(v1) < 2 || nrow(v2) < 2) return(NULL)

  stats <- sapply(seq_len(ncol(v1)), function(i) {
    a <- as.numeric(v1[, i]); b <- as.numeric(v2[, i])
    lfc <- mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)
    p <- tryCatch(wilcox.test(a, b)$p.value, error = function(e) 1)
    sign(lfc) * (-log10(p + 1e-300))
  })
  names(stats) <- colnames(v1)
  sort(stats[!is.na(stats)], decreasing = TRUE)
}

# --- Run ---------------------------------------------------------------------
all_res <- list()
for (nm in names(comparisons)) {
  comp <- comparisons[[nm]]
  message("\n-- GSEA: ", comp$label, " --")
  ranks <- make_ranked_list(comp$cells, comp)
  if (is.null(ranks) || length(ranks) < 100) {
    message("  skipped: too few genes or patients")
    next
  }

  pw <- lapply(pathway_sets, function(gs) intersect(gs, names(ranks)))
  pw <- Filter(function(gs) length(gs) >= GSEA_MIN_SIZE, pw)
  if (!length(pw)) {
    message("  no sets survived the gene filter")
    next
  }

  set.seed(42)
  res <- as.data.frame(fgseaMultilevel(pathways = pw, stats = ranks,
                                       minSize = GSEA_MIN_SIZE, maxSize = Inf, eps = 0))
  res$comparison <- comp$label
  res$lineage <- comp$lineage
  res$source <- pathway_source[res$pathway]
  res$display <- ifelse(res$pathway %in% names(display_names),
                           display_names[res$pathway], res$pathway)
  res <- res[order(res$NES, decreasing = TRUE), ]
  message("  ", nrow(res), " sets tested | FDR<0.05: ", sum(res$padj < 0.05, na.rm = TRUE))
  all_res[[nm]] <- res

  res$sig <- ifelse(res$padj < 0.05, "FDR < 0.05", "ns")
  res$display <- factor(res$display, levels = res$display[order(res$NES)])
  p_bar <- ggplot(res, aes(x = NES, y = display, fill = sig)) +
    geom_col(width = 0.68, alpha = 0.9) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    facet_grid(source ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("FDR < 0.05" = "#C0392B", "ns" = "grey75")) +
    labs(title = paste0("GSEA - ", comp$label),
         subtitle = paste0("pseudobulk-ranked fgsea (>= ", CAR_MIN, " CAR+ cells/patient) | red = FDR<0.05"),
         x = "NES", y = NULL, fill = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          strip.text.y = element_text(angle = 0, face = "bold", size = 9),
          strip.background = element_rect(fill = "grey92", color = NA))
  save_pdf(p_bar, file.path(FIG_DIR, paste0("gsea_", nm, "_bar.pdf")), w = 11, h = 9)

  # Enrichment curves for everything that clears FDR
  sig_paths <- res$pathway[which(res$padj < 0.05)]
  if (length(sig_paths) > 0) {
    pdf(file.path(FIG_DIR, paste0("gsea_", nm, "_enrichment.pdf")), width = 8, height = 5)
    for (pwn in sig_paths) {
      dn <- if (pwn %in% names(display_names)) display_names[[pwn]] else pwn
      tryCatch(print(plotEnrichment(pw[[pwn]], ranks) +
        labs(title = dn, subtitle = paste0(comp$label,
             " | NES=", round(res$NES[res$pathway == pwn], 2),
             " padj=", signif(res$padj[res$pathway == pwn], 2)))),
        error = function(e) message("  plotEnrichment skipped: ", pwn))
    }
    dev.off()
    message("  Enrichment curves: ", paste(sig_paths, collapse = ", "))
  }
}

if (length(all_res) > 0) {
  combined <- do.call(rbind, lapply(all_res, function(r)
    r[, c("comparison", "lineage", "source", "pathway", "display",
          "NES", "pval", "padj", "size")]))
  write.csv(combined, file.path(FIG_DIR, "gsea_all_results.csv"), row.names = FALSE)
  message("\nSaved: gsea_all_results.csv (", nrow(combined), " rows)")
}

# CD8 and CD4 cilta-vs-ide side by side, ordered on the CD8 NES
r8 <- all_res[["CD8_cilta_vs_ide"]]
r4 <- all_res[["CD4_cilta_vs_ide"]]
if (!is.null(r8) && !is.null(r4)) {
  rc <- rbind(r8, r4)
  rc$sig <- ifelse(rc$padj < 0.05, "FDR < 0.05", "ns")
  rc$display <- factor(rc$display, levels = r8$display[order(r8$NES)])
  rc$lineage <- factor(rc$lineage, levels = c("CD8", "CD4"))
  p_comb <- ggplot(rc[!is.na(rc$display), ], aes(x = NES, y = display, fill = sig)) +
    geom_col(width = 0.68, alpha = 0.9) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    facet_grid(source ~ lineage, scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("FDR < 0.05" = "#C0392B", "ns" = "grey75")) +
    labs(title = "GSEA - post-infusion cilta vs ide (CD8 and CD4)",
         subtitle = "pseudobulk-ranked fgsea | red = FDR<0.05",
         x = "NES", y = NULL, fill = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey92", color = NA))
  save_pdf(p_comb, file.path(FIG_DIR, "gsea_cilta_vs_ide_CD8_CD4.pdf"), w = 12, h = 10)
}

# -----------------------------------------------------------------------------
# Module scores for the curated sets  (Fig 4f)
# Methods: AddModuleScore, mean per patient, two-sided Wilcoxon
# -----------------------------------------------------------------------------
message("\n-- Module scoring, curated sets --")

MODULE_SETS <- names(pathway_source)[pathway_source == "Category"]
mod_label <- function(x) gsub("^Cat_", "", gsub("_", " ", x))

score_modules_lineage <- function(lineage) {
  md <- seu@meta.data[post_cells(lineage), ]
  md <- md[!is.na(md[[PRODUCT_COL]]) & md[[PRODUCT_COL]] %in% c(CILTA_VALUE, IDE_VALUE), ]
  md <- keep_min_carpos(md, patient_col = PATIENT_COL)
  if (nrow(md) < 4) {
    message("  ", lineage, ": too few cells or patients, skipped")
    return(NULL)
  }

  sub <- seu[, rownames(md)]
  DefaultAssay(sub) <- "RNA"
  feats <- setNames(lapply(MODULE_SETS, function(s) intersect(pathway_sets[[s]], rownames(sub))),
                    MODULE_SETS)
  for (s in MODULE_SETS)
    message(" ", lineage, " | ", s, ": ", length(feats[[s]]), "/",
            length(pathway_sets[[s]]), " genes present")

  set.seed(42)
  sub <- tryCatch(
    AddModuleScore(sub, features = feats, name = "MODscore_", assay = "RNA", search = FALSE),
    error = function(e) {
      message("... retrying AddModuleScore with nbin = 15 (", conditionMessage(e), ")")
      AddModuleScore(sub, features = feats, name = "MODscore_", assay = "RNA",
                     search = FALSE, nbin = 15)
    })
  score_cols <- paste0("MODscore_", seq_along(MODULE_SETS))

  sm <- sub@meta.data
  pts <- unique(sm[[PATIENT_COL]])
  do.call(rbind, lapply(pts, function(pt) {
    r <- sm[sm[[PATIENT_COL]] == pt, ]
    pm <- colMeans(r[, score_cols, drop = FALSE])
    row <- data.frame(patient = pt, product = unique(r[[PRODUCT_COL]])[1], lineage = lineage, stringsAsFactors = FALSE)
    for (i in seq_along(MODULE_SETS)) row[[MODULE_SETS[i]]] <- unname(pm[i])
    row
  }))
}

pb_mod <- do.call(rbind, lapply(LINEAGES, score_modules_lineage))

if (!is.null(pb_mod) && nrow(pb_mod) > 0) {

  mod_stats <- do.call(rbind, lapply(LINEAGES, function(l) {
    d <- pb_mod[pb_mod$lineage == l, ]
    if (!nrow(d)) return(NULL)
    do.call(rbind, lapply(MODULE_SETS, function(m) {
      a <- d[[m]][d$product == CILTA_VALUE]
      b <- d[[m]][d$product == IDE_VALUE]
      wt <- tryCatch(wilcox.test(a, b), error = function(e) NULL)
      data.frame(
        lineage = l,
        module = mod_label(m),
        n_cilta = length(a),
        n_ide = length(b),
        median_cilta = median(a, na.rm = TRUE),
        median_ide = median(b, na.rm = TRUE),
        direction = ifelse(median(a, na.rm = TRUE) >= median(b, na.rm = TRUE), "cilta>ide", "ide>cilta"),
        W = if (!is.null(wt)) unname(wt$statistic) else NA_real_,
        p = if (!is.null(wt)) wt$p.value else NA_real_,
        stringsAsFactors = FALSE)
    }))
  }))
  mod_stats$p_BH_within_lineage <- ave(mod_stats$p, mod_stats$lineage,
                                       FUN = function(p) p.adjust(p, method = "BH"))
  write.csv(mod_stats, file.path(FIG_DIR, "module_scores_stats.csv"), row.names = FALSE)
  message("\nSaved: module_scores_stats.csv")
  print(mod_stats[, c("lineage", "module", "n_cilta", "n_ide", "direction", "p", "p_BH_within_lineage")])
    long <- do.call(rbind, lapply(MODULE_SETS, function(m)
    data.frame(patient = pb_mod$patient, product = pb_mod$product,
               lineage = pb_mod$lineage, module = mod_label(m),
               score = pb_mod[[m]], stringsAsFactors = FALSE)))
  long$product <- factor(long$product, levels = c(IDE_VALUE, CILTA_VALUE),
                         labels = c("Ide", "Cilta"))
  long$lineage <- factor(long$lineage, levels = LINEAGES)
  long$module <- factor(long$module,  levels = mod_label(MODULE_SETS))

  lab <- mod_stats
  lab$lineage <- factor(lab$lineage, levels = LINEAGES)
  lab$module <- factor(lab$module,  levels = mod_label(MODULE_SETS))
  lab$label <- paste0("p=", signif(lab$p, 2))

  p_mod <- ggplot(long, aes(product, score, fill = product)) +
    geom_violin(trim = FALSE, alpha = 0.85, scale = "width") +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.08, height = 0, size = 0.7, alpha = 0.5) +
    geom_text(data = lab, aes(x = 1.5, y = Inf, label = label),
              vjust = 1.6, size = 3, inherit.aes = FALSE) +
    facet_grid(lineage ~ module, scales = "free_y") +
    scale_fill_manual(values = c(Ide = unname(product_colors["ide"]),
                                 Cilta = unname(product_colors["cilta"]))) +
    labs(title = "Module scores: post-infusion CAR+ T cells",
         subtitle = "Per-patient means | two-sided Wilcoxon, cilta vs ide",
         x = NULL, y = "Module score", fill = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = 8),
          strip.background = element_rect(fill = "grey92", color = NA))
  save_pdf(p_mod, file.path(FIG_DIR, "module_scores_all.pdf"), w = max(10, 2.1 * length(MODULE_SETS)), h = 7)

  # two panels for the manuscript
  SELECTED_PANELS <- list(c("CD8", "Cat_Effector_cytotoxic"),
                          c("CD8", "Cat_Memory"))

  make_mod_panel <- function(lin, set) {
    d <- pb_mod[pb_mod$lineage == lin, c("patient", "product", set)]
    names(d)[3] <- "score"
    d$product <- factor(d$product, levels = c(IDE_VALUE, CILTA_VALUE),
                        labels = c("Ide", "Cilta"))
    st <- mod_stats[mod_stats$lineage == lin & mod_stats$module == mod_label(set), ]
    pv <- if (nrow(st) > 0) st$p[1] else NA_real_
    ggplot(d, aes(product, score, fill = product)) +
      geom_violin(trim = FALSE, alpha = 0.85, scale = "width") +
      geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.85) +
      geom_jitter(width = 0.08, height = 0, size = 0.8, alpha = 0.6) +
      annotate("text", x = 1.5, y = Inf, vjust = 1.6, size = 3.2, label = paste0("p = ", signif(pv, 2))) +
      scale_fill_manual(values = c(Ide   = unname(product_colors["ide"]),  Cilta = unname(product_colors["cilta"]))) +
      labs(title = paste0(lin, " - ", mod_label(set)), x = NULL, y = "Module score") +
      theme_classic(base_size = 11) +
      theme(legend.position = "none", plot.title = element_text(face = "bold", size = 10))
  }

  sel_panels <- lapply(SELECTED_PANELS, function(x) make_mod_panel(x[[1]], x[[2]]))
  p_sel <- wrap_plots(sel_panels, nrow = 1) +
    plot_annotation(title    = "Selected module scores",
                    subtitle = "Per-patient means | two-sided Wilcoxon, cilta vs ide")
  save_pdf(p_sel, file.path(FIG_DIR, "module_scores_selected.pdf"),
           w = 4 * length(sel_panels), h = 4)
} else {
  message("No module-score results produced.")
}

message("\n-- 02_gsea_and_module_scores.R done --")
