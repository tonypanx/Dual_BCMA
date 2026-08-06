# -----------------------------------------------------------------------------
# Per-cell module scores for the signature panel, on CAR+ cells.
#
# Scores every signature with AddModuleScore and saves a slim metadata table
# (scores, a subset of gene expression, UMAP coordinates and the clinical columns) 
# so 03_regulon_pseudobulk.R can run off the scores without reloading the 1.8 GB object
#
# Input:  Rade et al. integrated T-cell object, Signatures
# Output: interim/rade_meta_scored_carpos.rds
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

source("config.R")
source(file.path("R", "common.R"))
source(file.path("R", "car_gating.R"))

N_GENES_TO_SLIM <- 20

# --- Signatures --------------------------------------------------------------
message("\n-- Loading gene signatures --")

if (!dir.exists(SIG_DIR)) stop("Signatures directory not found: ", SIG_DIR)

rds_files <- list.files(SIG_DIR, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
if (!length(rds_files)) stop("No .rds files in ", SIG_DIR)

raw_sigs <- setNames(lapply(rds_files, readRDS), tools::file_path_sans_ext(basename(rds_files)))
message("Signature files: ", paste(names(raw_sigs), collapse = ", "))

extract_genes <- function(obj, key = "") {
  if (is.null(obj))
    stop("Signature '", key, "' not found. Available: ", paste(names(raw_sigs), collapse = ", "))
  if (is.character(obj)) return(obj)
  if (is.data.frame(obj) && "gene" %in% colnames(obj)) return(obj$gene)
  if (is.list(obj) && length(obj) > 0) return(obj[[1]])
  stop("Unrecognised format for '", key, "': ", paste(class(obj), collapse = "/"))
}

gene_signatures <- list(
  CAR_TRM = extract_genes(raw_sigs[["CAR_TRM_gene_set"]], "CAR_TRM_gene_set"),
  RUNX3_regulon = extract_genes(raw_sigs[["RUNX3_Regulon"]], "RUNX3_Regulon"),
  KLF2_regulon = extract_genes(raw_sigs[["KLF2_Regulon"]], "KLF2_Regulon"),
  TGFb_signaling = extract_genes(raw_sigs[["TGFB"]], "TGFB"),
  Hallmark_IL2 = extract_genes(raw_sigs[["Hallmark_IL2"]], "Hallmark_IL2"),
  IL2_up = extract_genes(raw_sigs[["IL2_up"]], "IL2_up"),
  IL1_signaling = extract_genes(raw_sigs[["IL1A_up"]], "IL1A_up"),
  Naive_Tcell = extract_genes(raw_sigs[["Naive_CD8"]], "Naive_CD8"),
  Activation_Effector = c("B2M", "CD40LG", "CD81", "GZMA", "GZMB", "GNLY", "PRF1", "IFNG", 
                          "TNF", "LTA", "KLRG1", "CX3CR1", "FCGR3A", "BAX", "BAK1",
                          "FAS", "HAVCR2", "LAG3", "ENTPD1", "TBX21", "BHLHE40", "RORC")
)

CAR_TRM_CAP <- 500
gene_signatures[["CAR_TRM"]] <- head(gene_signatures[["CAR_TRM"]], CAR_TRM_CAP)

message("Signature sizes (CAR_TRM capped at ", CAR_TRM_CAP, "):")
for (nm in names(gene_signatures))
  message("  ", nm, ": ", length(gene_signatures[[nm]]), " genes")

# --- Score -------------------------------------------------------------------
message("\n-- Loading object and scoring --")
seu <- readRDS(EXT_RDS)
message("Loaded: ", ncol(seu), " cells | ", nrow(seu), " genes")

seu <- car_subset_obj(seu)

for (sig_name in names(gene_signatures)) {
  genes_present <- intersect(gene_signatures[[sig_name]], rownames(seu))
  message("  ", sig_name, ": ", length(genes_present), "/",
          length(gene_signatures[[sig_name]]), " genes in the dataset")
  if (length(genes_present) < 3) {
    warning("Too few genes for ", sig_name, ", skipping")
    next
  }
  seu <- AddModuleScore(seu, features = list(genes_present), name = paste0(sig_name, "_score"),
                        ctrl = 100, seed = 42)
  seu@meta.data[[paste0(sig_name, "_score")]] <- seu@meta.data[[paste0(sig_name, "_score1")]]
  seu@meta.data[[paste0(sig_name, "_score1")]] <- NULL
}
message("Scoring complete.")

# --- Slim metadata table -----------------------------------------------------
message("\n-- Building slim metadata --")

paper_genes <- list(
  TRM_memory = c("CD69", "CXCR4", "NR4A2", "IL7R", "CCR7", "LEF1",
                   "SELL", "FOXO1", "TGFBR2", "TGFBR3"),
  Effector_exh = c("GZMA", "GZMB", "PRF1", "GNLY", "NKG7",
                   "HAVCR2", "LAG3", "ENTPD1", "PDCD1"),
  Activation = c("TNFRSF4", "CD27", "TNFRSF9", "SLAMF7",
                   "CD40LG", "ICOS", "B2M"),
  Survival_TF = c("BCL2", "BAX", "BAK1", "TBX21", "BHLHE40",
                   "RORC", "MYB", "STAT1", "TCF7", "BATF"),
  IL_signaling = c("SATB1", "NFKBIA", "IKZF1", "STAT3",
                   "IL2RA", "CISH", "PIM1", "SOCS3")
)

sig_top_genes <- unique(unlist(lapply(gene_signatures, function(g)
  head(intersect(g, rownames(seu)), N_GENES_TO_SLIM))))
genes_to_slim <- unique(c(sig_top_genes,
                          intersect(unique(unlist(paper_genes)), rownames(seu))))
message("Genes carried into the slim table: ", length(genes_to_slim))

gene_expr_df <- as.data.frame(t(as.matrix(
  GetAssayData(seu, layer = "data")[genes_to_slim, , drop = FALSE])))
colnames(gene_expr_df) <- paste0("expr_", colnames(gene_expr_df))

umap_coords <- as.data.frame(Embeddings(seu, reduction = "umap"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

score_cols <- intersect(paste0(names(gene_signatures), "_score"), colnames(seu@meta.data))

keep_meta_cols <- intersect(
  c(PATIENT_COL, PRODUCT_COL, TIMEPOINT_COL, RESPONSE_COL,
    SUBTYPE_COL, "celltype_short_2", "celltype_short_3",
    T_LIN_COL, "CRS_GRADE", "CRS_GROUP", "SEX", "AGE_AT_CAR",
    CAR_COL, score_cols),
  colnames(seu@meta.data))

meta <- cbind(seu@meta.data[, keep_meta_cols], umap_coords, gene_expr_df)
meta$cell_barcode <- rownames(meta)

message("Slim metadata: ", nrow(meta), " cells x ", ncol(meta), " columns")
message("Score columns: ", paste(score_cols, collapse = ", "))

saveRDS(meta, SCORED_META_RDS)
message("Saved: ", SCORED_META_RDS,
        " (", round(file.size(SCORED_META_RDS) / 1e6, 1), " MB)")
