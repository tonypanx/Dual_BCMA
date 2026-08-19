library(Seurat)
library(plyr)
library(dplyr)

setwd("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/baseline_analysis")
source("/project/jhuangime/tony/Scripts/scRNA_core.R")

outdir <- "/project/jhuangime/tony/projects/BCMA_CAR/Pre_Infusion/Batch"

samples <- list.dirs(outdir, full.names = F, recursive = F)

mtxdirs <- paste0(outdir, '/', samples, '/outs/per_sample_outs/', samples, '/count/sample_filtered_feature_bc_matrix')

mtxdirs <- as.list(c(mtxdirs))
mtx_list <- lapply(mtxdirs, function(mm) return(Read10X(mm)))

seurat_list_gex <- lapply(mtx_list, function(mm) CreateSeuratObject(mm$'Gene Expression'))
seurat_list_adt <- lapply(mtx_list, function(mm) CreateAssayObject(mm$'Antibody Capture'))

for (i in 1:length(seurat_list_gex)){
  seurat_list_gex[[i]][["ADT"]] <- seurat_list_adt[[i]]
  seurat_list_gex[[i]]$sample <- samples[i]
}

obj_all <- merge(seurat_list_gex[[1]], seurat_list_gex[c(2:4)], merge.data = F, add.cell.ids = samples)
obj_all[['RNA']] <- JoinLayers(obj_all[['RNA']])
obj_all <- NormalizeData(obj_all)
obj_all[["percent.mt"]] <- PercentageFeatureSet(obj_all, pattern = "^MT-")

saveRDS(obj_all, "RDS/All_merged_noQC.rds")

VlnPlot(obj_all, features = c('percent.mt', 'nCount_RNA', 'nFeature_RNA'), pt.size = 0, group.by = 'sample')

obj_all <- subset(obj_all, percent.mt < 5)
obj_all <- subset(obj_all, nCount_RNA < 20000)
obj_all <- subset(obj_all, nFeature_RNA > 350)

saveRDS(obj_all, "RDS/All_merged_QC.rds")

obj_all <- Preprocess(obj_all, cellCycle = T, RunHarmony = F, nVgenes = 10000)

obj_all <- RunHarmony(obj_all, group.by.vars = "sample")
obj_all <- RunUMAP(obj_all, reduction = "harmony", dims = 1:50,
              min.dist = .3, n.neighbors = 30)
obj_all <- FindNeighbors(obj_all, dims = 1:50, reduction = 'harmony')
obj_all <- FindClusters(obj_all, resolution = .1)

saveRDS(obj_all, "RDS/All_integrated_QC.rds")
