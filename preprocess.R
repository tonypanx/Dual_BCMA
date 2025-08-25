library(Seurat)
library(plyr)
library(harmony)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(patchwork)
library(sctransform)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(SCpubr)
library(grid)
library(ggplotify)
library(cowplot)
library(ComplexHeatmap)
library(stats)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(RColorBrewer)
library(igraph)
require(tidygraph)
require(ggraph)
library(reshape2)
library(rlist)
library(scater)
library(ggrepel)
library(stringr)
library(scales)

Preprocess <- function(sobj, nVgenes=10000, minDistUmap=0.3, nNeighborUmap=30, cellCycle=F, harmVar, Human=T, RunHarmony = T, iter = 10, excludevars = NULL, 
                       isv5 = F, v5assay = NULL){
  
  
  sobj <- NormalizeData(sobj)
  sobj <- FindVariableFeatures(sobj, nfeatures = nVgenes, )
  
  
  if(isv5){
    DefaultLayer(sobj[[v5assay]]) <- 'data'
  }
  
  
  VariableFeatures(sobj) <- setdiff(VariableFeatures(sobj), excludevars)
  if (cellCycle){
    if (Human){
      s.genes <- cc.genes$s.genes
      g2m.genes <- cc.genes$g2m.genes
    } else {
      s.genes <- readRDS("/project2/jhuangime/tony/software/Mart/mouse_s_genes.rds")
      g2m.genes <- readRDS("/project2/jhuangime/tony/software/Mart/mouse_g2m_genes.rds")
    }
    sobj <- CellCycleScoring(sobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
    sobj$CC.Difference <- sobj$S.Score - sobj$G2M.Score
    sobj <- ScaleData(sobj, vars.to.regress = c("percent.mt", "CC.Difference", "nCount_RNA", "nFeature_RNA"))
  } else {
    sobj <- ScaleData(sobj, vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA"))
  }
  
  sobj <- RunPCA(sobj, npcs = 100, verbose = FALSE)
  if (RunHarmony){
    require(harmony)
    sobj <- RunHarmony(sobj, group.by.vars = harmVar, max_iter = iter)
    sobj <- RunUMAP(sobj, reduction = "harmony", dims = 1:50,
                    min.dist = minDistUmap, n.neighbors = nNeighborUmap)
  }else{
    sobj <- RunUMAP(sobj, reduction = "pca", dims = 1:50,
                    min.dist = minDistUmap, n.neighbors = nNeighborUmap)
  }
  return(sobj)
}
