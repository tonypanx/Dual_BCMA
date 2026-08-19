
library(Seurat)
library(plyr)
library(dplyr)
library(ggtext)
library(readxl)
library(fgsea)
library(clusterProfiler)
library(forcats)
library(Cairo)

options(max.print=9999)

c5 <- read.gmt("/project2/jhuangime/tony/software/gsea/c5.gmt")
c2 <- read.gmt("/project2/jhuangime/tony/software/gsea/c2.gmt")
c8 <- read.gmt("/project2/jhuangime/tony/software/gsea/c8.gmt")
c7 <- read.gmt("/project2/jhuangime/tony/software/gsea/c7.gmt")
cp <- read.gmt("/project2/jhuangime/tony/software/gsea/cp.gmt")
go <- read.gmt("/project2/jhuangime/tony/software/gsea/go.gmt")
hallmark <- read.gmt("/project2/jhuangime/tony/software/gsea/hallmark.gmt")

showtext::showtext_auto()

options(ggrepel.max.overlaps = Inf)

all.list <- rbind(c5, c2, c8, c7, cp, go, hallmark)
all.list$term <- as.character(all.list$term)

setwd("/project/jhuangime/tony/projects/BCMA_CAR/Pre_Infusion/Analysis")
source("/project/jhuangime/tony/Scripts/scRNA_core.R")

color_list <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting/colors.rds")

obj <- readRDS('RDS/obj_integrated_QC.rds')

Idents(obj) <- obj$Treatment_Time

cd4 <- subset(obj, cell_type == 'T CD4')
cd8 <- subset(obj, cell_type == 'T CD8')

t_obj <- subset(obj, cell_type %in% c('T CD4', 'T CD8'))
t_obj$ttype <- dplyr::recode(t_obj$cell_type, 'T CD4' = 'CD4', 'T CD8' = 'CD8')

cd4 <- Preprocess(cd4, cellCycle = T, RunHarmony = F, nVgenes = 6000)

cd4 <- RunHarmony(cd4, group.by.vars = "Sample")
cd4 <- RunUMAP(cd4, reduction = "harmony", dims = 1:50,
               min.dist = .3, n.neighbors = 30)
cd4 <- FindNeighbors(cd4, dims = 1:50, reduction = 'harmony')
cd4 <- FindClusters(cd4, resolution = .8)

cd4 <- subset(cd4, seurat_clusters != '5')

cd4 <- Preprocess(cd4, cellCycle = T, RunHarmony = F, nVgenes = 6000)

cd4 <- RunHarmony(cd4, group.by.vars = "Sample")
cd4 <- RunUMAP(cd4, reduction = "harmony", dims = 1:50,
               min.dist = .3, n.neighbors = 30)
cd4 <- FindNeighbors(cd4, dims = 1:50, reduction = 'harmony')
cd4 <- FindClusters(cd4, resolution = .5)

markers.cd4 <- clean_markers(FindAllMarkers(cd4, min.pct = .2))

###3 and 4 are myeloid and cd8

cd4 <- subset(cd4, seurat_clusters %!in% c('3', '4'))

cd4 <- Preprocess(cd4, cellCycle = T, RunHarmony = F, nVgenes = 6000)

cd4 <- RunHarmony(cd4, group.by.vars = "Sample")
cd4 <- RunUMAP(cd4, reduction = "harmony", dims = 1:50,
               min.dist = .3, n.neighbors = 30)
cd4 <- FindNeighbors(cd4, dims = 1:50, reduction = 'harmony')
cd4 <- FindClusters(cd4, resolution = .5)

markers.cd4 <- clean_markers(FindAllMarkers(cd4, min.pct = .2))

###2: Treg(Foxp3, IKZF2, Il2RA)
###1: CM (CCR7, TCF7, LEf1)
###0: EM (IL7R, GZMA, BHLHE40)

cd4$annofine <- dplyr::recode(cd4$seurat_clusters, '0' = 'EM', '1' = 'CM', '2' = 'Treg')

features <- c('CD4', 'HOPX', 'BHLHE40', 'CCR7', 'LEF1', "TCF7", 'FOXP3', 'IL2RA')

p <- DimPlot_better(cd4,  group.by = "annofine") + theme(plot.title = element_blank()) 

ggsave('Plots/new_plots/preinfusion_cd4_umap.pdf',  plot = p, width = 4200, height = 3200, dpi = 900, units = 'px')

p <- VlnPlot_better(cd4, genes = features, group.by = "annofine")

ggsave('Plots/new_plots/preinfusion_cd4_vln.pdf',  plot = p, width = 4200, height = 5400, dpi = 900, units = 'px')

meta <- cd4@meta.data %>%
  select(annofine, Treatment) %>%          
  mutate(Treatment = factor(Treatment, levels = c("Ide", "Cilta")))

p2 <- meta %>% 
  group_by(Treatment, annofine) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  group_by(Treatment) %>% 
  mutate(pct = n / sum(n)) %>% 
  ggplot(aes(
    x = fct_rev(fct_inorder(Treatment)),
    y = pct,
    fill = fct_rev(fct_inorder(annofine))
  )) +
  geom_bar(stat = "identity", color = 'black') +
  labs(
    x = NULL,
    y = "Proportion of Sample",
    fill = "CD4 Cluster"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12)
  ) +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE))

ggsave('Plots/new_plots/preinfusion_cd4_prop_bar.pdf',  plot = p2, width = 5200, height = 2400, dpi = 900, units = 'px')


cd8 <- Preprocess(cd8, cellCycle = T, RunHarmony = F, nVgenes = 6000)

cd8 <- RunHarmony(cd8, group.by.vars = "Sample")
cd8 <- RunUMAP(cd8, reduction = "harmony", dims = 1:50,
               min.dist = .3, n.neighbors = 30)
cd8 <- FindNeighbors(cd8, dims = 1:50, reduction = 'harmony')
cd8 <- FindClusters(cd8, resolution = .2)

markers.cd8 <- clean_markers(FindAllMarkers(cd8, min.pct = .2))

#0: TEm GZMK
#1: TE GNLY GZMB
#2: 
features <- c('GZMK', 'LEF1', 'IL7R', 'TCF7')

## 3 is contaminant

cd8 <- subset(cd8, seurat_clusters != '3')

cd8 <- Preprocess(cd8, cellCycle = T, RunHarmony = F, nVgenes = 6000)

cd8 <- RunHarmony(cd8, group.by.vars = "Sample")
cd8 <- RunUMAP(cd8, reduction = "harmony", dims = 1:50,
               min.dist = .3, n.neighbors = 30)
cd8 <- FindNeighbors(cd8, dims = 1:50, reduction = 'harmony')
cd8 <- FindClusters(cd8, resolution = .2)

markers.cd8 <- clean_markers(FindAllMarkers(cd8, min.pct = .2))


## remove 3 and 4, contaminant

cd8 <- subset(cd8, seurat_clusters %!in% c('3', '4'))

cd8 <- Preprocess(cd8, cellCycle = T, RunHarmony = F, nVgenes = 6000)

cd8 <- RunHarmony(cd8, group.by.vars = "Sample")
cd8 <- RunUMAP(cd8, reduction = "harmony", dims = 1:50,
               min.dist = .3, n.neighbors = 30)
cd8 <- FindNeighbors(cd8, dims = 1:50, reduction = 'harmony')
cd8 <- FindClusters(cd8, resolution = .3)

markers.cd8 <- clean_markers(FindAllMarkers(cd8, min.pct = .2))

features <- c('CD8A', 'GZMK','IL7R', 'TCF7', 'PRF1', 'GZMA', 'GZMB')

cd8$annofine <- dplyr::recode(cd8$seurat_clusters, '0' = 'EM', '1' = 'TE1', '2' = 'TE2')

p <- DimPlot_better(cd8,  group.by = "annofine") + theme(plot.title = element_blank()) 

ggsave('Plots/preinfusion_cd8_umap.pdf',  plot = p, width = 4200, height = 3200, dpi = 900, units = 'px')

p <- VlnPlot_better(cd8, genes = features, group.by = "annofine")

ggsave('Plots/new_plots/preinfusion_cd8_vln.pdf',  plot = p, width = 4200, height = 5800, dpi = 900, units = 'px')


meta <- cd8@meta.data %>%
  select(annofine, Treatment) %>%          
  mutate(Treatment = factor(Treatment, levels = c("Ide", "Cilta")))

p2 <- meta %>% 
  group_by(Treatment, annofine) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  group_by(Treatment) %>% 
  mutate(pct = n / sum(n)) %>% 
  ggplot(aes(
    x = fct_rev(fct_inorder(Treatment)),
    y = pct,
    fill = fct_rev(fct_inorder(annofine))
  )) +
  geom_bar(stat = "identity", color = 'black') +
  labs(
    x = NULL,
    y = "Proportion of Sample",
    fill = "CD8 Cluster"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12)
  ) +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE))

ggsave('Plots/preinfusion_cd8_prop_bar.pdf',  plot = p2, width = 5200, height = 2400, dpi = 900, units = 'px')

### testing DEGs for CAR-TRM sig ###

Idents(cd4) <- cd4$Treatment_Time

Idents(cd8) <- cd8$Treatment_Time

cd4.pre.contrast <- FindMarkers(cd4, ident.1 = 'Cilta_Pre', ident.2 = 'Ide_Pre', min.pct = .2, logfc.threshold = 0)
cd8.pre.contrast <- FindMarkers(cd8, ident.1 = 'Cilta_Pre', ident.2 = 'Ide_Pre', min.pct = .2, logfc.threshold = 0)

cd8.pre.contrast$gene <- rownames(cd8.pre.contrast)
cd4.pre.contrast$gene <- rownames(cd4.pre.contrast)

cd4.pre.contrast_filt <- clean_markers(cd4.pre.contrast)
cd8.pre.contrast_filt <- clean_markers(cd8.pre.contrast)

cd8.mem <- c('NELL2', 'LEF1', 'SELL', 'IL7R', 'TCF7', 'BACH2', 'CCR7', 'FOXO1')
cd8.exh <- c('HAVCR2', 'LAG3', 'ENTPD1', 'PDCD1', 'CTLA4',  'TIGIT', 'TOX')

cd8.pre.contrast$p_val_adj <- ifelse(cd8.pre.contrast$p_val_adj < 10^-10, 10^-10, cd8.pre.contrast$p_val_adj)
cd8.pre.contrast$avg_log2FC <- ifelse(cd8.pre.contrast$avg_log2FC > 2.5, 2.5, cd8.pre.contrast$avg_log2FC)
cd8.pre.contrast$avg_log2FC <- ifelse(cd8.pre.contrast$avg_log2FC < -2.5, -2.5, cd8.pre.contrast$avg_log2FC)

genes_validate <- c('RUNX3', 'RUNX2', 'TGFBR3', 'TGFBR2', 'STAT1', 'NCL', 'BCL2', 'FYN', 'LARP1', 'SATB1', 'NFKBIA', 'BATF', 'STAT3', 
                    'B2M', 'CD40LG', 'CD81', 'RAC2', 'GZMA', 'GZMB', 'GNLY', 'IFNG', 'CXCR3', 'LGALS3', 'BAX', 'BAK1', "FAS")


p1 <- volcanoplot_better(cd8.pre.contrast, genes.highlight = genes_validate) 
p2 <- VlnPlot(cd8, features = setdiff(genes_validate, rownames(cd8.pre.contrast)))


ggsave('Plots/cd8_preinfusion_deg_validate.pdf', plot = p1, height = 7600, width = 3600*3, dpi = 1200, units = 'px')
ggsave('Plots/cd8_preinfusion_deg_validate_missing.pdf', plot = p2, width = 3600, height = 3600, dpi = 300, units = 'px')

p3 <- volcanoplot_better(cd4.pre.contrast, genes.highlight = genes_validate) 
p4 <- VlnPlot(cd4, features = setdiff(genes_validate, rownames(cd4.pre.contrast)))

ggsave('Plots/cd4_preinfusion_deg_validate.pdf', plot = p3, height = 9600, width = 4800*3, dpi = 1200, units = 'px')
ggsave('Plots/cd4_preinfusion_deg_validate_missing.pdf', plot = p4, width = 3600, height = 3600, dpi = 300, units = 'px')

p5 <- volcanoplot_better(cd8.pre.contrast, genes.highlight = c(cd8.mem, cd8.exh)) 
p6 <- volcanoplot_better(cd4.pre.contrast, genes.highlight = c(cd8.mem, cd8.exh)) 
p7 <- VlnPlot(cd4, features = setdiff(c(cd8.mem, cd8.exh), rownames(cd4.pre.contrast)))

ggsave('Plots/cd8_preinfusion_deg_validate_mem_exh.pdf', plot = p5, height = 7600, width = 3600*3, dpi = 1200, units = 'px')
ggsave('Plots/cd4_preinfusion_deg_validate_mem_exh.pdf', plot = p6, height = 7600, width = 3600*3, dpi = 1200, units = 'px')
ggsave('Plots/cd4_preinfusion_deg_validate_mem_exh_missing.pdf', plot = p7, width = 3600, height = 3600, dpi = 300, units = 'px')

cd8.pre.contrast.ranks <- cd8.pre.contrast_filt$avg_log2FC
names(cd8.pre.contrast.ranks) <- rownames(cd8.pre.contrast_filt)

trm_car_data <- read_excel('/project/jhuangime/tony/projects/BCMA_CAR/external_data/trm_car_1.xlsx', col_names = T, skip = 0)
trm_car_data <- trm_car_data[trm_car_data$`P Value` < .05, ]

car_trm_sigs_filt2 <- list(trm_car_up =  trm_car_data[trm_car_data$`log FC` > .33, ]$`Gene ID`)

### cd8 path enrichment ##

cd8.paths <- gsea_simple(cd8.pre.contrast_filt, list(c5, c2, c7, hallmark), minGSSize = 1, pvalueCutoff = 1)
dplyr::filter(cd8.paths, Description %in% c("GSE9650_EFFECTOR_VS_MEMORY_CD8_TCELL_UP",  "GO_T_CELL_ACTIVATION", 
                                            "GO_CELL_KILLING",  "GO_T_CELL_APOPTOTIC_PROCESS" ,  "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP", "WP_TGFBETA_RECEPTOR_SIGNALING" ))

cd8_cyto <- readRDS('/project/jhuangime/tony/software/immune_Dict/T_cell_CD8_comb.rds')

cd8.csig <- gsea_simple(cd8.pre.contrast_filt, list(c5, c2, c7, hallmark), minGSSize = 1, pvalueCutoff = 1)


cd4.paths <- gsea_simple(cd4.pre.contrast_filt, list(c5, c2, c7, hallmark), minGSSize = 1, pvalueCutoff = 1)


## cd8 enrichment ##

out <- GSEA(cd8.pre.contrast.ranks[abs(cd8.pre.contrast.ranks) > .33], 
            TERM2GENE = data.frame(term= 'trm', gene=car_trm_sigs_filt2$trm_car_up), 
            pvalueCutoff = 1) ###NES = 1.43, p.adj = 0.0919


p1 <- plotEnrichment(car_trm_sigs_filt2$trm_car_up, cd8.pre.contrast.ranks[abs(cd8.pre.contrast.ranks) > .33]) + ggtitle('CD8 Pre-Infusion CAR-TRM vs CAR-TConv UP') +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text = element_text(size = 12), 
        plot.title = element_text(hjust = .5, face = 'bold', size = 14)) + ylab('Running Enrichment Score') + xlab('Rank') +
  annotate('text', x= 30, y = .055, label = 'NES=1.43', hjust = 0)+
  annotate('text', x= 30, y = .040, label = 'Adjusted p-value=0.0919', hjust = 0)

ggsave('Plots/new_plots/preinfusion_cd8_trm_gsea.pdf', plot = p1, height = 4000, width = 4400, dpi = 850, units = 'px')


## cd4 enrichment ##
cd4.pre.contrast <- clean_markers(cd4.pre.contrast)

cd4.pre.contrast.ranks <- cd4.pre.contrast$avg_log2FC
names(cd4.pre.contrast.ranks) <- rownames(cd4.pre.contrast)


out <- GSEA(cd4.pre.contrast.ranks[abs(cd4.pre.contrast.ranks) > .33], 
            TERM2GENE = data.frame(term= 'trm', gene=car_trm_sigs_filt2$trm_car_up), 
            pvalueCutoff = 1, minGSSize = 1) ###NES = -0.57, p.adj = .942

p1 <- plotEnrichment(car_trm_sigs_filt2$trm_car_up, cd4.pre.contrast.ranks[abs(cd4.pre.contrast.ranks) > .33]) + ggtitle('CD4 Pre-Infusion CAR-TRM vs CAR-TConv UP') +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text = element_text(size = 12), 
        plot.title = element_text(hjust = .5, face = 'bold', size = 14)) + ylab('Running Enrichment Score') + xlab('Rank') +
  annotate('text', x= 0, y = -.065, label = 'NES=-0.57', hjust = 0)+
  annotate('text', x= 00, y = -.08, label = 'Adjusted p-value=0.942', hjust = 0)

ggsave('Plots/new_plots/preinfusion_cd4_trm_gsea.pdf', plot = p1, height = 4000, width = 4400, dpi = 850, units = 'px')



### contrast scores ###


cd8$Treatment <- factor(cd8$Treatment, levels = c('Ide-pre', 'Cilta-pre'))

cd8 <- AddModuleScore(cd8, features = list(car_trm_sigs_filt2$trm_car_up), name = 'car_trm_up_')
cd4 <- AddModuleScore(cd4, features = list(car_trm_sigs_filt2$trm_car_up), name = 'car_trm_up_')

p1 <- vsplot(cd8, groups = 'Treatment', gene = c('car_trm_up_1'), pt.size = 0, colors = c('skyblue', 'red'), 
             feature_type = 'score') + ggtitle('Preinfusion CD8 CAR TRM Score')  + NoLegend() + 
  theme(axis.text.x = element_text(size = 15))

p2 <- vsplot(cd4, groups = 'Treatment', gene = c('car_trm_up_1'), pt.size = 0, colors = c('skyblue', 'red'), 
             feature_type = 'score') + ggtitle('Preinfusion CD4 CAR TRM Score')  + NoLegend() + 
  theme(axis.text.x = element_text(size = 15))

pa <- ggarrange(plotlist = list(p1, p2), nrow = 1)
ggsave('Plots/new_plots/preinfusion_cartrm_treatmentcomp.pdf', plot = pa, height = 3200, width = 7200, dpi = 850, units = 'px')


p1 <- vsplot(cd4, groups = 'Treatment', gene = c('car_trm_up_1'), pt.size = 0, colors = color_list$treat, 
             feature_type = 'score') + ggtitle('Preinfusion CD4 CAR TRM Score')


saveRDS(cd8, 'RDS/CD8.rds')
saveRDS(cd4, 'RDS/CD4.rds')
### gene contrast ###


### add regulon scores ###

regs_df <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/data/cd8_ip_reg_df.rds")
reg_runx3 <- dplyr::filter(regs_df, regulon %in% paste0(c('RUNX3'), '(+)')) %>% group_by(regulon) %>% 
  top_n(100, wt=importance)  
reg_klf2 <- dplyr::filter(regs_df, regulon %in% paste0(c('KLF2'), '(+)')) %>% group_by(regulon) %>% 
  top_n(100, wt=importance)  

cd8$Treatment <- dplyr::recode(cd8$Treatment, 'Ide' = 'Ide-pre', 'Cilta' = 'Cilta-pre')
cd4$Treatment <- dplyr::recode(cd4$Treatment, 'Ide' = 'Ide-pre', 'Cilta' = 'Cilta-pre')

cd8 <- AddModuleScore(cd8, features = list(reg_runx3$gene), name = 'RUNX3_Top_100_')
cd8 <- AddModuleScore(cd8, features = list(reg_klf2$gene), name = 'KLF2_Top_100_')

cd8$Treatment <- factor(cd8$Treatment, levels = c('Ide-pre', 'Cilta-pre'))
cd4$Treatment <- factor(cd4$Treatment, levels = c('Ide-pre', 'Cilta-pre'))

p1 <- vsplot(cd8, feature_type = 'score', assay5 = T, 
             gene = 'RUNX3_Top_100_1', groups = 'Treatment', pt.size = 0,  colors = c('skyblue', 'red')) + ggtitle('RUNX3 Regulon Score') +
  theme(axis.text.x = element_text(size = 15)) + NoLegend()
p2 <- vsplot(cd8, feature_type = 'score', assay5 = T, 
             gene = 'KLF2_Top_100_1', groups = 'Treatment', pt.size = 0,  colors = c('skyblue', 'red'))+ ggtitle('KLF2 Regulon Score')+
  theme(axis.text.x = element_text(size = 15)) + NoLegend()

pc <- ggarrange(plotlist = list(p1, p2), nrow = 1)

ggsave('Plots/preinfusion_cd8_runx3_klf2_comp.pdf', plot = pc, height = 3200, width = 7200, dpi = 850, units = 'px')

cd4 <- AddModuleScore(cd4, features = list(reg_runx3$gene), name = 'RUNX3_Top_100_')
cd4 <- AddModuleScore(cd4, features = list(reg_klf2$gene), name = 'KLF2_Top_100_')

p1 <- vsplot(cd4, feature_type = 'score', assay5 = T, 
             gene = 'RUNX3_Top_100_1', groups = 'Treatment', pt.size = 0,  colors = c('skyblue', 'red')) + ggtitle('RUNX3 Regulon Score') +
  theme(axis.text.x = element_text(size = 15)) + NoLegend()
p2 <- vsplot(cd4, feature_type = 'score', assay5 = T, 
             gene = 'KLF2_Top_100_1', groups = 'Treatment', pt.size = 0,  colors = c('skyblue', 'red'))+ ggtitle('KLF2 Regulon Score')+
  theme(axis.text.x = element_text(size = 15)) + NoLegend()

pc <- ggarrange(plotlist = list(p1, p2), nrow = 1)

ggsave('Plots/new_plots/preinfusion_cd4_runx3_klf2_comp.pdf', plot = pc, height = 3200, width = 7200, dpi = 850, units = 'px')



### cytokine signals ###

T_cell_CD8_comb <- readRDS("/project/jhuangime/tony/software/immune_Dict/T_cell_CD8_comb.rds")
colnames(T_cell_CD8_comb) <- c('term', 'gene')

T_cell_CD4_comb <- readRDS("/project/jhuangime/tony/software/immune_Dict/T_cell_CD4_comb.rds")
colnames(T_cell_CD4_comb) <- c('term', 'gene')

cd8 <- AddModuleScore(cd8, features = list(dplyr::filter(T_cell_CD8_comb, term == 'IL1a_up')$gene), name = 'IL1a_')
cd8 <- AddModuleScore(cd8, features = list(dplyr::filter(T_cell_CD8_comb, term == 'IL2_up')$gene), name = 'IL2_v2_')

cd8.top_terms <- c("HALLMARK_IL2_STAT5_SIGNALING" , 
                   "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP", "WP_TGFBETA_RECEPTOR_SIGNALING" )

cd4 <- AddModuleScore(cd4, features = list(dplyr::filter(T_cell_CD4_comb, term == 'IL1a_up')$gene), name = 'IL1a_')

paths_names <- c( 'IL2_', 'NAIVE_', "TGFB_")

for(i in 1:length(cd8.top_terms)){
  print(unique(dplyr::filter(all.list, term == cd8.top_terms[i])$gene))
  cd8 <- AddModuleScore(cd8, features = list(unique(dplyr::filter(all.list, term == cd8.top_terms[i])$gene)), name = paths_names[i])
}


p3 <- vsplot(cd8, feature_type = 'score', assay5 = T, 
             gene = 'IL1a_1', groups = 'Treatment', pt.size = 0,  colors = c('skyblue', 'red')) + ggtitle('IL-1a Score')+
  theme(axis.text.x = element_text(size = 15)) + NoLegend()
p4 <- vsplot(cd8, feature_type = 'score', assay5 = T, 
             gene = 'NAIVE_1', groups = 'Treatment', pt.size = 0,  colors = c('skyblue', 'red'))+ ggtitle('Naive CD8 Score')+
  theme(axis.text.x = element_text(size = 15)) + NoLegend()
p5 <- vsplot(cd8, feature_type = 'score', assay5 = T, 
             gene = 'TGFB_1', groups = 'Treatment', pt.size = 0,  colors = c('skyblue', 'red'))+ ggtitle('TGF-beta Score')+
  theme(axis.text.x = element_text(size = 15)) + NoLegend()

pe <- ggarrange(plotlist = list(p1, p2, p3, p4, p5), nrow = 1)

ggsave('Plots/new_plots/pre_cd8_cytosigs.pdf', plot = pe, height = 3200, width = 3000*5, dpi = 850, units = 'px')


for(i in 1:length(cd8.top_terms)){
  print(unique(dplyr::filter(all.list, term == cd8.top_terms[i])$gene))
  cd4 <- AddModuleScore(cd4, features = list(unique(dplyr::filter(all.list, term == cd8.top_terms[i])$gene)), name = paths_names[i])
}


p3 <- vsplot(cd4, feature_type = 'score', assay5 = T, 
             gene = 'IL1a_1', groups = 'Treatment', pt.size = 0,  colors = c('skyblue', 'red')) + ggtitle('IL-1a Score') +
  theme(axis.text.x = element_text(size = 15)) + NoLegend()
p4 <- vsplot(cd4, feature_type = 'score', assay5 = T, 
             gene = 'NAIVE_1', groups = 'Treatment', pt.size = 0,  colors = c('skyblue', 'red'))+ ggtitle('Naive CD8 Score') +
  theme(axis.text.x = element_text(size = 15)) + NoLegend()
p5 <- vsplot(cd4, feature_type = 'score', assay5 = T, 
             gene = 'TGFB_1', groups = 'Treatment', pt.size = 0,  colors = c('skyblue', 'red'))+ ggtitle('TGF-beta Score') +
  theme(axis.text.x = element_text(size = 15)) + NoLegend()

pe <- ggarrange(plotlist = list(p1, p2, p3, p4, p5), nrow = 1)

ggsave('Plots/new_plots/pre_cd4_cytosigs.pdf', plot = pe, height = 3200, width = 3000*5, dpi = 850, units = 'px')



### volcano plot ###

cd8.pre.contrast <- FindMarkers(cd8, ident.1 = 'Cilta_Pre', ident.2 = 'Ide_Pre', min.pct = .1, logfc.threshold = 0)

cd8.pre.contrast$gene <- rownames(cd8.pre.contrast)

cd8.pre.contrast$p_val_adj <- ifelse(cd8.pre.contrast$p_val_adj < 10^-10, 10^-10, cd8.pre.contrast$p_val_adj)
cd8.pre.contrast$avg_log2FC <- ifelse(cd8.pre.contrast$avg_log2FC > 2.5, 2.5, cd8.pre.contrast$avg_log2FC)
cd8.pre.contrast$avg_log2FC <- ifelse(cd8.pre.contrast$avg_log2FC < -2.5, -2.5, cd8.pre.contrast$avg_log2FC)

genes <- c('RUNX3', 'RUNX2', 'TGFBR3', 'TGFBR2', 'STAT1', 'NCL', 'BCL2', 'FYN', 'LARP1', 'SATB1', 'NFKBIA', 'BATF', 'STAT3', 
           'B2M', 'CD40LG', 'CD81', 'RAC2', 'GZMA', 'GZMB', 'GNLY', 'IFNG', 'CXCR3', 'LGALS3', 'BAX', 'BAK1', "FAS", 
           'SELL', 'LEF1', 'FOXO1', 'HAVCR2', 'LAG3', 'ENTPD1')

p1 <- volcanoplot_better(cd8.pre.contrast, genes.highlight = genes) 

ggsave('Plots/cd8_preinfusion_deg.pdf', plot = p1, height = 7600, width = 3600*3, dpi = 1200, units = 'px')





### heatmaps ###

t_obj <- subset(t_obj, cells = c(colnames(cd4), colnames(cd8)))

t_obj$tt <- paste0(t_obj$Treatment, '_', t_obj$ttype)

t_obj$tt <- factor(t_obj$tt, levels = c('Ide-pre_CD8', 'Ide-pre_CD4', 'Cilta-pre_CD8', 'Cilta-pre_CD4'))

t_obj$Treatment <- dplyr::recode(t_obj$Treatment, 'Ide' = 'Ide-pre', 'Cilta' = 'Cilta-pre')

saveRDS(t_obj, 'RDS/T_pre.rds')

cd8.core_top <- rbind(data.frame(pathway="T cell Stemness", gene = c('NELL2', 'LEF1', 'SELL', 'IL7R', 'TCF7', 'BACH2', 'CCR7', 'FOXO1')),
                      data.frame(pathway="T cell Exhaustion", gene = c('HAVCR2', 'LAG3', 'ENTPD1', 'PDCD1', 'CTLA4', 'TIGIT', 'TOX')),
                      data.frame(pathway="T cell Cytotoxicity", gene = c('GZMA', 'GZMB', 'GNLY', 'IFNG', 'CXCR3')),
                      data.frame(pathway="TGF-beta Signaling", gene = c('RUNX2', 'RUNX3', 'TGFBR3', 'TGFBR2', 'STAT1')), 
                      data.frame(pathway="IL-2/7/15 up", gene = c('NCL', 'BCL2', 'FYN', 'LARP1')), 
                      data.frame(pathway="IL-1a/b up", gene = c('SATB1', 'NFKBIA', 'BATF', 'STAT3')),
                      data.frame(pathway="T cell Activation", gene = c('B2M', 'CD40LG', 'CD81', 'RAC2')), 
                      data.frame(pathway="Apoptosis", gene = c('LGALS3', 'BAX',  'BAK1', 'FAS'))) 

CairoPDF('Plots/new_plots/cd8_pre_gene_heatmap_large.pdf', width = 5.2, height = 9)

t_obj$Treatment <- factor(t_obj$Treatment, levels = c('Ide-pre', 'Cilta-pre'))

enhancedHeatmap_with_grouping(t_obj, featuredf = cd8.core_top, group.by = 'tt', big.group = 'Treatment', scale = T, column_title = ' ', 
                              column_names_rot = 45, column_labels = c('CD8', 'CD4', 'CD8', 'CD4'), 
                              big.group_cols = c('skyblue','red'), specify_grouping = T, g1_length = 2, g2_length = 2)

dev.off()

trm_features <- c('ATF2', 'ATF3', 'ATF4', 'ATF6', 'BATF', 'BATF3', 'RUNX1', 'RUNX3', 'NFKB1', 'NFKB2', 'KLF2')


CairoPDF('Plots/new_plots/atf_heatmap.pdf', width = 3.5, height = 4.5)

t_obj$ttype <- factor(t_obj$ttype, levels = c('CD8', 'CD4'))

enhancedHeatmap_splitcat(t_obj, features = trm_features, group.by = 'ttype', splitcat = 'Treatment', scale = T,  
                         column_names_rot = 45, x_horizontal = F, cluster_rows = F, cluster_columns = F, cat_order = c('Ide-pre', 'Cilta-pre'))


dev.off()

saveRDS(cd8, 'RDS/CD8_pre.rds')
saveRDS(cd4, 'RDS/CD4_pre.rds')
saveRDS(t_obj, 'RDS/T_obj_pre.rds')
