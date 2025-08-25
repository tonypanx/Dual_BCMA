library(Seurat)
library(scRepertoire)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(EnhancedVolcano)
library(scplotter)
library(forcats)
library(WGCNA)
library(hdWGCNA)
library(RColorBrewer)
library(corrplot)
library(ggcorrplot)

setwd("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting")
source("/project/jhuangime/tony/Scripts/scRNA_core.R")

showtext::showtext_auto()

set.seed(1234)

early_post <- readRDS('/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/RDS/early_post_ide_cilta_carpos.rds')
cd8 <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/RDS/cilta_car_cd8_pb_bm_merge.rds")
cd8_te <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/RDS/cilta_car_te_cd8_pb_bm_merge.rds")

color_list <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting/colors.rds")


###---------------------------------------------Dimplot/Featureplot of Cilta D14 and Ide D9--------------------------------------####

p1 <- DimPlot_better(early_post, group.by = 'treatment', cols = color_list$treat, pt.size = .4) + theme(plot.title = element_blank())

p2 <- FeaturePlot_better(early_post, features = c('Ide-gene', 'Cilta-comb'), max.cutoff = 'q90', order = T, names = c('Ide CAR', 'Cilta CAR'))

ggsave('Fig5b/Fig5b_peak_CAR_UMAP.pdf', plot = p1, height = 2400, width = 3100, dpi = 850, units = 'px')
ggsave('Fig5b/Fig5b_peak_CAR_gene.pdf', plot = p2, height = 2800, width = 4500, dpi = 850, units = 'px')

###---------------------------------------------Volcano comparison of Cilta D14 vs Ide D9--------------------------------------####

markers.peak <- FindMarkers(early_post, ident.1 = 'Cilta', ident.2 = 'Ide', group.by = 'treatment', logfc.threshold = 0)
markers.peak_filt <- filter_markers(clean_markers(markers.peak))

markers.peak$p_val_adj[markers.peak$p_val_adj < 10^-50] <- 10^-50

early_treat.genes <- c('KLRC2', 'KLRC1', 'GZMK', 'FYN', 'CD27', 'FCGR3A', 'CD8A', 'CCL5', 'HAVCR2', 'CX3CR1', 'PRF1', 'NKG7', 'TIGIT', 'GZMA',
                       'FYN', 'NFKB1', 'CD69', 'JUNB', 'JUND', 'BHLHE40', 'JUN', 'ICOS', 'IL7R', 'FOS', 'TNF', 'KLRB1', 'CD4', 'IL2RA', 'TNFRSF4', 'CD40LG', 'FOXP3', 'TNFRSF18' )

p <- VolcanoPlot(markers.peak, x = 'avg_log2FC', y = 'p_val_adj', label_size = 4, flip_negatives = T, 
                 labels = early_treat.genes, highlight = early_treat.genes, y_cutoff = .05,
                 xlab =  bquote(paste('Average ', log[2], 'FC')), ylab = "-log<sub>10</sub>(adjusted p-value)")  + 
  theme(legend.position = "none", axis.title.y = ggtext::element_markdown()) 

ggsave('Fig5b/Fig5b_peak_volcano.pdf', plot = p, height = 5400, width = 5400, dpi = 850, units = 'px')

###---------------------------------------------CAR Fitness/Migration comparison--------------------------------------####

p1 <- vsplot_many(early_post, groups = 'treatment', gene = c('CAR_Fit_Up_1', 'CAR_Migratory_1'), nrows =  1, ncols = 2, pt.size = 0, colors = color_list$treat, 
                  feature_type = 'score', names = c('CAR Fit Score', 'CAR Migratory Score'))

ggsave('Fig5b/Fig5b_car_fit.pdf', plot = p1, height = 2200, width = 5200, dpi = 850, units = 'px')

###---------------------------------------------Cilta Post CAR UMAP--------------------------------------####
cd8$time <- revalue(cd8$sample, c('Cilta_BM_CAR_1' = 'D28_BM', 'Cilta_BM_CAR_2' = 'D28_BM', 'Ciltacel_D14_CAR' = 'D14', 
                                  'Ciltacel_D21_CAR' = 'D21', 'Ciltacel_D28_CAR' = 'D28'))

cd8$time <- factor(cd8$time, levels = c('D14', 'D21', 'D28', 'D28_BM'))

cd8$anno1 <- as.character(cd8$anno1)
cd8$anno1 <- revalue(cd8$anno1, c('CD8_Tem' = 'CD8_EM'))

cd8$anno1_short <- str_replace_all(cd8$anno1, 'CD8_', '')
cd8$anno1_short <- factor(cd8$anno1_short, levels = c('Prolif', 'EM', 'ISG', 'IFNG', 'TE', 'IL'))

cd8$anno2 <- as.character(cd8$anno2)
cd8$anno2 <- revalue(cd8$anno2, c('CD8_Tem' = 'CD8_EM'))

cd8$anno2_short <- str_replace_all(cd8$anno2, 'CD8_', '')
cd8$anno2_short <- factor(cd8$anno2_short, levels = c('Prolif', 'EM', 'ISG', 'IFNG', 'TE1', 'TE2', 'TE3', 'IL'))

saveRDS(cd8, "/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/RDS/cilta_car_cd8_pb_bm_merge.rds")

p1 <- DimPlot_better(cd8, group.by = 'anno1_short', cols = color_list$cd8_anno1[levels(cd8$anno1_short)], pt.size = .1) + theme(plot.title = element_blank())
ggsave('Fig5b/Fig5b_cilta_post_CAR_anno1_UMAP.pdf', plot = p1, height = 2400*1.2, width = 3000*1.2, dpi = 800, units = 'px')

p2 <- DimPlot_better(cd8, group.by = 'time', cols = color_list$time[levels(cd8$time)], pt.size = .1) + theme(plot.title = element_blank()) + 
  scale_color_manual(values = color_list$time[levels(cd8$time)], labels = c('D28_BM' = 'D28 BM'))
ggsave('Fig5b/Fig5b_cilta_post_CAR_time_UMAP.pdf', plot = p2, height = 2400*1.2, width = 3000*1.2, dpi = 800, units = 'px')

p3 <- DimPlot_better(cd8, group.by = 'anno2_short', cols = color_list$cd8_anno2[levels(cd8$anno2_short)], pt.size = .1) + theme(plot.title = element_blank())
ggsave('Fig5b/Fig5b_cilta_post_CAR_anno2_UMAP.pdf', plot = p3, height = 2400*1.2, width = 3000*1.2, dpi = 800, units = 'px')

###----------------------------------------------CD8 Post anno1 proportions--------------------------------------####
ftable <- FetchData(cd8, vars = c("anno1_short", "time")) 

dfx <- ftable %>% 
  group_by(time) %>%
  dplyr::count(anno1_short) %>%
  mutate(freq = n/sum(n)) %>%
  ungroup()

dfx$time <- as.character(dfx$time)
dfx$time <- revalue(dfx$time, c('D28_BM' = 'D28 BM'))
dfx$time <- factor(dfx$time, levels = c('D14', 'D21', 'D28', 'D28 BM'))

p <- ggplot(dfx, aes(x= fct_rev(time), y=freq, fill= anno1_short ))+
  geom_col(width = 0.85, color = 'black', position=position_fill(reverse=T))+ theme_pubr()+ scale_y_continuous(expand=c(0,0.01)) +
  guides(fill = guide_legend(override.aes = list(size=4))) + coord_flip() +
  theme_bw() + labs(y = "Proportion") + theme( axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm"), size=15, face = 'bold'),
                                               axis.title.y = element_blank(), legend.text = element_text(size = 13), axis.text.x = element_text(size = 12, color = 'black'),
                                               axis.text.y = element_text(size = 16, margin = margin(r=.1, unit = 'cm'), color = 'black'), legend.title = element_blank(), 
                                               panel.background = element_rect(color = 'black')) +
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, r=.5, unit = "cm"), legend.box.margin = margin(10, 10, 10, 10),
        plot.title = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(angle = 270, face = "bold", size=10), strip.placement = "outside",panel.grid.major.y = element_blank()) + scale_fill_manual(values = color_list$cd8_anno1[levels(cd8$anno1_short)]) 

ggsave('Fig5b/Fig5b_anno1_bar.pdf', plot = p, height = 2400, width = 6000, dpi = 850, units = 'px')

###----------------------------------------------CD8 Post anno2 proportions--------------------------------------####
ftable <- FetchData(cd8, vars = c("anno2_short", "time")) 

dfx <- ftable %>% 
  group_by(time) %>%
  dplyr::count(anno2_short) %>%
  mutate(freq = n/sum(n)) %>%
  ungroup()

dfx$time <- as.character(dfx$time)
dfx$time <- revalue(dfx$time, c('D28_BM' = 'D28 BM'))
dfx$time <- factor(dfx$time, levels = c('D14', 'D21', 'D28', 'D28 BM'))

p <- ggplot(dfx, aes(x= fct_rev(time), y=freq, fill= anno2_short ))+
  geom_col(width = 0.85, color = 'black', position=position_fill(reverse=T))+ theme_pubr()+ scale_y_continuous(expand=c(0,0.01)) +
  guides(fill = guide_legend(override.aes = list(size=4))) + coord_flip() +
  theme_bw() + labs(y = "Proportion") + theme( axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm"), size=15, face = 'bold'),
                                               axis.title.y = element_blank(), legend.text = element_text(size = 13), axis.text.x = element_text(size = 12, color = 'black'),
                                               axis.text.y = element_text(size = 16, margin = margin(r=.1, unit = 'cm'), color = 'black'), legend.title = element_blank(), 
                                               panel.background = element_rect(color = 'black')) +
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, r=.5, unit = "cm"), legend.box.margin = margin(10, 10, 10, 10),
        plot.title = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(angle = 270, face = "bold", size=10), strip.placement = "outside",panel.grid.major.y = element_blank()) + scale_fill_manual(values = color_list$cd8_anno2[levels(cd8$anno2_short)]) 

ggsave('Fig5b/Fig5b_anno2_bar.pdf', plot = p, height = 2400, width = 6000, dpi = 850, units = 'px')

###----------------------------------------------CD8 Post anno1 markers--------------------------------------####

rna.markers <- c('Cilta-comb', 'CD8A', 'MKI67', 'IL7R', "SELL", 'GZMK', 'ISG15', 'IFNG', 'CX3CR1')
adt.markers <- c('CD127', 'CD62L', 'CX3CR1.1')

p <- VlnPlot(cd8, features = c(rna.markers, adt.markers), group.by = 'anno1_short', stack = T, flip = T) + 
  NoLegend() + theme(strip.text = element_text(face = 'plain'), axis.title.x = element_blank(), axis.title.y = element_text(face = 'bold'))

p$data$feature <- paste0('*', p$data$feature, '*')

p$data$feature <- revalue(p$data$feature, c('*adt_CD127*' = 'CD127',
                                            '*adt_CD62L*' = 'CD62L',
                                            '*adt_CX3CR1.1*' = 'CX3CR1',
                                            '*adt_KLRG1.1*' = 'KLRG1', 
                                            '*Cilta-comb*' = '*Cilta CAR*'))

p$data$feature <- factor(p$data$feature, levels = c(paste0('*', c('Cilta CAR', rna.markers[2:length(rna.markers)]), '*'), c('CD127', 'CD62L', 'CX3CR1', 'KLRG1') ))

p <- p + theme(strip.text  =element_markdown())

ggsave('Fig5b/Fig5b_anno1_markers.pdf', plot = p, height = 5400, width = 3000, dpi = 850, units = 'px')

###----------------------------------------------CD8 Post anno2 markers--------------------------------------####
rna.markers <- c('KLRG1', 'CX3CR1', 'TBX21', 'SPON2', 'ZEB2', 'GZMK', 'CXCR4')

p <- VlnPlot(subset(cd8, anno1_short == 'TE'), features = c(rna.markers), group.by = 'anno2_short', stack = T, flip = T) + 
  NoLegend() + theme(strip.text = element_text(face = 'plain'), axis.title.x = element_blank(), axis.title.y = element_text(face = 'bold'))

p$data$feature <- paste0('*', p$data$feature, '*')


p$data$feature <- factor(p$data$feature, levels = paste0('*', rna.markers, '*'))

p <- p + theme(strip.text  =element_markdown())

ggsave('Fig5b/Fig5b_anno2_markers.pdf', plot = p, height = 4000, width = 3000, dpi = 850, units = 'px')

###----------------------------------------------CD8 Post time markers--------------------------------------####

rna.markers <- c('PRF1', 'CX3CR1', 'FCGR3A', 'CXCR3', 'TIGIT',  'NR4A2', 'CXCR4', 'CD69', 'RUNX3', 'NELL2', 'FOXO1', 'NFKB1', 'Cilta-comb')

p <- VlnPlot(cd8, features = c(rna.markers), group.by = 'time', stack = T, flip = T) + 
  NoLegend() + theme(strip.text = element_text(face = 'plain'), axis.title.x = element_blank(), axis.title.y = element_text(face = 'bold')) + 
  scale_x_discrete(labels = c('D14', 'D21', 'D28', 'D28 BM'))

p$data$feature <- paste0('*', p$data$feature, '*')

p$data$feature <- revalue(p$data$feature, c('*Cilta-comb*' = '*Cilta CAR*'))

p$data$feature <- factor(p$data$feature, levels = c(paste0('*', c( rna.markers[1:(length(rna.markers)-1)], 'Cilta CAR' ), '*')))


p <- p + theme(strip.text  =element_markdown())

ggsave('Fig5b/Fig5b_time_markers.pdf', plot = p, height = 5400, width = 3000, dpi = 850, units = 'px')

###----------------------------------------------wgcna modules for TE--------------------------------------####

MEs <- GetModuleScores(cd8_te, wgcna_name=cd8_te@misc$active_wgcna)
MEs <- as.matrix(MEs)
res <- Hmisc::rcorr(x = MEs, type = 'spearman')
rownames(res$r) <- str_replace_all(rownames(res$r), 'TE_', '') 
colnames(res$r) <- str_replace_all(colnames(res$r), 'TE_', '') 

p <- ggcorrplot(res$r, type = 'upper', show.diag = T, legend.title = 'Spearman') + scale_x_discrete(position='top') + 
  theme(axis.text.x= element_text(vjust = .1, hjust = -.1), panel.grid.major = element_blank()) 
                                                                                            
ggsave('Fig5b/Fig5b_te_module_corrplot.pdf', plot = p, height = 3600, width = 3600, dpi = 850, units = 'px')

###----------------------------------------------vln plot for key modules--------------------------------------####
cd8_te$condition2 <- as.character(cd8_te$condition)
cd8_te$condition2 <- str_replace_all(cd8_te$condition2, 'PB_', '')
cd8_te$condition2 <- str_replace_all(cd8_te$condition2, '_', ' ')
cd8_te$condition2 <- factor(cd8_te$condition2, levels = c('D14', 'D21', 'D28', 'BM D28'))

p <- vsplot_many(cd8_te, plot_signif=F, gene = c('TE_M1', 'TE_M2', 'TE_M5', 'TE_sig_1'), pt.size = 0,  names = c('M1 Score', 'M2 Score', 'M5 Score', 'Effector Score'),
            groups = 'condition2', isAssay5 = T, nrows = 1, ncols = 4, step = .1, feature_type = 'score') + theme(legend.position = 'none')

ggsave('Fig5b/Fig5b_te_mod_scores.pdf', plot = p, height = 2200, width = 7000, dpi = 850, units = 'px')



###----------------------------------------------modules top genes--------------------------------------####
te_modules <- GetModules(cd8_te)

p1 <- dplyr::arrange(dplyr::filter(te_modules, module == 'TE_M1'), desc(kME_TE_M1))[1:75, ] %>% 
  ggplot(aes(x=reorder(gene_name, kME_TE_M1), y = kME_TE_M1))+ geom_bar(stat = 'identity', fill = 'turquoise') + coord_flip() + theme_pubr() + 
  ylab('kME') + xlab('Gene') + theme(axis.text.y = element_text(size = 10), plot.title = element_text(size = 15, hjust = .5)) + ggtitle('M1 Genes')

p2 <- dplyr::arrange(dplyr::filter(te_modules, module == 'TE_M2'), desc(kME_TE_M2))[1:75, ] %>% 
  ggplot(aes(x=reorder(gene_name, kME_TE_M2), y = kME_TE_M2))+ geom_bar(stat = 'identity', fill = 'blue') + coord_flip() + theme_pubr() + 
  ylab('kME') + xlab('Gene') + theme(axis.text.y = element_text(size = 10), plot.title = element_text(size = 15, hjust = .5)) + ggtitle('M2 Genes')

p3 <- dplyr::arrange(dplyr::filter(te_modules, module == 'TE_M5'), desc(kME_TE_M5)) %>% 
  ggplot(aes(x=reorder(gene_name, kME_TE_M5), y = kME_TE_M5))+ geom_bar(stat = 'identity', fill = 'brown') + coord_flip() + theme_pubr() + 
  ylab('kME') + xlab('Gene') + theme(axis.text.y = element_text(size = 10), plot.title = element_text(size = 15, hjust = .5)) + ggtitle('M5 Genes')

ggsave('Fig5b/Fig5b_te_m1_genes.pdf', plot = p1, height = 8500, width = 6000, dpi = 850, units = 'px')
ggsave('Fig5b/Fig5b_te_m2_genes.pdf', plot = p2, height = 8500, width = 6000, dpi = 850, units = 'px')
ggsave('Fig5b/Fig5b_te_m5_genes.pdf', plot = p3, height = 7500, width = 6000, dpi = 850, units = 'px')

###----------------------------------------------top gene heat map--------------------------------------####

trm_genes <- c('TNFAIP3', 'NR4A2', 'CXCR4', 'CD69', 'RUNX3', 'SLC2A3', 'TCF7', 'NELL2')
mapk_genes <- c('PIK3R1', 'DUSP1', 'JUN', 'DUSP2', 'FYN', 'RASA2', 'DUSP5', 'RASA3')
actin_genes <- c('ACTB', 'ACTG1', 'PFN1', 'CORO1A', 'CFL1', 'MYL6')
effector_genes <- c('NKG7', 'PRF1', 'GZMA', 'GZMB', 'CTSW', 'GZMH', 'CX3CR1',  'EMP3',  'CXCR3')
effector_2_genes <- c('FGFBP2',   'FCGR3A', 'SPON2',  'TBX21',  'LAIR2')

cd8_te$anno2_short <- str_replace_all(cd8_te$anno2, 'CD8_', '')

cd8_te$condition2 <- as.character(cd8_te$condition2)
cd8_te$anno2_short <- as.character(cd8_te$anno2_short)

ht_opt$HEATMAP_LEGEND_PADDING =unit(.6, 'cm')

CairoPDF('Fig5b/Fig5b_m1_hm.pdf', width = 6, height = 4)

enhancedHeatmap_splitcat(cd8_te, features = c(trm_genes, mapk_genes), scale = T, group.by = 'anno2_short', 
                         cluster_rows = F, column_names_rot = 45, splitcat = 'condition2', cluster_columns = F, 
                         column_order = c(trm_genes, mapk_genes), cat_order  = c('D14', 'D21', 'D28', 'BM D28'))



dev.off()

CairoPDF('Fig5b/Fig5b_eff_hm.pdf', width = 5.5, height = 4.5)

enhancedHeatmap_splitcat(cd8_te, features = c(effector_genes, effector_2_genes), scale = T, group.by = 'anno2_short', 
                               cluster_rows = F, column_names_rot = 45, splitcat = 'condition2', cluster_columns = F, 
                               column_order = c(effector_genes, effector_2_genes), cat_order  = c('D14', 'D21', 'D28', 'BM D28'))

dev.off()



###----------------------------------------------IP enrichment--------------------------------------####

gsea_df <- readRDS('/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/data/post_cd8_ip_enrichment.rds')
gsea_df$condition <- revalue(gsea_df$condition, c('PB_D14' = 'D14', 'PB_D21' = 'D21', 'PB_D28' = 'D28', 'BM_D28' = 'BM D28'))

pa <- ggplot(dplyr::filter(gsea_df, pathway == 'IP_Pos'), aes(x = reorder(condition, NES), y = NES)) +
  geom_segment(aes(xend = condition, yend = 0), size = 1) +
  geom_point(aes(size = size, fill = padj), shape = 21, fill = color_list$treat['Cilta'], color = "black") +
  theme_minimal() + coord_flip() + 
  labs(x = "Condition", y = "NES", title = "IP Cilta vs Ide UP Enrichment", size = "Size") + theme_pubr()+
  theme(legend.position = "right", axis.title.y = element_blank(), plot.title = element_text(hjust = .5), axis.title = element_text(face = 'bold')) + 
  scale_size(range = c(3, 6)) 

ggsave('Fig5b/Fig5b_post_ip_up_enrich.pdf', plot = pa, height = 2400, width = 4000, dpi = 850, units = 'px')

pa <- ggplot(dplyr::filter(gsea_df, pathway == 'IP_Neg'), aes(x = reorder(condition, NES), y = NES)) +
  geom_segment(aes(xend = condition, yend = 0), size = 1) +
  geom_point(aes(size = size, fill = padj), shape = 21, fill = color_list$treat['Ide'], color = "black") +
  theme_minimal() + coord_flip() + 
  labs(x = "Condition", y = "NES", title = "IP Cilta vs Ide DN Enrichment", size = "Size") + theme_pubr()+
  theme(legend.position = "right", axis.title.y = element_blank(), plot.title = element_text(hjust = .5), axis.title = element_text(face = 'bold')) + 
  scale_size(range = c(3, 6)) 

ggsave('Fig5b/Fig5b_post_ip_down_enrich.pdf', plot = pa, height = 2400, width = 4000, dpi = 850, units = 'px')

###----------------------------------------------IP enrichment vlnplot--------------------------------------####
cd8_te$condition2 <- factor(cd8_te$condition2, levels = c('D14', 'D21', 'D28', 'BM D28'))
p1 <- vsplot(cd8_te, groups = 'condition2', gene = c('Cilta_IP_sig_1'), pt.size = 0, colors = as.character(color_list$time[c('D14', 'D21', 'D28', 'D28_BM')]), 
                  feature_type = 'score', comparisons = list(c('D28', 'BM D28'))) + ggtitle('Cilta IP Score') + theme(plot.title = element_text(hjust = .5))

ggsave('Fig5b/Fig5b_post_ip_score.pdf', plot = p1, height = 2200, width = 3200, dpi = 850, units = 'px')

####-------------------------saving----------------------------------------------##

saveRDS(cd8_te, "/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/RDS/cilta_car_te_cd8_pb_bm_merge.rds")

