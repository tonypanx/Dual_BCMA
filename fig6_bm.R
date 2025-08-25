library(Seurat)
library(plyr)
library(dplyr)
library(Azimuth)
library(data.table)
library(forcats)
library(xlsx)
library(EnhancedVolcano)
library(ggtext)
library(clusterProfiler)
library(Cairo)

showtext::showtext_auto()

setwd("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting")
source("/project/jhuangime/tony/Scripts/scRNA_core.R")

color_list <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting/colors.rds")

bm <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/BM_Analysis/RDS/integrated_bm_non_cleaned.rds")

bm$anno1_lab <- str_replace_all(bm$anno1, '_', ' ')

###---------------------------------------------Dimplot--------------------------------------####

p1 <- DimPlot_better(bm, group.by = 'anno1_lab' ,pt.size = .1) + theme(plot.title = element_blank())
ggsave('Fig9/Fig9_all_bm_UMAP.pdf', plot = p1, height = 2300*1.2, width = 3200*1.2, dpi = 800, units = 'px')


p2 <- DimPlot_better(subset(bm, treatment == 'Ide'), group.by = 'anno1_lab',  pt.size = .1) + theme(plot.title = element_blank()) + 
  ggtitle('Ide') + theme(plot.title = element_text(hjust = .5)) + NoLegend()

p3 <- DimPlot_better(subset(bm, treatment == 'Cilta'), group.by = 'anno1_lab',  pt.size = .1) + theme(plot.title = element_blank()) + 
  ggtitle('Cilta') + theme(plot.title = element_text(hjust = .5)) + NoLegend() 

p <- ggarrange(plotlist = list(p2, p3), ncol = 2, nrow = 1)

ggsave('Fig9/Fig9_bm_all_split_UMAP.pdf', plot = p, height = 3000, width = 5400, dpi = 800, units = 'px')

###---------------------------------------------Vlnplt--------------------------------------####

markers <- c('MS4A1', 'CD1C', 'ELANE', 'CD14', 'FCGR3A', 'NCAM1', 'LILRA4', 'MZB1', 'TNFRSF17',  'CD3G', 'CD4', 'CD8A')

p <- VlnPlot(bm, features = c(markers), group.by = 'anno1_lab', stack = T, flip = T) + 
  NoLegend() + theme(strip.text = element_text(face = 'plain'), axis.title.x = element_blank(), axis.title.y = element_text(face = 'bold'))

p$data$feature <- paste0('*', p$data$feature, '*')

p$data$feature <- factor(p$data$feature, levels = c(paste0('*', markers, "*" )))


p <- p + theme(strip.text  =element_markdown())

ggsave('Fig9/Fig9_anno1_markers.pdf', plot = p, height = 3100*1.3, width = 3400*1.3, dpi = 850, units = 'px')

###---------------------------------------------NK Analysis--------------------------------------####


nk <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/BM_Analysis/RDS/NK.rds")
nk$anno2_lab <- revalue(as.character(nk$anno2), c('NK_CD16' = 'CD16-hi', 'NK_CD56' = 'CD56-hi', 'NK_CD16_CX3CR1' = 'CD16-hi CX3CR1-hi', 
                                                  'NK_CD16_CD56' = 'CD56-hi CD16-hi'))
nk$anno2_lab <- factor(nk$anno2_lab, levels = c('CD56-hi', 'CD56-hi CD16-hi', 'CD16-hi', 'CD16-hi CX3CR1-hi'))

p1 <- DimPlot_better(nk, group.by = 'anno2_lab' ,pt.size = .5) + theme(plot.title = element_blank())
ggsave('Fig9/Fig9_nk_bm_UMAP.pdf', plot = p1, height = 2200*1.2, width = 3300*1.2, dpi = 800, units = 'px')

p2 <- FeaturePlot_better(nk, features = c('NCAM1', 'FCGR3A', 'CX3CR1'), ncol = 3, max.cutoff = 'q90', order = T, pt.size = .4)
ggsave('Fig9/Fig9_nk_feature_UMAP.pdf', plot = p2, height = 2200*1.2, width = 5400*1.2, dpi = 800, units = 'px')

###---------------------------------------------anno2 bar--------------------------------------####

ftable <- FetchData(nk, vars = c("anno2_lab", "treatment")) 

dfx <- ftable %>% 
  group_by(treatment) %>%
  dplyr::count(anno2_lab) %>%
  mutate(freq = n/sum(n)) %>%
  ungroup()

p <- ggplot(dfx, aes(x= fct_rev(treatment), y=freq, fill= anno2_lab ))+
  geom_col(width = 0.85, color = 'black', position=position_fill(reverse=T))+ theme_pubr()+ scale_y_continuous(expand=c(0,0.01)) +
  guides(fill = guide_legend(override.aes = list(size=4))) + coord_flip() +
  theme_bw() + labs(y = "Proportion") + theme( axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm"), size=15, face = 'bold'),
                                               axis.title.y = element_blank(), legend.text = element_text(size = 13), axis.text.x = element_text(size = 12, color = 'black'),
                                               axis.text.y = element_text(size = 16, margin = margin(r=.1, unit = 'cm'), color = 'black'), legend.title = element_blank(), 
                                               panel.background = element_rect(color = 'black')) +
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, r=.5, unit = "cm"), legend.box.margin = margin(10, 10, 10, 10),
        plot.title = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(angle = 270, face = "bold", size=10), strip.placement = "outside",panel.grid.major.y = element_blank()) 

ggsave('Fig9/Fig9_anno2_bar_bytreat.pdf', plot = p, height = 2400, width = 6000, dpi = 850, units = 'px')

###---------------------------------------------Signature plot--------------------------------------####


p2 <- FeaturePlot_better(nk, features = c('GO_TNF_Production_1', 'KEGG_NK_Cytotoxicity_1'), ncol = 2, max.cutoff = 'q90', order = T, pt.size = .5, 
                         names = c('TNF Score', 'NK Cytotoxicity Score'), face = 'bold')
ggsave('Fig9/Fig9_nk_sig_UMAP.pdf', plot = p2, height = 2200*1.2, width = 3800*1.2, dpi = 800, units = 'px')

###---------------------------------------------Vln plot--------------------------------------####

p3 <- vsplot_many(nk, gene = c('FCGR3A', 'NCAM1'), pt.size = 0, groups = 'treatment', colors = color_list$treat, nrows = 1, ncols = 2, isAssay5 = T)

ggsave('Fig9/Fig9_nk_genes.pdf', plot = p3, height = 2200, width = 4200, dpi = 850, units = 'px')

###---------------------------------------------heatmap--------------------------------------####
nk.features <- c('CD160', 'FCER1G', 'FCGR3A', 'CD247', 'CHST2', 'SPON2', 'FGFBP2', 'IGFBP7', 'MYOM2', 'CLIC3', 'AKR1C3', 'PRF1', 'KLRB1', 'GZMB')


CairoPDF('Fig9/Fig9_nk_hm.pdf', width = 6.7, height = 3.5)

enhancedHeatmap_splitcat(nk, features = nk.features, scale = T, group.by = 'anno2_lab', 
                         cluster_rows = F, column_names_rot = 45, splitcat = 'treatment', cluster_columns = F, 
                         column_order = nk.features)

dev.off()

###---------------------------------------------Mono Analysis--------------------------------------####
mono <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/BM_Analysis/RDS/mono.rds")

p1 <- DimPlot_better(mono, group.by = 'anno1' ,pt.size = .5) + theme(plot.title = element_blank())
ggsave('Fig9/Fig9_mono_bm_UMAP.pdf', plot = p1, height = 2200*1.2, width = 3000*1.2, dpi = 800, units = 'px')

p2 <- VlnPlot(mono, features = c('CD14', 'FCGR3A'), group.by = 'anno1', stack = T, flip = T) + 
  NoLegend() + theme(strip.text = element_text(face = 'plain'), axis.title.x = element_blank(), axis.title.y = element_text(face = 'bold'))

p2$data$feature <- paste0('*', p2$data$feature, '*')

p2 <- p2 + theme(strip.text  =element_markdown())

ggsave('Fig9/Fig9_mono_anno1_markers.pdf', plot = p2, height = 2000, width = 3000, dpi = 850, units = 'px')

###---------------------------------------------Signature plot--------------------------------------####


p1 <- FeaturePlot_better(mono, features = c('Hallmark_Inflam_1'),max.cutoff = 'q90', order = T, pt.size = .5, 
                         face = 'bold', split.by = 'treatment') & scale_color_viridis(
                           breaks = pretty_breaks(n = 4),
                           labels = function(x) signif(x, 1) 
                         )



ggsave('Fig9/Fig9_mono_sig_umap.pdf', plot = p1, height = 3000*1.1, width = 5200*1.1, dpi = 850, units = 'px')

###---------------------------------------------inflam vsplot--------------------------------------####

p2 <- vsplot(subset(mono, anno1 == 'CMono'), groups = 'treatment', colors = color_list$treat, gene = 'Hallmark_Inflam_1', pt.size = 0, feature_type = 'score') + 
  ggtitle('Inflammatory Response Score')

ggsave('Fig9/Fig9_cmono_inflam.pdf', plot = p2, height = 2200*1.1, width = 2800*1.1, dpi = 850, units = 'px')

###---------------------------------------------mono heatmap--------------------------------------####
mono.features <- c('C3AR1', 'NAMPT', 'PIK3R5', 'TLR2', 'HIF1A', 'IL4R', 'KLF6', 'NFKB1', 'CD69', 'HBEGF', 'MYC', 'LDLR', 'NLRP3', 'IL1B')

CairoPDF('Fig9/Fig9_mono_hm.pdf', width = 5.3, height = 2.6)

enhancedHeatmap_splitcat(mono, features = mono.features, scale = T, group.by = 'anno1', 
                         cluster_rows = F, column_names_rot = 45, splitcat = 'treatment', cluster_columns = F, 
                         column_order = mono.features)

dev.off()

###---------------------------------------------mono gt--------------------------------------####

gene_trajectory <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/BM_Analysis/Data/mono_gene_traj_res.rds")

early_genes <- c('CD14', 'SELL', 'CCR2')
late_genes <- c('FCGR3A', 'SPN', 'CX3CR1', 'VMO1', 'SH2D1B', 'CSF1R', 'MS4A7',  'PTGER4', 'C3AR1', 'CXCL16', 'NR4A1', 'FCGR3B', 'CDKN1C', 'CASP5' )

gene_trajectory$goi <- gene_trajectory$gene
gene_trajectory$goi <- case_when(gene_trajectory$goi %!in% c(early_genes, late_genes) ~ '', 
                                 TRUE ~ gene_trajectory$goi)

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(0.1, 'npc')
)


p <- ggplot(gene_trajectory, aes(x = DM_1, y = DM_2, color = lab, label = goi)) + geom_point() + 
  geom_text_repel(max.overlaps = 999) + xlim(c(-.02, .07)) + ylim(c(-.009, .016))+ theme_minimal()+
  theme(legend.position = 'top') + guides(color = guide_legend(title = 'Trajectory')) & 
  guides(x = axis, y = axis) &
  theme(axis.line = element_line(arrow=arrow(type = 'closed', length = unit(6, 'pt'))), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = 'right',
        axis.title = element_text(hjust = 0, size = 15), panel.grid = element_blank(),
        plot.title = element_text( margin = margin(25, 25, 25, 25))) & 
  ylab('DM2') & xlab('DM1')

ggsave('Fig9/Fig9_mono_traj.pdf', plot = p, height = 2200*2, width = 2800*2, dpi = 850, units = 'px')

p2 <- vsplot_many(subset(mono, anno1 == 'NMono'), gene = c('T2_1', 'T3_1'), groups = 'treatment', nrows = 1, ncols = 2, pt.size = 0, 
                  names = c('T2 Score', 'T3 Score'), feature_type = 'score', colors = color_list$treat)

ggsave('Fig9/Fig9_nomono_tscores.pdf', plot = p2, height = 2400, width = 5000, dpi = 850, units = 'px')


