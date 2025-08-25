library(Seurat)
library(plyr)
library(dplyr)
library(SeuratDisk)
library(ggplot2)
library(viridis)
library(scales)
library(stringr)
library(fgsea)
library(scplotter)
library(forcats)
library(extrafont)
library(ragg)
library(ggtext)
library(Cairo)

setwd("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting")
source("/project/jhuangime/tony/Scripts/scRNA_core.R")

font_import()

showtext::showtext_auto()

set.seed(1234)

cd8 <- readRDS('/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/RDS/KK_ip_combined_reprocessed_nomt_cd8.rds')
cd4 <- readRDS('/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/RDS/KK_ip_combined_reprocessed_nomt_cd4.rds')

color_list <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting/colors.rds")

markers <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/data/cd8_treatment_markers_nomt.rds") 
paths <- readRDS( "/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/data/cd8_treatment_paths_nomt.rds")
csig <- readRDS( "/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/data/cd8_csig.rds")

###----------------------------------------------Reharmonize labels--------------------------------------####

cd8$anno1 <- str_replace_all(cd8$anno1, 'Tem', 'EM')
cd8$anno1 <- str_replace_all(cd8$anno1, 'Teff', 'TE')

cd8$anno2 <- str_replace_all(cd8$anno2, 'Tem', 'EM')
cd8$anno2 <- str_replace_all(cd8$anno2, 'Teff', 'TE')

cd8$anno1 <- as.character(cd8$anno1)
cd8$anno1 <- str_replace_all(cd8$anno1, 'CD8_PEM', 'CD8_Prolif_EM')
cd8$anno1 <- str_replace_all(cd8$anno1, 'CD8_PTE', 'CD8_Prolif_TE')
cd8$anno1 <- factor(cd8$anno1, levels = c('CD8_Prolif_EM', 'CD8_Prolif_TE', 'CD8_TE'))

cd8$treatment <- cd8$treatment_short
cd8$treatment <- revalue(cd8$treatment, c('Cil' = 'Cilta'))

cd4$treatment <- cd4$treatment_short
cd4$treatment <- revalue(cd4$treatment, c('Cil' = 'Cilta'))

cd8$anno1_short <- str_replace_all(cd8$anno1, 'CD8_', '')
cd8$anno1_short <- factor(cd8$anno1_short, levels = c('Prolif_EM', 'Prolif_TE', 'TE'))

saveRDS(cd8, '/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/RDS/KK_ip_combined_reprocessed_nomt_cd8.rds')
saveRDS(cd4, '/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/RDS/KK_ip_combined_reprocessed_nomt_cd4.rds')

###----------------------------------------------CD4/CD8 ratio plot--------------------------------------####

t_props <- data.frame(treatment = c('Ide', 'Ide', 'Cilta', 'Cilta'), type = c('CD8', 'CD4', 'CD8', 'CD4') , 
                      counts = c(1458, 8211, 1169, 855))

dfx <- t_props %>% 
  group_by(treatment) %>%
  mutate(prop = counts/sum(counts))

dfx$type <- factor(dfx$type, levels = c('CD4', 'CD8'))

p <- ggplot(dfx, aes(x= treatment, y=prop, fill= type ))+
  geom_col(width = 0.85, color = 'black', position=position_stack(reverse=F))+ theme_pubr()+ scale_y_continuous(expand=c(0,0.01)) +
  guides(fill = guide_legend(override.aes = list(size=8)))  +
  theme_void() + labs(y = "Proportion") + theme( axis.title.y = element_text(margin = margin(t = 0.5, b = 0.5, r = .5, unit = "cm"), size=15, angle = 90),
                                                 axis.title.x = element_blank(), legend.text = element_text(size = 18), axis.text.x = element_text(size = 16),
                                                 axis.text.y = element_text(size = 12), legend.title = element_blank()) +
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, r=.5, unit = "cm"), legend.box.margin = margin(10, 10, 10, 10),
        plot.title = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(angle = 270, face = "bold", size=10), strip.placement = "outside",panel.grid.major.y = element_blank()) + 
  guides(fill = guide_legend(override.aes=list(size = 2))) + scale_fill_manual(values = color_list$ttype)

ggsave('Fig2/Fig2_ttype_bar.pdf', plot = p, height = 3600, width = 2400, dpi = 500, units = 'px')

###----------------------------------------------CD8 Dimplot--------------------------------------####\
cd8$anno1_short <- as.character(cd8$anno1_short)
p <- DimPlot_better(cd8,label = F, group.by = 'anno1_short', cols = color_list$cd8_anno1) + theme(plot.title = element_blank())

ggsave('Fig2/Fig2_anno1_UMAP.pdf', plot = p, height = 2400, width = 3300, dpi = 850, units = 'px')

###----------------------------------------------CD8 Dimplot by treatment--------------------------------------####
p1 <- DimPlot_better(subset(cd8, treatment == 'Cilta'), pt.size = .7, label = F, group.by = 'anno1_short', cols = color_list$cd8_anno1) + theme(plot.title = element_blank()) + NoLegend()
p2 <- DimPlot_better(subset(cd8, treatment == 'Ide'), pt.size = .7, label = F, group.by = 'anno1_short', cols = color_list$cd8_anno1) + theme(plot.title = element_blank()) + NoLegend()
p3 <- DimPlot_better(cd8,label = F, group.by = 'treatment', cols = color_list$treat) + theme(plot.title = element_blank())


p <- ggarrange(plotlist = list(p1, p2), ncol = 2, nrow = 1)
ggsave('Fig2/Fig2_anno1_UMAP_splittreat.pdf', plot = p, height = 2400, width = 5000, dpi = 850, units = 'px')

ggsave('Fig2/Fig2_UMAP_treatsplit.pdf', plot = p3, height = 2400, width = 3300, dpi = 850, units = 'px')

###----------------------------------------------CD8 proportions, anno1--------------------------------------####
ftable <- FetchData(cd8, vars = c("anno1_short", "treatment")) 

dfx <- ftable %>% 
  group_by(treatment) %>%
  dplyr::count(anno1_short) %>%
  mutate(freq = n/sum(n)) %>%
  ungroup()

p <- ggplot(dfx, aes(x= fct_rev(treatment), y=freq, fill= anno1_short ))+
  geom_col(width = 0.85, color = 'black', position=position_fill(reverse=T))+ theme_pubr()+ scale_y_continuous(expand=c(0,0.01)) +
  guides(fill = guide_legend(override.aes = list(size=4))) + coord_flip() +
  theme_bw() + labs(y = "Proportion") + theme( axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm"), size=15, face = 'bold'),
                                                 axis.title.y = element_blank(), legend.text = element_text(size = 13), axis.text.x = element_text(size = 12, color = 'black'),
                                                 axis.text.y = element_text(size = 16, margin = margin(), color = 'black'), legend.title = element_blank(), 
                                               panel.background = element_rect(color = 'black')) +
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, r=.5, unit = "cm"), legend.box.margin = margin(10, 10, 10, 10),
        plot.title = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(angle = 270, face = "bold", size=10), strip.placement = "outside",panel.grid.major.y = element_blank()) + scale_fill_manual(values = color_list$cd8_anno1)

ggsave('Fig2/Fig2_anno1_bar.pdf', plot = p, height = 2400, width = 6000, dpi = 850, units = 'px')

###----------------------------------------------Vln markers, anno1--------------------------------------####

rna.markers <- c('CD8A', 'MKI67', 'SELL', 'CCR7', 'LEF1', 'CD27', 'GZMA', 'TNF', 'GNLY', 'PRF1')
adt.markers <- c('CD127', 'CD62L', 'CD27.1', 'CD45RA', 'KLRG1.1', 'CD45RO')


p <- VlnPlot(cd8, features = c(rna.markers, adt.markers), group.by = 'anno1_short', stack = T, flip = T) + 
  NoLegend() + theme(strip.text = element_text(face = 'plain'), axis.title.x = element_blank(), axis.title.y = element_text(face = 'bold'))

p$data$feature <- paste0('*', p$data$feature, '*')

p$data$feature <- revalue(p$data$feature, c('*adt_KLRG1.1*' = 'KLRG1',
                                            '*adt_CD127*' = 'CD127',
                                            '*adt_CD62L*' = 'CD62L',
                                            '*adt_CD45RA*' = 'CD45RA',
                                            '*adt_CD45RO*' = 'CD45RO',
                                            '*adt_CD27.1*' = 'CD27'))

p$data$feature <- factor(p$data$feature, levels = c(paste0('*', rna.markers, '*'), c('CD127', 'CD62L', 'CD27', 'CD45RA', 'KLRG1', 'CD45RO') ))

p <- p + theme(strip.text  =element_markdown())

ggsave('Fig2/Fig2_anno1_markers.pdf', plot = p, height = 5400, width = 3000, dpi = 850, units = 'px')

###----------------------------------------------Vln markers, contrast--------------------------------------####

rna.eff <- c('GNLY', 'NKG7', 'GZMA',  'PRF1', 'CXCR3', 'CCL5', 'KLRG1', 'TBX21')
rna.mem <- c('SELL', 'LEF1', 'CD27','FOXO1')
rna.tex <- c('HAVCR2', 'LAG3', 'ENTPD1', 'TIGIT')
rna.per <- c('BATF3', 'BATF', 'BCL2', 'BAX')

adt.mem <- c( 'CD62L', 'CD27.1', 'CD127', 'CCR7.1', 'CD45RA' )
adt.eff <- c('CXCR3.1', 'KLRG1.1', 'CD57', 'CD45RO')
adt.tex <- c('PD-1', 'TIM-3', 'LAG-3', 'CD39', 'TIGIT.1')

p1 <- vsplot_many(cd8, groups = 'treatment', gene = rna.eff, nrows =  1, ncols = 8, pt.size = 0, colors = color_list$treat)
p2 <- vsplot_many(cd8, groups = 'treatment', gene = rna.mem, nrows =  1, ncols = 8, pt.size = 0, colors = color_list$treat)
p3 <- vsplot_many(cd8, groups = 'treatment', gene = rna.tex, nrows =  1, ncols = 8, pt.size = 0, colors = color_list$treat)
p4 <- vsplot_many(cd8, groups = 'treatment', gene = rna.per, nrows =  1, ncols = 8, pt.size = 0, colors = color_list$treat)

p5 <- vsplot_many(cd8, groups = 'treatment', gene = adt.mem, names = c('CD62L', 'CD27', 'CD127', 'CCR7', 'CD45RA'), 
                  assay = 'ADT', feature_type = 'adt', nrows =  1, ncols = 8, pt.size = 0, colors = color_list$treat)

p6 <- vsplot_many(cd8, groups = 'treatment', gene = adt.eff, names = c('CXCR3', 'KLRG1', 'CD57', 'CD45RO'), 
                  assay = 'ADT', feature_type = 'adt', nrows =  1, ncols = 8, pt.size = 0, colors = color_list$treat)

p7 <- vsplot_many(cd8, groups = 'treatment', gene = adt.tex, names = c('PD-1', 'TIM-3', 'LAG-3', 'CD39', 'TIGIT'), 
                  assay = 'ADT', feature_type = 'adt', nrows =  1, ncols = 8, pt.size = 0, colors = color_list$treat)


ggsave('Fig2/Fig2_tmarkers1.pdf', plot = p1, height = 2200, width = 12000, dpi = 850, units = 'px')
ggsave('Fig2/Fig2_tmarkers2.pdf', plot = p2, height = 2200, width = 12000, dpi = 850, units = 'px')
ggsave('Fig2/Fig2_tmarkers3.pdf', plot = p3, height = 2200, width = 12000, dpi = 850, units = 'px')
ggsave('Fig2/Fig2_tmarkers4.pdf', plot = p4, height = 2200, width = 12000, dpi = 850, units = 'px')
ggsave('Fig2/Fig2_tmarkers5.pdf', plot = p5, height = 2200, width = 12000, dpi = 850, units = 'px')
ggsave('Fig2/Fig2_tmarkers6.pdf', plot = p6, height = 2200, width = 12000, dpi = 850, units = 'px')
ggsave('Fig2/Fig2_tmarkers7.pdf', plot = p7, height = 2200, width = 12000, dpi = 850, units = 'px')

###----------------------------------------------GSEA --------------------------------------####

cd8.top_terms <- c("PID_NFAT_3PATHWAY", "GO_GTPASE_BINDING" , "HALLMARK_MTORC1_SIGNALING", "HALLMARK_IL2_STAT5_SIGNALING" , 
                   "GSE9650_EFFECTOR_VS_MEMORY_CD8_TCELL_UP", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "GO_T_CELL_ACTIVATION", "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
                   "GO_CELL_KILLING",  "GO_T_CELL_APOPTOTIC_PROCESS" ,  "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP", "WP_TGFBETA_RECEPTOR_SIGNALING" )

paths.filt <- dplyr::filter(paths, ID %in% cd8.top_terms)
paths.filt$label <- format_string(paths.filt$ID)

paths.filt$label[1] <- "PID: NFAT Pathway"
paths.filt$label[2] <- "GO: GTPase Binding" 
paths.filt$label[3] <- "WP: TGF-beta Receptor Signaling" 
paths.filt$label[4] <- "Hallmark: mTORC1 Signaling"
paths.filt$label[5] <- "Kaech: Naive Vs Day8 Eff CD8 T Cell Up"
paths.filt$label[6] <- "Hallmark: IL-2—STAT5 Signaling" 
paths.filt$label[8] <- 'Hallmark: Interferon Gamma Response'
paths.filt$label[11] <- "GSE9650: Effector Vs Memory CD8 T Cell Up"
paths.filt$label[12] <- 'Hallmark: Oxidative Phosphorylation'

paths.filt$enrichment <- ifelse(paths.filt$NES > 0, 'Cilta', 'Ide')

p <- ggplot(paths.filt, aes(x = reorder(label, NES), y = NES)) +
  geom_segment(aes(xend = label, yend = 0), size = 1) +
  geom_point(aes(size = setSize, fill = enrichment), shape = 21, color = "black") + guides(size = guide_legend(title = 'Size')) +
  scale_fill_manual(values = color_list$treat, guide = 'none') + xlab('Term') +
  theme_minimal() + coord_flip() +  theme_pubr()+ theme(legend.position = "right", axis.title = element_text(face = 'bold')) + 
  scale_size(range = c(3, 6)) 

ggsave('Fig2/Fig2_gsea.pdf', plot = p, height = 3600, width = 6400, dpi = 850, units = 'px')


###----------------------------------------------cytosig GSEA --------------------------------------####
csig$label <- str_replace_all(csig$ID, '_', ' ')
csig$enrichment <- ifelse(csig$NES > 0, 'Cilta', 'Ide')
csig_filt <- dplyr::filter(csig, ID %in% c('IL7_up', 'IL2_up', 'IL1b_up', 'IL1a_up', 'IL15_up', 'IL7_down', 'IL2_down', 'IL15_down', 'IL1b_down', 'IL1a_down'))
csig_filt$label <- str_replace_all(csig_filt$label, 'IL', 'IL-')
csig_filt$label <- str_replace_all(csig_filt$label, 'a', 'A')

p <- ggplot(csig_filt, aes(x = reorder(label, NES), y = NES)) +
  geom_segment(aes(xend = label, yend = 0), size = 1) +
  geom_point(aes(size = setSize, fill = enrichment), shape = 21, color = "black") + guides(size = guide_legend(title = 'Size')) +
  scale_fill_manual(values = color_list$treat, guide = 'none') + xlab('Cytokine') +
  theme_minimal() + coord_flip() +  theme_pubr()+ theme(legend.position = "right", axis.title = element_text(face = 'bold'))+ 
  scale_size(range = c(3, 6)) 

ggsave('Fig2/Fig2_cytosig.pdf', plot = p, height = 3600, width = 4500, dpi = 850, units = 'px')

###----------------------------------------------GSEA heatmap --------------------------------------####

cd8.core_top <- rbind(data.frame(pathway="GTPase Activity", gene = c('ELMO1', 'DOCK2', 'YBX1', 'DIAPH3', 'RAPGEF6', 'KPNB1')), 
                     data.frame(pathway="TGF-beta Signaling", gene = c('RUNX2', 'RUNX3', 'ZEB2', 'TGFBR3', 'TGFBR2', 'STAT1')), 
                     data.frame(pathway="T cell Memory", gene = c('RERE', 'SELL', 'SATB1', 'LEF1', 'CCR7', 'IL7R')), 
                     data.frame(pathway="IL2/STAT5 Signaling", gene = c('BATF3', 'BCL2', 'IL2RA', 'TNFSF10', 'BATF', 'IRF4')), 
                     data.frame(pathway="T cell Activation", gene = c('B2M', 'CCL5', 'CD40LG', 'CD81', 'RAC2')), 
                     data.frame(pathway="T cell Apoptosis", gene = c('LGALS3', 'BAX', 'KDELR1', 'BAK1', 'FAS')), 
                     data.frame(pathway="T cell Migration", gene = c('ITGA4', 'CD99', 'CCL3', 'CXCR3')), 
                     data.frame(pathway="T cell effectorness", gene = c('GZMA', 'GZMB', 'GNLY', 'IFNG')), 
                     data.frame(pathway="Oxidative Phosphorylation", gene = c('NDUFC2', 'NDUFA4', 'COX8A')))

cd8$ta <- paste0(cd8$treatment, '_', cd8$anno1_short)


CairoPDF('Fig2/Fig2_gene_hm.pdf', width = 6.5, height = 9.5)


enhancedHeatmap_with_grouping(cd8, featuredf = cd8.core_top, group.by = 'ta', big.group = 'treatment', scale = T, column_title = ' ', 
                              column_labels = c('Prolif_EM', 'Prolif_TE', 'TE', 'Prolif_EM', 'Prolif_TE', 'TE'), column_names_rot = 45, 
                              big.group_cols = color_list$treat)

dev.off()

###----------------------------------------------cytosig heatmap --------------------------------------####

cd8_csig_core_top <- rbind(data.frame(pathway = 'IL-2/IL-7/IL-15 up', 
                                      gene = c('NCL', 'BCL2', 'FYN', 'STAT1', 'YBX1', 'LARP1', 'RNF213', 'CDK6', 'CCND2')), 
                           data.frame(pathway = 'IL-1A/b up', 
                                      gene = c('SATB1', 'SELL', 'RUNX3', 'NFKBIA', 'STAT3', 'BATF', 'TGFBR2') ), 
                           data.frame(pathway = 'IL-2/IL-7/IL-15 down', 
                                      gene = c('EMP3', 'ITGA4', 'LTB', 'KLF2', 'CCL5')), 
                           data.frame(pathway = 'IL-1A/b down', 
                                      gene = c('NKG7', 'CXCR3', 'LCK', 'MXD4', 'RAC2', 'CTSW')) )

CairoPDF('Fig2/Fig2_cyto_gene_hm.pdf', width = 6, height = 7)


enhancedHeatmap_with_grouping(cd8, featuredf = cd8_csig_core_top, group.by = 'ta', big.group = 'treatment', scale = T, column_title = ' ', 
                              column_labels = c('Prolif_EM', 'Prolif_TE', 'TE', 'Prolif_EM', 'Prolif_TE', 'TE'), column_names_rot = 45, 
                              big.group_cols = color_list$treat)

dev.off()

###----------------------------------------------fitness sig --------------------------------------####

p1 <- vsplot_many(cd8, groups = 'treatment', gene = c('CAR_Fit_Up_1', 'CAR_Migratory_1'), nrows =  1, ncols = 2, pt.size = 0, colors = color_list$treat, 
                  feature_type = 'score', names = c('CAR Fit Score', 'CAR Migratory Score'))

ggsave('Fig2/Fig2_car_fit.pdf', plot = p1, height = 2200, width = 5200, dpi = 850, units = 'px')
