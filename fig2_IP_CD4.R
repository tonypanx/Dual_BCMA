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

markers <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/data/cd4_treatment_markers_nomt.rds")
paths <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/data/cd4_treatment_paths.rds")
csig <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/data/cd4_csig.rds")

###----------------------------------------------Reharmonize labels--------------------------------------####

cd4$anno1 <- str_replace_all(cd4$anno1, 'Tem', 'EM')
cd4$anno1 <- str_replace_all(cd4$anno1, 'TE', 'Th1')
cd4$anno1 <- str_replace_all(cd4$anno1, 'Tcm', 'CM')

cd4$anno2 <- str_replace_all(cd4$anno2, 'Tem', 'EM')
cd4$anno2 <- str_replace_all(cd4$anno2, 'TE', 'Th1')
cd4$anno2 <- str_replace_all(cd4$anno2, 'Tcm', 'CM')

cd4$anno1_short <- str_replace_all(cd4$anno1, 'CD4_', '')
cd4$anno1_short <- factor(cd4$anno1_short, levels = c('Prolif', 'CM', 'EM', 'Th1', 'Treg'))

cd4$treatment <- cd4$treatment_short
cd4$treatment <- revalue(cd4$treatment, c('Cil' = 'Cilta'))

saveRDS(cd8, '/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/RDS/KK_ip_combined_reprocessed_nomt_cd8.rds')
saveRDS(cd4, '/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/RDS/KK_ip_combined_reprocessed_nomt_cd4.rds')

###----------------------------------------------CD4 Dimplot--------------------------------------####

p <- DimPlot_better(cd4,label = F, group.by = 'anno1_short', cols = color_list$cd4_anno1[levels(cd4$anno1_short)]) + theme(plot.title = element_blank())

ggsave('Fig4/Fig4_anno1_UMAP.pdf', plot = p, height = 2400, width = 3300, dpi = 850, units = 'px')

###----------------------------------------------CD4 Dimplot by treatment--------------------------------------####
p1 <- DimPlot_better(subset(cd4, treatment == 'Cilta'), pt.size = .7, label = F, group.by = 'anno1_short', cols = color_list$cd4_anno1[levels(cd4$anno1_short)]) + theme(plot.title = element_blank()) + NoLegend()
p2 <- DimPlot_better(subset(cd4, treatment == 'Ide'), pt.size = .4, label = F, group.by = 'anno1_short', cols = color_list$cd4_anno1[levels(cd4$anno1_short)]) + theme(plot.title = element_blank()) + NoLegend()
p3 <- DimPlot_better(cd4,label = F, group.by = 'treatment', cols = color_list$treat) + theme(plot.title = element_blank())


p <- ggarrange(plotlist = list(p1, p2), ncol = 2, nrow = 1)
ggsave('Fig4/Fig4_anno1_UMAP_splittreat.pdf', plot = p, height = 2400, width = 5000, dpi = 850, units = 'px')

ggsave('Fig4/Fig4_UMAP_treatsplit.pdf', plot = p3, height = 2400, width = 3300, dpi = 850, units = 'px')

###----------------------------------------------CD4 proportions, anno1--------------------------------------####
ftable <- FetchData(cd4, vars = c("anno1_short", "treatment")) 

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
        strip.text.y = element_text(angle = 270, face = "bold", size=10), strip.placement = "outside",panel.grid.major.y = element_blank()) + scale_fill_manual(values = color_list$cd4_anno1[levels(cd4$anno1_short)])

ggsave('Fig4/Fig4_anno1_bar.pdf', plot = p, height = 2400, width = 6000, dpi = 850, units = 'px')

###----------------------------------------------Vln markers, anno1--------------------------------------####



rna.markers <- c('CD4', 'MKI67', 'IL7R', 'SELL', 'CCR7', 'TBX21', 'CXCR3', 'GZMA', 'CCL5','FOXP3', 'IL2RA')
adt.markers <- c('CD127', 'CD62L', 'CCR7.1', 'CD45RA', 'CXCR3.1', 'CD45RO')


p <- VlnPlot(cd4, features = c(rna.markers, adt.markers), group.by = 'anno1_short', stack = T, flip = T) + 
  NoLegend() + theme(strip.text = element_text(face = 'plain'), axis.title.x = element_blank(), axis.title.y = element_text(face = 'bold'))

p$data$feature <- paste0('*', p$data$feature, '*')

p$data$feature <- revalue(p$data$feature, c('*adt_CD127*' = 'CD127',
                                            '*adt_CD62L*' = 'CD62L',
                                            '*adt_CCR7.1*' = 'CCR7',
                                            '*adt_CD45RA*' = 'CD45RA',
                                            '*adt_CD45RO*' = 'CD45RO',
                                            '*adt_CXCR3.1*' = 'CXCR3'))

p$data$feature <- factor(p$data$feature, levels = c(paste0('*', rna.markers, '*'), c('CD127', 'CD62L', 'CCR7','CD45RA', 'CXCR3', 'CD45RO') ))

p <- p + theme(strip.text  =element_markdown())

ggsave('Fig4/Fig4_anno1_markers.pdf', plot = p, height = 5400, width = 3000, dpi = 850, units = 'px')

###----------------------------------------------Vln markers, contrast--------------------------------------####

rna.eff <- c('NKG7', 'GZMA',  'PRF1', 'CXCR3', 'IL17RB' )
rna.mem <- c('IL7R', 'CCR7', 'SELL', 'LEF1')
rna.tex <- c('HAVCR2', 'PDCD1', 'LAG3', 'TOX')
rna.per <- c('BATF3', 'BATF', 'BCL2', 'BAX')

adt.mem <- c( 'CD127', 'CCR7.1', 'CD62L', 'CD45RA' )
adt.eff <- c('CXCR3.1', 'CD57', 'CD45RO')
adt.tex <- c('PD-1', 'TIM-3', 'LAG-3', 'CD39', 'TIGIT.1')

p1 <- vsplot_many(cd4, groups = 'treatment', gene = rna.eff, nrows =  1, ncols = 5, pt.size = 0, colors = color_list$treat)
p2 <- vsplot_many(cd4, groups = 'treatment', gene = rna.mem, nrows =  1, ncols = 5, pt.size = 0, colors = color_list$treat)
p3 <- vsplot_many(cd4, groups = 'treatment', gene = rna.tex, nrows =  1, ncols = 5, pt.size = 0, colors = color_list$treat)
p4 <- vsplot_many(cd4, groups = 'treatment', gene = rna.per, nrows =  1, ncols = 5, pt.size = 0, colors = color_list$treat)

p5 <- vsplot_many(cd4, groups = 'treatment', gene = adt.mem, names = c('CD127', 'CCR7', 'CD62L', 'CD45RA'), 
                  assay = 'ADT', feature_type = 'adt', nrows =  1, ncols = 5, pt.size = 0, colors = color_list$treat)

p6 <- vsplot_many(cd4, groups = 'treatment', gene = adt.eff, names = c('CXCR3',  'CD57', 'CD45RO'), 
                  assay = 'ADT', feature_type = 'adt', nrows =  1, ncols = 5, pt.size = 0, colors = color_list$treat)

p7 <- vsplot_many(cd4, groups = 'treatment', gene = adt.tex, names = c('PD-1', 'TIM-3', 'LAG-3', 'CD39', 'TIGIT'), 
                  assay = 'ADT', feature_type = 'adt', nrows =  1, ncols = 5, pt.size = 0, colors = color_list$treat)


ggsave('Fig4/Fig4_tmarkers1.pdf', plot = p1, height = 2200, width = 7600, dpi = 850, units = 'px')
ggsave('Fig4/Fig4_tmarkers2.pdf', plot = p2, height = 2200, width = 7600, dpi = 850, units = 'px')
ggsave('Fig4/Fig4_tmarkers3.pdf', plot = p3, height = 2200, width = 7600, dpi = 850, units = 'px')
ggsave('Fig4/Fig4_tmarkers4.pdf', plot = p4, height = 2200, width = 7600, dpi = 850, units = 'px')
ggsave('Fig4/Fig4_tmarkers5.pdf', plot = p5, height = 2200, width = 7600, dpi = 850, units = 'px')
ggsave('Fig4/Fig4_tmarkers6.pdf', plot = p6, height = 2200, width = 7600, dpi = 850, units = 'px')
ggsave('Fig4/Fig4_tmarkers7.pdf', plot = p7, height = 2200, width = 7600, dpi = 850, units = 'px')

###----------------------------------------------GSEA --------------------------------------####

cd4.top_terms <- c("GSE11057_NAIVE_VS_EFF_MEMORY_CD4_TCELL_UP", "KEGG_MTOR_SIGNALING_PATHWAY", "WP_TGFBETA_RECEPTOR_SIGNALING", 
                   "GO_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "GO_T_CELL_MEDIATED_CYTOTOXICITY", "GSE11057_NAIVE_VS_EFF_MEMORY_CD4_TCELL_DN", "GO_APOPTOTIC_SIGNALING_PATHWAY", 
                   'GO_CELL_KILLING', "GO_GTPASE_BINDING" )


paths.filt <- dplyr::filter(paths, ID %in% cd4.top_terms)
paths.filt$label <- format_string(paths.filt$ID)

paths.filt$label[1] <- "WP: TGF-beta Receptor Signaling" 
paths.filt$label[2] <-"GSE11057: Naive Vs Eff Memory CD4 T Cell Up"
paths.filt$label[3] <- "KEGG: mTOR Signaling Pathway"  
paths.filt$label[4] <- "GO: GTPase Binding"
paths.filt$label[6] <- "Hallmark: Interferon Alpha Response"
paths.filt$label[8] <- 'Hallmark: Interferon Gamma Response'
paths.filt$label[9] <- "GSE11057: Naive Vs Eff Memory CD4 T Cell Dn"


paths.filt$enrichment <- ifelse(paths.filt$NES > 0, 'Cilta', 'Ide')

p <- ggplot(paths.filt, aes(x = reorder(label, NES), y = NES)) +
  geom_segment(aes(xend = label, yend = 0), size = 1) +
  geom_point(aes(size = setSize, fill = enrichment), shape = 21, color = "black") + guides(size = guide_legend(title = 'Size')) +
  scale_fill_manual(values = color_list$treat, guide = 'none') + xlab('Term') +
  theme_minimal() + coord_flip() +  theme_pubr()+ theme(legend.position = "right", axis.title = element_text(face = 'bold')) + 
  scale_size(range = c(3, 6)) 

ggsave('Fig4/Fig4_gsea.pdf', plot = p, height = 3600, width = 6400, dpi = 850, units = 'px')


###----------------------------------------------cytosig GSEA --------------------------------------####
csig$label <- str_replace_all(csig$ID, '_', ' ')
csig$label <- str_replace_all(csig$label, 'IL', 'IL-')
csig$label <- str_replace_all(csig$label, ':', '')
csig$enrichment <- ifelse(csig$NES > 0, 'Cilta', 'Ide')
csig_filt <- csig[ c('IL7_up', 'IL2_up', 'IL1a_up', 'IL7_down', 'IL1a_down'),] 

p <- ggplot(csig_filt, aes(x = reorder(label, NES), y = NES)) +
  geom_segment(aes(xend = label, yend = 0), size = 1) +
  geom_point(aes(size = setSize, fill = enrichment), shape = 21, color = "black") + guides(size = guide_legend(title = 'Size')) +
  scale_fill_manual(values = color_list$treat, guide = 'none') + xlab('Cytokine') +
  theme_minimal() + coord_flip() +  theme_pubr()+ theme(legend.position = "right", axis.title = element_text(face = 'bold')) +
  scale_size(range = c(3, 6)) 

ggsave('Fig4/Fig4_cytosig.pdf', plot = p, height = 2800, width = 4500, dpi = 850, units = 'px')

###----------------------------------------------GSEA heatmap --------------------------------------####

cd4.core_top <- rbind(
                      data.frame(pathway="GTPase Activity", gene = c('ELMO1', 'DOCK2', 'YBX1', 'DIAPH3', 'RAPGEF6', 'KPNB1')), 
                      data.frame(pathway="TGF-beta Signaling", gene = c('STAT1', 'FOS', 'TGFBR3', 'TGFBR2', 'RUNX2', 'RUNX3')) ,   
                      data.frame(pathway="T cell Memory", gene = c('RERE', 'CCR7', 'SELL', 'IL7R', 'TCF7')),
                      data.frame(pathway="Oxidative Phosphorylation", gene = c('NDUFA12', 'NDUFB9', 'UQCRB', "NDUFV2")),
                      data.frame(pathway="T cell Effectorness", gene = c('ITGA4', 'GZMA', 'GNLY', 'NKG7', 'IL4')), 
                      data.frame(pathway="T cell Apoptosis", gene = c('LGALS3', 'BAX', 'TRADD', 'BAK1', 'FAS')), 
                      data.frame(pathway="Interferon Alpha Response", gene = c('ISG15', 'IFITM3', 'IFITM2', 'ISG20'))
)


cd4$ta <- paste0(cd4$treatment, '_', cd4$anno1_short)


CairoPDF('Fig4/Fig4_gene_hm.pdf', width = 7, height = 8)


enhancedHeatmap_with_grouping(cd4, featuredf = cd4.core_top, group.by = 'ta', big.group = 'treatment', scale = T, column_title = ' ', 
                              column_labels = c('Prolif', 'CM', 'EM', 'Th1', 'Treg', 'Prolif', 'CM', 'EM', 'Th1', 'Treg'), column_names_rot = 45, 
                              big.group_cols = color_list$treat)

dev.off()

###----------------------------------------------cytosig heatmap --------------------------------------####

cd4_csig_core_top <- rbind(data.frame(pathway = 'IL-2/IL-7 up', 
                                      gene = c('NCL', 'BCL2', 'FYN', 'STAT1', 'YBX1', 'LARP1', 'RNF213', 'CDK6', 'CCND2')), 
                           data.frame(pathway = 'IL-1A up', 
                                      gene = c('SATB1', 'SELL', 'RUNX3', 'NFKBIA', 'STAT3', 'BATF', 'TGFBR2') ), 
                           data.frame(pathway = 'IL-7 down', 
                                      gene = c('ITGA4', 'CD28', 'MXD4')), 
                           data.frame(pathway = 'IL-1A down', 
                                      gene = c( 'EMP3', 'LTB', 'KLF2', 'RAC2')) )

CairoPDF('Fig4/Fig4_cyto_gene_hm.pdf', width = 6, height = 7)


enhancedHeatmap_with_grouping(cd4, featuredf = cd4_csig_core_top, group.by = 'ta', big.group = 'treatment', scale = T, column_title = ' ', 
                              column_labels = c('Prolif', 'CM', 'EM', 'Th1', 'Treg', 'Prolif', 'CM', 'EM', 'Th1', 'Treg'), column_names_rot = 45, 
                              big.group_cols = color_list$treat)

dev.off()

###----------------------------------------------fitness sig --------------------------------------####

p1 <- vsplot_many(cd4, groups = 'treatment', gene = c('CAR_Fit_Up_1', 'CAR_Migratory_1'), nrows =  1, ncols = 2, pt.size = 0, colors = color_list$treat, 
                  feature_type = 'score', names = c('CAR Fit Score', 'CAR Migratory Score'))

ggsave('Fig4/Fig4_car_fit.pdf', plot = p1, height = 2200, width = 5200, dpi = 850, units = 'px')
