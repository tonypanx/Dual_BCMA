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
library(EnhancedVolcano)
library(scplotter)
library(patchwork)
library(ggthemes)
library(RColorBrewer)
library(igraph)
require(tidygraph)
require(ggraph)
library(reshape2)
library(rlist)
library(scater)
library(ggrepel)
library(readxl)
library(fgsea)


setwd("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting")
source("/project/jhuangime/tony/Scripts/scRNA_core.R")

showtext::showtext_auto()

set.seed(1234)

cd8 <- readRDS('/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/RDS/KK_ip_combined_reprocessed_nomt_cd8.rds')

color_list <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting/colors.rds")
markers <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/data/cd8_treatment_markers_nomt.rds") 

cd8.reg_markers <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/data/cd8_ip_reg.rds")
cd8.reg_markers_all <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/data/cd8_ip_reg_all.rds")
regs_df <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/data/cd8_ip_reg_df.rds")

cd8$at <- paste0(cd8$anno1_short, '_', cd8$treatment)

regs_tallied <- readRDS("../pyscenic/regs_tallied_redo.rds")[[1]]
regs_imp_df <- readRDS("../pyscenic/regs_tallied_redo.rds")[[2]]




###----------------------------------------------Regulon Volcano--------------------------------------####

regs_trm <- paste0(c('BATF', 'BATF3', 'ATF3', 'ATF4', 'ATF5', 'ATF7', 'JUNB', 'JUND',  'FOSL2', 'FOSB', 'BACH2','IRF1', 'IRF2', 
                     'IRF3', 'IRF4', 'IRF5', 'IRF6', 'IRF8', 'IRF7', 'IRF8', 'IRF9', 'RUNX1', 'RUNX3', 'NFKB1', 'NFKB2', 'KLF2'), '(+)')

regs.highlight <- c('MYB', 'LEF1', 'RUNX1', 'STAT1', 'ETS1', 'NFKB2', 'BACH2','EOMES',  'RELB', 'BATF3', 'TCF7', 'MXD4', 'MAF', 'KLF2', 'TBX21', 
                    'JUND', 'JUNB', 'CREM', 'RORC', 'ETV1', 'BATF', 'ATF3', 'ATF4', 'ATF5', 'FOSL2', 'IRF4', 
                     'NFKB1')

cd8.reg_markers_all$p_val_adj[cd8.reg_markers_all$p_val_adj < 10^-150] <- 10^-150
cd8.reg_markers_all$reg <- rownames(cd8.reg_markers_all)
rownames(cd8.reg_markers_all) <- str_replace_all(rownames(cd8.reg_markers_all), "\\(\\+\\)$", "")

p <- VolcanoPlot(cd8.reg_markers_all, x = 'avg_log2FC', y = 'p_val_adj', label_size = 4, flip_negatives = T, 
            labels = regs.highlight, highlight = regs.highlight, y_cutoff = 10^-10, 
            xlab =  bquote(paste('Average ', log[2], 'FC')), ylab = "-log<sub>10</sub>(adjusted p-value)")  + 
  theme(legend.position = "none", axis.title.y = ggtext::element_markdown()) 


ggsave('Fig3/Fig3_cd8_reg_volcano.pdf', plot = p, height = 5400, width = 5400, dpi = 850, units = 'px')

###----------------------------------------------Top Regulon PLot--------------------------------------####


id_filt <- dplyr::arrange(dplyr::filter(regs_imp_df, Freq >=10), Freq)

id_filt$reglist <- factor(id_filt$reglist, levels = id_filt$reglist)


p <- ggplot(tail(id_filt, 60), aes(x = Freq, y = reglist, fill = Freq)) + geom_bar(stat='identity') + theme_void() + 
  theme(axis.text.y=element_text(size = 8), plot.margin = margin(t=0.5, b=0.5, l=0.5, r=.5, unit = "cm")) 

ggsave('Fig3/Fig3_cd8_reg_freq.pdf', plot = p, height = 6000, width = 4400, dpi = 850, units = 'px')

###----------------------------------------------Regulon Heatmap--------------------------------------####

top_genes <- as.character(tail(id_filt, 60)$reglist)

#enhancedHeatmap_splitcat(cd8, features = c(intersect(top_genes, rownames(cd8.reg_markers)), 'BATF(+)'), group.by = 'anno1', splitcat = 'treatment',  x_horizontal=T, show_row_names=T, assay = 'REG', 
#                         cluster_rows = F, column_names_rot = 45)

CairoPDF('Fig3/Fig3_regulon_hm.pdf', width = 4, height = 7)

hm_regs <- paste0(c('STAT1', 'RUNX3', 'NFKB2', 'MXD4', 'KLF2', 'JUND', 'IKZF1', 'GATA3', 'ETS1', 'MAF', 'BACH2', 'JUNB', 'EOMES', 
              'RUNX1', 'RORC', 'RELB', 'FOXO1', 'BHLHE40', 'NFKB1', 'MYB', 'LEF1', 'TCF7', 'TBX21', 'BATF', 'BATF3', 'ETV1', 'IRF9', 'GATA3'), '(+)')

enhancedHeatmap_splitcat(cd8, features = hm_regs, group.by = 'anno1_short', splitcat = 'treatment',  x_horizontal=F, assay = 'REG', 
                         cluster_rows = T, cluster_columns = F, column_names_rot = 45)

dev.off()

###----------------------------------------------Individual Regulon Plots--------------------------------------####
regs_df_cd8 <- regs_df
regs_df_cd8 <- dplyr::filter(regs_df_cd8, gene %in% rownames(markers))
regs_df_cd8$FC <- markers[regs_df_cd8$gene, ]$avg_log2FC
regs_df_cd8$absFC <- abs(regs_df_cd8$FC)

reglist <- vector(mode = 'list', length = 10)

regs_select <- c('RUNX3', 'RUNX1', 'BATF3', 'LEF1', 'STAT1', 'RORC', 'MXD4',  'ETV1', 'TBX21', 'BHLHE40')

for(i in 1:10){
  reg_sel <- dplyr::filter(regs_df_cd8, regulon %in% paste0(regs_select[i], '(+)')) %>% group_by(regulon) %>% 
    top_n(25, wt=importance) %>% top_n(25, wt = absFC)
  
  reg_sel$regulon <- str_replace_all(reg_sel$regulon, "\\(\\+\\)$", " ")
  
  TFgraph <- graph_from_data_frame(reg_sel, directed = F) 
  tidyTFgraph <- as_tbl_graph(TFgraph) 
  
  regname <- names(tidyTFgraph[1])[1]
  
  FCs <- markers[names(tidyTFgraph[1]), 'avg_log2FC']
  
  tidyTFgraph <- tidyTFgraph %>% activate(nodes) %>% mutate(FC = FCs, labels1 = ifelse(name == regname, regname, ''), labels2 = ifelse(name != regname, name, ''))

  reglist[[i]] <- ggraph(tidyTFgraph, 'stress') +
    geom_node_point(size = 2, aes(color = FC)) + 
    scale_color_gradient2(low = 'blue',  mid = 'white', high = 'red')+
    geom_edge_link(alpha=0.1,label_colour = "grey")+
    geom_node_text(aes(label=labels1), fontface = 'bold', repel = F,
                   size = 6, vjust = 1.5)+
    geom_node_text(aes(label=labels2), fontface = 'italic', repel = T,
                   size = 4)+
    theme_pubr()+
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.line=element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size=20, hjust = 0.5, face = 'bold'),
      legend.text = element_text(size=15),
      legend.title = element_text(size = 15, face = 'bold'),
      legend.position = 'right')+
    guides(shape = guide_legend(override.aes = list(size = 15)))
  
}

p <- ggarrange(plotlist=reglist, nrow = 2, ncol = 5)

ggsave('Fig3/Fig3_cd8_reg_sel.pdf', plot = p, height = 9000*.8, width = 27000*.8, dpi = 850, units = 'px')


###----------------------------------------------TRM signature vlnplt--------------------------------------####
  

p1 <- vsplot(cd8, groups = 'treatment', gene = c('car_trm_up_1'), pt.size = 0, colors = color_list$treat, 
                  feature_type = 'score') + ggtitle('CAR TRM Score')

ggsave('Fig3/Fig3_car_trm_vln.pdf', plot = p1, height = 2500, width = 2700, dpi = 850, units = 'px')

###----------------------------------------------TRM signature enrichment--------------------------------------####

cd8.ranks <- markers$avg_log2FC
names(cd8.ranks) <- rownames(markers)

trm_car_data <- read_excel('/project/jhuangime/tony/projects/BCMA_CAR/external_data/trm_car_1.xlsx', col_names = T, skip = 0)
trm_car_data <- trm_car_data[trm_car_data$`P Value` < .05, ]

car_trm_sigs_filt2 <- list(trm_car_up =  trm_car_data[trm_car_data$`log FC` > .5, ]$`Gene ID`)

p1 <- plotEnrichment(car_trm_sigs_filt2$trm_car_up, cd8.ranks[abs(cd8.ranks) > .33]) + ggtitle('CD8 CAR-TRM vs CAR-TConv UP') +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text = element_text(size = 12), 
        plot.title = element_text(hjust = .5, face = 'bold', size = 14)) + ylab('Running Enrichment Score') + xlab('Rank') +
  annotate('text', x= 50, y = .055, label = 'NES=1.57', hjust = 0)+
  annotate('text', x= 50, y = .040, label = 'Adjusted p-value=0.024', hjust = 0)

ggsave('Fig3/Fig3_car_trm_gsea.pdf', plot = p1, height = 4000, width = 3900, dpi = 850, units = 'px')

###----------------------------------------------TRM featureplots --------------------------------------####

p <- FeaturePlot_better(cd8, features = c('RUNX3', 'KLF2'), max.cutoff = 'q90', order = T, pt.size = .5) 


ggsave('Fig3/Fig3_trm_gene_UMAP.pdf', plot = p, height = 2800, width = 4400, dpi = 850, units = 'px')

p <- FeaturePlot_better(cd8, features = c('RUNX3(+)', 'KLF2(+)'), max.cutoff = 'q90', order = F, pt.size = .5, names = c("RUNX3", 'KLF2')) &
  theme(plot.title = element_text(face = 'bold')) & scale_color_viridis(
    breaks = pretty_breaks(n = 3),
    labels = function(x) signif(x, 1) 
  )

ggsave('Fig3/Fig3_trm_reg_UMAP.pdf', plot = p, height = 2800, width = 4400, dpi = 850, units = 'px')

###----------------------------------------------TRM genes --------------------------------------####

trm_features <- c('ATF2', 'ATF3', 'ATF4', 'ATF6', 'BATF', 'BATF3', 'RUNX1', 'RUNX3', 'NFKB1', 'NFKB2', 'KLF2')


p1 <- vsplot_many(cd8, groups = 'treatment', gene = trm_features, nrows =  2, ncols = 6, pt.size = 0, colors = color_list$treat)
ggsave('Fig3/Fig3_trm_markers.pdf', plot = p1, height = 4200, width = 10000, dpi = 850, units = 'px')
