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
library(numbat)

showtext::showtext_auto()

setwd("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting")
source("/project/jhuangime/tony/Scripts/scRNA_core.R")

color_list <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting/colors.rds")

###---------------------------------------------plasma --------------------------------------####

plasma <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/BM_Analysis/RDS/plasma_new.rds")

mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")

###---------------------------------------------Dimplot--------------------------------------####
plasma$tumor_meta <- revalue(plasma$tumor_meta, c('tumor' = 'Tumor', 'normal' = 'Normal'))

p1 <- DimPlot_better(plasma, group.by = 'tumor_meta' ,pt.size = .5, cols = c('Blue', 'Red')) + theme(plot.title = element_blank())

ggsave('Fig10/Fig10_plasma_tumor_UMAP.pdf', plot = p1, height = 2300, width = 3100, dpi = 800, units = 'px')

p2 <- DimPlot_better(subset(plasma, treatment == 'Ide'), group.by = 'tumor_meta',  pt.size = .5, cols = c('Blue', 'Red')) + theme(plot.title = element_blank()) + 
  ggtitle('Ide') + theme(plot.title = element_text(hjust = .5)) + NoLegend()

p3 <- DimPlot_better(subset(plasma, treatment == 'Cilta'), group.by = 'tumor_meta',  pt.size = .5, cols = c('Blue', 'Red')) + theme(plot.title = element_blank()) + 
  ggtitle('Cilta') + theme(plot.title = element_text(hjust = .5)) + NoLegend() 

p <- ggarrange(plotlist = list(p2, p3), ncol = 2, nrow = 1)

ggsave('Fig10/Fig10_plasma_split_UMAP.pdf', plot = p, height = 2300, width = 4300, dpi = 800, units = 'px')

###---------------------------------------------Featureplot--------------------------------------####

plasma.features <- c('MS4A1', 'CD19', 'TNFRSF17', 'SDC1', 'GPRC5D', 'CD38')

p3 <- FeaturePlot_better(plasma, features = plasma.features, max.cutoff= 'q90', order = T, ncol = 3, pt.size = .5)& scale_color_viridis(
  breaks = pretty_breaks(n = 4),
  labels = function(x) signif(x, 1) 
)

ggsave('Fig10/Fig10_plasma_features_UMAP.pdf', plot = p3, height = 4600, width = 5600, dpi = 800, units = 'px')

p4 <- FeaturePlot_better(plasma, features = 'tumor_cnv_post', max.cutoff= 'q90', order = F, pt.size = .5, face = 'bold') + ggtitle('Aneuploidy Probability') & scale_color_viridis(
  breaks = pretty_breaks(n = 3),
  labels = function(x) signif(x, 1) 
)

ggsave('Fig10/Fig10_aneuploidy_prob_UMAP.pdf', plot = p4, height = 2500, width = 2200, dpi = 800, units = 'px') 

p5 <- vsplot(plasma, gene = 'tumor_cnv_post', groups = 'treatment')+ ggtitle('Aneuploidy Probability') + 
  theme(plot.title = element_text(face = 'plain'), axis.text.x = element_text(size = 16)) + NoLegend()

ggsave('Fig10/Fig10_plasma_prob_comp.pdf', plot = p5, height = 2400, width = 2600, dpi = 850, units = 'px')

###---------------------------------------------anno tumor bar--------------------------------------####

ftable <- FetchData(plasma, vars = c("tumor_meta", "treatment")) 

dfx <- ftable %>% 
  group_by(treatment) %>%
  dplyr::count(tumor_meta) %>%
  mutate(freq = n/sum(n)) %>%
  ungroup()

p <- ggplot(dfx, aes(x= fct_rev(treatment), y=freq, fill= tumor_meta ))+
  geom_col(width = 0.85, color = 'black', position=position_fill(reverse=T))+ theme_pubr()+ scale_y_continuous(expand=c(0,0.01)) +
  guides(fill = guide_legend(override.aes = list(size=4))) + coord_flip() +
  theme_bw() + labs(y = "Proportion") + theme( axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm"), size=15, face = 'bold'),
                                               axis.title.y = element_blank(), legend.text = element_text(size = 13), axis.text.x = element_text(size = 12, color = 'black'),
                                               axis.text.y = element_text(size = 16, margin = margin(r=.1, unit = 'cm'), color = 'black'), legend.title = element_blank(), 
                                               panel.background = element_rect(color = 'black')) +
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, r=.5, unit = "cm"), legend.box.margin = margin(10, 10, 10, 10),
        plot.title = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(angle = 270, face = "bold", size=10), strip.placement = "outside",panel.grid.major.y = element_blank()) + 
  scale_fill_manual(values = c("Blue", 'Red'))


ggsave('Fig10/Fig10_annotumor_bar_bytreat.pdf', plot = p, height = 2400, width = 6000, dpi = 850, units = 'px')

###---------------------------------------------cnv map--------------------------------------####

nb = Numbat$new(out_dir = "/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/BM_Analysis/Numbat/out")

mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")

nb$cutree(n_cut = 3)
p <- nb$plot_phylo_heatmap() 

ggsave('Fig10/Fig10_sc_clones.pdf', plot = p, height = 4000, width = 6000, dpi = 850, units = 'px')

p2 <- nb$bulk_clones %>%
  plot_bulks(
    min_LLR = 10, # filtering CNVs by evidence
    legend = TRUE, 
  )

ggsave('Fig10/Fig10_bulk_clones.pdf', plot = p2, height = 6000, width = 6000, dpi = 850, units = 'px')

p1 <- DimPlot_better(plasma, group.by = 'tumor_clone' ,pt.size = .5, cols = mypal) + theme(plot.title = element_blank())

ggsave('Fig10/Fig10_plasma_clone_UMAP.pdf', plot = p1, height = 2300, width = 2800, dpi = 800, units = 'px')
