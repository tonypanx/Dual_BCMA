library(Seurat)
library(plyr)
library(dplyr)
library(ggtext)

setwd("/project/jhuangime/tony/projects/BCMA_CAR/Pre_Infusion/Analysis")
source("/project/jhuangime/tony/Scripts/scRNA_core.R")

showtext::showtext_auto()

theme_arial_all <- function(){a
  theme(text=element_text(family = "Arial"))
}

color_list <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting/colors.rds")

obj <- readRDS('RDS/obj_integrated_QC.rds')

obj <- FindClusters(obj, resolution = .3)

markers <- clean_markers(FindAllMarkers(obj, min.pct = .2))

markers %>% group_by(cluster) %>% top_n(20, avg_log2FC) %>% clean_markers() %>% print(n = 400)

obj$cell_type <- dplyr::recode(obj$seurat_clusters,
                               "0"  = "Mono Conv",
                               "1"  = "NK CD56",
                               "2"  = "T CD8",
                               "3"  = "NK CD16",
                               "4"  = "T CD4",
                               "5"  = "Mono NConv",
                               "6"  = "Mono Int",
                               "7"  = "B",
                               "8"  = "Cycling",
                               "9"  = "cDC",
                               "10" = "Megakaryocyte",
                               "11" = "Mono NConv",       
                               "12" = "T Innate-like",
                               "13" = "Mono Int",       
                               "14" = "Plasma",
                               "15" = "B",
                               "16" = "Neutrophil"
)

Idents(obj) <- obj$cell_type

markers <- clean_markers(FindAllMarkers(obj, min.pct = .2))

markers %>% group_by(cluster) %>% top_n(30, avg_log2FC) %>% clean_markers() %>% print(n = 420)

obj$cell_type <- as.character(obj$cell_type)

Idents(obj) <- factor(obj$cell_type, levels = sort(unique(obj$cell_type), decreasing = F))

marker.genes <-c('MS4A1', 'CD1C', 'MKI67', 'PPBP', 'CD14', 'FCGR3A', 'FCGR3B', 'NCAM1', 'MZB1', 'TNFRSF17', 'CD3G', 'CD4', 'CD8A', 'TRGV4')

cell_type_colors <- c(
  "B"              = "#E8836A",
  "cDC"            = "#E07B2A",
  "Mono Conv"      = "#6B7B2A",
  "Mono Int"       = "#4A9B3A",
  "Mono NConv"     = "#2A9B6A",
  "NK CD16"        = "#2A9B8B",
  "NK CD56"        = "#2ABBD4",
  "T Innate-like"  = "#1A6B9B",
  "Neutrophil"     = "#B8860B",
  "Megakaryocyte"  = "#D4A0C0",  
  "Plasma"         = "#C49BE0",
  "T CD4"          = "#E040A0",
  "T CD8"          = "#C0306A",
  "Cycling"        = "#A0A0A0"
)

p <- DimPlot_better(obj, cols = cell_type_colors, group.by = "cell_type") + theme(plot.title = element_blank()) 

ggsave('Plots/preinfusion_all_umap.pdf',  plot = p, width = 5000, height = 3200, dpi = 900, units = 'px')

p <- VlnPlot_better(obj, genes = marker.genes, group.by = "cell_type")

ggsave('Plots/preinfusion_all_vln.pdf',  plot = p, width = 5400, height = 5400, dpi = 900, units = 'px')

obj$Treatment_Time <- paste0(obj$Treatment, '_', obj$Time_Group)
### plot proportions ###

sample_colors <- c(
  "Cilta_Pre" = "red",  
  "Ide_Pre"    = "#87CEEB"  # skyblue
)


meta <- obj@meta.data %>%
  select(cell_type, Treatment_Time) %>%          
  mutate(Treatment_Time = factor(Treatment_Time, levels = c("Ide_Pre", "Cilta_Pre")))

p1 <- meta %>%
  group_by(Treatment_Time, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cell_type) %>%
  mutate(pct = n / sum(n) ) %>%
  ggplot(aes(x = Treatment_Time, y = pct, fill = Treatment_Time)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ cell_type, nrow = 2) +
  scale_fill_manual(values = sample_colors) +
  labs(x = NULL, y = "Proportion of Sample") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
    strip.text   = element_text(size = 12),
    legend.position = "none",
    plot.title   = element_text(size = 11)
  )

ggsave('Plots/preinfusion_all_props_by_celltype.pdf',  plot = p1, width = 8500, height = 5000, dpi = 900, units = 'px')

p2 <- meta %>%
  group_by(Treatment_Time, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Treatment_Time) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = Treatment_Time, y = pct, fill = cell_type)) +
  geom_bar(stat = "identity", position = 'dodge', color = 'black') +
  scale_fill_manual(values = cell_type_colors) +
  labs(x = NULL, y = "Proportion of Sample", fill = "Cell type") +
  theme_classic() + facet_wrap(~cell_type, scales = 'free', nrow =2)+
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12), 
    legend.text = element_text(size = 12),
    axis.title.y = element_text(size = 12), 
    strip.text = element_text(size = 12),
    legend.position = "none"
  )

ggsave('Plots/preinfusion_all_props_by_sample.pdf',  plot = p2, width = 9800, height = 5000, dpi = 900, units = 'px')

saveRDS(obj, 'RDS/obj_integrated_QC.rds')
