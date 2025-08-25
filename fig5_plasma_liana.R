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
library(circlize)
library(liana)
library(tibble)
library(tidyr)

showtext::showtext_auto()

setwd("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting")
source("/project/jhuangime/tony/Scripts/scRNA_core.R")

color_list <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting/colors.rds")

liana_anno1 <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/BM_Analysis/Data/liana_res_anno1.rds")

liana_anno1 <- liana_anno1 %>% bind_rows(.id = 'treatment')
liana_anno1$lr <- paste0(liana_anno1$ligand.complex, '^' , liana_anno1$receptor.complex)

plasma_interacts <- readRDS('/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/BM_Analysis/Data/plasma_interacts_strict.rds')

shared_types <- c('cDC', 'Mono_Conv', 'Mono_Int', 'Mono_NConv', 'NK_CD16', 'NK_CD56', 'T_CD4', 'T_CD8', 'Plasma')

liana_filt <- dplyr::filter(liana_anno1, (target %in% shared_types) & (source %in% shared_types) & treatment == 'Ide')

###---------------------------------------------heatmaps --------------------------------------####

plasma_sum <- plasma_interacts %>% group_by(source, target) %>% summarize(count=n_distinct(lr), .groups = 'drop')

liana_filt_mtx <- get_freq(liana_filt)
plasma_mtx <- get_freq(plasma_interacts)

rownames(liana_filt_mtx) <- str_replace_all(rownames(liana_filt_mtx), '_', ' ')
colnames(liana_filt_mtx) <- str_replace_all(colnames(liana_filt_mtx), '_', ' ')

CairoPDF('Fig10b/Fig10b_ide_filt_interacts.pdf', width = 5.5, height = 4.5)

draw(ComplexHeatmap::Heatmap(liana_filt_mtx, rect_gp = gpar(col = 'black', lwd = 1), name = "Interactions", column_title_rot=0, 
                             column_title = 'Target', row_title = 'Source', column_names_rot = 45, column_title_side = 'bottom'), 
     padding = unit(c(.6, .5, .2, .2), "cm"))

dev.off()

rownames(plasma_mtx) <- str_replace_all(rownames(plasma_mtx), '_', ' ')
colnames(plasma_mtx) <- str_replace_all(colnames(plasma_mtx), '_', ' ')

CairoPDF('Fig10b/Fig10b_ide_filt_plasma_interacts.pdf', width = 5, height = 4)

draw(ComplexHeatmap::Heatmap(plasma_mtx, rect_gp = gpar(col = 'black', lwd = 1), name = "Interactions", column_title_rot=0, 
                             column_title = 'Target', row_title = 'Source', column_names_rot = 45, column_title_side = 'bottom'), 
     padding = unit(c(.6, .5, .2, .2), "cm"))

dev.off()


###---------------------------------------------mono-plasma interacts plot --------------------------------------####

plasma_mono_ranks <- as.data.frame(dplyr::arrange(dplyr::filter(plasma.source_strict, target == 'Mono_Conv'), sca.LRscore))
plasma_mono_ranks$rank <- as.numeric(rownames(plasma_mono_ranks))

plasma_mono_ranks$lr <- str_replace_all(plasma_mono_ranks$lr, '\\^', '-')
plasma_mono_ranks$lr <- factor(plasma_mono_ranks$lr, levels = plasma_mono_ranks$lr)

p1 <- ggplot(data = plasma_mono_ranks, aes(y = lr, x = sca.LRscore))   + geom_point()  + xlab('LR Score') + ylab('Ligand-Receptor') + ggtitle("Mono-Plasma Specific Interactions")+
  theme_pubr() + theme(plot.title = element_text(hjust = .5), axis.title = element_text(face = 'bold', size = 14))

ggsave('Fig10b/Fig10b_mono_plasma_ranks.pdf', plot = p1, height = 6000, width = 6000, dpi = 850, units = 'px')

###---------------------------------------------enrichment --------------------------------------####

plasma.sourceinteract_sourceannot <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/BM_Analysis/Data/plasma_source_ligand_enrich.rds")

psi_sa_terms <- c('GOMF: Integrin Binding', 'GOMF: Cytokine Activity', 'GOMF: Growth Factor Activity', 'GOBP: Regulation Of Inflammatory Response', 
                  'GOBP: Positive Regulation Of Mapk Cascade')

filt_df <- dplyr::filter(plasma.sourceinteract_sourceannot[!duplicated(plasma.sourceinteract_sourceannot$ID), ], ID %in% psi_sa_terms)
filt_df$label <- filt_df$ID
filt_df$label[5] <- 'GOBP: Positive Regulation Of MAPK Cascade'

p2 <- ggplot(filt_df, aes(x = reorder(label, GeneRatio), y = GeneRatio)) +
  geom_segment(aes(xend = label, yend = 0), size = 1) +
  geom_point(aes(size = Count), shape = 21, color = "black", fill = 'orange') + scale_size_continuous(limits = c(1, 15))+
  theme_minimal() + coord_flip() +
  labs(x = "Pathway", y = "Gene Ratio", title = "Plasma-Mono Ligand Enrichment", size = "Size") + theme_pubr()+
  theme(legend.position = "right", axis.title.x = element_text(face = 'bold'), axis.title.y = element_blank()) 


ggsave('Fig10b/Fig10_plasma_source_plasma_enrich.pdf', plot = p2, height = 2400, width = 5800, dpi = 850, units = 'px')

plasma.source_receptor <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/BM_Analysis/Data/plasma_source_receptor_enrich.rds")

source.terms <- c('GOBP: Myeloid Leukocyte Activation', 'GOBP: Macrophage Activation', 'BROWN: Myeloid Cell Development Up')

filt_df <- dplyr::filter(plasma.source_receptor, ID %in% source.terms)

filt_df <- filt_df[!duplicated(filt_df$Description), ]

filt_df$label <- filt_df$ID
filt_df$label[3] <- 'Brown: Myeloid Cell Development Up'

p3 <- ggplot(filt_df, aes(x = reorder(label, GeneRatio), y = GeneRatio)) +
  geom_segment(aes(xend = label, yend = 0), size = 1) +
  geom_point(aes(size = Count), shape = 21, color = "black", fill = 'orange') + scale_size_continuous(limits = c(1, 15))+
  theme_minimal() + coord_flip() +
  labs(x = "Pathway", y = "Gene Ratio", title = "Plasma-Mono Receptor Enrichment", size = "Size") + theme_pubr()+
  theme(legend.position = "right", axis.title.x = element_text(face = 'bold'), axis.title.y = element_blank()) 

ggsave('Fig10b/Fig10_plasma_source_mono_enrich.pdf', plot = p3, height = 2400, width = 5800, dpi = 850, units = 'px')


