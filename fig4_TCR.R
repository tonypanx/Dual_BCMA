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
library(scRepertoire)

setwd("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting")
source("/project/jhuangime/tony/Scripts/scRNA_core.R")

font_import()

showtext::showtext_auto()

set.seed(1234)

cd8 <- readRDS('/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/RDS/KK_ip_combined_reprocessed_nomt_cd8.rds')
cd4 <- readRDS('/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/RDS/KK_ip_combined_reprocessed_nomt_cd4.rds')

color_list <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Plotting/colors.rds")
vdj <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/vdj_all/filt_vdj_obj.rds")
pheno_clono_table <- readRDS("/project/jhuangime/tony/projects/BCMA_CAR/Cilta_seq/Analysis/vdj_all/RDS_fixedbc_withvdj/pheno_clono_table.rds")

###----------------------------------------------Clone size distribution for CD8 CD4--------------------------------------####

clono_dfx_cd8 <- FetchData(cd8, c('treatment', 'cloneSize', 'anno1', 'clonalFrequency', 'CTstrict'))
clono_dfx_cd8$ttype = 'CD8'
clono_dfx_cd8 <- clono_dfx_cd8[complete.cases(clono_dfx_cd8$cloneSize), ]

clono_dfx_cd4 <- FetchData(cd4, c('treatment', 'cloneSize', 'anno1', 'clonalFrequency', 'CTstrict'))
clono_dfx_cd4$ttype = 'CD4'
clono_dfx_cd4 <- clono_dfx_cd4[complete.cases(clono_dfx_cd4$cloneSize), ]

clono_dfx <- rbind(clono_dfx_cd8, clono_dfx_cd4)
clono_dfx$cloneSize <- as.character(clono_dfx$cloneSize)

clono_dfx$cloneSize <- revalue(clono_dfx$cloneSize, c('Large (0.01 < X <= 0.1)' = 'Large (0.01 < X <= 1)' ,
                                                      'Hyperexpanded (0.1 < X <= 1)' = 'Large (0.01 < X <= 1)'))

dfx_cd8 <- dplyr::filter(clono_dfx, ttype == 'CD8') %>% 
  group_by(treatment, cloneSize) %>% summarize(count = n()) %>% 
  mutate(prop = count/sum(count))

dfx_cd4 <- dplyr::filter(clono_dfx, ttype == 'CD4') %>% 
  group_by(treatment, cloneSize) %>% summarize(count = n()) %>% 
  mutate(prop = count/sum(count))

p1 <- ggplot(dfx_cd8, aes(x= treatment, y=prop, fill= cloneSize ))+
  geom_col(width = 0.85, color = 'black', position=position_fill(reverse=F))+ theme_pubr()+ scale_y_continuous(expand=c(0,0.01)) +
  guides(fill = guide_legend(override.aes = list(size=4))) + 
  theme_bw() + labs(y = "Proportion") + theme( axis.title.y = element_text(margin = margin(l = 0.5, r = 0.5, unit = "cm"), size=15, face = 'bold'),
                                               axis.title.x = element_blank(), legend.text = element_text(size = 13), axis.text.x = element_text(size = 12, color = 'black'),
                                               axis.text.y = element_text(size = 16, margin = margin(r = .1, unit = 'cm'), color = 'black'), legend.title = element_blank(), 
                                               panel.background = element_rect(color = 'black')) + ggtitle('CD8')+
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, r=.5, unit = "cm"), legend.box.margin = margin(10, 10, 10, 10),
        plot.title = element_text(size = 15, face = "bold", hjust = .5),
        strip.text.y = element_text(angle = 270, face = "bold", size=10), strip.placement = "outside",panel.grid.major.y = element_blank()) + 
  scale_fill_manual(values = c('red', 'yellow', 'blue'))

p2 <- ggplot(dfx_cd4, aes(x= treatment, y=prop, fill= cloneSize ))+
  geom_col(width = 0.85, color = 'black', position=position_fill(reverse=F))+ theme_pubr()+ scale_y_continuous(expand=c(0,0.01)) +
  guides(fill = guide_legend(override.aes = list(size=4))) + 
  theme_bw() + labs(y = "Proportion") + theme( axis.title.y = element_text(margin = margin(l = 0.5, r = 0.5, unit = "cm"), size=15, face = 'bold'),
                                               axis.title.x = element_blank(), legend.text = element_text(size = 13), axis.text.x = element_text(size = 12, color = 'black'),
                                               axis.text.y = element_text(size = 16, margin = margin(r = .1, unit = 'cm'), color = 'black'), legend.title = element_blank(), 
                                               panel.background = element_rect(color = 'black')) + ggtitle('CD4')+
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, r=.5, unit = "cm"), legend.box.margin = margin(10, 10, 10, 10),
        plot.title = element_text(size = 15, face = "bold", hjust = .5),
        strip.text.y = element_text(angle = 270, face = "bold", size=10), strip.placement = "outside",panel.grid.major.y = element_blank()) + 
  scale_fill_manual(values = c('red', 'yellow', 'blue'))


p <- ggarrange(plotlist = list(p1, p2), nrow = 1, ncol = 2, common.legend = T, legend = 'right')

ggsave('Fig5a/Fig5a_ttype_clonesize_bar.pdf', plot = p, height = 3600, width = 5600, dpi = 500, units = 'px')

###----------------------------------------------Ttype and anno1 distribution of large clones--------------------------------------####

dfx_large <- dplyr::filter(clono_dfx, cloneSize == 'Large (0.01 < X <= 1)')
dfx_large_ttype <- dfx_large %>% 
  group_by(treatment, ttype) %>% summarize(count = n()) %>% 
  mutate(prop = count/sum(count))

dfx_large_anno1 <- dfx_large %>% 
  group_by(treatment, anno1) %>% summarize(count = n()) %>% 
  mutate(prop = count/sum(count))

anno1_order <- c('CD8_Prolif_EM', 'CD8_Prolif_TE', 'CD8_TE', 'CD4_Prolif', 'CD4_EM', 'CD4_Th1')

p1 <- ggplot(dfx_large_ttype, aes(x= treatment, y=prop, fill= ttype ))+
  geom_col(width = 0.85, color = 'black', position=position_fill(reverse=F))+ theme_pubr()+ scale_y_continuous(expand=c(0,0.01)) +
  guides(fill = guide_legend(override.aes = list(size=4))) + 
  theme_bw() + labs(y = "Proportion") + theme( axis.title.y = element_text(margin = margin(l = 0.5, r = 0.5, unit = "cm"), size=15, face = 'bold'),
                                               axis.title.x = element_blank(), legend.text = element_text(size = 13), axis.text.x = element_text(size = 12, color = 'black'),
                                               axis.text.y = element_text(size = 16, margin = margin(r = .1, unit = 'cm'), color = 'black'), legend.title = element_blank(), 
                                               panel.background = element_rect(color = 'black')) + ggtitle('Large Clones')+
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, r=.5, unit = "cm"), legend.box.margin = margin(10, 10, 10, 10),
        plot.title = element_text(size = 15, face = "bold", hjust = .5),
        strip.text.y = element_text(angle = 270, face = "bold", size=10), strip.placement = "outside",panel.grid.major.y = element_blank()) + 
  scale_fill_manual(values = c( "cornsilk3", 'deepskyblue3'))

p2 <- ggplot(dfx_large_anno1, aes(x= treatment, y=prop, fill= anno1 ))+
  geom_col(width = 0.85, color = 'black', position=position_fill(reverse=F))+ theme_pubr()+ scale_y_continuous(expand=c(0,0.01)) +
  guides(fill = guide_legend(override.aes = list(size=4))) + 
  theme_bw() + labs(y = "Proportion") + theme( axis.title.y = element_text(margin = margin(l = 0.5, r = 0.5, unit = "cm"), size=15, face = 'bold'),
                                               axis.title.x = element_blank(), legend.text = element_text(size = 13), axis.text.x = element_text(size = 12, color = 'black'),
                                               axis.text.y = element_text(size = 16, margin = margin(r = .1, unit = 'cm'), color = 'black'), legend.title = element_blank(), 
                                               panel.background = element_rect(color = 'black')) + ggtitle('Large Clones')+
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, r=.5, unit = "cm"), legend.box.margin = margin(10, 10, 10, 10),
        plot.title = element_text(size = 15, face = "bold", hjust = .5),
        strip.text.y = element_text(angle = 270, face = "bold", size=10), strip.placement = "outside",panel.grid.major.y = element_blank()) + 
  scale_fill_manual(values =  color_list$anno1[anno1_order])


ggsave('Fig5a/Fig5a_large_ttype.pdf', plot = p1, height = 3600, width = 2400, dpi = 500, units = 'px')
ggsave('Fig5a/Fig5a_large_anno1.pdf', plot = p2, height = 3600, width = 2800, dpi = 500, units = 'px')

###----------------------------------------------Tracing of top 10 CAR clones at each TP, Ide--------------------------------------####
p1 <- clonalCompare(vdj[c(2, 3, 4)], cloneCall = 'strict', relabel.clones = T, top.clones = 10,  group.by = 'sample', 
              order.by = c('Ide_IP', 'Ide_PB_D9_CAR', 'Ide_PB_D28_CAR')) + theme_pubr() + rotate_x_text(angle = 45) + 
  theme(axis.text.x = element_text(size = 14), legend.position = 'right', axis.title.y = element_text(face = 'bold', size = 15), axis.title.x = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = .5)) + scale_fill_manual(values = distinct_colors(32)) + ggtitle('Ide') +
  scale_x_discrete(labels = c('IP', 'D9', 'D28')) + labs(fill = 'Clone')

p1$data$clones <- str_replace_all(p1$data$clones, 'Clone: ', 'I')
p1$data$clones <- factor(p1$data$clones, levels = paste0('I', 1:25))

ggsave('Fig5a/Fig5a_ide_top10_eachTP.pdf', plot = p1, height = 2800, width = 2800, dpi = 500, units = 'px')

###----------------------------------------------Tracing of top 10 CAR clones at each TP, Cilta--------------------------------------####

p2 <- clonalCompare(vdj[c(11:14, 18)], cloneCall = 'strict', relabel.clones = T, top.clones = 10,  group.by = 'sample', 
              order.by = c('Cilta_IP', 'Cilta_PB_D14_CAR', 'Cilta_PB_D21_CAR', 'Cilta_PB_D28_CAR', "Cilta_BM_D28_CAR")) + theme_pubr() + rotate_x_text(angle = 45) + 
  theme(axis.text.x = element_text(size = 14), legend.position = 'right', axis.title.y = element_text(face = 'bold', size = 15), axis.title.x = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = .5)) + scale_fill_manual(values = distinct_colors(14)) + ggtitle('Cilta') +
  scale_x_discrete(labels = c('IP', 'D14', 'D21', 'D28', 'D28 BM')) + labs(fill = 'Clone')

p2$data$clones <- str_replace_all(p2$data$clones, 'Clone: ', 'C')
p2$data$clones <- factor(p2$data$clones, levels = paste0('C', 1:14))

ggsave('Fig5a/Fig5a_cilta_top10_eachTP.pdf', plot = p2, height = 2800, width = 2800, dpi = 500, units = 'px')

###----------------------------------------------Overlap of Cilta and Ide CAR clones, by TP--------------------------------------####

p3 <- clonalOverlap(vdj[c( 11, 12, 13, 14, 18, 2, 3, 4)], cloneCall = 'strict', method = 'morisita')+ theme_pubr() + rotate_x_text(angle = 45) + 
  theme(axis.text = element_text(size = 12), legend.position = 'right') + 
  scale_x_discrete(labels = c( 'Cilta IP', 'Cilta D14', 'Cilta D21', 'Cilta D28', 'Cilta D28 BM', 'Ide IP', 'Ide D9', 'Ide D28')) +
  scale_y_discrete(labels = c( 'Cilta IP', 'Cilta D14', 'Cilta D21', 'Cilta D28', 'Cilta D28 BM', 'Ide IP', 'Ide D9', 'Ide D28')) +
  theme(axis.title = element_blank())

ggsave('Fig5a/Fig5a_overlap.pdf', plot = p3, height = 3000, width = 3400, dpi = 500, units = 'px')

###----------------------------------------------Overlap of Cilta and Ide CAR clones, by treat/type--------------------------------------####
replab <- function(x){
  x$tt <- paste0(x$treatment, '_', x$type)
  return(x)
  }


vdj <- lapply(vdj, replab)

p4 <- clonalOverlap(vdj[c(2:9, 11:19)], cloneCall = 'strict', method = 'morisita', group.by = 'tt')+ theme_pubr() + rotate_x_text(angle = 45) + 
  theme(axis.text = element_text(size = 14), legend.position = 'right') +
  scale_x_discrete(labels = c( 'Cilta IP', expression("Cilta Post CAR"^" +"), expression("Cilta Post CAR"^" -"), 'Ide IP', expression("Ide Post CAR"^" +"), expression("Ide Post CAR"^" -"))) +
  scale_y_discrete(labels = c( 'Cilta IP', expression("Cilta Post CAR"^"+"), expression("Cilta Post CAR"^"-"), 'Ide IP', expression("Ide Post CAR"^"+"), expression("Ide Post CAR"^"-"))) +
  theme(axis.title = element_blank())
  
ggsave('Fig5a/Fig5a_overlap_bygroup.pdf', plot = p4, height = 3000, width = 3400, dpi = 500, units = 'px')

###----------------------------------------------Diversity Comparison--------------------------------------####
cd <- clonalDiversity(vdj[c(2, 3, 4, 11:14, 18)], cloneCall = 'strict', metrics = 'shannon', exportTable = T, group.by = 'tt', skip.boots = T )
cd$tt <- factor(cd$tt, levels = c('Cilta_IP','Ide_IP', 'Cilta_Post_CAR', 'Ide_Post_CAR'))
cd$treatment <- c('Cilta', 'Cilta', 'Ide', 'Ide')
cd$treatment <- factor(cd$treatment, levels = c('Cilta', 'Ide'))
cd$type <- c("IP", 'Post_CAR', 'IP', 'Post_CAR')


p5 <- ggplot(cd, aes(x=type, y = shannon, fill = treatment))+ geom_bar(stat = 'identity', position = 'dodge', color = 'black') + theme_pubr() + ylab('Shannon Diversity')+
  xlab('Sample Type') + theme(legend.position = 'right', legend.text = element_text(size = 12), 
                              axis.title.x = element_blank() , axis.title.y = element_text(face = 'bold', margin = margin(l = 0.3, r = 0.3, unit = "cm"), size=13)) + scale_x_discrete(labels = c('IP', 'Post CAR')) + 
  scale_fill_discrete(name = 'Treatment')+ scale_fill_manual(values = color_list$treat) + labs(fill = 'Treatment')


ggsave('Fig5a/Fig5a_diversity_bygroup.pdf', plot = p5, height = 2400, width = 2000, dpi = 500, units = 'px')

###----------------------------------------------Ide IP large clone trace--------------------------------------####

### create dummy Ide obj for combineExpression
ide_count <- matrix(rep(1, 2*nrow(vdj$Ide_IP)), nrow = 2)
colnames(ide_count) <- vdj$Ide_IP$barcode
rownames(ide_count) <- c('g1', 'g2')

ide_temp <- CreateSeuratObject(counts = ide_count)
ide_temp <- combineExpression(vdj, ide_temp, cloneCall = 'strict', chain = 'both')

ide_cc <- FetchData(ide_temp, c('CTstrict', 'cloneSize'))
ide_cc_large <- dplyr::filter(ide_cc, cloneSize %in% c("Large (0.01 < X <= 0.1)", 'Hyperexpanded (0.1 < X <= 1)'))
ide_cc_large <- ide_cc_large[!duplicated(ide_cc_large$CTstrict), ]

ide_large_trace <- lapply(vdj[c( 2, 3, 4)], function(x){return(nrow(dplyr::filter(x, CTstrict %in% ide_cc_large$CTstrict))/ nrow(x)) })
ide_large_trace <- as.data.frame(ide_large_trace)
ide_large_trace <- data.frame(Sample = names(ide_large_trace), Prop = as.numeric(ide_large_trace))

ide_large_trace$Sample <- factor(ide_large_trace$Sample, levels = ide_large_trace$Sample)

p6 <- ggplot(ide_large_trace, aes(x=Sample, y = Prop))+ geom_bar(stat = 'identity', color = 'black', fill = color_list$treat['Ide']) + theme_pubr() + ylab('Repertoire Occupancy')+
  theme(legend.position = 'right', legend.text = element_text(size = 12), plot.title = element_text(hjust = .5, face = 'bold'), 
                              axis.title.x = element_blank() , axis.title.y = element_text(face = 'bold', margin = margin(l = 0.3, r = 0.3, unit = "cm"), size=13)) + scale_x_discrete(labels = c('IP', 'D9', 'D28')) + 
  scale_fill_manual(values = color_list$treat) + ggtitle('Ide IP Large Clones')

###----------------------------------------------Cilta IP large clone trace--------------------------------------####

### create dummy cilta obj for combineExpression
cilta_count <- matrix(rep(1, 2*nrow(vdj$Cilta_IP)), nrow = 2)
colnames(cilta_count) <- vdj$Cilta_IP$barcode
rownames(cilta_count) <- c('g1', 'g2')

cilta_temp <- CreateSeuratObject(counts = cilta_count)
cilta_temp <- combineExpression(vdj, cilta_temp, cloneCall = 'strict', chain = 'both')

cilta_cc <- FetchData(cilta_temp, c('CTstrict', 'cloneSize'))
cilta_cc_large <- dplyr::filter(cilta_cc, cloneSize %in% c("Large (0.01 < X <= 0.1)", 'Hyperexpanded (0.1 < X <= 1)'))
cilta_cc_large <- cilta_cc_large[!duplicated(cilta_cc_large$CTstrict), ]

cilta_large_trace <- lapply(vdj[c( 11:14, 18)], function(x){return(nrow(dplyr::filter(x, CTstrict %in% cilta_cc_large$CTstrict))/ nrow(x)) })
cilta_large_trace <- as.data.frame(cilta_large_trace)
cilta_large_trace <- data.frame(Sample = names(cilta_large_trace), Prop = as.numeric(cilta_large_trace))

cilta_large_trace$Sample <- factor(cilta_large_trace$Sample, levels = cilta_large_trace$Sample)

p7 <- ggplot(cilta_large_trace, aes(x=Sample, y = Prop))+ geom_bar(stat = 'identity', color = 'black', fill = color_list$treat['Cilta']) + theme_pubr() + ylab('Repertoire Occupancy')+
  theme(legend.position = 'right', legend.text = element_text(size = 12), plot.title = element_text(hjust = .5, face = 'bold'), plot.margin = margin(t = .3, b = .3, r = .5, unit = "cm"),
        axis.title.x = element_blank() , axis.title.y = element_text(face = 'bold', margin = margin(l = 0.3, r = 0.3, unit = "cm"), size=13)) + scale_x_discrete(labels = c('IP', 'D14', 'D21', 'D28', 'D28 BM')) + 
  scale_fill_manual(values = color_list$treat) + ggtitle('Cilta IP Large Clones')

p8 <- ggarrange(plotlist = list(p6, p7), nrow = 1, ncol =2)

ggsave('Fig5a/Fig5a_ide_ip_large.pdf', plot = p6, height = 2400, width = 1600, dpi = 500, units = 'px')
ggsave('Fig5a/Fig5a_cilta_ip_large.pdf', plot = p7, height = 2400, width = 2000, dpi = 500, units = 'px')

###----------------------------------------------Ide IP large clone ttype trace--------------------------------------####

ide_large_pheno <- dplyr::filter(pheno_clono_table, CTstrict %in% ide_cc_large$CTstrict & treatment == 'Ide')

ide_large_pheno_ttype <- dplyr::filter(ide_large_pheno, type %in% c('IP', 'CAR')) %>% group_by(time, sample,  annotype) %>% summarize(count = n()) %>% mutate(Prop = count/sum(count))  %>% ungroup()
ide_large_pheno_ttype$annotype <- as.character(ide_large_pheno_ttype$annotype)
ide_large_pheno_ttype$annotype <- factor(ide_large_pheno_ttype$annotype, levels = c('CD4', 'CD8'))
ide_large_pheno_ttype$time <- factor(ide_large_pheno_ttype$time, levels = c('IP', 'D9', 'D28'))

p9 <- ggplot(ide_large_pheno_ttype, aes(x= time, y=Prop, fill= annotype )) + geom_bar(stat = 'identity', color = 'black') + theme_pubr() + ylab('Proportion')+
  theme(legend.position = 'right', legend.text = element_text(size = 12), plot.title = element_text(hjust = .5, face = 'bold'), plot.margin = margin(t = .3, b = .3, r = .5, unit = "cm"),
        axis.title.x = element_blank() , axis.title.y = element_text(face = 'bold', margin = margin(l = 0.3, r = 0.3, unit = "cm"), size=13)) + scale_x_discrete(labels = c('IP', 'D9', 'D28')) + 
  scale_fill_manual(values = color_list$ttype) + ggtitle('Ide IP Large Clones Type') + labs(fill = '')

ggsave('Fig5a/Fig5a_ide_ip_large_ttype_trace.pdf', plot = p9, height = 2400, width = 1800, dpi = 500, units = 'px')
  
###----------------------------------------------Ide IP large clone anno1 trace--------------------------------------####

ide_large_pheno_anno1 <- dplyr::filter(ide_large_pheno, type %in% c('IP', 'CAR')) %>% group_by(time, sample,  anno1) %>% summarize(count = n()) %>% mutate(Prop = count/sum(count))  %>% ungroup()
ide_large_pheno_anno1$anno1_lab <- str_replace_all(ide_large_pheno_anno1$anno1, 'CD4_', 'CD4 ')
ide_large_pheno_anno1$anno1_lab <- str_replace_all(ide_large_pheno_anno1$anno1, 'CD8_', 'CD8 ')
ide_large_pheno_anno1$anno1 <- factor(ide_large_pheno_anno1$anno1, levels = c("CD4_Prolif", "CD4_EM", "CD4_Th1", "CD8_Prolif_EM", "CD8_Prolif_TE", "CD8_TE"))
ide_large_pheno_anno1$time <- factor(ide_large_pheno_anno1$time, levels = c('IP', 'D9', 'D28'))

p10 <- ggplot(ide_large_pheno_anno1, aes(x= time, y=Prop, fill= anno1 )) + geom_bar(stat = 'identity', color = 'black') + theme_pubr() + ylab('Proportion')+
  theme(legend.position = 'right', legend.text = element_text(size = 12), plot.title = element_text(hjust = .5, face = 'bold'), plot.margin = margin(t = .3, b = .3, r = .5, unit = "cm"),
        axis.title.x = element_blank() , axis.title.y = element_text(face = 'bold', margin = margin(l = 0.3, r = 0.3, unit = "cm"), size=13)) + scale_x_discrete(labels = c('IP', 'D9', 'D28')) + 
  scale_fill_manual(values = color_list$anno1[levels(ide_large_pheno_anno1$anno1)], labels = c("CD4 Prolif", "CD4 EM", "CD4 Th1", "CD8 Prolif_EM", "CD8 Prolif_TE", "CD8 TE")) + ggtitle('Ide IP Large Clones Phenotype') + labs(fill = '') 

ggsave('Fig5a/Fig5a_ide_ip_large_anno1_trace.pdf', plot = p10, height = 2400, width = 2200, dpi = 500, units = 'px')
###----------------------------------------------Cilta IP large clone ttype trace--------------------------------------####

cilta_large_pheno <- dplyr::filter(pheno_clono_table, CTstrict %in% cilta_cc_large$CTstrict & treatment == 'Cilta')
cilta_large_pheno$sample2 <- revalue(cilta_large_pheno$sample, c('Cilta_BM_CAR_1' = 'BM', 'Cilta_BM_CAR_2' = 'BM'))

cilta_large_pheno_ttype <- dplyr::filter(cilta_large_pheno, type %in% c('IP', 'CAR')) %>% group_by(sample2,  annotype) %>% summarize(count = n()) %>% mutate(Prop = count/sum(count))  %>% ungroup()
cilta_large_pheno_ttype$annotype <- as.character(cilta_large_pheno_ttype$annotype)
cilta_large_pheno_ttype$annotype <- factor(cilta_large_pheno_ttype$annotype, levels = c('CD4', 'CD8'))
cilta_large_pheno_ttype$sample2 <- factor(cilta_large_pheno_ttype$sample2, levels = c('Cilta_IP', 'Ciltacel_D14_CAR', 'Ciltacel_D21_CAR', 'Ciltacel_D28_CAR', 'BM'))

p11 <- ggplot(cilta_large_pheno_ttype, aes(x= sample2, y=Prop, fill= annotype )) + geom_bar(stat = 'identity', color = 'black') + theme_pubr() + ylab('Proportion')+
  theme(legend.position = 'right', legend.text = element_text(size = 12), plot.title = element_text(hjust = .5, face = 'bold'), plot.margin = margin(t = .3, b = .3, r = .5, unit = "cm"),
        axis.title.x = element_blank() , axis.title.y = element_text(face = 'bold', margin = margin(l = 0.3, r = 0.3, unit = "cm"), size=13)) + scale_x_discrete(labels = c('IP', 'D14', 'D21', 'D28', 'D28 BM')) + 
  scale_fill_manual(values = color_list$ttype) + ggtitle('Cilta IP Large Clones Type') + labs(fill = '')

ggsave('Fig5a/Fig5a_cilta_ip_large_ttype_trace.pdf', plot = p11, height = 2400, width = 2300, dpi = 500, units = 'px')

###----------------------------------------------Cilta IP large clone pheno trace--------------------------------------####

cilta_large_pheno_anno1 <- dplyr::filter(cilta_large_pheno, type %in% c('IP', 'CAR')) %>% group_by(sample2,  anno1) %>% summarize(count = n()) %>% mutate(Prop = count/sum(count))  %>% ungroup()
cilta_large_pheno_anno1$anno1 <- factor(cilta_large_pheno_anno1$anno1, levels = c("CD4_Prolif", "CD4_EM", 'CD8_Prolif', 'CD8_Prolif_EM', 'CD8_Prolif_TE', 'CD8_EM', 'CD8_ISG', 'CD8_IFNG', 'CD8_TE', 'CD8_IL'))
cilta_large_pheno_anno1$sample2 <- factor(cilta_large_pheno_anno1$sample2, levels = c('Cilta_IP', 'Ciltacel_D14_CAR', 'Ciltacel_D21_CAR', 'Ciltacel_D28_CAR', 'BM'))


p12 <- ggplot(cilta_large_pheno_anno1, aes(x= sample2, y=Prop, fill= anno1 )) + geom_bar(stat = 'identity', color = 'black') + theme_pubr() + ylab('Proportion')+
  theme(legend.position = 'right', legend.text = element_text(size = 12), plot.title = element_text(hjust = .5, face = 'bold'), plot.margin = margin(t = .3, b = .3, r = .5, unit = "cm"),
        axis.title.x = element_blank() , axis.title.y = element_text(face = 'bold', margin = margin(l = 0.3, r = 0.3, unit = "cm"), size=13)) + scale_x_discrete(labels =c('IP', 'D14', 'D21', 'D28', 'D28 BM')) + 
  scale_fill_manual(values = color_list$anno1[levels(cilta_large_pheno_anno1$anno1)], labels = c("CD4 Prolif", "CD4 EM", 'CD8 Prolif', 'CD8 Prolif_EM', 'CD8 Prolif_TE', 'CD8 EM', 'CD8 ISG', 'CD8 IFNG', 'CD8 TE', 'CD8 IL')) + ggtitle('Cilta IP Large Clones Phenotype') + labs(fill = '') 

ggsave('Fig5a/Fig5a_cilta_ip_large_anno1_trace.pdf', plot = p12, height = 2400, width = 2700, dpi = 500, units = 'px')
