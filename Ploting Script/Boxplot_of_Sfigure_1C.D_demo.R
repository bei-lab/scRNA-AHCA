#####################-------------------1. suplementary UMI and Genes-------------###################
library(ggplot2)
library(stringr)
library(data.table)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(pheatmap)
library(ggthemes)
library(ggsci)

tissues_colors <- c(    "Common.bile.duct" = '#1f77b4',
                        Bladder = '#aec7e8',
                        Blood = '#ff7f0e',
                        Esophagus = '#ffbb78',
                        Heart = '#2ca02c',
                        # Kidney = '#98df8a',
                        Rectum = '#d62728',
                        Liver = '#9edae5',
                        # Lung = '#ff9896',
                        "Lymph.node" = '#9467bd',
                        Marrow = '#c5b0d5',
                        Muscle = '#8c564b',
                        # Pancreas = '#c49c94',
                        Skin = '#e377c2',
                        "Small.intestine" = '#f7b6d2',
                        Spleen = '#7f7f7f',
                        Stomach = '#c7c7c7',
                        # Testis = '#17becf',
                        Trachea = '#d62790'
)

summary_of_genes_and_UMI_after_filtering <- read.table("Reannotation_HCA_alltissues_meta.data_84363_cell.txt",
                                                       header = T,
                                                       row.names = 1,
                                                       sep = "\t",
                                                       stringsAsFactors = F,
                                                       comment.char = "")
                                                       
##-------------------------plot number of genes in each organ-----------------##

meta_data <- summary_of_genes_and_UMI_after_filtering
meta_data <- meta_data[, c(1, 2, 3)]
meta_data_melt <- as.data.frame(melt(meta_data))
meta_data_melt <- meta_data[, c(3,1)]
colnames(meta_data_melt) <- c("Genes", "Organ")

meta_data_melt$Organ <- factor(meta_data_melt$Organ, levels = sort(unique(meta_data_melt$Organ), decreasing = T))

png("Supplementary_figure1_1A.png", width = 10, height = 15, res = 400, units = "in")
ggplot(data = meta_data_melt, aes(x = Organ, y = (Genes))) + 
  geom_boxplot(aes(fill = Organ)) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"), axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
  ylab("Genes") +
  theme(legend.justification=c(0.5,0), legend.position=c(0.9, 0.2)) +
  scale_x_discrete(limits = (levels(meta_data_melt$Organ))) +
  scale_fill_manual(values = tissues_colors[tissues_colors %>% names %>% order()]) +
  scale_y_continuous(expand = c(0.03, 0.5), breaks = c(1000, 2000, 3000, 4000, 5000, 6000)) +
  coord_flip()
dev.off()

##-------------------------plot number of UMI in each organ-----------------##

meta_data <- summary_of_genes_and_UMI_after_filtering
meta_data <- meta_data[, c(1, 2, 3)]
meta_data_melt <- as.data.frame(melt(meta_data))
meta_data_melt <- meta_data[, c(2,1)]
colnames(meta_data_melt) <- c("UMI", "Organ")

meta_data_melt$Organ <- factor(meta_data_melt$Organ, levels = sort(unique(meta_data_melt$Organ), decreasing = T))

png("Supplementary_figure1_1B.png", width = 10, height = 15, res = 400, units = "in")
ggplot(data = meta_data_melt, aes(x = Organ, y = (UMI))) + 
  geom_boxplot(aes(fill = Organ)) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"), axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
  ylab("UMI") +
  theme(legend.justification=c(0.5,0), legend.position=c(0.9, 0.2)) +
  scale_x_discrete(limits = (levels(meta_data_melt$Organ))) +
  scale_fill_manual(values = tissues_colors[tissues_colors %>% names %>% order()]) +
  scale_y_continuous(expand = c(0.03, 0.5), breaks = c(1000, 2000, 5000, 10000, 20000, 30000, 40000, 50000, 60000)) +
  coord_flip()
dev.off()

