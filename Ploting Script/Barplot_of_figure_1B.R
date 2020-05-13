library(ggplot2)
library(stringr)
library(data.table)
library(dplyr)
library(plyr)
library(RColorBrewer)
# library(tidyverse)
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
                        # x1 = '#00cdaa',
                        # x2 = '#aacdaa',
                        # x3 = '#aaaa00',
                        # x4 = '#ff00ff',
                        # x5 = '#00acff',
                        # x6 = '#ccacff'
)

####################---------------read meta.data-----------------------################################
meta.data <- read.table(file = "GB_RV.cell_of_each_tissue.meta.data.txt", ##meta.data after filtering 
                        sep = "\t", 
                        header = T, 
                        col.names = c("Tissue", "Number of UMI", "Number of Genes", "Cluster", "Organ", "Cell_types", "cellbarcode"),
                        comment.char = "",
                        row.names = 1,
                        stringsAsFactors = F)

meta.data$Organ <- (meta.data$cellbarcode %>% as.character %>% str_split(., "_", simplify = T))[, 1]

#####################------------------- figure 1b -------------###################
#1. cells number of each tissue
cell_of_each_tissue <- unclass(table(meta.data$Organ))
dat <- data.frame(Organ = names(cell_of_each_tissue), Cell_numbers = unname(cell_of_each_tissue))
dat$Organ <- factor(dat$Organ, levels = sort(dat$Organ, decreasing = F))


# png("Figure1_B.png", width = 15, height = 15, res = 400, units = "in" )
pdf("Figure1_B.pdf", width = 15, height = 15)
ggplot(data = dat, aes(x = Organ, y = Cell_numbers, fill = Organ)) + 
  geom_bar(stat="identity", position = position_dodge(width = 0.5), width = 0.8) + 
  # geom_text(aes(y= Cell_numbers/2, label = Cell_numbers), hjust = 0,
  #           color="black", size=4) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.8, linetype = "solid")) +
  scale_x_discrete(name = "Organ") +
  scale_y_discrete(name ="Cells number", 
                   limits=c(0, 3000, 6000, 9000, 12000, 14000)) +
  theme(legend.background = element_rect(fill= NULL, size=.5, linetype= NULL)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0.001),
        axis.title.x = element_text(colour = "black", size = 12, family="Helvetica"),
        axis.title.y = element_text(colour = "black", size = 12, family="Helvetica"),
        # axis.ticks.length = unit(.25, "cm"),
        axis.text = element_text(colour = "black", face = "bold", family="Helvetica")) +
  scale_fill_manual(values = tissues_colors[order(names(tissues_colors), decreasing = T)]) +
  theme(legend.position = "none") +
  scale_x_discrete(limits = rev(levels(dat$Organ))) +
  coord_flip()
dev.off()
