###-------------------------1. suplementary feature pheatmap-------------###################
library(ggplot2)
library(stringr)
library(data.table)
library(dplyr)
library(plyr)
library(RColorBrewer)
# library(tidyverse)
library(pheatmap)
library(ggthemes)
library(ggsci)

meta_dat_all <- read.csv("BG.RV.meta_data_of_all_tissues_from_merged_data.txt",
                         header = T, stringsAsFactors = F,
                         sep = "\t",
                         row.names = 1
)

meta_dat_all$major_Cell_types_in_tissues <- meta_dat_all$Cell_types_combined_data %>%
  gsub(pattern = "T Cell .*", replacement = "T Cell") %>% 
  gsub(pattern = "B Cell .*", replacement = "B Cell") %>% 
  gsub(pattern = "Macrophage .*", replacement = "Macrophage") %>% 
  gsub(pattern = "Fibroblast .*", replacement = "Fibroblast") %>% 
  gsub(pattern = "FibSmo .*", replacement = "FibSmo") %>% 
  gsub(pattern = "Endothelial .*", replacement = "Endothelial") %>% 
  gsub(pattern = "Enterocyte .*", replacement = "Enterocyte") 

tmp <- c()
for (i in meta_dat_all$orig.ident %>% unique) {
  tmp1 <- meta_dat_all %>% subset(orig.ident == i) %>% dplyr::select(Cell_types_in_tissue) %>% unlist %>% 
    gsub(pattern = "T Cell .*", replacement = "T Cell") %>%
    gsub(pattern = "CD4 T Cell", replacement = "T Cell") %>% 
    gsub(pattern = "CD8 T Cell", replacement = "T Cell") %>% 
    gsub(pattern = "B Cell .*", replacement = "B Cell") %>% 
    gsub(pattern = "Plasma .*", replacement = "Plasma Cell") %>% 
    gsub(pattern = "Macrophage .*", replacement = "Macrophage") %>% 
    gsub(pattern = "Fibroblast .*", replacement = "Fibroblast") %>% 
    gsub(pattern = "FibSmo .*", replacement = "FibSmo") %>% 
    gsub(pattern = "Endothelial .*", replacement = "Endothelial") %>% 
    gsub(pattern = "Enterocyte .*", replacement = "Enterocyte") 
  tmp1 <- paste0(tmp1, " ", i)
  tmp <- c(tmp, tmp1)
}
  
meta_dat_all$Cell_types_in_tissue <- tmp


# meta_dat_all$major_Cell_types_all_cluster <- gsub(" [1-9]", "", meta_dat_all$cluster_all)

data_to_plot <- table(meta_dat_all[, c("Cell_types_in_tissue", "major_Cell_types_in_tissues")]) %>% unclass %>% as.data.frame()
data_to_plot <- data_to_plot[, c("B Cell",
                                 "Plasma Cell",
                                 "T Cell",
                                 "Fibroblast",
                                 "FibSmo",
                                 "Smooth Muscle Cell",
                                 "Endothelial",
                                 "Monocyte",
                                 "Macrophage",
                                 "Enterocyte",
                                 "Schwann Cell", colnames(data_to_plot) %>% setdiff(c("B Cell",
                                                                                       "Plasma Cell",
                                                                                       "T Cell",
                                                                                       "Fibroblast",
                                                                                       "FibSmo",
                                                                                       "Smooth Muscle Cell",
                                                                                       "Endothelial",
                                                                                       "Monocyte",
                                                                                       "Macrophage",
                                                                                       "Enterocyte",
                                                                                       "Schwann Cell")))]

# png("cell_type_cor.png", width = 10, height = 15, res = 500, units = "in")
pdf("Supplementary_Figure_1_G_cell_type_cor.pdf", height = 12, width = 7)
pheatmap::pheatmap(log2(data_to_plot[] + 1),
                   color = colorRampPalette(c("white", "red"))(5),
                   #color = colorRampPalette(brewer.pal(8, "PiYG"))(25),
                   #color = colorRampPalette(brewer.pal(8, "Blues"))(25),
                   #color = terrain.colors(256),
                   border_color = "white",
                   cluster_cols = F,
                   cluster_rows = T,
                   show_rownames = T,
                   show_colnames = T,
                   fontsize_row = 3,
                   fontsize_col = 5,
                   scale = "none",
                   silent = F,
                   treeheight_row = 0,
                   treeheight_col = 0
                   # border_color = "grey60"
                   # annotation_row = rep(c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
                   #                    '#98df8a', '#d62728', '#9edae5', '#ff9896', '#9467bd',
                   #                    '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2',
                   #                    '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d'),
                   #                    times = table(B_cell_clone_sort_h_2_l$Tissue)),
                   # annotation_col = rep(c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
                   #                        '#98df8a', '#d62728', '#9edae5', '#ff9896', '#9467bd',
                   #                        '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2',
                   #                        '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d'),
                   #                      times = table(B_cell_clone_sort_h_2_l$Tissue))
                   # cellheight = 3,
                   # cellwidth = 20
                   # filename = "Cell_type_cor.png"
)
dev.off()
