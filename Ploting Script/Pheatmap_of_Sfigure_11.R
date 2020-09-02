###-------------------------1. suplementary feature pheatmap-------------###################
library(ggplot2)
library(stringr)
library(data.table)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(tidyverse)
library(pheatmap)
library(ggthemes)
library(ggsci)

meta_dat_all <- read.csv("Reannotation_HCA_alltissues_meta.data_84363_cell.txt",
                         header = T, stringsAsFactors = F,
                         sep = "\t",
                         row.names = 1
)

meta_dat_all$major_Cell_types_in_tissues <- meta_dat_all$re_Cell_type_of_all_data %>%
  gsub(pattern = ".*T Cell .*", replacement = "T Cell") %>% 
  gsub(pattern = ".*B Cell .*", replacement = "B Cell") %>% 
  gsub(pattern = ".*Macrophage .*", replacement = "Macrophage") %>% 
  gsub(pattern = ".*Fibroblast .*", replacement = "Fibroblast") %>% 
  gsub(pattern = ".*FibSmo .*", replacement = "FibSmo") %>% 
  gsub(pattern = ".*Endothelial .*", replacement = "Endothelial") %>% 
  gsub(pattern = "Plasma.*", replacement = "Plasma") 


tmp <- c()
for (i in meta_dat_all$orig.ident %>% unique) {
  tmp1 <- meta_dat_all %>% subset(orig.ident == i) %>% dplyr::select(Cell_type_of_each_tissue) %>% unlist %>% 
    gsub(pattern = ".*T Cell .*", replacement = "T Cell") %>% 
    gsub(pattern = ".*T GD Cell .*", replacement = "T Cell") %>% 
    gsub(pattern = ".*T/NK.*", replacement = "NK/T") %>% 
    gsub(pattern = ".*NK/T.*", replacement = "NK/T") %>% 
    gsub(pattern = ".*CD4 T Cell", replacement = "T Cell") %>% 
    gsub(pattern = ".*CD8 T .*", replacement = "T Cell") %>% 
    gsub(pattern = ".*B Cell .*", replacement = "B Cell") %>% 
    gsub(pattern = "Plasma .*", replacement = "Plasma Cell") %>% 
    gsub(pattern = ".*Macrophage .*", replacement = "Macrophage") %>% 
    gsub(pattern = ".*Fibroblast .*", replacement = "Fibroblast") %>% 
    gsub(pattern = ".*FibSmo .*", replacement = "FibSmo") %>% 
    gsub(pattern = "Endothelial .*", replacement = "Endothelial") %>% 
    gsub(pattern = ".*Smooth Muscle Cell.*", replacement = "Smooth Muscle Cell") %>% 
    gsub(pattern = " .*_cDNA", replacement = "") #%>% 
    tmp1 <- paste0(tmp1, " ", i %>% gsub(pattern = "_cDNA", replacement = ""))
    tmp <- c(tmp, tmp1)
}
  
meta_dat_all$Cell_types_in_tissue <- tmp

data_to_plot <- table(meta_dat_all[, c("Cell_types_in_tissue", "major_Cell_types_in_tissues")]) %>% unclass %>% as.data.frame()

# png("cell_type_cor.png", width = 10, height = 15, res = 500, units = "in")
pdf("Supplementary_Figure_11_cell_type_cor.pdf", height = 10, width = 7)
pheatmap::pheatmap(log2(data_to_plot[] + 1),
                   color = colorRampPalette(c("white", "red"))(5),
                   #color = colorRampPalette(brewer.pal(8, "PiYG"))(25),
                   #color = colorRampPalette(brewer.pal(8, "Blues"))(25),
                   #color = terrain.colors(256),
                   border_color = "white",
                   cluster_cols = T,
                   cluster_rows = T,
                   show_rownames = T,
                   show_colnames = T,
                   fontsize_row = 3,
                   fontsize_col = 5,
                   scale = "none",
                   silent = F,
                   treeheight_row = 0,
                   treeheight_col = 0
)
dev.off()
