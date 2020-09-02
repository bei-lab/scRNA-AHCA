###-------------------Sfigure 2 bar ploting-------------------------

library(ggplot2)
library(stringr)
library(data.table)
library(dplyr)
library(plyr)
library(RColorBrewer)
# library(tidyverse)
library(ggthemes)
library(ggsci)

meta.data <- read.table(file = "Reannotation_HCA_alltissues_meta.data_84363_cell.txt", #"/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/8-reanalysis/Figure_1/meta.data_of_all_cells.txt",  ##meta.data before filtering
                        sep = "\t", 
                        header = T, 
                        # col.names = c("Tissue", "Number of UMI", "Number of Genes", "Cluster", "Organ", "Cell_types"),
                        comment.char = "",
                        row.names = 1,
                        stringsAsFactors = F)
meta.data$orig.ident <- factor(meta.data$orig.ident %>% as.character(),
                               levels = meta.data$orig.ident %>% as.character %>% unique %>% sort)


by(meta.data %>%
   dplyr::select(Cell_types),
   INDICES = meta.data$orig.ident,
   FUN = function(x) names(table(x))[order(table(x), decreasing = T)],
   simplify = FALSE) -> cell_types_order

celltypes_orders <- do.call(c, cell_types_order)

color_to_use <- c(pal_npg("nrc")(10)[-8], pal_aaas("default")(10), pal_nejm("default")(8))
cols <- do.call(c, lapply(do.call(c,lapply(cell_types_order, length)), FUN = function(x) color_to_use[1:x]))

cell_types_to_cols <- data.frame(cols = cols, celltypes = celltypes_orders)

meta.data$Cell_types <- factor(meta.data$Cell_types, levels = rev(celltypes_orders))
# meta.data$tissue <- factor(meta.data$orig.ident, levels = rev(meta.data %>% select(orig.ident) %>% unique() %>% `[`(,)))
meta.data %>% dplyr::select(orig.ident) %>% unique() %>% `[`(1:8,) -> Tentissues


###-------------8 tissues barplot-----------------


ggplot(meta.data %>% subset(orig.ident %in% Tentissues) %>% dplyr::select(Cell_types, orig.ident), aes(Cell_types, fill = Cell_types)) + geom_bar(width = 0.7) +
  coord_flip() + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 4)) +
  guides(fill = FALSE) +
  scale_fill_manual(values =  cell_types_to_cols %>%
                      subset(gsub(pattern = "\\d+", replacement = "", row.names(cell_types_to_cols)) %in% Tentissues) %>%
                      dplyr::select(cols) %>% `[`(,) %>% as.character %>% rev) +
  scale_y_continuous(expand = c(0.01, 0.01),  name="Number of cells", breaks = c(0, 300, 600, 900, 1200, 1500, 3000)) +
  theme(axis.ticks.y = element_blank()) + 
  geom_hline(yintercept = c(300, 600, 900, 1200, 1500), colour = "white")

###-------------7 tissues barplot-----------------

meta.data %>% dplyr::select(orig.ident) %>% unique() %>% `[`(9:15,) -> Ninetissues
ggplot(meta.data %>% subset(orig.ident %in% Ninetissues) %>% dplyr::select(Cell_types, orig.ident), aes(Cell_types, fill = Cell_types)) + geom_bar(width = 0.7) +
  coord_flip() + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 4)) +
  guides(fill = FALSE) +
  scale_fill_manual(values =  cell_types_to_cols %>%
                      subset(gsub(pattern = "\\d+", replacement = "", row.names(cell_types_to_cols)) %in% Ninetissues) %>%
                      dplyr::select(cols) %>% `[`(,) %>% as.character %>% rev) +
  scale_y_continuous(expand = c(0.01, 0.01), name="Number of cells", breaks = c(0, 300, 600, 900, 1200, 1500, 3000)) +
  theme(axis.ticks.y = element_blank()) + 
  geom_hline(yintercept = c(300, 600, 900, 1200, 1500), colour = "white") 

