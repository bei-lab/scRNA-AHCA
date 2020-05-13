##----------------This script was used to identify T cell with paired A and B chains.-------####
##---------------Take CD8 T cell for an example.--------------########
library(Seurat)
library(ggpubr)
library(plyr)
library(dplyr)
library(stringr)
library(ggthemes)
library(cowplot)
library(data.table)
library(RColorBrewer)
library(pheatmap)
library(tidyr)
library(Startrac)
library(ggplot2)
library(ggsci)
library(igraph)

tissues_colors <- c(    "Common.bile.duct" = '#1f77b4',
                        Bladder = '#aec7e8',
                        Blood = '#ff7f0e',
                        Esophagus = '#ffbb78',
                        Heart = '#2ca02c',
                        Kidney = '#98df8a',
                        Rectum = '#d62728',
                        Liver = '#9edae5',
                        Lung = '#ff9896',
                        "Lymph.node" = '#9467bd',
                        Marrow = '#c5b0d5',
                        Muscle = '#8c564b',
                        Pancreas = '#c49c94',
                        Skin = '#e377c2',
                        "Small.intestine" = '#f7b6d2',
                        Spleen = '#7f7f7f',
                        Stomach = '#c7c7c7',
                        Testis = '#17becf',
                        Trachea = '#d62790'
                        # x1 = '#00cdaa',
                        # x2 = '#aacdaa',
                        # x3 = '#aaaa00',
                        # x4 = '#ff00ff',
                        # x5 = '#00acff',
                        # x6 = '#ccacff'
)
tissus_to_numbers <- c( "Common.bile.duct" = 1,
                        Bladder = 1,
                        Blood = 1,
                        Esophagus = 1,
                        Heart = 1,
                        Kidney = 3,
                        Rectum = 1,
                        Liver = 1,
                        Lung = 2,
                        "Lymph.node" = 1,
                        Marrow = 1,
                        Muscle = 1,
                        Pancreas = 2,
                        Skin = 1,
                        "Small.intestine" = 1,
                        Spleen = 1,
                        Stomach = 1,
                        Testis = 1,
                        Trachea = 1
                        # x1 = '#00cdaa',
                        # x2 = '#aacdaa',
                        # x3 = '#aaaa00',
                        # x4 = '#ff00ff',
                        # x5 = '#00acff',
                        # x6 = '#ccacff'
)

#######################################--------------------------------------################################S
number_to_tissue <- read.table("number_corresponding_tissue_BCR.txt", header = F, col.names = c("Number", "Tissue"), stringsAsFactors = F)
T_cells_meta.data <- read.table("CD8_meta.data.csv", header = T, row.names = 1, sep = "\t", stringsAsFactors = F, comment.char = "")

#################
T_cell_clone_orig <- read.table("filtered_contig_annotations_T.csv", sep = ",", header = T, stringsAsFactors = F)
number_cor_tissue_T_TCR <- str_split(T_cell_clone_orig$barcode, pattern = "-", simplify = T)
number_cor_tissue_T_TCR[, 2] <- mapvalues(number_cor_tissue_T_TCR[, 2], from = number_to_tissue$Number, to = number_to_tissue$Tissue )
T_cell_clone_orig$Tissue <- number_cor_tissue_T_TCR[, 2]
T_cell_clone_orig$Barcode <- number_cor_tissue_T_TCR[, 1]
T_cell_clone_orig$cell_barcode <- (with(T_cell_clone_orig, paste0(Tissue, "_cDNA_", Barcode)))
################

T_cell_clone <- read.table("filtered_contig_annotations_T.csv", sep = ",", header = T, stringsAsFactors = F)

###-------------------remove the cells with only one chain---------#################
multiple_chains <- T_cell_clone[duplicated(T_cell_clone$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone <- T_cell_clone[T_cell_clone$barcode %in% multiple_chains, ]

###-------------------transform the barcode of cells---------------#################
number_cor_tissue_T <- str_split(T_cell_clone$barcode, pattern = "-", simplify = T)
number_cor_tissue_T[, 2] <- mapvalues(number_cor_tissue_T[, 2], from = number_to_tissue$Number, to = number_to_tissue$Tissue )
T_cell_clone$Tissue <- number_cor_tissue_T[, 2]
T_cell_clone$Barcode <- number_cor_tissue_T[, 1]

T_cell_clone$cell_barcode <- (with(T_cell_clone, paste0(Tissue, "_cDNA_", Barcode)))

###--------------------remove contaminated cells with BCR,D&G T cells, raw_clonetype_id None, cdr3 None, productive None--------------------------------------#######################
contamination_cells_T <- unique(T_cell_clone[(grepl(pattern = "IG|Multi|TRD|TRG",
                                                    x = with(T_cell_clone, paste0(chain, v_gene, d_gene, j_gene, c_gene)))),
                                             "cell_barcode"]) ## remove contaminated cells with BCR,D&G T cells

T_cell_clone <- T_cell_clone[!(T_cell_clone$cell_barcode %in% contamination_cells_T), ] ## remove contaminated cells

T_cell_clone_none_rm <- T_cell_clone[!T_cell_clone$raw_clonotype_id == "None", ] ## remove raw_clonotype_id None
T_cell_clone_none_rm <- T_cell_clone_none_rm[!T_cell_clone_none_rm$cdr3 == "None", ] ## remove the cdr3 None
T_cell_clone_none_rm <- T_cell_clone_none_rm[!T_cell_clone_none_rm$productive == "False", ] ## remove the productive False

###-------------------remove the cells with only one chain---------#################
multiple_chains <- T_cell_clone_none_rm[duplicated(T_cell_clone_none_rm$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone <- T_cell_clone_none_rm[T_cell_clone_none_rm$barcode %in% multiple_chains, ]
T_cell_clone_uniq <- T_cell_clone

###--------------------remove T cells only with A or B chain---------------------------------#############
barcodes_cells <- unique(T_cell_clone_uniq$cell_barcode)

qualited_T_cells <- lapply(barcodes_cells, FUN = function(x){
  tmp <- T_cell_clone[T_cell_clone_uniq$cell_barcode  %in% x, "chain"] %>% unlist
  if(all("TRA" %in% tmp, "TRB" %in% tmp)){
    return(x)
  }
}) %>% do.call(rbind, .)

T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$cell_barcode %in% qualited_T_cells, ]

###------------------remove non T cells TCR----------------------------#####################
T_cell_clone_uniq <- filter(T_cell_clone_uniq, cell_barcode %in% row.names(T_cells_meta.data)) ## 5560 T cells (2402 clones)  with certain clone type and paired A and B chain
T_cell_clone_uniq$annotation <- mapvalues(T_cell_clone_uniq$cell_barcode, 
                                          from = T_cells_meta.data %>% row.names,
                                          to = T_cells_meta.data$annotation)

