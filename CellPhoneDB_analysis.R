library(data.table)
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(plyr)
library(pheatmap)

cellphone_dat <- read.table("Receptor_and_ligation_significant.txt",
                            header = T,
                            row.names = NULL,
                            sep = "\t",
                            stringsAsFactors = F)

cellphone_dat$A_B <- cellphone_dat$interacting_pair
cellphone_dat <- cbind(cellphone_dat, str_split(cellphone_dat$interacting_pair, pattern = "_", simplify = T)[, 1:2] %>% as.data.frame %>% mutate(B_A = paste0(V2, "_", V1)))

number_fo_significant_RL <- table(cellphone_dat$interacting_pair) %>% sort(decreasing = T) %>% length() ##460

table(cellphone_dat$tissue) %>% sort(decreasing = T)

receptors_and_cell_types <- cellphone_dat[, 1:2]
receptors_and_cell_types$comb <- with(receptors_and_cell_types, paste(interacting_pair," ", variable))

unique_receptors_and_cell_types <- receptors_and_cell_types[!duplicated(receptors_and_cell_types$comb), ]

dat_for_plot <- table(unique_receptors_and_cell_types$interacting_pair) %>% unclass %>% as.data.frame() 
dat_for_plot <-  data.frame(name = row.names(dat_for_plot), frequence = dat_for_plot$.) %>% subset(frequence >= 50) 

##-----------------frequence of each paired ligation and receptors across different tissues---------------------### 

pdf("receptor_ligation_barplot.pdf", height = 10, width = 7)
ggplot(dat_for_plot, mapping = aes(x = reorder(name, frequence), y = frequence)) + 
  geom_bar(stat = "identity", fill = "#4DBBD5FF", width = 0.75) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, angle = 0, size = 8),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 6)) +
  scale_y_continuous(expand = c(0.01, 0.5)) + 
  coord_flip() +
  labs(x = "ligation_receptors")
dev.off()

##---------------number of paired ligation and receptors in each tissue---------#####

pdf("receptor_ligation_barplot_in_each_tissue.pdf", height = 10, width = 15)
cellphone_dat %>% ggplot(., aes(x = tissue)) +
  geom_bar(width = 0.7, fill = "#4DBBD5FF") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1)) +
  scale_y_continuous(expand = c(0.002, 0), limits = c(0, 3000), breaks = seq(0, 3000, 500))
dev.off()
  
###-----------------------------------------prepare the files for cytoscape---------------####
celltype_celltypes_annotation <- cellphone_dat %>% 
  mutate(Type1 = part1 %>% 
           gsub(pattern = ".*high|HPro_Esophagus|HPro_Skin|Stem|Absorptive|Follicular|Novel|Simple|Goblet|Pit|Secretory|Granular|Tuft|Ionocyte|Chief", replacement = "Epi"),
         Type2 = part2 %>%
           gsub(pattern = ".*high|HPro_Esophagus|HPro_Skin|Stem|Absorptive|Follicular|Novel|Simple|Goblet|Pit|Secretory|Granular|Tuft|Ionocyte|Chief", replacement = "Epi"))

celltype_celltypes_annotation$Type1 <- celltype_celltypes_annotation$Type1 %>%  gsub(., pattern = ".*CD4$", replacement = "CD4")
celltype_celltypes_annotation$Type2 <- celltype_celltypes_annotation$Type2 %>%  gsub(., pattern = ".*CD4$", replacement = "CD4")

celltype_celltypes_annotation$Type1 <- celltype_celltypes_annotation$Type1 %>%  gsub(., pattern = ".*CD8$", replacement = "CD8")
celltype_celltypes_annotation$Type2 <- celltype_celltypes_annotation$Type2 %>%  gsub(., pattern = ".*CD8$", replacement = "CD8")

celltype_celltypes_annotation$Type1 <- celltype_celltypes_annotation$Type1 %>%  gsub(., pattern = ".*Memory_B$", replacement = "B")
celltype_celltypes_annotation$Type1 <- celltype_celltypes_annotation$Type1 %>%  gsub(., pattern = ".*Naive_B$", replacement = "B")
celltype_celltypes_annotation$Type2 <- celltype_celltypes_annotation$Type2 %>%  gsub(., pattern = ".*Memory_B$", replacement = "B")
celltype_celltypes_annotation$Type2 <- celltype_celltypes_annotation$Type2 %>%  gsub(., pattern = ".*Naive_B$", replacement = "B")

celltype_celltypes_annotation$Type1 <- celltype_celltypes_annotation$Type1 %>%  gsub(., pattern = ".*LEC$", replacement = "Endo")
celltype_celltypes_annotation$Type1 <- celltype_celltypes_annotation$Type1 %>%  gsub(., pattern = ".*BEC$", replacement = "Endo")
celltype_celltypes_annotation$Type2 <- celltype_celltypes_annotation$Type2 %>%  gsub(., pattern = ".*LEC$", replacement = "Endo")
celltype_celltypes_annotation$Type2 <- celltype_celltypes_annotation$Type2 %>%  gsub(., pattern = ".*BEC$", replacement = "Endo")

celltype_celltypes_annotation$Type1 <- celltype_celltypes_annotation$Type1 %>%  gsub(., pattern = ".*Mon$", replacement = "Myeloid")
celltype_celltypes_annotation$Type1 <- celltype_celltypes_annotation$Type1 %>%  gsub(., pattern = ".*Mac$", replacement = "Myeloid")
celltype_celltypes_annotation$Type1 <- celltype_celltypes_annotation$Type1 %>%  gsub(., pattern = ".*DC$", replacement = "Myeloid")
celltype_celltypes_annotation$Type2 <- celltype_celltypes_annotation$Type2 %>%  gsub(., pattern = ".*Mon$", replacement = "Myeloid")
celltype_celltypes_annotation$Type2 <- celltype_celltypes_annotation$Type2 %>%  gsub(., pattern = ".*Mac$", replacement = "Myeloid")
celltype_celltypes_annotation$Type2 <- celltype_celltypes_annotation$Type2 %>%  gsub(., pattern = ".*DC$", replacement = "Myeloid")


celltype_celltypes_annotation <- data.frame(celltype_celltypes_annotation %>% dplyr::select(interacting_pair, A_B, V1, V2, Type1, Type2)) %>% mutate(Type = paste0(Type1, "_", Type2))
row.names(celltype_celltypes_annotation) <- 1:dim(celltype_celltypes_annotation)[1]

Tcells_CD4 <- celltype_celltypes_annotation[grep(celltype_celltypes_annotation$Type, pattern = "CD4", value = F), ]
Tcells_CD8 <- celltype_celltypes_annotation[grep(celltype_celltypes_annotation$Type, pattern = "CD8", value = F), ]
Bcells <- celltype_celltypes_annotation[grep(celltype_celltypes_annotation$Type, pattern = "B", value = F), ]
Plasma <- celltype_celltypes_annotation[grep(celltype_celltypes_annotation$Type, pattern = "Plasma", value = F), ]
Fib <- celltype_celltypes_annotation[grep(celltype_celltypes_annotation$Type, pattern = "^Fib_|_Fib$", value = F), ]
Smo <- celltype_celltypes_annotation[grep(celltype_celltypes_annotation$Type, pattern = "^Smo_|_Smo$", value = F), ]
FibSmo <- celltype_celltypes_annotation[grep(celltype_celltypes_annotation$Type, pattern = "FibSmo", value = F), ]
Endo <- celltype_celltypes_annotation[grep(celltype_celltypes_annotation$Type, pattern = "Endo", value = F), ]
Myeloid <- celltype_celltypes_annotation[grep(celltype_celltypes_annotation$Type, pattern = "Myeloid", value = F), ]
NK <- celltype_celltypes_annotation[grep(celltype_celltypes_annotation$Type, pattern = "NK", value = F), ]
Epi <- celltype_celltypes_annotation[grep(celltype_celltypes_annotation$Type, pattern = "Epi", value = F), ]


for (cell in c("Tcells_CD4", "Tcells_CD8", "Bcells", "Plasma", "Fib", "Smo", "FibSmo", "Myeloid", "NK", "Epi", "Endo")) {
  Tcells_CD4_for_cytoscape <- get(cell) %>% dplyr::select(Type) %>% table() %>% unclass %>% as.data.frame()
  Tcells_CD4_for_cytoscape <- Tcells_CD4_for_cytoscape %>% transmute(from = row.names(Tcells_CD4_for_cytoscape) %>% str_split(. ,pattern = "_", simplify = T) %>% `[`(, 1),
                                                                     to = row.names(Tcells_CD4_for_cytoscape) %>% str_split(. ,pattern = "_", simplify = T)%>% `[`(, 2),
                                                                     number = Tcells_CD4_for_cytoscape$.)
  
  write.table(Tcells_CD4_for_cytoscape, paste0(cell, "_for_cytoscape.txt"), sep = "\t", quote = F, row.names = F)
  
}

##----------summarize the number of interactions for each cell type
tmp <- data.frame()
for (cell in c("Tcells_CD4", "Tcells_CD8", "Bcells", "Plasma", "Fib", "Smo", "FibSmo", "Myeloid", "NK", "Epi", "Endo")) {
  Tcells_CD4_for_cytoscape <- get(cell) %>% dplyr::select(Type) %>% table() %>% unclass %>% as.data.frame()
  Tcells_CD4_for_cytoscape <- Tcells_CD4_for_cytoscape %>% transmute(from = row.names(Tcells_CD4_for_cytoscape) %>% str_split(. ,pattern = "_", simplify = T) %>% `[`(, 1),
                                                                     to = row.names(Tcells_CD4_for_cytoscape) %>% str_split(. ,pattern = "_", simplify = T)%>% `[`(, 2),
                                                                     comb = paste0(from, "_", to),
                                                                     number = Tcells_CD4_for_cytoscape$.)
  cells <- table(Tcells_CD4_for_cytoscape$from) %>% sort(decreasing = T) %>% `[`(1) %>% names
  row.names(Tcells_CD4_for_cytoscape) <- 1:dim(Tcells_CD4_for_cytoscape)[1]
  
  exchanged_position <- with(Tcells_CD4_for_cytoscape, grep(from, pattern = paste0("^", cells, "$"), value = F, fixed = F))
  from_1 <- Tcells_CD4_for_cytoscape[exchanged_position, "from"] 
  to_1 <- Tcells_CD4_for_cytoscape[exchanged_position, "to"] 
  
  ## begin to exchange
  Tcells_CD4_for_cytoscape[exchanged_position, "from"] <- to_1
  Tcells_CD4_for_cytoscape[exchanged_position, "to"] <- from_1
  
  result <- Tcells_CD4_for_cytoscape %>% dplyr::group_by(from, to) %>% dplyr::summarise(summary = sum(number), n = n()) %>% as.data.frame()
  tmp <- rbind(tmp, result)
}

write.table(tmp, "summary_dat_for_cytoscape.txt", sep = "\t", quote = F, row.names = F)

###------------------modify Receptor_and_ligation_all.txt data-----##############

all_cellphone_dat <-  read.table("Receptor_and_ligation_all_demo.txt",
                                    header = T,
                                    row.names = 1,
                                    sep = "\t",
                                    stringsAsFactors = F)
 
all_cellphone_dat <- all_cellphone_dat %>% 
   mutate(Type1 = part1 %>% 
           gsub(pattern = ".*high|HPro_Esophagus|HPro_Skin|Stem|Absorptive|Follicular|Novel|Simple|Goblet|Pit|Secretory|Granular|Tuft|Ionocyte|Chief", replacement = "Epi"),
           Type2 = part2 %>%
           gsub(pattern = ".*high|HPro_Esophagus|HPro_Skin|Stem|Absorptive|Follicular|Novel|Simple|Goblet|Pit|Secretory|Granular|Tuft|Ionocyte|Chief", replacement = "Epi"))


all_cellphone_dat$variable <- gsub(x = all_cellphone_dat$variable, pattern = ".Fib.Smo$", replacement = ".FibSmo")
all_cellphone_dat$variable <- gsub(x = all_cellphone_dat$variable, pattern = "^Fib.Smo.", replacement = "FibSmo.")
all_cellphone_dat$variable <- gsub(x = all_cellphone_dat$variable, pattern = "D.G_CD8", replacement = "DG_CD8")
all_cellphone_dat$variable <- gsub(x = all_cellphone_dat$variable, pattern = "Common.bile.duct", replacement = "Commonbileduct")
all_cellphone_dat$variable <- gsub(x = all_cellphone_dat$variable, pattern = "common.bile.duct", replacement = "Commonbileduct")
all_cellphone_dat$variable <- gsub(x = all_cellphone_dat$variable, pattern = "Small.intestine", replacement = "Smallintestine")

all_cellphone_dat$Type1 <- str_split(all_cellphone_dat$variable, "\\.", simplify = T)[, 1]
all_cellphone_dat$Type2 <- str_split(all_cellphone_dat$variable, "\\.", simplify = T)[, 2]

all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "Mac")] <- "Myeloid"
all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "Mon")] <- "Myeloid"
all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "Plasma")] <- "Plasma"
all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "Memory_B")] <- "B"
all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "Naive_B")] <- "B"
all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "LEC")] <- "Endo"
all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "BEC")] <- "Endo"
all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "^Smo$")] <- "Smo"
all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "FibSmo")] <- "FibSmo"
all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "^Fib$")] <- "Fib"
all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "DC")] <- "Myeloid"
all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "CD4")] <- "CD4"
all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "CD8")] <- "CD8"
all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "NK")] <- "NK"
all_cellphone_dat$cell1[grepl(all_cellphone_dat$Type1, pattern = "Epi")] <- "Epi"

all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "Mac")] <- "Myeloid" 
all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "Mon")] <- "Myeloid"
all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "Plasma")] <- "Plasma"
all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "Memory_B")] <- "B"
all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "Naive_B")] <- "B"
all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "LEC")] <- "Endo"
all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "BEC")] <- "Endo"
all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "^Smo$")] <- "Smo"
all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "FibSmo")] <- "FibSmo"
all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "^Fib$")] <- "Fib"
all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "DC")] <- "Myeloid"
all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "CD4")] <- "CD4"
all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "CD8")] <- "CD8"
all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "NK")] <- "NK"
all_cellphone_dat$cell2[grepl(all_cellphone_dat$Type2, pattern = "Epi")] <- "Epi"
 
all_cellphone_dat$Type <- paste0(all_cellphone_dat$cell1, "_", all_cellphone_dat$cell2)

write.table(all_cellphone_dat, "Receptor_and_ligation_all_modified.txt", sep = "\t", quote = F)

##########-------------------cell-cell interaction pheatmap----------------------------------#############

cellphone_dat <- read.table("Receptor_and_ligation_significant.txt",
                            header = T,
                            row.names = NULL,
                            sep = "\t",
                            stringsAsFactors = F)

cellphone_dat$A_B <- cellphone_dat$interacting_pair
cellphone_dat <- cbind(cellphone_dat, str_split(cellphone_dat$interacting_pair, pattern = "_", simplify = T)[, 1:2] %>% as.data.frame %>% mutate(B_A = paste0(V2, "_", V1)))

cellphone_dat <- cellphone_dat %>%
  mutate(Type1 = part1 %>%
           gsub(pattern = ".*high|HPro_Esophagus|HPro_Skin|Stem|Absorptive|Follicular|Novel|Simple|Goblet|Pit|Secretory|Granular|Tuft|Ionocyte|Chief", replacement = "Epi"),
         Type2 = part2 %>%
           gsub(pattern = ".*high|HPro_Esophagus|HPro_Skin|Stem|Absorptive|Follicular|Novel|Simple|Goblet|Pit|Secretory|Granular|Tuft|Ionocyte|Chief", replacement = "Epi"))

cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "Mac")] <- "Myeloid"
cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "Mon")] <- "Myeloid"
cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "Plasma")] <- "Plasma"
cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "Memory_B")] <- "B"
cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "Naive_B")] <- "B"
cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "LEC")] <- "Endo"
cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "BEC")] <- "Endo"
cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "^Smo$")] <- "Smo"
cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "FibSmo")] <- "FibSmo"
cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "^Fib$")] <- "Fib"
cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "DC")] <- "Myeloid"
cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "CD4")] <- "CD4"
cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "CD8")] <- "CD8"
cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "NK")] <- "NK"
cellphone_dat$cell1[grepl(cellphone_dat$Type1, pattern = "Epi")] <- "Epi"

cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "Mac")] <- "Myeloid"
cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "Mon")] <- "Myeloid"
cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "Plasma")] <- "Plasma"
cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "Memory_B")] <- "B"
cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "Naive_B")] <- "B"
cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "LEC")] <- "Endo"
cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "BEC")] <- "Endo"
cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "^Smo$")] <- "Smo"
cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "FibSmo")] <- "FibSmo"
cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "^Fib$")] <- "Fib"
cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "DC")] <- "Myeloid"
cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "CD4")] <- "CD4"
cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "CD8")] <- "CD8"
cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "NK")] <- "NK"
cellphone_dat$cell2[grepl(cellphone_dat$Type2, pattern = "Epi")] <- "Epi"

all_cellphone_dat <-  read.table("Receptor_and_ligation_all_modified.txt",
                                 header = T,
                                 row.names = 1,
                                 sep = "\t",
                                 stringsAsFactors = F)

###. find the significant paired RL-------------
T_epi <- data.frame(cellphone_dat %>% dplyr::select(interacting_pair, variable, A_B, V1, V2, cell1, cell2, pvalue, mean)) %>% mutate(Type = paste0(cell1, "_", cell2))
row.names(T_epi) <- 1:dim(T_epi)[1]
T_epi <- T_epi[grep(T_epi$Type, pattern = "CD8", value = F), ]
row.names(T_epi) <- 1:dim(T_epi)[1]
T_epi <- T_epi[grep(T_epi$Type, pattern = "Myeloid", value = F), ]
row.names(T_epi) <- 1:dim(T_epi)[1]

T_epi_pvalue <- T_epi %>% dplyr::select(interacting_pair, variable, pvalue)
paired_RL_names <- T_epi_pvalue %>% dplyr::select(interacting_pair) %>% unlist() %>% unname %>% unique %>% sort()

###.2 extracted the significant paired RL from raw data-------------
T_epi <- data.frame(all_cellphone_dat %>%
                      dplyr::select(interacting_pair, variable, cell1, cell2, Type, pvalue, mean, tissue)) %>%
  mutate(Types = paste0(variable, ".", tissue ))

row.names(T_epi) <- 1:dim(T_epi)[1]
T_epi <- T_epi[grep(T_epi$Type, pattern = "CD8", value = F), ]
row.names(T_epi) <- 1:dim(T_epi)[1]
T_epi <- T_epi[grep(T_epi$Type, pattern = "Myeloid", value = F), ]
row.names(T_epi) <- 1:dim(T_epi)[1]
T_epi_pvalue <- T_epi %>% dplyr::select(interacting_pair, Types, pvalue) %>% unique()
T_epi_mean <- T_epi %>% dplyr::select(interacting_pair, Types, mean) %>% unique()

TEpi_pvalue <- dcast(T_epi_pvalue, interacting_pair ~ Types) %>% melt()
TEpi_mean <- dcast(T_epi_mean, interacting_pair ~ Types) %>% melt()

merged_T_Epi <- merge(TEpi_pvalue, TEpi_mean, by = c("interacting_pair", "variable"))
setnames(merged_T_Epi, old = c("value.x", "value.y"), new = c("pvalue", "mean"))

merged_T_Epi$mean[is.na(merged_T_Epi$mean)] <- 0
merged_T_Epi$pvalue[is.na(merged_T_Epi$pvalue)] <- 1
row.names(merged_T_Epi) <- 1:dim(merged_T_Epi)[1]

merged_T_Epi <- merged_T_Epi[merged_T_Epi$interacting_pair %in% paired_RL_names, ]
merged_T_Epi$pvalue <- log10(merged_T_Epi$pvalue)
merged_T_Epi$pvalue[is.infinite(merged_T_Epi$pvalue)] <- -1.5
rownames(merged_T_Epi) <- 1:dim(merged_T_Epi)[1]

##---culster the RL using Pvlue; RL: Receptor and Ligation
dct_RL <- dcast(merged_T_Epi[, -4], interacting_pair ~ variable)
row.names(dct_RL) <- dct_RL$interacting_pair
dct_RL <- dct_RL[, -1] %>% as.matrix()
dct_RL[is.na(dct_RL)] <- 0
hdct_RL <- hclust(dist(dct_RL), method = "ward.D2")

##---culster the celltyeps using Pvlue; CT: Cell Type
dct_CT <- dcast(merged_T_Epi[, -4], interacting_pair ~ variable)
row.names(dct_CT) <- dct_CT$interacting_pair
dct_CT <- dct_CT[, -1] %>% as.matrix() %>% t
dct_CT[is.na(dct_CT)] <- 0
hdct_CT <- hclust(dist(dct_CT), method = "ward.D2")

merged_T_Epi$interacting_pair <- factor(merged_T_Epi$interacting_pair, levels = hdct_RL$labels[hdct_RL$order])
merged_T_Epi$variable <- factor(merged_T_Epi$variable, levels = hdct_CT$labels[hdct_CT$order])


png("CD8_VS_Myolid.png", width = 20, height = 20, units = "in", res = 600)
ggplot(merged_T_Epi, aes(x = variable, y = interacting_pair, size = -pvalue)) +
  geom_point(aes(colour = mean)) +
  theme(axis.text.x = element_text(angle = 90, size = 4, vjust = 0.2, hjust = 1),
        axis.text.y = element_text(angle = 0, size = 4, vjust = 0.2, hjust = 1),
        panel.background = element_rect(linetype = "solid", color = "black", size = 1, fill = "white"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.1, 0, 0.1, 0), "in"))  +
  scale_size_continuous(range = c(0, 2), name = "-log10(P-value)") +
  scale_color_gradientn(colours = colorRampPalette(c("#23405b", "#3e73a0", "yellow", "red", "#950000"))(10)) +
  labs(colour = "Mean",
       x = "Cell Type",
       y = "Interacting Pair")
dev.off()


##------epi interacting with Epi------##

T_epi <- data.frame(cellphone_dat %>% dplyr::select(interacting_pair, variable, A_B, V1, V2, cell1, cell2, pvalue, mean)) %>% mutate(Type = paste0(cell1, "_", cell2))
row.names(T_epi) <- 1:dim(T_epi)[1]
T_epi <- T_epi[grep(T_epi$cell1, pattern = "Epi", value = F), ]
row.names(T_epi) <- 1:dim(T_epi)[1]
T_epi <- T_epi[grep(T_epi$cell2, pattern = "Epi", value = F), ]
row.names(T_epi) <- 1:dim(T_epi)[1]

T_epi_pvalue <- T_epi %>% dplyr::select(interacting_pair, variable, pvalue)
paired_RL_names <- T_epi_pvalue %>% dplyr::select(interacting_pair) %>% unlist() %>% unname %>% unique %>% sort()

###.2 extracted the significant paired RL from raw data-------------
T_epi <- data.frame(all_cellphone_dat %>%
                      dplyr::select(interacting_pair, variable, cell1, cell2, Type, pvalue, mean, tissue)) %>%
  mutate(Types = paste0(variable, ".", tissue ))

row.names(T_epi) <- 1:dim(T_epi)[1]
T_epi <- T_epi[grep(T_epi$cell1, pattern = "Epi", value = F), ]
row.names(T_epi) <- 1:dim(T_epi)[1]
T_epi <- T_epi[grep(T_epi$cell2, pattern = "Epi", value = F), ]
row.names(T_epi) <- 1:dim(T_epi)[1]
T_epi_pvalue <- T_epi %>% dplyr::select(interacting_pair, Types, pvalue) %>% unique()
T_epi_mean <- T_epi %>% dplyr::select(interacting_pair, Types, mean) %>% unique()

TEpi_pvalue <- dcast(T_epi_pvalue, interacting_pair ~ Types) %>% melt()
TEpi_mean <- dcast(T_epi_mean, interacting_pair ~ Types) %>% melt()

merged_T_Epi <- merge(TEpi_pvalue, TEpi_mean, by = c("interacting_pair", "variable"))
setnames(merged_T_Epi, old = c("value.x", "value.y"), new = c("pvalue", "mean"))

merged_T_Epi$mean[is.na(merged_T_Epi$mean)] <- 0
merged_T_Epi$pvalue[is.na(merged_T_Epi$pvalue)] <- 1
row.names(merged_T_Epi) <- 1:dim(merged_T_Epi)[1]

merged_T_Epi <- merged_T_Epi[merged_T_Epi$interacting_pair %in% paired_RL_names, ]
merged_T_Epi$pvalue <- log10(merged_T_Epi$pvalue)
merged_T_Epi$pvalue[is.infinite(merged_T_Epi$pvalue)] <- -1.5
rownames(merged_T_Epi) <- 1:dim(merged_T_Epi)[1]

##---culster the RL using Pvlue; RL: Receptor and Ligation
dct_RL <- dcast(merged_T_Epi[, -4], interacting_pair ~ variable)
row.names(dct_RL) <- dct_RL$interacting_pair
dct_RL <- dct_RL[, -1] %>% as.matrix()
dct_RL[is.na(dct_RL)] <- 0
hdct_RL <- hclust(dist(dct_RL), method = "ward.D2")

##---culster the celltyeps using Pvlue; CT: Cell Type
dct_CT <- dcast(merged_T_Epi[, -4], interacting_pair ~ variable)
row.names(dct_CT) <- dct_CT$interacting_pair
dct_CT <- dct_CT[, -1] %>% as.matrix() %>% t
dct_CT[is.na(dct_CT)] <- 0
hdct_CT <- hclust(dist(dct_CT), method = "ward.D2")

merged_T_Epi$interacting_pair <- factor(merged_T_Epi$interacting_pair, levels = hdct_RL$labels[hdct_RL$order])
merged_T_Epi$variable <- factor(merged_T_Epi$variable, levels = hdct_CT$labels[hdct_CT$order])


png("Epi.png", width = 17, height = 14.5, res = 1000, units = "in")
ggplot(merged_T_Epi, aes(x = variable, y = interacting_pair, size = -pvalue)) +
  geom_point(aes(colour = mean)) +
  theme(axis.text.x = element_text(angle = 90, size = 4, vjust = 0.2, hjust = 1),
        axis.text.y = element_text(angle = 0, size = 4, vjust = 0.2, hjust = 1),
        panel.background = element_rect(linetype = "solid", color = "black", size = 1, fill = "white"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.1, 0, 0.1, 0), "in"))  +
  scale_size_continuous(range = c(0, 2), name = "-log10(P-value)") +
  scale_color_gradientn(colours = colorRampPalette(c("#23405b", "#3e73a0", "yellow", "red", "#950000"))(10)) +
  labs(colour = "Mean",
       x = "Cell Type",
       y = "Interacting Pair")
dev.off()


