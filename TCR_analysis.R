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
)

##----------------------------------Step 1. Identify T cell with paired A and B chains.-------------------------------####
number_to_tissue <- read.table("number_corresponding_tissue_TCR.txt", header = F, col.names = c("Number", "Tissue"), stringsAsFactors = F)
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

####----------the result of this step (T_cell_clone_uniq) was used as the input to identify sharing weight.------------######
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$cell_barcode %in% qualited_T_cells, ] 
#########################################################################################################

###------------------remove non T cells TCR----------------------------#####################
T_cell_clone_uniq <- filter(T_cell_clone_uniq, cell_barcode %in% row.names(T_cells_meta.data)) ## 5560 T cells (2402 clones)  with certain clone type and paired A and B chain
T_cell_clone_uniq$annotation <- mapvalues(T_cell_clone_uniq$cell_barcode, 
                                          from = T_cells_meta.data %>% row.names,
                                          to = T_cells_meta.data$annotation)

##--------------------------------------------------Step 2. TCR sharing across organs.--------------------------------####

##remove none clone cells
T_cell_clone_uniq_each <- T_cell_clone_uniq[ !duplicated(T_cell_clone_uniq$"cell_barcode"), ]
duplicated_clones <- names(table(T_cell_clone_uniq_each$raw_clonotype_id)[unname(table(T_cell_clone_uniq_each$raw_clonotype_id))>= 2])
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$raw_clonotype_id %in% duplicated_clones, ] ## 438 cells, 106 clones

T_cell_clone_uniq_sorted <- T_cell_clone_uniq %>% arrange(Tissue, raw_clonotype_id, Barcode) %>% dplyr::select(Barcode, raw_clonotype_id, Tissue)
T_cell_clone_uniq_sorted_h_2_low <- with(T_cell_clone_uniq_sorted, 
                                         lapply(unique(Tissue), function(x, dat){
                                           ord <- match(dat[dat$Tissue == x, "raw_clonotype_id"], 
                                                        names(sort(table(dat[dat$Tissue == x, "raw_clonotype_id"]), decreasing = T)))
                                           each_ordered <- dat[dat$Tissue == x, ][order(ord), ]
                                         }, dat = T_cell_clone_uniq_sorted))

T_cell_clone_sort_h_2_l <- do.call(rbind, T_cell_clone_uniq_sorted_h_2_low)
matrix_cor_expression <- with(T_cell_clone_sort_h_2_l, lapply(1:dim(T_cell_clone_sort_h_2_l)[1], 
                                                              function(x){(raw_clonotype_id %in% raw_clonotype_id[x])+0 }))
Cor_T_matrix <- do.call(cbind, matrix_cor_expression)

row.names(Cor_T_matrix) <- T_cell_clone_sort_h_2_l$Tissue
colnames(Cor_T_matrix) <- T_cell_clone_sort_h_2_l$Tissue

shared_matrix <- unclass(with(T_cell_clone_sort_h_2_l,table(raw_clonotype_id, Tissue)))
shared_matrix[ shared_matrix > 0 ] <- 1

shared_m <- t(shared_matrix) %*% (shared_matrix)
cluster.gr <- igraph::graph_from_adjacency_matrix(shared_m/sum(shared_m), 
                                                  mode="undirected", weighted=TRUE, diag=FALSE)

##------------------3. calculate the sharing weight (TCR tracking analysis across organs)----------------###

TCR <- T_cell_clone_uniq %>% ### the input "T_cell_clone_uniq" was generated from step 1 above. 
  select(c(cell_barcode, raw_clonotype_id, Tissue)) %>%
  mutate(patient = mapvalues(Tissue, from = names(tissus_to_numbers), to = unname(tissus_to_numbers)), loc = Tissue) %>% 
  setnames(old = c("cell_barcode", "raw_clonotype_id", "Tissue"), new = c("Cell_Name", "clone.id", "majorCluster"))

obj <- new("Startrac", TCR, aid = "HCA")
obj <- calIndex(obj)
# tic("pIndex")
obj <- pIndex(obj)

obj@pIndex.tran[is.na(obj@pIndex.tran)] <- 0

migration_across_tissue <- obj@pIndex.tran %>% select(-c(1:2)) 
row.names(migration_across_tissue) <- colnames(migration_across_tissue)
pheatmap::pheatmap(migration_across_tissue)

uppertri <- migration_across_tissue
uppertri[!upper.tri(uppertri)] <- 10
uppertri <- as.matrix(uppertri)
tmp <- as.vector(t(uppertri))
weight <- tmp[!((tmp == 10)|(tmp == 0)) ]

#########----------------------------4. ploting the TCR sharing across organs-----------------------------------------------------##

E(cluster.gr)$weight <- weight

E(cluster.gr)$width <- E(cluster.gr)$weight/6
E(cluster.gr)$width <- 1+E(cluster.gr)$weight/8

V(cluster.gr)$size <- c(T_cell_clone_uniq[ !duplicated(T_cell_clone_uniq$"cell_barcode"), ]$Tissue %>% table()) %>% log2
V(cluster.gr)$frame.color <- unname(tissues_colors)[match(x = names(V(cluster.gr)), table = names(tissues_colors))]

pdf("Tissue_TCR_sharing_network_CD8.pdf", height = 10, width = 10)
set.seed(2)
pt <- plot(cluster.gr,
           edge.width = igraph::E(cluster.gr)$weight*20,
           vertex.color = unname(tissues_colors)[match(x = names(V(cluster.gr)), table = names(tissues_colors))],
           vertex.label.dist = 0,
           vertex.label.color = unname(tissues_colors)[match(x = names(V(cluster.gr)), table = names(tissues_colors))],
           vertex.label.cex	= 1,
           layout = layout_in_circle(cluster.gr))
dev.off()

#####-----------------------------------------5. TCR tracking analysis across subclusters--------------------------------##

TCR <- T_cell_clone_uniq %>% ### the input "T_cell_clone_uniq" was generated from step 1 above.
  select(c(cell_barcode, raw_clonotype_id, Tissue)) %>%
  mutate(patient = Tissue, loc = Tissue,
         Tissue = mapvalues(cell_barcode, from = row.names(T_cells_meta.data), to = T_cells_meta.data$annotation) ) %>% 
  setnames(old = c("cell_barcode", "raw_clonotype_id", "Tissue"), new = c("Cell_Name", "clone.id", "majorCluster"))

obj <- new("Startrac", TCR, aid = "HCA")
obj <- calIndex(obj)
obj <- pIndex(obj)

obj@pIndex.tran[is.na(obj@pIndex.tran)] <- 0

migration_across_tissue <- obj@pIndex.tran %>% select(-c(1:2)) 
row.names(migration_across_tissue) <- colnames(migration_across_tissue)

pdf("heatmap_transition_across_clusters_CD8.pdf", width = 15, height = 15)
pheatmap::pheatmap(migration_across_tissue,
                   # color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                   # breaks = seq(0, 1, length.out = 100),
                   treeheight_row = 8,
                   treeheight_col = 8,
                   border_color = "white",
                   # cellwidth = 20,
                   # cellheight = 5,
                   scale = "none",
                   cluster_rows = T,
                   cluster_cols = T,
                   fontsize_row  = 5,
                   fontsize_col = 5,
                   # annotation_col = annotation_cols,
                   # annotation_colors = ann_colors,
                   show_colnames = T,
                   width = 10,
                   height = 14 
)
dev.off()

###-----------------------------6. distribution of T clonetyps across organs------------------##

TCR_clone_dat <- T_cell_clone_uniq %>% ### the input "T_cell_clone_uniq" was generated from step 1 above.
  select(c(Tissue, cell_barcode, raw_clonotype_id)) %>%
  mutate(cluster = mapvalues(cell_barcode, from = row.names(T_cells_meta.data), to = T_cells_meta.data$annotation) ) %>% unique

group_by_tissue <- by(TCR_clone_dat, TCR_clone_dat$Tissue, FUN = function(x) { `[`(x) })

result <- lapply(1:length(group_by_tissue), FUN = function(x, dat) {
  clone_number <- split(1:(dim(dat[[x]])[1]), dat[[x]]$raw_clonotype_id) %>%
    lapply(FUN = length) %>%
    do.call(what = rbind) %>%
    as.data.frame()
  clone_number[clone_number$V1 >= 3, ] <- 3
  dat[[x]]$clone_numbers <- mapvalues(dat[[x]]$raw_clonotype_id, from = row.names(clone_number), to = clone_number$V1)
  return(dat[[x]])
}, dat = group_by_tissue)

results <- do.call(result, what = rbind)
results$clone_numbers <- factor(results$clone_numbers, levels = c(3, 2, 1))

pdf("TCR_clone_structures_across_tissue_CD8.pdf", width = 15, height = 7)
ggplot(results, aes(x = Tissue,  fill = clone_numbers)) +
  geom_bar(stat = "count") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(expand = c(0.005, 0.1)) +
  scale_fill_discrete(name = "Clonality", labels = c("Clonal", "Duplicated", "Unique"))
dev.off()

###-----------------------------7. distribution of T clonetyps across clusters------------------###

TCR_clone_dat <- T_cell_clone_uniq %>% ### the input "T_cell_clone_uniq" was generated from step 1 above.
  select(c(Tissue, cell_barcode, raw_clonotype_id)) %>%
  mutate(cluster = mapvalues(cell_barcode, from = row.names(T_cells_meta.data), to = T_cells_meta.data$T_subtype) )  %>% unique

group_by_tissue <- by(TCR_clone_dat, TCR_clone_dat$cluster, FUN = function(x) { `[`(x) })

result <- lapply(1:length(group_by_tissue), FUN = function(x, dat) {
  clone_number <- split(1:(dim(dat[[x]])[1]), dat[[x]]$raw_clonotype_id) %>%
    lapply(FUN = length) %>%
    do.call(what = rbind) %>%
    as.data.frame()
  clone_number[clone_number$V1 >= 3, ] <- 3
  dat[[x]]$clone_numbers <- mapvalues(dat[[x]]$raw_clonotype_id, from = row.names(clone_number), to = clone_number$V1)
  return(dat[[x]])
}, dat = group_by_tissue)

results <- do.call(result, what = rbind)
results$clone_numbers <- factor(results$clone_numbers, levels = c(3, 2, 1))
results$cluster <- factor(results$cluster, levels = c('TN_SELL',
                                                      'TN_KLF2',
                                                      'TCM_LEF1',
                                                      'TCM_GADD45B',
                                                      'TEM_INFG',
                                                      'TEM_GZMK',
                                                      'TEM_GIMAP4',
                                                      'TEFF_TRBV4-2',
                                                      'TEFF_MT1E',
                                                      'TEFF_GNLY',
                                                      'TRM_HSPA1A',
                                                      'TRM_H2AFZ',
                                                      'TRM_GZMB',
                                                      'TRM_PRR4',
                                                      'TRM_RGS1',
                                                      'TRM_MT1X',
                                                      'TRM_NABP1',
                                                      'TRM_TYMS',
                                                      'IEL_TMIGD2',
                                                      'IEL_TRBV7-3',
                                                      'MAIT_SLC4A10'
))

pdf("TCR_clone_structures_across_clusters_CD8.pdf", width = 15, height = 7)
ggplot(results, aes(x = cluster,  fill = clone_numbers)) +
  geom_bar(stat = "count") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(expand = c(0.005, 0.1)) +
  scale_fill_discrete(name = "Clonality", labels = c("Clonal", "Duplicated", "Unique"))
dev.off()

#########------------------------8. TCR sharing across organs and clusters -----------------------------------------############

metdata <- T_cells_meta.data ## T_cell_clone_uniq: the input "T_cell_clone_uniq" was generated from step 1 above.
T_cell_clone_uniq <- T_cell_clone_uniq %>% subset(!(annotation %in% grep(T_cell_clone_uniq$annotation, pattern = "TGD", value = T)))
clo_cells <- T_cell_clone_uniq$cell_barcode
clone_metdat <- metdata[row.names(metdata) %in% clo_cells, ]

clone_metdat$raw_clonotype_id <- mapvalues(rownames(clone_metdat), from = T_cell_clone_uniq$cell_barcode, to = T_cell_clone_uniq$raw_clonotype_id) 
size_of_col_type <- sort(table(clone_metdat %>% subset$raw_clonotype_id), decreasing = T)

dat_to_plot <- clone_metdat[ ,c("orig.ident", "annotation", "raw_clonotype_id", "Color_of_tissues")]

clone_tissue <- unclass(t(table(dat_to_plot[, c(1, 3)])))
shared_clone_frequence <- lapply(1:(dim(clone_tissue)[1]), FUN = function(x) {result <- table(as.numeric(clone_tissue[x, ]) > 0)[2]})
shared_clone_frequences <- do.call(rbind, shared_clone_frequence)
row.names(shared_clone_frequences) <- row.names(clone_tissue)
shared_clone_frequences <- data.frame(clones = row.names(shared_clone_frequences), Freq = unname(shared_clone_frequences))
shared_clone_frequences <- shared_clone_frequences[order(shared_clone_frequences$Freq, decreasing = T), ]

dat <- data.frame(table(dat_to_plot[, c(1:3)]))
dat <- dat[!(dat$orig.ident == "Testis"), ]

col_fre_larger_than_x <- shared_clone_frequences[shared_clone_frequences$Freq >= 8, "clones"]
dat1 <- dat[dat$raw_clonotype_id %in% col_fre_larger_than_x, ]

dat1$cols_of_tissue <- mapvalues(as.character(dat1$orig.ident),
                                 from = unique(as.character(clone_metdat$orig.ident)),
                                 to = unique(as.character(clone_metdat$Color_of_tissues) ))

dat1$raw_clonotype_id <- factor(as.character(dat1$raw_clonotype_id),
                                levels = names(size_of_col_type)[sort(match(x = unique(as.character(dat1$raw_clonotype_id)),
                                                                            table = names(size_of_col_type)), decreasing = F)])
###---order the tissues
order_of_tissues <- dat1 %>% subset(Freq > 0) %>% group_by(orig.ident) %>% dplyr::summarise(., total = sum(Freq))
order_of_tissues <- order_of_tissues[order(order_of_tissues$total, decreasing = T), ] %>% `[`(1)  %>%  unlist %>% as.character
dat1 <- dat1 %>% subset(orig.ident %in% order_of_tissues)
dat1$orig.ident <- factor(dat1$orig.ident %>% as.character, levels = order_of_tissues)

###---order the Clusters
order_of_clusters <- dat1 %>% subset(Freq > 0) %>% group_by(annotation) %>% dplyr::summarise(., total = sum(Freq))
order_of_clusters <- order_of_clusters[order(order_of_clusters$total, decreasing = T), ] %>% `[`(1) %>%  unlist %>% as.character
dat1 <- dat1 %>% subset(annotation %in% order_of_clusters)
dat1$annotation <- factor(dat1$annotation %>% as.character, levels = order_of_clusters) 


png("TCR_share_across_organs_and_cell_types_CD8_with_legend.png", width = 20, height = 20, res = 500, units = "in")
ggplot(data = dat1[, ],
       aes(axis1 = raw_clonotype_id, axis2 = orig.ident, axis3 = annotation, y = Freq)) +
  scale_x_discrete(limits = c("raw_clonotype_id", "orig.ident", "annotation"), expand = c(.1, .05)) +
  xlab("Demographic") + 
  geom_alluvium(aes(fill = orig.ident)) +
  geom_stratum() +             
  geom_text(stat = "stratum", label.strata = T, check_overlap = T, size = 6) +
  theme_minimal() +
  # scale_fill_manual(values = as.character(unique(dat1$cols_to_use)[match(x = levels(dat1$T_subtype), table = unique(dat1$T_subtype))])) +
  scale_fill_manual(values = c(tissues_colors[match(levels(dat1$orig.ident), tissues_colors %>% names)]) %>% unname) +
  ggtitle("TCR_share_across_organs_and_cell_types") +
  theme(legend.position = "none")
dev.off()
