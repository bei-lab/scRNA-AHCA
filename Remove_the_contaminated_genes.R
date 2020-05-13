library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggthemes)
library(cowplot)
library(data.table)
library(parallel)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

## Take fibroblast and smooth muscle cell as an example.
## We remove the contaminated gene in each tissue before comparing the difference between different tissues.
## 1. 1) load the singlet data identfied in each tissue. 
##    2) select "the Fibroblast" and "Smooth Muscle Cell" clusters. 
##    3) Identfy the differently expressed genes between tissue.
## 2. calculate the top 1% highly expressed genes in each tissue.
## 3. identify the low-frequency genes in each tissue.
## 4. remove the contaminated gene and then normalize the data, run PCA, find the cluster and run tSNE
## 5. remove the contaminated clusters.


### load the HCA_data_singlet
load("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/4-Seurat_cluster/HCA_all_data_88638_cells.RData")

HCA_all_data_singlet <- HCA_data

##############--------------------------scale and normalize the data------------------------------------------#####################
Idents(HCA_data) <- factor(as.character(HCA_data$orig.ident))
celltypes <- c("Fib|Smo" )

for (celltype in celltypes) {
  subset_cells <- subset(x = HCA_data, cells = names(grep(celltype, HCA_data$Cell_types_in_tissue, value = T)))
  # subset_cells <- subset(x = subset_cells,
  #                        subset = nFeature_RNA >= 500 & nFeature_RNA < 4000 & percent.mito <= 0.25 & nCount_RNA >= 1000 & nCount_RNA <= 10000)
  print(table(subset_cells@active.ident))
  subset_cells <- NormalizeData(object = subset_cells, normalization.method = "LogNormalize", scale.factor = 1e4)
  subset_cells <- ScaleData(object = subset_cells, features = rownames(x = subset_cells))#, vars.to.regress = c("nCount_RNA", "percent.mito"))
  
  ###############-----------------------------1. Find all the marker genes of same cell types across the tissues-------------------####################
  ###subsetclusters markers
  result_subset_cells <- FindMarkers_parallel(object = subset_cells, mc.cores = 36)
  subset_cells_markers <- do.call(rbind, result_subset_cells)
  subset_cells_markers$gene <- unlist(mapply(rownames, result_subset_cells))
  subset_cells_markers$cluster <- rep(names(table(subset_cells@active.ident))[table(subset_cells@active.ident) > 30], 
                                      times = mapply(dim, result_subset_cells, SIMPLIFY = TRUE)[1,])
  
  subset_cells_markers %>% TOP_N(200) -> top200
  reordered_subset_cells_markers <- subset_cells_markers %>% TOP_N(4000)
  
  write.table(top200, paste0("HCA_top200_max300cells_", make.names(celltype), ".csv"), sep = ",", row.names = T, quote = F, col.names = NA)
  write.table(reordered_subset_cells_markers, paste0("HCA_allFC_max300cells_", make.names(celltype), ".csv"), sep = ",", row.names = T, quote = F, col.names = NA)
  rm(subset_cells_markers, top200, celltype, result_subset_cells)
  gc()
}

####------------------------------------2. calculate top 1% high expression genes in each organ-----------------------#############################################
##calculate the top 1% high expression genes
genes_of_2_percent_expression <- data.frame(Gene_names = NULL,
                                            UMI_sum = NULL,
                                            Tissues = NULL)

for (tissue in HCA_all_data_singlet$orig.ident %>% as.character %>% unique) {
  counts <- apply(as.matrix(GetAssayData(subset(HCA_data, idents = tissue) ,assay = "RNA", slot = "counts")),
                  MARGIN = 1,
                  FUN = function(x){sum(x, na.rm = TRUE)})
  counts <- sort(counts[counts > 0 ], decreasing = TRUE)
  top_2_percent_genes <- counts[counts > quantile(counts[counts > 1], seq(0, 1, 0.01), na.rm = TRUE)[99]]
  genes_of_2_percent_expression <- rbind(genes_of_2_percent_expression,
                                         data.frame(Gene_names = names(top_2_percent_genes),
                                                    UMI_sum = unname(top_2_percent_genes),
                                                    Tissues = rep(tissue, times = length(top_2_percent_genes))))
  print(length(top_2_percent_genes))
  print(c(tissue, quantile(counts[counts>1], seq(0, 1, 0.01), na.rm = T )))
  rm(counts, top_2_percent_genes)
  gc()
}
write.table(genes_of_2_percent_expression,
            file = "genes_of_2_percent_expression.txt",
            sep = ",",
            quote = F,
            row.names = T,
            col.names = NA)
####------------------------------------3. identify the low-frequency genes in each tissue-----------------######################################################
genes_of_2_percent_expression <- read.table("genes_of_2_percent_expression.txt",
                                            sep = ",",
                                            stringsAsFactors = F,
                                            header = T,
                                            row.names = 1)

for (celltype in celltypes) {
  DEG <- read.table(paste0("HCA_allFC_max300cells_", make.names(celltype), ".csv"),
                    sep = ",",
                    stringsAsFactors = F,
                    row.names = 1,
                    header = T)
  ###
  pct.2_less_than_0.05 <- unique(as.character(with(DEG, DEG[avg_logFC > 0 & pct.2 < 0.05, "gene"])))
  High_expression_genes <- unique(with(genes_of_2_percent_expression, genes_of_2_percent_expression[Tissues %in% DEG$cluster, "Gene_names"]))
  genes_to_removed <- intersect(pct.2_less_than_0.05, High_expression_genes)
  specific_genes_in_each_tissue <- with(DEG,
                                        data.frame(removed_or_not = (avg_logFC > 0 & pct.2 < 0.05) & (gene %in% High_expression_genes),
                                                   Tissues = cluster,
                                                   Genes = gene))
  write.table(specific_genes_in_each_tissue,
              paste0(make.names(celltype), "_specific_genes_in_each_tissue.txt"),
              sep = "\t",
              quote = FALSE,
              row.names = T,
              col.names = NA)
  print(table(specific_genes_in_each_tissue[, 1], specific_genes_in_each_tissue[, 2]))
  print(genes_to_removed)
  print(c(celltype, length(genes_to_removed)))
  write.table(data.frame(genes_to_removed = genes_to_removed), paste0(make.names(celltype), "_genes_to_removed.txt"),
              sep = "\t", row.names = T, col.names = NA, quote = FALSE)
  rm(celltype, genes_to_removed, DEG, pct.2_less_than_0.05, High_expression_genes, specific_genes_in_each_tissue)
  gc()
}

####------------------------------------4. remove the blacklist gene and then normalize the data, run PCA, find the cluster and run tSNE----------------------########################
############################--------------first round------------------------------------------###################################
##############----------------------------scale and normalize the data------------------------------------------#####################
celltypes <- c("Fib|Smo")

for (celltype in celltypes) {
  genes <- as.character(read.table(paste0(make.names(celltype), "_genes_to_removed.txt"),
                                   header = T,
                                   sep = "\t",
                                   stringsAsFactors = F,
                                   row.names = 1)[, 1])
  subset_cells <- subset(x = HCA_data, cells = names(grep(celltype, HCA_data$Cell_types_in_tissue, value = T)))
  subset_cells <- subset_cells[!(c(subset_cells %>% row.names()) %in% genes), ]
  # subset_cells <- subset(x = subset_cells,
  #                        subset = nFeature_RNA >= 500 & nFeature_RNA < 4000 & percent.mito <= 0.25 & nCount_RNA >= 1000 & nCount_RNA <= 10000)
  # print(table(subset_cells@active.ident))
  subset_cells <- NormalizeData(object = subset_cells, normalization.method = "LogNormalize", scale.factor = 1e4)
  subset_cells <- ScaleData(object = subset_cells, features = rownames(x = subset_cells))#, vars.to.regress = c("nCount_RNA", "percent.mito"))
  subset_cells <- FindVariableFeatures(object = subset_cells, selection.method = 'mean.var.plot', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf))
  subset_cells <- RunPCA(object = subset_cells, features = VariableFeatures(object = subset_cells), verbose = FALSE, npcs = 50)
  
  dim.use <- 25
  res.use <- 1.2
  
  subset_cells <- FindNeighbors(object = subset_cells, dims = 1:dim.use)
  subset_cells <- FindClusters(object = subset_cells, resolution = res.use)
  
  ###Run the TSNE
  subset_cells <- RunTSNE(object = subset_cells, dims = 1:dim.use)
  
  ### plot the tSNE
  png(paste0(make.names(celltype), "_tSNE_Round_1_", dim.use, "_", res.use, ".png"),
      width = 15, height = 15, units = "in", res = 300)
  p2 <- DimPlot(object = subset_cells, reduction = 'tsne', label = TRUE, pt.size = 1,) + NoLegend()
  print(p2)
  dev.off()
  
  png(paste0(make.names(celltype), "tSNE_Round_1_groupby_", dim.use, "_", res.use, ".png"),
      width = 15, height = 15, units = "in", res = 300)
  p3 <- DimPlot(object = subset_cells, reduction = 'tsne', group.by = "orig.ident", pt.size = 1,
                cols = unique(subset_cells$Color_of_tissues)[match(levels(subset_cells$orig.ident)[table(subset_cells@meta.data$orig.ident) > 0], unique(as.character(subset_cells$orig.ident)) )]) + NoLegend()
  print(p3)
  dev.off()
  
  png(paste0(make.names(celltype), "tSNE_Round_1_groupby_with_legend_", dim.use, "_", res.use, ".png"),
      width = 15, height = 15, units = "in", res = 300)
  p4 <- DimPlot(object = subset_cells, reduction = 'tsne', group.by = "orig.ident", pt.size = 1, 
                cols = unique(subset_cells$Color_of_tissues)[match(levels(subset_cells$orig.ident)[table(subset_cells@meta.data$orig.ident) > 0], unique(as.character(subset_cells$orig.ident)) )]) 
  print(p4)
  dev.off()
  
  ###-----------------------------Find all the marker genes-------------------------------------------
  ###subsetclusters markers
  ###find all the markers
  result <- FindMarkers_parallel(object = subset_cells, mc.cores = 36)
  all_markers <- do.call(rbind, result)
  all_markers$gene <- unlist(mapply(rownames, result))
  all_markers$cluster <- rep(levels(subset_cells@active.ident), times = mapply(dim, result, SIMPLIFY = TRUE)[1,])
  subset_cells.markers <- all_markers
  
  assign(make.names(celltype), subset_cells)
  save(list = make.names(celltype),
       file = paste0(
                     make.names(celltype),
                     "_before_filtering.RData"))
  
  subset_cells.markers %>% TOP_N(50) -> top50
  subset_cells.markers <- subset_cells.markers %>% TOP_N(5000)
  
  write.table(top50,
              file = paste0(
                            make.names(celltype), "_top50_", dim.use, "_", res.use, ".csv"),
              sep = ",",
              row.names = T,
              quote = F)
  
  write.table(subset_cells.markers,
              file = paste0(
                            make.names(celltype), "_all_DEGs_", dim.use, "_", res.use, ".csv"),
              sep = ",",
              row.names = T,
              quote = F)
  
  ### plot the heat map
  subset_cells.markers %>% TOP_N(10) -> top10
  png(paste0(
             make.names(celltype), "_", dim.use,"_",res.use, "_heatmap.png"),
      width = 15, height = 15, units = "in", res = 300)
  
  p5 <- Fixed_DoHeat_map(object = subset_cells, features = top10$gene, size = 2) + NoLegend() +
    theme(axis.text.x = element_text(size = 2),
          axis.text.y = element_text(size = 2),
          axis.title.x = element_text(colour = "red", size = 2),
          axis.title.y = element_text(colour = "black", size = 2))
  print(p5)
  dev.off()
}

###------------------------------------------use the feature plot and vlnplot to identify the cell types------------------------------###########################

png("vlnplot_of_FibSmo.png", height = 60, width = 20, units = "in", res = 400)
VlnPlot(subset_cells, features = c("CD8A", "CD4", "CD3E", "CD3D", "FCGR3A", "GNLY", "MMP2", "ACTA2", "CD14", "EPCAM", "PECAM1", "PTRRC"), pt.size = 0.1,
        # cols = c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10))[-8], 
        ncol = 1)  
dev.off()

png("feature_plot_of_FibSmo.png", height = 60, width = 30, units = "in", res = 400)
FeaturePlot(subset_cells, features = c("CD8A", "CD4", "CD3E", "CD3D", "FCGR3A", "GNLY", "MMP2", "ACTA2", "CD14", "EPCAM", "PECAM1", "PTRRC"), pt.size = 0.5,
            # cols = c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10))[-8], 
            ncol = 3)  
dev.off()

####-----------------------------------------5. removed non-FibSmo cell clusters---------------###################################
celltypes <- c("Fib|Smo")

for (celltype in celltypes) {
  all_certain_cells <- read.table("Fib.Smo_top50_25_1.2.csv",
                                  header = T,
                                  row.names = 1,
                                  stringsAsFactors = F,
                                  sep = ",")
  subset_cells$YES_OR_NO <- mapvalues(as.character(subset_cells@meta.data$seurat_clusters), from = as.character(all_certain_cells$cluster), to = as.character(all_certain_cells$cell.type))
  write.table(data.frame(table(subset_cells$YES_OR_NO)), paste0(make.names(celltype), "_cells_removed.txt"), sep = "\t", quote = F, row.names = T)
  subset_cells <- subset(x = subset_cells,
                         cells = row.names(subset_cells@meta.data)[as.character(subset_cells$YES_OR_NO) %in% "YES"])
  
  subset_cells <- NormalizeData(object = subset_cells, normalization.method = "LogNormalize", scale.factor = 1e4)
  subset_cells <- ScaleData(object = subset_cells, features = rownames(x = subset_cells))#, vars.to.regress = c("nCount_RNA", "percent.mito"))
  subset_cells <- FindVariableFeatures(object = subset_cells, selection.method = 'mean.var.plot', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf))
  subset_cells <- RunPCA(object = subset_cells, features = VariableFeatures(object = subset_cells), verbose = FALSE, npcs = 100)
  
  dim.use <- 20
  res.use <- 1
  
  subset_cells <- FindNeighbors(object = subset_cells, dims = 1:dim.use)
  subset_cells <- FindClusters(object = subset_cells, resolution = res.use)
  
  ### Run the TSNE
  subset_cells <- RunTSNE(object = subset_cells, dims = 1:dim.use)
  
  ### plot the tSNE
  png(paste0(
             make.names(celltype), "_tSNE_filtered_", dim.use, "_", res.use, ".png"),
      width = 15, height = 15, units = "in", res = 300)
  p2 <- DimPlot(object = subset_cells, reduction = 'tsne', label = TRUE, pt.size = 1, label.size = 5) + NoLegend()
  print(p2)
  dev.off()
  
  png(paste0(
             make.names(celltype), "tSNE_filtered_groupby_", dim.use, "_", res.use, ".png"),
      width = 15, height = 15, units = "in", res = 300)
  p3 <- DimPlot(object = subset_cells, reduction = 'tsne', group.by = "orig.ident", pt.size = 1, 
                cols = unique(subset_cells$Color_of_tissues)[match(levels(subset_cells$orig.ident)[table(subset_cells@meta.data$orig.ident) > 0], unique(as.character(subset_cells$orig.ident)) )]) + NoLegend()
  print(p3)
  dev.off()
  
  png(paste0(
             make.names(celltype), "tSNE_filtered_groupby_with_legend_", dim.use, "_", res.use, ".png"),
      width = 15, height = 15, units = "in", res = 300)
  p4 <- DimPlot(object = subset_cells, reduction = 'tsne', group.by = "orig.ident", pt.size = 1,
                cols = unique(subset_cells$Color_of_tissues)[match(levels(subset_cells$orig.ident)[table(subset_cells@meta.data$orig.ident) > 0], unique(as.character(subset_cells$orig.ident)) )])
  print(p4)
  dev.off()
  
  ###-----------------------------Find all the marker genes-------------------------------------------##########################
  ###subsetclusters markers
  ### find all the markers
  result <- FindMarkers_parallel(object = subset_cells, mc.cores = 36)
  all_markers <- do.call(rbind, result)
  all_markers$gene <- unlist(mapply(rownames, result))
  all_markers$cluster <- rep(levels(subset_cells@active.ident), times = mapply(dim, result, SIMPLIFY = TRUE)[1,])
  subset_cells.markers <- all_markers
  
  assign(paste0(make.names(celltype), "_filtered"), subset_cells)
  save(list = paste0(make.names(celltype), "_filtered"),
       file = paste0(make.names(celltype), "_filtered.RData"))
  
  subset_cells.markers %>% TOP_N(50) -> top50
  subset_cells.markers <-  subset_cells.markers %>% TOP_N(5000)
  write.table(top50,
              file = paste0(make.names(celltype), "_top50_", dim.use, "_", res.use, "_filtered.csv"),
              sep = ",",
              row.names = T,
              quote = F)
  
  write.table(subset_cells.markers,
              file = paste0(make.names(celltype), "_all_DEGs_", dim.use, "_", res.use, "_filtered.csv"),
              sep = ",",
              row.names = T,
              quote = F)
  
  subset_cells.markers %>% TOP_N(6000, fc.threshold = log(1.5, base = exp(1))) -> top_fc_lth_1.5 ## genes with fc > 1.5 
  png(paste0(make.names(celltype), "_", dim.use,"_",res.use, "_filtered_heatmap.png"),
      width = 15,
      height = 15,
      res = 300, units = "in")
  p5 <- DoHeatmap(object = subset_cells, features = top_fc_lth_1.5$gene, size = 2) + NoLegend() +
    theme(axis.text.x = element_text(size = 2),
          axis.text.y = element_text(size = 2),
          axis.title.x = element_text(colour = "red", size = 2),
          axis.title.y = element_text(colour = "black", size = 2))
  print(p5)
  dev.off()
  
}

