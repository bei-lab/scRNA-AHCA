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
## We remove the comtaminated gene in each tissue before comparing the difference between different tissues.
## 1. load the sinlet data identfied in each tissue.
## 2. select "the Fibroblast" and "Smooth Muscle Cell" clusters.


### load the HCA_data_singlet
load("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/4-Seurat_cluster/HCA_all_data_88638_cells.RData")

HCA_all_data_singlet <- HCA_data

##############--------------------------1. scale and normalize the data------------------------------------------#####################
Idents(HCA_data) <- factor(as.character(HCA_data$orig.ident))
celltypes <- c("Fib|Smo" )

for (celltype in celltypes) {
  subset_cells <- subset(x = HCA_data, cells = names(grep(celltype, HCA_data$Cell_types_in_tissue, value = T)))
  # subset_cells <- subset(x = subset_cells,
  #                        subset = nFeature_RNA >= 500 & nFeature_RNA < 4000 & percent.mito <= 0.25 & nCount_RNA >= 1000 & nCount_RNA <= 10000)
  print(table(subset_cells@active.ident))
  subset_cells <- NormalizeData(object = subset_cells, normalization.method = "LogNormalize", scale.factor = 1e4)
  subset_cells <- ScaleData(object = subset_cells, features = rownames(x = subset_cells))#, vars.to.regress = c("nCount_RNA", "percent.mito"))
  
  ###############-----------------------------2. Find all the marker genes of same cell types across the tissues-------------------####################
  ###subsetclusters markers
  result_subset_cells <- mclapply(names(table(subset_cells@active.ident))[table(subset_cells@active.ident) > 30],
                                  FUN =  function(x) {FindMarkers(subset_cells,
                                                                  ident.1 = x, 
                                                                  ident.2 = NULL,
                                                                  logfc.threshold = 0.25,
                                                                  only.pos = TRUE,
                                                                  max.cells.per.ident = 300)},
                                  mc.cores = 36)
  
  
  # RESULT <- result_subset_cells
  roundN <- 1
  while(any(mapply(length, result_subset_cells, SIMPLIFY = TRUE)!=5)){
    if(any(mapply(length, result_subset_cells, SIMPLIFY = TRUE)!=5)){
      recalculate_clusters <- which(mapply(length, result_subset_cells, SIMPLIFY = TRUE)!=5)-1
      print(recalculate_clusters)
      result1 <- mclapply(recalculate_clusters,
                          FUN =  function(x) {FindMarkers(subset_cells,
                                                          ident.1 = x,
                                                          ident.2 = NULL,
                                                          logfc.threshold = 0.25,
                                                          only.pos = TRUE,
                                                          max.cells.per.ident = 300)},
                          mc.cores = 35)
    }
    print(roundN + 1)
    for(i in 1:length(recalculate_clusters)){
      result_subset_cells[c(recalculate_clusters+1)[i]] <- result1[i]
    }
  }
  
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

####------------------------------------3. calculate top 2% high expression genes in each organ-----------------------#############################################
##loading the data
# laod("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/5-HCA-all/HCA_data_dims30_res1_only_singlet/HCA_all_data_singlet.RData")
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
####------------------------------------4. the low-frequency gene in each tissue-----------------######################################################
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

####------------------------------------5. remove the blacklist gene and then normalize the data, run PCA, find the cluster and run tSNE----------------------########################
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
  
  if(celltype %in% c("T |Immune", "Fib|Smo", "B Cells|Plasma Cells" )){
    dim.use <- 25
    res.use <- 1.2
  } else if(celltype %in% c("Monocytes|Macrophages|Langerhans|Myeloid", "Endothelial")){
    dim.use <- 20
    res.use <- 1
  } else if(celltype == "B Cells"){
    dim.use <- 15
    res.use <- 1
  } else {
    dim.use <- 25
    res.use <- 1 
  }
  
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
  date()
  system.time({
    date()
    result <- mclapply(as.numeric(levels(subset_cells@active.ident)),
                       FUN =  function(x) {FindMarkers(subset_cells, ident.1 = x, ident.2 = NULL)},
                       mc.cores = 35)
    date()
  })
  date()
  RESULT <- result
  
  roundN <- 1
  while(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
    if(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
      recalculate_clusters <- which(mapply(length, result, SIMPLIFY = TRUE)!=5)-1
      print(recalculate_clusters)
      result1 <- mclapply(recalculate_clusters,
                          FUN =  function(x) {FindMarkers(subset_cells, ident.1 = x, ident.2 = NULL)},
                          mc.cores = 35)
    }
    print(roundN + 1)
    for(i in 1:length(recalculate_clusters)){
      result[c(recalculate_clusters+1)[i]] <- result1[i]
    }
  }
  
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

###------------------------------------------use the feature plot to identify the cell types------------------------------###########################

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


####-----------------------------------------removed non-FibSmo cell clusters---------------###################################
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
  
  if(celltype %in% c("Fib|Smo")){
    dim.use <- 20
    res.use <- 1
  } else if(celltype %in% c("Monocytes ", "Macrophages ", "Smooth Muscle Cells")){
    dim.use <- 15
    res.use <- 1
  } else if(celltype == "B Cells"){
    dim.use <- 15
    res.use <- 1
  } else {
    dim.use <- 20
    res.use <- 1 
  }
  
  subset_cells <- FindNeighbors(object = subset_cells, dims = 1:dim.use)
  subset_cells <- FindClusters(object = subset_cells, resolution = res.use)
  
  ###4. Run the TSNE
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
  ###5. find all the markers
  date()
  system.time({
    date()
    result <- mclapply(as.numeric(levels(subset_cells@active.ident)),
                       FUN =  function(x) {FindMarkers(subset_cells, ident.1 = x, ident.2 = NULL)},
                       mc.cores = 35)
    date()
  })
  date()
  RESULT <- result
  
  roundN <- 1
  while(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
    if(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
      recalculate_clusters <- which(mapply(length, result, SIMPLIFY = TRUE)!=5)-1
      print(recalculate_clusters)
      result1 <- mclapply(recalculate_clusters,
                          FUN =  function(x) {FindMarkers(subset_cells, ident.1 = x, ident.2 = NULL)},
                          mc.cores = 35)
    }
    print(roundN + 1)
    for(i in 1:length(recalculate_clusters)){
      result[c(recalculate_clusters+1)[i]] <- result1[i]
    }
  }
  
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

###############################--------------------------------------organ-specific signature of same cell types----------------------------###################### 


cell_types <- c("B.Cells.Plasma.Cells_filtered.RData",
                "Endothelial_filtered.RData",
                "Fibroblasts.Smooth.Muscle.Cells_filtered.RData",
                "Monocytes.Macrophages.Langerhans_filtered.RData")

all_T <- read.table("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/7-HCA_final_results/Figure_2/T_Cells_cluster/T._top50_20_1.csv",
                    header = T,
                    row.names = 1,
                    stringsAsFactors = F,
                    sep = ",")

celltypes <- c("Monocytes ", "Macrophages |Langerhans|Kupffer", "T ", "Fibroblasts", "Smooth Muscle Cells",
               "Endothelial Cells", "B Cells" , "Schwann Cells", "Plasma Cells", "Macrophages/Monocytes")

celltype <- "T "
load("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/7-HCA_final_results/Figure_2/T..RData")
subset_cells <- T.
subset_cells$T_NK_or_not <- mapvalues(as.character(T.@meta.data$RNA_snn_res.1), from = as.character(all_T$cluster), to = as.character(all_T$Cluster))
subset_cells <- subset(x = subset_cells,
                       cells = row.names(subset_cells@meta.data)[as.character(subset_cells$T_NK_or_not) %in% "T Cells"])
Idents(subset_cells) <- subset_cells$orig.ident

subset_cells <- NormalizeData(object = subset_cells, normalization.method = "LogNormalize", scale.factor = 1e4)
subset_cells <- ScaleData(object = subset_cells, features = rownames(x = subset_cells))#, vars.to.regress = c("nCount_RNA", "percent.mito"))

###------------------------------------------Find all the marker genes between organs-------------------------------------------##########################
###subsetclusters markers
###5. find all the markers
result_subset_cells <- mclapply(names(table(subset_cells@active.ident))[table(subset_cells@active.ident) > 0],
                                FUN =  function(x) {FindMarkers(subset_cells,
                                                                ident.1 = x, 
                                                                ident.2 = NULL,
                                                                logfc.threshold = 0.25,
                                                                only.pos = FALSE,
                                                                max.cells.per.ident = Inf)},
                                mc.cores = 36)



subset_cells_markers <- do.call(rbind, result_subset_cells)
subset_cells_markers$gene <- unlist(mapply(rownames, result_subset_cells))
subset_cells_markers$cluster <- rep(names(table(subset_cells@active.ident))[table(subset_cells@active.ident) > 0], 
                                    times = mapply(dim, result_subset_cells, SIMPLIFY = TRUE)[1,])

subset_cells_markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC) -> top200

write.table(top200, paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/7-HCA_final_results/Figure_2/T_Cells_cluster/", "tissue_specific_genes_top200_removed_high_exp_genes_", make.names(celltype), ".csv"), sep = ",", row.names = T, quote = F, col.names = NA)
write.table(subset_cells_markers, paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/7-HCA_final_results/Figure_2/T_Cells_cluster/", "tissue_specific_genes_allFC_removed_high_exp_genes_", make.names(celltype), ".csv"), sep = ",", row.names = T, quote = F, col.names = NA)
# rm(subset_cells_markers, top200)
gc()

# assign(make.names(celltype), subset_cells)
# save(list = make.names(celltype), file = paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/7-HCA_final_results/Figure_2/T_Cells_cluster/", make.names(celltype), "_filtered.RData"))

subset_cells_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) -> top50
write.table(top50, paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/7-HCA_final_results/Figure_2/T_Cells_cluster/", make.names(celltype), "_top50_", dim.use, "_", res.use, "_organ_specfic.csv"),
            sep = ",",
            row.names = T,
            quote = F)
write.table(subset_cells_markers, paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/7-HCA_final_results/Figure_2/T_Cells_cluster/", make.names(celltype), "_all_DEGs_", dim.use, "_", res.use, "_organ_specfic.csv"),
            sep = ",",
            row.names = T,
            quote = F)

###7. plot the heat map
subset_cells_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
png(paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/7-HCA_final_results/Figure_2/T_Cells_cluster/", make.names(celltype), "_", dim.use,"_",res.use, "_organ_specfic_heatmap.png"),
    width = 15,
    height = 15,
    res = 300, units = "in")
p5 <- DoHeatmap(object = subset_cells, features = top50$gene, size = 2) + NoLegend() +
  theme(axis.text.x = element_text(size = 2),
        axis.text.y = element_text(size = 2),
        axis.title.x = element_text(colour = "red", size = 2),
        axis.title.y = element_text(colour = "black", size = 2))
print(p5)
dev.off()

#####################------------------------organ specific pheatmap----------------------------------------------------#################################
##calculate the average expression
cluster.averages <- AverageExpression(object = subset_cells)
cluster.averages <- as.data.frame(cluster.averages$RNA)

DEG_across_organ <- read.table(paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/7-HCA_final_results/Figure_2/T_Cells_cluster/",
                                      "tissue_specific_genes_top200_removed_high_exp_genes_", make.names(celltype), ".csv"),
                               sep = ",",
                               header = T,
                               row.names = 1,
                               stringsAsFactors = F)

label <- c("XCL1", "LGALS1", "XCL2", "NR4A1", "SQSTM1", "IFNG", "NR4A2", "TMIGD2", "S1PR1",
           "SELL", "CCR7", "LEF1", "TCF7", "IFNGR2", "TNF", "GADD45B", "KLRC1", "CAPG", "KLRB1",
           "STMN1", "TYMS", "MKI67", "LMNA", "MYADM", "NR4A3", "ITGAE", "CD69", "IKZF2", "CD44",
           "CTLA4", "GPR15","TRDV1", "IL7R", "KIR2DL4", "RORA", "TRDV2", "TRGV9", "RGS1", "ENTPD1")

# DEG_to_plot <- cluster.averages[unique(DEG_across_organ[(DEG_across_organ$avg_logFC >= 0.5 | DEG_across_organ$avg_logFC <= -0.6), "gene"]), ]
DEG_to_plot_1 <- cluster.averages[label, ]
DEG_to_plot <- apply(DEG_to_plot_1, MARGIN = 1, scale, center = TRUE)
DEG_to_plot[DEG_to_plot > 2] <- 2
DEG_to_plot[DEG_to_plot < -2] <- -2
row.names(DEG_to_plot) <- colnames(DEG_to_plot_1)
DEG_to_plot <- as.data.frame(t(DEG_to_plot))

# png("T_organ_specific_gene.png", height = 15, width = 15, res = 300, units = "in")
# pheatmap::pheatmap(t(DEG_to_plot),
#                    # color = colorRampPalette(c("green", "black", "red"))(50),
#                    # color = colorRampPalette(brewer.pal(8, "PiYG"))(25),
#                    # color = colorRampPalette(brewer.pal(8, "Blues"))(25),
#                    scale = "none",
#                    cluster_rows = T,
#                    cluster_cols = T,
#                    fontsize_row = 3)
# dev.off()
# 
# 
# x <- unlist(tissues_colors)
# # 
# # col_ano <-  data.frame(type = tissues_colors[c(1:17,19)])
# ha <-  HeatmapAnnotation(df = col_ano, col = list(type = "red"))
# draw(ha, 1:18)

label <- c("XCL1", "LGALS1", "XCL2", "NR4A1", "SQSTM1", "IFNG", "NR4A2", "TMIGD2", "S1PR1",
           "SELL", "CCR7", "LEF1", "TCF7", "IFNGR2", "TNF", "GADD45B", "KLRC1", "CAPG", "KLRB1",
           "STMN1", "TYMS", "MKI67", "LMNA", "MYADM", "NR4A3", "ITGAE", "CD69", "IKZF2", "CD44",
           "CTLA4", "GPR15","TRDV1", "IL7R", "KIR2DL4", "RORA", "TRDV2", "TRGV9", "RGS1", "ENTPD1")

subset_position <- match(label, label)

png(paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/7-HCA_final_results/Figure_2/T_Cells_cluster/T_organ_specific_gene.png"), height = 15, width = 15, res = 300, units = "in")
Heatmap(DEG_to_plot,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_row_dend = FALSE,
        show_column_dend = T) #+ 
# rowAnnotation(link =
#                 row_anno_link(at = subset_position, labels = label),
#               width = unit(0.5, "cm") + max_text_width(labels))
dev.off()

####################-------------------------pathway analysis----------------------------------------############################
library(clusterProfiler)
library(AnnotationHub)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(cowplot)
library(Seurat)
library(RColorBrewer)
library(ggthemes)
library(Seurattopython)
data("geneList")
library(monocle)

genes <- read.table("tissue_specific_genes_top200_removed_high_exp_genes_T..csv", sep = ",", header = T, row.names = 1)

gene.df <- bitr(with(genes, genes[cluster == "Marrow" & avg_logFC > 0, "gene"]), fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene          = gene.df$ENTREZID,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)
ego <- setReadable(ego, OrgDb = org.Hs.eg.db)
dotplot(ego, showCategory=30, font.size = 8)# +
#   theme(panel.grid.minor = element_blank(),
#         line = element_blank())

kk <- enrichKEGG(gene.df$ENTREZID, organism="hsa", pvalueCutoff=1, pAdjustMethod="BH", 
                 qvalueCutoff=1)

dotplot(kk, showCategory=30, font.size = 8)


# load("Monocytes._cells_data.RData")
load("T._cells_data.RData")
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
                        Trachea = '#d62790',
                        x1 = '#00cdaa',
                        x2 = '#aacdaa',
                        x3 = '#aaaa00',
                        x4 = '#ff00ff',
                        x5 = '#00acff',
                        x6 = '#ccacff'
)
celltypes <- c("Monocytes ", "Macrophages |Langerhans|Kupffer", "T ", "Fibroblasts", "Smooth Muscle Cells",
               "Endothelial Cells", "B Cells" , "Schwann Cells", "Plasma Cells", "Macrophages/Monocytes")





meta.data <- read.csv("E:\\8.Everyday_working\\R\\HCA_final_script-20190412-\\Figure_1\\meta_data_of_all_tissues_from_merged_data.txt",
                      sep = "\t",
                      header = T,
                      stringsAsFactors = F,
                      row.names = 1)
meta.data <- meta.data[meta.data$orig.ident != "Testis", ]
epith_cell_types <- unique(mapply(FUN = function(x) x[1:(length(x)-1)], str_split(string = meta.data$major_Cell_types_in_tissue, pattern = " ")))
epith_exp <- c("Epithelial|Esophageal Progenitor Cells|Proximal|Distal|AT2|High Proliferative Cells|AT1|Duct|Goblet|Merkel|Enterocytes|Absorptive|Stem Cells|Tuft|Enteroendocrine|Ionocytes")

table(meta.data$major_Cell_types_in_tissue)

epith_cells <- row.names(meta.data)[grep(pattern = epith_exp, x = meta.data$major_Cell_types_in_tissue)] 

pdf("KRT10.pdf") 
FeaturePlot(T., features = "KRT10")
dev.off()

######----------------------------------supplementary figures or tables---------------------------------#######################
all_DEG <- read.table("tissue_specific_genes_allFC_removed_high_exp_genes_T..csv",
                      sep = ",",
                      header = T,
                      stringsAsFactors = F,
                      row.names = 1)
top_200_DEGs <- read.table("tissue_specific_genes_top200_removed_high_exp_genes_T..csv",
                           sep = ",",
                           header = T,
                           stringsAsFactors = F,
                           row.names = 1)
fc_lth0.25  <- all_DEG %>% filter(p_val_adj < 0.05 & avg_logFC > 0.25)
sort(table((fc_lth0.25[, c(7)])))
length(unique(fc_lth0.25$gene))

DEGs_order_ac_fc <- data.frame()
for (tis in unique(all_DEG$cluster)) {
  DEGs_order_ac_fc_tmp <- all_DEG[all_DEG$cluster == tis, ] %>% arrange(desc(avg_logFC)) %>% filter(p_val < 0.05)
  DEGs_order_ac_fc <- rbind(DEGs_order_ac_fc, DEGs_order_ac_fc_tmp)
}

write.table(DEGs_order_ac_fc, "all_DEGs_order_ac_fc.csv", sep = ",", row.names = F, quote = F)

############---------------------------------FeaturePlot of interested genes---------------------------####################
png("featureplot_of_CD69_103.png", height = 7.5, width = 15, units = "in", res = 400)
FeaturePlot(T., features = c("ITGAE", "CD69"), pt.size = 0.5)  
dev.off()

png("featureplot_of_CD8_4_103_69.png", height = 15, width = 15, units = "in", res = 400)
FeaturePlot(subset_T, features = c("ITGAE", "CD69", "CD8A", "CD4"), pt.size = 0.5)  
dev.off()


pbmc.data <- Read10X(data.dir = getwd())
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")



subset_cells_markers %>% TOP_N(5000, sig.padj = 0.05) -> 
  write.table(subset_cells_markers, "T_organ_specfic.genes.csv", quote = F, sep = ",")


png(paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/7-HCA_final_results/Figure_2/T_organ_specfic_heatmap_to_20.png"),
    width = 15,
    height = 15,
    res = 300, units = "in")

pdf("T_organ_specfic_heatmap_to_20.pdf", height = 15, width = 10)
p5 <- DoHeatmap(object = subset_cells, features = all$gene, size = 2) + NoLegend() +
  theme(axis.text.x = element_text(size = 2),
        axis.text.y = element_text(size = 2),
        axis.title.x = element_text(colour = "red", size = 2),
        axis.title.y = element_text(colour = "black", size = 2))
print(p5)
dev.off()

###-------------------------------------part 2 T cell ratio in each tissues--------------
setwd("E:\\8.Everyday_working\\R\\HCA_final_script-20190412-\\Figure_1")
meta.dat <- read.table("meta_data_of_all_tissues_from_merged_data.txt",
                       sep = "\t",
                       header = T,
                       stringsAsFactors = F,
                       comment.char = "",
                       row.names = 1) %>% select(1, 9)
cell_retios <- data.frame()
for(tis in unique(meta.dat$orig.ident)){
  cellnumbers <- meta.dat[meta.dat$orig.ident == tis, 2]
  tmp <- data.frame(table(cellnumbers)/length(cellnumbers), tis)
  cell_retios <- rbind(cell_retios, tmp)
}

write.table(cell_retios, "cell_retios.csv", row.names = T, quote = F, sep = ",")




