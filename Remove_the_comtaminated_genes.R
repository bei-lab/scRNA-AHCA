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

## Take the fibroblast and smooth muscle cell as an example.

### load the HCA_data_singlet
load("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/8-reanalysis/All_singlet_data/HCA_all_data_singlet.RData")
HCA_data <- HCA_all_data_singlet

##############--------------------------1. scale and normalize the data------------------------------------------#####################
Idents(HCA_data) <- factor(as.character(HCA_data$orig.ident))
celltypes <- c("Fibroblasts|Smooth Muscle Cells")

for (celltype in celltypes) {
  subset_cells <- subset(x = HCA_data, cells = names(grep(celltype, HCA_data$cluster_in_each_tissue, value = T)))
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

####------------------------------------3. calculate top 2% high expression genes-----------------------#############################################
##1.loading the data
# laod("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/5-HCA-all/HCA_data_dims30_res1_only_singlet/HCA_all_data_singlet.RData")
##2.calculate the top 2% high expression genes
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
            file = "/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/8-reanalysis/All_singlet_data/genes_of_2_percent_expression_for_fib_soom.txt",
            sep = ",",
            quote = F,
            row.names = T,
            col.names = NA)
####------------------------------------4. remove the low-frequency genes in each tissue-----------------######################################################
genes_of_2_percent_expression <- read.table("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/8-reanalysis/All_singlet_data/genes_of_2_percent_expression_for_fib_soom.txt",
                                            sep = ",",
                                            stringsAsFactors = F,
                                            header = T,
                                            row.names = 1)

for (celltype in celltypes) {
  DEG <- read.table(paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/8-reanalysis/All_singlet_data/", "HCA_allFC_max300cells_", make.names(celltype), ".csv"),
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
HCA_data <- HCA_all_data_singlet
HCA_data <- subset(HCA_data, cells = grep(pattern = "Testis", colnames(HCA_data), invert = T, value = T))
HCA_data <- subset(x = HCA_data,
                   cells = setdiff(names(grep(".*", HCA_data$cluster_in_each_tissue, value = T)),
                                   names(grep("Unassigned|Mixed", HCA_data$cluster_all , value = T))))
# subset_cells 
##############---------------------------- scale and normalize the data------------------------------------------#####################
celltypes <- c("Fibroblasts|Smooth Muscle Cells")

for (celltype in celltypes) {
  genes <- as.character(read.table(paste0(make.names(celltype), "_genes_to_removed.txt"),
                                   header = T,
                                   sep = "\t",
                                   stringsAsFactors = F,
                                   row.names = 1)[, 1])
  subset_cells <- subset(x = HCA_data,
                         cells = names(grep(celltype, HCA_data$cluster_in_each_tissue, value = T)),
                         features = setdiff(row.names(HCA_data), genes))
  # subset_cells <- subset(x = subset_cells,
  #                        subset = nFeature_RNA >= 500 & nFeature_RNA < 4000 & percent.mito <= 0.25 & nCount_RNA >= 1000 & nCount_RNA <= 10000)
  # print(table(subset_cells@active.ident))
  subset_cells <- NormalizeData(object = subset_cells, normalization.method = "LogNormalize", scale.factor = 1e4)
  subset_cells <- ScaleData(object = subset_cells, features = rownames(x = subset_cells))#, vars.to.regress = c("nCount_RNA", "percent.mito"))
  subset_cells <- FindVariableFeatures(object = subset_cells, selection.method = 'mean.var.plot', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf))
  subset_cells <- RunPCA(object = subset_cells, features = VariableFeatures(object = subset_cells), verbose = FALSE, npcs = 50)
  
  if(celltype %in% c("Monocytes|Macrophages|Langerhans", "T |Immune", "Fibroblasts|Smooth Muscle Cells", "Endothelial", "B Cells|Plasma Cells" )){
    dim.use <- 15
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
  png(paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/8-reanalysis/All_singlet_data/", make.names(celltype), "_tSNE_Round_1_", dim.use, "_", res.use, ".png"),
      width = 15, height = 15, units = "in", res = 300)
  p2 <- DimPlot(object = subset_cells, reduction = 'tsne', label = TRUE, pt.size = 1,) + NoLegend()
  print(p2)
  dev.off()
  
  png(paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/8-reanalysis/All_singlet_data/", make.names(celltype), "tSNE_Round_1_groupby_", dim.use, "_", res.use, ".png"),
      width = 15, height = 15, units = "in", res = 300)
  p3 <- DimPlot(object = subset_cells, reduction = 'tsne', group.by = "orig.ident", pt.size = 1,
                cols = unique(subset_cells$Color_of_tissues)[match(levels(subset_cells$orig.ident)[table(subset_cells@meta.data$orig.ident) > 0], unique(as.character(subset_cells$orig.ident)) )]) + NoLegend()
  print(p3)
  dev.off()
  
  png(paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/8-reanalysis/All_singlet_data/", make.names(celltype), "tSNE_Round_1_groupby_with_legend_", dim.use, "_", res.use, ".png"),
      width = 15, height = 15, units = "in", res = 300)
  p4 <- DimPlot(object = subset_cells, reduction = 'tsne', group.by = "orig.ident", pt.size = 1, 
                cols = unique(subset_cells$Color_of_tissues)[match(levels(subset_cells$orig.ident)[table(subset_cells@meta.data$orig.ident) > 0], unique(as.character(subset_cells$orig.ident)) )]) 
  print(p4)
  dev.off()
  
  ###-----------------------------Find all the marker genes-------------------------------------------
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
  
  assign(make.names(celltype), subset_cells)
  save(list = make.names(celltype),
       file = paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/8-reanalysis/All_singlet_data/",
                     make.names(celltype),
                     "_before_filtering.RData"))
  
  subset_cells.markers %>% TOP_N(50) -> top50
  subset_cells.markers <- subset_cells.markers %>% TOP_N(5000)
  
  write.table(top50,
              file = paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/8-reanalysis/All_singlet_data/",
                            make.names(celltype), "_top50_", dim.use, "_", res.use, ".csv"),
              sep = ",",
              row.names = T,
              quote = F)
  
  write.table(subset_cells.markers,
              file = paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/8-reanalysis/All_singlet_data/",
                            make.names(celltype), "_all_DEGs_", dim.use, "_", res.use, ".csv"),
              sep = ",",
              row.names = T,
              quote = F)
  
  ###7. plot the heat map
  subset_cells.markers %>% TOP_N(10) -> top10
  png(paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/8-reanalysis/All_singlet_data/",
             make.names(celltype), "_", dim.use,"_",res.use, "_heatmap.png"),
      width = 15, height = 15, units = "in", res = 300)
  
  p5 <- DoHeatmap(object = subset_cells, features = top10$gene, size = 2) + NoLegend() +
    theme(axis.text.x = element_text(size = 2),
          axis.text.y = element_text(size = 2),
          axis.title.x = element_text(colour = "red", size = 2),
          axis.title.y = element_text(colour = "black", size = 2))
  print(p5)
  dev.off()
}



