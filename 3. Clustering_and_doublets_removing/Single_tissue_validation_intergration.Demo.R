library(Seurat)##Seurat V3.1.5 was used
library(stringr)
library(dplyr)
library(ggplot2)
library(parallel)
library(stringi)
library(plyr)
library(ggthemes)
library(data.table)
library(future)
library(ggsci)
library(patchwork)

options(future.globals.maxSize = 100*1000 * 1024^2)
plan("multiprocess", workers = 10)
plan()

color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10))[-8]

###------------------------1. read the data-----------------------------------###
sample_list <- list()

for (timepoint in dir() %>% grep(pattern = "Sample", value = T)) {
  subset_cell <- Read10X(data.dir = paste0("./", timepoint, "/"))
  colnames(subset_cell) <- paste0(colnames(subset_cell), "_", timepoint)
  subset_cell <- CreateSeuratObject(counts = subset_cell, project = timepoint, min.cells = 3, min.features = 200)
  subset_cell[["percent.mt"]] <- PercentageFeatureSet(subset_cell, pattern = "^MT-")
  sample_list[[timepoint]] <- subset_cell
}

# sample_list[["bladder1"]] <- bladder1 ## GX: 13,490 cells, 3 samples
# sample_list[["bladder2"]] <- bladder2 ## GGJ: 1,193 cells, 2 samples

for (i in 1:length(sample_list)) {
  sample_list[[i]] <- NormalizeData(sample_list[[i]], verbose = FALSE)
  sample_list[[i]] <- FindVariableFeatures(sample_list[[i]], selection.method = "vst", 
                                           nfeatures = 2000, verbose = FALSE)
}

##-------optional: remove sample with less than 2000 cells-------------
number <- length(sample_list)
for (i in 1:number) {
  if (sample_list[[i]] %>% dim %>% `[`(2)  %>% `>`(2000)) {
    print(c(i, dim(sample_list[[i]] %>% dim %>% `[`(2))))
  } else {
    sample_list[[i]] <- NULL
    print(c(i, "Changed !"))
  }
}

Bladder.anchors <- FindIntegrationAnchors(object.list = sample_list, dims = 1:30, k.filter = 200, verbose = T)
Bladder.anchors <- IntegrateData(anchorset = Bladder.anchors, dims = 1:30)

# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(Bladder.anchors) <- "integrated"

# Run the standard workflow for visualization and clustering
Bladder.anchors <- ScaleData(Bladder.anchors, verbose = T)
Bladder.anchors <- RunPCA(Bladder.anchors, npcs = 30, verbose = T)

Bladder.anchors <- FindNeighbors(Bladder.anchors, dims = 1:30)
Bladder.anchors <- FindClusters(Bladder.anchors, resolution = 1.5)

#Bladder.anchors <- RunUMAP(Bladder.anchors, reduction = "pca", dims = 1:30, reduction.model = "umap-learn")
Bladder.anchors <- RunTSNE(Bladder.anchors, reduction = "pca", dims = 1:30)

DefaultAssay(Bladder.anchors) <- "RNA" ##69,274 cells, 26 samples. 454 COCH cells
# Normalize RNA data for visualization purposes
Bladder.anchors <- NormalizeData(Bladder.anchors, verbose = FALSE)

png("Bladder_Gx_Vlnplot_of_markers.png", height = 7.5, width = 15, units = "in", res = 400)
VlnPlot(Bladder.anchors,
        features = c("MMP2", "ACTA2"), 
        pt.size = 0.0001,
        log = F,
        ncol = 1,
        cols = color_used)  
dev.off()

###---------------------------------------------dim plot-------------------------###
png("Bladder_tSNE_by_cluster_1.png",
    width = 15, height = 15, units = "in", res = 300)
p2 <- DimPlot(object = Bladder.anchors, reduction = 'tsne', label = F, pt.size = 1, cols = color_used) + NoLegend()
print(p2)
dev.off()

png("Bladder_GGJ_tSNE_by_cluster_with_legend_1.png",
    width = 15, height = 15, units = "in", res = 300)
p5 <- DimPlot(object = Bladder.anchors, reduction = 'tsne', label = TRUE, pt.size = 1, cols = color_used)
print(p5)
dev.off()

png("Bladder_tSNE_by_tissue_1.png",
    width = 15, height = 15, units = "in", res = 300)
p3 <- DimPlot(object = Bladder.anchors, reduction = 'tsne', group.by = "orig.ident", pt.size = 1,
              cols = color_used)
print(p3)
dev.off()

png("Bladder_tSNE_by_tissue_nolegend_1.png",
    width = 15, height = 15, units = "in", res = 300)
p4 <- DimPlot(object = Bladder.anchors, reduction = 'tsne', group.by = "orig.ident", pt.size = 1, 
              cols = color_used) +
  NoLegend()
print(p4)
dev.off()

###-------------2. find all markers-------------------
Bladder.anchors.markers <- FindMarkers_parallel(subset_cells, mc.cores = 10)

Bladder.anchors.markers %>% TOP_N(50, pct.1 = 0.2) -> top50
Bladder.anchors.markers <- Bladder.anchors.markers %>% TOP_N(5000)

write.table(top50,
            file = "Bladder.anchors_culster_top50_DEGs.csv",
            sep = ",",
            row.names = T,
            quote = F)

write.table(Bladder.anchors.markers,
            file = "Bladder.anchors_culster_all_DEGs.csv",
            sep = ",",
            row.names = T,
            quote = F)


png("subset_cells_fib_Vlnplot_of_markers.png", height = 15, width = 15, units = "in", res = 400)
VlnPlot(subset_cells_fib,
        features = c("PTPRC", "MMP2", "COCH", "ACTA2"),
        pt.size = 0.0001,
        log = F,
        ncol = 1,
        cols = color_used)  
dev.off()


