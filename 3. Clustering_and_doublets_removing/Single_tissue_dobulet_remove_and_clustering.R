library(Seurat)##Seurat V3.1.5 was used
library(stringr)
library(dplyr)
library(ggplot2)
library(parallel)
library(stringi)
library(data.table)
library(RColorBrewer)
library(DoubletFinder)
library(future)

##------optional------parallel environment seting--cores = 10, memory = 100G----##
options(future.globals.maxSize = 10*1000 * 1024^2)
plan("multiprocess", workers = 10)

color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10))[-8]

###-----------------------1. the first clustering-----------------------####
ags <- commandArgs(trailingOnly = T)
dir.create(paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/9-reanalysis/1-cluster_each_organ/", ags))
setwd(paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/9-reanalysis/1-cluster_each_organ/", ags))

##------------------------1). Read in all input expression matrices
# Create and setup Seurat objects for each dataset 
TenXdat1 <- Read10X(data.dir = paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/2-Mapping/", ags[1],"/outs/filtered_feature_bc_matrix"))
TenXdat2 <- CreateSeuratObject(counts = TenXdat1, min.cells = 3, min.features = 500, project = ags[1])

##------------------------2).remove the genes in minimal cells (>=0.1% total cells)
if(dim(TenXdat2@assays$RNA@data)[2] <= 3000){
  TenXdat <- TenXdat2
} else {
  TenXdat <- CreateSeuratObject(counts = TenXdat1, min.cells = round(dim(TenXdat2@assays$RNA@data)[2]/1000),
                                min.features = 500, project = ags[1])
}

##------------------------3).cells quality filtering
mito.features <- grep(pattern = "^MT-", x = rownames(x = TenXdat), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = TenXdat, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = TenXdat, slot = 'counts'))

TenXdat[["percent.mito"]] <- percent.mito

TenXdat <- subset(x = TenXdat, subset = nFeature_RNA >= 500 & nFeature_RNA < 25000 & percent.mito <= 0.25 & nCount_RNA >= 1000 & nCount_RNA <= 500000)
TenXdat <- NormalizeData(object = TenXdat, normalization.method = "LogNormalize", scale.factor = 1e4)
TenXdat <- FindVariableFeatures(object = TenXdat, selection.method = 'mean.var.plot', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf), nfeatures = 2000)
TenXdat <- ScaleData(object = TenXdat, features = rownames(x = TenXdat), vars.to.regress = c("nCount_RNA", "percent.mito"))

TenXdat <- RunPCA(object = TenXdat, features = VariableFeatures(object = TenXdat), verbose = FALSE)

##------------------------4). dims and resulations used
subset_cells <- TenXdat
all_tissues <- c('Bladder_cDNA', 'Blood_cDNA', 'Common.bile.duct_cDNA', 'Esophagus_cDNA', 'Heart_cDNA',
                 'Liver_cDNA', 'Lymph.node_cDNA', 'Marrow_cDNA', 'Muscle_cDNA', 'Rectum_cDNA', 'Skin_cDNA',
                 'Small.intestine_cDNA', 'Spleen_cDNA', 'Stomach_cDNA', 'Trachea_cDNA')

TS <- ags
if(TS == 'Blood_cDNA' ){
  dim.use <- 15; res.use <- 0.4
} else if(TS == 'Liver_cDNA'){
  dim.use <- 30; res.use <- 0.8
} else if(TS == 'Lymph.node_cDNA'){
  dim.use <- 25; res.use <- 1
} else if(TS == 'Marrow_cDNA'){
  dim.use <- 25; res.use <- 0.6
} else if(TS %in% c('Small.intestine_cDNA', 'Spleen_cDNA')){
  dim.use <- 25; res.use <- 0.8
} else if(TS == 'Trachea_cDNA'){
  dim.use <- 35; res.use <- 1
} else {
  dim.use <- 30; res.use <- 1
}

subset_cells <- FindNeighbors(object = subset_cells, dims = 1:dim.use)
subset_cells <- FindClusters(object = subset_cells, resolution = res.use)

#ags <- ags %>% str_split(pattern = "_", simplify = T) %>% `[`(1)

##------------------------5). Run the TSNE
subset_cells <- RunTSNE(object = subset_cells, dims = 1:dim.use)

##------------------------## plot the first run tSNE
pdf(paste0(ags, "_", dim.use, "_", res.use, "_First_Run_tSNE.pdf"))
DimPlot(object = subset_cells, reduction = 'tsne', label = TRUE, cols = color_used)
dev.off()

write.table(data.frame(Tissue = ags[1], genes = dim(subset_cells@assays$RNA@data)[1], cells = dim(subset_cells@assays$RNA@data)[2]),
            file = paste0(ags[1], "_before_dobuletfinder.txt"),
            sep = "\t", row.names = F, quote = F)

##------------------------2.doublet finder processing---------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(subset_cells, PCs = subset_cells@commands$FindNeighbors.RNA.pca$dims)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- bcmvn[which(max(bcmvn$BCmetric) == bcmvn$BCmetric),2] %>% as.character %>% as.numeric

##------------------------1). Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(subset_cells@active.ident)
nExp_poi <- round(0.05*length(colnames(subset_cells)))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

##------------------------2). Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
subset_cells <- doubletFinder_v3(subset_cells, PCs = subset_cells@commands$FindNeighbors.RNA.pca$dims, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE)
subset_cells <- doubletFinder_v3(subset_cells, PCs = subset_cells@commands$FindNeighbors.RNA.pca$dims, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = grep("pANN", names(subset_cells@meta.data), value = T))

##------------------------3). Plot results --------------------------------------------------------------------------------------------------------------
high_of_low <- (subset_cells@meta.data[, grep("^DF\\.classifications", names(subset_cells@meta.data), value = T)] == "Singlet")+0
subset_cells@meta.data[high_of_low[, 1] + high_of_low[, 2] == 2, "DF_hi.lo"] <- "Singlet"
subset_cells@meta.data[high_of_low[, 1] + high_of_low[, 2] == 1, "DF_hi.lo"] <- "Doublet_lo"
subset_cells@meta.data[high_of_low[, 1] + high_of_low[, 2] == 0, "DF_hi.lo"] <- "Doublet_hi"

##------------------------3.remove the doublets and recluster-----------------#######
subset_cells <- subset_cells[, subset_cells@meta.data[grepl(subset_cells$DF_hi.lo, pattern = "Singlet"), ] %>% row.names()]

##------------------------. recluster-----------------#######
subset_cells <- NormalizeData(object = subset_cells, normalization.method = "LogNormalize", scale.factor = 1e4)
subset_cells <- FindVariableFeatures(object = subset_cells, selection.method = 'mean.var.plot', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf), nfeatures = 2000)
subset_cells <- ScaleData(object = subset_cells, features = rownames(x = subset_cells), vars.to.regress = c("nCount_RNA", "percent.mito"))

subset_cells <- RunPCA(object = subset_cells, features = VariableFeatures(object = subset_cells), verbose = FALSE)

subset_cells <- FindNeighbors(object = subset_cells, dims = 1:dim.use)
subset_cells <- FindClusters(object = subset_cells, resolution = res.use)

##------------------------ Run the TSNE
subset_cells <- RunTSNE(object = subset_cells, dims = 1:dim.use)

##------------------------ plot the tSNE
pdf(paste0(ags, "_", dim.use, "_", res.use, "_doublets_removed_tSNE.pdf"))
DimPlot(object = subset_cells, reduction = 'tsne', label = TRUE, cols = color_used)
dev.off()

write.table(data.frame(Tissue = ags, genes = dim(subset_cells@assays$RNA@data)[1], cells = dim(subset_cells@assays$RNA@data)[2]),
            file = paste0(ags[1], "_after_dobuletfinder.txt"),
            sep = "\t", row.names = F, quote = F)

plan("multiprocess", workers = 1)
plan()

##------------------------ find all the markers
tissue.markers <- FindMarkers_parallel(subset_cells, mc.cores = 10)

tissue.markers %>% TOP_N(n = 50) -> top50
write.table(top50, paste0(ags, "top50.csv"), sep = ",", row.names = T, quote = F)
write.table(tissue.markers %>% TOP_N(n = 2000), paste0(ags, ".csv"), sep = ",", row.names = T, quote = F)

##------------------------ plot the heat map

tissue.markers %>% TOP_N(n = 10) -> top10
pdf(paste0(ags, "_", dim.use, "_", res.use,"_doublets_removed_heatmap.pdf"), width = 14, height = 7)
Fixed_DoHeatmap(object = subset_cells, features = top10$gene, size = 2) + NoLegend() +
  theme(axis.text.x = element_text(size = 0), ##control the x label of cell barcodes
        axis.text.y = element_text(size = 0) ##control the gene label
  )
dev.off()

##------------------------ plot the vlnplot

png(paste0(ags, "_Vlnplot_of_markers.png"), height = 25, width = 15, units = "in", res = 400)
VlnPlot(subset_cells,
        features = c("CD3E", "CD3D", "FCGR3A", "MS4A1", "XBP1", "CD14", "MMP2", "ACTA2", "PECAM1", "CD8A", "CD8B", "CD4", "TRDV1", "TRDV2", "PAX7", "EPCAM"), ## active and repression marker genes
        # features = c("MKI67", "PCNA", "TYMS"),
        pt.size = 0.01,
        log = F,
        ncol = 1,
        cols = color_used)  
dev.off()

subset_cells[["tissue"]] <- ags[1]
assign(ags[1], subset_cells)

###it is so import that "list" args was used in save function.
save(list = ags[1], file = paste0(ags[1], ".RData"))

## package information
sessionInfo()
R version 3.6.3 (2020-02-29)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux Server release 6.6 (Santiago)

Matrix products: default
BLAS/LAPACK: /data/home/heshuai/Miniconda3/envs/RV_3.6/lib/libopenblasp-r0.3.10.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] future_1.18.0       DoubletFinder_2.0.3 RColorBrewer_1.1-2  data.table_1.13.0   stringi_1.4.6       ggplot2_3.3.2       dplyr_1.0.1         stringr_1.4.0       Seurat_3.1.5       

loaded via a namespace (and not attached):
 [1] ggrepel_0.8.2      Rcpp_1.0.5         rsvd_1.0.3         ape_5.4            lattice_0.20-41    ica_1.0-2          tidyr_1.1.1        listenv_0.8.0      png_0.1-7          zoo_1.8-8          digest_0.6.25      lmtest_0.9-37     
[13] R6_2.4.1           plyr_1.8.6         ggridges_0.5.2     httr_1.4.2         pillar_1.4.6       rlang_0.4.7        lazyeval_0.2.2     irlba_2.3.3        leiden_0.3.3       Matrix_1.2-18      reticulate_1.16    splines_3.6.3     
[25] Rtsne_0.15         htmlwidgets_1.5.1  uwot_0.1.8         igraph_1.2.5       munsell_0.5.0      sctransform_0.2.1  compiler_3.6.3     pkgconfig_2.0.3    globals_0.12.5     htmltools_0.5.0    tidyselect_1.1.0   gridExtra_2.3     
[37] tibble_3.0.3       RANN_2.6.1         codetools_0.2-16   fitdistrplus_1.1-1 viridisLite_0.3.0  withr_2.2.0        crayon_1.3.4       MASS_7.3-51.6      rappdirs_0.3.1     grid_3.6.3         tsne_0.1-3         nlme_3.1-148      
[49] jsonlite_1.7.0     gtable_0.3.0       lifecycle_0.2.0    magrittr_1.5       scales_1.1.1       KernSmooth_2.23-17 future.apply_1.6.0 pbapply_1.4-2      reshape2_1.4.4     ROCR_1.0-11        ellipsis_0.3.1     generics_0.0.2    
[61] vctrs_0.3.2        cowplot_1.0.0      RcppAnnoy_0.0.16   tools_3.6.3        glue_1.4.1         purrr_0.3.4        survival_3.2-3     colorspace_1.4-1   cluster_2.1.0      plotly_4.9.2.1     patchwork_1.0.1   

