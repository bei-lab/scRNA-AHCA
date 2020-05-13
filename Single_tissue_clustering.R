
library(Seurat)##Seurat V3.0 was used
library(Matrix)
library(stringr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(parallel)

ags <- commandArgs(trailingOnly = T)
dir.create(paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/4-Seurat_cluster/", ags))
setwd(paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/4-Seurat_cluster/", ags))

#1. Read in all input expression matrices
# Create and setup Seurat objects for each dataset 
TenXdat1 <- Read10X(data.dir = paste0("/data4/heshuai/RAW_data/1-SingleCell/3-HCA/2-Mapping/", ags[1],"/outs/filtered_feature_bc_matrix"))
TenXdat2 <- CreateSeuratObject(counts = TenXdat1, min.cells = 3, min.features = 500, project = ags[1])

#2.remove the genes in minimal cells (>=0.1% total cells)
if(dim(TenXdat2@assays$RNA@data)[2] <= 3000){
  TenXdat <- TenXdat2
} else {
  TenXdat <- CreateSeuratObject(counts = TenXdat1, min.cells = round(dim(TenXdat2@assays$RNA@data)[2]/1000),
                                min.features = 500, project = ags[1])
}

#3. cells quality filtering
mito.features <- grep(pattern = "^MT-", x = rownames(x = TenXdat), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = TenXdat, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = TenXdat, slot = 'counts'))

# The [[ operator can add columns to object metadata, and is a great place to stash QC stats
# TenXdat[['percent.ribo']] <- percent.ribo
TenXdat <- subset(x = TenXdat, subset = nFeature_RNA >= 500 & nFeature_RNA < 25000 & percent.mito <= 0.25 & nCount_RNA >= 1000 & nCount_RNA <= 500000)
TenXdat <- NormalizeData(object = TenXdat, normalization.method = "LogNormalize", scale.factor = 1e4)
TenXdat <- ScaleData(object = TenXdat, features = rownames(x = TenXdat), vars.to.regress = c("nCount_RNA", "percent.mito"))
TenXdat <- FindVariableFeatures(object = TenXdat, selection.method = 'mean.var.plot', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf), nfeatures = 1000)

TenXdat <- RunPCA(object = TenXdat, features = VariableFeatures(object = TenXdat), verbose = FALSE)
# TenXdat <- ProjectDim(object = TenXdat)

#4. dims and resulation choose
tissuename <- ags
tissue <- TenXdat
all_tissues <- c('bile_cDNA', 'bladder_cDNA', 'blood_cDNA', 'esophagus_cDNA', 'heart_cDNA', 'kidney_cDNA', 
                 'largeintestine_cDNA', 'liver_cDNA', 'lung_cDNA', 'lymphnode_cDNA', 'marrow_cDNA', 
                 'muscle_cDNA', 'pancreas_cDNA', 'skin_cDNA', 'smallintestine_cDNA', 'spleen_cDNA', 
                 'stomach_cDNA', 'testis_cDNA', 'Trachea_cDNA')
TS <- ags
if(TS == 'blood_cDNA' ){
    dim.use <- 15; res.use <- 0.4
  } else if(TS == 'liver_cDNA'){
    dim.use <- 30; res.use <- 0.8
  } else if(TS == 'lymphnode_cDNA'){
    dim.use <- 25; res.use <- 1
  } else if(TS == 'marrow_cDNA'){
    dim.use <- 25; res.use <- 0.6
  } else if(TS == 'pancreas_cDNA'){
    dim.use <- 35; res.use <- 1
  } else if(TS %in% c('smallintestine_cDNA', 'spleen_cDNA')){
    dim.use <- 25; res.use <- 0.8
  } else if(TS == 'testis_cDNA'){
    dim.use <- 34; res.use <- 1
  } else if(TS == 'Trachea_cDNA'){
    dim.use <- 35; res.use <- 1
  } else {
    dim.use <- 30; res.use <- 1
  }


tissue <- FindNeighbors(object = tissue, dims = 1:dim.use)
tissue <- FindClusters(object = tissue, resolution = res.use)

tissuename <- ags
#5. Run the TSNE
tissue <- RunTSNE(object = tissue, dims = 1:dim.use)
### plot the tSNE
pdf(paste0(tissuename, "_", dim.use, "_", res.use, "_tSNE.pdf"))
DimPlot(object = tissue, reduction = 'tsne', label = TRUE)
dev.off()

#6. find all the markers
result <- mclapply(as.numeric(levels(tissue@active.ident)),
                   FUN =  function(x) {FindMarkers(tissue, ident.1 = x, ident.2 = NULL)},
                   mc.cores = 16)

all_markers <- do.call(rbind, result)
all_markers$gene <- unlist(mapply(rownames, result))
all_markers$cluster <- rep(levels(tissue@active.ident), times = mapply(dim, result, SIMPLIFY = TRUE)[1,])
tissue.markers <- all_markers

tissue.markers %>% group_by(cluster) %>% TOP_N(n = 50, wt = avg_logFC) -> top50
write.table(top50, paste0(ags, "top50.csv"), sep = ",", row.names = T, quote = F)
write.table(tissue.markers, paste0(ags, ".csv"), sep = ",", row.names = T, quote = F)

#7. plot the heat map
tissue.markers %>% group_by(cluster) %>% TOP_N(n = 10, wt = avg_logFC) -> top10
pdf(paste0(tissuename, "_", dim.use, "_", res.use,"_heatmap.pdf"), width = 14, height = 7)
DoHeatmap(object = tissue, features = top10$gene, size = 2) + NoLegend() +
  theme(axis.text.x = element_text(size = 3),
        axis.text.y = element_text(size = 3),
        axis.title.x = element_text(colour = "red", size = 3),
        axis.title.y = element_text(colour = "black", size = 3))
dev.off()

tissue[["tissue"]] <- ags[1]
assign(ags[1], tissue)

###it is so import that "list" args was used in save function.
save(list = ags[1], file = paste0(ags[1], ".RData"))
write.table(data.frame(tissue = ags[1], genes = dim(TenXdat@assays$RNA@data)[1], cells = dim(TenXdat@assays$RNA@data)[2]),
            file = paste0(ags[1], ".txt"),
            sep = "\t", row.names = F, quote = F)

sessionInfo()
R version 3.6.1 (2019-07-05)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux Server release 6.6 (Santiago)

Matrix products: default
BLAS/LAPACK: /data/home/heshuai/anaconda3/lib/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] gridExtra_2.3 stringr_1.4.0 Matrix_1.2-18 Seurat_3.1.2  rlang_0.4.5   scales_1.1.0  dplyr_0.8.3   ggplot2_3.2.1 ggsci_2.9    

loaded via a namespace (and not attached):
 [1] tsne_0.1-3          nlme_3.1-143        bitops_1.0-6        RcppAnnoy_0.0.14    RColorBrewer_1.1-2  httr_1.4.1          numDeriv_2016.8-1.1 sctransform_0.2.1   tools_3.6.1         R6_2.4.1            irlba_2.3.3         KernSmooth_2.23-16 
[13] uwot_0.1.5          lazyeval_0.2.2      BiocGenerics_0.32.0 colorspace_1.4-1    sn_1.5-4            npsurv_0.4-0        withr_2.1.2         tidyselect_0.2.5    mnormt_1.5-5        compiler_3.6.1      Biobase_2.46.0      TFisher_0.2.0      
[25] plotly_4.9.1        sandwich_2.5-1      caTools_1.17.1.3    lmtest_0.9-37       mvtnorm_1.0-11      ggridges_0.5.1      pbapply_1.4-2       rappdirs_0.3.1      digest_0.6.23       R.utils_2.9.2       htmltools_0.4.0     pkgconfig_2.0.3    
[37] bibtex_0.4.2.1      plotrix_3.7-7       htmlwidgets_1.5.1   zoo_1.8-6           jsonlite_1.6        ica_1.0-2           gtools_3.8.1        R.oo_1.23.0         magrittr_1.5        Rcpp_1.0.3          munsell_0.5.0       ape_5.3            
[49] reticulate_1.14     lifecycle_0.2.0     R.methodsS3_1.7.1   stringi_1.4.3       multcomp_1.4-11     gbRd_0.4-11         MASS_7.3-51.5       gplots_3.0.1.1      Rtsne_0.15          plyr_1.8.5          grid_3.6.1          gdata_2.18.0       
[61] listenv_0.8.0       ggrepel_0.8.1       crayon_1.3.4        lattice_0.20-38     cowplot_1.0.0       splines_3.6.1       multtest_2.42.0     SDMTools_1.1-221.2  pillar_1.4.3        igraph_1.2.4.2      reshape2_1.4.3      future.apply_1.3.0 
[73] codetools_0.2-16    stats4_3.6.1        leiden_0.3.1        mutoss_0.1-12       glue_1.3.1          lsei_1.2-0          metap_1.2           RcppParallel_4.4.4  data.table_1.12.8   png_0.1-7           vctrs_0.2.4         Rdpack_0.11-1      
[85] tidyr_1.0.0         gtable_0.3.0        RANN_2.6.1          purrr_0.3.3         future_1.15.1       assertthat_0.2.1    rsvd_1.0.2          viridisLite_0.3.0   survival_3.1-12     tibble_2.1.3        cluster_2.1.0       globals_0.12.5     
[97] fitdistrplus_1.0-14 TH.data_1.0-10      ROCR_1.0-7         


