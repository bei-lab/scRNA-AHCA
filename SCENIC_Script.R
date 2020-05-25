library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(ggthemes)
library(parallel)
library(dplyr)
library(tidyverse)
library(reshape2)
library(SCENIC)
library(ggsci)
library(AUCell)
library(ComplexHeatmap)
library(pheatmap)


###-------1. basical settings--------------------
org <- "hgnc" # or hgnc, or dmel
dbDir <- "/data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/8-reanalysis/All_singlet_data/" # RcisTarget databases location
myDatasetTitle <- "HCL" # choose a name for your analysis
data(defaultDbNames)
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org = org, dbDir = dbDir, dbs = dbs, datasetTitle = myDatasetTitle, nCores = 40) 
saveRDS(scenicOptions, file = paste0("int/scenicOptions.Rds"))

##----2. input file----------------------------------
## only 250 cells were selected in each cluster

meta.dat <- subset_cells@meta.data
tmp.dat <- data.frame()
selected_cells <- c()

for (cluster in unique(meta.dat$seurat_clusters %>% levels)) {
  tmp <- meta.dat %>% subset(seurat_clusters %in% cluster)
  if(dim(tmp)[1] < 250){
    cat("do nothig \n")
    selected_cells <- c(selected_cells, row.names(tmp))
    tmp.dat <- rbind(tmp.dat, tmp)
  } else {
    set.seed(1)
    tmp <- tmp[sample(1:dim(tmp)[1], replace = FALSE, size = 250), ]
    tmp.dat <- rbind(tmp.dat, tmp)
    selected_cells <- c(selected_cells, row.names(tmp))
  }
}

subset_cells_selected <- subset_cells[, selected_cells]

cellInfo <- data.frame(seuratCluster = Idents(subset_cells_selected))
singleCellMatrix <- subset_cells_selected@assays$RNA@counts %>% as.matrix

#-------3. remove low-frequency genes------###
genesKept <- geneFiltering(singleCellMatrix, scenicOptions = scenicOptions,
                           minCountsPerGene = 3*.01*ncol(singleCellMatrix),
                           minSamples = ncol(singleCellMatrix)*.05)

exprMat_filtered <- singleCellMatrix[genesKept, ]
dim(exprMat_filtered)

runCorrelation(exprMat_filtered, scenicOptions)
exprMat_log <- log2(singleCellMatrix + 1) 

#-------4. Run GENIE3-----------#####
exprMat_filtered_log <- log2(exprMat_filtered + 1) 
runGenie3(exprMat_filtered_log, scenicOptions)

#############-----------------------------###################
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 40
scenicOptions@settings$seed <- 123
#--------5. Run SCENIC---------########
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions) #
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

######---------to find the most significantly enriched TF in each cluster--------------######

SUBSET_CELLS <- subset_cells
subset_cells <- subset_cells_selected

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

subset_cells@assays$RNA@data <- regulonAUC@assays$data$AUC
genes_modified <- subset_cells %>% row.names %>% gsub(., pattern="_extended*", replacement = "")
row.names(subset_cells@assays$RNA@data) <- genes_modified

# Idents(subset_cells) <- subset_cells$orig.ident

result <- mclapply(levels(subset_cells@active.ident),
                   FUN =  function(x) {FindMarkers(subset_cells, ident.1 = x, ident.2 = NULL, logfc.threshold = 0)},
                   mc.cores = 36)
RESULT <- result

all_markers <- do.call(rbind, result)
all_markers$gene <- unlist(mapply(rownames, result, SIMPLIFY = F))
all_markers$cluster <- rep(levels(subset_cells@active.ident), times = mapply(dim, result, SIMPLIFY = TRUE)[1,])
subset_cells.markers <- all_markers

row.names(regulonAUC@assays@data$AUC) <- regulonAUC@assays$data$AUC %>% row.names %>% gsub(., pattern = "_extended*", replacement = "")
tissue_specific_gene <- subset_cells.markers %>% TOP_N(10, pct.1 = 0.2, sig.padj = 0.05, fc.threshold = 0.01) %>% `[`("gene") %>% unlist %>% unique()

Cells <- "Epi_HCA"

write.table(subset_cells.markers %>% TOP_N(1000, pct.1 = 0.2, sig.padj = 0.05, fc.threshold = 0.01), paste0(Cells, "_order_of_TFs.txt"), row.names = T, col.names = T, sep = "\t", quote = F)

##---plot with average value---------------####
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[, cells]))

# TFs <- grep(row.names(regulonActivity_byCellType), pattern = "extended", value = T, invert = T)
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType[, ]), center = T, scale=T))
row.names(regulonActivity_byCellType_Scaled) <- regulonActivity_byCellType_Scaled %>% row.names %>% gsub(., pattern="_extended*", replacement = "")

###---preparate the annotation materials-------------- 
anto <- data.frame(cluster = colnames(regulonActivity_byCellType_Scaled))
row.names(anto) <- colnames(regulonActivity_byCellType_Scaled)

cluster <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10))[-8]
# names(cluster) <- 0:((cellInfo$seuratCluster %>% unique %>% length) - 1)
names(cluster) <- anto$cluster %>% as.character
# cluster <- cluster[1:17]

ancols <- list(cluster = cluster)

# png("Epi_to10_significant_cells.reg.png", height = 15, width = 15, units = "in", res = 400)
pdf("Epi_to10_significant_cells.reg.pdf", height = 15, width = 10)
pheatmap::pheatmap(regulonActivity_byCellType_Scaled[tissue_specific_gene %>% unique, ], 
                   color = colorRampPalette(c("#181495", "blue", "white", "red", "#930107"))(100),
                   breaks = seq(-3, 3, length.out = 100),
                   treeheight_row = 10,
                   treeheight_col = 10,
                   border_color = "white",
                   cellwidth = 12,
                   cellheight = 3,
                   fontsize_row  = 3,
                   fontsize_col = 6,
                   cluster_rows = T,
                   cluster_cols = F,
                   width = 15,
                   height = 7.5,
                   annotation_col = anto,
                   annotation_colors = ancols,
                   )
dev.off()

                                     
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
[1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] pheatmap_1.0.12      ComplexHeatmap_2.2.0 AUCell_1.8.0         SCENIC_1.1.2-2       reshape2_1.4.3       forcats_0.5.0        stringr_1.4.0        purrr_0.3.3          readr_1.3.1          tidyr_1.0.0          tibble_2.1.3        
[12] tidyverse_1.3.0      ggthemes_4.2.0       RColorBrewer_1.1-2   Seurat_3.1.2         rlang_0.4.5          scales_1.1.0         dplyr_0.8.3          ggplot2_3.2.1        ggsci_2.9           

loaded via a namespace (and not attached):
  [1] reticulate_1.14             R.utils_2.9.2               tidyselect_0.2.5            RSQLite_2.1.5               AnnotationDbi_1.48.0        htmlwidgets_1.5.1           BiocParallel_1.20.0         Rtsne_0.15                 
  [9] munsell_0.5.0               codetools_0.2-16            mutoss_0.1-12               ica_1.0-2                   future_1.15.1               withr_2.1.2                 colorspace_1.4-1            Biobase_2.46.0             
 [17] rstudioapi_0.10             stats4_3.6.1                SingleCellExperiment_1.8.0  ROCR_1.0-7                  gbRd_0.4-11                 listenv_0.8.0               Rdpack_0.11-1               GenomeInfoDbData_1.2.2     
 [25] mnormt_1.5-5                bit64_0.9-7                 vctrs_0.2.4                 generics_0.0.2              TH.data_1.0-10              R6_2.4.1                    GenomeInfoDb_1.22.0         clue_0.3-57                
 [33] rsvd_1.0.2                  bitops_1.0-6                DelayedArray_0.12.1         assertthat_0.2.1            promises_1.1.0              SDMTools_1.1-221.2          multcomp_1.4-11             gtable_0.3.0               
 [41] npsurv_0.4-0                globals_0.12.5              sandwich_2.5-1              GlobalOptions_0.1.1         splines_3.6.1               lazyeval_0.2.2              broom_0.5.5                 modelr_0.1.6               
 [49] backports_1.1.6             httpuv_1.5.2                tools_3.6.1                 gplots_3.0.1.1              BiocGenerics_0.32.0         ggridges_0.5.1              TFisher_0.2.0               Rcpp_1.0.3                 
 [57] plyr_1.8.5                  zlibbioc_1.32.0             RCurl_1.95-4.12             pbapply_1.4-2               GetoptLong_0.1.7            cowplot_1.0.0               S4Vectors_0.24.1            zoo_1.8-6                  
 [65] SummarizedExperiment_1.16.1 haven_2.2.0                 ggrepel_0.8.1               cluster_2.1.0               fs_1.3.1                    magrittr_1.5                data.table_1.12.8           circlize_0.4.8             
 [73] lmtest_0.9-37               reprex_0.3.0                RANN_2.6.1                  mvtnorm_1.0-11              fitdistrplus_1.0-14         matrixStats_0.55.0          hms_0.5.3                   lsei_1.2-0                 
 [81] mime_0.8                    xtable_1.8-4                XML_3.98-1.20               readxl_1.3.1                IRanges_2.20.1              gridExtra_2.3               shape_1.4.4                 compiler_3.6.1             
 [89] KernSmooth_2.23-16          crayon_1.3.4                R.oo_1.23.0                 htmltools_0.4.0             later_1.0.0                 RcppParallel_4.4.4          lubridate_1.7.8             DBI_1.1.0                  
 [97] dbplyr_1.4.2                MASS_7.3-51.5               rappdirs_0.3.1              Matrix_1.2-18               cli_2.0.2                   R.methodsS3_1.7.1           gdata_2.18.0                metap_1.2                  
[105] igraph_1.2.4.2              GenomicRanges_1.38.0        pkgconfig_2.0.3             sn_1.5-4                    numDeriv_2016.8-1.1         plotly_4.9.1                xml2_1.2.2                  annotate_1.64.0            
[113] multtest_2.42.0             XVector_0.26.0              bibtex_0.4.2.1              rvest_0.3.5                 digest_0.6.23               sctransform_0.2.1           RcppAnnoy_0.0.14            tsne_0.1-3                 
[121] graph_1.64.0                cellranger_1.1.0            leiden_0.3.1                uwot_0.1.5                  GSEABase_1.48.0             shiny_1.4.0.2               gtools_3.8.1                rjson_0.2.20               
[129] lifecycle_0.2.0             nlme_3.1-143                jsonlite_1.6                viridisLite_0.3.0           fansi_0.4.1                 pillar_1.4.3                lattice_0.20-38             fastmap_1.0.1              
[137] httr_1.4.1                  plotrix_3.7-7               survival_3.1-12             glue_1.3.1                  png_0.1-7                   bit_1.1-14                  stringi_1.4.3               blob_1.2.0                 
[145] caTools_1.17.1.3            memoise_1.1.0               irlba_2.3.3                 future.apply_1.3.0          ape_5.3                    

