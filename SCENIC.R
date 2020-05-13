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

subset_cells@assays$RNA@data <- regulonAUC@assays@data$AUC
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

row.names(regulonAUC@assays@data$AUC) <- regulonAUC@assays@data$AUC %>% row.names %>% gsub(., pattern = "_extended*", replacement = "")
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

