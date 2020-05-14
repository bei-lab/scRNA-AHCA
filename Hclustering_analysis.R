library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(data.table)
library(RColorBrewer)
library(reshape2)
library(dendextend)

### Take epithelial cell as an example.
subset_cells <- get(load("Epithelial_Cell_annotation.RData"))
cluster.averages <- AverageExpression(object = subset_cells, return.seurat = TRUE)
cluster.average <- as.matrix(cluster.averages@assays$RNA@data)

##------------------read the differently expressed genes list across different epithelial cell clusters------------------####
DEGs <- read.table("Epi_Cells_culster_all_DEGs_clean.csv", sep = ",", header = T, row.names = 1, stringsAsFactors = F) %>%
  TOP_N(1000, pct.1 = 0.2, sig.padj = 0.05, fc.threshold = 0.25)

spearman_cor <- cor(cluster.averages@assays$RNA@data %>% as.matrix() %>% `[`(DEGs$gene %>% unique(), ), method = "spearman")
row.names(spearman_cor) <- mapvalues(row.names(spearman_cor), from = subset_cells$seurat_clusters %>% as.character, to = subset_cells$annotation %>% as.character)
colnames(spearman_cor) <- mapvalues(colnames(spearman_cor), from = subset_cells$seurat_clusters %>% as.character, to = subset_cells$annotation %>% as.character)

distance_of_cluster <- as.dist(1-spearman_cor)

hc <- hclust(distance_of_cluster, method = "complete")
hc <- as.dendrogram(hc)
asggden <- as.ggdend(hc)

pdf("HCA_HCluster_Epi_tissue.pdf", width = 15, height = 15)
ggplot(asggden, horiz = TRUE, theme = NULL, offset_labels = -0.01) +
  ylim(0.8, -0.5) +
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 

dev.off()

##--------------or-----------##
pdf("HCA_HCluster_Epi_tissue.pdf", width = 15, height = 15)
plot(as.phylo(hc), type = "cladogram", cex = 0.6,
     edge.color = "steelblue", edge.width = 2, edge.lty = 2,
     tip.color = "steelblue")
dev.off()




