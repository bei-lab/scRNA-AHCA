library(monocle)
library(Seurat)
library(plyr)
library(pheatmap)
library(ggsci)
library(ggplo2)

### Take CD8 cells as an example.
subset_cells <- get(load("subset_cells.RData"))

expression_matrix <- subset_cells@assays$RNA@counts %>% as.matrix()
cell_metadata <- subset_cells@meta.data
gene_annotation <- data.frame(gene_short_name = row.names(subset_cells@assays$RNA@counts), row.names = row.names(subset_cells@assays$RNA@counts), stringsAsFactors = F)

pd <- new("AnnotatedDataFrame", data = cell_metadata)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
cds <- newCellDataSet(expression_matrix, phenoData = pd, featureData = fd)

DelayedArray:::set_verbose_block_processing(TRUE)

# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size = 1000e6)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
disp_table <- dispersionTable(cds)

## use different mean_expression and dispersion_empirical 
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.2 &
                           dispersion_empirical >= 0.1 * dispersion_fit)$gene_id ##

# ordering_genes <- ordering_genes[!ordering_genes %in% grep("MT-|^RPL|^RPS", ordering_genes, value = T)] ## optional
cds <- setOrderingFilter(cds, ordering_genes)
cds <- preprocessCDS(cds, num_dim = 20)
cds <- fixed_reduceDimension(cds, reduction_method = 'UMAP')
cds <- partitionCells(cds)
cds <- learnGraph(cds,  RGE_method = 'SimplePPT')

cols = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', 
         '#98df8a', '#d62728', '#d62728', '#ff9896', '#9467bd', 
         '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', 
         '#7f7f7f', '#c7c7c7', '#17becf', '#9edae5',  '#00cdaa',
         '#aacdaa')

# a helper function to identify the root principal points:
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  
  closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
                                           (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

MPP_node_ids = get_correct_root_state(cds,
                                      cell_phenotype = "annotation",
                                      root_type = "CCR7_CD8_Naive")
cds <- orderCells(cds, root_pr_nodes = MPP_node_ids)

png("CD8_monocle_monocyte_all_cells_by_clusters_facet_Pseudotime.png", height = 10, width = 10, res = 400, units = "in")
plot_cell_trajectory(cds)
dev.off()

# names(cols) <- cell_type_color
cols <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12))[-8]

png("CD8_trajectory_by_cluster.png", height = 15, width = 15, res = 400, units = "in")
plot_cell_trajectory(cds, 
                     color_by = "annotation", cell_size = 1, cell_name_size = 6)  + # facet_wrap(~orig.ident, nrow = 2) +
  scale_color_manual(values = cols) + theme(legend.position = c(0.7, 0.8),
                                            legend.background = element_rect(fill= "transparent"),
                                            legend.key.size = unit(0.5, 'lines')) +
  # # legend.direction = "horizontal")# +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)),
         col = guide_legend(ncol = 1, byrow = T)) +
  labs(colour= "clusters") +
  theme(legend.position = "none")
dev.off()

tissue_levels <- pData(cds) %>% select("orig.ident", "Color_of_tissues") %>% `[`(, "orig.ident") %>% levels()
cols <- mapvalues(tissue_levels, from = pData(cds) %>% select(orig.ident) %>% unlist %>% as.character,
                  to = pData(cds) %>% select(Color_of_tissues) %>% unlist %>% as.character)

png("CD8_trajectory_by_tissue_facet.png", height = 10, width = 20, res = 400, units = "in")
plot_cell_trajectory(cds, 
                     color_by = "orig.ident", cell_size = 0.5, cell_name_size = 6)  +  facet_wrap(~orig.ident, nrow = 2) +
  scale_color_manual(values = cols) + theme(legend.position = c(0.7, 0.8),
                                            legend.background = element_rect(fill= "transparent"),
                                            legend.key.size = unit(0.5, 'lines')) +
  # # legend.direction = "horizontal")# +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)),
         col = guide_legend(ncol = 1, byrow = T)) +
  labs(colour= "clusters") +
  theme(legend.position = "none")
dev.off()

cols <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12))[-8]

png("CD8_trajectory_by_cluster_facet.png", height = 10, width = 20, res = 400, units = "in")
plot_cell_trajectory(cds, 
                     color_by = "annotation", cell_size = 0.5, cell_name_size = 6)  +  facet_wrap(~annotation, nrow = 3) +
  scale_color_manual(values = cols) + theme(legend.position = c(0.7, 0.8),
                                            legend.background = element_rect(fill= "transparent"),
                                            legend.key.size = unit(0.5, 'lines')) +
  # # legend.direction = "horizontal")# +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)),
         col = guide_legend(ncol = 1, byrow = T)) +
  labs(colour= "clusters") +
  theme(legend.position = "none")
dev.off()

