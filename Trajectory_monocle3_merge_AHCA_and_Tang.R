library(Seurat)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(data.table)
library(ggsci)
library(monocle3)

load("Tang_large_intestine.RData")
orig <- subset_cells

Rectum_cells <- subset_cells@meta.data %>% subset(orig.ident == "Rectum_cDNA") %>% row.names()#%>% subset(! annotation %in% c("cDC1", "cDC2") ) %>%  row.names()

AHCA <- subset_cells[, Rectum_cells]
TANG$seurat_clusters <- factor(19, levels = 19)
merged_seurat_object <- merge(AHCA, TANG)

merged_seurat_object$annotation[is.na(merged_seurat_object$annotation)] <- "Tang"
subset_cells <- merged_seurat_object

Idents(subset_cells) <- factor(subset_cells$seurat_clusters, levels = subset_cells$seurat_clusters %>% unique)
subset_cells <- subset(subset_cells, ident = c(2, 4, 10, 12, 19))

expression_matrix <- subset_cells@assays$RNA@counts %>% as.matrix()
cell_metadata <- subset_cells@meta.data
gene_annotation <- data.frame(gene_short_name = row.names(subset_cells@assays$RNA@counts), row.names = row.names(subset_cells@assays$RNA@counts), stringsAsFactors = F)

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds$batch <- ifelse(c(cds$seurat_clusters %>% as.numeric) < 19, "AHCA", "Tang")
cds <- preprocess_cds(cds, num_dim = 25) ## previous 30
cds <- align_cds(cds, alignment_group = "batch")

cds <- reduce_dimension(cds, umap.min_dist = 0.1, umap.n_neighbors = 15L)

# png("Trajectory_monocle3.png", height = 10, width = 10, res = 400, units = "in")
# plot_cells(cds, label_groups_by_cluster = FALSE,  color_cells_by = "annotation", cell_size = 1)
# dev.off()

cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = T)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="Intermediate_Mon_CCL20"){
  cell_ids <- which(colData(cds)[, "annotation"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))


png("Rectum_psedutime_test.png", height = 10, width = 10, res = 400, units = "in")
# pdf("monocyte3_Rectum_psedutime.pdf", height = 10, width = 10)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5,
           cell_size = 2)
dev.off()


png("Rectum_trajectory_annotation_test.png", height = 10, width = 10, res = 400, units = "in")
# pdf("monocyte3__Rectum_trajectory_annotation.pdf", height = 10, width = 10)
plot_cells(cds,
           color_cells_by = "annotation",
           label_cell_groups = F,
           label_leaves = T,
           label_branch_points = T,
           cell_size = 2,
           graph_label_size = 2
           ) +
  scale_color_manual(values = eval(parse(text = color_used)))
dev.off()

png("Rectum_by_clusters_genes.png", height = 10, width = 20, res = 400, units = "in")
plot_cells(cds,
           genes=c("MKI67", "PCNA"),
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           cell_size = 2)
dev.off()
