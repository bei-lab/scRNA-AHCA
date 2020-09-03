library(Seurat)
library(dplyr)

##----------This code were from the "Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris" with a few modification. These human dissociation-related genes 
##----------are the homolog of mouse provied by the article.  

dissociation_genes <-  read.table("Dissociation_induced_gene_human.txt", sep = "\t", header = F,row.names = NULL, stringsAsFactors = F) %>% unlist()

tiss <-  RunPCA(HCA_data, features  = dissociation_genes)
tiss <- ProjectDim(object = tiss)

tiss$plate.barcode.x <- tiss@meta.data %>% row.names()

dissociation_genes_pca_medians.pdf <- FetchData(tiss, vars = c('orig.ident', 'PC_1', 'PC_2', 'Color_of_tissues')) %>% 
  group_by(orig.ident, Color_of_tissues) %>%
  summarize(median_PC1 = median(PC_1),
            median_PC2 = median(PC_2)) %>%
  ggplot(aes(x = median_PC1, y = median_PC2, color = Color_of_tissues)) + geom_point(aes(color = Color_of_tissues)) +
  scale_color_identity(breaks = tiss$orig.ident %>% unique, 
                       labels = tiss$Color_of_tissues %>% unique, 
                       guide = "legend") + 
  guides(colour = guide_legend(override.aes = list(size=2)))

ggsave('dissociation_genes_pca_medians.pdf')

dissociation_genes_in_data <-  dissociation_genes[dissociation_genes %in% rownames(tiss@assays$RNA@scale.data)]
length(dissociation_genes_in_data)
tiss@meta.data$total.dissociation <- Matrix::colMeans(tiss@assays$RNA@counts[dissociation_genes_in_data,])


dissociation_genes_tissue_distributions.pdf <- FetchData(tiss, vars = c('plate.barcode.x', 'orig.ident', 'Color_of_tissues', 'total.dissociation'))%>%
  ggplot(aes(total.dissociation, ..density.., colour = Color_of_tissues)) +
  geom_freqpoly(binwidth = .05) +
  facet_wrap(~ orig.ident, ncol = 5) +
  scale_color_identity(breaks = tiss$Color_of_tissues %>% unique,
                       labels = tiss$orig.ident %>% unique,
                       guide = "legend") +
  ggtitle("Total dissociation")# + 
#theme(legend.position = "none")
ggsave('dissociation_genes_tissue_distributions.pdf')

