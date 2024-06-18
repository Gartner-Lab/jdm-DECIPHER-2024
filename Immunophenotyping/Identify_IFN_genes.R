#Loading packages
library(xlsx)
library(dendsort)
library(dittoSeq)
library(Seurat)
library(tidyverse)
library(dplyr)
library(vctrs)
library(ggpubr)
library(readr)
library(RColorBrewer)
library(ComplexHeatmap)
library(tidyr)
library(circlize)
library(readr)


#Define functions to sort heatmap
std.error <- function(x) sd(x)/sqrt(length(x))
callback = function(hc, ...){dendsort(hc)}
fun <- function(names){
  t <- sample(names)
  t <- which(t == 'TNJDM')
  t <- std.error(t)}
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

#Load data
markers_rna <- xlsx::read.xlsx("/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/markers_rna_tnjdm_hc_min_lfc_1_updated.xlsx", sheetIndex =1)
sobj_final <- readRDS(file = "/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/sobj_final.rds")


#Combine disease state and cluster
Idents(sobj_final) <- sobj_final$new_diseasegroups
disease_groups <- sobj_final$new_diseasegroups
cluster <- sobj_final$new_labels
new <- data.frame(cluster, disease_groups)
new$new <- paste0(new$cluster, "&", new$disease_groups)
cluster_group <- new$new
Idents(sobj_final) <- cluster_group
sobj_final[["cluster_group"]] <- Idents(sobj_final)

#Pseudobulk data
mean_exp <- AverageExpression(sobj_final, assays = "RNA", return.seurat = T, verbose = T)

#Split into disease group and cluster
mean_exp[["disease_group"]] <- Idents(mean_exp)
mean_exp[["cluster"]] <- Idents(mean_exp)

disease_group <- mean_exp$disease_group
levels(disease_group) <- sub("(.*)&","", levels(disease_group))
mean_exp$disease_group <- disease_group
Idents(mean_exp) <- mean_exp$disease_group

cluster <- mean_exp$cluster
levels(cluster) <- sub("&.*","", levels(cluster))
mean_exp$cluster <- cluster
Idents(mean_exp) <- mean_exp$cluster

#Remove cell types with less than 90 cells in active, HC and TNJDM
table(sobj_final$new_labels, sobj_final$new_diseasegroups)
mean_exp <- subset(mean_exp, cluster %in% c('B_naive4', 'PBs', 'CD8+ mem', 'CD56bright NK', 'CDCs'), invert = T)
mean_exp$cluster <- droplevels(mean_exp$cluster)

#Filter DEGs. Exclude groups with less than 90 cells in HC TNJDM (the compared groups) and include only genes
#that are differentially expressed in min. 2 cell types
genes <- markers_rna %>%
  dplyr::filter(! celltype %in% c('Plasmablasts', 'CD8+ memory_resting', 'CD56bright NK', 'Classic dendritic')) %>%
  dplyr::select(gene) %>%
  dplyr::group_by(gene) %>%
  dplyr::count()%>%
  dplyr::filter(n > 1)
genes <- genes$gene


#Create heatmap with genes to extract IFN module
t <- dittoHeatmap(mean_exp, genes = genes,  annot.by = c('cluster', 'disease_group'), order.by = c('cluster', 'disease_group'),
                  cluster_rows = T, cluster_cols = T, show_rownames = T, 
                  show_colnames = T, fontsize = 7, data.out = T, cutree_row = 8, highlight.features = c('A1BG', 'ADK', 'AC007952.4', 'AIM2', 'ADAR', 'ADPRHL2', 'BST2', 'CCR7'))
t$clustering_callback <- callback
t <- do.call(pheatmap::pheatmap, t)
sort <- sort(cutree(t$tree_row, k = 8))
gene_list <- list()
for(i in 1:8){
  markers <-  as.vector(names(sort[sort == i]))
  gene_list[[paste0("genes_",i)]] <- markers
}
genes <- c(gene_list$genes_4)
genes
gene_list <- list('genes' = genes)

saveRDS(gene_list, "/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/IFN_genes.RDS")


