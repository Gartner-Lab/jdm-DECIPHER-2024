#Loading packages
library(Seurat)
library(tidyverse)
library(dittoSeq)
library(scales)
library(cowplot)
library(openxlsx)
library(MAST)
library(Seurat)
library(tidyverse)
library(xlsx)
library(scales)
library(vctrs)

#Laoding data
sobj_final <- readRDS(file = "/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/sobj_final.rds")
markers <- xlsx::read.xlsx("/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/markers_adt_tnjdm_hc.xlsx", sheetIndex = 1)

#Fig 2D
#Filter markers to seleceted cell types and more stringent LFC and p-value cutoffs
markers <- filter(markers, abs(log2FoldChange) > 0.5 & padj < 0.05 & celltype %in% c('CD4+ eff', 'CD4+ Tregs', 'B_naive1'))

#Combine disease state and cluster
Idents(sobj_final) <- sobj_final$disease_group
sobj <- subset(sobj_final$new_labels %in% c('CD4+ Tregs', 'CD4+ eff', 'B_naive1'))
disease_groups <- sobj$new_diseasegroups
cluster <- sobj$new_labels
new <- data.frame(cluster, disease_groups)
new$new <- paste0(new$cluster, "&", new$disease_groups)
cluster_group <- new$new
Idents(sobj) <- cluster_group
sobj[["cluster_group"]] <- Idents(sobj)

#Pseudobulk data
mean_exp <- AverageExpression(sobj, assays = "integrated.ADT", return.seurat = T, verbose = T)
mean_exp[["disease_group"]] <- Idents(mean_exp)
mean_exp[["cluster"]] <- Idents(mean_exp)
mean_exp[['disease_group_cluster']] <- Idents(mean_exp)
disease_group <- mean_exp$disease_group
levels(disease_group) <- sub("(.*)&","", levels(disease_group))
mean_exp$disease_group <- disease_group
Idents(mean_exp) <- mean_exp$disease_group
levels(mean_exp) <- c('TNJDM','Active', 'Inactive', 'HC')
mean_exp$disease_group <- Idents(mean_exp)

cluster <- mean_exp$cluster
levels(cluster) <- sub("&[^_]+$","", levels(cluster))
mean_exp$cluster <- cluster
Idents(mean_exp) <- mean_exp$cluster
levels <- c('CD4+ Tregs', 'CD4+ eff', 'B_naive1')
levels(mean_exp) <- levels
mean_exp$cluster <- Idents(mean_exp)

#To split plot by cell type visually
split <- as.factor(levels)
split <- factor(split, levels = levels)
split <- vec_rep_each(split, 4)

#Select ADTs to show
adts <- c('MICA-MICB', 'CD1c', 'CD268--BAFF-R', 'CD274--B7-H1-PD-L1', 'CD366--Tim-3', 
          'CD278--ICOS', 'CD164', 'CD38', 'CD101--BB27', 'KLRG1--MAFA', 'CD279--PD-1')

#Draw heatmap
dittoHeatmap(mean_exp, genes = adts, show_colnames = F, 
             annot.by = c('disease_group'), order.by = c('cluster', 'disease_group'),
             show_rownames = T, complex = T, column_split = split, cluster_column_slices = F,
             fontsize = 8, row_names_max_width = unit(12, 'cm'), 
             annot.colors = c('#F3756D','#A680BA', '#7BAF41','#88CDEA'))

