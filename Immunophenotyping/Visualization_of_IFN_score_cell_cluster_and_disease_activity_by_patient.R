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

#Load data
sobj_final <- readRDS(file = "/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/sobj_final.rds")

#Fig3a
#Pseudubulk data by study_id and cluster
studyid <- sobj_final$paper_id
cluster <- sobj_final$new_labels
new <- data.frame(cluster, studyid)
new$new <- paste0(new$cluster, "&", new$studyid)
cluster_studyid <- new$new
Idents(sobj_final) <- cluster_studyid
sobj_final[["cluster_studyid"]] <- Idents(sobj_final)
Idents(sobj_final) <- sobj_final$cluster_studyid
mean_exp <- AverageExpression(sobj_final, assays = "RNA", return.seurat = T, verbose = T)

#Add IFN score
gene_list <- readRDS("/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/IFN_genes.RDS")
score <- AddModuleScore(mean_exp, features = gene_list)

#Create dataframe with study_id, diseasegroup, cluster and IFN score, filtering out clusters with few cells
meta <- score@meta.data
meta$disease_group <- NA
meta$patient <- rownames(meta)
meta$patient <- sub("(.*)&","", meta$patient)
meta$cluster <- sub("&.*","", rownames(meta))
for (i in 1:nrow(meta)) {
  meta$disease_group[i] <- sobj_final$new_diseasegroups[match(meta$patient[i], sobj_final$paper_id)]
}
mat <- meta[,c('disease_group', 'Cluster1', 'cluster', 'patient')]
mat <- filter(mat, ! cluster %in% c('B_naive4', 'PBs', 'CD8+ mem', 'CD56bright NK', 'CDCs'))
mat$new <- paste0(mat$patient, "&", mat$disease_group)
mat <- mat[,c('new', 'Cluster1', 'cluster')]
mat <- pivot_wider(mat, names_from = cluster, values_from = Cluster1)
rownames <- mat$new
mat <- mat[,c(2:17)]
mat <- data.matrix(mat)
rownames <- sub("(.*)&","", rownames)
rownames(mat) <- rownames
mat <- t(mat)

#Add additional parameters to heatmap
col_fun = colorRamp2(c(0, 10), c("white", "darkgoldenrod1"))
meta <- sobj_final@meta.data %>%
  group_by(paper_id) %>%
  filter(row_number() == 1) %>%
  select(vasglobal, paper_id)
vas <- meta$vasglobal
names(vas) <- meta$paper_id
column_ha = HeatmapAnnotation(t = rownames, score = vas,
                              col = list(t = c("TNJDM" = "#F3756D", "Active" = "#A680BA", "Inactive" = "#7BAF42", 'HC' = '#88CDEA'),
                                         score = col_fun),
                              na_col = "black")

#Draw heatmap
t <- Heatmap(mat, top_annotation = column_ha)
draw(t)


#Fig3b
#Pseudobulk data pr study_id
Idents(sobj_final) <- sobj_final$paper_id
mean_exp <- AverageExpression(sobj_final, assays = "RNA", return.seurat = T, verbose = T)

#Add IFN score
score <- AddModuleScore(mean_exp, features = gene_list)

#Create dataframe with study_id, disease group, vas global and ifn score
meta1 <- score@meta.data
meta1$paper_id <- row.names(meta1)
meta2 <- sobj_final@meta.data   %>%
  group_by(paper_id) %>%
  filter(row_number() == 1) %>%
  select(vasglobal, new_diseasegroups, paper_id) %>%
  na.omit()
t <- left_join(meta2, meta1)

ggplot(t, aes(Cluster1, vasglobal)) +
  geom_smooth(method=lm) + 
  theme_bw() +
  stat_cor(method = 'spearman') +
  geom_point(aes(color = new_diseasegroups)) +
  ylim(0,10) +
  scale_color_manual(values = c('Inactive' = 'gold1',
                                'TNJDM' = 'red','Active' = 'orange'))
