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

#Fig3c-e
#Create data frame for CD14+ mono with information on CD169, IFN-score and VAS global
Idents(sobj_final) <- sobj_final$cluster_group
mean_exp <- AverageExpression(sobj_final, assays = "integrated.ADT", return.seurat = T, verbose = T)
meta <- mean_exp@meta.data
meta$vasglobal <- NA
meta$patient <- rownames(meta)
meta$patient <- sub("(.*)&","", meta$patient)
meta$cluster <- sub("&.*","", rownames(meta))
for (i in 1:nrow(meta)) {
  meta$vasglobal[i] <- sobj_final$vasglobal[match(meta$patient[i], sobj_final$paper_id)]
}

expression <- GetAssayData(mean_exp, assay = 'integrated.ADT')
expression <- t(expression)
expression <- as.data.frame(expression)
expression <-  dplyr::select(expression, 'CD169--Sialoadhesin-Siglec-1')
expression$disease_group <- rownames(expression)
expression$patient <- sub("(.*)&","", rownames(expression))
expression$cluster <- sub("&.*","", rownames(expression))

collected <- left_join(meta, expression)
collected <- na.omit(collected)

mean_exp <- AverageExpression(sobj_final, assays = "RNA", return.seurat = T, verbose = T)
gene_list <- readRDS("/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/IFN_genes.RDS")
mean_exp <- AddModuleScore(mean_exp, features = gene_list)
meta <- mean_exp@meta.data
meta$disease_group <- rownames(meta)
meta$patient <- sub("(.*)&","", rownames(meta))
meta$cluster <- sub("&.*","", rownames(meta))

collected <- left_join(collected, meta)
collected <- na.omit(collected)

collected$disease_group <- NA
for (i in 1:nrow(collected)) {
  collected$disease_group[i] <- sobj_final$new_diseasegroups[match(collected$patient[i], sobj_final$paper_id)]
}
mono <- filter(collected, cluster == 'CD14+ mono')


#Fig3c
ggplot(mono, aes(`CD169--Sialoadhesin-Siglec-1`, vasglobal)) +
  geom_point(aes(color = disease_group)) +
  geom_smooth(method=lm, color = 'black') + 
  theme_bw() +
  stat_cor(method = 'spearman') +
  xlab('CD169 - Siglec-1') +
  ylab('VAS global') +
  ylim(0,10) +
  scale_color_manual(values = c(
    'Inactive' = 'gold1',
    'TNJDM' = 'red','Active' = 'orange'))  +
  theme(legend.position = 'none')


#Fig3d
ggplot(mono, aes(Cluster1, vasglobal)) +
  geom_point(aes(color = disease_group)) +
  geom_smooth(method=lm, color = 'black') + 
  theme_bw() +
  stat_cor(method = 'spearman') +
  xlab('IFN score') +
  ylab('VAS global') + ylim(0,10) +
  scale_color_manual(values = c(
    'Inactive' = 'gold1',
    'TNJDM' = 'red','Active' = 'orange'))  +
  theme(legend.position = 'none')

#Fig3e
ggplot(mono, aes(`CD169--Sialoadhesin-Siglec-1`, Cluster1)) +
  geom_point(aes(color = disease_group)) +
  geom_smooth(method=lm, color = 'black') + 
  theme_bw() +
  stat_cor(method = 'spearman') +
  xlab('CD169 - Siglec-1') +
  ylab('IFN score') +
  scale_color_manual(values = c(
    'Inactive' = 'gold1',
    'TNJDM' = 'red','Active' = 'orange')) +
  theme(legend.position = 'none')
