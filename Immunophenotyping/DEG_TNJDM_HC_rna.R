#Differential expression
library(DESeq2)
library(tidyverse)
library(Seurat)
library(dittoSeq)
library(scran)
library(glmGamPoi)
library(xlsx)

#Loading data
sobj_final <- readRDS(file = "/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/sobj_final.rds")
load("/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/sobj_w_raw_counts.RData")
sobj_raw_rna <- sobj

#Adding metadata to objects
sobj_raw_rna <- AddMetaData(sobj_raw_rna, sobj_final@meta.data)

#Setting default assays
DefaultAssay(sobj_raw_rna) <- 'RNA'

#Setting idents
Idents(sobj_final) <- sobj_final$new_labels
Idents(sobj_raw_rna) <- sobj_raw_rna$new_labels

#Picking labels to perform differential expression analysis for
labels <- levels(sobj_final$new_labels)

#Making lists for differential expression analysis
markers_rna <- list()
m <- 0
n <- 0

#### DESEQ2 5% for RNA
for(label in labels){
  print(label)
  
  #Subset object to relevant label
  sobj_subset <- subset(sobj_raw_rna, idents = label)
  
  #Set ident to disease_group
  Idents(sobj_subset) <- sobj_subset$disease_group
  
  #Subset to HC and TNJDM
  sobj_subset <- subset(sobj_subset, idents = c('HC', 'TNJDM'))
  
  #Get counts
  counts <- Seurat::GetAssayData(sobj_subset, assay = 'RNA', slot = "counts")
  
  #Get relevant meta data
  metadata <- sobj_subset@meta.data[c('donor', 'disease_group', 'well')]
  
  # Trim counts to just genes that meet the min.pct cutoff
  counts.1 <- counts[,metadata[,'disease_group'] %in% 'HC']
  counts.2 <- counts[,metadata[,'disease_group'] %in% 'TNJDM']
  pct.1 <- apply(counts.1, 1, function(x) mean(x>0))
  pct.2 <- apply(counts.2, 1, function(x) mean(x>0))
  keep <- pct.1 >= 0.05 | pct.2 >= 0.05
  
  counts <- counts[keep, metadata[,'disease_group'] %in% c('HC', 'TNJDM')]

  # Create DESeq2 data set and add it to the DESeq.vars
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~well + disease_group)
  
  #Set HC as reference
  dds$disease_group <- relevel(dds$disease_group, ref = "HC")
  
  #Compute sumfactors
  dds <- computeSumFactors(dds)
  
  #Run  DESeq using suggested parameters for single cell data
  dds <- DESeq(dds,test = 'LRT', fitType = 'glmGamPoi', 
               useT = T, minmu = 1e-6, minReplicatesForReplace = Inf, reduced = ~well)
  
  #Save results
  res <- as.data.frame(results(dds))
  markers_rna[[label]] <- res
  save(markers_rna, file = "/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/deseq2_rna_5%_tnjdm_hc.RData") 
  n <- n+1
  print(n)
}

#Load DeSeq2/results
load(file = "/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/deseq2_rna_5%_final_tnjdm_hc.RData") 

#Make dataframe for differentially expressed ADTs
rna <- as.data.frame(matrix(ncol = 8, nrow = 0))
res <- markers_rna[[1]]
colnames(rna) <- c(colnames(res), 'celltype', 'gene')

#Filter results using an adjusted p-value of < 0.05 and a Logfold change of 1 and create one collected dataframe for all cells
for(n in 1:length(markers_rna)){
  t <- markers_rna[[n]]
  t <- t %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
  t$celltype <- names(markers_rna[n])
  t$gene <- rownames(t)
  rna <- rbind(rna, t)
}

#Save results
write.xlsx(rna, "/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/markers_rna_tnjdm_hc_min_lfc_1.xlsx")
