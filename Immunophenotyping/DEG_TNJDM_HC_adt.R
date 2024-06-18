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
load("/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/sobj_w_raw_adt.RData")
sobj_raw <- sobj

#Adding metadata to objects
sobj_raw <- AddMetaData(sobj_raw, sobj_final@meta.data)

#Setting default assays
DefaultAssay(sobj_raw) <= 'integrated.ADT'

#Setting idents
Idents(sobj_final) <- sobj_final$new_labels
Idents(sobj_raw) <- sobj_raw$new_labels

#Picking labels to perform differential expression analysis for
labels <- levels(sobj_final$new_labels)

#Making lists for differential expression analysis
markers_rna <- list()
m <- 0
n <- 0

#### DESEQ2 5% for ADTs
for(label in labels){
  print(label)
  #Subset object to relevant label
  sobj_subset <- subset(sobj_raw, idents = label)
  
  #Set ident to disease_group
  Idents(sobj_subset) <- sobj_subset$disease_group
  
  #Subset to HC and TNJDM
  sobj_subset <- subset(sobj_subset, idents = c('HC', 'TNJDM'))
  
  #Get counts
  counts <- Seurat::GetAssayData(sobj_subset, assay = 'integrated.ADT', slot = "counts")
  
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
  markers_adt[[label]] <- res
  save(markers_adt, file = "/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/deseq2_adt_5%_final_tnjdm_hc.RData") 
  n <- n+1
  print(n)
}

#Load DeSeq2-results
load(file = "/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/deseq2_adt_5%_final_tnjdm_hc.RData") 

#Make dataframe for differentially expressed ADTs
adt <- as.data.frame(matrix(ncol = 8, nrow = 0))
res <- markers_adt[[1]]
colnames(adt) <- c(colnames(res), 'celltype', 'gene')

#Filter results using an adjusted p-value of < 0.05 and create one collected dataframe for all cells
for(n in 1:length(markers_adt)){
  t <- markers_adt[[n]]
  t <- t %>% filter(padj < 0.05)
  if(nrow(t) < 1){
    next
  } else{
    t$celltype <- names(markers_adt[n])
    t$gene <- rownames(t)
    adt <- rbind(adt, t)
  }}

#Save results
write.xlsx(adt, "/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/markers_adt_tnjdm_hc.xlsx")