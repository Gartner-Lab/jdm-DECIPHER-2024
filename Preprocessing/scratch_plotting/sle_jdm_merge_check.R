library(Seurat)
library(tidyverse)
library(dittoSeq)
library(cowplot)
library(harmony)


DATA_DIR = "/krummellab/data1/erflynn/premier/jdm/data/"
jdm_data = readRDS(sprintf("%s/jdm/analyzed_seurat_object (1).rds", DATA_DIR))

load(sprintf("%s/sobj_GSE135779/merged_df_proc2/GSE135779_merged_processed.RData", DATA_DIR)) 
# --> sobj
csle_data = sobj
rm(sobj)

# check overlap in genes

c_genes = rownames(csle_data@assays[["RNA"]]) # 20,903
j_genes = rownames(jdm_data@assays[["RNA"]]) # 18,465 ß
comb_genes = intersect(c_genes, j_genes) # 13,411 genes intersect
# ... not great but not terrible?

jdm_data@reductions


