# adds additional filters to the data
# - nCount_ADT upper bound
# - gene filter for genes present in at least three cells

library(Seurat)
library(tidyverse)
source(sprintf("%s/scripts/utils/generate_profile_plot.R", CODE_DIR))


jdm_obj = readRDS( sprintf("%s/seurat_filt_postdemuxDF_w_meta.rds", JDM_DIR))

## all singlets - already filtered
# table(jdm_obj@meta.data$DROPLET.TYPE)

# min count of 3 cells with gene to keep genes
count_mat = jdm_obj@assays[["RNA"]]@counts
nnz = tabulate(count_mat@i + 1)
keep.features = rownames(count_mat)[nnz>=3]
adt.features = rownames(jdm_obj@assays$ADT@counts)

print(sprintf("Started with %s genes, filtered to at least 3 counts per gene --> %s genes", nrow(count_mat), length(keep.features)))

rm(count_mat)
gc()
jdm_obj = subset(jdm_obj, features=c(keep.features, adt.features))

# now filter based on metadata QC features
jdm_obj@meta.data %>% 
  as_tibble(rownames="cell_id") %>%
  select(nCount_RNA:nFeature_ADT, percent.mt, percent.ribo) %>%
  summary()

generate_all_profile_plots(sprintf("%s/seurat_obj", JDM_DIR))

cells.keep = jdm_obj@meta.data %>%
  as_tibble(rownames="cell_id") %>%
  filter(nCount_ADT <= 5000 & percent.ribo <= 60 & percent.mt <= 15) %>%
  pull(cell_id)

print(sprintf("Started with %s cells, now have %s", nrow(jdm_obj@meta.data), length(cells.keep)))

jdm_obj = subset(jdm_obj, cells=cells.keep )
saveRDS(jdm_obj, sprintf("%s/seurat_filt_postdemuxDF_w_meta_filt2.rds", JDM_DIR), compress=F)
