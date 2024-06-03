library(Seurat)

setwd(JDM_DIR)
load("wnn/reclus/sobj_w_lab.RData")

filtered = Read10X("raw/filtered_feature_bc_matrix/")

# check cell names, genes
gex = filtered$`Gene Expression`
# make sure all cell names, genes are PRESENT
cells = colnames(sobj)
genes = rownames(sobj@assays$RNA)
stopifnot(length(setdiff(cells, colnames(gex)))==0)
stopifnot(length(setdiff(genes, rownames(gex)))==0)

sobj@assays$RNA@counts = gex[genes, cells]
save(sobj, file="wnn/reclus/sobj_w_raw_counts.RData")


# check ordering

# add back in and save