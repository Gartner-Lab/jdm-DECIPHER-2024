# based on this tutorial https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html

library('SoupX')
library('Seurat')
library('tidyverse')

args = commandArgs(trailingOnly=T)
idx = as.numeric(args[1])
raw_well = Read10X(sprintf("%s/raw/crosslong_raw/jdm_crosslong%s/outs/raw_feature_bc_matrix/", JDM_DIR, idx))
filtered = Read10X(sprintf("%s/raw/filtered_feature_bc_matrix/", JDM_DIR))

raw_rna_well = raw_well$`Gene Expression`
filtered_rna = filtered$`Gene Expression`

# double check on row IDs
stopifnot(all(rownames(filtered_rna) == rownames(raw_rna_well)))

# fix and check the column names 
colnames(raw_rna_well) = str_replace_all(colnames(raw_rna_well), "-1", sprintf("-%s", idx))
filtered_rna_well = filtered_rna[,str_detect(colnames(filtered_rna), sprintf("-%s", idx))]

stopifnot(length(setdiff(colnames(filtered_rna_well), colnames(raw_rna_well)))==0)

sc = SoupChannel(raw_rna_well, filtered_rna_well)

### note you want to read in your seurat object here
# instead of running this section - this was just to sanity check it works
sobj <- CreateSeuratObject(counts = filtered_rna_well,
                           project = "well1")
sobj = NormalizeData(sobj) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims=1:30)
sobj = FindNeighbors(sobj)
sobj = FindClusters(sobj, res=0.8)
###

# double check that there is overlap
stopifnot(length(setdiff(colnames(sc$toc), rownames(sobj@meta.data)))==0)
sobj = subset(sobj, cells=colnames(sc$toc))

# fill this in with your clusters
sc = setClusters(sc, setNames(sobj@meta.data$seurat_clusters,
                              rownames(sobj@meta.data)))
sc = setDR(sc, sobj@reductions$umap@cell.embeddings[colnames(sc$toc),])
sc = autoEstCont(sc)
out = adjustCounts(sc)
sobj_adj = CreateSeuratObject(out)
save(sc, file=sprintf("%s/soupx/soupx_out_well%s.RData", JDM_DIR, idx))

# load our filtered JDM object
jdm_obj = readRDS(sprintf("%s/seurat_filt_postdemuxDF_w_meta_filt2.rds", JDM_DIR))
DefaultAssay(jdm_obj) = "RNA"
jdm_well = subset(jdm_obj, well==idx)

# subset the soupx object to the same cells and features as the filtered JDM object
stopifnot(length(setdiff(colnames(jdm_well), colnames(sobj_adj))) == 0)
stopifnot(length(setdiff(rownames(jdm_well), rownames(sobj_adj))) == 0)

sobj_adj = subset(sobj_adj, cells=colnames(jdm_well))
sobj_adj = subset(sobj_adj, features=rownames(jdm_well))

stopifnot(all(colnames(jdm_well)==colnames(sobj_adj)))
stopifnot(all(rownames(jdm_well)==rownames(sobj_adj)))

# replace the original jdm object counts with the adjusted soupx counts
stopifnot(all(sobj_adj@assays[["RNA"]]@counts==sobj_adj@assays[["RNA"]]@data))
jdm_well@assays[["RNA"]]@counts=sobj_adj@assays[["RNA"]]@counts
jdm_well@assays[["RNA"]]@data=sobj_adj@assays[["RNA"]]@counts

saveRDS(jdm_well, file=sprintf("%s/soupx/soupx_sobj_adj_well%s.RDS", JDM_DIR, idx), compress=F)



# head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
# cntSoggy = rowSums(sc$toc > 0)
# cntStrained = rowSums(out > 0)
# mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
# mostZeroed
# 
# 
# plotChangeMap(sc, out, "GNLY")+
#   ggtitle("GNLY change in expression")
