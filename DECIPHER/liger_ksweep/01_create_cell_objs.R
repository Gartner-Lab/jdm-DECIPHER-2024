# Liger workflow, adapted for use on C4

#### SET UP LIGER OBJECTS ####

# Create Liger objects on server

# CHANGE THIS TO OWN DIRECTORY

JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm"
setwd(JDM_DIR)
library(rliger)
library(Matrix)
library(Seurat)
library(matrixStats)
library(parallel)
dir.create('k_sweep')
dir.create('k_sweep/Liger_objects')

## load the data (make it a diet seurat object)
#RMs.seu = readRDS("wnn_diet_obj.RDS")
RMs.seu = readRDS("wnn_v2/diet_sobj.RDS")
DefaultAssay(RMs.seu)="RNA"
#if ("integrated" %in% Assays(RMs.seu)){
#  RMs.seu[['integrated']] <- NULL
#}
  

RMs.seu = AddMetaData(RMs.seu, RMs.seu$broad_ids, "type_spec")

# Identify cell types with at least 10 samples with greater than 100 cells
# Only run LIGER/network analyses on these cell types, since others won't have enough statistical power
cell.types=names(which(colSums(table(RMs.seu@meta.data[,c("study_id_visit","type_spec")])>100)>10))
save(cell.types, file="k_sweep/Liger_objects/cell.types.robj")
Idents(RMs.seu) <- RMs.seu$type_spec


# Determine minibatch sizes - set to minimium dataset size if smaller than default
# online iNMF uses "minibatches" to process large datasets more quickly, see the documentation in rliger package
# if you haven't read in a separate metadata object, need to create (metadata = seurat.object@meta.data)
metadata = RMs.seu@meta.data
saveRDS(metadata, file = "RMs.metadata.rds")
names(cell.types) = cell.types
minibatch_size = lapply(cell.types, function(x) {
  size = min(table(subset(metadata, type_spec==x)$well))
  if (size > 1000) {
    size = floor(size/1000)*1000
  } else if (size > 500) {
    size = floor(size/500)*500
  } else if (size > 100){
    size = floor(size/100)*100
  } else size = NA 
})

saveRDS(minibatch_size, 
        file = "k_sweep/Liger_objects/minibatch_size.rds", version = 2)


# Create basic Liger objects for each cell type
mclapply(cell.types, function(x){
  seu_sub <- subset(RMs.seu, type_spec==x)
  if (length(unique(RMs.seu@meta.data$well))>1){
    Batch.list <- SplitObject(seu_sub, split.by="well")
    Liger.list=list()
    for (j in 1:length(Batch.list)){
      Liger.list[[j]]=Batch.list[[j]]@assays$RNA@counts
    }
    names(Liger.list)=names(Batch.list)
    Liger <- createLiger(Liger.list)
  } else {
    Liger <- createLiger(seu_sub@assays$RNA@counts)
  }
  # normalize so all cells have same total counts
  Liger <- rliger::normalize(Liger)
  Liger <- selectGenes(Liger)
  # log normalize - this normalization worked best on my dataset, but the original LIGER paper doesn't log normalize
  # tried out both with some simulated data and felt like log normalization was best, but probably worth re-doing this for a couple different datasets to test
  for (k in 1:length(Liger@norm.data)){
    Liger@norm.data[[k]]=log1p(Liger@norm.data[[k]]*10000)
  }
  Liger <- rliger::scaleNotCenter(Liger)
  saveRDS(Liger, paste("k_sweep/Liger_objects/Liger", x, "rds", sep="."))
  return(NULL)
}, mc.cores = min(detectCores(), length(cell.types)))



