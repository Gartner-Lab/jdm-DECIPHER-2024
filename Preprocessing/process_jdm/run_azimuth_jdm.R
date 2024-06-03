library(Seurat)
library(tidyverse)
library(dittoSeq)
library(Azimuth)

reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")

load("/krummellab/data1/erflynn/premier/jdm/data/jdm/merged_jdm_soupx/merged_processed.RData")
query=merged_data
rm(merged_data); gc()

# Preprocess with SCTransform
query <- SCTransform(
  object = query,
  assay = "RNA",
  new.assay.name = "refAssay",
  residual.features = rownames(x = reference$map),
  reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
  method = 'glmGamPoi',
  ncells = 2000,
  n_genes = 2000,
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)

# Find anchors between query and reference
anchors <- FindTransferAnchors(
  reference = reference$map,
  query = query,
  k.filter = NA,
  reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "refAssay",
  query.assay = "refAssay",
  reference.reduction = "refDR",
  normalization.method = "SCT",
  features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
  dims = 1:50,
  n.trees = 20,
  mapping.score.k = 100
)

# Transfer cell type labels and impute protein expression
#
# Transferred labels are in metadata columns named "predicted.*"
# The maximum prediction score is in a metadata column named "predicted.*.score"
# The prediction scores for each class are in an assay named "prediction.score.*"
# The imputed assay is named "impADT" if computed

refdata <- lapply(X = "celltype.l2", function(x) {
  reference$map[[x, drop = TRUE]]
})
names(x = refdata) <- "celltype.l2"
if (TRUE) {
  refdata[["impADT"]] <- GetAssayData(
    object = reference$map[['ADT']],
    slot = 'data'
  )
}
query <- TransferData(
  reference = reference$map,
  query = query,
  dims = 1:50,
  anchorset = anchors,
  refdata = refdata,
  n.trees = 20,
  store.weights = TRUE
)

# Calculate the embeddings of the query data on the reference SPCA
query <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = reference$map,
  query = query,
  reductions = "pcaproject",
  reuse.weights.matrix = TRUE
)

# Calculate the query neighbors in the reference
# with respect to the integrated embeddings
query[["query_ref.nn"]] <- FindNeighbors(
  object = Embeddings(reference$map[["refDR"]]),
  query = Embeddings(query[["integrated_dr"]]),
  return.neighbor = TRUE,
  l2.norm = TRUE
)

# The reference used in the app is downsampled compared to the reference on which
# the UMAP model was computed. This step, using the helper function NNTransform,
# corrects the Neighbors to account for the downsampling.
query <- Azimuth:::NNTransform(
  object = query,
  meta.data = reference$map[[]]
)

# Project the query to the reference UMAP.
query[["proj.umap"]] <- RunUMAP(
  object = query[["query_ref.nn"]],
  reduction.model = reference$map[["refUMAP"]],
  reduction.key = 'UMAP_'
)


# Calculate mapping score and add to metadata
query <- AddMetaData(
  object = query,
  metadata = MappingScore(anchors = anchors),
  col.name = "mapping.score"
)
azimuth_meta = query@meta.data
JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm/"
save(azimuth_meta, file=sprintf("%s/azimuth_jdm.RData", JDM_DIR))
