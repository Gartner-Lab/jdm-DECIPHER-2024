library(Seurat)
library(tidyverse)
library(dittoSeq)

norm.method = "DSB"
marker.res = 1.4

#args=commandArgs(trailingOnly=T)
#soupx=as.logical(args[[1]]) 
JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm"

#if (soupx){
  IN_DIR=sprintf("%s/merged_jdm_soupx/", JDM_DIR)
  OUT_DIR=sprintf("%s/wnn_v2/", JDM_DIR)
#} else {
#  IN_DIR=sprintf("%s/merged_jdm/", JDM_DIR)
#  OUT_DIR=sprintf("%s/wnn2/", JDM_DIR)
#}


dir.create(OUT_DIR, showWarnings=T)



if (!file.exists(sprintf("%s/wnn_obj.RData", OUT_DIR))){
  
  # load the ADT
  integ_adt = readRDS(sprintf("%s/rpca_v2/integ_ADT_%s.rds", JDM_DIR, norm.method))
  
  # load the harmonized RNA
  load(sprintf("%s/merged_processed.RData", IN_DIR))
  sobj = merged_data
  rm(merged_data)
  gc()
  
  # subset 
  sobj = subset(sobj, cells=colnames(integ_adt))
  #stopifnot(all(colnames(sobj)==colnames(integ_adt)))
  
  # add the ADT
  DefaultAssay(sobj) = "RNA"
  sobj[['ADT']] = NULL
  sobj[["integrated.ADT"]] = integ_adt@assays[["integrated.ADT"]]
  DefaultAssay(sobj) = "integrated.ADT"
  
  sobj[["a.pca"]] = integ_adt@reductions$pca
  sobj[["a.umap"]] = integ_adt@reductions$umap
  
  print("Running WNN")
  
  # run WNN on harmony + a.harmony reduction
  sobj = FindMultiModalNeighbors(
    sobj, reduction.list=list("harmony", "a.pca"),
    dims.list=list(1:30, 1:18), # TODO: should this be different?
    prune.SNN = 1/20, 
    modality.weight.name="RNA.weight"
  )
  
  sobj = RunUMAP(sobj, nn.name = "weighted.nn", 
                 reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  
  save(sobj, file=sprintf("%s/wnn_obj.RData", OUT_DIR))
} else {
  load(file=sprintf("%s/wnn_obj.RData", OUT_DIR))
}


print("running clustering")
for (res in c(0.8, 1.0, 1.2, 1.4, 1.6)){
  
  sobj = FindClusters(sobj, graph.name = "wsnn", 
                      algorithm = 4 , method="igraph", resolution = res, verbose = FALSE)
  sobj@meta.data[[paste0('wsnn_res.', res)]] <- sobj@meta.data$seurat_clusters
  
  png(filename=sprintf('%s/wsnn_%s.png', OUT_DIR, res), 
      width = 5, height = 5, units = "in", 
      res = 300)
  print(DimPlot(sobj, group.by=paste0('wsnn_res.', res), label=T, 
                reduction="wnn.umap") + NoLegend())
  dev.off()
}

save(sobj, file=sprintf("%s/wnn_obj2.RData", OUT_DIR))

Idents(sobj) = sprintf("wsnn_res.%s", marker.res)

#print("finding markers")
#DefaultAssay(sobj) = "integrated.ADT"
#wnn_markers_adt = FindAllMarkers(sobj, only.pos=T)
#save(wnn_markers_adt, file=sprintf("%s/wnn_adt_markers_%s.RData", OUT_DIR, marker.res))

#DefaultAssay(sobj) = "RNA"
#wnn_markers = FindAllMarkers(sobj, test.use="MAST",only.pos=T)
#save(wnn_markers, file=sprintf("%s/wnn_rna_markers_%s.RData", OUT_DIR, marker.res))


