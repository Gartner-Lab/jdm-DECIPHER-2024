library(Seurat)
library(tidyverse)

args = commandArgs(trailingOnly=T)
# fname = args[[1]]
# outdir = args[[2]]
# assay = args[[3]]
# idents = args[[4]]
# cluster_idx = as.numeric(args[[5]])

cluster_idx = as.numeric(args[[1]])
CS_DIR= "/krummellab/data1/erflynn/premier/jdm/data/jdm/"
IN_DIR=sprintf("%s/wnn_v2/reclus", CS_DIR)
OUT_DIR=sprintf("%s/wnn_v2/reclus", CS_DIR)
fname=sprintf("%s/wnn_reclus.RData", OUT_DIR)
outdir=OUT_DIR
assay="RNA"
idents="wsnn_leiden_reclus_res.1.4"
load(fname) # --> sobj
DefaultAssay(sobj) = assay
#if (!idents %in% colnames(sobj@meta.data )){
#  print(sprintf("Error. Ident %s specified is not in metadata", idents))
#  exit()
#}
Idents(sobj) = idents

#nclus = length(unique(sobj@active.ident))
#print(sprintf("Object has %s clusters", nclus))
#if (cluster_idx > (nclus-1)){
#  print(sprintf("Max number of clusters is %s, %s exceeds that", nclus, cluster_idx))
#  exit()
#}
markers = FindMarkers(sobj, ident.1=cluster_idx, test.use="negbinom", only.pos=T, min.pct=0.3)

saveRDS(markers, file=sprintf("%s/markers/%s_%s_markers_%s.RDS", outdir, assay, idents, cluster_idx))






