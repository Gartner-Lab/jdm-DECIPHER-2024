

library(tidyverse)
library(Seurat)

CODE_DIR="/krummellab/data1/erflynn/premier/jdm/"
JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm"

#source(sprintf("%s/scripts/utils/label_transfer.R", CODE_DIR))
#source(sprintf("%s/scripts/utils/plot_wnn.R", CODE_DIR))


IN_DIR=sprintf("%s/wnn_v2/", JDM_DIR)
OUT_DIR=sprintf("%s/reclus", IN_DIR)
dir.create(OUT_DIR)
load(sprintf("%s/sobj_w_lab.RData", IN_DIR))
unique(sobj@meta.data$wnn_cluster_label)
meta_keep = sobj@meta.data %>%
  as_tibble(rownames="cell_id") %>%
  filter(wnn_cluster_label != "T4_Naive* 1" &
           wnn_cluster_label != "T4_Naive* 2") 

meta_keep
sobj = subset(sobj, cells=meta_keep$cell_id)

sobj = FindMultiModalNeighbors(
  sobj, reduction.list=list("harmony", "a.pca"),
  dims.list=list(1:30, 1:18), # TODO: should this be different?
  prune.SNN = 1/20, 
  modality.weight.name="RNA.weight"
)

sobj = RunUMAP(sobj, nn.name = "weighted.nn", 
               reduction.name = "wnn.umap", 
               reduction.key = "wnnUMAP_")


for (res in c( 1.0, 1.2, 1.4, 1.6)){
  
  sobj = FindClusters(sobj, graph.name = "wsnn", 
                      algorithm = 4 , resolution = res, verbose = FALSE,
                      method="igraph")
  sobj@meta.data[[paste0('wsnn_leiden_reclus_res.', res)]] <- sobj@meta.data$seurat_clusters
  
  png(filename=sprintf('%s/wsnn_leiden_%s.png', OUT_DIR, res), 
      width = 5, height = 5, units = "in", 
      res = 300)
  print(DimPlot(sobj, group.by=paste0('wsnn_leiden_reclus_res.', res), label=T, 
                reduction="wnn.umap") + NoLegend())
  dev.off()
}

save(sobj, file=sprintf("%s/wnn_reclus.RData", OUT_DIR))

