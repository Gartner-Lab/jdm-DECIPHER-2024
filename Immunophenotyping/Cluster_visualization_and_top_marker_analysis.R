#Loading packages
library(Seurat)
library(tidyverse)
library(dittoSeq)
library(scales)
library(cowplot)
library(openxlsx)
library(MAST)
library(Seurat)
library(tidyverse)
library(xlsx)
library(scales)
library(vctrs)

#Laoding data
load("/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/sobj_w_raw_counts.RData")
sobj_final <- readRDS(file = "/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/sobj_final.rds")

#Figure 2A
th <-   theme(axis.text.x = element_text(color = "black", size = 18, face = 'bold'), 
              strip.text = element_text(color = "black", size = 10),
              axis.text.y.left = element_text(size = 18, color = 'black', face = 'bold'),
              plot.title = element_text(size = 17, face = 'bold'),
              legend.title = element_text(size = 18, face = 'bold'),
              plot.subtitle = element_text(size = 16, face = 'bold'),
              legend.text = element_text(size = 18))
dittoDimPlot(sobj_final, var = "new_labels", reduction.use = "wnn.umap", size = 0.3, do.label = TRUE, labels.size = 4, legend.show = T) +
  guides(colour=guide_legend(ncol =6, override.aes = list(size=7))) +
  ggtitle('') +
  theme(plot.caption = element_text(hjust = 0, size = 10),
        axis.title = element_text(size = 12, face = 'bold'),
        legend.position = 'none',
        axis.text = element_text(size = 12))

### Figure 2B
#Changing name
sobj_raw_rna <- sobj

#Adding meta-data
sobj_raw_rna <- AddMetaData(sobj_raw_rna, sobj_final@meta.data)

#Setting idents
Idents(sobj_raw_rna) <- sobj_raw_rna$new_labels

#Setting default assay
DefaultAssay(sobj_raw_rna) <- 'RNA'

#Running FindAllMarkers between cell clusters
markers <- FindAllMarkers(object = sobj_raw_rna, 
                          only.pos = T, 
                          verbose = T, 
                          test.use = 'MAST')

#Select top2 markers pr cell type based on LFC
top <- markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) %>%
  select(gene) %>%
  ungroup()
genes <- unique(top$gene)

#Setting idents and assay
Idents(sobj_final) <- 'new_labels'
DefaultAssay(sobj_final) <- 'RNA'

#Pseudobulking
sobj_avg <- AverageExpression(sobj_final, assays = 'RNA', return.seurat = T)
sobj_avg$new_labels <- Idents(sobj_avg)

#Making heatmap
DoHeatmap(sobj_avg, features = genes, draw.lines = F, size = 3) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red')