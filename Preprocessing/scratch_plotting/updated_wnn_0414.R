library(Seurat)
library(tidyverse)
library(dittoSeq)
JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm"
norm.method="DSB"

OUT_DIR = sprintf("%s/wnn_out_rpca_%s/", JDM_DIR, norm.method)

load(file=sprintf("%s/wnn_obj.RData", OUT_DIR))
Idents(sobj) = "wsnn_res.1.4"

VlnPlot(sobj, features = "RNA.weight", 
        group.by = 'wsnn_res.1.4', sort = TRUE, pt.size = 0) +
  NoLegend()
ggsave(sprintf("%s/wnn_weight1.4.png", OUT_DIR))


## add labels
## label the clusters
jdm_obj_initial= readRDS(sprintf("%s/analyzed_seurat_object (1).rds", JDM_DIR))
meta = jdm_obj_initial@meta.data %>% 
  as_tibble(rownames="cell_id")
rm(jdm_obj_initial)
length(unique(sobj@meta.data$wsnn_res.1.4)) # 36
clus_lab = sobj@meta.data %>% 
  as_tibble(rownames="cell_id") %>%
  dplyr::select(cell_id, wsnn_res.1.4) %>%
  left_join(  
    meta %>% 
      dplyr::select(cell_id, manual_labels))
length(unique(clus_lab$wsnn_res.1.4))
clus_to_lab = clus_lab %>%
  dplyr::rename(cluster="wsnn_res.1.4") %>%
  group_by(cluster) %>%
  mutate(tot=n()) %>%
  group_by(cluster, manual_labels, tot) %>%
  dplyr::count() %>%
  arrange(cluster, desc(n)) %>%
  mutate(frac=n/tot)
length(unique(clus_to_lab$cluster)) # 36

clus_to_lab %>%
  ungroup() %>%
  group_by(cluster) %>%
  filter(frac > 0.4) %>%
  View()



clus_to_lab2 = clus_to_lab %>%
  ungroup() %>%
  filter(tot > 10) %>%
  group_by(cluster) %>%
  slice_max(1)  # 22

clus_to_lab3 = clus_to_lab2 %>%
  mutate(manual_labels=as.character(manual_labels)) %>%
  mutate(assigned_label=ifelse(frac < 0.4 | is.na(manual_labels), "unknown", manual_labels)) %>%
  arrange(assigned_label) %>%
  ungroup() %>%
  group_by(assigned_label) %>%
  mutate(idx=1:n(),
         nlab=n()) %>%
  mutate(cluster_label=ifelse(nlab ==1, assigned_label,
                              paste(assigned_label, idx))) %>%
  ungroup()
clus_map = clus_to_lab3 %>% distinct(cluster, cluster_label)

new_data = sobj@meta.data %>% 
  as_tibble(rownames="cell_id") %>%
  select(cell_id, wsnn_res.1.4) %>%
  left_join(clus_to_lab3 %>% select(cluster, cluster_label) %>%
              mutate(cluster_label=as.factor(cluster_label)), by=c("wsnn_res.1.4"="cluster"))

manual_lab_df = clus_lab %>% 
  mutate(manual_labels=as.character(manual_labels),
         manual_labels=ifelse(is.na(manual_labels), "unknown",
                              manual_labels)) 

sobj = AddMetaData(sobj, new_data$cluster_label, col.name="cluster_label")
sobj = AddMetaData(sobj, manual_lab_df$manual_labels, col.name="manual_labels")


col_vec = dittoColors()[1:36]
names(col_vec)= unique(sobj@meta.data$cluster_label)

#TODO: fix color vec for manual labels


Idents(sobj) = "cluster_label"
## remake the plots
DimPlot(sobj, reduction = 'wnn.umap', 
        label = TRUE, cols=col_vec,
        label.size=3) + NoLegend()
ggsave(sprintf("%s/wnn_umap_clus_lab.pdf", OUT_DIR), height=5, width=5)


DimPlot(sobj, reduction = 'wnn.umap', group.by="manual_labels",
        label = TRUE, cols=dittoColors(),
        label.size=3) 
ggsave(sprintf("%s/wnn_umap_manual_lab.pdf", OUT_DIR), height=5,
       width=9)



DimPlot(sobj, reduction = 'umap', 
        label = TRUE, cols=col_vec,
        label.size=3) + NoLegend()+
  ylab("rnaUMAP_2")+
  xlab("rnaUMAP_1")
ggsave(sprintf("%s/wnn_umap_RNA.pdf", OUT_DIR), height=5, width=5)
DimPlot(sobj, reduction = 'umap', group.by = "manual_labels",
        label = TRUE, cols=dittoColors(),
        label.size=3)+
  ylab("rnaUMAP_2")+
  xlab("rnaUMAP_1")
ggsave(sprintf("%s/wnn_umap_manual_RNA.pdf", OUT_DIR), height=5, width=5)


DimPlot(sobj, reduction = 'a.umap', 
        label = TRUE, cols=col_vec,
        label.size=3) + NoLegend()
ggsave(sprintf("%s/wnn_umap_ADT.pdf", OUT_DIR), height=5, width=5)


DimPlot(sobj, reduction = 'a.umap', group.by = "manual_labels",
        label = TRUE, cols=dittoColors(),
        label.size=3)
ggsave(sprintf("%s/wnn_umap_manual_ADT.pdf", OUT_DIR), height=5, width=5)


VlnPlot(sobj, features = "RNA.weight", cols=col_vec, 
        group.by = 'cluster_label', sort = TRUE, pt.size = 0) +
  NoLegend()+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
        axis.title.x = element_blank())
ggsave(sprintf("%s/wnn_vln_clus_lab.png", OUT_DIR), height=5, width=8)


### check on doublets
DefaultAssay(sobj) = "integrated.ADT"
FeaturePlot(sobj, features=c("CD4", "CD8", "CD16"), reduction = 'wnn.umap')
p1 = VlnPlot(sobj, features=c("CD4"), cols=col_vec, 
        group.by = 'cluster_label', sort = FALSE, pt.size = 0) +
  NoLegend()+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
        axis.title.x = element_blank())
p2 = VlnPlot(sobj, features=c("CD8"), cols=col_vec, 
        group.by = 'cluster_label', sort = FALSE, pt.size = 0) +
  NoLegend()+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
        axis.title.x = element_blank())
p3 = VlnPlot(sobj, features=c("CD16"), cols=col_vec, 
        group.by = 'cluster_label', sort = FALSE, pt.size = 0) +
  NoLegend()+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
        axis.title.x = element_blank())

cowplot::plot_grid(p1, p2, p3, ncol=1)
ggsave(sprintf("%s/cd_doublet_check.png", OUT_DIR), height=20, width=20)


p1 = VlnPlot(sobj, features=c("nFeature_RNA"), cols=col_vec, 
        group.by = 'cluster_label', sort = FALSE, pt.size = 0) +
  NoLegend()+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
        axis.title.x = element_blank())
p2 = VlnPlot(sobj, features=c("nCount_RNA"), cols=col_vec, 
        group.by = 'cluster_label', sort = FALSE, pt.size = 0) +
  NoLegend()+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
        axis.title.x = element_blank())

p4= VlnPlot(sobj, features=c("nFeature_ADT"), cols=col_vec, 
        group.by = 'cluster_label', sort = FALSE, pt.size = 0) +
  NoLegend()+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
        axis.title.x = element_blank())
p5 = VlnPlot(sobj, features=c("nCount_ADT"), cols=col_vec, 
        group.by = 'cluster_label', sort = FALSE, pt.size = 0) +
  NoLegend()+
  ylim(0,20000)+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
        axis.title.x = element_blank())


p3 =VlnPlot(sobj, features=c("percent.mt"), cols=col_vec, 
        group.by = 'cluster_label', sort = FALSE, pt.size = 0) +
  NoLegend()+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
        axis.title.x = element_blank())
p6 = VlnPlot(sobj, features=c("percent.ribo"), cols=col_vec, 
        group.by = 'cluster_label', sort = FALSE, pt.size = 0) +
  NoLegend()+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
        axis.title.x = element_blank())

cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=3)
ggsave(sprintf("%s/qc_vln.png", OUT_DIR), height=10, width=20)

# dotplots
load(sprintf("%s/wnn_rna_1.4.RData", OUT_DIR)) # wnn_markers
load(sprintf("%s/wnn_adt_1.4.RData", OUT_DIR)) # wnn_markers_adt



#### ADT PLOT
Idents(sobj) = "cluster_label"
d = cbind(sobj@meta.data, as.data.frame(t(as.matrix(sobj@assays$integrated.ADT@data))))

# calculate the median protein expression separately for each cluster 
adt_plot = d %>%
  dplyr::group_by(cluster_label) %>% 
  dplyr::summarize_at(.vars = rownames(sobj@assays[["integrated.ADT"]]), 
                      .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("cluster_label") 


# plot a heatmap of the average dsb normalized values for each cluster
adt_plot2 = adt_plot
adt_plot2[adt_plot2>100]=100
colmax = apply(adt_plot2, 2, max)
adt_plot3=adt_plot2[,colmax >10]
png(filename=sprintf("%s/wnn_adt_heatmap.png", OUT_DIR), width = 12, height = 23, units = "in", res = 300)
print(pheatmap::pheatmap(t(adt_plot2), 
                         color = viridis::viridis(25, option = "B"), 
                         fontsize_row = 8, border_color = NA))
dev.off()


png(filename=sprintf("%s/wnn_adt_heatmap2.png", OUT_DIR), width = 12, height = 15, units = "in", res = 300)
print(pheatmap::pheatmap(t(adt_plot3), 
                         color = viridis::viridis(25, option = "B"), 
                         fontsize_row = 8, border_color = NA))
dev.off()

DefaultAssay(sobj) = "integrated.ADT"

bcells_n_unknown = subset(sobj, 
                          cells=(sobj@meta.data %>% as_tibble(rownames="cell_id") %>% 
                                   filter(str_detect(cluster_label, "B_") | cluster_label=="unknown 2", cluster_label!="B_Naive 3") %>% pull(cell_id)))

adt_bcells = FindAllMarkers(bcells_n_unknown, test.use="MAST", only.pos=T)
save(adt_bcells, file=sprintf("%s/adt_bcells_markers.RData", OUT_DIR))


top5_bcell = adt_bcells %>% 
  group_by(cluster) %>%
  slice_max(n=5, order_by = avg_log2FC) %>%
  pull(gene)
source("/krummellab/data1/erflynn/scratch_scripts/FuncPctDotPlot.R")
FuncPctDotPlot(bcells_n_unknown, assay="integrated.ADT", features = unique(top5_bcell),
               cols="RdYlBu")+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave(sprintf("%s/bcell_dotplot.png", OUT_DIR))

tcells = subset(sobj, cells=(sobj@meta.data %>%
                               as_tibble(rownames="cell_id") %>% filter(str_detect(cluster_label, "^T") |
                                                                          cluster_label=="unknown 1") %>% pull(cell_id)))
adt_tcells = FindAllMarkers(tcells, test.use="MAST", only.pos=T)
save(adt_tcells, file=sprintf("%s/adt_tcells_markers.RData", OUT_DIR))

top5_tcell = adt_tcells %>% 
  group_by(cluster) %>%
  slice_max(n=5, order_by = avg_log2FC) %>%
  pull(gene)
FuncPctDotPlot(tcells, assay="integrated.ADT", features = unique(c(top5_tcell), c("CD4", "CD8")),
               cols="RdYlBu")+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave(sprintf("%s/tcell_dotplot.png", OUT_DIR), height=12, width=8)

### PLOTS



###



top5_adt_genes = wnn_markers_adt %>% 
  left_join(clus_map) %>%
  group_by(cluster_label) %>%
  slice_max(n=5, order_by = avg_log2FC) %>%
  arrange(cluster_label)
top5_adt_genes$gene = factor(top5_adt_genes$gene, levels=unique(top5_adt_genes$gene))
source("/krummellab/data1/erflynn/scratch_scripts/FuncPctDotPlot.R")
FuncPctDotPlot(sobj, assay="integrated.ADT", features = unique(top5_adt_genes$gene),
               cols="RdYlBu")+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave(sprintf("%s/adt_dotplot.pdf", OUT_DIR), width=14, height=17)

top5_genes = wnn_markers %>% 
  left_join(clus_map) %>%
  group_by(cluster_label) %>%
  slice_max(n=5, order_by = avg_log2FC) %>%
  pull(gene)
DotPlot(sobj, assay="RNA", features=unique(top5_genes),
        cols="RdYlBu")+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave(sprintf("%s/rna_dotplot.pdf", OUT_DIR), width=11, height=20)

### ALLUVIAL


###

### check mixing

DimPlot(sobj, reduction="wnn.umap", group.by="well")
ggsave(sprintf("%s/wnn_umap_well.pdf", OUT_DIR), height=5, width=5)


sobj = AddMetaData(sobj, sobj@meta.data %>% 
                     as_tibble(rownames="cell_id") %>%
                     dplyr::select(cell_id) %>%
                     left_join(  
                       meta %>% 
                         select(cell_id, donor)) %>% 
                     pull(donor), "donor")

DimPlot(sobj, reduction="wnn.umap", group.by="donor")
ggsave(sprintf("%s/wnn_umap_donor.pdf", OUT_DIR))


sobj = AddMetaData(sobj, sobj@meta.data %>% 
                     as_tibble(rownames="cell_id") %>%
                     dplyr::select(cell_id) %>%
                     left_join(  
                       meta %>% 
                         select(cell_id, study_id_visit)) %>% 
                     pull(study_id_visit), "study_id_visit")

DimPlot(sobj, reduction="wnn.umap", group.by="study_id_visit")
ggsave(sprintf("%s/wnn_umap_study_id_visit.pdf", OUT_DIR))

sobj@meta.data %>%
  dplyr::select(study_id_visit, cluster_label) %>%
  ggplot(aes(x=cluster_label, fill=study_id_visit))+
  geom_bar(position="fill")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))+
  ylab("fraction of cells")
ggsave(sprintf("%s/wnn_clus_study_id_visit.pdf", OUT_DIR))

