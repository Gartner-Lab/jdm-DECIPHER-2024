library(Seurat)
library(tidyverse)
library(dittoSeq)
## 1) plot louvain clusters without the tiny clusters
JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm"

WNN_DIR = sprintf("%s/wnn_out", JDM_DIR)
load(file=sprintf("%s/wnn_obj.RData", WNN_DIR))

### to re-load wsnn1
if (reload == T){
  sobj =  subset(sobj, wsnn_res.1 %in% (sobj@meta.data %>% 
                                          group_by(wsnn_res.1 ) %>% 
                                          dplyr::count() %>% filter(n>=10) %>%
                                          pull(wsnn_res.1)))
  
  load(file=sprintf("%s/wsnn_res.1_clus_meta.RData", WNN_DIR))
  sobj@meta.data = meta.data
  Idents(sobj) = "cluster_label"
  DimPlot(sobj, reduction = 'wnn.umap',
          label = TRUE, cols=dittoColors())+NoLegend()
}




#  --- 0.4 --- #
sobj_no_tiny =  subset(sobj, wsnn_res.0.4 %in% (sobj@meta.data %>% 
                                                  group_by(wsnn_res.0.4 ) %>% 
                                                  dplyr::count() %>% filter(n>=10) %>%
                                                  pull(wsnn_res.0.4)))

DimPlot(sobj_no_tiny, group.by="wsnn_res.0.4", label=T, reduction="wnn.umap",
        cols=dittoColors()) + NoLegend()
ggsave(sprintf("%s/wsnn_0.4_no_tiny.png", WNN_DIR), width=5, height=5, units="in")



#  --- 0.6 --- #
sobj_no_tiny =  subset(sobj, wsnn_res.0.6 %in% (sobj@meta.data %>% 
                                                  group_by(wsnn_res.0.6 ) %>% 
                                                  dplyr::count() %>% filter(n>=10) %>%
                                                  pull(wsnn_res.0.6)))

DimPlot(sobj_no_tiny, group.by="wsnn_res.0.6", label=T, reduction="wnn.umap",
        cols=dittoColors()) + NoLegend()
ggsave(sprintf("%s/wsnn_0.6_no_tiny.png", WNN_DIR), width=5, height=5, units="in")





#  --- 0.6 --- #
sobj_no_tiny =  subset(sobj, wsnn_res.0.6 %in% (sobj@meta.data %>% 
                                                  group_by(wsnn_res.0.6 ) %>% 
                                                  dplyr::count() %>% filter(n>=10) %>%
                                                  pull(wsnn_res.0.6)))

DimPlot(sobj_no_tiny, group.by="wsnn_res.0.6", label=T, reduction="wnn.umap",
        cols=dittoColors()) + NoLegend()
ggsave(sprintf("%s/wsnn_0.6_no_tiny.png", WNN_DIR), width=5, height=5, units="in")



#  --- 0.8 --- #
sobj_no_tiny =  subset(sobj, wsnn_res.0.8 %in% (sobj@meta.data %>% 
                                 group_by(wsnn_res.0.8 ) %>% 
                                 dplyr::count() %>% filter(n>=10) %>% pull(wsnn_res.0.8)))

DimPlot(sobj_no_tiny, group.by="wsnn_res.0.8", label=T, reduction="wnn.umap",
        cols=dittoColors()) + NoLegend()
ggsave(sprintf("%s/wsnn_0.8_no_tiny.png", WNN_DIR), width=5, height=5, units="in")


#  --- 1.0 --- #
sobj_no_tiny =  subset(sobj, wsnn_res.1 %in% (sobj@meta.data %>% 
                                                  group_by(wsnn_res.1 ) %>% 
                                                  dplyr::count() %>% filter(n>=10) %>%
                                                  pull(wsnn_res.1)))

DimPlot(sobj_no_tiny, group.by="wsnn_res.1", label=T, reduction="wnn.umap",
        cols=dittoColors()) + NoLegend()
ggsave(sprintf("%s/wsnn_1_no_tiny.png", WNN_DIR), width=5, height=5, units="in")


#  --- 1.2 --- #
load(sprintf("%s/wnn_obj_extra_graph_clus.RData", WNN_DIR))
sobj@meta.data = meta_data
sobj_no_tiny =  subset(sobj, wsnn_graph.1.2 %in% (sobj@meta.data %>% 
                                                group_by(wsnn_graph.1.2 ) %>% 
                                                dplyr::count() %>% filter(n>=10) %>%
                                                pull(wsnn_graph.1.2)))

DimPlot(sobj_no_tiny, group.by="wsnn_graph.1.2", label=T, reduction="wnn.umap",
        cols=dittoColors()) + NoLegend()
ggsave(sprintf("%s/wsnn_1.2_no_tiny.png", WNN_DIR), width=5, height=5, units="in")


#  --- 1.4 --- #
sobj_no_tiny =  subset(sobj, wsnn_graph.1.4 %in% (sobj@meta.data %>% 
                                                    group_by(wsnn_graph.1.4 ) %>% 
                                                    dplyr::count() %>% filter(n>=10) %>%
                                                    pull(wsnn_graph.1.4)))

DimPlot(sobj_no_tiny, group.by="wsnn_graph.1.4", label=T, reduction="wnn.umap",
        cols=dittoColors()) + NoLegend()
ggsave(sprintf("%s/wsnn_1.4_no_tiny.png", WNN_DIR), width=5, height=5, units="in")


#  --- 1.6 --- #
sobj_no_tiny =  subset(sobj, wsnn_graph.1.6 %in% (sobj@meta.data %>% 
                                                    group_by(wsnn_graph.1.6 ) %>% 
                                                    dplyr::count() %>% filter(n>=10) %>%
                                                    pull(wsnn_graph.1.6)))

DimPlot(sobj_no_tiny, group.by="wsnn_graph.1.6", label=T, reduction="wnn.umap",
        cols=dittoColors()) + NoLegend()
ggsave(sprintf("%s/wsnn_1.6_no_tiny.png", WNN_DIR), width=5, height=5, units="in")


#  --- 2.0 --- #
sobj_no_tiny =  subset(sobj, wsnn_graph.2 %in% (sobj@meta.data %>% 
                                                    group_by(wsnn_graph.2 ) %>% 
                                                    dplyr::count() %>% filter(n>=10) %>%
                                                    pull(wsnn_graph.2)))

DimPlot(sobj_no_tiny, group.by="wsnn_graph.2", label=T, reduction="wnn.umap",
        cols=dittoColors()) + NoLegend()
ggsave(sprintf("%s/wsnn_2_no_tiny.png", WNN_DIR), width=5, height=5, units="in")





#########

## keep going with 1.0 ##
sobj =  subset(sobj, wsnn_res.1 %in% (sobj@meta.data %>% 
                                                group_by(wsnn_res.1 ) %>% 
                                                dplyr::count() %>% filter(n>=10) %>%
                                                pull(wsnn_res.1)))

#
sobj@meta.data  = sobj@meta.data %>% mutate(wsnn_res.1=as.numeric(as.character(wsnn_res.1)))
Idents(sobj) = "wsnn_res.1"
DimPlot(sobj, group.by="wsnn_res.1", label=T, reduction="wnn.umap", cols=dittoColors())+
  NoLegend()
ggsave(sprintf("%s/wsnn_res1.0/wnn_umap.pdf", WNN_DIR))

VlnPlot(sobj, features = "RNA.weight", 
        group.by = 'wsnn_res.1', sort = TRUE, pt.size = 0) +
  NoLegend()
ggsave(sprintf("%s/wsnn_res1.0/wnn_weight.pdf", WNN_DIR))

VlnPlot(sobj, features = "RNA.weight", 
        group.by = 'wsnn_res.1', sort = FALSE, pt.size = 0) +
  NoLegend()
ggsave(sprintf("%s/wsnn_res1.0/wnn_weight_unsrt.pdf", WNN_DIR))


#######
## label the clusters
jdm_obj_initial = readRDS(sprintf("%s/analyzed_seurat_object (1).rds", JDM_DIR))
meta = jdm_obj_initial@meta.data %>% 
  as_tibble(rownames="cell_id")
rm(jdm_obj_initial)
length(unique(sobj@meta.data$wsnn1)) # 35
clus_lab = sobj@meta.data %>% 
  as_tibble(rownames="cell_id") %>%
  dplyr::select(cell_id, wsnn_res.1) %>%
  left_join(  
    meta %>% 
      dplyr::select(cell_id, manual_labels))
length(unique(clus_lab$wsnn_res.1))
clus_to_lab = clus_lab %>%
  dplyr::rename(cluster="wsnn_res.1") %>%
  group_by(cluster) %>%
  mutate(tot=n()) %>%
  group_by(cluster, manual_labels, tot) %>%
  dplyr::count() %>%
  arrange(cluster, desc(n)) %>%
  mutate(frac=n/tot)
length(unique(clus_to_lab$cluster)) # 35

clus_to_lab %>%
  ungroup() %>%
  group_by(cluster) %>%
  filter(frac > 0.4) %>%
  View()
clus_to_lab2 = clus_to_lab %>%
  ungroup() %>%
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


new_data = sobj@meta.data %>% 
  as_tibble(rownames="cell_id") %>%
  dplyr::select(cell_id, wsnn_res.1) %>%
  left_join(clus_to_lab3 %>% dplyr::select(cluster, cluster_label) %>%
              mutate(cluster_label=as.factor(cluster_label)), by=c("wsnn_res.1"="cluster")) 

manual_lab_df = clus_lab %>% 
  mutate(manual_labels=as.character(manual_labels),
         manual_labels=ifelse(is.na(manual_labels), "unknown",
                              manual_labels)) 


sobj = AddMetaData(sobj, new_data$cluster_label, col.name="cluster_label")
sobj = AddMetaData(sobj, manual_lab_df$manual_labels, col.name="manual_labels")
meta.data = sobj@meta.data
save(meta.data, file=sprintf("%s/wsnn_res.1_clus_meta.RData", WNN_DIR))

Idents(sobj) = "cluster_label"
## remake the plots
DimPlot(sobj, reduction = 'wnn.umap', 
        label = TRUE, cols=dittoColors(),
        label.size=3) + NoLegend()
ggsave(sprintf("%s/wsnn_res1.0/wnn_umap_clus_lab.pdf", WNN_DIR))


VlnPlot(sobj, features = "RNA.weight", 
        group.by = 'cluster_label', cols=dittoColors(),
        sort = TRUE, pt.size = 0) +
  NoLegend()
ggsave(sprintf("%s/wsnn_res1.0/wnn_umap_clus_lab_weight.pdf", WNN_DIR))

## ones with high ADT
# look at unknown 2, 3, 6, 7, 8
# look at cM 6, T4 Naive 3


# compare to UMAP of RNA, ADT
sobj <- RunUMAP(sobj, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
sobj <- RunUMAP(sobj, reduction = 'a.pca', dims = 1:18, assay = 'ADT', 
                reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

Idents(sobj) = "manual_labels"

p1 = DimPlot(sobj, reduction = 'wnn.umap', 
             label = TRUE, cols=dittoColors()[1:30],
             repel=T,
             label.size=3) + NoLegend()
p2 = DimPlot(sobj, reduction = 'rna.umap', cols=dittoColors()[1:30], label = TRUE, 
             repel = TRUE, label.size = 3) + NoLegend()
p3 = DimPlot(sobj, reduction = 'adt.umap', cols=dittoColors()[1:30], label = TRUE, 
             repel = TRUE, label.size = 3) + NoLegend()
cowplot::plot_grid(p1,p2,p3, ncol=3)
ggsave(sprintf("%s/wsnn_res1.0/wnn_umap_compare.pdf", WNN_DIR))

Idents(sobj) = "cluster_label"

# make results dataframe 
d = cbind(sobj@meta.data, as.data.frame(t(sobj@assays$ADT@data)))

# calculate the median protein expression separately for each cluster 
adt_plot = d %>%
  dplyr::group_by(cluster_label) %>% 
  dplyr::summarize_at(.vars = rownames(sobj@assays[["ADT"]]), 
                      .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("cluster_label") 


# plot a heatmap of the average dsb normalized values for each cluster
adt_plot2 = adt_plot
adt_plot2[adt_plot2>100]=100
colmax = apply(adt_plot2, 2, max)
adt_plot3=adt_plot2[,colmax >30]
png(filename=sprintf("%s/wsnn_res1.0/wnn_adt_heatmap.png", WNN_DIR), width = 12, height = 23, units = "in", res = 300)
print(pheatmap::pheatmap(t(adt_plot2), 
                         color = viridis::viridis(25, option = "B"), 
                         fontsize_row = 8, border_color = NA))
dev.off()


png(filename=sprintf("%s/wsnn_res1.0/wnn_adt_heatmap_expr.png", WNN_DIR), width = 12, height = 15, units = "in", res = 300)
print(pheatmap::pheatmap(t(adt_plot3), 
                         color = viridis::viridis(25, option = "B"), 
                         fontsize_row = 8, border_color = NA))
dev.off()

load(sprintf("%s/wnn_adt_markers_res_1.RData", WNN_DIR))
sobj@meta.data %>% distinct(wsnn_res.1, cluster_label)
top5_adt_genes = wnn_adt_markers %>% 
  mutate(cluster=as.numeric(as.character(cluster))) %>%
  left_join(sobj@meta.data %>% 
              distinct(wsnn_res.1, cluster_label) %>%
              mutate(wsnn_res.1 = as.numeric(as.character(wsnn_res.1))), 
            by=c("cluster"="wsnn_res.1")) %>%
  group_by(cluster_label) %>%
  slice_max(n=5, order_by = avg_log2FC)  %>%
  ungroup() %>%
  distinct(gene)

source("/krummellab/data1/erflynn/scratch_scripts/FuncPctDotPlot.R")
FuncPctDotPlot(sobj, assay="ADT", features = top5_adt_genes$gene,
               cols="RdYlBu")+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave(sprintf("%s/wsnn_res1.0/adt_dotplot.pdf", WNN_DIR), width=10, height=16)


# look at the RNA markers
load(sprintf("%s/wnn_markers_res_1.RData", WNN_DIR)) # wnn_markers
top5_gene_l = wnn_markers %>% 
  mutate(cluster=as.numeric(as.character(cluster))) %>%
  left_join(sobj@meta.data %>% 
              distinct(wsnn_res.1, cluster_label) %>%
              mutate(wsnn_res.1 = as.numeric(as.character(wsnn_res.1))), 
            by=c("cluster"="wsnn_res.1")) %>%
  group_by(cluster_label) %>%
  slice_max(n=5, order_by = avg_log2FC) %>%
  ungroup() %>%
  distinct(gene)

top5_gene_l$gene = factor(top5_gene_l$gene, levels=top5_gene_l$gene )

DotPlot(sobj, 
        assay = "RNA", 
        features = top5_gene_l$gene,
        cols="RdYlBu")+
  coord_flip()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.1))
ggsave(sprintf("%s/wsnn_res1.0/rna_dotplot.pdf", 
               WNN_DIR),
       height=19, width=11.5)

## QC plot
p1 = VlnPlot(sobj, "percent.mt", cols=dittoColors(),sort=F, pt.size = 0) + NoLegend()
p2 = VlnPlot(sobj, "percent.ribo", cols=dittoColors(),sort=F, pt.size = 0) + NoLegend()
p3 = VlnPlot(sobj, "nCount_RNA", cols=dittoColors(),sort=F, pt.size = 0) + NoLegend()
p4 = VlnPlot(sobj, "nFeature_RNA", cols=dittoColors(),sort=F, pt.size = 0) + NoLegend()
p5 = VlnPlot(sobj, "nFeature_RNA", cols=dittoColors(),sort=F, pt.size = 0) + NoLegend()+ylim(0, 1000)
p6 = VlnPlot(sobj, "nCount_ADT", cols=dittoColors(),sort=F, pt.size = 0) + NoLegend()
p7 = VlnPlot(sobj, "nCount_ADT", cols=dittoColors(),sort=F, pt.size = 0) + NoLegend()+ylim(0, 1000)
p8 = VlnPlot(sobj, "nFeature_ADT", cols=dittoColors(),sort=F, pt.size = 0) + NoLegend()
cowplot::plot_grid( p3, p4, p5, p8, p6, p7, p1, p2, ncol=4)
ggsave(sprintf("%s/wsnn_res1.0/qc_vln.pdf", WNN_DIR), width=24, height=9)

# It would be helpful to see a UMAP (WNN) of the clusters themselves side-by-side 
# with the manual labels, and then we can assign them together? 

p1 = DimPlot(sobj, reduction="wnn.umap", group.by="wsnn_res.1", label=T, cols=dittoColors())+NoLegend()
p2 = DimPlot(sobj, reduction="wnn.umap", group.by="manual_labels", 
             label=T, cols=c(dittoColors()[1:26],"light gray"),
             label.size=2.5)+NoLegend()
plot_grid(p1, p2, ncol=2)
ggsave(sprintf("%s/wsnn_res1.0/wnn_umap_clusters_labels.pdf", WNN_DIR), height=5, width=10)

p1.1 = DimPlot(sobj, reduction="wnn.umap", group.by="cluster_label", 
               label=T, cols=dittoColors(), label.size=2.5)+NoLegend()
plot_grid(p1.1, p2, ncol=2)
ggsave(sprintf("%s/wsnn_res1.0/wnn_umap_clusters_labeled2_labels.pdf", WNN_DIR), height=5, width=10)


# It would also be helpful to see umap colored by batch and study_id_visit 
# variables to see if we still have those clusters separating out by patient 
# with this method.

DimPlot(sobj, reduction="wnn.umap", group.by="well")
ggsave(sprintf("%s/wsnn_res1.0/wnn_umap_well.pdf", WNN_DIR), height=5, width=5)

#meta %>% 
head(sobj@meta.data)
sobj = AddMetaData(sobj, sobj@meta.data %>% 
                     as_tibble(rownames="cell_id") %>%
                     dplyr::select(cell_id) %>%
                     left_join(  
                       meta %>% 
                         dplyr::select(cell_id, donor)) %>% 
                     pull(donor), "donor")

DimPlot(sobj, reduction="wnn.umap", group.by="donor")
ggsave(sprintf("%s/wsnn_res1.0/wnn_umap_donor.pdf", WNN_DIR), height=5, width=5)


sobj = AddMetaData(sobj, sobj@meta.data %>% 
                     as_tibble(rownames="cell_id") %>%
                     dplyr::select(cell_id) %>%
                     left_join(  
                       meta %>% 
                         dplyr::select(cell_id, study_id_visit)) %>% 
                     pull(study_id_visit), "study_id_visit")

DimPlot(sobj, reduction="wnn.umap", group.by="study_id_visit")
ggsave(sprintf("%s/wsnn_res1.0/wnn_umap_study_id_visit.pdf", WNN_DIR), height=5, width=5)

sobj@meta.data %>%
  dplyr::select(study_id_visit, cluster_label) %>%
  ggplot(aes(x=cluster_label, fill=study_id_visit))+
  geom_bar(position="fill")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))+
  ylab("fraction of cells")
ggsave(sprintf("%s/wsnn_res1.0/wnn_clus_study_id_visit.pdf", WNN_DIR), height=5, width=10)


### TODO: QC plots


##### LOOK AT THE PATH OF SOME OF THE high ADT ones ####

# esp compare to their neighbors

# unknown 7 + unknown 4 are near the B cells



## ones with high ADT
# look at unknown 2, 3, 6, 7, 8
# look at cM 6, T4 Naive 3

# do the cM cells actually look different?

cM6_adt_markers = FindMarkers(sobj, ident.1="cM 6", ident.2=c("cM 1", "cM 2"), test.use="MAST", only.pos=T, assay="ADT")
cM6_rna_markers = FindMarkers(sobj, ident.1="cM 6", ident.2=c("cM 1", "cM 2"),  test.use="MAST", only.pos=T, assay="RNA")

# look at unknown 2,4,5,6



# what was the location of these cells in the RNA?

p0 = DimPlot(sobj, reduction = 'umap', group.by="manual_labels", label=T, label.size = 3) + 
  NoLegend()+ylab("rna.umap2")+xlab("rna.umap1")

p9= DimPlot(sobj, reduction = 'umap', 
             cells.highlight=rownames(sobj@meta.data %>% filter(cluster_label=="cM 6")), 
             label.size = 3) + NoLegend()+
  ylab("rna.umap2")+xlab("rna.umap1")+
  ggtitle("cM 6")
 
p10= DimPlot(sobj, reduction = 'umap', 
         cells.highlight=rownames(sobj@meta.data %>% filter(cluster_label=="T4_Naive 3")), 
         label.size = 3) + NoLegend()+
   ylab("rna.umap2")+xlab("rna.umap1")+
   ggtitle("T4_Naive 3")


p1 = DimPlot(sobj, reduction = 'umap', 
        cells.highlight=rownames(sobj@meta.data %>% filter(cluster_label=="unknown 1")), 
        label.size = 3) + NoLegend()+
  ylab("rna.umap2")+xlab("rna.umap1")+
  ggtitle("unknown 1")

p2 = DimPlot(sobj, reduction = 'umap', 
        cells.highlight=rownames(sobj@meta.data %>% filter(cluster_label=="unknown 2")), 
        label.size = 3) + NoLegend()+
  ggtitle("unknown 2")+ylab("rna.umap2")+xlab("rna.umap1")

p3 = DimPlot(sobj, reduction = 'umap', 
        cells.highlight=rownames(sobj@meta.data %>% filter(cluster_label=="unknown 3")), 
        label.size = 3) + NoLegend()+
  ggtitle("unknown 3")+ylab("rna.umap2")+xlab("rna.umap1")

p4 = DimPlot(sobj, reduction = 'umap', 
        cells.highlight=rownames(sobj@meta.data %>% filter(cluster_label=="unknown 4")), 
        label.size = 3) + NoLegend()+
  ggtitle("unknown 4")+ylab("rna.umap2")+xlab("rna.umap1")

p5 = DimPlot(sobj, reduction = 'umap', 
        cells.highlight=rownames(sobj@meta.data %>% filter(cluster_label=="unknown 5")), 
        label.size = 3) + NoLegend()+
  ggtitle("unknown 5")+ylab("rna.umap2")+xlab("rna.umap1")

p6 = DimPlot(sobj, reduction = 'umap', 
        cells.highlight=rownames(sobj@meta.data %>% filter(cluster_label=="unknown 6")), 
        label.size = 3) + NoLegend()+
  ggtitle("unknown 6")+ylab("rna.umap2")+xlab("rna.umap1")

p7 = DimPlot(sobj, reduction = 'umap', 
        cells.highlight=rownames(sobj@meta.data %>% filter(cluster_label=="unknown 7")), 
        label.size = 3) + NoLegend()+
  ggtitle("unknown 7")+ylab("rna.umap2")+xlab("rna.umap1")

p8 = DimPlot(sobj, reduction = 'umap', 
        cells.highlight=rownames(sobj@meta.data %>% filter(cluster_label=="unknown 8")), 
        label.size = 3) + NoLegend()+
  ggtitle("unknown 8")+ylab("rna.umap2")+xlab("rna.umap1")


cowplot::plot_grid( p0, p9, p10, p1, p2, p3, p4, p5, p6, p7, p8, ncol=4)
ggsave(sprintf("%s/wsnn_res1.0/high_adt_in_rna.pdf", WNN_DIR), width=16, height=12)

## 
Idents(sobj) = "wsnn_res.1"
high_adt_cts = sobj@meta.data %>%
  as_tibble(rownames="cell_id") %>%
  filter(nCount_ADT > 30000) %>%
  pull(cell_id)
# only 43

VlnPlot(sobj, "nFeature_ADT", pt.size=0)
VlnPlot(sobj, "nCount_ADT", pt.size=0)
DimPlot(sobj, cells.highlight=high_adt_cts, reduction="wnn.umap")
FeaturePlot(subset(sobj, nCount_ADT < 30000), "nCount_ADT",
            reduction="wnn.umap")






