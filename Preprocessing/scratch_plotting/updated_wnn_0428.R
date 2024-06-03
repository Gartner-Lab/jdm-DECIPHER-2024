library(Seurat)
library(tidyverse)
library(dittoSeq)
DATA_DIR = "/krummellab/data1/erflynn/premier/jdm/data/"
setwd(sprintf("%s/jdm/wnn_out_rpca_DSB3/", DATA_DIR))

source("/krummellab/data1/erflynn/scratch_scripts/FuncPctDotPlot.R")


load(file=sprintf("%s/jdm/wnn_out_rpca_DSB3/sobj_filt.RData", DATA_DIR)) # --> sobj_filt
jdm_data = sobj_filt
DefaultAssay(jdm_data)= "RNA"
rm(sobj_filt)

## add the JDM metadata
jdm_obj_meta = readRDS( sprintf("%s/jdm/seurat_filt_postdemuxDF_w_meta.rds", DATA_DIR))
jdm_meta = jdm_obj_meta@meta.data
rm(jdm_obj_meta)


new_meta = jdm_data@meta.data %>%
  as_tibble(rownames="cell_id") %>%
  select(-study_id_visit, -donor) %>% # they match and the other meta is more complete
  left_join(jdm_meta %>% 
              as_tibble(rownames="cell_id") %>% 
              select(cell_id, study_id_visit:age),
            by="cell_id")
all(new_meta$cluster_label==jdm_data@meta.data$cluster_label)
jdm_data@meta.data = new_meta %>% 
  column_to_rownames("cell_id")

# plot the RNA at a bunch of resolutions
DimPlot(jdm_data, group.by="RNA_snn_res.0.8", 
        label=T, cols=dittoColors(),
        reduction="umap")+
  NoLegend()
ggsave("dimplot_RNA_res.0.8.png", height=5, width=5)
DimPlot(jdm_data, group.by="RNA_snn_res.1", 
        label=T, cols=dittoColors(),
        reduction="umap")+
  NoLegend()
ggsave("dimplot_RNA_res.1.png", height=5, width=5)

DimPlot(jdm_data, group.by="wsnn_res.1.4", 
        label=T, cols=dittoColors(),
        reduction="wnn.umap")+
  NoLegend()
ggsave("dimplot_wnn_res.1.4.png", height=5, width=5)



# pick one and run find all markers


# to compare RNA and WNN:
# NEED:
# [x] RNA clustered at a similar resolution
# [ ] RNA dotplots
# [ ] redo WNN DSB dotplots w/ wilcoxon


### TODO:
# re-run everything removing uninformative genes
# min.cell = 3
rs = rowSums(jdm_data@assays[["RNA"]]@data)
summary(rs)
#cs = colSums(jdm_data@assays[["RNA"]]@data)
#summary(cs) # min is 261, so we're fine


### plot WNN vs RNA

# [ ] add labels + redo dimplots
# "RNA_snn_res.1"
# cluster_label



clus_to_lab = jdm_data@meta.data  %>% 
  as_tibble(rownames="cell_id") %>%
  dplyr::select(cell_id, RNA_snn_res.1, manual_labels) %>%
  dplyr::rename(cluster="RNA_snn_res.1") %>%
  group_by(cluster) %>%
  mutate(tot=n()) %>%
  group_by(cluster, manual_labels, tot) %>%
  dplyr::count() %>%
  arrange(cluster, desc(n)) %>%
  mutate(frac=n/tot)

clus_to_lab2 = clus_to_lab %>%
  ungroup() %>%
  filter(tot > 10) %>%
  group_by(cluster) %>%
  slice_max(1)

clus_to_lab3 = clus_to_lab2 %>%
  mutate(manual_labels=as.character(manual_labels)) %>%
  mutate(assigned_label=ifelse(frac < 0.4 | is.na(manual_labels), "unknown", manual_labels)) %>%
  arrange(assigned_label) %>%
  ungroup() %>%
  group_by(assigned_label) %>%
  mutate(idx=1:n(),
         nlab=n()) %>%
  mutate(rna_cluster_label=ifelse(nlab ==1, assigned_label,
                              paste(assigned_label, idx))) %>%
  ungroup()
clus_map = clus_to_lab3 %>% distinct(cluster, rna_cluster_label)

new_data = jdm_data@meta.data %>% 
  as_tibble(rownames="cell_id") %>%
  select(cell_id, RNA_snn_res.1) %>%
  left_join(clus_to_lab3 %>% select(cluster, rna_cluster_label) %>%
              mutate(rna_cluster_label=as.factor(rna_cluster_label)), by=c("RNA_snn_res.1"="cluster"))

jdm_data = AddMetaData(jdm_data, new_data$rna_cluster_label, col.name="rna_cluster_label")


wnn_clus_to_lab = jdm_data@meta.data  %>% 
  as_tibble(rownames="cell_id") %>%
  dplyr::select(cell_id, wsnn_res.1.4, manual_labels) %>%
  dplyr::rename(cluster="wsnn_res.1.4") %>%
  group_by(cluster) %>%
  mutate(tot=n()) %>%
  group_by(cluster, manual_labels, tot) %>%
  dplyr::count() %>%
  arrange(cluster, desc(n)) %>%
  mutate(frac=n/tot)

wnn_clus_to_lab2 = wnn_clus_to_lab %>%
  ungroup() %>%
  filter(tot > 10) %>%
  group_by(cluster) %>%
  slice_max(1)

wnn_clus_to_lab3 = wnn_clus_to_lab2 %>%
  mutate(manual_labels=as.character(manual_labels)) %>%
  mutate(assigned_label=ifelse(frac < 0.4 | is.na(manual_labels), "unknown", manual_labels)) %>%
  arrange(assigned_label) %>%
  ungroup() %>%
  group_by(assigned_label) %>%
  mutate(idx=1:n(),
         nlab=n()) %>%
  mutate(wnn_cluster_label=ifelse(nlab ==1, assigned_label,
                                  paste(assigned_label, idx))) %>%
  ungroup()
wnn_clus_map = wnn_clus_to_lab3 %>% distinct(cluster, wnn_cluster_label)

wnn_new_data = jdm_data@meta.data %>% 
  as_tibble(rownames="cell_id") %>%
  select(cell_id, wsnn_res.1.4) %>%
  left_join(wnn_clus_to_lab3 %>% select(cluster, wnn_cluster_label) %>%
              mutate(wnn_cluster_label=as.factor(wnn_cluster_label)), 
            by=c("wsnn_res.1.4"="cluster"))
jdm_data = AddMetaData(jdm_data, wnn_new_data$wnn_cluster_label, col.name="wnn_cluster_label")



Idents(jdm_data) = "rna_cluster_label"
DimPlot(jdm_data, reduction = 'umap', 
        label = TRUE, cols=dittoColors(),
        label.size=3) + NoLegend()+
  ylab("rnaUMAP_2")+
  xlab("rnaUMAP_1")
ggsave("rna_umap_clus_lab.pdf", height=5, width=5)


#unknown_jdm = subset(jdm_data, manual_labels=="unknown")
#save(unknown_jdm, file="unknown_jdm.RData")
load("azimuth_unknown_jdm.RData")
unknown_rna_map = unknown_meta %>% 
  group_by(rna_cluster_label) %>%
  mutate(tot=n()) %>%
  group_by(rna_cluster_label, predicted.celltype.l2, tot) %>%
  dplyr::count() %>%
  arrange(rna_cluster_label, desc(n)) %>%
  mutate(frac=n/tot) %>%
  ungroup() %>%
  group_by(rna_cluster_label) %>%
  slice_max(1) %>%
  filter(str_detect(rna_cluster_label, "unknown")) %>%
  filter(tot > 10, frac > 0.4)

unknown_wnn_map = unknown_meta %>% 
  left_join(jdm_data@meta.data %>% distinct(wsnn_res.1.4, wnn_cluster_label)) %>%
  group_by(wnn_cluster_label) %>%
  mutate(tot=n()) %>%
  group_by(wnn_cluster_label, predicted.celltype.l2, tot) %>%
  dplyr::count() %>%
  arrange(wnn_cluster_label, desc(n)) %>%
  mutate(frac=n/tot) %>%
  ungroup() %>%
  group_by(wnn_cluster_label) %>%
  slice_max(1) %>%
  filter(str_detect(wnn_cluster_label, "unknown")) %>%
  filter(tot > 10, frac > 0.4)

unknown_rna_map2 = tibble(
  rna_cluster_label = unknown_rna_map$rna_cluster_label,
  rna_cluster_label2=c("B_Mem_Intermediate 2*", "T8_TEMRA 2*", "Eryth*", 
              "B_Mem_Intermediate 3*", "T4_Mem 3*")
)

unknown_wnn_map2 = tibble(
  wnn_cluster_label = unknown_wnn_map$wnn_cluster_label,
  wnn_cluster_label2=c( "B_Mem_Intermediate 2*", "T8_TEMRA 2*", "T4_Mem 2*")
)
jdm_data@meta.data$rna_cluster_label2=NULL
jdm_data = AddMetaData(jdm_data, jdm_data@meta.data %>%
  left_join(unknown_rna_map2, by="rna_cluster_label") %>%
    mutate(rna_cluster_label2=ifelse(is.na(rna_cluster_label2), 
                                     as.character(rna_cluster_label), rna_cluster_label2)) %>%
    pull(rna_cluster_label2), "rna_cluster_label2")
jdm_data@meta.data %>% distinct(RNA_snn_res.1, rna_cluster_label2)

jdm_data@meta.data$wnn_cluster_label2=NULL

jdm_data = AddMetaData(jdm_data, jdm_data@meta.data %>%
                         left_join(unknown_wnn_map2, by="wnn_cluster_label") %>%
                         mutate(wnn_cluster_label2=ifelse(is.na(wnn_cluster_label2), 
                                                          as.character(wnn_cluster_label), 
                                                          wnn_cluster_label2)) %>%
                         pull(wnn_cluster_label2), "wnn_cluster_label2")
jdm_data@meta.data %>% distinct(wsnn_res.1.4, wnn_cluster_label2)

labs = sort(union(jdm_data@meta.data$rna_cluster_label2, 
                  jdm_data@meta.data$wnn_cluster_label2))
col_vec2 = dittoColors()[1:length(labs)]
names(col_vec2) = labs
Idents(jdm_data) = "wnn_cluster_label2"
VlnPlot(jdm_data, sort=T, features="RNA.weight", cols=col_vec2, pt.size=0)+NoLegend()
ggsave("filt_rna_weight2.pdf")

DimPlot(jdm_data, reduction = 'wnn.umap', 
        label = TRUE, cols=col_vec2,
        label.size=3) + NoLegend()
ggsave("wnn_umap_clus_lab2.pdf", height=5, width=5)

Idents(jdm_data) = "rna_cluster_label2"
DimPlot(jdm_data, reduction = 'umap', 
        label = TRUE, cols=col_vec2,
        label.size=3) + NoLegend()+
  ylab("rnaUMAP_2")+
  xlab("rnaUMAP_1")
ggsave("rna_umap_clus_lab2.pdf", height=5, width=5)

# [ ] ALLUVIAL
library(ggalluvial)

meta_alluv =jdm_data@meta.data  %>%
  as_tibble(rownames="cell_id")  %>%
  select(cell_id, rna_cluster_label2, wnn_cluster_label2) %>%
  dplyr::rename(RNA=rna_cluster_label2,
                WNN=wnn_cluster_label2) %>%
  pivot_longer(c("RNA",  "WNN"),
               names_to="annot", values_to="cell_type") %>%
  mutate(freq=1) %>%
  mutate(annot=factor(annot, levels=c("RNA", "WNN")))


ggplot(meta_alluv,
       aes(x = annot, stratum = cell_type, alluvium = cell_id,
           y = freq,
           fill = cell_type, label = cell_type)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "none") +
  theme_bw()+
  theme(panel.grid=element_blank())+
  ylab("number of cells")+
  xlab("annotation")
ggsave("alluv_annot2.pdf", height=15, width=10)

# [ ] dotplot top RNA genes per cluster for WNN, RNA clusters
load("rna_only_markers.RData") # oops --> wnn_markers
rna_markers = wnn_markers
rm(wnn_markers)
top5_genes_rna_only = rna_markers %>% 
  group_by(cluster) %>%
  slice_max(n=5, order_by = avg_log2FC) %>%
  pull(gene)
Idents(jdm_data) = "rna_cluster_label2"
DotPlot(jdm_data, assay="RNA", features=unique(top5_genes_rna_only),
        cols="RdYlBu", cluster.idents = T)+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("rna_only_clus_dotplot2.pdf", width=11, height=20)

load("wnn_rna_filt_1.4.RData")
top5_genes = wnn_markers %>% 
  group_by(cluster) %>%
  slice_max(n=5, order_by = avg_log2FC) %>%
  pull(gene)
Idents(jdm_data) = "wnn_cluster_label2"
DotPlot(jdm_data, assay="RNA", features=unique(top5_genes),
        cols="RdYlBu", cluster.idents = T)+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("wnn_rna_dotplot2.pdf", width=11, height=20)


load("wnn_adt_wilcoxon.RData")
top5_adt_genes = wnn_markers_adt %>% 
  left_join(jdm_data@meta.data %>%
              distinct(wsnn_res.1.4, wnn_cluster_label2),
            by=c("cluster"="wsnn_res.1.4")) %>%
  slice_max(n=5, order_by = avg_log2FC) %>%
  arrange(wnn_cluster_label2) %>%
  pull(gene)
Idents(jdm_data) = "wnn_cluster_label2"

FuncPctDotPlot(jdm_data, assay="integrated.ADT", features = unique(top5_adt_genes),
               cols="RdYlBu",
               cluster.idents = T)+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("adt_dotplot2.pdf", width=12, height=13)


## try dotplots together??


# [ ] heatmap

sobj=jdm_data
Idents(sobj) = wnn_cluster_label2

d = cbind(sobj@meta.data, as.data.frame(t(as.matrix(sobj@assays$integrated.ADT@data))))

# calculate the median protein expression separately for each cluster 
adt_plot = d %>%
  dplyr::group_by(wnn_cluster_label2) %>% 
  dplyr::summarize_at(.vars = rownames(sobj@assays[["integrated.ADT"]]), 
                      .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("wnn_cluster_label2") 


# plot a heatmap of the average dsb normalized values for each cluster
adt_plot2 = adt_plot
adt_plot2[adt_plot2>100]=100
colmax = apply(adt_plot2, 2, max)
adt_plot3=adt_plot2[,colmax >10]


png(filename="adt_heatmap_0429.png", width = 12, height = 15, units = "in", res = 300)
print(pheatmap::pheatmap(t(adt_plot3), 
                         color = viridis::viridis(25, option = "B"), 
                         fontsize_row = 8, border_color = NA))
dev.off()


adt_plot4=adt_plot2[,colmax >30]
png(filename="adt_heatmap_0429_v2.png", width = 12, height = 13, units = "in", res = 300)
print(pheatmap::pheatmap(t(adt_plot4), 
                         color = viridis::viridis(25, option = "B"), 
                         fontsize_row = 8, border_color = NA))
dev.off()

t4_only = adt_plot2[str_detect(rownames(adt_plot2), "T4"),]
adt_plot5=t4_only[,  apply(t4_only, 2, max) >15]
png(filename="adt_heatmap_0429_t4.png", width = 7, height = 7, units = "in", res = 300)
print(pheatmap::pheatmap(t(adt_plot5), 
                         color = viridis::viridis(25, option = "B"), 
                         fontsize_row = 8, border_color = NA))
dev.off()


t8_only = adt_plot2[str_detect(rownames(adt_plot2), "T8"),]
adt_plot5=t8_only[,  apply(t8_only, 2, max) >15]
png(filename="adt_heatmap_0429_t8.png", width = 7, height = 7, units = "in", res = 300)
print(pheatmap::pheatmap(t(adt_plot5), 
                         color = viridis::viridis(25, option = "B"), 
                         fontsize_row = 8, border_color = NA))
dev.off()


b_only = adt_plot2[str_detect(rownames(adt_plot2), "^B"),]
adt_plot5=b_only[,  apply(b_only, 2, max) >15]
png(filename="adt_heatmap_0429_b.png", width = 7, height = 8, units = "in", res = 300)
print(pheatmap::pheatmap(t(adt_plot5), 
                         color = viridis::viridis(25, option = "B"), 
                         fontsize_row = 8, border_color = NA))
dev.off()


mono_only = adt_plot2[str_detect(rownames(adt_plot2), "^cM"),]
adt_plot5=mono_only[,  apply(mono_only, 2, max) >20]
png(filename="adt_heatmap_0429_mono.png", width = 7, height = 9, units = "in", res = 300)
print(pheatmap::pheatmap(t(adt_plot5), 
                         color = viridis::viridis(25, option = "B"), 
                         fontsize_row = 8, border_color = NA))
dev.off()
# [ ] redo subset DSB dotplots

cd4_tcells = subset(sobj, cells=rownames(sobj@meta.data)[
  str_detect(sobj@meta.data$wnn_cluster_label2, "T4")])
DefaultAssay(cd4_tcells) = "integrated.ADT"
cd4_tcell_markers = FindAllMarkers(cd4_tcells)
save(cd4_tcell_markers, file="cd4_tcell_markers.RData")

cd4_tcell_gene = cd4_tcell_markers %>%
  filter(!str_detect(gene, "isotype")) %>%
  group_by(cluster) %>%
  slice_max(n=8, order_by = avg_log2FC) %>%
  pull(gene) 

cd4_tcells2 = subset(cd4_tcells, cells =
                       rownames(cd4_tcells@meta.data)[cd4_tcells@meta.data$wnn_cluster_label2!= "T4_Naive 2"] )
FuncPctDotPlot(cd4_tcells2,
               features=unique(cd4_tcell_gene),
               cols="RdYlBu", 
               cluster.idents=T)+
  coord_flip()+
  ylab("WNN cluster max annot")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("cd4_tcell_dotplot.pdf", height=7, width=7)



cd4_tnaive_cells = subset(sobj, cells=rownames(sobj@meta.data)[
  str_detect(sobj@meta.data$wnn_cluster_label2, "T4_Naive") & sobj@meta.data$wnn_cluster_label2!="T4_Naive 2"])
DefaultAssay(cd4_tnaive_cells) = "integrated.ADT"
cd4_tnaive_markers = FindAllMarkers(cd4_tnaive_cells)

cd4_tnaive_gene = cd4_tnaive_markers %>%
  filter(!str_detect(gene, "isotype")) %>%
  group_by(cluster) %>%
  slice_max(n=8, order_by = avg_log2FC) %>%
  arrange(cluster) %>%
  pull(gene) 

FuncPctDotPlot(cd4_tnaive_cells,
               features=unique(cd4_tnaive_gene),
               cols="RdYlBu", 
               cluster.idents=T)+
  coord_flip()+
  ylab("WNN cluster max annot")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("cd4_tnaive_dotplot.pdf", height=7, width=7)


cd8_tcells = subset(sobj, cells=rownames(sobj@meta.data)[
  str_detect(sobj@meta.data$wnn_cluster_label2, "T8")])
DefaultAssay(cd8_tcells) = "integrated.ADT"
cd8_tcell_markers = FindAllMarkers(cd8_tcells, only.pos=T)
save(cd8_tcell_markers, file="cd8_tcell_markers.RData")

cd8_tcell_gene = cd8_tcell_markers %>%
  filter(!str_detect(gene, "isotype")) %>%
  group_by(cluster) %>%
  slice_max(n=8, order_by = avg_log2FC) %>%
  pull(gene) 

cd8_tcells2 = subset(cd8_tcells, cells=rownames(cd8_tcells@meta.data)[cd8_tcells@meta.data$wnn_cluster_label2!="T8_MAIT"])
FuncPctDotPlot(cd8_tcells2,
               features=unique(cd8_tcell_gene),
               cols="RdYlBu", 
               cluster.idents=T)+
  coord_flip()+
  ylab("WNN cluster max annot")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("cd8_tcell_dotplot.pdf", height=7, width=7)


cd8_tnaive = subset(sobj, cells=rownames(sobj@meta.data)[
  str_detect(sobj@meta.data$wnn_cluster_label2, "T8_Naive")])
DefaultAssay(cd8_tnaive) = "integrated.ADT"
cd8_tnaive_markers = FindAllMarkers(cd8_tnaive, only.pos=T)

cd8_tnaive_gene =cd8_tnaive_markers %>%
  filter(!str_detect(gene, "isotype")) %>%
  group_by(cluster) %>%
  slice_max(n=8, order_by = avg_log2FC) %>%
  pull(gene) 

FuncPctDotPlot(cd8_tnaive,
               features=unique(cd8_tnaive_gene),
               cols="RdYlBu", 
               cluster.idents=T)+
  coord_flip()+
  ylab("WNN cluster max annot")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("cd8_tnaivecell_dotplot.pdf", height=7, width=7)


bcells = subset(sobj, cells=rownames(sobj@meta.data)[
  str_detect(sobj@meta.data$wnn_cluster_label2, "^B") & 
    sobj@meta.data$wnn_cluster_label2 != "B_Naive 3"])
DefaultAssay(bcells) = "integrated.ADT"
bcell_markers = FindAllMarkers(bcells, only.pos=T)
save(bcell_markers, file="bcell_markers.RData")

bcell_gene = bcell_markers %>%
  filter(!str_detect(gene, "isotype")) %>%
  group_by(cluster) %>%
  slice_max(n=8, order_by = avg_log2FC) %>%
  pull(gene) 


FuncPctDotPlot(bcells,
               features=unique(bcell_gene),
               cols="RdYlBu", 
               cluster.idents=T)+
  coord_flip()+
  ylab("WNN cluster max annot")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("bcell_dotplot.pdf", height=7, width=7)


mono = subset(sobj, cells=rownames(sobj@meta.data)[
  str_detect(sobj@meta.data$wnn_cluster_label2, "cM")])
DefaultAssay(mono) = "integrated.ADT"
mono_markers = FindAllMarkers(mono, only.pos=T)
save(mono_markers, file="mono_markers.RData")

mono_gene = mono_markers %>%
  filter(!str_detect(gene, "isotype")) %>%
  group_by(cluster) %>%
  slice_max(n=8, order_by = avg_log2FC) %>%
  pull(gene) 


FuncPctDotPlot(mono,
               features=unique(mono_gene),
               cols="RdYlBu", 
               cluster.idents=T)+
  coord_flip()+
  ylab("WNN cluster max annot")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("mono_dotplot.pdf", height=7, width=7)



cd4_tcells = subset(jdm_data, cells=rownames(jdm_data@meta.data)[
  str_detect(jdm_data@meta.data$wnn_cluster_label2, "T4")] )
DefaultAssay(cd4_tcells) = "integrated.ADT"
FuncPctDotPlot(cd4_tcells,
               features=c("CD4", "CD62L","CD45RA",
                          "CD194--CCR4","CD278--ICOS","CD122--IL-2Rbeta","CD127--IL-7Ralpha",
                          "rna_CCR7", "rna_SLC2A3" ),
               cols="RdYlBu", 
               cluster.idents=T)+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))+
  xlab("")+
  ylab("")
ggsave("cd4_tcell_selected_marker_dotplot.pdf")
rownames(cd4_tcells)





# [ ] try subclustering

# [ ] comparison dotplots for subsets RNA


### try soupX on one well

### update the DSB dotplots to wilcoxon



# submit JDM / cSLE integration
#  [x] LIGER
load("../csle_jdm/liger_out/int_data_seurat.RData")
load("../csle_jdm/liger_out/int_data.RData")

int.data <- runUMAP(int.data, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)
all.plots <- plotByDatasetAndCluster(int.data, 
                                     axis.labels = c('UMAP 1', 'UMAP 2'), 
                                     return.plots = T)
plot_grid(all.plots[[1]], all.plots[[2]])
ggsave("../csle_jdm/liger_out/post_int.png", height=5, width=11)


p1 = DimPlot(int.data.seurat, cells.highlight=
               rownames(int.data.seurat@meta.data)[int.data.seurat@meta.data$orig.ident=="csle"])+NoLegend()+
  ggtitle("csle")
p2 = DimPlot(int.data.seurat, cells.highlight=
               rownames(int.data.seurat@meta.data)[int.data.seurat@meta.data$orig.ident=="jdm"])+NoLegend()+
  ggtitle("jdm")
cowplot::plot_grid(p1, p2)
ggsave("../csle_jdm/liger_out/post_int_sep.png", height=5, width=11)


gene_loadings <- plotGeneLoadings(int.data, 
                                  do.spec.plot = FALSE, 
                                  return.plots = TRUE)
pdf("../csle_jdm/liger_out/gene_loadings.pdf")
gene_loadings
dev.off()


# tentative cluster labels
DATA_DIR = "/krummellab/data1/erflynn/premier/jdm/data/"

IN_DIR=sprintf("%s/csle_jdm/refPCA_seurat/", DATA_DIR)
load(sprintf("%s/processed_jdm.RData", IN_DIR))
load(sprintf("%s/processed_csle.RData", IN_DIR))

jdm_data = subset(jdm_data, case_control=="HC")
csle_data = subset(csle_data, Groups=="cHD")
csle_data@meta.data$orig.ident="csle"
jdm_data@meta.data$orig.ident="jdm"

csle_meta = csle_data@meta.data
jdm_meta = jdm_data@meta.data
rm(jdm_data)
rm(csle_data)
gc()

int_meta = int.data.seurat@meta.data %>%
  as_tibble(rownames="cell_lab") %>%
  select(cell_lab, cluster) %>%
  separate(cell_lab, into=c("orig.ident", "cell_id"), sep="_", extra="merge")

add_meta = jdm_meta %>% 
  as_tibble(rownames="cell_id") %>%
  select(cell_id, orig.ident, manual_labels) %>%
  rename(prev_label=manual_labels) %>%
  bind_rows(
    csle_meta %>% 
      as_tibble(rownames="cell_id") %>%
      select(cell_id, orig.ident, predicted.celltype.l2) %>%
      rename(prev_label=predicted.celltype.l2)
  )


int_meta2 = int_meta %>% 
  left_join(add_meta)

int_meta3 = int_meta2 %>%
  mutate(prev_label=case_when(
    prev_label == "NK_CD56bright" ~ "NK_CD56++", 
    prev_label == "CD14 Mono" ~ "cM",
    prev_label== "CD16 Mono" ~ "nCM",
    prev_label == "B memory" ~ "B_Mem",
    prev_label == "B naive" ~ "B_Naive",
    prev_label == "dnT" ~ "DN_T+",
    prev_label=="gdT" ~ "Tgd",
    prev_label=="Tgd1" ~ "Tgd",
    prev_label=="Tgd2" ~ "Tgd",
    prev_label=="ILC" ~ "ILC1_2_3",
    prev_label == "CD4 Naive" ~ "T4_Naive",
    prev_label=="CD8 Naive" ~ "T8_Naive",
    prev_label=="MAIT" ~ "T8_MAIT",
    prev_label == "CD8 TEM" ~ "T8_TEMRA",
    prev_label == "CD4 TCM" ~ "T4_Mem",
    prev_label == "CD4 CTL" ~ "T_Tox",
    prev_label == "CD4 TEM" ~ "T4_Mem",
    prev_label == "CD8 TCM" ~  "T8_Mem",
#    prev_label == "CD8 Proliferating" ~,
#    prev_label == "CD4 Proliferating" ~,
#    prev_label == "NK Proliferating" ~,
    prev_label == "Plasmablast" ~ "PB",
    prev_label == "HSPC" ~ "HSC" ,
 #   prev_label == "ASDC" ~,
    prev_label == "B intermediate" ~ "B_Mem_intermediate",
TRUE ~ prev_label
  ))

clus_to_lab = int_meta3 %>%
  group_by(cluster) %>%
  mutate(tot=n()) %>%
  group_by(cluster, prev_label, tot) %>%
  count() %>%
  ungroup() %>%
  group_by(cluster, tot) %>%
  mutate(frac=n/tot) %>%
  arrange(cluster, desc(n)) %>%
  slice_max(1)  %>%
  mutate(prev_label=case_when(cluster == 3 ~ "T4_Mem 2",
                   cluster==8 ~ "B_Naive 2",
                   TRUE ~ prev_label)) %>%
  ungroup()

int.data.seurat = AddMetaData(int.data.seurat, int_meta3 %>% pull(prev_label), "prev_label")


int.data.seurat = AddMetaData(int.data.seurat, int.data.seurat@meta.data %>%
  left_join(clus_to_lab %>% select(cluster, prev_label) %>% rename(cluster_label=prev_label)) %>%
    pull(cluster_label), "cluster_label")

int_meta4 = int.data.seurat@meta.data
clus_fracs = int_meta4 %>%
  group_by(orig.ident) %>%
  mutate(ntot=n()) %>%
  group_by(cluster_label, orig.ident, ntot) %>%
  count() %>%
  mutate(frac=n/ntot) %>%
  arrange(desc(frac)) 

clus_fracs$cluster_label = factor(clus_fracs$cluster_label, levels=unique(clus_fracs$cluster_label))
clus_fracs %>%
  ggplot(aes(x=cluster_label, y=frac, fill=orig.ident))+
  geom_bar(stat="identity", position="dodge")+
  theme_bw()+
  ylab("fraction")+
  labs(fill="dataset")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.1))
ggsave("../csle_jdm/liger_out/clus_fracs.png")

Idents(int.data.seurat) = "cluster_label"
p0 = DimPlot(int.data.seurat, cols=dittoColors(), label=T)+NoLegend()
p1 =FeaturePlot(int.data.seurat, features="CD3D")
p2 = FeaturePlot(int.data.seurat, features="CD14")
p3 = FeaturePlot(int.data.seurat, features="CD19")
plot_grid(p0, p1, p2, p3)
ggsave("../csle_jdm/liger_out/clus_lab.png")

#  [x] REF - looks awful, but mb gene overlap
#  [x] RPCA
load("../csle_jdm/rpca_seurat/csle_jdm_integ.RData")
load("../csle_jdm/rpca_seurat/integ_clus.RData")
sobj_integ = AddMetaData(sobj_integ, 
                         sobj_integ@meta.data$RNA_snn_res.1, 
                         "integ_clus")
unique(sobj_integ@meta.data$RNA_snn_res.1)

Idents(sobj_integ) = "integ_clus"
DimPlot(sobj_integ)

p1 = DimPlot(sobj_integ, cells.highlight=
          rownames(sobj_integ@meta.data)[sobj_integ@meta.data$orig.ident=="csle"])+NoLegend()+
  ggtitle("csle")
p2 = DimPlot(sobj_integ, cells.highlight=
               rownames(sobj_integ@meta.data)[sobj_integ@meta.data$orig.ident=="jdm"])+NoLegend()+
  ggtitle("jdm")
cowplot::plot_grid(p1, p2)
ggsave("integ_rpca_after.png", height=5, width=10)
