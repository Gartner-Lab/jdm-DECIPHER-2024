library(Seurat)
library(tidyverse)
library(dittoSeq)
JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm"
norm.method="DSB"

OUT_DIR = sprintf("%s/wnn_out_rpca_%s3/", JDM_DIR, norm.method)
setwd(OUT_DIR)
load("wnn_obj.RData")
Idents(sobj) = "wsnn_res.1.4"

VlnPlot(sobj, features = "RNA.weight", 
        group.by = 'wsnn_res.1.4', sort = TRUE, pt.size = 0) +
  NoLegend()
ggsave("wnn_weight1.4.png")


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


col_vec = dittoColors()[1:length(unique(sobj@meta.data$cluster_label))]
names(col_vec)= unique(sobj@meta.data$cluster_label)

Idents(sobj) = "cluster_label"
## remake t
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



load(sprintf("%s/wnn_rna_1.4.RData", OUT_DIR)) # wnn_markers
load(sprintf("%s/wnn_adt_1.4.RData", OUT_DIR)) # wnn_markers_adt


top5_adt_genes = wnn_markers_adt %>% 
  left_join(clus_map) %>%
  group_by(cluster_label) %>%
  slice_max(n=5, order_by = avg_log2FC) %>%
  arrange(cluster_label)
top5_adt_genes$gene = factor(top5_adt_genes$gene, levels=unique(top5_adt_genes$gene))
source("/krummellab/data1/erflynn/scratch_scripts/FuncPctDotPlot.R")
FuncPctDotPlot(sobj, assay="integrated.ADT", features = unique(top5_adt_genes$gene),
               cols="RdYlBu",
               cluster.idents = T)+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave(sprintf("%s/adt_dotplot.pdf", OUT_DIR), width=14, height=17)

top5_genes = wnn_markers %>% 
  left_join(clus_map) %>%
  group_by(cluster_label) %>%
  slice_max(n=5, order_by = avg_log2FC) %>%
  pull(gene)
DotPlot(sobj, assay="RNA", features=unique(top5_genes),
        cols="RdYlBu", cluster.idents = T)+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave(sprintf("%s/rna_dotplot.pdf", OUT_DIR), width=11, height=20)




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

## ALLUVIAL
# alluvial for data changes
library(ggalluvial)

meta = sobj@meta.data  %>%
  as_tibble(rownames="cell_id") 

meta_alluv = meta %>% 
  select(cell_id, manual_labels, cluster_label) %>%
  dplyr::rename(RNA=manual_labels,
                WNN=cluster_label) %>%
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
ggsave("alluv_annot.pdf", height=7, width=10)



## SMALL DOTPLOT

cd4_tcells = subset(sobj, cells=rownames(sobj@meta.data)[
  str_detect(sobj@meta.data$cluster_label, "T4")])
DefaultAssay(cd4_tcells) = "integrated.ADT"
cd4_tcell_markers = FindAllMarkers(cd4_tcells, test.use="MAST", only.pos=T)
save(cd4_tcell_markers, file="cd4_tcell_markers.RData")

cd4_tcell_gene = cd4_tcell_markers %>%
  filter(!str_detect(gene, "isotype")) %>%
  group_by(cluster) %>%
  slice_max(n=8, order_by = avg_log2FC) %>%
  pull(gene) 


FuncPctDotPlot(cd4_tcells,
               features=unique(cd4_tcell_gene),
               cols="RdYlBu", 
               cluster.idents=T)+
  coord_flip()+
  ylab("WNN cluster max annot")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("cd4_tcell_dotplot.pdf", height=7, width=7)


cd8_tcells = subset(sobj, cells=rownames(sobj@meta.data)[
  str_detect(sobj@meta.data$cluster_label, "T8")])
DefaultAssay(cd8_tcells) = "integrated.ADT"
cd8_tcell_markers = FindAllMarkers(cd8_tcells, test.use="MAST", only.pos=T)
save(cd8_tcell_markers, file="cd8_tcell_markers.RData")

cd8_tcell_gene = cd8_tcell_markers %>%
  filter(!str_detect(gene, "isotype")) %>%
  group_by(cluster) %>%
  slice_max(n=8, order_by = avg_log2FC) %>%
  pull(gene) 


FuncPctDotPlot(cd8_tcells,
               features=unique(cd8_tcell_gene),
               cols="RdYlBu", 
               cluster.idents=T)+
  coord_flip()+
  ylab("WNN cluster max annot")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("cd8_tcell_dotplot.pdf", height=7, width=7)


bcells = subset(sobj, cells=rownames(sobj@meta.data)[
  str_detect(sobj@meta.data$cluster_label, "^B") & sobj@meta.data$cluster_label != "B_Naive 3"])
DefaultAssay(bcells) = "integrated.ADT"
bcell_markers = FindAllMarkers(bcells, test.use="MAST", only.pos=T)
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
  str_detect(sobj@meta.data$cluster_label, "cM")])
DefaultAssay(mono) = "integrated.ADT"
mono_markers = FindAllMarkers(mono, test.use="MAST", only.pos=T)
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


## REMOVE + RECLUSTER
cells.keep = sobj@meta.data %>%
  as_tibble(rownames="cell_id") %>%
  filter(!cluster_label %in% c("unknown 3",
    "unknown 5",
    "unknown 7",
    "NKs_mix_prolif",
    "B_Naive 3")) %>%
  pull(cell_id)
sobj_filt = subset(sobj, cells=cells.keep)

DimPlot(sobj_filt, reduction = 'wnn.umap', 
        label = TRUE, cols=col_vec,
        label.size=3) + NoLegend()
ggsave(sprintf("%s/wnn_umap_clus_lab_filt.pdf", OUT_DIR), 
       height=5, width=5)

DefaultAssay(sobj_filt) = "RNA"
sobj_filt = FindMultiModalNeighbors(
  sobj_filt, reduction.list=list("harmony", "a.pca"), ## TODO rename these and switch to harmony, aharmony
  dims.list=list(1:30, 1:18), # TODO: should this be different?
  prune.SNN = 1/20, 
  modality.weight.name="RNA.weight"
)

sobj_filt = RunUMAP(sobj_filt, nn.name = "weighted.nn", 
               reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

sobj_filt = FindClusters(sobj_filt, graph.name = "wsnn", 
                    algorithm = 3 , resolution = 1.4, verbose = FALSE)
save(sobj_filt, file="sobj_filt.RData")

## PRE/POST HARMONY MIXING
