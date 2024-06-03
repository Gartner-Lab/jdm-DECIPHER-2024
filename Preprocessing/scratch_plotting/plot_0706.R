library(tidyverse)
library(Seurat)

JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm/"
CODE_DIR="/krummellab/data1/erflynn/premier/jdm"
source(sprintf("%s/scripts/utils/label_transfer.R", CODE_DIR))
source(sprintf("%s/scripts/utils/plot_wnn.R", CODE_DIR))

# load the data
OUT_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm/wnn_v2/"
setwd(OUT_DIR)

load("wnn_obj2.RData")
head(sobj@meta.data)


marker.res=1.4
load(sprintf("wnn_rna_markers_%s.RData", marker.res))
load(sprintf("wnn_adt_markers_%s.RData",  marker.res))

wnn_markers %>%
  group_by(cluster) %>%
  arrange(cluster, p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n=50) %>%
  write_csv(sprintf("top50_rna_markers_per_clus_%s.csv", marker.res))

wnn_markers_adt %>%
  group_by(cluster) %>%
  arrange(cluster, p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n=50) %>% 
  rename(wilcoxon_effect_size=avg_log2FC)  %>% 
  write_csv(sprintf("top50_adt_markers_per_clus_%s.csv", marker.res))


## add the labels
# manual_labels
sobj = sobj %>% add_manual_labels() 


# azimuth labels
load("../azimuth_jdm.RData")
azimuth_meta2 = azimuth_meta %>%
    map_azimuth_label(predicted.celltype.l2, azimuth_label_mapped) %>%
    as_tibble(rownames="cell_id") %>%
    select(cell_id, azimuth_label_mapped)
sobj = addMetaJoin(sobj, azimuth_meta2, "azimuth_label_mapped")
sobj@meta.data

sobj = transfer_annot(sobj, "wsnn_res.1.4", "wnn_cluster_label")

sobj@meta.data$manual_labels0 = sobj@meta.data$manual_labels
sobj@meta.data$manual_labels = sobj@meta.data$azimuth_label_mapped

sobj = transfer_annot(sobj, "wsnn_res.1.4", "azimuth_cluster_label")
sobj@meta.data$manual_labels = sobj@meta.data$manual_labels0
sobj@meta.data$manual_labels0 = NULL

Idents(sobj) = "wnn_cluster_label"

## save
save(sobj, file="sobj_w_lab.RData")

## remove a few clusters
unique(sobj@active.ident)
dim(sobj)
cells.keep = sobj@meta.data %>%
  as_tibble(rownames="cell_id" ) %>%
  filter(!wnn_cluster_label %in%  c("T4_Naive* 1", "T4_Naive* 2")) %>%
  pull(cell_id)
length(cells.keep)

sobj_filt = subset(sobj, cells=cells.keep)

## add "broad_ids"
unique(RMs.seu@meta.data$broad_ids)
# "B"     "myeloid"   "CD4T"   "NK"    "gdT"   "CD8T"   "platelets"

broad_ids = RMs.seu@meta.data %>%
  as_tibble(rownames="cell_id") %>%
  select(cell_id, wnn_cluster_label, broad_ids)

broad_ids %>% 
  group_by(broad_ids, wnn_cluster_label) %>% 
  count() %>%
  arrange(wnn_cluster_label, desc(n)) %>%
  View()

unique(sobj_filt@active.ident)

sobj_filt@meta.data = sobj_filt@meta.data %>%
  mutate(broad_ids=case_when(
    str_detect(wnn_cluster_label, "^B_") ~ "B",
    wnn_cluster_label=="Tgd1" | wnn_cluster_label=="T8_MAIT" ~ "gdT",
    str_detect(wnn_cluster_label, "^T8_") ~ "CD8T",
    str_detect(wnn_cluster_label, "^T4_") ~ "CD4T",
    str_detect(wnn_cluster_label, "^NK") ~ "NK",
    str_detect(wnn_cluster_label, "cM") ~ "myeloid",
    wnn_cluster_label == "PB" ~ "B",
    wnn_cluster_label == "pDC" ~ "myeloid",
    wnn_cluster_label == "cDC2" ~ "myeloid"
  )) 

## save diet version
DefaultAssay(sobj_filt) = "RNA"
dsobj = DietSeurat(sobj_filt, assay="RNA")

saveRDS(dsobj, file="diet_sobj.RDS", compress=F)

## standard plots
col_vec = dittoColors()[1:length(unique(sobj@meta.data$wnn_cluster_label))]
names(col_vec) = unique(sobj@meta.data$wnn_cluster_label)


p1 = DimPlot(sobj, group.by="wsnn_res.1.4",reduction="wnn.umap", cols=unname(col_vec), label=T,
        label.size=3, repel=T)+NoLegend()
p2 = DimPlot(sobj, group.by="manual_labels", reduction="wnn.umap", cols=unname(col_vec), label=T,
        label.size=3, repel=T)+NoLegend()
p3 = DimPlot(sobj, group.by="azimuth_label_mapped", reduction="wnn.umap", cols=unname(col_vec), label=T,
        label.size=3, repel=T)+NoLegend()+ggtitle("azimuth_label")
p4 = DimPlot(sobj, group.by="wnn_cluster_label", reduction="wnn.umap", cols=col_vec, label=T,
        label.size=3, repel=T)+NoLegend()+ggtitle("wnn_max_cluster")
plot_grid(p1, p2, p3, p4, ncol=2)
ggsave(sprintf("%s/dimplot_labels.pdf", OUT_DIR), height=13, width=12)

DimPlot(sobj, group.by="azimuth_cluster_label", reduction="wnn.umap", cols=unname(col_vec), label=T,
        label.size=3, repel=T)+NoLegend()

# azimuth max per cluster geom_tile
sobj@meta.data %>% 
  group_by(azimuth_cluster_label) %>% 
  mutate(tot=n()) %>%
  group_by(azimuth_cluster_label, wnn_cluster_label, tot) %>%
  count() %>%
  mutate(frac=n/tot) %>%
  ggplot(aes(x=azimuth_cluster_label, y=wnn_cluster_label, fill=frac))+
  geom_tile()+
  theme_bw()+
  scale_color_gradient(low="white", high="blue")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.1))

sobj@meta.data %>% 
  group_by(manual_labels) %>% 
  mutate(tot=n()) %>%
  group_by(manual_labels, azimuth_label_mapped, tot) %>%
  count() %>%
  mutate(frac=n/tot) %>%
  ungroup() %>%
  select(-tot, -n) %>%
  mutate(manual_labels=ifelse(is.na(manual_labels), "unlabeled", manual_labels)) %>%
  pivot_wider(names_from="manual_labels", values_from="frac", values_fill=0) %>%
  pivot_longer(-azimuth_label_mapped, names_to="manual_labels", 
               values_to="frac") %>%
  ggplot(aes(x=azimuth_label_mapped, y=manual_labels, fill=frac))+
  labs(fill="fraction")+
  geom_tile()+
  theme_bw()+
  scale_fill_gradient(low="white", high="blue")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.3))
ggsave(sprintf("%s/azimuth_vs_manual.pdf", OUT_DIR))


VlnPlot(sobj, features=c( "CD4", "CD8", "CD16", "CD45RO", "CD45RA", "CD56--NCAM",
                          "CD14", "rna_CD14", "rna_FCGR3A"), pt.size=0,
        cols=col_vec)
ggsave(sprintf("%s/vln_plots.pdf", OUT_DIR), height=15, width=19)

make_wnn_plots(sobj, col_vec)
ggsave(sprintf("%s/wnn_plots.pdf", OUT_DIR), width=15, height=5)

make_wnn_wt_plot(sobj, col_vec)
ggsave(sprintf("%s/wnn_wt_plot.pdf", OUT_DIR), width=9, height=5)

# dimplot for RNA

# dimplot

# donor/well plots
DimPlot(sobj, group.by="donor", cols=dittoColors(), reduction="wnn.umap")
ggsave(sprintf("%s/donor_dimplot.pdf", OUT_DIR), width=6, height=5)

DimPlot(sobj, group.by="well", cols=dittoColors(), reduction="wnn.umap")
ggsave(sprintf("%s/well_dimplot.pdf", OUT_DIR), width=6, height=5)

make_meta_barplot(sobj, donor, "wnn_cluster_label")
ggsave(sprintf("%s/donor_barplot.pdf", OUT_DIR), width=11, height=8)

make_meta_barplot(sobj, well, "wnn_cluster_label")
ggsave(sprintf("%s/well_barplot.pdf", OUT_DIR), width=10, height=5)


# vln plot qc + doublet check
make_doublet_vln_plots(sobj, col_vec, "wnn_cluster_label")
ggsave(sprintf("%s/cd_doublet_check.pdf", OUT_DIR), height=20, width=15)

make_qc_vln_plots(sobj, col_vec, "wnn_cluster_label")
ggsave(sprintf("%s/qc_vln.png", OUT_DIR), height=15, width=20)

# cluster size
plot_clus_sizes(sobj, col_vec, "wnn_cluster_label")
ggsave(sprintf("%s/wnn_clus_size.pdf", OUT_DIR), width=8, height=5)


# dotplots
gen_dotplot(sobj, "integrated.ADT", wnn_markers_adt)
ggsave(sprintf("%s/adt_dotplot.pdf", OUT_DIR), width=11, height=23)

gen_dotplot(sobj, "RNA", wnn_markers)
ggsave(sprintf("%s/rna_dotplot.pdf", OUT_DIR), width=11, height=20)



# cluster diffuse-ness?
split_clus = sobj@meta.data %>%
  as_tibble(rownames="cell_id") %>%
  select(cell_id, wnn_cluster_label) %>%
  group_split(wnn_cluster_label)
plts = lapply(split_clus, function(x){
  DimPlot(sobj, cells.highlight = x$cell_id, reduction="wnn.umap")+
    NoLegend()+
    ggtitle(sprintf("%s (n=%s)",x$wnn_cluster_label[[1]], nrow(x)))
})
plot_grid(plotlist=plts)
ggsave("cluster_diffuse.pdf", height=20, width=20)

# subcluster plots
cd4_tcells = subset(sobj, cells=rownames(sobj@meta.data)[
  str_detect(sobj@meta.data$wnn_cluster_label, "T4")])
DefaultAssay(cd4_tcells) = "integrated.ADT"
cd4_tcell_markers = FindAllMarkers(cd4_tcells)
save(cd4_tcell_markers, file="cd4_tcell_markers.RData")

cd4_tcell_gene = cd4_tcell_markers %>%
  filter(!str_detect(gene, "isotype")) %>%
  group_by(cluster) %>%
  slice_max(n=8, order_by = avg_log2FC) %>%
  pull(gene) 

#cd4_tcells2 = subset(cd4_tcells, cells =
#                       rownames(cd4_tcells@meta.data)[cd4_tcells@meta.data$wnn_cluster_label!= "T4_Naive 2"] )
FuncPctDotPlot(cd4_tcells,
               features=unique(cd4_tcell_gene),
               cols="RdYlBu", 
               cluster.idents=T)+
  coord_flip()+
  ylab("WNN cluster max annot")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("cd4_tcell_dotplot.pdf", height=10, width=10)



#cd4_tnaive_cells = subset(sobj, cells=rownames(sobj@meta.data)[
#  str_detect(sobj@meta.data$wnn_cluster_label, "T4_Naive") & 
#    sobj@meta.data$wnn_cluster_label!="T4_Naive 2"])
#DefaultAssay(cd4_tnaive_cells) = "integrated.ADT"
#cd4_tnaive_markers = FindAllMarkers(cd4_tnaive_cells)

#cd4_tnaive_gene = cd4_tnaive_markers %>%
#  filter(!str_detect(gene, "isotype")) %>%
#  group_by(cluster) %>%
#  slice_max(n=8, order_by = avg_log2FC) %>%
#  arrange(cluster) %>%
#  pull(gene) 

#FuncPctDotPlot(cd4_tnaive_cells,
#               features=unique(cd4_tnaive_gene),
#               cols="RdYlBu", 
#               cluster.idents=T)+
#  coord_flip()+
#  ylab("WNN cluster max annot")+
#  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
#ggsave("cd4_tnaive_dotplot.pdf", height=7, width=7)


cd8_tcells = subset(sobj, cells=rownames(sobj@meta.data)[
  str_detect(sobj@meta.data$wnn_cluster_label, "T8")])
DefaultAssay(cd8_tcells) = "integrated.ADT"
cd8_tcell_markers = FindAllMarkers(cd8_tcells, only.pos=T)
save(cd8_tcell_markers, file="cd8_tcell_markers.RData")

cd8_tcell_gene = cd8_tcell_markers %>%
  filter(!str_detect(gene, "isotype")) %>%
  group_by(cluster) %>%
  slice_max(n=8, order_by = avg_log2FC) %>%
  pull(gene) 

#cd8_tcells2 = subset(cd8_tcells, cells=rownames(cd8_tcells@meta.data)[cd8_tcells@meta.data$wnn_cluster_label!="T8_MAIT"])
FuncPctDotPlot(cd8_tcells,
               features=unique(cd8_tcell_gene),
               cols="RdYlBu", 
               cluster.idents=T)+
  coord_flip()+
  ylab("WNN cluster max annot")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("cd8_tcell_dotplot.pdf", height=7, width=7)


# cd8_tnaive = subset(sobj, cells=rownames(sobj@meta.data)[
#   str_detect(sobj@meta.data$wnn_cluster_label, "T8_Naive")])
# DefaultAssay(cd8_tnaive) = "integrated.ADT"
# cd8_tnaive_markers = FindAllMarkers(cd8_tnaive, only.pos=T)
# 
# cd8_tnaive_gene =cd8_tnaive_markers %>%
#   filter(!str_detect(gene, "isotype")) %>%
#   group_by(cluster) %>%
#   slice_max(n=8, order_by = avg_log2FC) %>%
#   pull(gene) 
# 
# FuncPctDotPlot(cd8_tnaive,
#                features=unique(cd8_tnaive_gene),
#                cols="RdYlBu", 
#                cluster.idents=T)+
#   coord_flip()+
#   ylab("WNN cluster max annot")+
#   theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
# ggsave("cd8_tnaivecell_dotplot.pdf", height=7, width=7)


bcells = subset(sobj, cells=rownames(sobj@meta.data)[
  str_detect(sobj@meta.data$wnn_cluster_label, "^B")])
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
  str_detect(sobj@meta.data$wnn_cluster_label, "cM")])
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


