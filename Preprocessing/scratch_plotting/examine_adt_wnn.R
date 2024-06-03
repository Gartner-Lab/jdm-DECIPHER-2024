# 6/3/2022
library(Seurat)
library(tidyverse)
library(dittoSeq)
# Look at the ADT DSB data
# -- can we resolve the problem of some clusters being low ADT?

CODE_DIR="/krummellab/data1/erflynn/premier/jdm/"
JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm"

source(sprintf("%s/scripts/utils/label_transfer.R", CODE_DIR))
source(sprintf("%s/scripts/utils/plot_wnn.R", CODE_DIR))

# -  original
OUT_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm/wnn_soupx/"
load(sprintf("%s/wnn_obj.RData", OUT_DIR))

# add labels
sobj = sobj %>% add_manual_labels() %>% add_azimuth_labels()
sobj = transfer_annot(sobj, "wsnn_res.1.4", "wnn_cluster_label")
meta = sobj@meta.data

Idents(sobj) = "wnn_cluster_label"
DimPlot(sobj, label=T, label.size=3, repel=T, cols=dittoColors(), reduction="wnn.umap")+NoLegend()
ggsave(sprintf("%s/orig_labeled_dimplot.pdf", OUT_DIR), height=5, width=5)

VlnPlot(sobj, "RNA.weight", pt.size=0, sort=T)+NoLegend()

VlnPlot(sobj, "nCount_ADT", pt.size=0, sort=T)+NoLegend()
VlnPlot(sobj, "nFeature_ADT", pt.size=0, sort=T)+NoLegend()

t4_naive_cells = rownames(sobj@meta.data)[str_detect(sobj@meta.data$wnn_cluster_label, "^T4_Naive")]
t4_naive = subset(sobj, cells=t4_naive_cells)

col_vec = dittoColors(length(unique(t4_naive$wnn_cluster_label)))
names(col_vec) = unique(t4_naive$wnn_cluster_label)
DimPlot(t4_naive, label=T, reduction="wnn.umap", repel=T, 
        cols=col_vec)+NoLegend()
ggsave(sprintf("%s/orig_t4_naive_dimplot.pdf", OUT_DIR), height=5, width=5)

# cells highlight 
# T4_Naive 3 is the biggest problem cluster
meta %>% select(wnn_cluster_label, well, nFeature_ADT) %>% 
  filter(wnn_cluster_label=="T4_Naive 3") %>%
  group_by(well) %>% 
  summarize(m=median(nFeature_ADT),
            max=max(nFeature_ADT),
            min=min(nFeature_ADT),
            n=n()) #%>%
 # count() # pretty evenly distributed in terms of count, and all have low mins


p1 = VlnPlot(t4_naive, "RNA.weight", pt.size=0, sort=T, cols=col_vec)+NoLegend()
p2 = VlnPlot(t4_naive, "nCount_ADT", pt.size=0, sort=T, cols=col_vec)+NoLegend()
p3 = VlnPlot(t4_naive, "nFeature_ADT", pt.size=0, sort=T, cols=col_vec)+NoLegend()
plot_grid(p3, p2, p1)
ggsave(sprintf("%s/orig_t4_naive_vln.pdf", OUT_DIR), height=6.5, width=5)

t4_naive_markers = FindAllMarkers(t4_naive, assay="integrated.ADT", only.pos = T)
gen_dotplot(t4_naive, "integrated.ADT", t4_naive_markers)
ggsave(sprintf("%s/orig_t4_naive_adt_dotplot.pdf", OUT_DIR), height=9, width=10)


#### look at before & after DSB, & in CLR

# could we simplify / replicate by just looking at one well? e.g. well 3


t4_naive3_cells = meta %>% 
  as_tibble(rownames="cell_id") %>%
  filter(wnn_cluster_label=="T4_Naive 3") %>%
  pull(cell_id)

other_t4_cells = meta %>%
  as_tibble(rownames="cell_id") %>%
  filter(str_detect(wnn_cluster_label, "^T4_Naive")) %>%
  filter(wnn_cluster_label!= "T4_Naive 3") %>%
  pull(cell_id)

other_cells = setdiff(rownames(meta), union(other_t4_cells, t4_naive3_cells))

# merged_processed.RData - has levels before DSB?
load("/krummellab/data1/erflynn/premier/jdm/data/jdm/merged_jdm_soupx/merged_processed.RData",
     verbose=T) # --> merged_data
t4_naive3 = subset(merged_data, cells=t4_naive3_cells)
other_t4 = subset(merged_data, cells=other_t4_cells)
other_sobj = subset(merged_data, cells=other_cells)

rm(merged_data); gc()

t4_naive3_adt = t4_naive3@assays$ADT@counts
other_t4_adt = other_t4@assays$ADT@counts
other_adt = other_sobj@assays$ADT@counts

isotype_ctls = rownames(t4_naive3_adt)[str_detect(rownames(t4_naive3_adt), "Iso")]

t4_naive3_isotype = as.matrix(t4_naive3_adt[isotype_ctls,])
other_t4_isotype = as.matrix(other_t4_adt[isotype_ctls,])
other_adt_isotype = as.matrix(other_adt[isotype_ctls,])



apply(t4_naive3_isotype, 1, max)
quantile(as.vector(t4_naive3_isotype), c(0.99, 0.995, 0.999, 0.9995, 0.9999, 1))
apply(t4_naive3_isotype, 1, mean)
dim(t4_naive3_isotype)

dim(t4_naive3_isotype[,apply(t4_naive3_isotype, 2, max) < 100])
t4_naive3_isotype[,apply(t4_naive3_isotype, 2, max) > 100] # only nine cells

ncol(t4_naive3_isotype[,apply(t4_naive3_isotype, 2, max) > 50]) # 20

apply(other_t4_isotype, 1, max) 
apply(other_t4_isotype, 1, mean)
dim(other_t4_isotype)
dim(other_t4_isotype[,apply(other_t4_isotype, 2, max) < 100])
other_t4_isotype[,apply(other_t4_isotype, 2, max) > 100] # only three cells
ncol(other_t4_isotype[,apply(other_t4_isotype, 2, max) > 50]) # 16
quantile(as.vector(other_t4_isotype), c(0.99, 0.995, 0.999, 0.9995, 0.9999, 1))


apply(other_adt_isotype, 1, max) 
apply(other_adt_isotype, 1, mean)
dim(other_adt_isotype)
quantile(as.vector(other_adt_isotype), c(0.99, 0.995, 0.999, 0.9995, 0.9999, 1))
         
dim(other_adt_isotype[,apply(other_adt_isotype, 2, max) < 100])
other_adt_isotype[,apply(other_adt_isotype, 2, max) > 100] # only seven cells
ncol(other_adt_isotype[,apply(other_adt_isotype, 2, max) > 50]) # 37


### what are the cutoffs the DSB paper uses for ADT filtering
# prot.size = log10(Matrix::colSums(prot))
# do we need to filter the background for ADT?


# DSB: levels post-DSB
# rpca_orig: levels post RPCA




# - after regressing out
OUT_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm/wnn_soupx2"
load(sprintf("%s/wnn_obj2.RData", OUT_DIR))

# add labels
sobj = sobj %>% add_manual_labels() %>% add_azimuth_labels()
sobj = transfer_annot(sobj, "wsnn_res.1.2", "wnn_cluster_label")
meta = sobj@meta.data

Idents(sobj) = "wnn_cluster_label"
DimPlot(sobj, label=T, label.size=3, repel=T, cols=dittoColors(), reduction="wnn.umap")+NoLegend()
ggsave(sprintf("%s/labeled_dimplot.pdf", OUT_DIR), height=5, width=5)

VlnPlot(sobj, "RNA.weight", pt.size=0, sort=T)+NoLegend()

VlnPlot(sobj, "nCount_ADT", pt.size=0, sort=T)+NoLegend()
VlnPlot(sobj, "nFeature_ADT", pt.size=0, sort=T)+NoLegend()

t4_naive_cells = rownames(sobj@meta.data)[str_detect(sobj@meta.data$wnn_cluster_label, "^T4_Naive")]
t4_naive = subset(sobj, cells=t4_naive_cells)
library(dittoSeq)
col_vec = dittoColors(7)
names(col_vec) = unique(t4_naive$wnn_cluster_label)
DimPlot(t4_naive, label=T, reduction="wnn.umap", repel=T, 
        cols=col_vec)+NoLegend()
ggsave(sprintf("%s/t4_naive_dimplot.pdf", OUT_DIR), height=5, width=5)

# cells highlight


p1 = VlnPlot(t4_naive, "RNA.weight", pt.size=0, sort=T, cols=col_vec)+NoLegend()
p2 = VlnPlot(t4_naive, "nCount_ADT", pt.size=0, sort=T, cols=col_vec)+NoLegend()
p3 = VlnPlot(t4_naive, "nFeature_ADT", pt.size=0, sort=T, cols=col_vec)+NoLegend()
plot_grid(p3, p2, p1)
ggsave(sprintf("%s/t4_naive_vln.pdf", OUT_DIR))

t4_naive_markers = FindAllMarkers(t4_naive, assay="integrated.ADT", only.pos = T)
gen_dotplot(t4_naive, "integrated.ADT", t4_naive_markers)
ggsave(sprintf("%s/t4_naive_adt_dotplot.pdf", OUT_DIR), height=9, width=10)



# dotplot


# - after lower bound nFeature --> reclustering
OUT_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm/wnn_soupx/"
load(sprintf("%s/wnn_reclus.RData", OUT_DIR))

# add labels
sobj = sobj %>% add_manual_labels() %>% add_azimuth_labels()
sobj = transfer_annot(sobj, "wsnn_res.1.4", "wnn_cluster_label")
meta = sobj@meta.data

Idents(sobj) = "wnn_cluster_label"
DimPlot(sobj, label=T, label.size=3, repel=T, cols=dittoColors(), reduction="wnn.umap")+NoLegend()
ggsave(sprintf("%s/reclus_labeled_dimplot.pdf", OUT_DIR), height=5, width=5)

VlnPlot(sobj, "RNA.weight", pt.size=0, sort=T)+NoLegend()

VlnPlot(sobj, "nCount_ADT", pt.size=0, sort=T)+NoLegend()
VlnPlot(sobj, "nFeature_ADT", pt.size=0, sort=T)+NoLegend()

t4_naive_cells = rownames(sobj@meta.data)[str_detect(sobj@meta.data$wnn_cluster_label, "^T4_Naive")]
t4_naive = subset(sobj, cells=t4_naive_cells)

col_vec = dittoColors(length(unique(t4_naive$wnn_cluster_label)))
names(col_vec) = unique(t4_naive$wnn_cluster_label)
DimPlot(t4_naive, label=T, reduction="wnn.umap", repel=T, 
        cols=col_vec)+NoLegend()
ggsave(sprintf("%s/reclus_t4_naive_dimplot.pdf", OUT_DIR), height=5, width=5)

# cells highlight


p1 = VlnPlot(t4_naive, "RNA.weight", pt.size=0, sort=T, cols=col_vec)+NoLegend()
p2 = VlnPlot(t4_naive, "nCount_ADT", pt.size=0, sort=T, cols=col_vec)+NoLegend()
p3 = VlnPlot(t4_naive, "nFeature_ADT", pt.size=0, sort=T, cols=col_vec)+NoLegend()
plot_grid(p3, p2, p1)
ggsave(sprintf("%s/reclus_t4_naive_vln.pdf", OUT_DIR), height=6.5, width=5)

t4_naive_markers = FindAllMarkers(t4_naive, assay="integrated.ADT", only.pos = T)
gen_dotplot(t4_naive, "integrated.ADT", t4_naive_markers)
ggsave(sprintf("%s/reclus_t4_naive_adt_dotplot.pdf", OUT_DIR), height=9, width=10)




