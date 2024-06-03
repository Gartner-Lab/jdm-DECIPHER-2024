library(tidyverse)
library(Seurat)

# write out the markers files

JDM_DIR = "/krummellab/data1/erflynn/premier/jdm/data/jdm"

marker.res = 1.4
soupx=F
if (soupx){
  OUT_DIR=sprintf("%s/wnn_soupx/", JDM_DIR)
} else {
  OUT_DIR=sprintf("%s/wnn/", JDM_DIR)
}


load(sprintf("%s/wnn_rna_markers_%s.RData", OUT_DIR, marker.res))
load(sprintf("%s/wnn_adt_markers_%s.RData", OUT_DIR, marker.res))

wnn_markers %>%
  group_by(cluster) %>%
  arrange(cluster, p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n=50) %>%
  write_csv(sprintf("%s/top50_rna_markers_per_clus_%s.csv", OUT_DIR, marker.res))
  

wnn_markers_adt %>%
  group_by(cluster) %>%
  arrange(cluster, p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n=50) %>% 
  rename(wilcoxon_effect_size=avg_log2FC)  %>% 
  write_csv(sprintf("%s/top50_adt_markers_per_clus_%s.csv", OUT_DIR, marker.res))


###############################
# try to do soupx comparisons #
###############################

# (1) RNA dotplot -- soupx vs non-soupx
# top 5 markers from BOTH -- put next to each other
#

# soupx object
OUT_DIR=sprintf("%s/wnn_soupx/", JDM_DIR)
load(sprintf("%s/wnn_obj.RData", OUT_DIR)) # sobj
load(file=sprintf("%s/meta_w_labels.RData", OUT_DIR)) # meta
sobj@meta.data = meta
load(sprintf("%s/wnn_rna_markers_%s.RData", OUT_DIR, marker.res)) # wnn_markers
wnn_markers_soupx = wnn_markers
Idents(sobj) = "wnn_cluster_label"

# create a DietSeurat version that we can actually look at 
DefaultAssay(sobj) = "RNA"
sobj_soupx = DietSeurat(sobj, assays="RNA")
rm(sobj, meta, wnn_markers)
gc()

# non-soupx
OUT_DIR=sprintf("%s/wnn/", JDM_DIR)
load(sprintf("%s/wnn_obj.RData", OUT_DIR)) # sobj
load(file=sprintf("%s/meta_w_labels.RData", OUT_DIR)) # meta
load(sprintf("%s/wnn_rna_markers_%s.RData", OUT_DIR, marker.res)) # wnn_markers
wnn_markers_no_soupx = wnn_markers
sobj@meta.data = meta
Idents(sobj) = "wnn_cluster_label"

# create a DietSeurat version that we can actually look at 
DefaultAssay(sobj) = "RNA"
sobj_no_soupx = DietSeurat(sobj, assays="RNA")
rm(sobj, meta, wnn_markers)
gc()

# now start the comparison 
OUT_DIR=sprintf("%s/compare_soupx/", JDM_DIR)
#dir.create(OUT_DIR)
sobj_soupx
sobj_no_soupx
head(wnn_markers_no_soupx)
head(wnn_markers_soupx)

# intersection of labels?
top5_soupx = wnn_markers_soupx %>%
  group_by(cluster) %>%
  arrange(cluster) %>%
  slice_max(n=5, order_by=avg_log2FC) 

top5_no_soupx = wnn_markers_no_soupx %>%
  group_by(cluster) %>%
  arrange(cluster) %>%
  slice_max(n=5, order_by=avg_log2FC) 


both_top = union(top5_soupx$gene, top5_no_soupx$gene) # 133, 123 of which match
setdiff(top5_soupx$gene, top5_no_soupx$gene) #  "FCN1"   "MT-CO3" "IFI6"   "BCL11A"
setdiff(top5_no_soupx$gene, top5_soupx$gene) # "MT-ND2" "MX1"    "TTN"    "IL1R2"  "MNDA"   "MEF2C" 

# add cluster labels to the markers -- then attempt to compare?
DotPlot(sobj_soupx, features=both_top, cols="RdYlBu", cluster.idents = T)+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave(sprintf("%s/soupx_dotplot.pdf", OUT_DIR), height=20, width=12)
DotPlot(sobj_no_soupx, features=both_top, cols="RdYlBu", cluster.idents = T)+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave(sprintf("%s/no_soupx_dotplot.pdf", OUT_DIR), height=20, width=12)

# filter
# sobj_soupx = subset(sobj_soupx, cells=(
#   rownames(sobj_soupx@meta.data)[!sobj_soupx@meta.data$wnn_cluster_label %in%
#                                    c("T4_Naive 6", "T4_Naive* 2", "B_Naive*")]
# ))
# 
# sobj_no_soupx = subset(sobj_no_soupx, cells=(
#   rownames(sobj_no_soupx@meta.data)[!sobj_no_soupx@meta.data$wnn_cluster_label %in%
#                                    c("T4_Naive 5", "T4_Naive* 2", "T4_Naive* 3")]
# ))
# 
# DotPlot(sobj_soupx, features=both_top, cols="RdYlBu", cluster.idents = T)+
#   coord_flip()+
#   theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
# ggsave(sprintf("%s/soupx_dotplot2.pdf", OUT_DIR), height=20, width=12)
# DotPlot(sobj_no_soupx, features=both_top, cols="RdYlBu", cluster.idents = T)+
#   coord_flip()+
#   theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
# ggsave(sprintf("%s/no_soupx_dotplot2.pdf", OUT_DIR), height=20, width=12)


# little to no difference

### look into the soupx output itself -- which genes do we expect to see difference on?
# sanity check this -- make sure there is not a bug                                                       
load("/krummellab/data1/erflynn/premier/jdm/data/jdm/soupx/soupx_out_well1.RData", verbose=T)
# --> sc

# the counts should be different!
count_mat_soupx = sobj_soupx@assays$RNA@counts
count_mat_no_soupx = sobj_no_soupx@assays$RNA@counts
all(colnames(count_mat_soupx)==colnames(count_mat_no_soupx))
rand_cols = sample(ncol(count_mat_soupx), 1000)
all(count_mat_soupx[,rand_cols]==count_mat_no_soupx[,rand_cols])
data_mat_soupx = sobj_soupx@assays$RNA@data
data_mat_no_soupx = sobj_no_soupx@assays$RNA@data
all(data_mat_soupx[,rand_cols]==data_mat_no_soupx[,rand_cols])

diff_mat = count_mat_soupx-count_mat_no_soupx
diff_mat_data = data_mat_soupx-data_mat_no_soupx

rs = rowSums(diff_mat)
rs_data = rowSums(diff_mat_data)
cs = colSums(diff_mat)
summary(rs)
most_diff_data = rownames(diff_mat_data)[rs_data < -2000]


most_diff = rownames(diff_mat)[rs < -20000]
sobj_soupx@active.ident = factor(sobj_soupx@active.ident,
                                    levels=sort(as.character(unique(sobj_soupx@active.ident))))
DotPlot(sobj_soupx, features=most_diff, cols="RdYlBu")+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave(sprintf("%s/soupx_dotplot_max_diff.pdf", OUT_DIR), height=24, width=12)

sobj_no_soupx@active.ident = factor(sobj_no_soupx@active.ident,
                                 levels=sort(as.character(unique(sobj_no_soupx@active.ident))))

DotPlot(sobj_no_soupx, features=most_diff, cols="RdYlBu")+
  coord_flip()+  
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))

ggsave(sprintf("%s/no_soupx_dotplot_max_diff.pdf", OUT_DIR), height=24, width=12)

DotPlot(sobj_soupx, features=setdiff(most_diff, c("MALAT1", "IGKC")), cols="RdYlBu",scale=F)+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave(sprintf("%s/soupx_dotplot_max_diff_ns.pdf", OUT_DIR), height=24, width=12)


DotPlot(sobj_no_soupx, features=setdiff(most_diff, c("MALAT1", "IGKC")), cols="RdYlBu", scale=F)+
  coord_flip()+  
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))

ggsave(sprintf("%s/no_soupx_dotplot_max_diff_ns.pdf", OUT_DIR), height=24, width=12)

# max median difference
rand_cols = sample(ncol(count_mat_soupx), 10000)
row_meds = apply(diff_mat[,rand_cols], 1, median)

head(rownames(diff_mat)[order(row_meds)], 15)
# plot raw counts

library(SoupX)
out = adjustCounts(sc)
plotChangeMap(sc, out, "IGKC")

