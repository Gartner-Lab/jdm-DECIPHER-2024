
library(Seurat)
library(tidyverse)
library(dittoSeq)
library(cowplot)
DATA_DIR = "/krummellab/data1/erflynn/premier/jdm/data/"

source("/krummellab/data1/erflynn/scratch_scripts/FuncPctDotPlot.R")


# --- add in metadata from JDM, cSLE --- #
jdm_obj_initial = readRDS(sprintf("%s/jdm/analyzed_seurat_object (1).rds", DATA_DIR))
jdm_meta = jdm_obj_initial@meta.data %>% 
  as_tibble(rownames="cell_id")
rm(jdm_obj_initial)

jdm_obj2 = readRDS(sprintf("%s/jdm/seurat_filt_postdemuxDF.rds", DATA_DIR))
jdm_meta2 = jdm_obj2@meta.data %>%
  as_tibble(rownames="cell_id")
rm(jdm_obj2)

# fix var coding
tn_jdm <- c("A1-019_V1", "A1-030_V1", "A1-031_V1", "A1-026_V1", "A1-034_V1", "A1-036_V1", "A1-037_V1", "A1-032_V1", "A1-025_V1")
inact_off_meds <- c("A1-009_V6", "A2-005_V3", "A1-019_V7", "A2-018_V1", "A2-028_V1", "A2-011_V3")
inact_on_meds <- c("A4-027_V4", "A1-026_V3")
flare <- c("A1-009_V8", "A2-018_V2", "A4-027_V1", "A1-026_V5", "A1-034_V2")

jdm_meta$disease_group <- ifelse(jdm_meta$study_id_visit %in% tn_jdm, "TNJDM", ifelse(
  jdm_meta$study_id_visit %in% inact_off_meds, "InactOffMeds", ifelse(
    jdm_meta$study_id_visit %in% inact_on_meds, "InactOnMeds", ifelse(
      jdm_meta$study_id_visit %in% flare, "Flare", "HC"
    )
  )
))

# fix the na/nan study_id/donor

lapply(jdm_meta2 %>% 
         distinct(well, BEST.GUESS) %>% 
         arrange(well, BEST.GUESS) %>% group_split(well), function(x) x %>%  pull(BEST.GUESS))


lapply(jdm_meta2 %>% 
         distinct(well, BEST.GUESS) %>% 
         arrange(well, BEST.GUESS) %>% group_split(well), function(x) nrow(x))



# how many per pool-well?
# some wells are missing certain samples:
#  well 1: sample 3,13
#  well 2: sample 3,13
#  well 3: sample 6,9
# well 4: sample 6,9
# well 5: sample 3,7,11
# well 6: sample 3,7,11

pt_data = jdm_meta %>% 
  distinct(study_id_visit, donor, study_id, disease_group,
           on_meds, sex, age, group, visit, group_pool, well, BEST.GUESS)
pt_data0 = pt_data %>%
  mutate(study_id=as.character(study_id)) %>%
  mutate(study_id=case_when(
    study_id != "na" ~ study_id,
    TRUE ~ str_extract(study_id_visit, "A[0-9]+-[0-9]+") 
  )) %>%
  mutate(donor= str_extract(study_id, "A[0-9]+-[0-9]+") ) 

pt_data0 %>% filter(is.nan(age)) %>%
  dplyr::select(-sex, -age) %>%
  left_join(pt_data0 %>% distinct(donor, sex, age) %>%
              filter(!is.nan(age)), by="donor")

pt_data0 %>% 
  distinct(donor, disease_group, on_meds, sex, age) %>%
  arrange(donor) %>% View()
# "A5-039" # we don't have info on
# "A1-032", "034", "036", "037" also no info

# what is the difference btw study_id, study_id_visit, and donor? some study_id have visit numbers

pt_data0 %>% 
  distinct(well, BEST.GUESS, study_id) %>%
  arrange(well, BEST.GUESS) %>% View()


### load the data and add the info
JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm"
OUT_DIR = sprintf("%s/wnn_out_rpca_DSB3/", JDM_DIR)
setwd(OUT_DIR)
load("sobj_filt.RData")

# add the pt metadata
meta = sobj_filt@meta.data %>%
  as_tibble(rownames="cell_id") %>%
  mutate(BEST.GUESS=str_extract(as.character(BEST.GUESS), "^[0-9]+"),
         well=as.character(well)) 
meta2 = meta %>%
  select(-donor, -study_id_visit) %>%
  left_join(pt_data0, by=c("well", "BEST.GUESS")) %>%
  select(-cell_id)

sobj_filt@meta.data = data.frame(meta2)
rownames(sobj_filt@meta.data ) = meta$cell_id
length(unique(sobj_filt@meta.data$seurat_clusters)) # 28
all((sobj_filt@meta.data$seurat_clusters)==(sobj_filt@meta.data$wsnn_res.1.4)) # 28

clus_lab = sobj_filt@meta.data  %>% 
  as_tibble(rownames="cell_id") %>%
  dplyr::select(cell_id, wsnn_res.1.4) %>%
  left_join(  
    jdm_meta %>% 
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
length(unique(clus_to_lab$cluster)) # 

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

new_data = sobj_filt@meta.data %>% 
  as_tibble(rownames="cell_id") %>%
  select(cell_id, wsnn_res.1.4) %>%
  left_join(clus_to_lab3 %>% select(cluster, cluster_label) %>%
              mutate(cluster_label=as.factor(cluster_label)), by=c("wsnn_res.1.4"="cluster"))

manual_lab_df = clus_lab %>% 
  mutate(manual_labels=as.character(manual_labels),
         manual_labels=ifelse(is.na(manual_labels), "unknown",
                              manual_labels)) 

sobj_filt = AddMetaData(sobj_filt, new_data$cluster_label, col.name="cluster_label")
sobj_filt = AddMetaData(sobj_filt, manual_lab_df$manual_labels, col.name="manual_labels")


col_vec = dittoColors()[1:length(unique(sobj_filt@meta.data$cluster_label))]
names(col_vec)= unique(sobj_filt@meta.data$cluster_label)
dir.create("reclus")
setwd("reclus")


Idents(sobj_filt) = "cluster_label"
## remake the plots
DimPlot(sobj_filt, reduction = 'wnn.umap', 
        label = TRUE, cols=col_vec,
        label.size=3) + NoLegend()
ggsave("wnn_umap_clus_lab.pdf", height=5, width=5)


DimPlot(sobj_filt, reduction = 'wnn.umap', group.by="manual_labels",
        label = TRUE, cols=dittoColors(),
        label.size=3) 
ggsave("wnn_umap_manual_lab.pdf", height=5,
       width=9)

VlnPlot(sobj_filt, features = "RNA.weight", cols=col_vec, 
        group.by = 'cluster_label', sort = TRUE, pt.size = 0) +
  NoLegend()+
  theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
        axis.title.x = element_blank())
ggsave("wnn_vln_clus_lab.png", height=5, width=8)




## redo the plots
DimPlot(sobj_filt, reduction="wnn.umap", group.by="well")
ggsave("wnn_umap_well.pdf", height=5, width=5)

DimPlot(sobj_filt, reduction="wnn.umap", group.by="donor")
ggsave("wnn_umap_donor.pdf", height=5, width=5)

sobj_filt@meta.data %>%
  dplyr::select(study_id_visit, cluster_label) %>%
  ggplot(aes(x=cluster_label, fill=study_id_visit))+
  geom_bar(position="fill")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))+
  ylab("fraction of cells")
ggsave("wnn_clus_study_id_visit.pdf")


sobj_filt@meta.data %>%
  dplyr::select(well, cluster_label) %>%
  ggplot(aes(x=cluster_label, fill=well))+
  geom_bar(position="fill")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))+
  ylab("fraction of cells")
ggsave("wnn_clus_well.pdf")


## ALLUVIAL
# alluvial for data changes
library(ggalluvial)

meta_alluv =sobj_filt@meta.data  %>%
  as_tibble(rownames="cell_id")  %>%
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

sobj = sobj_filt
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





# dotplots
load("../wnn_rna_filt_1.4.RData")
load("../wnn_adt_filt_1.4.RData")


top5_adt_genes = wnn_markers_adt %>% 
  left_join(clus_map) %>%
  group_by(cluster_label) %>%
  slice_max(n=5, order_by = avg_log2FC) %>%
  arrange(cluster_label)
top5_adt_genes$gene = factor(top5_adt_genes$gene, levels=unique(top5_adt_genes$gene))

FuncPctDotPlot(sobj, assay="integrated.ADT", features = unique(top5_adt_genes$gene),
               cols="RdYlBu",
               cluster.idents = T)+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("adt_dotplot.pdf", width=14, height=17)

top5_genes = wnn_markers %>% 
  left_join(clus_map) %>%
  group_by(cluster_label) %>%
  slice_max(n=5, order_by = avg_log2FC) %>%
  pull(gene)
DotPlot(sobj, assay="RNA", features=unique(top5_genes),
        cols="RdYlBu", cluster.idents = T)+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave("rna_dotplot.pdf", width=11, height=20)

