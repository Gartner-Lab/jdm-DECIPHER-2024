
library(Seurat)
library(tidyverse)
library(dittoSeq)
library(cowplot)

DATA_DIR = "/krummellab/data1/erflynn/premier/jdm/data/"
setwd("/krummellab/data1/erflynn/premier/jdm/data")
prefix='merged_harmony_w_wnn/merged_'

load(file=paste0(prefix, '_merged_w_clus.RData'))




# --- add in metadata from JDM, cSLE --- #
jdm_obj_initial = readRDS(sprintf("%s/jdm/analyzed_seurat_object (1).rds", DATA_DIR))
jdm_meta = jdm_obj_initial@meta.data %>% 
  as_tibble(rownames="cell_id")
rm(jdm_obj_initial)



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


pt_data = jdm_meta %>% 
  distinct(study_id_visit, donor, study_id, disease_group, on_meds, sex, age, group, visit, group_pool)
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


jdm_cl_to_study_id = meta.data %>%
  as_tibble(rownames="cell_id") %>%
  mutate(BEST.GUESS=str_extract(as.character(BEST.GUESS), "^[0-9]+"),
         well=as.character(well)) %>%
  left_join(jdm_meta %>% 
              distinct(BEST.GUESS, well, study_id_visit, donor, case_control, group, visit, age, sex, 
                       disease_group, on_meds),
              mutate(BEST.GUESS=as.character(BEST.GUESS),
                     well=as.character(well)), by=c("BEST.GUESS", "well")) %>%
  left_join(jdm_meta %>% dplyr::select(cell_id, manual_labels),
            by="cell_id")

# clean up the metadata
meta.data2 = jdm_cl_to_study_id %>% 
  mutate(orig.ident = ifelse(orig.ident=="SeuratProject", "jdm", orig.ident)) %>%
  dplyr::select(cell_id, orig.ident, LIBRARY, nCount_RNA, nFeature_RNA,
                predicted.celltype.l2, Groups, Batch, Age, Gender, Race,
                harmony_res0.8, donor:study_id_visit, case_control, manual_labels) %>%
  dplyr::rename(csle_labels=predicted.celltype.l2, jdm_labels=manual_labels) %>%
  mutate(sex=as.character(sex)) %>%
  mutate(age = ifelse(is.na(age), age, Age),
         sex= ifelse(is.na(Gender), sex, Gender)) %>%
  dplyr::select(-Gender, -Age) %>%
  mutate(case_control=case_when(
    case_control == "JDM" ~ "JDM",
    case_control == "HC" ~ "JDM - HC",
    Groups =="cHD" ~ "cSLE - HC",
    Groups=="cSLE" ~ "cSLE"
  )) 
rm(jdm_meta)
rm(jdm_cl_to_study_id)
rm(meta.data)

# add the metadata
load(file=paste0(prefix, '_merged_processed0.RData'))
meta.data2 %>% mutate()
meta.data3 = data.frame(meta.data2 %>% dplyr::select(-cell_id))
rownames(meta.data3) = meta.data2$cell_id
merged_data@meta.data = meta.data3


# plot with the different groups colored?
DimPlot(merged_data, split.by='case_control', cols=dittoColors())
p1 = DimPlot(merged_data, cells.highlight=meta.data2 %>% filter(case_control=="JDM - HC") %>% pull(cell_id))+
  NoLegend()+ggtitle("JDM - HC")

p2 = DimPlot(merged_data, cells.highlight=meta.data2 %>% filter(case_control=="cSLE - HC") %>% pull(cell_id))+
  NoLegend()+ggtitle("cSLE - HC")

p3 = DimPlot(merged_data, cells.highlight=meta.data2 %>% filter(case_control=="JDM") %>% pull(cell_id))+
  NoLegend()+ggtitle("JDM")

p4 = DimPlot(merged_data, cells.highlight=meta.data2 %>% filter(case_control=="cSLE") %>% pull(cell_id))+
  NoLegend()+ggtitle("cSLE")
cowplot::plot_grid(p1, p2, p3, p4, ncol=4)
ggsave(sprintf("%s_case_control.png", prefix), height=5, width=16)

Idents(merged_data) = "harmony_res0.8"
DimPlot(merged_data, label=T, cols=dittoColors())+NoLegend()
ggsave(sprintf("%s_clusters0.8.png", prefix), height=5, width=5)

DimPlot(merged_data, group.by="csle_labels", 
        split.by="orig.ident", label=T, cols=dittoColors(), ncol=3)+NoLegend()
ggsave(sprintf("%s_csle_labels.png", prefix), height=6, width=18)
DimPlot(merged_data, group.by="jdm_labels",
        split.by="orig.ident", label=T, cols=dittoColors(), ncol=3)+NoLegend()
ggsave(sprintf("%s_jdm_labels.png", prefix), height=6, width=18)


### look at the meta



#load(sprintf("%s/sobj_GSE135779/merged_df_proc2/GSE135779_merged_processed.RData", DATA_DIR)) # --> sobj
#csle_meta = sobj@meta.data %>% 
#  as_tibble(rownames="cell_id")
#rm(sobj)

clus_lab = meta.data2 %>% 
  dplyr::select(cell_id, orig.ident, harmony_res0.8, jdm_labels, csle_labels)

clus_to_lab = clus_lab %>%
  dplyr::rename(cluster="harmony_res0.8") %>%
  group_by(cluster ) %>%
  mutate(tot=n()) %>%
  group_by(cluster, orig.ident, jdm_labels, csle_labels, tot) %>%
  dplyr::count() %>%
  arrange(cluster, desc(n)) %>%
  mutate(frac=n/tot)
length(unique(clus_to_lab$cluster)) # 35

clus_to_lab %>%
  ungroup() %>%
  group_by(cluster, orig.ident) %>%
  filter(frac > 0.4) %>%
  View()
clus_lab_by_ident = clus_to_lab %>%
  ungroup() %>%
  group_by(cluster, orig.ident) %>%
  slice_max(1) %>%
  ungroup() %>%
  mutate(jdm_labels=as.character(jdm_labels),
         csle_labels=as.character(csle_labels)) %>%
  mutate(prev_label=ifelse(is.na(csle_labels), jdm_labels, csle_labels)) %>%
  dplyr::select(-csle_labels, -jdm_labels) %>%
  dplyr::select(cluster, orig.ident, prev_label, everything()) %>%
  pivot_wider(id_cols=c("cluster", "tot"), names_from="orig.ident",
              values_from=c("prev_label", "n", "frac")) %>%
  dplyr::select(everything(), tot) 

clus_lab_by_ident %>% write_csv(sprintf("%s_clus_lab_ident.csv", prefix))

clus_one_lab = clus_to_lab %>%
  ungroup() %>%
  group_by(cluster) %>%
  slice_max(1) %>%
  ungroup() %>%
  mutate(jdm_labels=as.character(jdm_labels),
         csle_labels=as.character(csle_labels)) %>%
  mutate(prev_label=ifelse(is.na(csle_labels), jdm_labels, csle_labels))

clus_one_lab2 = clus_one_lab %>% distinct(cluster, prev_label) %>%
  arrange(prev_label, cluster) %>%
  group_by(prev_label) %>%
  mutate(idx=1:n(),
         nlab=n()) %>%
  mutate(cluster_label=ifelse(nlab ==1, prev_label,
                              paste(prev_label, idx))) %>%
  ungroup()

meta.data4 = meta.data3 %>% left_join(clus_one_lab2 %>% 
                           distinct(cluster_label, cluster), 
                         by=c("harmony_res0.8"="cluster") )
# now plot with this?
meta.data4 %>%
  ggplot(aes(x=cluster_label, fill=orig.ident)) +
  geom_bar(position="fill")+
  theme_bw()+
  scale_fill_manual(values=dittoColors())+
  ylab("fraction")+
  xlab("max cluster label")+
  labs(fill="dataset")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave(sprintf("%s_fraction_dataset.png", prefix))  


meta.data4 %>%  
  ggplot(aes(x=cluster_label, fill=case_control)) +
  geom_bar(position="fill")+
  theme_bw()+
  scale_fill_manual(values=dittoColors())+
  ylab("fraction")+
  xlab("max cluster label")+
  labs(fill="group")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
ggsave(sprintf("%s_fraction_groups.png", prefix)) 



load(file=paste0(prefix, '_merged_processed0.RData'))
merged_data@meta.data = meta.data4
#Idents(merged_data) = "cluster_label"
unique(merged_data@meta.data$cluster_label)
DimPlot(merged_data, group.by="cluster_label", label=T, cols=dittoColors())+NoLegend()
ggsave(sprintf("%s_clusters_lab.png", prefix), height=5, width=5)


### STOP

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



# cluster labels

# cluster fractions

# healthy control locations