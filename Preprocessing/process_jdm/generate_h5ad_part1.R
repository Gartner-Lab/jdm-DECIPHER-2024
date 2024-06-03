# this version *does not* concatenate the RNA and ADT
# instead it writes out just the RNA and saves the ADT as a matrix

library(Seurat)
library(tidyverse)
library(SeuratDisk)
library(readxl)
library(dittoSeq)

setwd("/krummellab/data1/erflynn/premier/jdm/data/jdm/h5ad")
# three parts: (1) metadata, (2) rna, (3) adt

dat = Read10X("../raw/filtered_feature_bc_matrix/")
raw_rna = dat$`Gene Expression`
raw_adt = dat$`Antibody Capture`
save(raw_rna, raw_adt, file="raw_mats.RData")
load("raw_mats.RData")

load("../wnn/reclus/sobj_w_raw_adt.RData", verbose=T)

## load the new labels, disease activity from Camilla
labels = read_tsv("../annot/labels_07-22-23.tsv")
patient_meta2 = read_excel("../annot/metadata_updated_08-09-23.xlsx") # use this for age
re_info = read_csv("../annot/JDMcross_long_raceethnicity_08-11-23.csv")

# check that we have the exact same set of cells
stopifnot(length(labels$barcodes)==length(intersect(labels$barcodes, colnames(sobj))))

stopifnot(length(intersect(unique(patient_meta2$donor), 
                           unique(re_info$`Study ID`)))==min(nrow(re_info), nrow(patient_meta2)))

# subset the object to remove final set of cells
sobj = subset(sobj, cells=labels$barcodes)


## as much metadata as possible
cols_to_keep =  c("orig.ident" , "nCount_RNA" , "nFeature_RNA", "nCount_ADT", "nFeature_ADT",
                  "well",  "percent.mt", "percent.ribo", "S.Score", "G2M.Score", "Phase",
                  "study_id_visit", "donor","visit" , "case_control",
                  "sex", "RNA.weight", "wsnn_leiden_reclus_res.1.4")

meta = sobj@meta.data %>%
  as_tibble(rownames="cell_id") %>%
  select(cell_id, cols_to_keep) %>%
  left_join(patient_meta2 %>% select(-well, -donor, -case_control, -sex, -visit, 
                                     -medritux, -medcyc, -p155, -nxp2, -mda5), 
            by=c("study_id_visit")) %>%
  rename(disease_group=new_diseasegroups,
         age=agey,
         cluster_number=wsnn_leiden_reclus_res.1.4) %>%
  left_join(re_info %>% select(`Study ID`, Race), by=c("donor"="Study ID")) %>%
 # left_join(dz_activity %>% rename("cell_id"="...1"), 
#            by="cell_id") %>%
  left_join(labels %>% select(barcodes, new_labels_long) %>%
              rename(cluster_label=new_labels_long), by=c("cell_id"="barcodes")) %>%
  mutate(orig.ident=well) %>%
  select(-well) %>%
  select(orig.ident:Phase, RNA.weight, cluster_number, cluster_label, study_id_visit:visit, 
         disease_group, case_control:age, Race, on_meds, steroids,  
         other_immune_suppressant, MSA, DA_cat, everything() ) %>%
  column_to_rownames("cell_id") 
stopifnot(nrow(meta)==nrow(sobj@meta.data))



# check we selected correct cluster_number -- is this right?
meta %>% distinct(cluster_number, cluster_label) %>%  nrow()
meta %>% distinct(cluster_label) %>% nrow()

#title: CITEseq of JDM PBMCs
#description: Mapping immune dysregulation in juvenile dermatomyositis at single cell resolution
#For the four VAS scores, these can be continuous variables; same with CDASIact scores;
#"menzymeyn" is binary (y/n), "mmt8score" continuous, "chaqscore" continuous; remove "p155", "nxp2", 
#"mda5" as this is all captured in "MSA" which is categorical;  "DA_cat" keep categorical
#For the IDs, everything looks right, but I am wondering if it would make sense to indicate the HC patients in the ID somehow? 
#like HC instead of JDM1 at the start? 
#For the metadata file, many of the numbers have multiple decimal places--would it make sense to round to two decimals? 
# remove donor


sample_id_mapping = meta %>% 
  distinct(donor, case_control) %>%
  arrange(donor) %>% 
  mutate(ind_idx = 1:n()) %>%
  mutate(donor_id=sprintf("%s%s", case_control, ind_idx),
         colabs_id=sprintf("JDM1-HS%s", ind_idx)) %>%
  ungroup() %>%
  left_join(meta %>% distinct(donor, visit)) %>%
  mutate(visit=ifelse(is.na(visit), "V1", visit)) %>% 
  mutate(visit_id=sprintf("%s-E%sPB1", colabs_id, str_extract(visit, "[0-9]+")))
sample_id_mapping %>% write_csv("sample_id_mapping.csv")


meta2 = meta %>%
  as_tibble(rownames="cell_id") %>%
  mutate(visit=ifelse(is.na(visit), "V1", visit)) %>% 
  left_join(sample_id_mapping %>% distinct(donor, donor_id), by="donor") %>%
  select(-donor, -study_id_visit) %>%
  mutate(cluster_label=case_when(
    cluster_label=="B_cluster5" ~ "B mem",
    str_detect(cluster_label, "^gdT") ~ str_replace_all(cluster_label, "_cluster", " "),
    str_detect(cluster_label, "^B") ~ str_replace_all(cluster_label, "_cluster", " naive "), 
    cluster_label=="CD8+ memory_resting" ~ "CD8+ mem",
    TRUE ~ cluster_label
  )) %>% # fix the cell type labels
  column_to_rownames("cell_id") 
stopifnot(nrow(meta2)==nrow(meta))
meta2 %>% as_tibble(rownames="cell_id") %>% write_csv("metadata_09-18-23.csv")
meta2 = read_csv("metadata_09-18-23.csv") %>%
  column_to_rownames("cell_id")

# required by CZI:
# options: under 2 (HsapDv:0000083), 2-5 (HsapDv:0000084), 6-12 (HsapDv:0000085), 13-18 (HsapDv:0000086)
# -or- by year ontology

meta_w_czi_lab = meta2 %>% 
  mutate(
    organism_ontology_term_id="NCBITaxon:9606",
    assay_ontology_term="EFO:0009294", # CITE-seq
    tissue_ontology_term_id="CL:2000001", # UBERON PBMC
    suspension_type="cell",
    developmental_stage_ontology_term_id= case_when(
      age < 2 ~ "HsapDv:0000083",
      age <= 5 ~ "HsapDv:0000084",
      age <= 12 ~ "HsapDv:0000085",
      age <= 18 ~ "HsapDv:0000086"
    ),
    sex_ontology_term_id=case_when(sex=="M" ~ "PATO:0000384", 
                                   sex=="F" ~ "PATO:0000383"),
    disease_ontology_term_id=case_when(case_control=="HC" ~ "PATO:0000461", # normal
                                         case_control=="JDM" ~ "MONDO:0008054"), # JDM
    self_reported_ethnicity_ontology_term_id=case_when(
      Race == "White" ~ "HANCESTRO:0005",
      Race == "Latino, or Spanish origin" ~ "HANCESTRO:0014",
      Race == "Asian" ~ "HANCESTRO:0008",
      Race == "Middle Eastern/North African"  ~ "HANCESTRO:0015",
      Race != "Other" ~ "multiethnic"
    ),
    cell_type_ontology_term_id=case_when(
      str_detect(cluster_label, "^B naive") ~ "CL:0000788", # naive B cell
      cluster_label == "B mem" ~ "CL:0000787",
      str_detect(cluster_label, "^gdT") ~"CL:0000798", 
      cluster_label=="CD14+ monocytes" ~ "CL:0002057", # CD14-positive, CD16-negative classical monocyte
      cluster_label=="CD4+ Tregs" ~ "CL:0000792", # CD4-positive, CD25-positive, alpha-beta regulatory T cell 
      cluster_label == "CD8+ naive T" ~ "CL:0000900 ", # naive thymus-derived CD8-positive, alpha-beta T cell
      cluster_label== "CD56dim NK"  ~ "CL:0000939" , # CD16-negative, CD56-dim natural killer cell, human
      cluster_label=="CD4+ effector T" ~ "CL:0001044", 
      cluster_label=="CD4+ naive T" ~ "CL:0000895", # naive thymus-derived CD4-positive, alpha-beta T cell
      cluster_label=="CD56bright NK" ~ "CL:0000938", # CD16-negative, CD56-bright natural killer cell, human
      cluster_label=="Classic dendritic" ~ "CL:0000990", # Conventional dendritic cell 
      cluster_label=="CD8+ GZMA/Bhi" ~ "CL:0001050", # effector CD8-positive, alpha-beta T cell
      cluster_label == "PDCs" ~ "CL:0000784",
      cluster_label=="CD16+ monocytes" ~ "CL:0002396", # CD14-low, CD16-positive monocyte
      cluster_label=="CD8+ GZMKhi"  ~ "CL:0001050" ,   # effector CD8-positive, alpha-beta T cell    
      cluster_label=="CD8+ mem" ~ "CL:0000909",
      cluster_label=="Plasmablasts" ~ "CL:0000980"
    )
  ) 
meta_w_czi_lab2 = meta_w_czi_lab %>%
  select(-sex, -Race) 

# expand so that the norm data has *ALL* genes, just with zeros for the missing genes
removed_genes = setdiff(rownames(raw_rna), rownames(sobj@assays$RNA@data))

tibble("gene"=removed_genes) %>% write_csv("filtered_genes.csv")


rna_zero_mat = matrix(0, ncol=ncol(sobj@assays$RNA@data), nrow=length(removed_genes))
rownames(rna_zero_mat) = removed_genes
colnames(rna_zero_mat) = colnames(sobj@assays$RNA@data)
norm_rna = rbind(sobj@assays$RNA@data, rna_zero_mat)
norm_rna2 = norm_rna[rownames(raw_rna),]



sobj_new = CreateSeuratObject(counts=raw_rna[, rownames(meta_w_czi_lab2)], meta.data=meta_w_czi_lab2)
sobj_new@assays$RNA@data = norm_rna2



Idents(sobj_new) = "cluster_label"
sobj_new@reductions[["wnn.umap"]] = CreateDimReducObject(sobj@reductions$wnn.umap@cell.embeddings,
                                                    key="wnnUMAP", assay="RNA")
sobj_new@reductions[["rna.umap"]] = CreateDimReducObject(sobj@reductions$umap@cell.embeddings,
                                                    key="rnaUMAP", assay="RNA")
sobj_new@reductions[["adt.umap"]] = CreateDimReducObject(sobj@reductions$a.umap@cell.embeddings,
                                                    key="adtUMAP", assay="RNA")
DimPlot(sobj_new, label=T, reduction="wnn.umap", repel=T, size=3,
        cols=dittoColors())

sobj_new@meta.data$orig.ident = sapply(sobj_new@meta.data$orig.ident, function(x) sprintf("library%s", x))
save(sobj_new, file="sobj_before_conversion_09-20-23.RData")


SaveH5Seurat(sobj_new, filename = "jdm_obj_09-20-23.h5Seurat")
Convert("jdm_obj_09-20-23.h5Seurat", dest = "h5ad")

### ADT
# fix the rownames for the raw_adt vs integrated.ADT
norm_adt_feats = rownames(sobj@assays$integrated.ADT@data)
raw_adt_feats = str_replace_all(str_replace_all(rownames(raw_adt), "\\.1", ""), "_", "-")
rownames(raw_adt) = raw_adt_feats
stopifnot(length(setdiff(norm_adt_feats, raw_adt_feats))==0)

removed_adts = setdiff(raw_adt_feats, norm_adt_feats)
tibble("adt"=removed_adts) %>% write_csv("filtered_adts.csv")
adt_zero_mat = matrix(0, ncol=ncol(sobj@assays$integrated.ADT@data), nrow=length(removed_adts))
rownames(adt_zero_mat) = removed_adts
colnames(adt_zero_mat) = colnames(sobj@assays$integrated.ADT@data)
norm_adt = rbind(sobj@assays$integrated.ADT@data,adt_zero_mat)
norm_adt2 = norm_adt[rownames(raw_adt),]
sobj_adt = CreateSeuratObject(counts=raw_adt[, rownames(meta_w_czi_lab2)], meta.data=meta_w_czi_lab2)
sobj_adt@assays$RNA@data = norm_adt2

SaveH5Seurat(sobj_adt, filename = "jdm_adt_09-20-23.h5Seurat")
Convert("jdm_adt_09-20-23.h5Seurat", dest = "h5ad")

# format the feature info
feat_info = read_excel("../annot/human_TSA_99786_B307739_Barcode List.xlsx")
stopifnot(nrow(feat_info)==length(raw_adt_feats))
feat_info$feature_id = raw_adt_feats
feat_info %>% View() # check to make sure these line up!
feat_info %>% write_csv("adt_feature_df.csv")
