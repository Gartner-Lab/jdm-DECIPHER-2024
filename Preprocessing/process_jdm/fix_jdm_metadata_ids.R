library(Seurat)
library(tidyverse)

## code for fixing the metadata ids
jdm_obj_initial = readRDS(sprintf("%s/analyzed_seurat_object.rds", JDM_DIR))
jdm_meta0 = jdm_obj_initial@meta.data %>% 
  as_tibble(rownames="cell_id")
rm(jdm_obj_initial)  

jdm_obj = readRDS(sprintf("%s/seurat_filt_postdemuxDF.rds", JDM_DIR))

# add percent.ribo b/c it's not there
jdm_obj = jdm_obj %>% 
  PercentageFeatureSet(pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", 
                       col.name = "percent.ribo")


jdm_meta = jdm_obj@meta.data  %>%
  mutate(best.guess.well = paste(BEST.GUESS, well, sep = "_")) %>%
  mutate(study_id_visit=case_when(
    best.guess.well %in% c("7,7_1", "7,7_2") ~ "A1-019_V1",
    best.guess.well %in% c("7,7_3", "7,7_4") ~ "A1-019_V7",
    best.guess.well %in% c("6,6_1", "6,6_2") ~ "A1-009_V6",
    best.guess.well %in% c("6,6_5", "6,6_6") ~ "A1-009_V8",
    best.guess.well %in% c("9,9_1", "9,9_2") ~ "A2-018_V1",
    best.guess.well %in% c("9,9_5",  "9,9_6") ~ "A2-018_V2", 
    best.guess.well %in% c("11,11_1", "11,11_2") ~ "A4-027_V1",
    best.guess.well %in% c("11,11_3", "11,11_4") ~ "A4-027_V4",
    best.guess.well %in% c("5,5_1", "5,5_2") ~ "A1-026_V1",
    best.guess.well %in% c("5,5_3", "5,5_4") ~ "A1-026_V5",
    best.guess.well %in% c("5,5_5",  "5,5_6") ~ "A1-026_V3",
    best.guess.well %in% c("3,3_3", "3,3_4") ~ "A1-031_V1",
    best.guess.well %in% c("13,13_5", "13,13_6") ~ "A1-034_V1", 
    best.guess.well %in% c("13,13_3","13,13_4") ~ "A1-034_V2", 
    BEST.GUESS=="19,19" ~ "A1-030_V1", 
    BEST.GUESS=="10,10" ~ "A1-025_V1", 
    BEST.GUESS=="18,18" ~ "A2-028_V1", 
    BEST.GUESS=="14,14" ~ "A2-005_V3", 
    BEST.GUESS=="8,8"   ~ "A2-011_V3", 
    BEST.GUESS=="17,17" ~ "A5-006", 
    BEST.GUESS=="16,16" ~ "A5-021", 
    BEST.GUESS=="15,15" ~ "A5-022", 
    BEST.GUESS=="12,12" ~ "A5-024", 
    BEST.GUESS=="0,0"   ~ "A5-039", 
    BEST.GUESS=="1,1"   ~ "A1-037_V1", 
    BEST.GUESS=="2,2"   ~ "A1-036_V1", 
    TRUE ~"A1-032_V1")
    )


# fix var coding
tn_jdm <- c("A1-019_V1", "A1-030_V1", "A1-031_V1", "A1-026_V1", "A1-034_V1", 
            "A1-036_V1", "A1-037_V1", "A1-032_V1", "A1-025_V1")
inact_off_meds <- c("A1-009_V6", "A2-005_V3", "A1-019_V7", "A2-018_V1", 
                    "A2-028_V1", "A2-011_V3")
inact_on_meds <- c("A4-027_V4", "A1-026_V3")
flare <- c("A1-009_V8", "A2-018_V2", "A4-027_V1", "A1-026_V5", "A1-034_V2")

jdm_meta = jdm_meta %>%
  mutate(disease_group = case_when(
    study_id_visit %in% tn_jdm ~ "TNJDM", 
    study_id_visit %in% inact_off_meds ~ "InactOffMeds",
    study_id_visit %in% inact_on_meds ~ "InactOnMeds",
    study_id_visit %in% flare ~ "Flare", 
    TRUE ~ "HC")
    ) %>%
  mutate(case_control= case_when(
    disease_group=="HC" ~ "HC",
    TRUE ~ "JDM"))
  
  
table(jdm_meta$case_control, jdm_meta$study_id_visit)

jdm_meta = jdm_meta %>%
  select(-best.guess.well) %>%
  separate(study_id_visit, 
           into=c("donor", "visit"), 
           sep = "_", remove=F) %>%
  left_join(jdm_meta0 %>%  # add age/sex info
              distinct(donor, age, sex) %>%
              filter(donor!="nan"))

# note: 5 donors are missing age/sex info

# add the cell_ids before replacing the metadata
rownames(jdm_meta) = rownames(jdm_obj@meta.data) 
jdm_obj@meta.data = jdm_meta

saveRDS(jdm_obj, sprintf("%s/seurat_filt_postdemuxDF_w_meta.rds", JDM_DIR),
        compress = F)

