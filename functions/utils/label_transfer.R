
# map azimuth label to Jess's labels
map_azimuth_label = function(ds, azimuth_label, new_label){
  ds %>%
    mutate({{new_label}}:=case_when(
      {{azimuth_label}} == "NK_CD56bright" ~ "NK_CD56++", 
      {{azimuth_label}} == "CD14 Mono" ~ "cM",
      {{azimuth_label}} == "CD16 Mono" ~ "nCM",
      {{azimuth_label}} == "B memory" ~ "B_Mem",
      {{azimuth_label}} == "B naive" ~ "B_Naive",
      {{azimuth_label}} == "dnT" ~ "DN_T+",
      {{azimuth_label}} =="gdT" ~ "Tgd",
      {{azimuth_label}} =="ILC" ~ "ILC1_2_3",
      {{azimuth_label}} == "CD4 Naive" ~ "T4_Naive",
      {{azimuth_label}} =="CD8 Naive" ~ "T8_Naive",
      {{azimuth_label}} =="MAIT" ~ "T8_MAIT",
      {{azimuth_label}} == "CD8 TEM" ~ "T8_TEMRA",
      {{azimuth_label}} == "CD4 TCM" ~ "T4_Mem",
      {{azimuth_label}} == "CD4 CTL" ~ "T_Tox",
      {{azimuth_label}} == "CD4 TEM" ~ "T4_Mem",
      {{azimuth_label}} == "CD8 TCM" ~  "T8_Mem",
      #    azimuth_label == "CD8 Proliferating" ~,
      #    azimuth_label == "CD4 Proliferating" ~,
      {{azimuth_label}} == "NK Proliferating" ~ "NKs_mix_prolif",
      {{azimuth_label}} == "Plasmablast" ~ "PB",
      {{azimuth_label}} == "HSPC" ~ "HSC" ,
      {{azimuth_label}} == "B intermediate" ~ "B_Mem_intermediate",
      TRUE ~ {{azimuth_label}}
    ))
}
 

# wrapper around Seurat AddMetaData to add metadata through a join
addMetaJoin = function(sobj, new_data, column_name, overwrite=F){
  meta = sobj@meta.data %>% 
    as_tibble(rownames="cell_id") 
  
  if (column_name %in% colnames(meta)){
    if (!overwrite){
      print(sprintf("Column %s already present in the metadata. Switch overwrite to TRUE to add anyway.", column_name))
      return(sobj)
    }
    print(sprintf("Column %s already present in the metadata. Overwriting.", column_name))
    meta = meta %>%
      select(-{{column_name}})
  }
  sobj_annot = meta %>%
    left_join(new_data, by=c("cell_id")) %>%
    select(cell_id, {{column_name}})


  
  sobj = AddMetaData(sobj, sobj_annot %>% pull({{column_name}}), column_name)
  return(sobj)
}

reform_meta = function(meta, my_cluster){
 meta  %>% 
    as_tibble(rownames="cell_id") %>%
    select(cell_id, {{my_cluster}}, manual_labels, azimuth_label_mapped) %>%
    rename(cluster={{my_cluster}}) %>%
    mutate(manual_labels=ifelse(is.na(manual_labels), azimuth_label_mapped, manual_labels),
           label_src=ifelse(is.na(manual_labels), "azimuth", "manual")) %>%
    select(-azimuth_label_mapped)
}

reform_meta0 = function(meta, my_cluster){
  meta  %>% 
    as_tibble(rownames="cell_id") %>%
    select(cell_id, {{my_cluster}}, manual_labels) %>%
    rename(cluster={{my_cluster}}) %>%
    mutate(label_src="manual")
}


add_manual_labels = function(sobj){
  orig = readRDS(sprintf("%s/analyzed_seurat_object.rds", JDM_DIR))
  orig_meta = orig@meta.data 
  cl_to_manual = orig_meta %>% 
    as_tibble(rownames="cell_id") %>%
    select(cell_id, manual_labels) %>%
    mutate(manual_labels=as.character(manual_labels))
  
  sobj = addMetaJoin(sobj, cl_to_manual, "manual_labels")
  return(sobj)
}

add_azimuth_labels = function(sobj){
  load(sprintf("%s/azimuth_unknown_jdm.RData", JDM_DIR))
  unknown_meta2 = unknown_meta %>%
    map_azimuth_label(predicted.celltype.l2, azimuth_label_mapped) %>%
    as_tibble(rownames="cell_id") %>%
    select(cell_id, azimuth_label_mapped)
  sobj = addMetaJoin(sobj, unknown_meta2, "azimuth_label_mapped")
  return(sobj)
}


transfer_annot = function(sobj, my_cluster, new_cluster_label, min.frac=0.4,
                          overwrite=F){
  meta0 = reform_meta(sobj@meta.data, my_cluster)
  
  clus_to_lab = meta0 %>%
    group_by(cluster) %>%
    mutate(tot=n()) %>%
    group_by(cluster, manual_labels, tot) %>%
    mutate(n_manual=sum(label_src=="manual")) %>%
    group_by(cluster, manual_labels, tot, n_manual) %>%
    dplyr::count() %>%
    mutate(frac_manual=n_manual/n) %>%
    arrange(cluster, desc(n)) %>%
    mutate(frac=n/tot)
  
  clus_to_lab2 = clus_to_lab %>%
    ungroup() %>%
    group_by(cluster) %>%
    slice_max(order_by=frac, n=1)
  
  clus_to_lab3 = clus_to_lab2 %>%
    mutate(manual_labels=as.character(manual_labels)) %>%
    mutate(assigned_label=ifelse(frac < min.frac | is.na(manual_labels), sprintf("%s*", manual_labels), 
                                 manual_labels)) %>%
    arrange(assigned_label) %>%
    ungroup() %>%
    group_by(assigned_label) %>%
    mutate(idx=1:n(),
           nlab=n()) %>%
    mutate({{new_cluster_label}}:=ifelse(nlab ==1, assigned_label,
                                    paste(assigned_label, idx))) %>%
    ungroup()

  clus_map = clus_to_lab3 %>% select(cluster, {{new_cluster_label}}) %>% distinct()
  stopifnot(nrow(clus_map)==length(unique(meta0$cluster)))
  
  sobj = addMetaJoin(sobj, meta0 %>% left_join(clus_map, by="cluster"), new_cluster_label,
                     overwrite)
  return(sobj)
}


