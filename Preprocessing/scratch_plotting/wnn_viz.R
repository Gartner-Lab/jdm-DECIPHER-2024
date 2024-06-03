
load(sprintf("%s/wnn_obj.RData", OUT_DIR)) # --> sobj

VlnPlot(sobj, features = "RNA.weight", 
        group.by = 'wsnn0.4', sort = TRUE, pt.size = 0) +
  NoLegend()
ggsave(sprintf("%s/wnn_weight0.4.png", OUT_DIR))
       
       VlnPlot(sobj, features = "RNA.weight", 
               group.by = 'wsnn0.4', sort = FALSE, pt.size = 0) +
         NoLegend()
       ggsave(sprintf("%s/wnn_weight0.4_unsrt.png", OUT_DIR))
       
       VlnPlot(sobj, features = "PF4", 
               group.by = 'wsnn0.4', sort = FALSE, pt.size = 0) +
         #  ylim(0, 1000)+
         NoLegend()
       
       sobj@meta.data %>%
         group_by(wsnn_res.1) %>%
         dplyr::count() %>%
         filter(n>=10) %>% 
         nrow()
       
       
       #######
       ## label the clusters
       jdm_obj_initial = readRDS(sprintf("%s/analyzed_seurat_object (1).rds", JDM_DIR))
       meta = jdm_obj_initial@meta.data %>% 
         as_tibble(rownames="cell_id")
       rm(jdm_obj_initial)
       length(unique(sobj@meta.data$wsnn0.4)) # 33
       clus_lab = sobj@meta.data %>% 
         as_tibble(rownames="cell_id") %>%
         dplyr::select(cell_id, wsnn0.4) %>%
         left_join(  
           meta %>% 
             dplyr::select(cell_id, manual_labels))
       length(unique(clus_lab$wsnn0.4))
       clus_to_lab = clus_lab %>%
         dplyr::rename(cluster="wsnn0.4") %>%
         group_by(cluster) %>%
         mutate(tot=n()) %>%
         group_by(cluster, manual_labels, tot) %>%
         dplyr::count() %>%
         arrange(cluster, desc(n)) %>%
         mutate(frac=n/tot)
       length(unique(clus_to_lab$cluster)) # 33
       
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
       
       
       new_data = sobj@meta.data %>% 
         as_tibble(rownames="cell_id") %>%
         select(cell_id, wsnn0.4) %>%
         left_join(clus_to_lab3 %>% select(cluster, cluster_label) %>%
                     mutate(cluster_label=as.factor(cluster_label)), by=c("wsnn0.4"="cluster"))
       
       manual_lab_df = clus_lab %>% 
         mutate(manual_labels=as.character(manual_labels),
                manual_labels=ifelse(is.na(manual_labels), "unknown",
                                     manual_labels)) 
       
       sobj = AddMetaData(sobj, new_data$cluster_label, col.name="cluster_label")
       sobj = AddMetaData(sobj, manual_lab_df$manual_labels, col.name="manual_labels")
       tiny_clusters = clus_to_lab %>%
         ungroup() %>%
         filter(tot < 10) %>% distinct(cluster) %>% pull(cluster)
       keep.cls = rownames(sobj@meta.data %>% filter(!wsnn0.4 %in% tiny_clusters))
       dim(sobj); length(keep.cls) # removing 22 cells
       sobj = subset(sobj, cells=keep.cls)
       
       Idents(sobj) = "cluster_label"
       ## remake the plots
       DimPlot(sobj, reduction = 'wnn.umap', 
               label = TRUE, cols=dittoColors(),
               label.size=3) + NoLegend()
       ggsave(sprintf("%s/wnn_umap_clus_lab.pdf", OUT_DIR))
       
       
       VlnPlot(sobj, features = "RNA.weight", 
               group.by = 'cluster_label', cols=dittoColors(),
               sort = TRUE, pt.size = 0) +
         NoLegend()
       ggsave(sprintf("%s/wnn_umap_clus_lab_weight.pdf", OUT_DIR))
       
       
       # look at unknown 2,4,5,6
       
       
       # compare to UMAP of RNA, ADT
       sobj <- RunUMAP(sobj, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                       reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
       sobj <- RunUMAP(sobj, reduction = 'a.pca', dims = 1:18, assay = 'ADT', 
                       reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
       
       Idents(sobj) = "manual_labels"
       
       p1 = DimPlot(sobj, reduction = 'wnn.umap', 
                    label = TRUE, cols=dittoColors()[1:30],
                    repel=T,
                    label.size=3) + NoLegend()
       p2 = DimPlot(sobj, reduction = 'rna.umap', cols=dittoColors()[1:30], label = TRUE, 
                    repel = TRUE, label.size = 3) + NoLegend()
       p3 = DimPlot(sobj, reduction = 'adt.umap', cols=dittoColors()[1:30], label = TRUE, 
                    repel = TRUE, label.size = 3) + NoLegend()
       cowplot::plot_grid(p1,p2,p3, ncol=3)
       ggsave(sprintf("%s/wnn_umap_compare.pdf", OUT_DIR))
       
       Idents(sobj) = "cluster_label"
       
       # make results dataframe 
       d = cbind(sobj@meta.data, as.data.frame(t(sobj@assays$ADT@data)))
       
       # calculate the median protein expression separately for each cluster 
       adt_plot = d %>%
         dplyr::group_by(cluster_label) %>% 
         dplyr::summarize_at(.vars = rownames(sobj@assays[["ADT"]]), 
                             .funs = median) %>% 
         tibble::remove_rownames() %>% 
         tibble::column_to_rownames("cluster_label") 
       
       
       # plot a heatmap of the average dsb normalized values for each cluster
       adt_plot2 = adt_plot
       adt_plot2[adt_plot2>100]=100
       colmax = apply(adt_plot2, 2, max)
       adt_plot3=adt_plot2[,colmax >10]
       png(filename=sprintf("%s/wnn_adt_heatmap.png", OUT_DIR), width = 12, height = 23, units = "in", res = 300)
       print(pheatmap::pheatmap(t(adt_plot2), 
                                color = viridis::viridis(25, option = "B"), 
                                fontsize_row = 8, border_color = NA))
       dev.off()
       
       
       png(filename=sprintf("%s/wnn_adt_heatmap2.png", OUT_DIR), width = 12, height = 15, units = "in", res = 300)
       print(pheatmap::pheatmap(t(adt_plot3), 
                                color = viridis::viridis(25, option = "B"), 
                                fontsize_row = 8, border_color = NA))
       dev.off()
       
       
       uk6_markers = FindMarkers(sobj, ident.1="unknown 6", test.use="MAST", only.pos=T, assay="ADT")
       uk6_rna_markers = FindMarkers(sobj, ident.1="unknown 6", test.use="MAST", only.pos=T, assay="RNA")
       
       # look at unknown 2,4,5,6
       
       adt_markers = FindAllMarkers(sobj, test.use="MAST", only.pos=T, assay="ADT")
       top5_adt_genes = adt_markers %>% 
         group_by(cluster) %>%
         slice_max(n=5, order_by = avg_log2FC) %>%
         pull(gene) 
       source("/krummellab/data1/erflynn/scratch_scripts/FuncPctDotPlot.R")
       FuncPctDotPlot(sobj, assay="ADT", features = unique(top5_adt_genes),
                      cols="RdYlBu")+
         coord_flip()+
         theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
       ggsave(sprintf("%s/adt_dotplot.pdf", OUT_DIR), width=8, height=15)
       
       # what was the location of these cells in the RNA?
       DimPlot(sobj, reduction = 'rna.umap', 
               cells.highlight=rownames(sobj@meta.data %>% filter(cluster_label=="unknown 2")), 
               label.size = 3) + NoLegend()
       
       DimPlot(sobj, reduction = 'adt.umap', 
               cells.highlight=rownames(sobj@meta.data %>% filter(cluster_label=="unknown 2")), 
               label.size = 3) + NoLegend()
       
       DimPlot(sobj, reduction = 'adt.umap', 
               cells.highlight=rownames(sobj@meta.data %>% filter(cluster_label=="unknown 6")), 
               label.size = 3) + NoLegend()
       
       DimPlot(sobj, reduction = 'rna.umap', 
               cells.highlight=rownames(sobj@meta.data %>% filter(cluster_label=="unknown 6")), 
               label.size = 3) + NoLegend()
       
       # It would be helpful to see a UMAP (WNN) of the clusters themselves side-by-side 
       # with the manual labels, and then we can assign them together? 
       
       p1 = DimPlot(sobj, reduction="wnn.umap", group.by="wsnn0.4", label=T, cols=dittoColors())+NoLegend()
       p2 = DimPlot(sobj, reduction="wnn.umap", group.by="manual_labels", 
                    label=T, cols=c(dittoColors()[1:26],"light gray"),
                    label.size=2.5)+NoLegend()
       plot_grid(p1, p2, ncol=2)
       ggsave(sprintf("%s/wnn_umap_clusters_labels.pdf", OUT_DIR))
       
       p1.1 = DimPlot(sobj, reduction="wnn.umap", group.by="cluster_label", 
                      label=T, cols=dittoColors(), label.size=2.5)+NoLegend()
       plot_grid(p1.1, p2, ncol=2)
       ggsave(sprintf("%s/wnn_umap_clusters_labeled2_labels.pdf", OUT_DIR))
       
       
       # It would also be helpful to see umap colored by batch and study_id_visit 
       # variables to see if we still have those clusters separating out by patient 
       # with this method.
       
       DimPlot(sobj, reduction="wnn.umap", group.by="well")
       ggsave(sprintf("%s/wnn_umap_well.pdf", OUT_DIR))
       
       #meta %>% 
       head(sobj@meta.data)
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
       
       ## 
       Idents(sobj) = "wsnn0.4"
       high_adt_cts = sobj@meta.data %>%
         as_tibble(rownames="cell_id") %>%
         filter(nCount_ADT > 30000) %>%
         pull(cell_id)
       # only 43
       
       VlnPlot(sobj, "nFeature_ADT", pt.size=0)
       VlnPlot(sobj, "nCount_ADT", pt.size=0)
       DimPlot(sobj, cells.highlight=high_adt_cts, reduction="wnn.umap")
       FeaturePlot(subset(sobj, nCount_ADT < 30000), "nCount_ADT",
                   reduction="wnn.umap")
       