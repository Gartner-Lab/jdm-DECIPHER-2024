source(sprintf("%s/scripts/utils/FuncPctDotPlot.R", CODE_DIR))



gen_dotplot = function(sobj, my_assay, markers, nmarkers=5, cluster=T){

  top5_genes = markers %>%
    group_by(cluster) %>%
    arrange(cluster) %>%
    slice_max(n=nmarkers, order_by=avg_log2FC) 
  top5_genes$gene = factor(top5_genes$gene, levels=unique(top5_genes$gene ))
  print(length(unique(top5_genes$gene)))
  if (my_assay %in% c("RNA", "SCT")){
    DotPlot(sobj, assay=my_assay, features = unique(top5_genes$gene),
                   cols="RdYlBu",
                   cluster.idents = cluster)+
      coord_flip()+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
  } else{
    FuncPctDotPlot(sobj, assay=my_assay, features = unique(top5_genes$gene),
                   cols="RdYlBu",
                   cluster.idents = cluster)+
      coord_flip()+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))
  }
}

make_wnn_plots = function(sobj, col_vec){
  p1 = DimPlot(sobj, reduction="wnn.umap", cols=col_vec, repel=T, label=T, label.size=3)+NoLegend()
  p2 = DimPlot(sobj, reduction="umap", cols=col_vec, label=T, repel=T, label.size=3)+
    xlab("rnaUMAP_1")+ylab("rnaUMAP_2")+NoLegend()
  p3 = DimPlot(sobj, reduction="a.umap", cols=col_vec, label=T, repel=T, label.size=3)+
    xlab("adtUMAP_1")+ylab("adtUMAP_2")+NoLegend()

  plot_grid(p1, p2, p3, ncol=3)
}

make_wnn_wt_plot = function(sobj, col_vec){
  VlnPlot(sobj, features="RNA.weight", sort=T, pt.size=0,
               cols=col_vec)+NoLegend()+
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2))+
    theme(axis.title.x = element_blank())
}


make_doublet_vln_plots = function(sobj, col_vec, cluster){
  DefaultAssay(sobj) = "integrated.ADT"
  p1 = VlnPlot(sobj, features=c("CD4"), cols=col_vec, 
               group.by = cluster, sort = FALSE, pt.size = 0) +
    NoLegend()+
    theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  p2 = VlnPlot(sobj, features=c("CD8"), cols=col_vec, 
               group.by = cluster, sort = FALSE, pt.size = 0) +
    NoLegend()+
    theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  p3 = VlnPlot(sobj, features=c("CD16"), cols=col_vec, 
               group.by = cluster, sort = FALSE, pt.size = 0) +
    NoLegend()+
    theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  cowplot::plot_grid(p1, p2, p3, ncol=1)
}



make_qc_vln_plots = function(sobj, col_vec, cluster){
  p1 = VlnPlot(sobj, features=c("nFeature_RNA"), cols=col_vec, 
               group.by = cluster, sort = FALSE, pt.size = 0) +
    NoLegend()+
    theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  p2 = VlnPlot(sobj, features=c("nCount_RNA"), cols=col_vec, 
               group.by = cluster, sort = FALSE, pt.size = 0) +
    NoLegend()+
    theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  p4= VlnPlot(sobj, features=c("nFeature_ADT"), cols=col_vec, 
              group.by = cluster , sort = FALSE, pt.size = 0) +
    NoLegend()+
    theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  p5 = VlnPlot(sobj, features=c("nCount_ADT"), cols=col_vec, 
               group.by = cluster, sort = FALSE, pt.size = 0) +
    NoLegend()+
    theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  p3 =VlnPlot(sobj, features=c("percent.mt"), cols=col_vec, 
              group.by = cluster, sort = FALSE, pt.size = 0) +
    NoLegend()+
    theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  p6 = VlnPlot(sobj, features=c("percent.ribo"), cols=col_vec, 
               group.by = cluster, sort = FALSE, pt.size = 0) +
    NoLegend()+
    theme(axis.text.x = element_text(angle=90, vjust=0.2, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=3)
  
}

make_meta_barplot = function(sobj, feature, cluster, col_vec){
  sobj@meta.data %>%
    select({{feature}}, {{cluster}}) %>%
    rename(cluster={{cluster}}) %>%
    ggplot(aes(x=cluster, fill={{feature}}))+
    geom_bar(position="fill")+
    theme_bw()+
    scale_fill_manual(values=dittoColors(length(unique(sobj@meta.data %>% pull({{feature}})))))+
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2),
          axis.title.x=element_blank())+
    ylab("fraction of cells")
}


plot_clus_sizes = function(sobj,  col_vec, cluster){
  meta = sobj@meta.data %>%
    rename(cluster={{cluster}}) %>%
    group_by(cluster) %>%
    count() %>%
    arrange(desc(n))
  meta$cluster = factor(meta$cluster, levels=meta$cluster)
  
  ggplot(meta, aes(x=cluster, y=n, fill=cluster))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=col_vec)+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2),
          axis.title.x=element_blank(),
          legend.position = "none")+
    ylab("number of cells")
}
