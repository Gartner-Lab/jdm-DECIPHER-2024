# Code from Arjun for generating profile plots for QC
library(tidyverse)
library(assertthat)  # For when we want to make sanity checks
library(grid)        # For plotting multiple plots in one frame
library(gridExtra)   # For plotting multiple plots in one frame
library(scales)      # To access break formatting functions


generate_profile_plot <- function(sobj, feature1, feature2, feature1_binwidth=100,
                                  feature2_binwidth=100, visual_outlier_cutoff1=0.999,
                                  visual_outlier_cutoff2=0.999) {
  suppmsg <- assert_that(feature1 %in% colnames(sobj@meta.data), 
                         msg=paste0(feature1, " was not present in the metadata of sobj"))
  suppmsg <- assert_that(feature2 %in% colnames(sobj@meta.data), 
                         msg=paste0(feature2, " was not present in the metadata of sobj"))
  suppmsg <- assert_that(0 < visual_outlier_cutoff1 && visual_outlier_cutoff1 <=1.0, 
                         msg="visual_outlier_cutoff1 must be in the range (0,1]")
  suppmsg <- assert_that(0 < visual_outlier_cutoff2 && visual_outlier_cutoff2 <=1.0, 
                         msg="visual_outlier_cutoff2 must be in the range (0,1]")
  
  lay <- rbind(c(1,  1,2,2,2,2),
               c(1,  1,2,2,2,2),
               c(1,  1,2,2,2,2),
               c(NA,NA,3,3,3,3),
               c(NA,NA,3,3,3,3))
  
  lims = as.vector(
    c(quantile(sobj@meta.data[[feature1]], visual_outlier_cutoff1), 
      quantile(sobj@meta.data[[feature2]], visual_outlier_cutoff2)))
  xticks <- as.vector(quantile(sobj@meta.data[[feature1]], seq(0, max(0.9, visual_outlier_cutoff1), 0.1)))
  if (xticks[length(xticks)] != lims[1]) {
    xticks <- c(xticks, lims[1])
  }
  yticks <- as.vector(quantile(sobj@meta.data[[feature2]], seq(0, max(0.9, visual_outlier_cutoff2), 0.1)))
  if (yticks[length(yticks)] != lims[2]) {
    yticks <- c(yticks, lims[2])
  }
  main <- ggplot(sobj@meta.data, aes_string(x=feature1, y=feature2)) + 
    geom_point(aes(col="red"), size=0.5) + 
    xlim(NA, lims[1]) +
    ylim(NA, lims[2]) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()
    ) +
    NoLegend()
  
  y_hist <- ggplot(sobj@meta.data, aes_string(x=feature2)) + 
    geom_histogram(aes(col="red", fill="red"), binwidth=feature2_binwidth) +
    theme_bw() +
    theme(axis.title.x=element_blank(), 
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank()) +
    scale_x_continuous(limits=c(NA, lims[2]),
                       sec.axis = sec_axis(trans = ~., 
                                           breaks=yticks, labels=NULL)) +
    NoLegend() + 
    coord_flip() + 
    scale_y_reverse()
  
  
  x_hist <- ggplot(sobj@meta.data, aes_string(x=feature1)) + 
    geom_histogram(aes(col="red", fill="red"), binwidth=feature1_binwidth) +
    theme_bw() +
    theme(axis.title.y=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) +
    scale_x_continuous(limits=c(NA, lims[1]),
                       sec.axis = sec_axis(trans = ~., 
                                           breaks=xticks, labels=NULL)) +
    
    NoLegend() +
    scale_y_reverse()
  empty_plot <- ggplot() + theme_void()
  grid.arrange(grobs=list(y_hist, main, x_hist), layout_matrix = lay)
}



generate_all_profile_plots = function(sobj, prefix){
  plots <- list(p1=generate_profile_plot(sobj,
                                         feature1 = "nCount_RNA",
                                         feature2 = "percent.mt",
                                         feature1_binwidth=100,
                                         feature2_binwidth=0.1),
                p2=generate_profile_plot(sobj,
                                         feature1 = "nCount_RNA",
                                         feature2 = "percent.ribo",
                                         feature1_binwidth=100,
                                         feature2_binwidth=0.1),
                p3=generate_profile_plot(sobj,
                                         feature1 = "nCount_RNA",
                                         feature2 = "nFeature_RNA",
                                         feature1_binwidth=100,
                                         feature2_binwidth=100),
                p4=generate_profile_plot(sobj,
                                         feature1 = "percent.ribo",
                                         feature2 = "percent.mt",
                                         feature1_binwidth=0.1,
                                         feature2_binwidth=0.1),
                p5=generate_profile_plot(sobj,
                                         feature1 = "nFeature_RNA",
                                         feature2 = "percent.mt",
                                         feature1_binwidth=100,
                                         feature2_binwidth=0.1),
                p6=generate_profile_plot(sobj,
                                         feature1 = "percent.ribo",
                                         feature2 = "nFeature_RNA",
                                         feature1_binwidth=0.1,
                                         feature2_binwidth=100),
                p7=generate_profile_plot(sobj,
                                        feature1 = "nCount_ADT",
                                        feature2 = "percent.mt",
                                        feature1_binwidth=100,
                                        feature2_binwidth=0.1),
                p8=generate_profile_plot(sobj,
                                         feature1 = "nCount_ADT",
                                         feature2 = "percent.ribo",
                                         feature1_binwidth=100,
                                         feature2_binwidth=0.1),
                p9=generate_profile_plot(sobj,
                                         feature1 = "nCount_ADT",
                                         feature2 = "nFeature_RNA",
                                         feature1_binwidth=100,
                                         feature2_binwidth=100),
                p10=generate_profile_plot(sobj,
                                         feature1 = "nCount_ADT",
                                         feature2 = "nCount_RNA",
                                         feature1_binwidth=100,
                                         feature2_binwidth=100),
                p10=generate_profile_plot(sobj,
                                         feature1 = "nCount_ADT",
                                         feature2 = "nFeature_ADT",
                                         feature1_binwidth=100,
                                         feature2_binwidth=1)
  )
  df <- data.frame(
    cell_counts=seq(0, 1.01, 0.1)*dim(sobj@meta.data)[1],
    percent.mt=quantile(sobj@meta.data[["percent.mt"]], seq(0, 1.01, 0.1)),
    percent.ribo=quantile(sobj@meta.data[["percent.ribo"]], seq(0, 1.01, 0.1)),
    nCount_RNA=quantile(sobj@meta.data[["nCount_RNA"]], seq(0, 1.01, 0.1)),
    nFeature_RNA=quantile(sobj@meta.data[["nFeature_RNA"]], seq(0, 1.01, 0.1)),
    nCount_ADT=quantile(sobj@meta.data[["nCount_ADT"]], seq(0, 1.01, 0.1)),
    nFeature_ADT=quantile(sobj@meta.data[["nFeature_ADT"]], seq(0, 1.01, 0.1)),
    row.names=seq(0, 1.01, 0.1)
  )
  write.table(format(df, digits=2), file=paste0(prefix, "_dualscatter_pre.tsv"), 
              row.names=T, col.names=T, quote=F, sep="\t")
  pdf(paste0(prefix, "_dualscatter_pre.pdf"), 
      width = 21, height = 14)
  print(CombinePlots(plots = plots, ncol=3))
  dev.off()
  
}