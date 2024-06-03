#function to run proxy module score to validate programs in external data

proxy_AUCscore <- function(sobj, type, nmf_type, markers, program, ...){
  require(Seurat); require(AUCell); require(tidyverse); require(dplyr)
  suppressWarnings({
  cells <- subset(sobj, broad_label == type)
  counts <- cells@assays$RNA@counts
  rankings <- AUCell_buildRankings(counts, plotStats=F)
  markers_tbl <- rownames_to_column(markers[[nmf_type]], var = 'gene') %>% as_tibble()
  program_markers <- markers_tbl[, c('gene', program)]
  program_markers[,"quantile"] <- percent_rank(program_markers[,program])
  program_markers <- subset(program_markers, subset = quantile > 0.95)
  genes = program_markers[["gene"]]
  
  #calc AUC using input genes and cell type rankings
  AUC <- AUCell_calcAUC(genes, rankings, aucMaxRank = ceiling(0.05 * nrow(rankings)), verbose = F)
  
  # Extract the scores and add to metadata
  AUCell_res <- as.data.frame(t(getAUC(AUC)))
  auc_label <- paste0(nmf_type, regmatches(program,regexpr('Program[0-9]{1,2}', program)))
  colnames(AUCell_res) <- c(auc_label)
  cells <- AddMetaData(cells, AUCell_res)
  
  #plot umap
  directory = '~/Library/CloudStorage/Box-Box/neely.jdm/NMF/validation_analyses/JDM_validation_cohort/AUCell_Results/'
  pdf(paste0(directory, 'jdm_val_', auc_label,'_AUC_umap.pdf'), height = 14, width =8)
  DefaultAssay(cells) <- 'RNA'
  print(FeaturePlot(cells, features = auc_label, reduction = 'umap') + DimPlot(cells, group.by = 'case_ctrl', reduction = 'umap'))
  dev.off()
  
  #plot boxplots
  means <-cells@meta.data %>% group_by(individual)  %>% mutate(AUC = mean(!!rlang::sym(auc_label)))
  scores <- dplyr::distinct(means[, c("individual", "case_ctrl", "AUC")])
  scores <- drop_na(scores) %>% select(individual, case_ctrl, AUC) %>% setNames(c('individual', 'case_ctrl', auc_label))
  
  pdf(paste0(directory, 'jdm_val_', auc_label, '_boxplot.pdf'), height = 8, width =8)
  plot <- ggplot(transform(scores, xjit = jitter(as.numeric(as.factor(case_ctrl)))), 
         aes(x = case_ctrl, y = rlang::sym(auc_label), color = case_ctrl)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(x = xjit)) +
    theme_classic() +
    theme(legend.position = 'none') +
    ggpubr::stat_compare_means(label = "p.signif", method = "t.test",
                       ref.group = "HC") 
  print(plot)
  dev.off()
  })
  return(scores)
}
