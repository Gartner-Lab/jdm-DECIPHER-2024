#Loading packages
library(Seurat)
library(tidyverse)
library(dittoSeq)
library(scales)
library(cowplot)
library(openxlsx)
library(MAST)
library(Seurat)
library(tidyverse)
library(xlsx)
library(scales)
library(vctrs)

#Laoding data
sobj_final <- readRDS(file = "/home/ubuntu/citeseq_21/clustered_seurat_object_8_11_22/sobj_final.rds")

## Figure 2C
th <-   theme(axis.text.x = element_text(color = "black", size = 18, face = 'bold'), 
              strip.text = element_text(color = "black", size = 10),
              axis.text.y.left = element_text(size = 18, color = 'black', face = 'bold'),
              plot.title = element_text(size = 17, face = 'bold'),
              legend.title = element_text(size = 18, face = 'bold'),
              plot.subtitle = element_text(size = 16, face = 'bold'),
              legend.text = element_text(size = 18))
#Setting idents
Idents(sobj_final) <- sobj_final$new_labels

#Making combined meta-info using diseasegroups and study_id
disease_groups <- sobj_final$new_diseasegroups
donors <- sobj_final$paper_id
new <- data.frame(disease_groups, donors)
new$new <- paste0(new$disease_groups, "_", new$donors)
donor_group <- new$new
Idents(sobj_final) <- donor_group
sobj_final[["donor_group"]] <- Idents(sobj_final)

#Getting information on % cell type pr study_id
dittoBarPlot(sobj_final, "new_labels", group.by = "donor_group")
cell_prop_sample <- dittoBarPlot(sobj_final, "new_labels", group.by = "donor_group",data.out = TRUE)$data
cell_prop_sample[c("disease_group", "donor")] <- str_split_fixed(cell_prop_sample$grouping, "_" , 2)

# Add VAS global by looping through celltypes and ids, and grabbing scores for all scores per iteration
cell_prop_sample$vasglobal <- NA
for (i in 1:nrow(cell_prop_sample)) {
  cell_prop_sample$vasglobal[i] <- sobj_final$vasglobal[match(cell_prop_sample$grouping[i], sobj_final$donor_group)]
}

#Add vascutaneous
cell_prop_sample$vascutaneous <- NA
for (i in 1:nrow(cell_prop_sample)) {
  cell_prop_sample$vascutaneous[i] <- sobj_final$vascutaneous[match(cell_prop_sample$grouping[i], sobj_final$donor_group)]
}

#Add vasmuscle
cell_prop_sample$vasmuscle <- NA
for (i in 1:nrow(cell_prop_sample)) {
  cell_prop_sample$vasmuscle[i] <- sobj_final$vasmuscle[match(cell_prop_sample$grouping[i], sobj_final$donor_group)]
}


# Add rough annotation
cell_prop_sample$rough <- NA
for (i in 1:nrow(cell_prop_sample)) {
  cell_prop_sample$rough[i] <- as.character(sobj_final$rough_annotation[match(cell_prop_sample$label[i], sobj_final$new_labels)])
}
cell_prop_sample$rough <- as.factor(cell_prop_sample$rough)

levels <- levels(sobj_final$new_labels)
cell_prop_sample$label <- factor(cell_prop_sample$label, levels = levels)


#Remove B_naive4 and HC
cell_prop_sample <- filter(cell_prop_sample, label != 'B_naive4')
cell_prop_sample2 <- subset(cell_prop_sample, disease_group != "HC")


#Calculate correlation between cell type proportion and different disease measures
types <- c('vasglobal', 'vascutaneous', 'vasmuscle')
for_plot <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(for_plot) <- c("prop_cor", 'prop_p','prop_p_adj', 'cell', 'disease measure')
for(n in types){
  prop_cor  <- unlist(sapply(
    unique(cell_prop_sample2$label),
    function(cell_type) {
      t <- subset(cell_prop_sample2, label==cell_type)
      t <- na.omit(t)
      out <- cor(t[[n]],
                 t$percent,
                 method = "spearman",
                 use = "complete.obs")
      names(out) <- paste(cell_type)
      out
    },
    simplify = FALSE))
  
  prop_p<- unlist(sapply(
    unique(cell_prop_sample2$label),
    function(cell_type) {
      t <- subset(cell_prop_sample2, label==cell_type)
      t <- na.omit(t)
      out <- round(cor.test(t[[n]],
                            t$percent,
                            method = "spearman", use = "complete.obs", exact = FALSE)$p.value, 3)
      names(out) <- paste(cell_type)
      out
    },
    simplify = FALSE))
  
  corr_comp <- data.frame(prop_cor, prop_p, check.names = TRUE)
  corr_comp$cell <- rownames(corr_comp)
  corr_comp$prop_p_adj <- p.adjust(prop_p, method = 'fdr')
  rownames(corr_comp) <- NULL
  corr_comp['scale'] <- n
  for_plot <- rbind(for_plot, corr_comp)
}

#Prepare dataframe for visualization
for_plot$Correlation <- for_plot$prop_cor
for_plot$Significant <- as.factor(ifelse(for_plot$prop_p_adj < 0.05, "TRUE", "FALSE"))
for_plot['Absolut correlation'] <- abs(for_plot$prop_cor)
for_plot$scale[for_plot$scale == 'vascutaneous'] <- 'Skin'
for_plot$scale[for_plot$scale == 'vasglobal'] <- 'Global'
for_plot$scale[for_plot$scale == 'vasmuscle'] <- 'Muscle'
for_plot$scale <- as.factor(for_plot$scale)
for_plot$cell <- as.factor(for_plot$cell)
levels <- levels(cell_prop_sample$label)
for_plot$cell <- factor(for_plot$cell, levels = levels)

# Add rough annotation
for_plot$rough <- NA
for (i in 1:nrow(for_plot)) {
  for_plot$rough[i] <- as.character(sobj_final$rough_annotation[match(for_plot$cell[i], sobj_final$new_labels)])
}
cell_prop_sample$rough <- as.factor(cell_prop_sample$rough)
levels <- levels(sobj_final$new_labels)
cell_prop_sample$label <- factor(cell_prop_sample$label, levels = levels)


#Filter cells without significant differences
cell_prop_sample3 <- filter(cell_prop_sample, ! label %in% c('gdT_c1', 'CD4+ naive', 
                                                             'PBs',
                                                             'CD8+ GZMKhi', 'CD8+ GZMA/Bhi', 
                                                             'CD8+ mem', 'CD14+ mono', 
                                                             'CD16+ mono', 'PDCs', 'CD8+ naive', 'B_mem'))
for_plot3 <- filter(for_plot, ! cell %in% c('gdT_c1', 'CD4+ naive', 
                                            'PBs',
                                            'CD8+ GZMKhi', 'CD8+ GZMA/Bhi', 
                                            'CD8+ mem', 'CD14+ mono',
                                            'CD16+ mono', 'PDCs', 'CD8+ naive', 'B_mem'))

#Add overall cell categoreis
cell_prop_sample3$rough <- case_when(cell_prop_sample3$rough == 'CD4 T cells' ~ 'T cells',
                                     cell_prop_sample3$rough == 'CD8 T cells' ~ 'T cells', 
                                     cell_prop_sample3$rough == 'gdT' ~ 'T cells',
                                     cell_prop_sample3$rough == 'NK' ~ 'Myeloid',
                                     cell_prop_sample3$rough == 'Monocytes' ~ 'Myeloid',
                                     .default = as.character(cell_prop_sample3$rough))
for_plot3$rough <- case_when(for_plot3$rough == 'CD4 T cells' ~ 'T cells',
                             for_plot3$rough == 'CD8 T cells' ~ 'T cells', 
                             for_plot3$rough == 'gdT' ~ 'T cells',
                             for_plot3$rough == 'NK' ~ 'Myeloid',
                             for_plot3$rough == 'Monocytes' ~ 'Myeloid',
                             .default = as.character(for_plot3$rough))

#Prepare dataframes for visualization
cell_prop_sample3$disease_group <- factor(cell_prop_sample3$disease_group, levels = c('TNJDM', 'Active', 'Inactive', 'HC'))
for_plot3$rough <- as.factor(for_plot3$rough)
for_plot3$scale <- factor(for_plot3$scale, levels = c('Global', 'Muscle', 'Skin'))
for_plot3$cell <- droplevels(for_plot3$cell)
for_plot3$cell <- factor(for_plot3$cell, levels = c('B_naive1', 'B_naive2', 'B_naive3','B_mem',
                                                    'CD56dim NK', 'CD56bright NK', 'CDCs',
                                                    'gdT_c2', 'CD4+ eff', 'CD4+ Tregs', 'CD8+ naive'))
cell_prop_sample3$rough <- as.factor(cell_prop_sample3$rough)
cell_prop_sample3$label <- droplevels(cell_prop_sample3$label)
cell_prop_sample3$label <- factor(cell_prop_sample3$label, levels = c('B_naive1', 'B_naive2', 'B_naive3','B_mem',
                                                                      'CD56dim NK', 'CD56bright NK', 'CDCs',
                                                                      'gdT_c2', 'CD4+ eff', 'CD4+ Tregs', 'CD8+ naive'))

#Make plots
plot1 <- ggplot(cell_prop_sample3, aes(label, percent, fill = disease_group)) +
  geom_boxplot(lwd = 0.2) +
  geom_point(position = position_dodge(width = 0.75)) +
  theme_bw() +
  theme(axis.text.x = element_text(face = 'bold', size = 18, angle = 25, hjust = 1, vjust = 1),
        axis.title.y.left = element_text(size = 15, colour = 'black', face = 'bold'),
        axis.title = element_blank(),
        legend.position = 'top') +
  ylab('log(percent)') +
  th +
  guides(fill = guide_legend(title = 'Disease groups')) +
  scale_y_log10() +
  scale_fill_manual(values = c('HC' = '#5CB4E5',
                               'Inactive' = '#7BAF41',
                               'TNJDM' = '#F3756D','Active' = '#A680BA'))


plot2 <- ggplot(for_plot3, aes(scale, cell, size = `Absolut correlation`, fill = Correlation, colour = Significant)) +
  geom_point(shape = 21, stroke = 1) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.9,0.9)) +
  scale_colour_manual(values = c("grey", "black")) +
  scale_size_continuous(range = c(1,7), breaks = seq(0,1,0.2)) +
  theme_bw() +  
  #th +
  theme(legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.position = 'top',
        axis.text.x.bottom = element_blank()) +
  labs(size = 'Absolut\ncorrelation') +
  th +
  ylab('') +
  xlab('') +
  coord_flip()

#Plot plots
plot_grid(plot2, plot1,
          align = 'hv', axis = 'l', nrow = 2, ncol = 1, rel_heights = c(1,2))

#Calculating Dunn's post test, comparing between HC, TNJDM and inactive
cell_prop_sample$label <- droplevels(cell_prop_sample$label) 
cell_prop_sample <- cell_prop_sample %>% filter(disease_group != 'Active') 
cell_prop_sample$disease_group <- as.factor(cell_prop_sample$disease_group)
for(n in levels(cell_prop_sample$label)){
  t <- cell_prop_sample %>% filter(label == n) %>% ungroup()
  print(1)
  t$disease_group <- droplevels(t$disease_group)
  print(2)
  t$disease_group <- as.character(t$disease_group)
  print(2)
  comparisons <- list()
  print(4)
  means <- rstatix::dunn_test(t, formula = percent ~ disease_group)
  print(5)
  print(n)
  print(means)
}