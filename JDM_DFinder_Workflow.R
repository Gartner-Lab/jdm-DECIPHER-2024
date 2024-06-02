
Idents(data) <- data$well
split <- SplitObject(data, split.by = 'well')
saveRDS(split, file = '~/citeseq_21/split_object.rds')

#run the below lines to save each lane as a separate seurat object (easier for downstream functions instead of reading entire dataset in at the beginning)
lane_IDs <- as.factor(unique(data@metadata$well)) #save lane/well IDs as list of strings
saveRDS(lane_IDs, file = '~/citeseq_21/10X_laneIDs.rds')
data_split <- readRDS('~/citeseq_21/split_object.rds')
index <- list(1:length(data_split))
lapply(index, function(x){
  seu <- data_split[[x]]
  saveRDS(seu, file = paste0('~/citeseq_21/seu_well', lane_IDs[[x]]), '.rds')
  return(NULL)
})
rm(data_split) #remove list of subsetted objects as well as object with entire dataset from environment before running the below lines
#every other line of code reads in the corresponding subsetted object using the filename as the input (this preserves more memory since R automatically wipes the environment of each function's internal objects)
subsetted_objects <- grep('seu_well', list.files('~/citeseq_21/'), value = TRUE) #list with subsetted object filenames
dir.create('~/citeseq_21/doublet_finder_outputs/') #directory to save intermediate and final objects generated during doubletfinder workflow

#source custom doubletFinder function in accompanying script
source('~/citeseq_21/doublet_finder_outputs/conserve_mem_DoubletFinder.R') #I just added the 'conserve.memory = TRUE' option to the SCT processing line of the doublet finder function

## pK Identification step
lane_pKs <- lapply(subsetted_objects, function(x){
  seu <- readRDS(paste0('~/citeseq_21/', x))
  sweep.res.list <- paramSweep_v3(seu, PCs = 1:30, sct = TRUE) #may need to add 'conserve.memory' option to paramSweep function as well
  saveRDS(paste0('~/citeseq_21/doublet_finder_outputs/sweep.res_', x ))
  gt.calls <- data@meta.data[rownames(sweep.res.list[[1]]), 'DROPLET.TYPE']   ## GT is a vector containing “Singlet” and “Doublet” calls recorded using sample multiplexing classification and/or in silico geneotyping results
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE, GT.calls = gt.calls)
  bcmvn <- find.pK(sweep.stats)
  pk_choose <- bcmvn[which.max(bcmvn$BCmetric), 'pK']
})
names(lane_pKs) <-lane_IDs 
saveRDS(lane_pKs, file = '~/citeseq_21/doublet_finder_outputs/lane_pKs.rds') 

## Homotypic Doublet Proportion Estimate
lanes_nExp_poi.adj <- lapply(subsetted_objects, function(x){
  data <- readRDS(paste0('~/citeseq_21/', x))
  homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters) 
  nExp_poi <- round(homotypic.prop*nrow(data@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
})
names(lanes_nExp_poi.adj) <- lane_IDs
saveRDS(lanes_nExp_poi.adj, file = '~/citeseq_21/doublet_finder_outputs/lane_nExp.adj.rds')

#run doublet finder function on each subsetted object
lapply(index, function(x){
  data <- readRDS(paste0('~/citeseq_21/', subsetted_objects[[x]]))
  data <- doubletFinder_v3(data, PCs = 1:30, pN = 0.25, pK = lane_pKs[[x]], nExp = lanes_nExp_poi.adj[[x]], reuse.pANN = FALSE, sct = TRUE)
  saveRDS(paste0('~/citeseq_21/doublet_finder_outputs/doub_found_lane', lane_IDs[[x]], '_seu.sub.rds'))
  return(NULL)
})

