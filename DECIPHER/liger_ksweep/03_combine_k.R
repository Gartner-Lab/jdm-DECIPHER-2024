
library(tidyverse)

# COLLATE RESULTS FOR EASIER DATA GRABBING #### - do this later once KL calculated for all cell types at all ranks
# (probably could have set this up from the beginning as a list, but I think it was faster to run in parallel and collate later)
JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm"
setwd(JDM_DIR)

K.MIN=2
K.MAX=40

load("k_sweep/Liger_objects/cell.types.robj") # --> cell.types
#cell.types = setdiff(cell.types, "CD4T")
error_list = lapply(cell.types, function(cell_type){
  readRDS(sprintf("k_sweep/%s/consensus_res/error_list.rds", cell_type))
})
names(error_list) = cell.types
#error_list = list()
#error_list[["B"]] = readRDS('k_sweep/B/consensus_res/error_list.rds')
#error_list[["myeloid"]] = readRDS('k_sweep/myeloid/consensus_res/error_list.rds')
# error_list[["LP"]] = readRDS('k_sweep/LP/consensus_res/error_list.rds')
# error_list[["FB"]] = readRDS('k_sweep/FB/consensus_res/error_list.rds')
# error_list[["Blood_Endo"]] = readRDS('k_sweep/Blood_Endo/consensus_res/error_list.rds')

cell.types2 = cell.types[sapply(cell.types, function(cell_type) file.exists(sprintf("k_sweep/%s/H_metrics.rds", cell_type)))]

H_metrics_list = lapply(cell.types2, function(cell_type){
  print(cell_type)
  readRDS(sprintf("k_sweep/%s/H_metrics.rds", cell_type))
})
names(H_metrics_list) = cell.types2

#H_metrics_list = list()
#H_metrics_list[["B"]] = readRDS('k_sweep/B/H_metrics.rds')
#H_metrics_list[["myeloid"]] = readRDS('k_sweep/myeloid/H_metrics.rds')
# H_metrics_list[["LP"]] = readRDS('k_sweep/LP/H_metrics.rds')
# H_metrics_list[["FB"]] = readRDS('k_sweep/FB/H_metrics.rds')
# H_metrics_list[["Blood_Endo"]] = readRDS('k_sweep/Blood_Endo/H_metrics.rds')

dir.create('k_sweep/results_summary')
saveRDS(error_list, file = "k_sweep/results_summary/error_list.rds")
saveRDS(H_metrics_list, file = "k_sweep/results_summary/H_metrics_list.rds")


rank = K.MIN:K.MAX
names(rank) = paste0("R", rank)

H_consensus_list = list()
for (i in cell.types){
  H_consensus_list[[i]] = lapply(rank, function(k) {
    H = t(readRDS(file = paste0('k_sweep/', i, '/consensus_res/Liger.', i, '.consensus.R', k, '.rds'))$H)
  })
}
format(object.size(H_consensus_list), units = "Mb")
saveRDS(H_consensus_list, file = "k_sweep/results_summary/H_consensus_list.rds")

W_consensus_list = list()
for (i in cell.types){
  W_consensus_list[[i]] = lapply(rank, function(k) {
    W = t(readRDS(paste0('k_sweep/', i, '/consensus_res/Liger.', i, '.consensus.R', k, '.rds'))$W)
  })
}
format(object.size(W_consensus_list), units = "Mb")
saveRDS(W_consensus_list, file = "k_sweep/results_summary/W_consensus_list.rds")


# CALCULATE H SCORE (AVERAGE "EXPRESSION" ACROSS EACH SAMPLE) ####
H_consensus_list = readRDS('k_sweep/results_summary/H_consensus_list.rds')
RMs.metadata = readRDS('RMs.metadata.rds')

H.score = lapply(H_consensus_list, function(H){
  lapply(H, function(x){
    df = t(x)
    res = data.frame(matrix(nrow = length(unique(RMs.metadata$study_id_visit)), ncol = ncol(df)))
    rownames(res) = unique(RMs.metadata$study_id_visit)
    colnames(res) = colnames(df)
    for (i in unique(RMs.metadata$study_id_visit)){
      res[i,] = colMeans(df[intersect(rownames(df), rownames(subset(RMs.metadata, study_id_visit==i))),])
    }
    return(res)
  })
})

saveRDS(H.score, file = "k_sweep/results_summary/H.score.list.rds")


# PLOT SOME RESULTS ####
col2=c('#33A02B')
par(mfrow=c(1,1), mar=c(2.5,4.5,1,1))

list_lab = seq(5, K.MAX, 5)

kl_plot = function(cell_type){
  H_metrics <- readRDS(file = sprintf('k_sweep/%s/H_metrics.rds', cell_type))
  x = K.MIN:K.MAX
  y= unlist(lapply(H_metrics, `[[`, "KL"))
  plot(x, y, pch = 16, 
       type = "b", col = col2[1], lwd = 2, ylab = "KL divergence", xlab = "", 
       main = sprintf('%s KSweep KL Divergence', cell_type))
  points(list_lab, y[list_lab-K.MIN+1], pch = 16, col = "black", ylab = "KL divergence", xlab = "", 
       main = sprintf('%s KSweep KL Divergence', cell_type))
  text(list_lab+0.5, y[list_lab-K.MIN+1]-0.08, labels=list_lab)
}

pdf("k_sweep/results_summary/kl_divergence_plots.pdf", width=12, height=8)
par(mfrow=c(2,3))
lapply(cell.types, kl_plot)
dev.off()


# get.elbow.points.index <- function(x, y) {
#   d1 <- diff(y) / diff(x) # first derivative
#   d2 <- diff(d1) / diff(x[-1]) # second derivative
#   elbow <- d2[which(abs(diff(d2))==max(abs(diff(d2)), na.rm=T))]
#   index <- which(x==x[-c(1:2)][which(d2 == elbow)])
#   return(index)
# }
# 
# k.use=NULL
# for (i in cell.types){
#   h_metrics<- readRDS(file=paste0(i, '/H_metrics.rds'))
#   rank = sapply(strsplit(names(h_metrics), 'R'), function(x){y <- as.numeric(x[2])})
#   df=as.data.frame(cbind(rank, sapply(h_metrics, '[[', 'KL')))
#   colnames(df) <- c('k', 'kl')
#   ind=get.elbow.points.index(df$k,df$kl)
#   plot(df$k, df$kl, pch=19, cex=1, main=i)
#   points(df$k, df$kl_loess, pch=19, cex=1, col="blue")
#   points(df$k[ind], df$kl[ind], pch=19, cex=2,  col="red")
#   k.use[i]=df$k[ind]
# }
#names(k.use)=cell.types
#saveRDS(k.use, file="k.use.rds")


