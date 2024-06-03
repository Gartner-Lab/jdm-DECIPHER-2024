#### RUN K SWEEP (SERVER) ####

library(rliger)
library(parallel)
library(cluster)
library(liger)
library(stats)
library(RANN)
library(matrixStats)
library(parallel)
library(mixtools)
library(rlist)


##### separate script #####
args=commandArgs(trailingOnly=T)
cell_type=args[1]
K.MIN=2
K.MAX=40

#conda activate R4.0LM
JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm"
setwd(JDM_DIR)

dir.create(sprintf('k_sweep/%s', cell_type))
Liger = readRDS(sprintf('k_sweep/Liger_objects/Liger.%s.rds', cell_type))

# see batch size calculation above
minibatch_size = readRDS('k_sweep/Liger_objects/minibatch_size.rds')

LIGER.rep.function <- function(x){
  # select random seeds
  seed = sample(1:10000, 1)
  Liger <- online_iNMF(Liger, k = x, max.epochs=10, seed = seed, miniBatch_size = 
                         minibatch_size[[cell_type]])
  res = t(Liger@W)
  colnames(res) = paste0("R", x, "_Program", 1:x)
  H_res = list.rbind(Liger@H)
  colnames(H_res) = paste0("R", x, "_Program", 1:x)
  V_res = Liger@V
  V_res = lapply(V_res, function(y){
    rownames(y) = paste0("R", x, "_Program", 1:x)
    y = t(y)
    return(y)
  })
  return(list(W = res, H = H_res, V = V_res, seed = seed))
}

# Depending on biological variation expected across patients, can do a much more restricted sweep depending on the number of programs you input
k = c(K.MIN:K.MAX)
for (j in k){
  reps.k = rep(j, 10)
  names(reps.k) = paste0("rep", 1:length(reps.k))
  Liger_list = mclapply(reps.k, FUN = LIGER.rep.function, mc.cores = 10)
  for (i in 1:length(Liger_list)){
    colnames(Liger_list[[i]][["W"]]) =  paste0("Rep", i, "_", colnames(Liger_list[[i]][["W"]])) 
    colnames(Liger_list[[i]][["H"]]) = paste0("Rep", i, "_", colnames(Liger_list[[i]][["H"]]))
    Liger_list[[i]][["V"]] = lapply(Liger_list[[i]][["V"]], function(y) {
      colnames(y) = paste0("Rep", i, "_", colnames(y))
      return(y)})
  }
  saveRDS(Liger_list, file = 
            sprintf("k_sweep/%s/Liger.%s.list.R%s.rds", cell_type, cell_type, j))
}
rm(Liger_list)
rm(Liger)
gc()

#### CALCULATE CONSENSUS VALUES ####  run replicates and find consensus

# IDENTIFY OUTLIERS
# determine distance threshold - find nearest neighbors on server, then move to local computer to inspect
files = grep(sprintf("Liger.%s.list", cell_type), list.files(sprintf('k_sweep/%s/', cell_type)), 
             value = T) #grabbing file names so that file names are passed to below functions rather than objects themselves
files = files[order(as.numeric(gsub("R", "", unlist(lapply(strsplit(files, '\\.'), `[`, 4)))))]

# find average distance of each program/rep's two nearest neighbors
nn.dist.full = mclapply(files, function(x){
  Liger_list = readRDS(sprintf('k_sweep/%s/%s', cell_type, x))
  W_list = lapply(Liger_list, `[[`, "W")
  W_list = lapply(W_list, function(x) {apply(x, MARGIN = 2, FUN = function(y) {y/norm(y, type = "2")})})
  # L2 normalize (as in cNMF paper)
  W.mat = list.rbind(lapply(W_list, FUN = t))
  nn.dist = rowMeans(nn2(W.mat)$nn.dists[,2:3])
  # the [,2:3] here is the two nearest neighbors
}, mc.cores = 4)
names(nn.dist.full) = unlist(lapply(strsplit(files, '\\.'), `[`, 4))
saveRDS(nn.dist.full, file = sprintf("k_sweep/%s/nn.dist.full.rds", cell_type))
rm(nn.dist.full, files)
gc()

# for some reason this function isn't automatically found unless I call it to the environment
solveNNLS <- function(C, B) {
  .Call('_rliger_solveNNLS', PACKAGE = 'rliger', C, B)
}

# FIND CONSENSUS RESULTS #
dir.create(sprintf('k_sweep/%s/consensus_res', cell_type))
# rank should be over the same sweep as k above
rank = K.MIN:K.MAX #CD4T starts at 5, all others start at 2
names(rank) = paste0("R", rank)

Liger = readRDS(sprintf('k_sweep/Liger_objects/Liger.%s.rds', cell_type))

#cell types to run consensus results section on: CD8T, NK, expanded myeloid, expanded B cells, expanded CD4T cells
# THIS BOTH SAVES THE CONSENSUS RESULTS AND OUTPUTS A SUMMARY OF THE KMEANS RESULTS AS KMEANS_LIST
kmeans_list = mclapply(rank, function(k){
  Liger_list = readRDS(sprintf('k_sweep/%s/Liger.%s.list.R%s.rds', cell_type, cell_type, k))
  W_list = lapply(Liger_list, `[[`, "W")
  W_2norm = lapply(W_list, function(x) {apply(x, MARGIN = 2, FUN = function(y){norm(y, type ="2")})}) #why do we need to store the norm of y?
  W_list = lapply(W_list, function(x) {apply(x, MARGIN = 2, FUN = function(y) {y/norm(y, type = "2")})}) #normalizing over columns instead of rows bc in this liger package, columns are the activity programs (K)
  W.mat = list.rbind(lapply(W_list, FUN = t))
  km.res = kmeans(W.mat, centers = k, nstart = 100, iter.max = 100) #kmeans clustering (1 cluster for each program)
  nn.dist = rowMeans(nn2(W.mat)$nn.dists[,2:3]) #stores 2 nearest neighbors for each replicate of each program at rank k, computes mean distance
  names(nn.dist) = rownames(W.mat) 
  min = diff(quantile(nn.dist, probs = c(0.25, 0.75)))*0.5 + quantile(nn.dist, probs = c(0.75)) 
  if (length(which(nn.dist>min))>0) {
    W.mat.filt = W.mat[-which(nn.dist>min),] 
    km.res.filt = kmeans(W.mat.filt, centers = k, nstart = 100, iter.max = 100) #recalculate kmeans without those outliers
  } else {
    W.mat.filt = W.mat
    km.res.filt = km.res
  } #this filtering step removes nn's that are more than half the iq range from 75th quantile
  table(km.res.filt$cluster)
  W_consensus = matrix(nrow = k, ncol = ncol(W.mat.filt))
  for (i in seq(k)){
    row.ind = which(km.res.filt$cluster==i) 
    if (length(row.ind) > 1){
      W_consensus[i,] = colMedians(W.mat.filt[row.ind,]) #median expression across programs if multiple clusters assigned to 1 program
    } else W_consensus[i,] = W.mat.filt[row.ind,]
  }
  rownames(W_consensus) = paste0("R", k, "_Program", seq(k))
  colnames(W_consensus) = colnames(W.mat.filt)
  V_list = lapply(Liger_list, `[[`, "V") #batch effects 
  for (i in names(V_list)){
    V_list[[i]] = lapply(V_list[[i]], function(x){t(t(x)*(1/W_2norm[[i]]))})
  }
  V.mat = t(list.cbind(lapply(V_list, list.rbind)))
  if (length(which(nn.dist>min))>0) {
    V.mat.filt = V.mat[-which(nn.dist>min),]
  } else {
    V.mat.filt = V.mat
  }
  V_consensus = matrix(nrow = k, ncol = ncol(V.mat.filt))
  for (i in seq(k)){
    row.ind = which(km.res.filt$cluster==i)
    if (length(row.ind) > 1){
      V_consensus[i,] = colMedians(V.mat.filt[row.ind,])
    } else V_consensus[i,] = V.mat.filt[row.ind,]
  }
  rownames(V_consensus) = paste0("R", k, "_Program", seq(k))
  colnames(V_consensus) = colnames(V.mat.filt)
  Batch.info = 1:length(V_list$rep1)
  names(Batch.info) = names(V_list$rep1)
  V_consensus = lapply(Batch.info, function(x){
    V = V_consensus[,((x-1)*length(colnames(W_consensus))+1):(x*length(colnames(W_consensus)))]
  })
  # If no batch effects, can just use predict
  # Liger_predict = online_iNMF(Liger, k = 9, max.epochs=20, miniBatch_size = 1000, projection = T, W.init = t(W_consensus))
  # H_predict = list.rbind(Liger_predict@H)
  # Otherwise, use solveNNLS with W and V initializations
  Batch.info = names(V_list$rep1)
  names(Batch.info) = names(V_list$rep1)
  H_consensus_list = lapply(Batch.info, function(x){
    H = solveNNLS(rbind(t(W_consensus) + t(V_consensus[[x]]), 
                        sqrt(10) * t(V_consensus[[x]])), rbind(t(Liger@scale.data[[x]]), matrix(0, 
                                                                                                dim(W_consensus)[2], dim(Liger@raw.data[[x]])[2])))
  })#input to (sqrt) should be number of replicates of NMF (to account for random seeds)
  H_consensus_list = lapply(H_consensus_list, t)
  H_consensus = list.rbind(H_consensus_list)
  rownames(H_consensus) = rownames((lapply(Liger_list, `[[`, "H"))[[1]])
  colnames(H_consensus) = paste0("R", k, "_Program", seq(k))
  consensus_res = list(H = H_consensus, W = W_consensus, V = V_consensus)
  # This saves the consensus matrices for each rank separately
  saveRDS(consensus_res, file = sprintf("k_sweep/%s/consensus_res/Liger.%s.consensus.R%s.rds", cell_type, cell_type, k))
  kmeans_result = list(km = km.res, km_filt = km.res.filt)
  return(kmeans_result)
}, mc.cores = 10)
saveRDS(kmeans_list, file = sprintf("k_sweep/%s/kmeans_list.rds", cell_type))
rm(kmeans_list)
gc()


# METRICS FOR CHOOSING K - all calculated on consensus results (generated in chunk above)#### 


# 2. H_metrics: KL divergence/JSD of cell loadings from uniform (as in original LIGER paper)  ####
# create new kl divergence function that works on the consensus H matrix instead of needing the full Liger object
median_kl_divergence_uniform = function(Hs){
  n_cells = lapply(Hs, nrow)
  n_factors = ncol(Hs)
  dataset_list = list()
  scaled = scale(Hs, center=FALSE, scale=TRUE)
  inflated = t(apply(scaled, 1, function(x) {
    replace(x, x == 0, 1e-20)
  }))
  inflated = inflated/rowSums(inflated)
  divs = apply(inflated, 1, function(x) {log2(n_factors) + sum(log2(x) * x)})
  res = median(divs)
  return(res)
}

H_metrics = lapply(rank, function(k){
  H = readRDS(sprintf('k_sweep/%s/consensus_res/Liger.%s.consensus.R%s.rds',cell_type, cell_type, k))$H 
  KL = median_kl_divergence_uniform(H)
  H_list = split(H, seq(nrow(H)))
  #ignoring JSD for rank optimization; just use KL divergence
  # JSD uses the KL divergence to calculate a normalized score that is symmetrical (and ranges from 0 to 1)
  #JSD = median(unlist(mclapply(H_list, function(x) {
  #JSD(as.data.frame(rbind(x, rep(1, length(x)))), est.prob = "empirical")
  #3}, mc.cores = 4)))
  res = list(KL = KL)
})
#JSD_max = lapply(rank, function(x){
#JSD(rbind(c(rep(0, x-1), 1), rep(1, x)), est.prob = "empirical")
#})

saveRDS(H_metrics, file = sprintf('k_sweep/%s/H_metrics.rds', cell_type)) #H_metrics didn't fully finish for CD8T


# 1. error of result - run for each cell type####
Liger.data = readRDS(sprintf('k_sweep/Liger_objects/Liger.%s.rds', cell_type))@scale.data
error_list = mclapply(rank, function(k){ #error is optional metric after looking at KL divergence
  consensus = readRDS(sprintf('k_sweep/%s/consensus_res/Liger.%s.consensus.R%s.rds', cell_type, cell_type, k))
  data = list()
  for (i in names(consensus$V)){
    data[[i]] = consensus$H[rownames(Liger.data[[i]]),] %*% (consensus$W + consensus$V[[i]])
  }
  data = list.rbind(data)
  Liger.data = list.rbind(Liger.data)
  error = abs(Liger.data - data)
  error_list[[paste0("R", k)]] = norm(error, type = "F")
}, mc.cores = 10)

saveRDS(error_list, file=sprintf('k_sweep/%s/consensus_res/error_list.rds', cell_type))
rm(error_list)
gc()



# 3. stability across replicates (JSD from uniform distribution) #### -->skip this step (metric #3)
#library(philentropy)
#library(diceR)
# kmeans_list = list()
# replicate_stability_metrics = list()
# for (i in cell.types){
#   kmeans_list[[i]] = readRDS(paste0('k_sweep/', i, '/kmeans_list.rds'))
#   kmeans_clusters = lapply(lapply(kmeans_list[[i]], '[[', "km"), '[[', 'cluster')
#   JSD_res = lapply(kmeans_clusters, function(x){
#     JSD = JSD(rbind(table(x), rep(5, length(table(x)))), est.prob = 'empirical')
#     # maximum JSD (i.e. worst performance) will depend on the rank, normalize so we can compare across ranks
#     k_length = length(table(x))
#     JSD_max = JSD(rbind(c(rep(0, k_length-1), 5*k_length), rep(5, k_length)), est.prob = 'empirical')
#     JSD_norm =  JSD/JSD_max
#     res = list(JSD = JSD, JSD_norm = JSD_norm)
#     return(res)
#   })
#   replicate_stability_metrics[[i]] = list(clusters = kmeans_clusters, JSD = lapply(JSD_res, '[[', "JSD"),
#                                           JSD_norm = lapply(JSD_res, '[[', "JSD_norm"))
# }
# 
# saveRDS(replicate_stability_metrics, file = "k_sweep/results_summary/replicate_stability_metrics.rds")

