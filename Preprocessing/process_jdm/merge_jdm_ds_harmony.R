library(cowplot)
library(RColorBrewer)
library(Seurat)
library(harmony)
library(tidyverse)
require(gridExtra)


args=commandArgs(trailingOnly=T)
soupx=as.logical(args[[1]]) 
if (soupx){
  OUT_MERGE_DIR=sprintf("%s/merged_jdm_soupx/", JDM_DIR)
} else {
  OUT_MERGE_DIR=sprintf("%s/merged_jdm/", JDM_DIR)
  jdm0 = readRDS(sprintf("%s/seurat_filt_postdemuxDF_w_meta_filt2.rds", JDM_DIR))
}

dir.create(OUT_MERGE_DIR, showWarnings=T)


# all singlets

MTPATTERN="^MT-"
RIBOPATTERN="^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sobjs = list()
for (idx in 1:6){
  sample_name=sprintf("well%s", idx)
  if (soupx){
    sobj = readRDS(sprintf("%s/soupx/soupx_sobj_adj_well%s.RDS", JDM_DIR, idx))
  } else {
    sobj = subset(jdm0, well==idx)
  }
  sobj@active.assay = 'RNA'
  sobj = NormalizeData(sobj, 
                       assay='RNA', 
                       normalization.method = "LogNormalize",
                       scale.factor = 10000)
  sobj = sobj %>%
    PercentageFeatureSet(pattern = RIBOPATTERN, col.name = "percent.ribo") %>%
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE) 
  
  sobj@meta.data$LIBRARY = sample_name
  sobj@meta.data %>%
    write_csv(sprintf("%s/%s_metadata_temp.csv", OUT_MERGE_DIR, sample_name))
  
  # now add it to the list
  sobjs[[sample_name]] = sobj 
}



first_sobj <- sobjs[[names(sobjs)[1]]]
sobjs[[names(sobjs)[1]]] <- NULL
merged_data <- merge(x = first_sobj, y = unname(sobjs)) 

merged_data <- FindVariableFeatures(merged_data, selection.method = "vst", nfeatures = 3000)
merged_data <- ScaleData(merged_data, vars.to.regress=c('percent.mt', 'percent.ribo', 
                                                        'S.Score', 'G2M.Score'), verbose = FALSE)
merged_data <- RunPCA(merged_data, npcs = 30, verbose = FALSE)
merged_data <- RunUMAP(merged_data, dims = 1:30)

save(merged_data, file=sprintf('%s/merged_temp.RData', OUT_MERGE_DIR))
raw_metadata <- merged_data@meta.data
for (redn in c('pca', 'umap')){
  cell_embeddings <- as.data.frame(merged_data@reductions[[redn]]@cell.embeddings[, c(1: min(5, dim(merged_data@reductions[[redn]]@cell.embeddings)[2]))])
  raw_metadata <- merge(raw_metadata,
                        cell_embeddings,
                        by=0)
  rownames(raw_metadata) <- raw_metadata$Row.names
  raw_metadata$Row.names <- NULL
}

write.table(raw_metadata,
            file=sprintf('%s/merged_metadata_raw.tsv', OUT_MERGE_DIR),
            sep='\t',
            row.names=T,
            col.names=T,
            quote=F)


png(filename=sprintf('%s/merged_raw_library_umap.png', OUT_MERGE_DIR), 
    width = 5, height = 5, units = "in", res = 300)
print(DimPlot(merged_data, group.by='LIBRARY') + NoLegend())
dev.off()

png(filename=sprintf('%s/merged_raw_split_library_umap.png', OUT_MERGE_DIR), 
    width = 15, height = 30, units = "in", res = 300)
print(DimPlot(merged_data, split.by='LIBRARY', ncol=4, raster=F) + 
        theme(legend.position="none", axis.title=element_blank(), 
              axis.text=element_blank()))
dev.off()

pdf(sprintf('%s/merged_harmony_convergence.pdf', OUT_MERGE_DIR))
merged_data <- RunHarmony(merged_data,
                          "LIBRARY",
                          assay.use='RNA',
                          plot_convergence = TRUE,
                          max.iter.harmony=20,
                          max.iter.cluster=30)
dev.off()

merged_data <- RunUMAP(merged_data,
                       dims = 1:30,  # Num PCs to use
                       reduction='harmony',
                       n.neighbors = 30,  # Default. Controls how UMAP balances local (low) versus global (large) structure in the data
                       min.dist = 0.3,   # Default. Controls the size of the clusters. Should be smaller than spread
                       spread = 1,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                       a = NULL,  # Default. Can be used with b instead of using min.dist/spread
                       b = NULL,  # Default. Can be used with a instead of using min.dist/spread
                       verbose = FALSE)

# Calculate the neighborhood graph
merged_data <- FindNeighbors(merged_data,
                             dims = 1:30,  # Num PCs to use
                             reduction='harmony',
                             k.param = 20,  # k for the knn algorithm
                             verbose = FALSE
)

png(filename=sprintf('%s/merged_harmony_library_umap.png', OUT_MERGE_DIR), 
    width = 5, height = 5, units = "in", res = 300)
print(DimPlot(merged_data, group.by='LIBRARY') + NoLegend())
dev.off()

png(filename=sprintf('%s/merged_harmony_split_library_umap.png', OUT_MERGE_DIR), 
    width = 15, height = 30, units = "in", res = 300)
print(DimPlot(merged_data, split.by='LIBRARY', ncol=4, raster=F) + 
        theme(legend.position="none", axis.title=element_blank(), 
              axis.text=element_blank()))
dev.off()

# do it by highlighting cells - easier to viz
libraries = unique(merged_data@meta.data$LIBRARY)
list_cl = sapply(libraries, function(my_library) rownames(merged_data@meta.data %>% 
                                                            filter(LIBRARY==my_library)))
png(filename=sprintf('%s/merged_harmony_split_library_umap_hl.png', OUT_MERGE_DIR), 
    width = 15, height = 30, units = "in", res = 300)
plots = lapply(list_cl, function(.x) DimPlot(merged_data, cells.highlight=.x, combine=T)+NoLegend())
do.call(grid.arrange,  plots)
dev.off()

for (res in c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)){
  if (paste0('louvain_res', res) %in% colnames(merged_data@meta.data)){
    next
  }
  merged_data <- FindClusters(merged_data, verbose = TRUE,
                              algorithm = 1,
                              resolution = res)
  merged_data@meta.data[[paste0('louvain_res', res)]] <- merged_data@meta.data$seurat_clusters
  
  png(filename=sprintf('%s/merged_louvain_res_%s.png', OUT_MERGE_DIR, res), width = 5, height = 5, units = "in", res = 300)
  print(DimPlot(merged_data, group.by=paste0('louvain_res', res), label=T) + NoLegend())
  dev.off()
}

png(filename=sprintf('%s/merged_harmony_split_library_umap_clus.png', OUT_MERGE_DIR), width = 15, height = 12, units = "in", res = 300)
print(DimPlot(merged_data, split.by='LIBRARY', ncol=4) + theme(legend.position="none", axis.title=element_blank(), axis.text=element_blank()))
dev.off()

save(merged_data, file=sprintf('%s/merged_processed.RData', OUT_MERGE_DIR))

metadata <- merged_data@meta.data
for (redn in c('pca', 'umap')){
  cell_embeddings <- as.data.frame(merged_data@reductions[[redn]]@cell.embeddings[, c(1: min(5, dim(merged_data@reductions[[redn]]@cell.embeddings)[2]))])
  metadata <- merge(metadata,
                    cell_embeddings,
                    by=0)
  rownames(metadata) <- metadata$Row.names
  metadata$Row.names <- NULL
}


write.table(metadata,
            file=sprintf('%s/merged_metadata.tsv', OUT_MERGE_DIR),
            sep='\t',
            row.names=T,
            col.names=T,
            quote=F)

### CLUSTERS
Idents(merged_data) = merged_data$RNA_snn_res.1.2
all_markers = FindAllMarkers(merged_data, test.use="MAST", only.pos=T)
save(all_markers, file=sprintf("%s/all_markers_res1.2.RData", OUT_MERGE_DIR))

