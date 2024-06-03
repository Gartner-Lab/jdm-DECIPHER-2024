library(Seurat)
library(tidyverse)
library(dittoSeq)


args = commandArgs(trailingOnly=T)
norm.method = args[1]
JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm"


OUT_DIR = sprintf("%s/rpca_v2/", JDM_DIR)
dir.create(OUT_DIR, showWarnings=T)

jdm0 = readRDS(sprintf("%s/seurat_filt_postdemuxDF_w_meta_filt2.rds", JDM_DIR))
jdm0 = AddMetaData(jdm0, sapply(jdm0$well, function(x) paste0("well", x)), "well_name")
DefaultAssay(jdm0) = "ADT"

if (norm.method == "CLR"){
  sobjs.list <- SplitObject(jdm0, split.by = "well_name")
  sobjs.list = lapply(sobjs.list, function(sobj) NormalizeData(sobj,
                                                               normalization.method = 'CLR'))
}
if (norm.method == "DSB"){
  sobjs.list = list()
  for (idx in 1:6){
    sample_name=sprintf("well%s", idx)
    sobj = subset(jdm0, well==idx)
    adt_dat = readRDS(sprintf("%s/dsb_v2/adt_dsb_%s.RDS", JDM_DIR, idx))
    sobj=subset(sobj, cells=colnames(adt_dat))
   # stopifnot(all(colnames(adt_dat)==colnames(sobj)))
    DefaultAssay(sobj) = "RNA"
    sobj[["ADT"]] = NULL
    sobj[["ADT"]] = CreateAssayObject(data=adt_dat)
    DefaultAssay(sobj) = "ADT"
    # now add it to the list
    sobjs.list[[sample_name]] = sobj 
  }
}

list_wells = unique(jdm0@meta.data$well_name)
names(sobjs.list) = list_wells


#Using all non-isotype markers for the following analyses.
features <- grep("isotype", rownames(sobjs.list[[1]][["ADT"]]), ignore.case = T, value = T, invert = T)
# Remove the isotype controls
features <- features[ features %in% grep("isotype", 
                                         rownames(sobjs.list[[1]][["ADT"]]), ignore.case = T, value = T, invert = T) ]
print(features)

sobjs.list <- lapply(list_wells, FUN = function(SAMPLE_NAME) {
  x = sobjs.list[[SAMPLE_NAME]]
  x <- ScaleData(x, features = features, verbose = TRUE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  sobj_before_file = sprintf("%s/%s_ADT_%s.png", OUT_DIR, SAMPLE_NAME, norm.method)
  if (!file.exists(sobj_before_file)){
    x <- RunUMAP(x, dims=1:18)
    DimPlot(x)+NoLegend()
    ggsave(sobj_before_file, height=5, width=5)
  }
  return(x)
  
})
names(sobjs.list) = list_wells

# RPCA
immune.anchors <- FindIntegrationAnchors(object.list = sobjs.list, 
                                         anchor.features = features, 
                                         reference = 1, 
                                         reduction = "rpca" )

integ_data <- IntegrateData(anchorset = immune.anchors, new.assay.name = "integrated.ADT")
DefaultAssay(integ_data) <- "integrated.ADT"
integ_data <- ScaleData(integ_data, vars.to.regress=c('nCount_ADT', 'nFeature_ADT', 
                                                      'S.Score', 'G2M.Score'),  verbose = FALSE)
integ_data <- RunPCA(integ_data, npcs = 30, verbose = FALSE)
integ_data <- RunUMAP(integ_data, reduction = "pca", dims = 1:30)
DimPlot(integ_data, split.by="orig.ident")+NoLegend()
ggsave(sprintf("%s/post_%s.pdf", OUT_DIR, norm.method), height=10, width=15)


saveRDS(integ_data, file = sprintf("%s/integ_ADT_%s.rds", OUT_DIR, norm.method))