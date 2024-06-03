library(dsb)
library(Seurat)
library(Matrix)
library(tidyverse)

JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm/"
CODE_DIR="/krummellab/data1/erflynn/premier/jdm/"
source(sprintf("%s/scripts/utils/dsb_utils.R", CODE_DIR))


# usage "run_dsb_jdm_v2.R $WELL"
args = commandArgs(trailingOnly=T)
WELL= as.numeric(args[1])



nfeature.adt.filter = TRUE
nfeature.adt.min = 70
isotype.ctl.filter = TRUE
isotype.ctl.max = 50
background.filter = TRUE


jdm_obj = readRDS(sprintf("%s/seurat_filt_postdemuxDF_w_meta_filt2.rds", JDM_DIR))
sobj = subset(jdm_obj, well==WELL) # grab only the well we want
rm(jdm_obj); gc()

DefaultAssay(sobj) = 'RNA'

ADT_counts=sobj@assays$ADT@counts
isotype_ctls = rownames(ADT_counts)[str_detect(rownames(ADT_counts), "Iso")]

isotype_ctl_data = as.matrix(ADT_counts[isotype_ctls,])
for(i in 1:length(rownames(isotype_ctl_data))){
  isotype_name = rownames(isotype_ctl_data)[[i]]
  sobj=AddMetaData(sobj, 
                   unlist(isotype_ctl_data[isotype_name,]), 
                   sprintf("%s",str_replace_all(isotype_name, "-", "_")))
}
sobj = AddMetaData(sobj, apply(isotype_ctl_data, 2, max), "isotype_ctl_max")

#### ---- APPLY PARAMETERS ---- ####
## 1) nfeature_ADT filter
if (nfeature.adt.filter){
  print(sprintf("FILTER: nFeatureADT > %s: Removing %s cells out of %s", 
                nfeature.adt.min,
                sobj@meta.data %>% 
                  filter(nFeature_ADT <= nfeature.adt.min) %>% 
                  nrow(),
                nrow(sobj@meta.data)))
  sobj = subset(sobj, nFeature_ADT > nfeature.adt.min)
}


if (isotype.ctl.filter){
  print(sprintf("FILTER: isotype control max < %s: Removing %s cells out of %s", 
                isotype.ctl.max,
                sobj@meta.data %>% 
                  filter(isotype_ctl_max >= isotype.ctl.max) %>% 
                  nrow(),
                nrow(sobj@meta.data)))
  sobj = subset(sobj, isotype_ctl_max < isotype.ctl.max)
}

ADT_counts = sobj[['ADT']]@counts 


### DSB ###
my_path = sprintf("%s/raw/crosslong_raw/jdm_crosslong%s/outs/raw_feature_bc_matrix/", JDM_DIR, WELL)
mtx = read_format_matrix(my_path)
colnames(mtx) = str_replace_all(colnames(mtx), "-1", sprintf("-%s", WELL))
stopifnot(length(setdiff(colnames(ADT_counts), colnames(mtx)))==0)
cor_names = get_correct_ab_names(mtx)


## 3a) DSB try different empty droplet filter
if (background.filter){
  substr_row = unname(sapply(rownames(mtx), function(x) substr(x, 1,2)))
  prot_mtx = mtx[substr_row!="EN",]
  gene_mtx = mtx[substr_row=="EN",]
  md = data.frame(
    prot.size = log10(colSums(prot_mtx)), # -Inf from 0 ones
    rna.size = log10(colSums(gene_mtx)) # -Inf from 0 ones
  )
  nrow(md)
  
  md = md[md$rna.size > 0 & md$prot.size > 0, ]
  nrow(md)
  
  background_drops = rownames(
    md[ md$prot.size > 1.5 & 
          md$prot.size < 3 & 
          md$rna.size < 2.5, ]
  ) 
  print(sprintf("FILTER background with %s empty drops",
                length(background_drops)))
  empty_ab_mtx = as.matrix(prot_mtx[ , background_drops])
  
} else {
  empty_droplets = colnames(mtx)[colSums(mtx) < 100] # TODO should we define these differently
  substr_row = unname(sapply(rownames(mtx), function(x) substr(x, 1,2)))
  empty_ab_mtx = mtx[substr_row!="EN", empty_droplets]
  print(sprintf("%s empty droplets", length(empty_droplets)))
}

# rename rows and subset for the background
rownames(empty_ab_mtx) = cor_names$mtx_names
ADT_background = empty_ab_mtx[cor_names$sobj_names,]
rownames(ADT_counts) = cor_names$sobj_names

# make sure these weren't included in the actual data (they should all be good)
stopifnot(length(intersect(colnames(ADT_background), colnames(ADT_counts)))==0)
stopifnot(rownames(ADT_counts)==rownames(ADT_background))


adt_norm_out = DSBNormalizeProtein(
  # cell-containing droplet raw protein count matrix
  cell_protein_matrix = ADT_counts, 
  # empty/background droplet raw protein counts
  empty_drop_matrix = ADT_background,
  # recommended step: model + remove the technical component of each cell's protein library
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  quantile.clipping = TRUE, # adding this in to clip extreme outliers
  quantile.clip = c(0.01, 0.99), # slightly more stringent than the default to clip large outliers
  return.stats = TRUE,
  isotype.control.name.vec = isotype_ctls
)
adt_norm = adt_norm_out$dsb_normalized_matrix
print(quantile(adt_norm, c(0, 0.001, 0.01, 0.1, 0.9, 0.99, 0.999, 1)))

non_isotype = setdiff(rownames(adt_norm), isotype_ctls)
adt_norm2=adt_norm[non_isotype,]

saveRDS(adt_norm2, file=sprintf("%s/dsb_v2/adt_dsb_%s.RDS", JDM_DIR, WELL))
