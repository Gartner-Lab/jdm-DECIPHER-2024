library(tidyverse)
library(Matrix)

read_format_matrix = function(my_path){
  mtx = readMM(sprintf("%s/matrix.mtx.gz", my_path))
  features = read_tsv(sprintf("%s/features.tsv.gz", my_path), col_names=F)
  cells = read_tsv(sprintf("%s/barcodes.tsv.gz", my_path), col_names=F)
  rownames(mtx) = features$X1
  colnames(mtx) = cells$X1
  return(mtx)
}

get_correct_ab_names = function(mtx){
  # fix the formatting
  substr_row = unname(sapply(rownames(mtx), function(x) substr(x, 1,2)))
  ab_names = rownames(mtx)[substr_row!="EN"]
  ab_names2 = ab_names %>% 
    str_replace_all( "^[A-z\\_]+;", "") %>% 
    str_replace_all("^;", "") %>%
    str_replace_all(";;1$", "") %>% 
    str_replace_all( ";1$", "") %>%
    str_replace_all(";", "--") %>%
    str_replace_all("\\_", "\\-") %>%
    str_replace_all("\\.1$", "") %>%
    str_replace_all("----2", "\\.2")
  orig_ab_names = rownames(ADT_counts) %>%
    str_replace_all("\\.1$", "") %>%
    str_replace_all("\\_", "\\-")
  
  ## isotype control duplicates
  orig_ab_names[duplicated(orig_ab_names)]=unname(sapply(orig_ab_names[duplicated(orig_ab_names)], function(x) paste(x, 1, sep=".")))
  ab_names2[duplicated(ab_names2)]=unname(sapply(ab_names2[duplicated(ab_names2)], function(x) paste(x, 1, sep=".")))
  stopifnot(length(setdiff(ab_names2, orig_ab_names)) ==0)
  stopifnot(length(setdiff(orig_ab_names, ab_names2)) == 0)
  
  return(list("mtx_names"=ab_names2, "sobj_names"=orig_ab_names))
}
