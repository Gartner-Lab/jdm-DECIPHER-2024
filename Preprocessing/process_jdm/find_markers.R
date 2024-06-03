library(Seurat) 
JDM_DIR="/krummellab/data1/erflynn/premier/jdm/data/jdm"
OUT_DIR=sprintf("%s/wnn_v2/reclus", JDM_DIR)

load(sprintf("%s/wnn_reclus.RData", OUT_DIR))
marker.res=1.4
Idents(sobj) = sprintf("wsnn_leiden_reclus_res.%s", marker.res)

#print("finding markers")
#DefaultAssay(sobj) = "integrated.ADT"
#wnn_markers_adt = FindAllMarkers(sobj, only.pos=T)
#save(wnn_markers_adt, file=sprintf("%s/wnn_adt_markers_%s.RData", OUT_DIR, marker.res))

DefaultAssay(sobj) = "RNA"
wnn_markers = FindAllMarkers(sobj, test.use="negbinom", only.pos=T, min.pct=0.3)
save(wnn_markers, file=sprintf("%s/wnn_rna_markers_%s.RData", OUT_DIR, marker.res))
