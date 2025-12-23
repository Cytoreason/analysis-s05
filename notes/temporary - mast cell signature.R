library(cytoreason.deconvolution)

need = c(mast = "CRCL_0000009", ilc2 = "CRCL_0000163", ilc3 = "CRCL_0000164", basophil = "CRCL_0000357", eosinophil = "CRCL_0000358")

gut <- SignatureCollection("gut_v9")
  gut_markers = get_markers(gut)
  gut_markers = split(gut_markers$SYMBOL, gut_markers$celltype)
lung <- SignatureCollection("lung_v4")
  lung_markers = get_markers(lung)
  lung_markers = split(lung_markers$SYMBOL, lung_markers$celltype)
skin <- SignatureCollection("skin_v12")
  skin_markers = get_markers(skin)
  skin_markers$cellName = cytoreason.datasource::translate_ids2names(skin_markers$celltype)
  skin_markers = split(skin_markers$SYMBOL, skin_markers$celltype)
skeletal <- SignatureCollection("skeletal_v2")
  skeletal_markers = get_markers(skeletal)
  skeletal_markers = split(skeletal_markers$SYMBOL, skeletal_markers$celltype)
synovium <- SignatureCollection("synovium_v3_hybrid")
  synovium_markers = get_markers(synovium)
  synovium_markers = split(synovium_markers$SYMBOL, synovium_markers$celltype)
  
  
overlap_gut_skin = sapply(skin_markers, function(x) sapply(gut_markers[which(names(gut_markers) %in% need)], function(y) length(intersect(x, y))))
  rownames(overlap_gut_skin) = names(need)[match(rownames(overlap_gut_skin), need)]
  overlap_gut_skin = t(overlap_gut_skin)
  rownames(overlap_gut_skin) = skin_mapping$celltype[match(rownames(overlap_gut_skin), skin_mapping$id)]
  

overlap_lung_skin = sapply(skin_markers, function(x) sapply(lung_markers[which(names(lung_markers) %in% need)], function(y) length(intersect(x, y))))
  rownames(overlap_lung_skin) = names(need)[match(rownames(overlap_lung_skin), need)]
  overlap_lung_skin = t(overlap_lung_skin)
  rownames(overlap_lung_skin) = skin_mapping$celltype[match(rownames(overlap_lung_skin), skin_mapping$id)]
  
## Using information from Ben to choose the right signature
cell_mapping = read.csv("~/exportedFiles/cells_mapping.csv")

lung <- SignatureCollection("lung_v4")
  lung = get_markers(lung)
  lung = split(lung$SYMBOL, lung$celltype)
  lung = lung$CRCL_0000009

ileum <- SignatureCollection("ileum_v2", type = "cyto_marker")
  ileum = get_markers(ileum)
  ileum = split(ileum$SYMBOL, ileum$celltype)
  ileum = ileum$CRCL_0000009
  
pan_cancer <- SignatureCollection("pan_cancer_v1")
  pan_cancer = get_markers(pan_cancer)
  pan_cancer = split(pan_cancer$SYMBOL, pan_cancer$celltype)
  pan_cancer = pan_cancer$CRCL_0000009
  
gut <- SignatureCollection("gut_v9")
  gut = get_markers(gut)
  gut = split(gut$SYMBOL, gut$celltype)
  gut = gut$CRCL_0000009
  
myocardium <- SignatureCollection("myocardium_v1")
  myocardium = get_markers(myocardium)
  myocardium = split(myocardium$SYMBOL, myocardium$celltype)
  myocardium = myocardium$CRCL_0000009
  
all(gut == ileum)

allSigs = list(lung, ileum, pan_cancer, myocardium)
allSigs <- do.call(cbind, lapply(allSigs, function(v) c(v, rep(NA, max(lengths(allSigs)) - length(v)))))
write.csv(allSigs,"~/exportedFiles/mast_cell_signatures.csv")
