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
  