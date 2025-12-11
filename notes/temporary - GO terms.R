library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(GO.db)

ids = c("GO:0042092","GO:0002828","GO:0002829","GO:0002830",
        "GO:0045064","GO:0002297","GO:0045628","GO:0045629",
        "GO:0045630","GO:0035712","GO:2000569","GO:2000570",
        "GO:0035745","GO:2000551","GO:2000552","GO:2000553",
        "GO:0035708","GO:2000571","GO:2000572","GO:0048295",
        "GO:0005136","GO:0005137","GO:0005144","GO:0004913",
        "GO:0004914","GO:0016515","GO:0019767","GO:0019863",
        "GO:0016516","GO:0005895","GO:0005898","GO:0070450",
        "GO:0032998")

genes_for_go <- select(org.Hs.eg.db,
                       keys = ids,
                       columns = c("ENTREZID", "SYMBOL", "GENENAME"),
                       keytype = "GO")

go_term <- select(GO.db,
                  keys = ids,
                  columns = "TERM",
                  keytype = "GOID")

go = merge(go_term, genes_for_go, by.x = "GOID", by.y = "GO")
go_list = split(go$SYMBOL, go$GOID)

go_term = genes_for_go %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  group_by(GO) %>%
  summarise(nGenes = n(),
            Genes = paste0(SYMBOL,collapse = ", ")) %>%
  merge(go_term, ., by.x = "GOID", by.y = "GO", all = T)
