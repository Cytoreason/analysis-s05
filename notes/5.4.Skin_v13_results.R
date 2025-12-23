devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)

ccm_v13 = as_ccm_fit("wf-d4626ee8f7")
ccm_original = as_ccm_fit("wf-08a6a0a503")

skin = read.csv("~/data/skin_modified.csv", row.names = NULL, sep = "\t", header = F)
skin = rbind(skin, c("basophil","CRCL_0000357"))

signatureMapping = readRDS(get_workflow_outputs("wf-0cb82886cc"))

## Extract ct_test
## ============================
ct_test = rbind(statistic_table(ccm_v13$meta$AD$ct_test),
                statistic_table(ccm_v13$meta$L_vs_HC$ct_test),
                statistic_table(ccm_v13$meta$L_vs_NL$ct_test),
                statistic_table(ccm_v13$meta$NL_vs_HC$ct_test))
ct_test$feature_id = skin$V1[match(ct_test$feature_id, skin$V2)]
ct_test$skin_signature = "v13"
uploadToBQ(ct_test, bqdataset = "s05_atopic_dermatitis", tableName = "skinV13_ct_test")

ct_test_previous = rbind(statistic_table(ccm_original$meta$AD$ct_test),
                         statistic_table(ccm_original$meta$L_vs_HC$ct_test),
                         statistic_table(ccm_original$meta$L_vs_NL$ct_test),
                         statistic_table(ccm_original$meta$NL_vs_HC$ct_test))
ct_test_previous$skin_signature = "v12"


ct_tests = rbind(ct_test_previous, ct_test)
or = ct_test_previous %>% 
  dplyr::filter(term == "L_vs_NL") %>%
  arrange(log10(fdr) * sign(estimate)) %>%
  dplyr::select(feature_id) %>%
  rbind("basophils","innate lymphoid cell 2","mast cell")

ggplot(ct_tests, aes(x = -log10(fdr) * sign(estimate), y = factor(feature_id, rev(or$feature_id)))) +
  geom_col() +
  facet_grid(cols = vars(term), rows = vars(skin_signature)) +
  theme_minimal() +
  ggpubr::border() +
  labs(y = "Cell")
ggsave("~/analysis-s05/figures/skin_v13/ct_tests.png", height = 2000, units = "px", bg = "white")


## Extract ct_test per dataset
## =========================================
ct_tests = lapply(names(ccm_v13$datasets), function(d){
  lapply(names(ccm_v13$datasets[[d]]$model), function(model){
    statistic_table(analysisResultElement(ccm_v13$datasets[[d]]$model[[model]],"ct_test"))
  }) %>% bind_rows()
}) %>% bind_rows()

ct_tests$feature_id = skin$V1[match(ct_tests$feature_id, skin$V2)]
ct_tests = ct_tests[which(ct_tests$term %in% c("DZ_vs_HC","L_vs_HC","L_vs_NL","NL_vs_HC")),]

ggplot(ct_tests[which(ct_tests$feature_id %in% c("basophil","innate lymphoid cell 2","mast cell")),], aes(x = -log10(fdr) * sign(estimate), y = experiment_id)) +
  geom_col(aes(fill = ifelse(-log10(fdr) * sign(estimate) > 0, "#58508d", "#ffa600"))) +
  facet_grid(cols = vars(model), rows = vars(feature_id), scales = "free_x") +
  scale_fill_identity()+
  theme_minimal() +
  ggpubr::border() +
  labs(y = "Cell")



## extract geneset-cell correlations
## ======================================
geneset_cell_corr = statistic_table(ccm_v13$meta$L_vs_NL$feature_cell_correlations, term = "L")
  geneset_cell_corr = geneset_cell_corr[which(geneset_cell_corr$feature_type1 == "gene_set"),]
  geneset_cell_corr$feature_id_2 = skin$V1[match(geneset_cell_corr$feature_id_2, skin$V2)]
  geneset_cell_corr = geneset_cell_corr[,c(15,12,13,30,11,5,20)]
  colnames(geneset_cell_corr)[3] = "cell"
  colnames(geneset_cell_corr)[5] = "correlation"
  geneset_cell_corr$log10_fdr = -log10(geneset_cell_corr$fdr)

idx = which(geneset_cell_corr$collection == "X2")
geneset_cell_corr$collection[idx] = signatureMapping$collection[match(geneset_cell_corr$feature_id_1[idx], signatureMapping$previousID)]
geneset_cell_corr$feature_id_1[idx] = signatureMapping$signature[match(geneset_cell_corr$feature_id_1[idx], signatureMapping$previousID)]

# Summarize random and shuffled random
for(id in c("random","top50","bottom50")){
  idx = which(str_detect(geneset_cell_corr$feature_id_1,id))
  if(id != "random") { id = paste0("smoothedRandom_",id)}
  
  tmp = geneset_cell_corr[idx,] %>%
    dplyr::group_by(submodel, cell, collection, signature_collection) %>%
    summarise(fdr = mean(fdr),
              correlation = mean(correlation)) %>%
    mutate(log10_fdr = -log10(fdr)) %>%
    mutate(feature_id_1 = id) %>%
    dplyr::select(colnames(geneset_cell_corr))
  
  geneset_cell_corr = geneset_cell_corr[-idx,]
  geneset_cell_corr = rbind(geneset_cell_corr, tmp)
  rm(idx, tmp)
}

uploadToBQ(geneset_cell_corr, bqdataset = "s05_atopic_dermatitis", tableName = "skinV13_geneset_cell_corr")

x2 = geneset_cell_corr[which(geneset_cell_corr$collection == "X2"),]
x2$feature_id_1 = signatureMapping$New_identifier[match(x2$feature_id_1, signatureMapping$signature)]
x2 = x2[which(x2$feature_id_1 %in% c("aIgE late activation","X2 early activation refined","X2 early general inhibition","X2 early general inhibition refined",
                                     "X2 late inhibition","X2 late inhibition refined","X2 early activated inhibition","X2 early activated inhibition refined")),]
x2 = rbind(x2, geneset_cell_corr[which(geneset_cell_corr$collection %in% c("th2","Mast","Neuronal","neuroinflammation","epidermis")),])

ggplot(x2[which(x2$submodel == "bulk"),], aes(x = correlation, y = log10_fdr, color = collection)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
  facet_wrap(~cell, scales = "free", ncol = 6) +
  theme_minimal() +
  ggpubr::border()
ggsave("~/analysis-s05/figures/skin_v13/geneset_cell_correlation.png", width = 5000, height = 2500, units = "px", bg = "white")


## extract ct_test post Dupilumab
## ======================================
ct_test_post_dupi = statistic_table(analysisResultElement(ccm_v13$datasets$GSE130588__GPL570$model$Treatment__GSE130588__GPL570__Dupilumab,"ct_test"),
                                    term = c("W16_vs_W0:DupilumabL","W16_vs_W0:DupilumabNL","W16_vs_W0:PlaceboL",
                                             "W4_vs_W0:DupilumabL","W4_vs_W0:PlaceboL"))
ct_test_post_dupi$feature_id = skin$V1[match(ct_test_post_dupi$feature_id, skin$V2)]
ct_test_post_dupi = ct_test_post_dupi[,c("experiment_id", "term", "feature_id", "signature_collection", "estimate", "fdr")]
ct_test_post_dupi$log10_fdr = -log10(ct_test_post_dupi$fdr)
colnames(ct_test_post_dupi)[3] = "cell"
uploadToBQ(geneset_cell_corr, bqdataset = "s05_atopic_dermatitis", tableName = "skinV13_ct_test_post_dupi")

ggplot(ct_test_post_dupi, aes(x = log10_fdr*sign(estimate), y = reorder(cell, -log10_fdr*sign(estimate)))) +
  geom_col() +
  facet_wrap(~term, scales = "free")
ggsave("~/analysis-s05/figures/skin_v13/ct_test_post_dupi.png", width = 5000, height = 2500, units = "px", bg = "white")


## cell loadings
## ======================================
cellLoadings = ccm_v13$multi$meta_pca$`1`$v[,1:3]
cellLoadings[,1] = (-1) * cellLoadings[,1]
cellLoadings = reshape2::melt(cellLoadings)
colnames(cellLoadings) = c("Cell","PC","Loading")
cellLoadings$PC = paste0("PC",cellLoadings$PC)
cellLoadings$Cell = skin$V1[match(cellLoadings$Cell, skin$V2)]
uploadToBQ(geneset_cell_corr, bqdataset = "s05_atopic_dermatitis", tableName = "skinV13_cellLoadings")

ggplot(cellLoadings, aes(x = Loading, y = Cell)) +
  geom_col() +
  facet_wrap(~PC, scales = "free")
ggsave("~/analysis-s05/figures/skin_v13/cellLoadings.png", bg = "white")


## cell loadings per dataset
## ======================================
cellLoadings = lapply(names(ccm_v13$datasets), function(d) {return(data.frame(dataset = d, ccm_v13$datasets[[d]]$cell_pca$rotation[,1:2]))})
cellLoadings = lapply(cellLoadings, function(x) {
  x$Cell = skin$V1[match(rownames(x),skin$V2)]
  x = reshape2::melt(x)
  return(x)
}) %>% do.call(rbind,.)
colnames(cellLoadings)[3:4] = c("PC","Loading")
cellLoadings$dataset = sapply(cellLoadings$dataset, function(x) strsplit(x,"__")[[1]][1])

pc1 = cellLoadings[which(cellLoadings$PC == "PC1"),]
pc1$Cell_Dataset = cytoreason.gx::reorder_within(pc1$dataset, pc1$Loading, pc1$Cell)

ggplot(pc1, aes(x = Loading, y = Cell_Dataset)) +
  geom_col() +
  facet_wrap(~Cell, scales="free") +
  cytoreason.gx::scale_y_reordered()+
  theme_minimal() +
  ggpubr::border() +
  theme(axis.text.y = element_text(size = 6))
ggsave("~/analysis-s05/figures/skin_v13/cellLoadings_perDataset.png", height = 4000, width = 3000, units = "px", bg = "white")

