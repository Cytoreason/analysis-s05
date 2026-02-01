devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)

# In this part, we use hybrid skin signatures that includes mast cells, basophils and ILC2 from different tissues
# in order to mitigate some concerns by the client that we are not capturing type 2 immunity
config = read.csv(get_task_inputs("wf-08a6a0a503","0", files_names_grepl_pattern = "ccm-metadata.csv"))
config$asset_id[which(config$experiment_id == "GSE153007")] <- 'wf-2ac040bba1:0:GSE153007.RDS'
config$asset_id[which(config$experiment_id == "GSE60709")] <- 'wf-923925db68:0:GSE60709.RDS'
config$ccm_annotate = "feature"

sample_ann = readRDS(get_workflow_outputs("wf-558ccafbfe"))
pheno_vars = readRDS(get_workflow_outputs("wf-4c71fb2030"))

ccm_stage_prepare_dataset_collection(config, annotate = list(sample_annotation_table = sample_ann))

IMAGE <- ccm_cyto_cc("master_1.0.1",
                     dataset.cell_contribution.service = "eu.gcr.io/cytoreason/ci-cytoreason.deconvolution-package:SUP_6578_add_skin_v14_v15_v16_2.2.8",
                     cell_specific_differences = "eu.gcr.io/cytoreason/ci-cytoreason.deconvolution-package:SUP_6578_add_skin_v14_v15_v16_2.2.8")

tissue_signature = "skin_v13:SUP-6496-skin_v13"
tissue_signature = "skin_v14:SUP-6578-add-skin_v14-v15-v16"


## First run: integrate the skin signatures into the disease model
## ======================================================================================
ccm_api_generate_ccm(config,
                     qc = FALSE,
                     term_metadata = FALSE, # fails for some reason
                     prepare_data = list(annotate = list(sample_annotation_table = sample_ann)),
                     dataset = list(cell_contribution = list(signature_collection = tissue_signature), 
                                    services = .skip("feature_pca")),
                     dataset2 = list(gene_set_activity = list(collection = c("h", "kegg", "reactome", "btm"))),
                     model = list(services = .skip("gene_set_activity_differences","feature_geneset_correlations")),
                     image = IMAGE,
                     tags = list(tissue="skin_v14", condition="AD", project="evo", analysis="AD_skin_v14"))
# generate_ccm -- Fri Dec 19 10:48:35 2025: wf-e82c5ab874 []
# generate_ccm -- Fri Dec 19 15:06:38 2025: wf-5f05a83bc6 [] - tag with skin_v13
# generate_ccm -- Fri Dec 19 15:06:51 2025: wf-dd8c8843e6 [] - tag with skin_v13:SUP-6496-skin_v13
# generate_ccm -- Fri Dec 19 16:14:48 2025: wf-c82f4f91e9 [] - change deconvolution image
# generate_ccm -- Thu Dec 25 12:57:59 2025: wf-6f00628526 [] - skin_v14
# generate_ccm -- Thu Dec 25 12:57:12 2025: wf-52a7393eb5 [] - skin_v15
# generate_ccm -- Thu Dec 25 12:57:28 2025: wf-6df6f6bdca [] - skin_v16


## Second run: analysis of the custom gene sets with the new cell signatures
## ======================================================================================
signatures = readRDS(get_workflow_outputs("wf-d22894f90b"))[-3]
th2 = list(th2 = readRDS(get_workflow_outputs("wf-410536ebd3")),
           neuroinflammation = readRDS(get_workflow_outputs("wf-4ce41a599a")),
           epidermis = readRDS(get_workflow_outputs("wf-d3b177cd82")))

# reorder
th2$epidermis = c(th2$epidermis, th2$th2[c("lichenification","Terminal_Differentiation_and_Lipids")])
th2$th2 = th2$th2[-which(names(th2$th2) %in% c("lichenification","Terminal_Differentiation_and_Lipids","Th1_Related","Th17_Related","TH22_IL22_Related"))]

allSignatures = c(signatures, th2)
allSignatures = lapply(allSignatures, function(x) lapply(x, unique))
rm(th2, signatures)

gene_set_limits <- setNames(replicate(length(allSignatures), c(1L, Inf), simplify = FALSE), names(allSignatures)) # Define gene set size per signature

ccm_api_run_custom_gene_set_analysis(cytoreason.assets::AssetData("wf-6df6f6bdca") , # when using AssetData it will pull relevant tags automatically, and makes the relationship traceable
                                     custom_gene_set_collections = allSignatures,
                                     prepare_data = list(annotate = list(sample_annotation_table = sample_ann)),
                                     dataset = list(cell_contribution = list(signature_collection = tissue_signature),
                                                    services = .skip("feature_pca")),
                                     dataset2 = list(gene_set_activity = list(collection = c("h", "kegg", "reactome", "btm"))),
                                     model = list(pheno_feature_correlations = list(phenotypic_variables = pheno_vars)),
                                     meta = list(gx_diff = list(collection_size_limits = gene_set_limits)),
                                     submodel = c("bulk","exprs_adjusted__1__1"),
                                     workflow_overrides=list(activeDeadlineSeconds=172800),
                                     image = IMAGE,
                                     tags = list(tissue="skin_v16", condition="AD", project="evo", analysis = "X2_skin_v16"),
                                     data_access = "s05")
# generate_ccm -- Thu Dec 18 15:52:10 2025: wf-a88918a66c [] - failed
# generate_ccm -- Fri Dec 19 18:41:16 2025: wf-2ec72d4434 [] - adjusted failed
# generate_ccm -- Fri Dec 19 18:42:40 2025: wf-c94ee2ad07 [] - no adjusted
# generate_ccm -- Fri Dec 19 18:43:35 2025: wf-d4626ee8f7 [] - with adjusted
# generate_ccm -- Sat Dec 20 19:12:00 2025: wf-496715f457 [] - with cells as pheno feature
# generate_ccm -- Sun Dec 28 09:20:45 2025: wf-882a48484e [] - skin_v16 (final version)


## Extract results
## ==============================
ccm = as_ccm_fit("wf-882a48484e")
skin = read.csv("~/data/skin_modified.csv", row.names = NULL, sep = "\t", header = F)
skin = rbind(skin, c("basophil","CRCL_0000357"))
signatureMapping = readRDS(get_workflow_outputs("wf-aa75ed069b"))

## Extract ct_test
## ----------------------
ct_test = rbind(statistic_table(ccm$meta$AD$ct_test),
                statistic_table(ccm$meta$L_vs_HC$ct_test),
                statistic_table(ccm$meta$L_vs_NL$ct_test),
                statistic_table(ccm$meta$NL_vs_HC$ct_test))

ct_test$feature_id = skin$V1[match(ct_test$feature_id, skin$V2)]
ct_test$skin_signature = "v16"
uploadToBQ(ct_test, bqdataset = "s05_atopic_dermatitis", tableName = "skinV16_ct_test")

ct_test$Cell = cytoreason.gx::reorder_within(ct_test$feature_id, -log10(ct_test$fdr) * sign(ct_test$estimate), ct_test$term)
ggplot(ct_test, aes(x = -log10(fdr) * sign(estimate), y = Cell)) +
  geom_col(aes(fill = ifelse(feature_id %in% c("mast cell","basophil","innate lymphoid cell 2"), "#bc5090", "grey30"))) +
  cytoreason.gx::scale_y_reordered() +
  facet_wrap(~term, scales = "free") +
  scale_fill_identity() +
  theme_minimal() +
  ggpubr::border() +
  labs(y = NULL)
ggsave("~/analysis-s05/figures/skin_v16/ct_test.png", width = 2000, units = "px", bg = "white", scale = 2)


## Extract geneset-cell correlations
## ---------------------------------------
geneset_cell_corr = statistic_table(ccm$meta$L_vs_NL$feature_cell_correlations, term = "L", filter = ~(feature_type1 == "gene_set"))
geneset_cell_corr$feature_id_2 = skin$V1[match(geneset_cell_corr$feature_id_2, skin$V2)]
geneset_cell_corr = geneset_cell_corr[,c("submodel","feature_id_1","feature_id_2","collection","value","fdr","signature_collection")]
colnames(geneset_cell_corr)[c(2:3,5)] = c("pathway","cell", "correlation")
geneset_cell_corr$log10_fdr = -log10(geneset_cell_corr$fdr)

idx = which(geneset_cell_corr$collection == "X2")
geneset_cell_corr$collection[idx] = signatureMapping$collection[match(geneset_cell_corr$pathway[idx], signatureMapping$previousID)]
geneset_cell_corr$pathway[idx] = signatureMapping$New_identifier[match(geneset_cell_corr$pathway[idx], signatureMapping$previousID)]

# Summarize random and shuffled random
for(id in c("random","top50","bottom50")){
  idx = which(str_detect(geneset_cell_corr$pathway,id))
  if(id != "random") { id = paste0("smoothedRandom_",id)}
  
  tmp = geneset_cell_corr[idx,] %>%
    dplyr::group_by(submodel, cell, collection, signature_collection) %>%
    summarise(fdr = mean(fdr),
              correlation = mean(correlation)) %>%
    mutate(log10_fdr = -log10(fdr)) %>%
    mutate(pathway = id) %>%
    dplyr::select(colnames(geneset_cell_corr))
  
  geneset_cell_corr = geneset_cell_corr[-idx,]
  geneset_cell_corr = rbind(geneset_cell_corr, tmp)
  rm(idx, tmp)
}

uploadToBQ(geneset_cell_corr, bqdataset = "s05_atopic_dermatitis", tableName = "skinV16_geneset_cell_corr")

x2 = geneset_cell_corr[which(geneset_cell_corr$collection %in% c("neuroinflammation","th2","X2","Ligands","Positives","Mast","negativeControls")),]
x2 = x2[which(x2$submodel == "bulk"),]
x2 = x2[!str_detect(x2$pathway, "CST14|PAMP12|SP|Icatibant|insilico|x2_activated|x2_general|x2_activation"),]
x2 = x2[-which(x2$pathway %in% c("X2 late activation","X2 late inhibition")),]
x2$pathway = str_remove(x2$pathway, "neuroinflammation:|th2:")

ggplot(x2, aes(x = correlation, y = log10_fdr, color = collection)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  facet_wrap(~cell, scales = "free") +
  theme_minimal() +
  ggpubr::border()
ggsave("~/analysis-s05/figures/skin_v16/geneset_cell_correlation.png", width = 5000, height = 2500, units = "px", bg = "white")



## Extract ct_test post Dupilumab
## ======================================
ct_test_post_dupi = rbind(statistic_table(analysisResultElement(ccm$datasets$GSE130588__GPL570$model$Treatment__GSE130588__GPL570__Dupilumab,"ct_test"),
                                    term = c("W16_vs_W0:DupilumabL","W16_vs_W0:DupilumabNL","W16_vs_W0:PlaceboL", "W4_vs_W0:DupilumabL")),
                          statistic_table(analysisResultElement(ccm$datasets$GSE130588__GPL570$model$Dupilumab_response__GSE130588__GPL570__R_vs_NR,"ct_test"),
                                          term = c("W4_vs_W0:NR_L", "W4_vs_W0:R_L", "W16_vs_W0:NR_L", "W4_vs_W0:R_L")))

ct_test_post_dupi$feature_id = skin$V1[match(ct_test_post_dupi$feature_id, skin$V2)]
ct_test_post_dupi = ct_test_post_dupi[,c("experiment_id", "term", "feature_id", "signature_collection", "estimate", "fdr")]
ct_test_post_dupi$log10_fdr = -log10(ct_test_post_dupi$fdr)
colnames(ct_test_post_dupi)[3] = "cell"
uploadToBQ(ct_test_post_dupi, bqdataset = "s05_atopic_dermatitis", tableName = "skinV16_ct_test_post_dupi")

ct_test_post_dupi$Cell = cytoreason.gx::reorder_within(ct_test_post_dupi$cell, -log10(ct_test_post_dupi$fdr) * sign(ct_test_post_dupi$estimate), ct_test_post_dupi$term)

ggplot(ct_test_post_dupi, aes(x = log10_fdr*sign(estimate), y = Cell)) +
  geom_col(aes(fill = ifelse(cell %in% c("mast cell","basophil","innate lymphoid cell 2"), "#bc5090", "grey30"))) +
  cytoreason.gx::scale_y_reordered() +
  scale_fill_identity() +
  facet_wrap(~term, scales = "free")
ggsave("~/analysis-s05/figures/skin_v16/ct_test_post_dupi.png", width = 5000, height = 2500, units = "px", bg = "white")


## cell loadings
## ======================================
cellLoadings = ccm$multi$meta_pca$`1`$v[,1:3]
cellLoadings[,1] = (-1) * cellLoadings[,1]
cellLoadings = reshape2::melt(cellLoadings)
colnames(cellLoadings) = c("Cell","PC","Loading")
cellLoadings$PC = paste0("PC",cellLoadings$PC)
cellLoadings$Cell = skin$V1[match(cellLoadings$Cell, skin$V2)]
uploadToBQ(cellLoadings, bqdataset = "s05_atopic_dermatitis", tableName = "skinV16_cellLoadings")


cellLoadings$cell = cytoreason.gx::reorder_within(cellLoadings$Cell, cellLoadings$Loading, cellLoadings$PC)
ggplot(cellLoadings, aes(x = Loading, y = cell)) +
  geom_col(aes(fill = ifelse(Cell %in% c("mast cell","basophil","innate lymphoid cell 2"), "#bc5090", "grey30"))) +
  scale_fill_identity() +
  cytoreason.gx::scale_y_reordered() +
  facet_wrap(~PC, scales = "free")
ggsave("~/analysis-s05/figures/skin_v16/cellLoadings.png", width = 5000, units = "px", bg = "white")
