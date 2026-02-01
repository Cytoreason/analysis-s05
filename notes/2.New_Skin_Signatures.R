## This script incorportes two previous scripts: Skin_v13_results.R and CCM_run_new_skin_signatures.R

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

## 1. Generating CCM runs with different skin signatures
## ======================================================================================
# 1.1. First run: integrate the skin signatures into the disease model
# ------------------------------------------------------------------
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
# generate_ccm -- Fri Dec 19 16:14:48 2025: wf-c82f4f91e9 [] - skin_v13
# generate_ccm -- Thu Dec 25 12:57:59 2025: wf-6f00628526 [] - skin_v14
# generate_ccm -- Thu Dec 25 12:57:12 2025: wf-52a7393eb5 [] - skin_v15
# generate_ccm -- Thu Dec 25 12:57:28 2025: wf-6df6f6bdca [] - skin_v16


# 1.2. Second run: analysis of the custom gene sets with the new cell signatures
# -----------------------------------------------------------------------------
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
# generate_ccm -- Sat Dec 20 19:12:00 2025: wf-496715f457 [] - skin_v13
# generate_ccm -- Sun Dec 28 09:20:45 2025: wf-882a48484e [] - skin_v16


## 2. Exploring skin_v13
## ======================================================================================
ccm_v13 = as_ccm_fit("wf-d4626ee8f7")
ccm_original = as_ccm_fit("wf-08a6a0a503")

skin = read.csv("~/data/skin_modified.csv", row.names = NULL, sep = "\t", header = F)
skin = rbind(skin, c("basophil","CRCL_0000357"))

signatureMapping = readRDS(get_workflow_outputs("wf-0cb82886cc"))

# 2.1. meta ct_test
# ----------------------
ct_test = rbind(statistic_table(ccm_v13$meta$AD$ct_test),
                statistic_table(ccm_v13$meta$L_vs_HC$ct_test),
                statistic_table(ccm_v13$meta$L_vs_NL$ct_test),
                statistic_table(ccm_v13$meta$NL_vs_HC$ct_test))
ct_test$feature_id = skin$V1[match(ct_test$feature_id, skin$V2)]
ct_test$skin_signature = "v13"
uploadToBQ(ct_test, bqdataset = "s05_atopic_dermatitis", tableName = "skinV13_ct_test")
# wf-63bb56b0aa

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


# 2.2. ct_test per dataset
# -----------------------------
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

# 2.3. geneset-cell correlations
# ----------------------------------
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
# wf-5d9977e1ed


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


# 2.4. ct_test post Dupilumab
# ----------------------------------
ct_test_post_dupi = statistic_table(analysisResultElement(ccm_v13$datasets$GSE130588__GPL570$model$Treatment__GSE130588__GPL570__Dupilumab,"ct_test"),
                                    term = c("W16_vs_W0:DupilumabL","W16_vs_W0:DupilumabNL","W16_vs_W0:PlaceboL",
                                             "W4_vs_W0:DupilumabL","W4_vs_W0:PlaceboL"))
ct_test_post_dupi$feature_id = skin$V1[match(ct_test_post_dupi$feature_id, skin$V2)]
ct_test_post_dupi = ct_test_post_dupi[,c("experiment_id", "term", "feature_id", "signature_collection", "estimate", "fdr")]
ct_test_post_dupi$log10_fdr = -log10(ct_test_post_dupi$fdr)
colnames(ct_test_post_dupi)[3] = "cell"
uploadToBQ(ct_test_post_dupi, bqdataset = "s05_atopic_dermatitis", tableName = "skinV13_ct_test_post_dupi")
# wf-0ae36f7272

ggplot(ct_test_post_dupi, aes(x = log10_fdr*sign(estimate), y = reorder(cell, -log10_fdr*sign(estimate)))) +
  geom_col() +
  facet_wrap(~term, scales = "free")
ggsave("~/analysis-s05/figures/skin_v13/ct_test_post_dupi.png", width = 5000, height = 2500, units = "px", bg = "white")


# 2.5. cell loadings
# -------------------
cellLoadings = ccm_v13$multi$meta_pca$`1`$v[,1:3]
cellLoadings[,1] = (-1) * cellLoadings[,1]
cellLoadings = reshape2::melt(cellLoadings)
colnames(cellLoadings) = c("Cell","PC","Loading")
cellLoadings$PC = paste0("PC",cellLoadings$PC)
cellLoadings$Cell = skin$V1[match(cellLoadings$Cell, skin$V2)]
uploadToBQ(cellLoadings, bqdataset = "s05_atopic_dermatitis", tableName = "skinV13_cellLoadings")
# wf-4f055620ac

ggplot(cellLoadings, aes(x = Loading, y = Cell)) +
  geom_col() +
  facet_wrap(~PC, scales = "free")
ggsave("~/analysis-s05/figures/skin_v13/cellLoadings.png", bg = "white", width = 2000, height = 900, units = "px", scale = 3)


# 2.6. cell loadings per dataset
# --------------------------------------
cellLoadings = lapply(names(ccm_v13$datasets), function(d) {
  return(
    data.frame(dataset = d, ccm_v13$datasets[[d]]$cell_pca$rotation[,1:2]))
})
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


## 3. Comparing skin signatures (Gil Eshel's analysis)
## ======================================================================================
# Function to extract ct-test results from CCM wf-id:
get_meta_ct_test <- function(ccm_wf_id="wf-08a6a0a503",model=c("AD","L_vs_HC","L_vs_NL","NL_vs_HC")){
  #' Decide whether feature_id_1 contains Cytoreason cell ids or cell-types
  #' @param x Character vector of cell identifiers from feature_id_1
  #' @return Logical scalar; TRUE means translate_ids2names should be attempted
  should_translate_cells <- function(x){
    #
    x <- unique(as.character(x))
    x <- x[!is.na(x) & x!=""]
    if(length(x)==0)return(FALSE)
    has_spaces <- mean(grepl("\\s",x))>0.1
    has_letters <- mean(grepl("[A-Za-z]",x))>0.2
    has_punct <- mean(grepl("[/_:]",x))>0.2
    all_numeric <- all(grepl("^[0-9]+$",x))
    likely_name <- has_spaces || (has_letters && !has_punct)
    likely_id <- all_numeric || (has_punct && !has_spaces)
    if(likely_name)return(FALSE)
    if(likely_id)return(TRUE)
    return(FALSE)
  }
  #
  ccm_fit <- as_ccm_fit(ccm_wf_id)
  cell_ids <- rownames(ccm_fit[["multi"]][["meta_pca"]][["1"]][["v"]])
  if(should_translate_cells(cell_ids)){
    cell_id_name <- cytoreason.datasource::translate_ids2names(cell_ids)
  }else{
    cell_id_name <- cytoreason.datasource::translate_names2ids(cell_ids)
  }
  #
  ct_test_res <- do.call('rbind',lapply(model,function(m){
    ct_test <- statistic_table(ccm_fit[["meta"]][[m]][["ct_test"]])
    ct_test[["model"]] <- m
    if(should_translate_cells(cell_ids)){
      ct_test[["cell_id"]] <- ct_test[["feature_id"]]
      ct_test[["cell_name"]] <- cell_id_name[["cell_type_name"]][match(ct_test[["feature_id"]],cell_id_name[["cell_type_id"]])]
    }else{
      ct_test[["cell_id"]] <- cell_id_name[["cell_type_id"]][match(ct_test[["feature_id"]],cell_id_name[["queried_cell_type_name"]])]
      ct_test[["cell_name"]] <- ct_test[["feature_id"]]
    }
    ct_test
  }))
  ct_test_res
}

ct_test_original <- get_meta_ct_test(ccm_wf_id="wf-08a6a0a503",model=c("AD","L_vs_HC","L_vs_NL","NL_vs_HC")); ct_test_original[["run"]] <- "original"
ct_test_skin_v13 <- get_meta_ct_test(ccm_wf_id="wf-d4626ee8f7",model=c("AD","L_vs_HC","L_vs_NL","NL_vs_HC")); ct_test_skin_v13[["run"]] <- "skin_v13"
ct_test_skin_v14 <- get_meta_ct_test(ccm_wf_id="wf-6f00628526",model=c("AD","L_vs_HC","L_vs_NL","NL_vs_HC")); ct_test_skin_v14[["run"]] <- "skin_v14"
ct_test_skin_v15 <- get_meta_ct_test(ccm_wf_id="wf-52a7393eb5",model=c("AD","L_vs_HC","L_vs_NL","NL_vs_HC")); ct_test_skin_v15[["run"]] <- "skin_v15"
ct_test_skin_v16 <- get_meta_ct_test(ccm_wf_id="wf-6df6f6bdca",model=c("AD","L_vs_HC","L_vs_NL","NL_vs_HC")); ct_test_skin_v16[["run"]] <- "skin_v16"

ct_test_res <- rbind(ct_test_original,ct_test_skin_v13,ct_test_skin_v14,ct_test_skin_v15,ct_test_skin_v16)

(wf_ct_test_res <- run_function_dist(FUN=function(my_object){return (my_object)},
                                     my_object=ct_test_res,
                                     image="eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_0.55.1")) # wf-9b5a81b0b3
ct_test_res <- readRDS(get_task_outputs('wf-9b5a81b0b3',"0")) # ct_test_res

# Save ct-test results on BQ:
bq_perform_upload(bq_table("cytoreason","s05_atopic_dermatitis",table="skinSignatures_ct_test"),ct_test_res, disposition = "WRITE_TRUNCATE")



## 4. Exploring skin_v16
## ======================================================================================
ccm_v16 = as_ccm_fit("wf-882a48484e")
skin = read.csv("~/data/skin_modified.csv", row.names = NULL, sep = "\t", header = F)
skin = rbind(skin, c("basophil","CRCL_0000357"))
signatureMapping = readRDS(get_workflow_outputs("wf-aa75ed069b"))

# 4.1. meta ct_test
# ----------------------
ct_test = rbind(statistic_table(ccm_v16$meta$AD$ct_test),
                statistic_table(ccm_v16$meta$L_vs_HC$ct_test),
                statistic_table(ccm_v16$meta$L_vs_NL$ct_test),
                statistic_table(ccm_v16$meta$NL_vs_HC$ct_test))

ct_test$feature_id = skin$V1[match(ct_test$feature_id, skin$V2)]
ct_test$skin_signature = "v16"
uploadToBQ(ct_test, bqdataset = "s05_atopic_dermatitis", tableName = "skinV16_ct_test")
# wf-28a0a166df

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


# 4.2. geneset-cell correlations
# -------------------------------------
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
# wf-3d3dc01301

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



# 4.3. Extract ct_test post Dupilumab
# -----------------------------------------
ct_test_post_dupi = rbind(statistic_table(analysisResultElement(ccm$datasets$GSE130588__GPL570$model$Treatment__GSE130588__GPL570__Dupilumab,"ct_test"),
                                          term = c("W16_vs_W0:DupilumabL","W16_vs_W0:DupilumabNL","W16_vs_W0:PlaceboL", "W4_vs_W0:DupilumabL")),
                          statistic_table(analysisResultElement(ccm$datasets$GSE130588__GPL570$model$Dupilumab_response__GSE130588__GPL570__R_vs_NR,"ct_test"),
                                          term = c("W4_vs_W0:NR_L", "W4_vs_W0:R_L", "W16_vs_W0:NR_L", "W4_vs_W0:R_L")))

ct_test_post_dupi$feature_id = skin$V1[match(ct_test_post_dupi$feature_id, skin$V2)]
ct_test_post_dupi = ct_test_post_dupi[,c("experiment_id", "term", "feature_id", "signature_collection", "estimate", "fdr")]
ct_test_post_dupi$log10_fdr = -log10(ct_test_post_dupi$fdr)
colnames(ct_test_post_dupi)[3] = "cell"
uploadToBQ(ct_test_post_dupi, bqdataset = "s05_atopic_dermatitis", tableName = "skinV16_ct_test_post_dupi")
# wf-a7d4dd6d7f

ct_test_post_dupi$Cell = cytoreason.gx::reorder_within(ct_test_post_dupi$cell, -log10(ct_test_post_dupi$fdr) * sign(ct_test_post_dupi$estimate), ct_test_post_dupi$term)

ggplot(ct_test_post_dupi, aes(x = log10_fdr*sign(estimate), y = Cell)) +
  geom_col(aes(fill = ifelse(cell %in% c("mast cell","basophil","innate lymphoid cell 2"), "#bc5090", "grey30"))) +
  cytoreason.gx::scale_y_reordered() +
  scale_fill_identity() +
  facet_wrap(~term, scales = "free")
ggsave("~/analysis-s05/figures/skin_v16/ct_test_post_dupi.png", width = 5000, height = 2500, units = "px", bg = "white")


# 4.4 cell loadings
# --------------------------
cellLoadings = ccm$multi$meta_pca$`1`$v[,1:3]
cellLoadings[,1] = (-1) * cellLoadings[,1]
cellLoadings = reshape2::melt(cellLoadings)
colnames(cellLoadings) = c("Cell","PC","Loading")
cellLoadings$PC = paste0("PC",cellLoadings$PC)
cellLoadings$Cell = skin$V1[match(cellLoadings$Cell, skin$V2)]
uploadToBQ(cellLoadings, bqdataset = "s05_atopic_dermatitis", tableName = "skinV16_cellLoadings")
# wf-fb527cb845

cellLoadings$cell = cytoreason.gx::reorder_within(cellLoadings$Cell, cellLoadings$Loading, cellLoadings$PC)
ggplot(cellLoadings, aes(x = Loading, y = cell)) +
  geom_col(aes(fill = ifelse(Cell %in% c("mast cell","basophil","innate lymphoid cell 2"), "#bc5090", "grey30"))) +
  scale_fill_identity() +
  cytoreason.gx::scale_y_reordered() +
  facet_wrap(~PC, scales = "free")
ggsave("~/analysis-s05/figures/skin_v16/cellLoadings.png", width = 5000, units = "px", bg = "white")
