devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)

## 1. Inputs
## ================================
ccm = as_ccm_fit("wf-08a6a0a503")
## removing the extra data set isn't feasible within the constraints of the project - need to re-run the disease model, causing compatibility issues

# preparing the sample annotation table we will pass, and integrate MS/CS/pheno into it
sample_ann = lapply(ccm$datasets, function(d){
  p = pData(assayDataExpression(d))
  if("week:ch1" %in% colnames(p)) {
    p[,"week:ch1"] = as.character(p[,"week:ch1"])
  }
  if("Week" %in% colnames(p)) {
    p[,"Week"] = as.character(p[,"Week"])
  }
  if("batch_date:ch1" %in% colnames(p)) {
    p[,"batch_date:ch1"] = as.character(p[,"batch_date:ch1"])
  }
  if("Batch_date" %in% colnames(p)) {
    p[,"Batch_date"] = as.character(p[,"Batch_date"])
  }
  if("channel_count" %in% colnames(p)) {
    p[,"channel_count"] = as.character(p[,"channel_count"])
  }
  if("taxid_ch1" %in% colnames(p)) {
    p[,"taxid_ch1"] = as.character(p[,"taxid_ch1"])
  }
  if("data_row_count" %in% colnames(p)) {
    p[,"data_row_count"] = as.character(p[,"data_row_count"])
  }
  if("patient_id:ch1" %in% colnames(p)) {
    p[,"patient_id:ch1"] = as.character(p[,"patient_id:ch1"])
  }
  if("Patient_id" %in% colnames(p)) {
    p[,"Patient_id"] = as.character(p[,"Patient_id"])
  }
  if("description.1" %in% colnames(p)) {
    p[,"description.1"] = as.character(p[,"description.1"])
  }
  return(p)
})
sample_ann = bind_rows(sample_ann)


## 2. Integrate molecular score
## ================================
MS = readRDS(get_workflow_outputs("wf-db0c88cca6"))
# calculated by the data science group https://cytoreason.atlassian.net/wiki/spaces/DSG/pages/4154261510/New+Molecular+Score+-+Documentation#Molecular-score-results-on-different-disease-models

MS = pData(MS$ms_results$combat_ms) # MS based on L vs HC
MS = MS %>%
  dplyr::select(condition, sample_classification, geo_accession,molecular_score) %>%
  mutate(condition = str_replace(condition, "control", "healthy")) %>%
  mutate(condition = factor(condition, ordered = T, levels = c("atopic dermatitis","contact dermatitis","psoriasis","healthy")))

sample_ann$MolecularScore = MS$molecular_score[match(sample_ann$geo_accession, MS$geo_accession)]

pheno_vars = "MolecularScore"
               
## 3. Unify SCORAD/EASI scores
## ================================
sample_ann <- sample_ann %>% 
  dplyr::select(-c("SCORAD","Scorad","easi:ch1")) %>% # redundant columns
  dplyr::rename(SCORAD = "scorad:ch1")

pheno_vars = c(pheno_vars, "SCORAD","EASI")

## 4. Integrate with pathways meta pcs
## ======================================
# 4.1. All meta PCs
# --------------------
pathwayPCs = readRDS(get_workflow_outputs("wf-39ae62c8d2"))
pathwayPCs = pathwayPCs %>%
  dplyr::filter(collection == "all" & term %in% c("AD","L_vs_HC")) %>%
  mutate(term = ifelse(term == "AD", "DZ_vs_HC", term)) %>%
  mutate(submodel = ifelse(submodel != "bulk", "adjusted", "bulk")) %>%
  mutate(suffix = paste0(term,"_",submodel,"_",collection)) %>%
  dplyr::select(-c(term:collection))

pathwayPCs = reshape2::melt(pathwayPCs, id.vars = c(1,2,6), measured.vals = 3:5, variable.name = "prefix", value.name = "sampleScore")
pathwayPCs = pathwayPCs %>%
  mutate(label = paste0(prefix,"_",suffix)) %>%
  dplyr::select(-c(suffix:prefix))

pathwayPCs = reshape2::dcast(pathwayPCs, dataset + sample_id ~ label, value.var = "sampleScore")
colnames(pathwayPCs) = str_replace(colnames(pathwayPCs),"_all","_key")

sample_ann = merge(sample_ann, pathwayPCs, by.x = c("sample_id", "dataset_id"), by.y = c("sample_id", "dataset"), all = T)
pheno_vars = c(pheno_vars, colnames(pathwayPCs[,-c(1:2)]))
               
# 4.2. Th2 meta PCs
# ----------------------
pathwayPCs_th2_bulk = readRDS(get_workflow_outputs("wf-a7ac64d7d6"))
pathwayPCs_th2_bulk = lapply(names(pathwayPCs_th2_bulk), function(x) {
  y = pathwayPCs_th2_bulk[[x]]
  colnames(y) = c("sample_id",paste0("pathway_meta_pc",1:3))
  colnames(y)[-1] = paste0(colnames(y)[-1],"_",x,"_bulk")
  return(y)
}) %>% purrr::reduce(., full_join, by = "sample_id")

pathwayPCs_th2_adj = readRDS(get_workflow_outputs("wf-f76769a6b6"))
pathwayPCs_th2_adj = lapply(names(pathwayPCs_th2_adj), function(x) {
  y = pathwayPCs_th2_adj[[x]]
  colnames(y) = c("sample_id",paste0("pathway_meta_pc",1:3))
  colnames(y)[-1] = paste0(colnames(y)[-1],"_",x,"_adj")
  return(y)
}) %>% purrr::reduce(., full_join, by = "sample_id")

sample_ann = merge(sample_ann, pathwayPCs_th2_bulk, by = "sample_id", all = T)
sample_ann = merge(sample_ann, pathwayPCs_th2_adj, by = "sample_id", all = T)

pheno_vars = c(pheno_vars, colnames(pathwayPCs_th2_bulk[,-1]), colnames(pathwayPCs_th2_adj[,-1]))


# 4.3. Key + Th2 meta PCs
# --------------------------
pathwayPCs = readRDS(get_workflow_outputs("wf-6dba8d483f"))
pathwayPCs = pathwayPCs %>%
  mutate(term = ifelse(term == "AD", "DZ_vs_HC", term)) %>%
  mutate(submodel = ifelse(submodel != "bulk", "adjusted", "bulk")) %>%
  mutate(suffix = paste0(term,"_",submodel,"_",collection)) %>%
  dplyr::select(-c(term:collection))

pathwayPCs = reshape2::melt(pathwayPCs, id.vars = c(1,2,6), measured.vals = 3:5, variable.name = "prefix", value.name = "sampleScore")
pathwayPCs = pathwayPCs %>%
  mutate(label = paste0(prefix,"_",suffix)) %>%
  dplyr::select(-c(suffix:prefix))

pathwayPCs = reshape2::dcast(pathwayPCs, dataset + sample_id ~ label, value.var = "sampleScore")

sample_ann = merge(sample_ann, pathwayPCs, by = "sample_id", all = T)
pheno_vars = c(pheno_vars, colnames(pathwayPCs[,-c(1:2)]))


## 5. Custom gene sets
## ============================
signatures = readRDS(get_workflow_outputs("wf-d22894f90b"))[-3]
th2 = list(th2 = readRDS(get_workflow_outputs("wf-410536ebd3")),
           neuroinflammation = readRDS(get_workflow_outputs("wf-4ce41a599a")),
           epidermis = readRDS(get_workflow_outputs("wf-d3b177cd82")))

allSignatures = c(signatures, th2)
allSignatures = lapply(allSignatures, function(x) lapply(x, unique))

# 
# ## 6. Add canonical pathways to phenoFeature
# ## ===================================================
# canonical = cytoreason.gx::load_gene_set_collections(collection = c("kegg","btm","reactome","h"))
# names(canonical) = c("KEGG","BTM","REACTOME","HALLMARK")

## 6. Add all phenoFeatures to the CCM object (I tried using sample_annotation_table but it didn't work)
## ===========================================================================================================
sample_ann$meta1_pc1 = (-1) * sample_ann$meta1_pc1 # flipped the sign of the cell meta-PC1 (because HC>DZ)
for (x in names(ccm$datasets)){
  cat("\n",x)
  pData(ccm$datasets[[x]]) <- sample_ann[which(sample_ann$geo_accession %in% pData(ccm$datasets[[x]])$geo_accession),]
}


pushToCC(sample_ann, tagsToPass = list(list(name="object", value="AD_sample_annotation")))
# wf-558ccafbfe

pheno_vars = c(pheno_vars, "meta1_pc1","meta1_pc2","meta1_pc3")
pushToCC(pheno_vars, tagsToPass = list(list(name="object", value="AD_pheno_vars")))
# wf-4c71fb2030

## 7. Running the CCM
## ===================================
# Set gene-set size limit for ssgsea and GSEA for our custom gene lists:
gene_set_limits <- setNames(replicate(length(allSignatures), c(1L, Inf), simplify = FALSE), names(allSignatures)) # Define gene set size per signature

IMAGE <- ccm_cyto_cc_image(image = "master_1.0.1",
                           "model/feature_geneset_correlations" = "[100Gi]")

ccm_api_run_custom_gene_set_analysis(ccm , # when using AssetData it will pull relevant tags automatically, and makes the relationship traceable
                                     custom_gene_set_collections = allSignatures,
                                     dataset2 = list(gene_set_activity = list(collection = c("h", "kegg", "reactome", "btm"))),
                                     model = list(gx_gsa = list(collection_size_limits = gene_set_limits),
                                                  pheno_feature_correlations = list(phenotypic_variables = pheno_vars)),
                                     meta = list(gx_diff = list(collection_size_limits = gene_set_limits)),
                                     submodel = c("bulk", "adjusted__1__1"),
                                     image = IMAGE,
                                     workflow_overrides=list(activeDeadlineSeconds=172800),
                                     tags = list(tissue="skin", condition="AD", project="evo", analysis = "X2_V8"),
                                     data_access = "s05")
# generate_ccm -- Mon Nov 24 19:05:53 2025: wf-ef57aebb52 [] - without adj pathway meta
# generate_ccm -- Tue Nov 25 09:47:02 2025: wf-8e948630d7 []
# generate_ccm -- Wed Dec 17 05:24:53 2025: wf-30b44d3c6f [] - some tasks failed. integration of many meta PCs incl. Th2
# generate_ccm -- Thu Dec 18 10:53:58 2025: wf-4ac16e6903 [] - sample_annotation_table, didn't work
# generate_ccm -- Thu Dec 18 11:06:02 2025: wf-3a2dc09ccc [] - workflow_overrides=list(activeDeadlineSeconds=172800) - sample_annotation_table, didn't work
# generate_ccm -- Fri Dec 19 10:05:40 2025: wf-9cbfd36e26 [] - using ccm_fit instead of sample_annotation_table. Time extended.
# generate_ccm -- Sat Dec 20 19:08:52 2025: wf-8d7f6929ca [] - adding cell meta pcs to pheno features
# generate_ccm -- Sat Dec 20 19:11:05 2025: wf-3e419ff83b [] - with cells as pheno feature

## 7. Running the CCM with skin_v13
## ===================================
# In this part, we use a hybrid skin signature that includes mast cells, basophils and ILC2 from lung_v4
# in order to mitigate some concerns by the client that we are not capturing type 2 immunity
config = read.csv(get_task_inputs("wf-08a6a0a503","0", files_names_grepl_pattern = "ccm-metadata.csv"))
config$asset_id[which(config$experiment_id == "GSE153007")] <- 'wf-2ac040bba1:0:GSE153007.RDS'
config$asset_id[which(config$experiment_id == "GSE60709")] <- 'wf-923925db68:0:GSE60709.RDS'
config$ccm_annotate = "feature"

sample_ann = readRDS(get_workflow_outputs("wf-558ccafbfe"))

ccm_stage_prepare_dataset_collection(config, annotate = list(sample_annotation_table = sample_ann))

IMAGE <- ccm_cyto_cc("master_1.0.1",
                     dataset.cell_contribution.service = "eu.gcr.io/cytoreason/ci-cytoreason.deconvolution-package:develop_2.2.7",
                     cell_specific_differences = "eu.gcr.io/cytoreason/ci-cytoreason.deconvolution-package:develop_2.2.7")

ccm_api_generate_ccm(config,
                     qc = FALSE,
                     term_metadata = FALSE, # fails for some reason
                     prepare_data = list(annotate = list(sample_annotation_table = sample_ann)),
                     dataset = list(cell_contribution = list(signature_collection = "skin_v13:SUP-6496-skin_v13"), 
                                    services = .skip("feature_pca")),
                     dataset2 = list(gene_set_activity = list(collection = c("h", "kegg", "reactome", "btm"))),
                     model = list(services = .skip("gene_set_activity_differences","feature_geneset_correlations")),
                     workflow_overrides=list(activeDeadlineSeconds=172800),
                     image = IMAGE,
                     tags = list(tissue="skin_v13:SUP-6496-skin_v13", condition="AD", project="evo", analysis="AD_skin_v13_reducedAnalyses"))
# generate_ccm -- Fri Dec 19 10:48:35 2025: wf-e82c5ab874 []
# generate_ccm -- Fri Dec 19 15:06:38 2025: wf-5f05a83bc6 [] - tag with skin_v13
# generate_ccm -- Fri Dec 19 15:06:51 2025: wf-dd8c8843e6 [] - tag with skin_v13:SUP-6496-skin_v13
# generate_ccm -- Fri Dec 19 16:14:48 2025: wf-c82f4f91e9 [] - change deconvolution image

ccm_api_run_custom_gene_set_analysis(cytoreason.assets::AssetData("wf-c82f4f91e9") , # when using AssetData it will pull relevant tags automatically, and makes the relationship traceable
                                     custom_gene_set_collections = allSignatures,
                                     prepare_data = list(annotate = list(sample_annotation_table = sample_ann)),
                                     dataset = list(cell_contribution = list(signature_collection = "skin_v13:SUP-6496-skin_v13"),
                                                    services = .skip("feature_pca")),
                                     dataset2 = list(gene_set_activity = list(collection = c("h", "kegg", "reactome", "btm"))),
                                     model = list(pheno_feature_correlations = list(phenotypic_variables = pheno_vars)),
                                     meta = list(gx_diff = list(collection_size_limits = gene_set_limits)),
                                     submodel = c("bulk","adjusted__1__meta1_pc1"), # adjusted failed
                                     workflow_overrides=list(activeDeadlineSeconds=172800),
                                     image = IMAGE,
                                     tags = list(tissue="skin_v13", condition="AD", project="evo", analysis = "X2_V7_skin_v13"),
                                     data_access = "s05")
# generate_ccm -- Thu Dec 18 15:52:10 2025: wf-a88918a66c [] - failed
# generate_ccm -- Fri Dec 19 18:41:16 2025: wf-2ec72d4434 [] - adjusted failed
# generate_ccm -- Fri Dec 19 18:42:40 2025: wf-c94ee2ad07 [] - no adjusted
# generate_ccm -- Fri Dec 19 18:43:35 2025: wf-d4626ee8f7 [] - with adjusted
# generate_ccm -- Sat Dec 20 19:12:00 2025: wf-496715f457 [] - with cells as pheno feature

## 7. Only correlations to SCORAD
## ===================================
# Pull the config
config = read.csv(get_task_inputs("wf-08a6a0a503","0", files_names_grepl_pattern = "ccm-metadata.csv"))

# Pull the sample annotation
sample_ann = readRDS(get_workflow_outputs("wf-558ccafbfe"))
scorad = unique(sample_ann$experiment_id[which(!is.na(sample_ann$SCORAD))])

sample_ann_sub = sample_ann[which(sample_ann$experiment_id %in% scorad),]
sample_ann_sub$SCORAD2 = sample_ann_sub$SCORAD # to be certain we are calculating on the right samples

# for GSE27887 inputting SCORAD from the paper
GSE27887 = read.csv("~/analysis-s05/data/GSE27887_SCORAD.csv", sep = "\t")
idx = which(sample_ann_sub$experiment_id == "GSE27887")
sample_ann_sub$SCORAD2[idx] = NA # I don't trust the annotations there

idx = which(sample_ann_sub$experiment_id == "GSE27887" & sample_ann_sub$sample_classification == "Lesion" & sample_ann_sub$time == "Pre")
sample_ann_sub$SCORAD2[idx] = GSE27887$Pre[match(sample_ann_sub$Subject_ID[idx], GSE27887$Patient)]

idx = which(sample_ann_sub$experiment_id == "GSE27887" & sample_ann_sub$sample_classification == "Lesion" & sample_ann_sub$time == "Post")
sample_ann_sub$SCORAD2[idx] = GSE27887$Post[match(sample_ann_sub$Subject_ID[idx], GSE27887$Patient)]

# Change the config to calculate the wanted correlation
config = config[which(config$experiment_id %in% scorad & config$effect_id == "L_vs_NL"),]
config$model_pairing[1] = "Subject.Id"
config$ccm_exclude = NA
config$dataset_id = paste0(config$experiment_id, "__", config$platform_id)
config$comparison_id = paste0("L_vs_NL__",config$dataset_id,"__Lesion_vs_non_lesion")
config$ccm_meta_pca = NA

# signatures
signatures = readRDS(get_workflow_outputs("wf-d22894f90b"))[-3]
th2 = list(th2 = readRDS(get_workflow_outputs("wf-410536ebd3")),
           neuroinflammation = readRDS(get_workflow_outputs("wf-4ce41a599a")),
           epidermis = readRDS(get_workflow_outputs("wf-d3b177cd82")))
th2$th2$McCluskey = unique(th2$th2$McCluskey)
allSignatures = c(signatures, th2)

ccm_stage_prepare_dataset_collection(config, annotate = list(sample_annotation_table = sample_ann_sub))

gene_set_limits <- setNames(replicate(length(allSignatures), c(1L, Inf), simplify = FALSE), names(allSignatures)) # Define gene set size per signature

ccm_api_run_custom_gene_set_analysis(config, 
                                     custom_gene_set_collections = allSignatures,
                                     prepare_data = list(annotate = list(sample_annotation_table = sample_ann_sub)),
                                     multi = list(services = .skip("meta_pca")),
                                     dataset2 = list(services=c("adjust_data","gene_set_activity","annotate")),
                                     model = list(services = c("pheno_feature_correlations"),
                                                  pheno_feature_correlations = list(phenotypic_variables = "SCORAD2")),
                                     meta = list(services = "pheno_feature_correlations"),
                                     image = "master_1.0.1",
                                     tags = list(tissue="skin", condition="AD", project="evo", analysis = "SCORAD_meta_correlation"),
                                     data_access = "s05")
# generate_ccm -- Wed Dec 17 15:44:45 2025: wf-0e2c39f807 [] - succeeded but no meta :(
# generate_ccm -- Thu Dec 18 14:05:23 2025: wf-47f31aef48 [] - succeeded but forgot to remove the data from the weird dataset
# generate_ccm -- Fri Dec 19 11:20:40 2025: wf-0c357e0a58 [] - using SCOARD2 (without the dataset) - worked but without the weird dataset I didn't get meta
# generate_ccm -- Fri Dec 19 14:58:05 2025: wf-b1d96b1738 [] - added SCORAD from the paper


################### Testing with Renaud and Matan
config = read.csv(get_task_inputs("wf-08a6a0a503","0", files_names_grepl_pattern = "ccm-metadata.csv"))
config$asset_id[which(config$experiment_id == "GSE153007")] <- 'wf-2ac040bba1:0:GSE153007.RDS'
config$asset_id[which(config$experiment_id == "GSE60709")] <- 'wf-923925db68:0:GSE60709.RDS'
sample_ann = lapply(ccm$datasets, function(d){
  p = pData(assayDataExpression(d))
  if("week:ch1" %in% colnames(p)) {
    p[,"week:ch1"] = as.character(p[,"week:ch1"])
  }
  if("Week" %in% colnames(p)) {
    p[,"Week"] = as.character(p[,"Week"])
  }
  if("batch_date:ch1" %in% colnames(p)) {
    p[,"batch_date:ch1"] = as.character(p[,"batch_date:ch1"])
  }
  if("Batch_date" %in% colnames(p)) {
    p[,"Batch_date"] = as.character(p[,"Batch_date"])
  }
  if("channel_count" %in% colnames(p)) {
    p[,"channel_count"] = as.character(p[,"channel_count"])
  }
  if("taxid_ch1" %in% colnames(p)) {
    p[,"taxid_ch1"] = as.character(p[,"taxid_ch1"])
  }
  if("data_row_count" %in% colnames(p)) {
    p[,"data_row_count"] = as.character(p[,"data_row_count"])
  }
  if("patient_id:ch1" %in% colnames(p)) {
    p[,"patient_id:ch1"] = as.character(p[,"patient_id:ch1"])
  }
  if("Patient_id" %in% colnames(p)) {
    p[,"Patient_id"] = as.character(p[,"Patient_id"])
  }
  if("description.1" %in% colnames(p)) {
    p[,"description.1"] = as.character(p[,"description.1"])
  }
  return(p)
})
sample_ann = bind_rows(sample_ann)
config$ccm_annotate = "feature"

ccm_stage_prepare_dataset_collection(config, annotate = list(sample_annotation_table = sample_ann))

ccm_api_generate_ccm(config,
                     prepare_data = list(annotate = list(sample_annotation_table = sample_ann)),
                     qc = FALSE,
                     term_metadata = FALSE, # fails for some reason
                     adjustment_models = list("1" = list(c("meta1_pc1"), c("CRCL_0000348"),c("meta1_pc1","CRCL_0000348"))),
                     model = .skip("cell_specific_differences"),
                     image = "master_1.0.1",
                     tags = list(tissue="skin", condition="AD", project="evo", analysis="test"))
# generate_ccm -- Tue Dec 16 12:24:40 2025: wf-d9aec2329a []