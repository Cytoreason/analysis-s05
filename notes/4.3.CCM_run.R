devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)

ccm = as_ccm_fit("wf-08a6a0a503")
sample_ann = readRDS(get_workflow_outputs("wf-adbac4aaf3"))
pheno_vars = readRDS(get_workflow_outputs("wf-d83dfe1bc9"))

all(pheno_vars %in% colnames(sample_ann)) # needs to be TRUE


## 1. Add all phenoFeatures to the CCM object (I tried using sample_annotation_table but it didn't work)
## ===========================================================================================================
for (x in names(ccm$datasets)){
  cat("\n",x)
  pData(ccm$datasets[[x]]) <- sample_ann[which(sample_ann$geo_accession %in% pData(ccm$datasets[[x]])$geo_accession),]
}

## 2. Custom gene sets
## ============================
signatures = readRDS(get_workflow_outputs("wf-d22894f90b"))[-3]
th2 = list(th2 = readRDS(get_workflow_outputs("wf-410536ebd3")),
           neuroinflammation = readRDS(get_workflow_outputs("wf-4ce41a599a")),
           epidermis = readRDS(get_workflow_outputs("wf-d3b177cd82")))

# reorder
th2$epidermis = c(th2$epidermis, th2$th2[c("lichenification","Terminal_Differentiation_and_Lipids")])
th2$th2 = th2$th2[-which(names(th2$th2) %in% c("lichenification","Terminal_Differentiation_and_Lipids","Th1_Related","Th17_Related","TH22_IL22_Related"))]

allSignatures = c(signatures, th2)
allSignatures = lapply(allSignatures, function(x) lapply(x, unique))


## 3. Running the CCM
## ===================================
# Set gene-set size limit for ssgsea and GSEA for our custom gene lists:
gene_set_limits <- setNames(replicate(length(allSignatures), c(1L, Inf), simplify = FALSE), names(allSignatures)) # Define gene set size per signature

IMAGE <- ccm_cyto_cc_image(image = "master_1.0.1",
                           "model/feature_geneset_correlations" = "[100Gi]",
                           "meta/feature_geneset_correlations" = "[100Gi]")

ccm_api_run_custom_gene_set_analysis(ccm , # when using AssetData it will pull relevant tags automatically, and makes the relationship traceable
                                     custom_gene_set_collections = allSignatures,
                                     dataset2 = list(gene_set_activity = list(collection = c("h", "kegg", "reactome", "btm"))),
                                     model = list(gx_gsa = list(collection_size_limits = gene_set_limits),
                                                  pheno_feature_correlations = list(phenotypic_variables = pheno_vars)),
                                     meta = list(gx_diff = list(collection_size_limits = gene_set_limits)),
                                     submodel = c("bulk", "adjusted__1__1"),
                                     image = IMAGE,
                                     workflow_overrides=list(activeDeadlineSeconds=172800),
                                     tags = list(tissue="skin", condition="AD", project="evo", analysis = "X2_V9"),
                                     data_access = "s05")
# generate_ccm -- Mon Nov 24 19:05:53 2025: wf-ef57aebb52 [] - without adj pathway meta
# generate_ccm -- Tue Nov 25 09:47:02 2025: wf-8e948630d7 [] - pre th2
# generate_ccm -- Thu Dec 18 10:53:58 2025: wf-4ac16e6903 [] - sample_annotation_table, without time extension
# generate_ccm -- Thu Dec 18 11:06:02 2025: wf-3a2dc09ccc [] - workflow_overrides=list(activeDeadlineSeconds=172800) - sample_annotation_table
# generate_ccm -- Fri Dec 19 10:05:40 2025: wf-9cbfd36e26 [] - using ccm_fit instead of sample_annotation_table. Time extended.
# generate_ccm -- Sat Dec 20 19:08:52 2025: wf-8d7f6929ca [] - adding cell meta pcs to pheno features
# generate_ccm -- Sat Dec 20 19:11:05 2025: wf-3e419ff83b [] - with cells as pheno feature
# generate_ccm -- Sun Dec 28 20:58:28 2025: wf-abde4bfab0 [] - treatment meta pcs (final version)



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