library(cytoreason.cc.client)
devtools::load_all("~/cytoreason.ccm.pipeline/")
library(tidyverse)

# 1. Inputs
# -----------------------------
ccm = as_ccm_fit("wf-08a6a0a503")
<<<<<<< Updated upstream
signatures = readRDS(get_workflow_outputs("wf-b33649db09"))
## make sure:
## A. removing the extra dataset
## all ligands is in entrez
=======
signatures = readRDS(get_workflow_outputs("wf-54d564beb9"))
## removing the extra dataset isn't feasable within the contstraints of the project - need to re-run the disease model, causing compatibility issues
>>>>>>> Stashed changes


# 2. Molecular and clinical scores
# ------------------------------------
MS = cytoreason.assets::read_asset("wf-666afb07dc") # generated in a p13 project: https://github.com/Cytoreason/analysis-p13/blob/master/notes/MS/AD/MS_for_AD_v2.rmd
MS$gene_set_res$calculated_score$clinical_score = MS$res$calculated_score$clinical_score[match(MS$gene_set_res$calculated_score$sample_id, MS$res$calculated_score$sample_id)]
MS$gene_set_res$calculated_score$clinical_score[which(MS$gene_set_res$calculated_score$clinical_score == -1)] <- NA

MS.old = MS$gene_set_res$calculated_score # Custom signature MS, based on DZ vs HC
MS.old = MS.old %>%
  dplyr::select(sample_id,dataset,norm_molecular_score) %>%
  dplyr::mutate(sample_id = as.character(sample_id)) %>%
  dplyr::mutate(dataset = as.character(dataset)) %>%
  dplyr::mutate(MS.old = 1-norm_molecular_score) # score is inverted

MS.new = MS$res$calculated_score # MS based on L vs NL
MS.new = MS.new %>%
  dplyr::select(sample_id,dataset,norm_molecular_score) %>%
  dplyr::mutate(sample_id = as.character(sample_id)) %>%
  dplyr::mutate(dataset = as.character(dataset)) %>%
  dplyr::mutate(MS.new = 1-norm_molecular_score) # score is inverted

MS = merge(MS.old[,c(1,2,4)], MS.new[,c(1,2,4)])
rm(MS.new, MS.old)

# add molecular score & clinical information per sample in the pData
for (x in names(ccm$datasets)){
  sample_ids <- sampleNames(ccm$datasets[[x]])
  MS.old <- MS$MS.old[match(sample_ids, MS$sample_id)]
  MS.new <- MS$MS.new[match(sample_ids, MS$sample_id)]
  pData(ccm$datasets[[x]])$MS.old <- MS.old
  pData(ccm$datasets[[x]])$MS.new <- MS.new
  
  if(x %in% c("GSE130588__GPL570","GSE27887__GPL570","GSE58558__GPL570")) { # there we have SCORAD
    idx = which(str_detect(colnames(pData(ccm$datasets[[x]])),"SCORAD|Scorad"))[1]
    pData(ccm$datasets[[x]])$SCORAD = as.numeric(pData(ccm$datasets[[x]])[,idx])
    rm(idx)
  }
  if(x %in% c("GSE130588__GPL570","GSE59294__GPL570")) { # there we have EASI
    idx = which(str_detect(colnames(pData(ccm$datasets[[x]])),"EASI"))[1]
    pData(ccm$datasets[[x]])$EASI = as.numeric(pData(ccm$datasets[[x]])[,idx])
    rm(idx)
  }
  
  rm(sample_ids, MS.new, MS.old)
}


## General
## ------------------
# Set gene-set size limit for ssgsea and GSEA for our custom gene lists:
gene_set_limits <- setNames(replicate(length(signatures), c(1L, Inf), simplify = FALSE), names(signatures)) # Define gene set size per signature


## Running the CCM
## ----------------------
ccm_api_run_custom_gene_set_analysis(ccm, # when using AssetData it will pull relevant tags automatically, and makes the relationship traceable
                                     custom_gene_set_collections = signatures,
                                     model = list(gx_gsa = list(collection_size_limits = gene_set_limits),
                                                  pheno_feature_correlations = list(phenotypic_variables = c("MS.old","MS.new","SCORAD","EASI"))), # add the new phenotypic variables
                                     meta = list(gx_diff = list(collection_size_limits = gene_set_limits)),
                                     submodel = c("bulk", "adjusted__1__1"),
                                     image = "master_1.0.1",
                                     tags = list(list(tissue="skin"),
                                                 list(condition="AD"),
                                                 list(project="evo"),
<<<<<<< Updated upstream
                                                 list(name="analysis",value="X2_V1")),
                                     data_access = "s05")
# generate_ccm -- Fri Nov 14 15:07:26 2025: wf-eafe4014ef []
=======
                                                 list(name="analysis",value="X2_V3")),
                                     data_access = "s05")
# generate_ccm -- Mon Nov 17 13:12:37 2025: wf-c022438e0b []
# generate_ccm -- Mon Nov 17 15:31:47 2025: wf-ab6dd9441b [] - fixed ep
# generate_ccm -- Tue Nov 18 09:12:32 2025: wf-67e687bf9f []
# generate_ccm -- Thu Nov 20 07:54:00 2025: wf-5d25fa144d []

########### This did not work due to using a legacy model
## 2. Modify CCM
## ===========================
config = read.csv(get_task_inputs("wf-08a6a0a503", task_id = 0, files_names_grepl_pattern = "ccm-metadata"))
config$ccm_exclude[which(config$experiment_id == "GSE147424")] <- 1
config$asset_id[which(config$asset_id == "18057:adhoc:GSE153007.RDS")] <- "wf-381431a592:0:GSE153007.RDS"
config$asset_id[which(config$asset_id == "18057:adhoc:GSE60709.RDS")] <- "wf-ff39d033df:0:GSE60709.RDS"

ccm_stage_prepare_dataset_collection(config)
ccm_api_generate_ccm(config = config,
                     submodel = c("bulk", "adjusted__1__1"),
                     image = "master_1.0.1",
                     tags = list(list(tissue="skin"),
                                 list(condition="AD"),
                                 list(name="analysis",value="DiseaseModel_sansGSE147424")))
# generate_ccm -- Mon Nov 17 08:03:09 2025: wf-ab6dd9441b []


########### This did not work because it ran a network stage, and upon skipping threw a bunch of errors
## Running the CCM
## ----------------------
UnifiedMetadata = readRDS(get_workflow_outputs("wf-e82a5ab3b6")) # from p13

sample_annotation = merge(MS, UnifiedMetadata[,c("sample_id", "clinical_score_value_SCORAD", "clinical_score_value_EASI")], by = "sample_id", all = T)
sample_annotation = sample_annotation %>%
  rename(SCORAD = clinical_score_value_SCORAD, EASI = clinical_score_value_EASI)


ccm_api_run_phenofeature_correlations(ccm, # when using AssetData it will pull relevant tags automatically, and makes the relationship traceable
                                      custom_gene_set_collections = signatures,
                                      stages = .skip("network"),
                                      prepare_data = list(annotate = list(sample_annotation_table = sample_annotation)),
                                      model = list(gx_gsa = list(collection_size_limits = gene_set_limits)),
                                      meta = list(gx_diff = list(collection_size_limits = gene_set_limits)),
                                      submodel = c("bulk", "adjusted__1__1"),
                                      image = "master_1.0.1", 
                                      tags = list(list(tissue="skin"),
                                                  list(condition="AD"),
                                                  list(project="evo"),
                                                  list(name="analysis",value="X2_V2")),
                                      data_access = "s05")
# generate_ccm -- Mon Nov 17 11:45:22 2025: wf-789d36544f []
>>>>>>> Stashed changes
