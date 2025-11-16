library(cytoreason.cc.client)
devtools::load_all("~/cytoreason.ccm.pipeline/")
library(tidyverse)

# 1. Inputs
# -----------------------------
ccm = as_ccm_fit("wf-08a6a0a503")
signatures = readRDS(get_workflow_outputs("wf-b33649db09"))
## make sure:
## A. removing the extra dataset
## all ligands is in entrez


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
                                                 list(name="analysis",value="X2_V1")),
                                     data_access = "s05")
# generate_ccm -- Fri Nov 14 15:07:26 2025: wf-eafe4014ef []