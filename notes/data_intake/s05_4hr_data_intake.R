library(bigrquery)
devtools::load_all("~/code/cytoreason.ccm.pipeline")
library("cytoreason.validator.apps.client")
# library(cytoreason.deconvolution)
# library(cytoreason.xdesign)
library(dplyr)
library(ggplot2)

#####################################################################################################################################
#####     This code includes all pre-processing of Evommune 4hr experiment from validator stage until CCM services analysis     #####
#####################################################################################################################################

retrieve_me <- function(my_object = my_object ){
  return(my_object)
}

###############################################
#                   Import                    #
###############################################

# Pull data from validator
s05_4hr_raw_eset = get_eset_from_dataset("s05_X2_RNAseq_4hr_single_2SlaF", "rnaseq")
experimentID(s05_4hr_raw_eset) <- "s05_4hr"
s05_4hr_raw_pData = pData(s05_4hr_raw_eset)
cytoreason.cc.client::run_function_dist(retrieve_me, my_object = s05_4hr_raw_eset) # wf-cb183f8dad # SUCCEEDED
s05_4hr_raw_eset = readRDS(get_workflow_outputs("wf-cb183f8dad"))


###############################################
#             Meta-data handling              #
###############################################


###############################################
#                     QC                      #
###############################################

# Generate cell contribution eset for cells PCA
s05_4hr_ct_service = ccm_service_cell_contribution(expression_data = s05_4hr_raw_eset, signature_collection = "skin_v12")
s05_4hr_ct_eset = s05_4hr_ct_service$eset
cytoreason.cc.client::run_function_dist(retrieve_me, my_object = s05_4hr_ct_eset) # wf-34a4a0fa5d # SUCCEEDED

# Train the cells PCA
cells_PCA <- cytoreason.analytics::service_pca(object = s05_4hr_ct_eset, n_pc = 8)
cytoreason.cc.client::run_function_dist(retrieve_me, my_object = cells_PCA) # wf-2964588404 # SUCCEEDED

# Extract cells PCA sample loadings
s05_4hr_sample_loadings_cells_PCA = statistic_table(cells_PCA)$sample_loadings %>% 
  tidyr::pivot_wider(names_from = pc, values_from = value)
  # dplyr::mutate(Sample=sample_id) %>% 
  # left_join(pData(s05_4hr_raw_eset), by = "Sample")
bq_job <- bq_perform_upload(bq_table("cytoreason","s05_atopic_dermatitis",table="s05_4hr_cells_PCA_sample_loadings"),
                  s05_4hr_sample_loadings_cells_PCA,write_disposition = "WRITE_TRUNCATE")
bq_job_status(bq_job)

# Top 1000 variable features PCA
# Find top 1000 variable features
expr_matrix = exprs(s05_4hr_raw_eset)
gene_variability <- apply(expr_matrix, 1, var)
top_genes <- names(sort(gene_variability, decreasing = TRUE))[1:1000]

# Filter eset for top variable genes
s05_4hr_bulk_eset_top1000_var_genes = s05_4hr_raw_eset[top_genes,]

# Train the top 1000 features PCA
features_top1000_PCA <- cytoreason.analytics::service_pca(object = s05_4hr_bulk_eset_top1000_var_genes, n_pc = 8)

# Extract top 1000 features PCA sample loadings
s05_4hr_sample_loadings_top1000_features_PCA = statistic_table(features_top1000_PCA)$sample_loadings %>% 
  tidyr::pivot_wider(names_from = pc, values_from = value)
  # dplyr::mutate(Sample=sample_id) %>% 
  # left_join(pData(s05_4hr_raw_eset), by = "Sample")
bq_job <- bq_perform_upload(bq_table("cytoreason","s05_atopic_dermatitis",table="s05_4hr_top1000_features_PCA_sample_loadings"),
                  s05_4hr_sample_loadings_top1000_features_PCA,write_disposition = "WRITE_TRUNCATE")
bq_job_status(bq_job)
