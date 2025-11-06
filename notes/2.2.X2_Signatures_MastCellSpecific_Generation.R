library(cytoreason.ccm.pipeline)
library(cytoreason.cc.client)
library(tidyverse)
library(reshape2)
devtools::load_all("~/analysis-s05/R/utils.R")


# Extract top 50/100 based on the estimate
# --------------------------------------------------
gxdiff = readRDS(get_workflow_outputs("wf-edcade5c80"))

activation = gxdiff %>%
  dplyr::filter(term == "activated_vs_unactivated") %>%
  # dplyr::filter(agonist %in% c("All_4hr","All_24hr")) %>%
  group_by(agonist) %>%
  dplyr::slice_max(order_by = estimate, n = 100) %>%
  mutate(in_top50 = rank(-estimate) <= 50) %>%
  dplyr::select(agonist, comparison, feature_id, estimate, pvalue, fdr, in_top50) %>%
  ungroup()

inhibition <- gxdiff %>%
  dplyr::filter(term == "inhibited_vs_uninhibited") %>%
  # dplyr::filter(agonist %in% c("All_4hr", "All_24hr")) %>%
  group_by(agonist) %>%
  slice_min(order_by = estimate, n = 100) %>%
  mutate(in_bottom50 = rank(estimate) <= 50) %>%
  dplyr::select(agonist, comparison, feature_id, estimate, pvalue, fdr, in_bottom50) %>%
  ungroup()

allTopGenes = bind_rows(activation, inhibition) %>%
  mutate(gx_diff_wfid = "wf-edcade5c80")


# Generate signatures
# ------------------------------------------
get_signature <- function(df, agonist_value, filter_col = NULL) {
  df_filtered <- df %>%
    dplyr::filter(agonist == agonist_value)
  
  if (!is.null(filter_col)) {
    df_filtered <- df_filtered %>%
      dplyr::filter(!!sym(filter_col))
  }
  
  df_filtered %>%
    dplyr::select(feature_id, estimate)
}


X2_Signatures <- list(
  genes_and_estimates = list(
    x2_activation_early_50 = deframe(get_signature(activation, "All_4hr", "in_top50")),
    x2_activation_late_50 = deframe(get_signature(activation, "All_24hr", "in_top50")),
    x2_activation_early_100 = deframe(get_signature(activation, "All_4hr")),
    x2_activation_late_100 = deframe(get_signature(activation, "All_24hr")),
    x2_inhibition_early_50 = deframe(get_signature(inhibition, "All_4hr", "in_bottom50")),
    x2_inhibition_late_50 = deframe(get_signature(inhibition, "All_24hr", "in_bottom50")),
    x2_inhibition_early_100 = deframe(get_signature(inhibition, "All_4hr")),
    x2_inhibition_late_100 = deframe(get_signature(inhibition, "All_24hr")),
    
    SP_activation_early_50 = deframe(get_signature(activation, "SP_4hr", "in_top50")),
    SP_activation_late_50 = deframe(get_signature(activation, "SP_24hr", "in_top50")),
    SP_activation_early_100 = deframe(get_signature(activation, "SP_4hr")),
    SP_activation_late_100 = deframe(get_signature(activation, "SP_24hr")),
    SP_inhibition_early_50 = deframe(get_signature(inhibition, "SP_4hr", "in_bottom50")),
    SP_inhibition_late_50 = deframe(get_signature(inhibition, "SP_24hr", "in_bottom50")),
    SP_inhibition_early_100 = deframe(get_signature(inhibition, "SP_4hr")),
    SP_inhibition_late_100 = deframe(get_signature(inhibition, "SP_24hr")),

    CST14_activation_early_50 = deframe(get_signature(activation, "CST14", "in_top50")),
    CST14_activation_early_100 = deframe(get_signature(activation, "CST14")),
    CST14_inhibition_early_50 = deframe(get_signature(inhibition, "CST14", "in_bottom50")),
    CST14_inhibition_early_100 = deframe(get_signature(inhibition, "CST14")),
    
    Icatibant_activation_early_50 = deframe(get_signature(activation, "Icatibant", "in_top50")),
    Icatibant_activation_early_100 = deframe(get_signature(activation, "Icatibant")),
    Icatibant_inhibition_early_50 = deframe(get_signature(inhibition, "Icatibant", "in_bottom50")),
    Icatibant_inhibition_early_100 = deframe(get_signature(inhibition, "Icatibant")),
    
    PAMP12_activation_early_50 = deframe(get_signature(activation, "PAMP12", "in_top50")),
    PAMP12_activation_early_100 = deframe(get_signature(activation, "PAMP12")),
    PAMP12_inhibition_early_50 = deframe(get_signature(inhibition, "PAMP12", "in_bottom50")),
    PAMP12_inhibition_early_100 = deframe(get_signature(inhibition, "PAMP12")),
    
    Untreated_inhibition_early_50 = deframe(get_signature(inhibition, "Untreated", "in_bottom50")),
    Untreated_inhibition_early_100 = deframe(get_signature(inhibition, "Untreated")),
    
    aIgE_activation_late_50 = deframe(get_signature(activation, "aIgE", "in_top50")),
    aIgE_activation_late_100 = deframe(get_signature(activation, "aIgE")),
    aIgE_inhibition_late_50 = deframe(get_signature(inhibition, "aIgE", "in_bottom50")),
    aIgE_inhibition_late_100 = deframe(get_signature(inhibition, "aIgE"))
  ),
  
  genes_only = list(
    x2_activation_early_50 = get_signature(activation, "All_4hr", "in_top50")$feature_id,
    x2_activation_late_50 = get_signature(activation, "All_24hr", "in_top50")$feature_id,
    x2_activation_early_100 = get_signature(activation, "All_4hr")$feature_id,
    x2_activation_late_100 = get_signature(activation, "All_24hr")$feature_id,
    x2_inhibition_early_50 = get_signature(inhibition, "All_4hr", "in_bottom50")$feature_id,
    x2_inhibition_late_50 = get_signature(inhibition, "All_24hr", "in_bottom50")$feature_id,
    x2_inhibition_early_100 = get_signature(inhibition, "All_4hr")$feature_id,
    x2_inhibition_late_100 = get_signature(inhibition, "All_24hr")$feature_id,
    
    SP_activation_early_50 = get_signature(activation, "SP_4hr", "in_top50")$feature_id,
    SP_activation_late_50 = get_signature(activation, "SP_24hr", "in_top50")$feature_id,
    SP_activation_early_100 = get_signature(activation, "SP_4hr")$feature_id,
    SP_activation_late_100 = get_signature(activation, "SP_24hr")$feature_id,
    SP_inhibition_early_50 = get_signature(inhibition, "SP_4hr", "in_bottom50")$feature_id,
    SP_inhibition_late_50 = get_signature(inhibition, "SP_24hr", "in_bottom50")$feature_id,
    SP_inhibition_early_100 = get_signature(inhibition, "SP_4hr")$feature_id,
    SP_inhibition_late_100 = get_signature(inhibition, "SP_24hr")$feature_id,
    
    CST14_activation_early_50 = get_signature(activation, "CST14", "in_top50")$feature_id,
    CST14_activation_early_100 = get_signature(activation, "CST14")$feature_id,
    CST14_inhibition_early_50 = get_signature(inhibition, "CST14", "in_bottom50")$feature_id,
    CST14_inhibition_early_100 = get_signature(inhibition, "CST14")$feature_id,
    
    Icatibant_activation_early_50 = get_signature(activation, "Icatibant", "in_top50")$feature_id,
    Icatibant_activation_early_100 = get_signature(activation, "Icatibant")$feature_id,
    Icatibant_inhibition_early_50 = get_signature(inhibition, "Icatibant", "in_bottom50")$feature_id,
    Icatibant_inhibition_early_100 = get_signature(inhibition, "Icatibant")$feature_id,
    
    PAMP12_activation_early_50 = get_signature(activation, "PAMP12", "in_top50")$feature_id,
    PAMP12_activation_early_100 = get_signature(activation, "PAMP12")$feature_id,
    PAMP12_inhibition_early_50 = get_signature(inhibition, "PAMP12", "in_bottom50")$feature_id,
    PAMP12_inhibition_early_100 = get_signature(inhibition, "PAMP12")$feature_id,
    
    Untreated_inhibition_early_50 = get_signature(inhibition, "Untreated", "in_bottom50")$feature_id,
    Untreated_inhibition_early_100 = get_signature(inhibition, "Untreated")$feature_id,
    
    aIgE_activation_late_50 = get_signature(activation, "aIgE", "in_top50")$feature_id,
    aIgE_activation_late_100 = get_signature(activation, "aIgE")$feature_id,
    aIgE_inhibition_late_50 = get_signature(inhibition, "aIgE", "in_bottom50")$feature_id,
    aIgE_inhibition_late_100 = get_signature(inhibition, "aIgE")$feature_id
  )
)
pushToCC(X2_Signatures, tagsToPass = list(list(name="object",value="X2Signatures_MCS")))
# wf-ba31d91a68

signature_df = lapply(X2_Signatures$genes_and_estimates, function(x){
  enframe(x, name = "feature_id", value = "estimate")
}) %>% bind_rows(.id = "signature")

signatures = merge(signature_df, allTopGenes, by.x = c("feature_id", "estimate"), by.y = c("feature_id","estimate"))
uploadToBQ(signatures, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures")
pushToCC(signatures, tagsToPass = list(list(name="object",value="X2Signatures_all_top_genes")))
# wf-99209e4581


## Hypergeometric test
## ======================================
library(cytoreason.gx)

# 1. use all genes as background
# --------------------------------
background_genes = unique(gxdiff$feature_id)

hypergeometric = lapply(X2_Signatures$genes_and_estimates, function(genes){
  data.frame(statistic_table(gx_gsa(
    x = genes, method = "hypergeometric", nPerm = 10000, background = background_genes)))
})

hypergeometric = bind_rows(hypergeometric, .id = "signature")
hypergeometric$background_genes = "all genes"
uploadToBQ(hypergeometric, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures_hypergeometric")
pushToCC(hypergeometric, tagsToPass = list(list(name="object",value="X2Signatures_hypergeometric_bg_allgenes")))
# wf-abc1528954

# 2. background genes as the up/downregulated genes
# ---------------------------------------------------
mapSignatureToGXDIFF = unique(signatures[,c("signature","comparison","agonist")]) %>%
  mutate(term = ifelse(str_detect(comparison, "activated"),"activated_vs_unactivated","inhibited_vs_uninhibited"))
mapSignatureToGXDIFF = merge(mapSignatureToGXDIFF, gxdiff, by = c("comparison","agonist","term"), all = T)

hypergeometric = lapply(names(X2_Signatures$genes_and_estimates), function(signature){
  genes = X2_Signatures$genes_and_estimates[[signature]]
  if(str_detect(signature, "activation")) { # for activation we take the upregulated genes
    background_genes = mapSignatureToGXDIFF$feature_id[which(mapSignatureToGXDIFF$signature == signature & mapSignatureToGXDIFF$estimate > 0)]
    dir = "upregulated genes"
  } else if(str_detect(signature, "inhibition")) { # for inhibition we take the downregulated genes
    background_genes = mapSignatureToGXDIFF$feature_id[which(mapSignatureToGXDIFF$signature == signature & mapSignatureToGXDIFF$estimate < 0)]
    dir = "downregulated genes"
  }
  data.frame(statistic_table(gx_gsa(
    x = genes, method = "hypergeometric", nPerm = 10000, background = background_genes)), background_genes = dir)
})
names(hypergeometric) = names(X2_Signatures$genes_and_estimates)
hypergeometric = bind_rows(hypergeometric, .id = "signature")

uploadToBQ(hypergeometric, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures_hypergeometric", disposition = "WRITE_APPEND")
pushToCC(hypergeometric, tagsToPass = list(list(name="object",value="X2Signatures_hypergeometric_bg_up/down regulated")))
# wf-b7fddcf37a