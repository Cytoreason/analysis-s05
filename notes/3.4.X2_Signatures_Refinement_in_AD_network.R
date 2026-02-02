devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.cc.client)
library(tidyverse)

bulk_ad_network = readRDS(get_workflow_outputs("wf-06da2e64ae", files_names_grepl_pattern = "igraph_NW_sub_0.6.rds")) # ad_nw_bulk_pruned_0.6
bulk_ad_network = igraph::as_data_frame(bulk_ad_network)
allgenes_network = c(bulk_ad_network$from, bulk_ad_network$to)
allgenes_network = data.frame(table(allgenes_network))
colnames(allgenes_network)[1] = "gene"
allgenes_network$gene = as.character(allgenes_network$gene)

refine_signatures = function(signatures) {
  signatures = lapply(signatures, function(x2){
    allgenes_network$gene[which(allgenes_network$gene %in% x2)]
  })
  names(signatures) = paste0(names(signatures),"_refined")
  return(signatures)
}

allSignatures = bq_table_download(x = bq_table(project = "cytoreason", dataset = "s05_atopic_dermatitis", table="X2Signatures"))
allSignatures = allSignatures[!str_detect(allSignatures$signature,"mast_tryptase"),]
allSignatures = split(allSignatures$feature_id, allSignatures$signature)                                  

x2_inNet = refine_signatures(allSignatures)
pushToCC(x2_inNet, tagsToPass = list(list(name = "object",value="X2_refined")))
# wf-3495fcab6f


# Append signatures to BQ
# ------------------------------------
append_signatures = function(wfid, rankingMetric){
  X2_Signatures = readRDS(get_workflow_outputs(wfid))
  signature_df = lapply(names(X2_Signatures), function(signature){
    x = X2_Signatures[[signature]]
    x = enframe(x, name = "agonist", value = "feature_id")
    agonist = str_split(signature,"_")[[1]][1]
    if(agonist == "x2") { agonist = "All" }
    if(agonist == "x2+cov") { agonist = "All+cov" }
    x$agonist = agonist
    x$signature = signature
    return(x)
  }) %>% bind_rows()
  signature_df$rankingMetric = rankingMetric
  new_cols <- setNames(rep(list(NA), 5), c("comparison", "log10_fdr", "log10_pvalue", "in50","estimate"))
  
  signature_df = signature_df %>%
    mutate(!!!new_cols) %>%
    mutate(wfid = wfid)
}

x2_refined = append_signatures("wf-3495fcab6f", "refined")

uploadToBQ(x2_refined, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures", disposition = "WRITE_APPEND")

