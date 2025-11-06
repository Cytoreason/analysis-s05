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


# 1. Top 50/100
# -----------------------------
x2_raw = readRDS(get_workflow_outputs("wf-ac3be6cf37"))$genes_only
x2_inNet = refine_signatures(x2_raw)

pushToCC(x2_inNet, tagsToPass = list(list(name = "object",value="X2_refined"),
                                      list(name="network_propagation",value="none"),
                                      list(name="origin",value="top")))
# wf-3fe1dce17d


# 2. Propagated in ARCHS
# -----------------------------
x2_archs = readRDS(get_workflow_outputs("wf-f5930577e9"))$genes
x2_inNet = refine_signatures(x2_archs)

pushToCC(x2_inNet, tagsToPass = list(list(name = "object",value="X2_refined"),
                                     list(name="network_propagation",value="archs"),
                                     list(name="origin",value="gxdiff")))
# wf-3b31fba621


# 3. Propagated in STRING
# -----------------------------
x2_string = readRDS(get_workflow_outputs("wf-ad41f5ce55"))$genes
x2_inNet = refine_signatures(x2_string)

pushToCC(x2_inNet, tagsToPass = list(list(name = "object",value="X2_refined"),
                                     list(name="network_propagation",value="string"),
                                     list(name="origin",value="gxdiff")))
# wf-6fff51f9f7


# 4. in-silico, propagated in ARCHS
# -----------------------------
x2_insilico_archs = readRDS(get_workflow_outputs("wf-6028a1f83b"))$genes
x2_inNet = refine_signatures(x2_insilico_archs)

pushToCC(x2_inNet, tagsToPass = list(list(name = "object",value="X2_refined"),
                                     list(name="network_propagation",value="archs"),
                                     list(name="origin",value="insilico")))
# wf-7305f8adfa


# 5. in-silico, propagated in STRING
# -----------------------------
x2_insilico_string = readRDS(get_workflow_outputs("wf-b6999b00da"))$genes
x2_inNet = refine_signatures(x2_insilico_string)

pushToCC(x2_inNet, tagsToPass = list(list(name = "object",value="X2_refined"),
                                     list(name="network_propagation",value="string"),
                                     list(name="origin",value="insilico")))
# wf-843f67359a





# 5. Append signatures to BQ
# ------------------------------------
append_signatures = function(wfid){
  X2_Signatures = readRDS(get_workflow_outputs(wfid))
  signature_df = lapply(names(X2_Signatures), function(signature){
    x = X2_Signatures[[signature]]
    x = enframe(x, name = "agonist", value = "feature_id")
    agonist = str_split(signature,"_")[[1]][1]
    if(agonist == "x2") { agonist = "All" }
    x$agonist = agonist
    x$signature = signature
    return(x)
  }) %>% bind_rows()
  
  new_cols <- setNames(rep(list(NA), 6), c("comparison", "fdr", "pvalue", "in_top50", "in_bottom50","estimate"))
  
  signature_df = signature_df %>%
    mutate(!!!new_cols) %>%
    mutate(gx_diff_wfid = wfid)
}
archs_refined = append_signatures("wf-3b31fba621")

uploadToBQ(archs_refined, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures", disposition = "WRITE_APPEND")
