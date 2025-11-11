library(cytoreason.cc.client)
library(tidyverse)
library(purrr)
devtools::load_all("~/analysis-s05/R/utils.R")

### Prep
### ===========================
gxdiff = readRDS(get_workflow_outputs("wf-edcade5c80"))
gxdiff <- gxdiff %>%
  mutate(signedP = sign(estimate) * log10_pvalue) %>%
  mutate(estimate_logp = estimate * log10_pvalue) %>%
  rowwise() %>%
  mutate(agonist = case_when(
    !agonist %in% c("All_4hr", "All_24hr", "SP_4hr", "SP_24hr") ~ 
      paste0(agonist, "_", tail(strsplit(comparison, "_")[[1]], 1)),
    TRUE ~ agonist
  )) %>%
  ungroup()


ranking_methods <- list(
  estimate = list(col = "estimate", suffix = ""),
  pvalue = list(col = "signedP", suffix = "_p"),
  estimate_logp = list(col = "estimate_logp", suffix = "_ep")
)

agonists <- unique(gxdiff$agonist)

### Helper functions
### ==================================
get_top_genes <- function(df, term, direction, rank_by, label_prefix) {
  label_col <- paste0("in_", label_prefix, "50")
  df = df[which(df$term == term), ]
  
  if (label_prefix == "top") {
    df <- dplyr::filter(df, estimate > 0)
  } else if (label_prefix == "bottom") {
    df <- dplyr::filter(df, estimate < 0)
  }
  
  df <- df %>%
    group_by(agonist) %>%
    group_modify(~ {
      if (direction == "descending") {
        ranked_df <- dplyr::slice_max(.x, order_by = .x[[rank_by]], n = 100)
        ranked_df[[label_col]] <- rank(-ranked_df[[rank_by]]) <= 50
      } else {
        ranked_df <- dplyr::slice_min(.x, order_by = .x[[rank_by]], n = 100)
        ranked_df[[label_col]] <- rank(ranked_df[[rank_by]]) <= 50
      }
      
      ranked_df
    }) %>%
    ungroup()
  
  dplyr::select(df, agonist, comparison, feature_id, estimate, log10_pvalue, log10_fdr, dplyr::all_of(label_col))
}

# Signature extraction
get_signature <- function(df, agonist_value, filter_col = NULL) {
  df_filtered <- df[df$agonist == agonist_value, ]
  if (!is.null(filter_col)) {
    df_filtered <- df_filtered[df_filtered[[filter_col]], ]
  }
  df_filtered[, c("feature_id", "estimate")]
}


# Generate signatures
generate_signatures <- function(activation_list, inhibition_list, agonists, ranking_methods) {
  genes_and_estimates <- list()
  genes_only <- list()
  
  for (method_name in names(ranking_methods)) {
    suffix <- ranking_methods[[method_name]]$suffix
    act_df <- activation_list[[method_name]]
    inh_df <- inhibition_list[[method_name]]
    
    for (agonist in agonists) {
      # Activation
      genes_and_estimates[[paste0(agonist, "_activation_50", suffix)]] <- deframe(get_signature(act_df, agonist, "in_top50"))
      genes_and_estimates[[paste0(agonist, "_activation_100", suffix)]] <- deframe(get_signature(act_df, agonist))
      genes_only[[paste0(agonist, "_activation_50", suffix)]] <- get_signature(act_df, agonist, "in_top50")$feature_id
      genes_only[[paste0(agonist, "_activation_100", suffix)]] <- get_signature(act_df, agonist)$feature_id
      
      # Inhibition
      genes_and_estimates[[paste0(agonist, "_inhibition_50", suffix)]] <- deframe(get_signature(inh_df, agonist, "in_bottom50"))
      genes_and_estimates[[paste0(agonist, "_inhibition_100", suffix)]] <- deframe(get_signature(inh_df, agonist))
      genes_only[[paste0(agonist, "_inhibition_50", suffix)]] <- get_signature(inh_df, agonist, "in_bottom50")$feature_id
      genes_only[[paste0(agonist, "_inhibition_100", suffix)]] <- get_signature(inh_df, agonist)$feature_id
    }
  }
  
  list(genes_and_estimates = genes_and_estimates, genes_only = genes_only)
}


### Generate signatures 
### ==============================================
# Generate all activation/inhibition datasets
# ---------------------------------------------
activation_list <- map(ranking_methods, ~get_top_genes(gxdiff, "activated_vs_unactivated", "descending", .x$col, "top"))
inhibition_list <- map(ranking_methods, ~get_top_genes(gxdiff, "inhibited_vs_uninhibited", "ascending", .x$col, "bottom"))
activation <- dplyr::bind_rows(activation_list, .id = "rankingMetric") %>%
  rename(in50 = in_top50)
inhibition <- dplyr::bind_rows(inhibition_list, .id = "rankingMetric") %>%
  rename(in50 = in_bottom50)
allTopGenes <- rbind(activation, inhibition)
                          
X2_Signatures <- generate_signatures(activation_list, inhibition_list, agonists, ranking_methods)

rename_signature_keys <- function(sig_list) {
  new_names <- names(sig_list) %>%
    gsub("_4hr", "_early", .) %>%
    gsub("_24hr", "_late", .) %>%
    gsub("^All", "x2", .) %>%
    gsub("^(x2|SP|CST14|Icatibant|PAMP12|Untreated|aIgE)_(early|late)?_(activation|inhibition)_(\\d+)(.*)$",
         "\\1_\\3_\\2_\\4\\5", .)
  
  
  names(sig_list) <- new_names
  sig_list
}

X2_Signatures$genes_and_estimates <- rename_signature_keys(X2_Signatures$genes_and_estimates)
X2_Signatures$genes_only <- rename_signature_keys(X2_Signatures$genes_only)

# remove untreated activation
X2_Signatures <- lapply(X2_Signatures, function(x) x[lengths(x) > 0])


# Final signatures
# -----------------------
pushToCC(X2_Signatures, tagsToPass = list(list(name="object",value="X2Signatures")))
# wf-712302a4c4

signature_df = lapply(X2_Signatures$genes_and_estimates, function(x){
  enframe(x, name = "feature_id", value = "estimate")
}) %>% bind_rows(.id = "signature") %>%
  mutate(rankingMetric = case_when(str_detect(signature, "_p") ~ "pvalue",
                                   str_detect(signature, "_ep") ~ "estimate_logp",
                                   .default = "estimate")) %>%
  mutate(wfid = "wf-712302a4c4")

signatures = merge(signature_df, allTopGenes, by.x = c("feature_id", "estimate","rankingMetric"), by.y = c("feature_id","estimate","rankingMetric"))
uploadToBQ(signatures, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures")
pushToCC(signatures, tagsToPass = list(list(name="object",value="X2Signatures_all_top_genes")))
# wf-5b45a8e1bc

