library(cytoreason.cc.client)
library(tidyverse)
library(purrr)
devtools::load_all("~/analysis-s05/R/utils.R")

### Prep
### ===========================
gxdiff = readRDS(get_workflow_outputs("wf-fe0c7701a0"))
gxdiff <- gxdiff %>%
  mutate(signedP = sign(estimate) * log10_pvalue) %>%
  mutate(estimate_logp = estimate * log10_pvalue) %>%
  mutate(agonist = str_replace(agonist, "All","Pooled")) %>%
  mutate(agonist = ifelse(str_detect(comparison, "_cov"), "Pooled+cov",agonist))
gxdiff$time = sapply(gxdiff$comparison, function(x) tail(strsplit(x, "_")[[1]], 1))
gxdiff$agonist = paste0(gxdiff$agonist, "_", gxdiff$time)

ranking_methods <- list(
  pvalue = list(col = "signedP", suffix = "_p")
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
generate_signatures <- function(activation_list, general_inhibition_list, activated_inhibition_list, agonists, ranking_methods) {
  genes_and_estimates <- list()
  genes_only <- list()
  
  for (method_name in names(ranking_methods)) {
    # suffix <- ranking_methods[[method_name]]$suffix
    suffix = NULL
    act_df <- activation_list[[method_name]]
    gen_inh_df <- general_inhibition_list[[method_name]]
    act_inh_df <- activated_inhibition_list[[method_name]]
    
    for (agonist in agonists) {
      # Activation
      genes_and_estimates[[paste0(agonist, "_activation_50", suffix)]] <- deframe(get_signature(act_df, agonist, "in_top50"))
      genes_and_estimates[[paste0(agonist, "_activation_100", suffix)]] <- deframe(get_signature(act_df, agonist))
      genes_only[[paste0(agonist, "_activation_50", suffix)]] <- get_signature(act_df, agonist, "in_top50")$feature_id
      genes_only[[paste0(agonist, "_activation_100", suffix)]] <- get_signature(act_df, agonist)$feature_id
      
      # General Inhibition
      genes_and_estimates[[paste0(agonist, "_general_inhibition_50", suffix)]] <- deframe(get_signature(gen_inh_df, agonist, "in_bottom50"))
      genes_and_estimates[[paste0(agonist, "_general_inhibition_100", suffix)]] <- deframe(get_signature(gen_inh_df, agonist))
      genes_only[[paste0(agonist, "_general_inhibition_50", suffix)]] <- get_signature(gen_inh_df, agonist, "in_bottom50")$feature_id
      genes_only[[paste0(agonist, "_general_inhibition_100", suffix)]] <- get_signature(gen_inh_df, agonist)$feature_id
      
      # Activated Inhibition
      genes_and_estimates[[paste0(agonist, "_activated_inhibition_50", suffix)]] <- deframe(get_signature(act_inh_df, agonist, "in_bottom50"))
      genes_and_estimates[[paste0(agonist, "_activated_inhibition_100", suffix)]] <- deframe(get_signature(act_inh_df, agonist))
      genes_only[[paste0(agonist, "_activated_inhibition_50", suffix)]] <- get_signature(act_inh_df, agonist, "in_bottom50")$feature_id
      genes_only[[paste0(agonist, "_activated_inhibition_100", suffix)]] <- get_signature(act_inh_df, agonist)$feature_id
    }
  }
  
  list(genes_and_estimates = genes_and_estimates, genes_only = genes_only)
}


### Generate signatures 
### ==============================================
# Generate all activation/inhibition datasets
# ---------------------------------------------
activation_list <- map(ranking_methods, ~get_top_genes(gxdiff, "activated_vs_unactivated", "descending", .x$col, "top"))
general_inhibition_list <- map(ranking_methods, ~get_top_genes(gxdiff, "inhibited_vs_uninhibited", "ascending", .x$col, "bottom"))
activated_inhibition_list <- map(ranking_methods, ~get_top_genes(gxdiff, "inhibited_vs_activated", "ascending", .x$col, "bottom"))

activation <- dplyr::bind_rows(activation_list, .id = "rankingMetric") %>%
  dplyr::rename(in50 = in_top50)
general_inhibition <- dplyr::bind_rows(general_inhibition_list, .id = "rankingMetric") %>%
  dplyr::rename(in50 = in_bottom50)
activated_inhibition <- dplyr::bind_rows(activated_inhibition_list, .id = "rankingMetric") %>%
  dplyr::rename(in50 = in_bottom50)
allTopGenes <- rbind(activation, general_inhibition, activated_inhibition)
allTopGenes$comparison = str_remove(allTopGenes$comparison, "_cov")
# wf-6e5b903b12

X2_Signatures <- generate_signatures(activation_list, general_inhibition_list, activated_inhibition_list, agonists, ranking_methods)

rename_signature_keys <- function(sig_list) {
  new_names <- names(sig_list) %>%
    gsub("_4hr", "_early", .) %>%
    gsub("_24hr", "_late", .) %>%
    gsub("^Pooled", "x2", .) %>%
    gsub("^(x2|x2\\+cov|SP|CST14|Icatibant|PAMP12|Untreated|aIgE)_(early|late)?_(activation|general_inhibition|activated_inhibition)_(\\d+)(.*)$",
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
# wf-45c79e4b82

signature_df = lapply(X2_Signatures$genes_and_estimates, function(x){
  enframe(x, name = "feature_id", value = "estimate")
}) %>% bind_rows(.id = "signature") %>%
  mutate(rankingMetric = "pvalue") %>%
  mutate(wfid = "wf-45c79e4b82")
signature_df$agonist = sapply(signature_df$signature, function(x) {
  pre = str_split(x,"_")[[1]][1]
  post = ifelse(str_detect(x,"late"),"_24hr","_4hr")
  return(paste0(pre,post))
})
allTopGenes$agonist = str_replace(allTopGenes$agonist, "Pooled","x2")

signatures = merge(signature_df, allTopGenes, by.x = c("agonist","feature_id", "estimate","rankingMetric"), 
                   by.y = c("agonist","feature_id","estimate","rankingMetric"), all = T)
  signatures = signatures[-which(str_detect(signatures$signature,"general_inhibition") & str_detect(signatures$comparison, "inhibited_vs_activated")),]
  signatures = signatures[-which(str_detect(signatures$signature,"activated_inhibition|activate_inhibition") & 
                                   str_detect(signatures$comparison, "inhibited_vs_uninhibited")),]
signatures$signature = str_replace(signatures$signature,"\\+","")
uploadToBQ(signatures, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures")
pushToCC(signatures, tagsToPass = list(list(name="object",value="X2Signatures_all_top_genes")))
# wf-932d6ea274