devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.cc.client)
library(cytoreason.patern)
library(tidyverse)


# Vector prep
# =================================
gxdiff = readRDS(get_workflow_outputs("wf-fe0c7701a0"))
gxdiff <- gxdiff %>%
  mutate(signedP = sign(estimate) * log10_pvalue) %>%
  mutate(estimate_logp = estimate * log10_pvalue)

fc_mat_p = reshape2::dcast(gxdiff, feature_id ~ comparison, value.var = "signedP")
fc_mat_p = column_to_rownames(fc_mat_p, var = "feature_id")
fc_mat_p = as.matrix(fc_mat_p)

insilico = matrix(data = 0, nrow = nrow(fc_mat_p), dimnames = list(rownames(fc_mat_p),"X2"))
insilico["117194",] <- 1

tryptase = matrix(data = 0, nrow = nrow(fc_mat_p), dimnames = list(rownames(fc_mat_p),"mast_tryptase"))
tryptase["7177",] <- 1

alphas <- 0.75

# Propagation in ARCHS
# =================================
# adj_matrix <- readRDS(get_task_outputs('wf-3917825f8b','0'))
adj_matrix = consume_workflow_task_output("wf-3917825f8b", task_id = 0)
tags <- list(list(name="service", value="patern"),
             list(name="network", value="ARCHS"),
             list(name="network_wfid", value="wf-3917825f8b"),
             list(name="alphas", value=paste0(alphas, collapse = '_')),
             list(name="signatures",value="X2"))


# Run Function
# =================================
sapply_dist(X = alphas, function(alpha,
                                 cor_mat,
                                 fc_mat,
                                 edge_cutoff,
                                 edge_cutoff_abs,
                                 stop_val,
                                 direction_split){
  library(cytoreason.patern)
  service_patern(cor_mat = cor_mat,
                 fc_mat = fc_mat,
                 alpha = alpha,
                 edge_cutoff = edge_cutoff,
                 edge_cutoff_abs = edge_cutoff_abs,
                 stop_val = stop_val,
                 direction_split = direction_split)
},
cor_mat = adj_matrix,
fc_mat = fc_mat_p,
# fc_mat = tryptase,
edge_cutoff = NA,
edge_cutoff_abs = FALSE,
stop_val = 1e-03,
direction_split = T,
image = 'eu.gcr.io/cytoreason/ci-cytoreason.patern-package:ori_latest',
cpu_request = '2000m',
memory_request = '25Gi',
tags = tags)

# wf-48a7d142f5 - ARCHS - in silico
# wf-e2d2817c9f - ARCHS - pValue
# wf-785dc86b82 - ARCHS - tryptase


# Useful functions
# =================================
alphas <- 0.75

extract_signatures = function(wfid, gxdiff, extract_top_50 = TRUE){
  smoothed_archs = lapply(get_workflow_outputs(wfid), function(x) { # extract the smoothed matrices
      x = readRDS(x)
      x = x[["pos_smooth_mat"]]
    })
    names(smoothed_archs) = paste0('alpha=', alphas)

    extract_set = function(mat, n_genes) { # extracting top 50/100
      result = list()
      for (col in colnames(mat)) {
        vals <- mat[, col]
        suffix <- paste0("_", n_genes)
        new_name <- paste0(col, suffix)

        if (gxdiff) {
          if (grepl("activated_vs_unactivated", col)) {
            top <- order(vals, decreasing = TRUE)[1:n_genes]
            result[[new_name]] <- list(gene_and_smoothing_value = vals[top], genes = rownames(mat)[top])
          } else if (grepl("inhibited", col)) {
            bottom <- order(vals)[1:n_genes]
            result[[new_name]] <- list(gene_and_smoothing_value = vals[bottom], genes = rownames(mat)[bottom])
          }
        } else {
          top <- order(vals, decreasing = TRUE)[1:n_genes]
          result[[new_name]] <- list(gene_and_smoothing_value = vals[top], genes = rownames(mat)[top])
        }
      }
      return(result)
    }

    result <- lapply(smoothed_archs, function(mat) {
      res <- extract_set(mat, 100)
      if (extract_top_50) {
        res_50 <- extract_set(mat, 50)
        res <- c(res, res_50)
      }
      return(res)
    })

  X2_Smoothed <- list( # reordering list
    gene_and_smoothing_value = unlist(lapply(result, function(mat_result) {
      lapply(mat_result, function(x) x$gene_and_smoothing_value)
    }), recursive = FALSE),

    genes = unlist(lapply(result, function(mat_result) {
      lapply(mat_result, function(x) x$genes)
    }), recursive = FALSE)
  )

}

nameChange = function(X2_Smoothed, suffix, insilico = F) {
  lapply(X2_Smoothed, function(x) {
    x = x[str_detect(names(x),"alpha=0.75")]
    signatureNames = names(x)
    signatureNames = signatureNames %>%
      str_remove("alpha=0.75.") %>%
      str_replace("24hr","late") %>%
      str_replace("4hr","early") %>%
      str_replace("_activated_vs_unactivated","_activation") %>%
      str_replace("_inhibited_vs_uninhibited","_general_inhibition") %>%
      str_replace("activated_vs_unactivated_cov","x2cov_activation") %>%
      str_replace("activated_vs_unactivated","x2_activation") %>%
      str_replace("inhibited_vs_uninhibited_cov","x2cov_general_inhibition") %>%
      str_replace("inhibited_vs_uninhibited","x2_general_inhibition") %>%
      str_replace("inhibited_vs_activated_cov","x2cov_activated_inhibition") %>%
      str_replace("inhibited_vs_activated","x2_activated_inhibition")
    if(!is.null(suffix)) {
      signatureNames = paste0(signatureNames,"_",suffix)
    }
    if(insilico) {
      signatureNames = str_replace(signatureNames, "X2_","x2_insilico_")
    }
    names(x) = signatureNames
    return(x)
  })
}


# 1. ARCHS - gxdiff
# --------------------------------
X2_Smoothed = extract_signatures("wf-e2d2817c9f", gxdiff = T, extract_top_50 = T)
pushToCC(X2_Smoothed, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_withAlphas")))
# wf-c547cb8388

x2_archs = nameChange(X2_Smoothed, "archs")
pushToCC(x2_archs, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_alpha0.75"),
                                     list(name="network",value="archs"),
                                     list(name="origin",value="gx_diff")))
# wf-d9a757211a


# 2. ARCHS - in silico
# ---------------------------------------------
X2_Smoothed = extract_signatures("wf-48a7d142f5", gxdiff = F, extract_top_50 = T)
pushToCC(X2_Smoothed, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_withAlphas")))
# wf-54c7f42333

x2_archs_insilico = nameChange(X2_Smoothed, "archs", insilico = T)
pushToCC(x2_archs_insilico, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_alpha0.75"),
                                               list(name="network",value="archs"),
                                               list(name="origin",value="in_silico")))
# wf-cc74409fc1


# 3. ARCHS - tryptase
# ---------------------------------------------
X2_Smoothed = extract_signatures("wf-785dc86b82", gxdiff = F, extract_top_50 = T)
pushToCC(X2_Smoothed, tagsToPass = list(list(name = "object",value="Tryptase_smoothed_signatures_withAlphas")))
# wf-f8c8f62b97

tryptase = nameChange(X2_Smoothed, NULL, insilico = T)
pushToCC(tryptase, tagsToPass = list(list(name = "object",value="tryptase_smoothed_signatures_alpha0.75"),
                                              list(name="network",value="archs"),
                                              list(name="origin",value="in_silico")))
# wf-acea788f44


# 5. Append signatures to BQ
# -------------------------------------------
append_signatures = function(wfid, rankingMetric){
  X2_Signatures = readRDS(get_workflow_outputs(wfid))
  signature_df = lapply(names(X2_Signatures$gene_and_smoothing_value), function(signature){
    x = X2_Signatures$gene_and_smoothing_value[[signature]]
    x = enframe(x, name = "feature_id", value = "estimate")
    agonist = str_split(signature,"_")[[1]][1]
    if(agonist == "x2") { agonist = "All" }
    x$agonist = agonist
    x$signature = signature
    return(x)
  }) %>% bind_rows()
  signature_df$rankingMetric = rankingMetric
  new_cols <- setNames(rep(list(NA), 4), c("comparison", "log10_fdr", "log10_pvalue", "in50"))
  
  signature_df = signature_df %>%
    mutate(!!!new_cols) %>%
    mutate(wfid = wfid)
}

archs = append_signatures("wf-d9a757211a","pvalue")
insilico_archs = append_signatures("wf-6028a1f83b","insilico")
tryptase = append_signatures("wf-acea788f44","tryptase")

all = rbind(archs, insilico_archs, tryptase)

uploadToBQ(all, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures", disposition = "WRITE_APPEND")



# 6. Add in-silico from string
# -------------------------------------
x2_fromString <- cytoreason.assets::read_asset("ccw://wf-9329a6bc38@0:signatures_18803")
x2_fromString = x2_fromString[["MRGPRX2"]]
x2_fromString = data.frame(feature_id = x2_fromString, estimate = NA, agonist = NA, signature = "MRGPRX2", comparison = NA, log10_fdr = NA,
                           log10_pvalue = NA, in50 = NA, wfid = "wf-9329a6bc38", rankingMetric = "insilico_Shiran")
uploadToBQ(x2_fromString, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures", disposition = "WRITE_APPEND")
