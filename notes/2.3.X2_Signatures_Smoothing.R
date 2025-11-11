devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.cc.client)
library(cytoreason.patern)
library(tidyverse)


# Vector prep
# =================================
gxdiff = readRDS(get_workflow_outputs("wf-edcade5c80"))
gxdiff <- gxdiff %>%
  mutate(signedP = sign(estimate) * log10_pvalue) %>%
  mutate(estimate_logp = estimate * log10_pvalue)

fc_mat = reshape2::dcast(gxdiff[-which(gxdiff$term == "unactivated"),],
                         feature_id ~ comparison, value.var = "estimate")
fc_mat = column_to_rownames(fc_mat, var = "feature_id")
fc_mat = as.matrix(fc_mat)

fc_mat_p = reshape2::dcast(gxdiff[-which(gxdiff$term == "unactivated"),],
                           feature_id ~ comparison, value.var = "signedP")
fc_mat_p = column_to_rownames(fc_mat_p, var = "feature_id")
fc_mat_p = as.matrix(fc_mat_p)

fc_mat_ep = reshape2::dcast(gxdiff[-which(gxdiff$term == "unactivated"),],
                           feature_id ~ comparison, value.var = "estimate_logp")
fc_mat_ep = column_to_rownames(fc_mat_ep, var = "feature_id")
fc_mat_ep = as.matrix(fc_mat_ep)

insilico = matrix(data = NA, nrow = nrow(fc_mat), dimnames = list(rownames(fc_mat),"X2"))
insilico["117194",] <- 1


# Propagation in ARCHS
# =================================
# adj_matrix <- readRDS(get_task_outputs('wf-3917825f8b','0'))
adj_matrix = consume_workflow_task_output("wf-3917825f8b", task_id = 0)
tags <- list(list(name="service", value="patern"),
             list(name="network", value="ARCHS"),
             list(name="network_wfid", value="wf-3917825f8b"),
             list(name="alphas", value=paste0(alphas, collapse = '_')),
             list(name="signatures",value="X2"))



# Propagation in STRINGdb
# =================================
net_file <- get_task_inputs(res = 'wf-690e1e4e36', task_id = 0, files_names_grepl_pattern = '^g$')
adj_matrix <- readRDS(net_file)
adj_matrix = as.matrix(adj_matrix)
tags <- list(list(name="service", value="patern"),
             list(name="network", value="STRINGdb"),
             list(name="network_wfid", value="wf-690e1e4e36"),
             list(name="alphas", value=paste0(alphas, collapse = '_')),
             list(name="signatures",value="X2"))


# Run Function
# =================================
alphas <- c(0.25, 0.5, 0.75)
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
# fc_mat = fc_mat,
# fc_mat = insilico,
# fc_mat = fc_mat_p,
fc_mat = fc_mat_ep,
edge_cutoff = NA,
edge_cutoff_abs = FALSE,
stop_val = 1e-03,
direction_split = T,
image = 'eu.gcr.io/cytoreason/ci-cytoreason.patern-package:ori_latest',
cpu_request = '2000m',
memory_request = '25Gi',
tags = tags)

# wf-57d094c630 - ARCHS - comparisons based on estimate
# wf-48a7d142f5 - ARCHS - in silico
# wf-0e724fe2fe - ARCHS - comparisons based on signed -log10(pvalue)
# wf-c09ac29d65 - ARCHS - comparisons based on logP * estimate
# wf-b06f33a12c - STRING - in silico
# wf-77d7f5bbee - STRING - comparisons based on estimate
# wf-420015561b - STRING - comparisons based on signed -log10(pvalue)
# wf-01758afcc6 - STRING - comparisons based on logP * estimate

# Useful functions
# =================================
alphas <- c(0.25, 0.5, 0.75)

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
          if (grepl("activated", col)) {
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
      str_replace("_inhibited_vs_uninhibited","_inhibition") %>%
      str_replace("activated_vs_unactivated","x2_activation") %>%
      str_replace("inhibited_vs_uninhibited","x2_inhibition") %>%
      paste0(names(signatureNames),"_",suffix)
    if(insilico) {
      signatureNames = str_replace(signatureNames, "X2_","x2_insilico_")
    }
    names(x) = signatureNames
    return(x)
  })
}


# 1. ARCHS - gxdiff
# --------------------------------
# 1.1. based on estimates
# --------------------------
X2_Smoothed = extract_signatures("wf-57d094c630", gxdiff = T, extract_top_50 = T)
pushToCC(X2_Smoothed, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_withAlphas")))
# wf-f37d18b864

x2_archs = nameChange(X2_Smoothed, "archs")
pushToCC(x2_archs, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_alpha0.75"),
                                      list(name="network",value="archs"),
                                      list(name="origin",value="gx_diff")))
# wf-f5930577e9


# 1.2. based on pvalues
# ------------------------
X2_Smoothed = extract_signatures("wf-0e724fe2fe", gxdiff = T, extract_top_50 = T)
pushToCC(X2_Smoothed, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_withAlphas")))
# wf-34e747bb2d

x2_archs = nameChange(X2_Smoothed, "p_archs")
pushToCC(x2_archs, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_alpha0.75"),
                                     list(name="network",value="archs"),
                                     list(name="origin",value="gx_diff")))
# wf-d8c2843d95


# 1.3. based on logP * estimate
# ---------------------------------
X2_Smoothed = extract_signatures("wf-c09ac29d65", gxdiff = T, extract_top_50 = T)
pushToCC(X2_Smoothed, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_withAlphas")))
# wf-b0443d3a6c

x2_archs = nameChange(X2_Smoothed, "ep_archs")
pushToCC(x2_archs, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_alpha0.75"),
                                     list(name="network",value="archs"),
                                     list(name="origin",value="gx_diff")))
# wf-36b72dc9e9



# 2. STRING - gxdiff
# ------------------------------------
# 2.1. based on estimates
# ---------------------------
X2_Smoothed = extract_signatures("wf-77d7f5bbee", gxdiff = T, extract_top_50 = T)
pushToCC(X2_Smoothed, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_withAlphas"),
                                        list(name="network",value="string"),
                                        list(name="origin",value="gx_diff")))
# wf-32818a296c

x2_string = nameChange(X2_Smoothed, "string")
pushToCC(x2_string, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_alpha0.75"),
                                      list(name="network",value="string"),
                                      list(name="origin",value="gx_diff")))
# wf-ad41f5ce55


# 2.2. based on pvalues
# ------------------------
X2_Smoothed = extract_signatures("wf-420015561b", gxdiff = T, extract_top_50 = T)
pushToCC(X2_Smoothed, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_withAlphas")))
# wf-34e747bb2d

x2_string = nameChange(X2_Smoothed, "p_string")
pushToCC(x2_string, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_alpha0.75"),
                                     list(name="network",value="string"),
                                     list(name="origin",value="gx_diff")))
# wf-9e4ec09432


# 2.3. based on logP * estimate
# ---------------------------------
X2_Smoothed = extract_signatures("wf-01758afcc6", gxdiff = T, extract_top_50 = T)
pushToCC(X2_Smoothed, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_withAlphas")))
# wf-c797aa790a

x2_string = nameChange(X2_Smoothed, "ep_string")
pushToCC(x2_string, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_alpha0.75"),
                                      list(name="network",value="string"),
                                      list(name="origin",value="gx_diff")))
# wf-bbda22b623



# 3. ARCHS - in silico
# ---------------------------------------------
X2_Smoothed = extract_signatures("wf-48a7d142f5", gxdiff = F, extract_top_50 = T)
pushToCC(X2_Smoothed, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_withAlphas")))
# wf-54c7f42333

x2_archs_insilico = nameChange(X2_Smoothed, "archs", insilico = T)
pushToCC(x2_archs_insilico, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_alpha0.75"),
                                               list(name="network",value="archs"),
                                               list(name="origin",value="in_silico")))
# wf-6028a1f83b



# 4. STRING - in silico
# ---------------------------------------------
X2_Smoothed = extract_signatures("wf-b06f33a12c", gxdiff = F, extract_top_50 = T)
pushToCC(X2_Smoothed, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_withAlphas")))
# wf-c3a57f4af6

x2_string_insilico = nameChange(X2_Smoothed, "string", insilico = T)
pushToCC(x2_string_insilico, tagsToPass = list(list(name = "object",value="X2_smoothed_signatures_alpha0.75"),
                                               list(name="network",value="string"),
                                               list(name="origin",value="in_silico")))
# wf-b6999b00da


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
archs = append_signatures("wf-f5930577e9","estimate")
archs_p = append_signatures("wf-d8c2843d95","pvalue")
archs_ep = append_signatures("wf-36b72dc9e9","estimate_logp")
string = append_signatures("wf-ad41f5ce55","estimate")
string_p = append_signatures("wf-9e4ec09432","pvalue")
string_ep = append_signatures("wf-bbda22b623","estimate_logp")
insilico_archs = append_signatures("wf-6028a1f83b","insilico")
insilico_string = append_signatures("wf-b6999b00da","insilico")

all = rbind(archs, archs_p, archs_ep, string, string_p, string_ep, insilico_archs, insilico_string)

uploadToBQ(all, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures", disposition = "WRITE_APPEND")



# 6. Add in-silico from string
# -------------------------------------
x2_fromString <- cytoreason.assets::read_asset("ccw://wf-9329a6bc38@0:signatures_18803")
x2_fromString = x2_fromString[["MRGPRX2"]]
x2_fromString = data.frame(feature_id = x2_fromString, estimate = NA, agonist = NA, signature = "MRGPRX2", comparison = NA, log10_fdr = NA,
                           log10_pvalue = NA, in50 = NA, wfid = "wf-9329a6bc38", rankingMetric = "insilico_Shiran")
uploadToBQ(x2_fromString, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures", disposition = "WRITE_APPEND")
