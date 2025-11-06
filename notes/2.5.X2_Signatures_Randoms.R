devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.cc.client)
library(tidyverse)
library(ComplexHeatmap)

gxdiff = readRDS(get_workflow_outputs("wf-edcade5c80"))
nPerm = 100

# 1. Negative Control - pure random
# =========================================
sample_gene_sets <- function(genes, sigSize, nPerm, prefix = "random") {
  sets <- lapply(1:nPerm, function(i) sample(genes, size = sigSize, replace = FALSE))
  names(sets) <- paste0(prefix, 1:nPerm)
  return(sets)
}

# Usage
genes <- unique(gxdiff$feature_id)

nc = list(nc_50 = sample_gene_sets(genes, sigSize = 50, nPerm = nPerm),
          nc_100 = sample_gene_sets(genes, sigSize = 100, nPerm = nPerm))
                           
pushToCC(nc, tagsToPass = list(list(name="object",value="randomNC"))) # wf-

overlap_50 = sapply(nc$nc_50, function(x) sapply(nc$nc_50, function(y) length(intersect(x,y))))
overlap_100 = sapply(nc$nc_100, function(x) sapply(nc$nc_100, function(y) length(intersect(x,y))))

png("~/analysis-s05/figures/X2_Signature/overlap_nc_50.png", res = "100", bg = "transparent", width = 1200, height = 900)
Heatmap(overlap_50, name = "overlap", row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize=5), cell_fun = function(j, i, x, y, width, height, fill) {
          if(overlap_50[i, j] >= 25) grid.text(overlap_50[i, j], x, y, gp = gpar(fontsize = 6, col = "white"))})
dev.off()

png("~/analysis-s05/figures/X2_Signature/overlap_nc_100.png", res = "100", bg = "transparent", width = 1200, height = 900)
Heatmap(overlap_100, name = "overlap", row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize=5), cell_fun = function(j, i, x, y, width, height, fill) {
          if(overlap_100[i, j] >= 25) grid.text(overlap_100[i, j], x, y, gp = gpar(fontsize = 6, col = "white"))})
dev.off()



# 2. Negative Control - random that passes smoothing
# ========================================================
# Strategy: from L/NL meta analysis, shuffle FCs and send to smoothing. use resulting signatures

ad.ccm = as_ccm_fit("wf-08a6a0a503")
genePool = build_service_result_tables(ad.ccm$meta$L_vs_HC$gx_diff$gx_diff)
genePool = genePool$bulk[which(genePool$bulk$term == "L_vs_HC"),]

nPerm = 100
ncFCs = lapply(1:nPerm, function(i){
  fcs = genePool[, c("feature_id","estimate")]
  fcs$feature_id = sample(fcs$feature_id, replace = F)
  return(data.frame(list = paste0("smoothedRandom",i),fcs,row.names = NULL))
}) %>% do.call(rbind,.)
ncFCs = dcast(ncFCs, feature_id ~ list, value.var = "estimate", fill = 0, drop = F)
rownames(ncFCs) = ncFCs$feature_id
ncFCs = ncFCs[,-1]

run_function_dist(function(alpha,
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
cor_mat = consume_workflow_task_output("wf-3917825f8b", task_id = 0),
alpha = 0.75,
fc_mat = ncFCs,
edge_cutoff = NA,
edge_cutoff_abs = FALSE,
stop_val = 1e-03,
direction_split = T,
image = 'eu.gcr.io/cytoreason/ci-cytoreason.patern-package:ori_latest',
cpu_request = '2000m',
memory_request = '25Gi',
tags = list(list(name="service", value="patern"),
            list(name="network", value="ARCHS"),
            list(name="network_wfid", value="wf-3917825f8b"),
            list(name="alpha", value="0.75"),
            list(name="signatures",value="randomControls")))

# wf-ee3aef5da2

ncPATERN = readRDS(get_workflow_outputs("wf-ee3aef5da2"))
ncPATERN = ncPATERN[["pos_smooth_mat"]]
ncPATERN = lapply(seq_len(ncol(ncPATERN)), function(i) {
    col_vals <- ncPATERN[, i]
    r = list(
      top50 = sort(col_vals, decreasing = TRUE)[1:50],
      bottom50 = sort(col_vals)[1:50],
      top100 = sort(col_vals, decreasing = TRUE)[1:100],
      bottom100 = sort(col_vals)[1:100]
    )
    names(r) = paste0("smoothedRandom", i,"_",names(r))
    return(r)
  })
ncPATERN = unlist(ncPATERN, recursive = F)                                                            

# Combine negative controls
# =======================================
negativeControls = list(random = nc,
                        smoothedRandom = ncPATERN)
pushToCC(negativeControls, tagsToPass = list(list(name="object",value="randomNC"))) # wf-905e97fa64
