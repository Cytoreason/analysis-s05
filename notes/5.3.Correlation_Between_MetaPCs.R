devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(ComplexHeatmap)
library(tidyverse)


Results = readRDS(get_workflow_outputs("wf-798227871e"))
Results = Results[c("Target_Cell_PCA", "Target_Pathway_PCA")] %>% bind_rows()

Results.corr = reshape2::dcast(Results, Target.Identifier + Type ~ Criteria.Identifier, value.var = "metricValue") %>%
  split(., .[,"Type"])
Results.corr = lapply(Results.corr, function(submodel){
  submodel = submodel %>%
    dplyr::select(-Type) %>%
    data.frame(., row.names = NULL) %>%
    column_to_rownames(var = "Target.Identifier") %>%
    as.matrix
})
Results.fdr = reshape2::dcast(Results, Target.Identifier + Type ~ Criteria.Identifier, value.var = "log10_fdr") %>%
  split(., .[,"Type"])
Results.fdr = lapply(Results.fdr, function(submodel){
  submodel = submodel %>%
    dplyr::select(-Type) %>%
    data.frame(., row.names = NULL) %>%
    column_to_rownames(var = "Target.Identifier") %>%
    as.matrix
})

corr_corr = lapply(Results.corr, cor, method = "spearman")

png("~/analysis-s05/figures/Results/metaPCs_correlation_bulk.png", res = "100", bg = "white", width = 600, height = 500)
Heatmap(corr_corr$bulk, cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(round(corr_corr$bulk[i, j],1), x, y, gp = gpar(fontsize = 10))},
  name = "Spearman\nCorrelation", row_labels = str_replace_all(rownames(corr_corr$bulk),"_"," "), 
  column_labels = str_replace_all(rownames(corr_corr$bulk),"_"," "), column_title = "Bulk")
dev.off()

png("~/analysis-s05/figures/Results/metaPCs_correlation_adj.png", res = "100", bg = "white", width = 600, height = 500)
Heatmap(corr_corr$adjusted__1__1, cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(round(corr_corr$adjusted__1__1[i, j],1), x, y, gp = gpar(fontsize = 10))},
  name = "Spearman\nCorrelation", row_labels = str_replace_all(rownames(corr_corr$adjusted__1__1),"_"," "), 
  column_labels = str_replace_all(rownames(corr_corr$adjusted__1__1),"_"," "), column_title = "Adjusted")
dev.off()

