library(reshape2)
library(ComplexHeatmap)
library(cytoreason.cc.client)
library(tidyverse)

enrichment_results = apply(get_workflow_outputs("wf-5cb62acda1", files_names_grepl_pattern = ".csv"), 1, read.csv, sep = ",")
enrichment_results = lapply(enrichment_results, function(x) {colnames(x)[which(colnames(x) == "geneset")] <- "signature" ; return(x)})
enrichedindisease =enrichment_results$statistics_X2_enrichedInDisease.csv

enrichment_results = apply(get_workflow_outputs("wf-0fdb40ddc1", files_names_grepl_pattern = ".csv"), 1, read.csv, sep = ",")
enrichment_results = lapply(enrichment_results, function(x) {colnames(x)[which(colnames(x) == "geneset")] <- "signature" ; return(x)})
enrichedintreatment =enrichment_results$treatments_X2_enrichments.csv

rm(enrichment_results)

# sig.act = c("x2_activation_early_50_p","x2_activation_early_50_p_archs_refined")
sig.inh = c("x2_inhibition_early_50_p","x2_inhibition_early_50_p_archs_refined","x2_inhibition_late_50_p","x2_inhibition_late_50_p_archs_refined")

enrichedindisease = enrichedindisease[which(enrichedindisease$signature %in% c(sig.inh)),]
enrichedintreatment = enrichedintreatment[which(enrichedintreatment$signature %in% c(sig.inh)),]

enrichedindisease = enrichedindisease[which(enrichedindisease$submodel == "bulk"),c(1,4,7,9,10,12,15,16)]
enrichedintreatment = enrichedintreatment[,c(2,4,5,6,8:10,13,14,16,19,20)]



dis = dcast(enrichedindisease, disease ~ signature, value.var = "NES") %>% column_to_rownames(var = "disease") %>% as.matrix
dis.st = dcast(enrichedindisease, disease ~ signature, value.var = "star") %>% column_to_rownames(var = "disease") %>% as.matrix

Heatmap(dis, cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(dis.st[i, j], x, y, gp = gpar(fontsize = 12))}, 
  column_labels = c("early inhibition","early inhibition\nrefined","late inhibition","late inhibition\nrefined"),
  name = "NES", column_title = "Signature enrichment across diseases")


enrichedintreatment = enrichedintreatment %>%
  rowwise() %>%
  mutate(treatment = paste(Disease,Treatment,Dosage,Timepoint,dataset,sep = "_"))

tre = enrichedintreatment %>%
  dplyr::filter(Disease == "AD") %>%
  mutate(treatment = str_remove(treatment,"_NA"))

tre.mat = dcast(tre, treatment ~ signature, value.var = "NES") %>% column_to_rownames(var = "treatment") %>% as.matrix
tre.st = dcast(tre, treatment ~ signature, value.var = "star") %>% column_to_rownames(var = "treatment") %>% as.matrix

Heatmap(tre.mat, cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(tre.st[i, j], x, y, gp = gpar(fontsize = 8))}, 
  column_labels = c("early inhibition","early inhibition\nrefined","late inhibition","late inhibition\nrefined"),
  name = "NES", column_title = "Signature enrichment across treatments", row_names_gp = gpar(fontsize = 8))




# gene overlap
# ---------------------
overlap = get_workflow_outputs("wf-76f2324cdf", files_names_grepl_pattern = ".csv", should_download_files = T)
overlap = read.csv(overlap["geneOverlap_X2.csv",])
overlap = overlap[which(overlap$genelist_1 %in% sig.inh & overlap$genelist_2 %in% sig.inh) ,]

overlap.mat = dcast(overlap, genelist_1 ~ genelist_2, value.var = "overlap") %>% column_to_rownames(var = "genelist_1") %>% as.matrix

Heatmap(overlap.mat, cell_fun = function(j, i, x, y, width, height, fill) {
  if(overlap.mat[i,j] > 0) grid.text(overlap.mat[i, j], x, y, gp = gpar(fontsize = 10))}, 
  column_labels = c("early inhibition","early inhibition\nrefined","late inhibition","late inhibition\nrefined"),
  row_labels = c("early inhibition","early inhibition\nrefined","late inhibition","late inhibition\nrefined"),
  name = "overlapping\ngenes", column_title = "Signature overlap")


