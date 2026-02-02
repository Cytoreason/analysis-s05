library(cytoreason.ccm.pipeline)
library(cytoreason.cc.client)
library(tidyverse)
library(reshape2)
devtools::load_all("~/analysis-s05/R/utils.R")

### X2 signatures signal - narrowing the feature space with the AD network
### -------------------------------------------------------------------------------------
bulk_ad_network = readRDS(get_workflow_outputs("wf-06da2e64ae", files_names_grepl_pattern = "igraph_NW_sub_0.6.rds")) # ad_nw_bulk_pruned_0.6
bulk_ad_network = igraph::as_data_frame(bulk_ad_network)

allgenes_network = c(bulk_ad_network$from, bulk_ad_network$to)
allgenes_network = as.data.frame(table(allgenes_network))
allgenes_keep = as.character(allgenes_network$allgenes_network[which(allgenes_network$Freq > 4)])

X2_4hr = readRDS(get_workflow_outputs("wf-09b09047b2"))
X2_24hr = readRDS(get_workflow_outputs("wf-84a90d45bb"))

X2_4hr = X2_4hr[which(rownames(X2_4hr) %in% allgenes_keep),]
# wf-fd7f5efe8d - Freq >1, 9801 genes
# wf-eeeba8d2cd - Freq >2, 7542 genes
# wf-7f39d71564 - Freq >4, 4948 genes
X2_24hr = X2_24hr[which(rownames(X2_24hr) %in% allgenes_keep),]
# wf-f54269d9d2 - Freq >1
# wf-0befa8710e - Freq >2
# wf-5b067e0f12 - Freq >4

config = readRDS(get_workflow_outputs("wf-a399720b79"))
config$asset_id[which(config$experiment_id == "X2_4hr")] <- "ccw://wf-7f39d71564:0:output.rds"
config$asset_id[which(config$experiment_id == "X2_24hr")] <- "ccw://wf-5b067e0f12:0:output.rds"

ccm_stage_prepare_dataset_collection(config) # you need to have a folder called 'output' in the current wd

ccm_api_generate_ccm(config = config,
                     adjustment_models = list(),
                     stages = .skip("meta", "dataset2","dataset"),
                     model = .skip('gene_cell_correlations','cell_cell_correlations','feature_cell_correlations','pheno_feature_correlations',
                                   'cell_specific_differences','gene_set_activity_differences','survival_analysis_tme'),
                     image = "master@0.72.0")
# generate_ccm -- Sun Nov  2 12:29:04 2025: wf-ed79c8d52c [] - reduction to ~9.8k genes, Freq>1
# generate_ccm -- Sun Nov  2 12:31:01 2025: wf-97a7fccee9 [] - reduction to ~7.5k genes, Freq>2
# generate_ccm -- Sun Nov  2 12:33:49 2025: wf-f9d7e2846b [] - reduction to ~5k genes, Freq>4


wfids = c("Freq>1" = "wf-ed79c8d52c","Freq>2" = "wf-97a7fccee9", "Freq>4" = "wf-f9d7e2846b")

gx = list()
lapply(wfids, function(wfid){
  ccm = as_ccm_fit(wfid)
  allOptions = lapply(names(ccm$datasets), function(dataset){
    data.frame(experiment = dataset,
               comparison = names(ccm$datasets[[dataset]]$model))
  }) %>% do.call(rbind,.)
  allOptions = allOptions[-which(str_detect(allOptions$comparison,"_cov")),]
  
  apply(allOptions, 1, function(x){
    gx_tmp = statistic_table(analysisResultElement(ccm$datasets[[x[1]]]$model[[x[2]]],"gx_diff"))
    gx_tmp = gx_tmp[which(gx_tmp$term %in% c("activated_vs_unactivated","inhibited_vs_uninhibited")),]
    agonist = strsplit(x[2],"_")[[1]][1]
    if(agonist %in% c("inhibited","activated")) { agonist = "All" }
    if(agonist %in% c("SP","All")) { agonist = paste0(agonist,"_", strsplit(x[1],"_")[[1]][2]) }
    
    gx_tmp = data.frame(wfid = names(wfids)[which(wfids == wfid)], dataset_id = x[1], comparison = x[2], agonist = agonist, gx_tmp)
    gx <<- append(gx, list(gx_tmp))
    NULL
  })
})

gx_binded = do.call(rbind, gx)
gx_binded$term = droplevels(gx_binded$term)
sigGenes.counts = gx_binded %>%
  dplyr::filter(fdr <= 0.05) %>%
  group_by(wfid, agonist, term, comparison) %>%
  summarise(nSigGenes = n())

# Graph to see if there's a signal
ggplot(sigGenes.counts, aes(x = agonist, y = nSigGenes, fill = wfid)) +
  geom_col(position = "dodge") +
  facet_grid(cols = vars(agonist), rows = vars(term), scales = "free_x")
