### Using the results from the custom CCM run, we will extract the result tables we are interested 
### in, and customize the tables to have a uniform built. 

# 1. Libraries and data loading
# ------------------------------------------
devtools::load_all("~/analysis-p03-poc1/R/mayan.R")
library(cytoreason.ccm.pipeline) # make sure to load version >= 1.0.1
library(tidyverse)


ccm = as_ccm_fit("wf-eafe4014ef")
targetMapping = read_asset("wf-01ed02b9dd") # taken from "2.2. Generating Signatures - Bispecifics.R"
geneCollections = unique(targetMapping$Target.Collection)
geneMapping = toTable(org.Hs.eg.db::org.Hs.egSYMBOL)

# 2. Extraction
# ----------------------------
# We use the L_vs_NL effect to extract results because we need only the information in within lesional samples
# and the L_vs_NL effect has the highest number of datasets
# this stage needs a lot of RAM because of the large number of signatures, make sure your capsule has ~100GB
run_function_dist(FUN = function(ccm_wfid){
  library(cytoreason.ccm.pipeline)
  
  ccm = as_ccm_fit(ccm_wfid)
  geneMapping = toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
  
  MS.CS = statistic_table(ccm$meta$L_vs_NL$pheno_feature_correlations) # MS/CS - L
    MS.CS = MS.CS$gene_set_phenotypic_variable_correlations
    
  TargetFeatures = statistic_table(ccm$meta$L_vs_NL$feature_geneset_correlations) # genes/cells/pathways/metaPCs - L
  DZEnrichment = statistic_table(ccm$meta$L_vs_NL$gx_diff$gx_gsa) # enrichment in disease - L
  
  leadingEdge = DZEnrichment
  leadingEdge$ID = paste0(leadingEdge$collection, "__", leadingEdge$pathway)
  leadingEdge$leadingEdge = sapply(leadingEdge$hit, function(x) {
    x = strsplit(x, "\\;")[[1]]
    x = geneMapping$symbol[match(x, geneMapping$gene_id)]
    x = paste(x, collapse = ";")
  })
  leadingEdge = leadingEdge[,c("ID", "leadingEdge")] # wf-18e19958fc
  
  Results = list(DZEnrichment = DZEnrichment, TargetFeatures = TargetFeatures, MS.CS = MS.CS)
  return(Results)
}, 
ccm_wfid = "wf-eafe4014ef", 
memory_request = "100Gi",
image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:master_latest")
# wf-67799934cd

Results = readRDS(get_workflow_outputs("wf-70a31a8336"))

# 3. Editing to keep only what we're interested in
# ----------------------------------------------------
# This function filters the tables to the terms we want to keep:
# L vs NL for enrichment in disease but in calculations within lesional samples for the rest of the criteria
# And unifies the column names - what's related to the signatures is under 'Target'
# and what's related to the criteria is under 'Criteria'.

processResults = function(data, tableName) {
  if(tableName == "DZEnrichment") {
    targetID = "pathway"
    criteriaID = "term"
    colnames(data)[which(colnames(data) == "nes")] <- "value"
    data$DataType = paste0("Enrichment in ",data$term)
    data$pathway = paste0(data$collection, "__",data$pathway)
    terms = c("L_vs_NL","L_vs_HC","NL_vs_HC","AD","DZ_vs_HC")
  } else if(tableName == "TargetFeatures") {
    targetID = "feature_id_2"
    criteriaID = "feature_id_1"
    data$feature_type1 = str_replace(data$feature_type1, "gene_set", "pathway")
    data$DataType = paste0("Target_", str_to_title(data$feature_type1))
    data$collection = str_split(data$feature_id_2,"__",simplify = TRUE)[,1]
    terms = "L"
  } else {
    targetID = "feature_id_1"
    criteriaID = "feature_id_2"
    data$DataType = paste0("Target_", sapply(data$feature_id_2, function(x) ifelse(strsplit(x,"\\.")[[1]][1] == "MS", "MS", "CS")))
    data$collection = str_split(data$feature_id_1,"__",simplify = TRUE)[,1]
    terms = "L"
  }
  
  newtable = data %>%
    dplyr::filter(submodel == "bulk", term %in% terms) %>%
    dplyr::rename(Target.Identifier = targetID) %>%
    mutate(Target.ID = Target.Identifier) %>%
    mutate(Target.Identifier = str_split(Target.Identifier,"__",simplify = TRUE)[,2]) %>%
    dplyr::rename(Target.Collection = collection) %>%
    dplyr::filter(Target.Collection %in% geneCollections) %>%
    mutate(Criteria.Identifier = get(criteriaID)) %>%
    mutate(Criteria.Collection = str_split(Criteria.Identifier,"__",simplify = TRUE)[,1]) %>%
    mutate(log10_fdr = -log10(fdr)) %>%
    mutate(metricType = case_when(tableName == "DZEnrichment" ~ "NES", .default = "bi-weight mid-correlation")) %>%
    dplyr::rename(metricValue = value, Type = submodel) %>%
    dplyr::select(DataType, Target.ID, Target.Identifier, Target.Collection, Type, Criteria.Identifier, Criteria.Collection, metricValue, metricType, log10_fdr)
  
  return(newtable)
}


Results = lapply(names(Results), function(dataName) processResults(data = Results[[dataName]], tableName = dataName))
Results = list(DZEnrichment = Results[[1]],
               Target_Cell = Results[[2]][which(Results[[2]]$DataType == "Target_Cell"),],
               Target_Gene = Results[[2]][which(Results[[2]]$DataType == "Target_Gene"),],
               Target_Pathway = Results[[2]][which(Results[[2]]$DataType == "Target_Pathway"),],
               Target_MS = Results[[3]][which(Results[[3]]$DataType == "Target_MS"),],
               Target_CS = Results[[3]][which(Results[[3]]$DataType == "Target_CS"),])

rm(DZEnrichment, MS.CS, TargetFeatures)
