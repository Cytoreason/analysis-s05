devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)

# ccm = as_ccm_fit("wf-08a6a0a503")
wantedTerms = read.delim("~/analysis-s05/data/wantedTerms.tsv")

# 1. Dupi pathway white space - gx_gsa
# ==============================================
# 1.1. Extract data
# ------------------------
run_function_dist(FUN = function(ccm_wfid, wantedTerms){
  library(cytoreason.ccm.pipeline)
  
  ccm = as_ccm_fit(ccm_wfid)
  
  Treatments = lapply(unique(wantedTerms$data_set), function(datasetID){
    res=lapply(unique(wantedTerms$model[which(wantedTerms$data_set == datasetID)]), function(model){
      cat("\n",datasetID,model,"\n")
      terms = unique(wantedTerms$Terms[which(wantedTerms$data_set == datasetID & wantedTerms$model == model)])
      res = statistic_table(analysisResultElement(ccm$datasets[[datasetID]]$model[[model]],"gx_gsa"),
                            submodel = c("bulk","adjusted__1__1"), term = terms)
    })
    names(res) = unique(wantedTerms$model[which(wantedTerms$data_set == datasetID)])
    return(res)
  })
  names(Treatments) = unique(wantedTerms$data_set)
  return(Treatments)
},
ccm_wfid = "wf-08a6a0a503",
wantedTerms = wantedTerms,
memory_request = "5Gi",
image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:master_latest")
# wf-ca1574909f

Treatment_path = readRDS(get_workflow_outputs("wf-ca1574909f"))
Treatment_path = lapply(names(Treatment_path), function(x) {
  y = bind_rows(Treatment_path[[x]]) %>%
    dplyr::filter(collection %in% c("h","kegg","reactome"))
  y$dataset = x
  return(y)
}) %>% bind_rows()

uploadToBQ(Treatment_path, bqdataset = "s05_atopic_dermatitis", tableName = "treatment_PrePostDupi")
pushToCC(Treatment_path, tagsToPass = list(list(name="object",value="treatment_PrePostDupi")))
# wf-eb2ac52a50

# 1.2. Define white space
# ---------------------------
pathways_in_dz = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "AD_gx_gsa")
pathways_in_dz = pathways_in_dz %>%
  group_by(term, submodel) %>%
  dplyr::filter(submodel %in% c("bulk","adjusted__1__1") & collection %in% c("h","kegg","reactome")) %>%
  dplyr::select(term, submodel, collection, pathway, NES, FDR, log10_fdr) %>%
  mutate(dir = ifelse(FDR > 0.05, "unchanged_in_dz",
                      ifelse(NES >0, "up_in_dz","down_in_dz"))) %>%
  dplyr::rename(term_disease = term) %>%
  ungroup()
  

pathway_ws = Treatment_path %>%
  group_by(dataset, term, submodel) %>%
  dplyr::filter(str_detect(term,"Dupilumab") & collection %in% c("h","kegg","reactome")) %>%
  dplyr::select(dataset, term, submodel, collection, pathway, NES, FDR, log10_fdr) %>%
  mutate(dir = ifelse(FDR > 0.05, "unchanged_in_dupi",
                      ifelse(NES >0, "up_in_dupi","down_in_dupi"))) %>%
  mutate(term = as.character(term)) %>%
  dplyr::rename(term_treatment = term) %>%
  ungroup()


pathways_ws = merge(pathways_in_dz, pathway_ws, by = c("submodel","collection","pathway"), suffix = c("_disease","_treatment"))
pushToCC(pathways_ws, tagsToPass = list(list(name="object",value="pathways_ws_all")))
# wf-152a08ff2d

pathways_ws_filter = pathways_ws %>%
  dplyr::filter(dir_disease != "unchanged_in_dz" & dir_treatment == "unchanged_in_dupi") %>%
  dplyr::filter(term_treatment %in% c("W16_vs_W0:DupilumabL","W12_vs_W0:DupilumabL","W4_vs_W0:DupilumabL"))

uploadToBQ(pathways_ws_filter, bqdataset = "s05_atopic_dermatitis", tableName = "whiteSpace_pathways")
pushToCC(pathways_ws_filter, tagsToPass = list(list(name="object",value="pathways_ws_filter")))
# wf-84f066601f



# 2. Dupi cell white space - ct_test
# --------------------------------------------------------------
run_function_dist(FUN = function(ccm_wfid, wantedTerms){
  library(cytoreason.ccm.pipeline)
  
  ccm = as_ccm_fit(ccm_wfid)
  
  Treatments = lapply(unique(wantedTerms$data_set), function(datasetID){
    res=lapply(unique(wantedTerms$model[which(wantedTerms$data_set == datasetID)]), function(model){
      cat("\n",datasetID,model,"\n")
      terms = unique(wantedTerms$Terms[which(wantedTerms$data_set == datasetID & wantedTerms$model == model)])
      res = statistic_table(analysisResultElement(ccm$datasets[[datasetID]]$model[[model]],"ct_test"), term = terms)
    })
    names(res) = unique(wantedTerms$model[which(wantedTerms$data_set == datasetID)])
    return(res)
  })
  names(Treatments) = unique(wantedTerms$data_set)
  return(Treatments)
},
ccm_wfid = "wf-08a6a0a503",
wantedTerms = wantedTerms[which(wantedTerms$drug == "Dupilumab"),],
memory_request = "5Gi",
image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:master_latest")
# wf-945dca171d

Treatment_cells = readRDS(get_workflow_outputs("wf-945dca171d"))
Treatment_cells = lapply(names(Treatment_cells), function(x) {
  y = bind_rows(Treatment_cells[[x]])
  y$dataset = x
  return(y)
}) %>% bind_rows()


pushToCC(Treatment_cells, tagsToPass = list(list(name="object",value="treatment_PrePostDupi_cells")))
# wf-bac4719c73


# 2.2. Define white space
# ---------------------------
cells_in_dz = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "AD_ct_test")
cells_in_dz = cells_in_dz %>%
  group_by(term) %>%
  dplyr::select(term, feature_id, estimate, fdr) %>%
  mutate(dir = ifelse(fdr > 0.05, "unchanged_in_dz",
                      ifelse(estimate >0, "up_in_dz","down_in_dz"))) %>%
  mutate(log10_fdr = -log10(fdr)) %>%
  dplyr::rename(term_disease = term, cell = feature_id) %>%
  ungroup()


cells_ws_all = Treatment_cells %>%
  group_by(dataset, term) %>%
  dplyr::filter(str_detect(term,"Dupilumab")) %>%
  dplyr::select(dataset, term, feature_id, estimate, FDR, log10_fdr) %>%
  mutate(dir = ifelse(FDR > 0.05, "unchanged_in_dupi",
                      ifelse(estimate >0, "up_in_dupi","down_in_dupi"))) %>%
  mutate(term = as.character(term)) %>%
  dplyr::rename(term_treatment = term, cell = feature_id) %>%
  ungroup()


cells_ws = merge(cells_in_dz, cells_ws_all, by = c("cell"), suffix = c("_disease","_treatment"))
pushToCC(cells_ws, tagsToPass = list(list(name="object",value="cells_ws_all")))
# wf-2e7a260c22

cells_ws_filter = cells_ws %>%
  dplyr::filter(dir_disease != "unchanged_in_dz" & dir_treatment == "unchanged_in_dupi") %>%
  dplyr::filter(term_treatment %in% c("W16_vs_W0:DupilumabL","W12_vs_W0:DupilumabL","W4_vs_W0:DupilumabL"))

uploadToBQ(cells_ws_filter, bqdataset = "s05_atopic_dermatitis", tableName = "whiteSpace_cells")
pushToCC(cells_ws_filter, tagsToPass = list(list(name="object",value="cells_ws_filtered")))
# wf-7b6f3f86f7
