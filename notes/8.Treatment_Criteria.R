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
ccm_wfid = "wf-8e948630d7",
wantedTerms = wantedTerms,
memory_request = "5Gi",
image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:master_latest")
# wf-ca1574909f

Treatment_path = readRDS(get_workflow_outputs("wf-ca1574909f"))
Treatment_path = lapply(names(Treatment_path), function(x) {
  y = bind_rows(Treatment_path[[x]]) %>%
    dplyr::filter(collection %in% c("h","kegg","reactome","btm"))
  y$dataset = x
  return(y)
}) %>% bind_rows()

uploadToBQ(Treatment_path, bqdataset = "s05_atopic_dermatitis", tableName = "treatment_PrePostDupi")
pushToCC(Treatment_path, tagsToPass = list(list(name="object",value="treatment_PrePostDupi")))
# wf-88d79c9c27

# 1.2. Define white space
# ---------------------------
pathways_in_dz = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "AD_gx_gsa")
pathways_in_dz = pathways_in_dz %>%
  group_by(term, submodel) %>%
  dplyr::filter(submodel %in% c("bulk","adjusted__1__1") & collection %in% c("h","kegg","reactome","btm")) %>%
  dplyr::select(term, submodel, collection, pathway, NES, FDR, log10_fdr) %>%
  mutate(dir = ifelse(FDR > 0.05, "unchanged",
                      ifelse(NES >0, "up","down"))) %>%
  ungroup()
  
dz_wide = pathways_in_dz %>%
  pivot_wider(id_cols = c(pathway,collection, submodel),
              names_from  = term,
              values_from = dir)

pathway_treatment = Treatment_path %>%
  group_by(dataset, term, submodel) %>%
  dplyr::filter(str_detect(term,"Dupilumab") & collection %in% c("h","kegg","reactome","btm")) %>%
  dplyr::select(dataset, term, submodel, collection, pathway, NES, FDR, log10_fdr) %>%
  mutate(term = as.character(term)) %>%
  ungroup()


pathway_space = merge(pathway_treatment, dz_wide, all = T)
pathway_space = pathway_space %>%
  mutate(in_white_space_L_vs_NL = case_when(NES > 0 & L_vs_NL == "up" ~ "yesUp",
                                            NES < 0 & FDR > 0.05 & L_vs_NL == "up" ~ "yesUp",
                                            NES < 0 & L_vs_NL == "down" ~ "yesDown",
                                            NES > 0 & FDR > 0.05 & L_vs_NL == "down" ~ "yesDown",
                                            .default = "no")) %>%
  mutate(in_white_space_L_vs_HC = case_when(NES > 0 & L_vs_HC == "up" ~ "yesUp",
                                            NES < 0 & FDR > 0.05 & L_vs_HC == "up" ~ "yesUp",
                                            NES < 0 & L_vs_HC == "down" ~ "yesDown",
                                            NES > 0 & FDR > 0.05 & L_vs_HC == "down" ~ "yesDown",
                                            .default = "no")) %>%
  mutate(in_white_space_NL_vs_HC = case_when(NES > 0 & NL_vs_HC == "up" ~ "yesUp",
                                             NES < 0 & FDR > 0.05 & NL_vs_HC == "up" ~ "yesUp",
                                             NES < 0 & NL_vs_HC == "down" ~ "yesDown",
                                             NES > 0 & FDR > 0.05 & NL_vs_HC == "down" ~ "yesDown",
                                             .default = "no")) %>%
  mutate(in_white_space_DZ_vs_HC = case_when(NES > 0 & DZ_vs_HC == "up" ~ "yesUp",
                                             NES < 0 & FDR > 0.05 & DZ_vs_HC == "up" ~ "yesUp",
                                             NES < 0 & DZ_vs_HC == "down" ~ "yesDown",
                                             NES > 0 & FDR > 0.05 & DZ_vs_HC == "down" ~ "yesDown",
                                             .default = "no")) 

pushToCC(pathway_space, tagsToPass = list(list(name="object",value="pathways_ws_all")))
# wf-8fc42d8d1b
uploadToBQ(pathways_ws_filter, bqdataset = "s05_atopic_dermatitis", tableName = "whiteSpace_pathways")



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
  mutate(dir = ifelse(fdr > 0.05, "unchanged",
                      ifelse(estimate >0, "up","down"))) %>%
  mutate(log10_fdr = -log10(fdr)) %>%
  dplyr::rename(cell = feature_id) %>%
  ungroup()

dz_wide = cells_in_dz %>%
  pivot_wider(id_cols = cell,
              names_from  = term,
              values_from = dir)

cells_treatment = Treatment_cells %>%
  group_by(dataset, term) %>%
  dplyr::filter(str_detect(term,"Dupilumab")) %>%
  dplyr::select(dataset, term, feature_id, estimate, FDR, log10_fdr) %>%
  mutate(term = as.character(term)) %>%
  dplyr::rename(cell = feature_id) %>%
  ungroup()


cells_space = merge(cells_treatment, dz_wide, by = "cell", all = T)
cells_space = pathway_space %>%
  mutate(in_white_space_L_vs_NL = case_when(NES > 0 & L_vs_NL == "up" ~ "yesUp",
                                            NES < 0 & FDR > 0.05 & L_vs_NL == "up" ~ "yesUp",
                                            NES < 0 & L_vs_NL == "down" ~ "yesDown",
                                            NES > 0 & FDR > 0.05 & L_vs_NL == "down" ~ "yesDown",
                                            .default = "no")) %>%
  mutate(in_white_space_L_vs_HC = case_when(NES > 0 & L_vs_HC == "up" ~ "yesUp",
                                            NES < 0 & FDR > 0.05 & L_vs_HC == "up" ~ "yesUp",
                                            NES < 0 & L_vs_HC == "down" ~ "yesDown",
                                            NES > 0 & FDR > 0.05 & L_vs_HC == "down" ~ "yesDown",
                                            .default = "no")) %>%
  mutate(in_white_space_NL_vs_HC = case_when(NES > 0 & NL_vs_HC == "up" ~ "yesUp",
                                             NES < 0 & FDR > 0.05 & NL_vs_HC == "up" ~ "yesUp",
                                             NES < 0 & NL_vs_HC == "down" ~ "yesDown",
                                             NES > 0 & FDR > 0.05 & NL_vs_HC == "down" ~ "yesDown",
                                             .default = "no")) %>%
  mutate(in_white_space_DZ_vs_HC = case_when(NES > 0 & DZ_vs_HC == "up" ~ "yesUp",
                                             NES < 0 & FDR > 0.05 & DZ_vs_HC == "up" ~ "yesUp",
                                             NES < 0 & DZ_vs_HC == "down" ~ "yesDown",
                                             NES > 0 & FDR > 0.05 & DZ_vs_HC == "down" ~ "yesDown",
                                             .default = "no")) 

uploadToBQ(cells_space, bqdataset = "s05_atopic_dermatitis", tableName = "whiteSpace_cells")
pushToCC(cells_space, tagsToPass = list(list(name="object",value="cells_space")))
# wf-82b0b996ed



# 1. Dupi pathway white space - gx_diff
# ==============================================
# 1.1. Extract data
# ------------------------
run_function_dist(FUN = function(ccm_wfid, wantedTerms){
  library(cytoreason.ccm.pipeline)
  
  ccm = as_ccm_fit(ccm_wfid)
  
  Treatments = lapply(unique(wantedTerms$data_set), function(datasetID){
    res=lapply(unique(wantedTerms$model[which(wantedTerms$data_set == datasetID)]), function(model){
      cat("\n",datasetID,model,"\n")
      res = statistic_table(analysisResultElement(ccm$datasets[[datasetID]]$model[[model]],"gx_diff"),
                            submodel = c("bulk","adjusted__1__1"))
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
# wf-0ca32e1ef5

Treatment_genes = readRDS(get_workflow_outputs("wf-0ca32e1ef5"))
Treatment_genes = lapply(names(Treatment_genes), function(x) {
  lapply(Treatment_genes[[x]], function(y){
    y = y[which(y$term %in% wantedTerms$Terms),]
    y$dataset = x
    return(y)
  }) %>% bind_rows()
}) %>% bind_rows()

ccm = as_ccm_fit("wf-08a6a0a503")
dz = rbind(statistic_table(ccm$meta$L_vs_NL$gx_diff$gx_diff),
           statistic_table(ccm$meta$AD$gx_diff$gx_diff),
           statistic_table(ccm$meta$L_vs_HC$gx_diff$gx_diff),
           statistic_table(ccm$meta$NL_vs_HC$gx_diff$gx_diff)) %>%
  dplyr::filter(term %in% c("L_vs_NL","L_vs_HC","NL_vs_HC","DZ_vs_HC")) %>%
  dplyr::filter(submodel %in% c("bulk","adjusted__1__1"))
  
dz = dz %>%
  mutate(dir = case_when(fdr <= 0.05 & effect_size > 0 ~ "up",
                         fdr <= 0.05 & effect_size < 0 ~ "down",
                         .default = "unchanged"))
dz_wide = dz %>%
  pivot_wider(id_cols = c(feature_id, submodel),
              names_from  = term,
              values_from = dir)

gene_space = merge(Treatment_genes, dz_wide, all = T)
gene_space = gene_space %>%
  mutate(in_white_space_L_vs_NL = case_when(estimate > 0 & L_vs_NL == "up" ~ "yesUp",
                                            estimate < 0 & fdr > 0.05 & L_vs_NL == "up" ~ "yesUp",
                                            estimate < 0 & L_vs_NL == "down" ~ "yesDown",
                                            estimate > 0 & fdr > 0.05 & L_vs_NL == "down" ~ "yesDown",
                                            .default = "no")) %>%
  mutate(in_white_space_L_vs_HC = case_when(estimate > 0 & L_vs_HC == "up" ~ "yesUp",
                                            estimate < 0 & fdr > 0.05 & L_vs_HC == "up" ~ "yesUp",
                                            estimate < 0 & L_vs_HC == "down" ~ "yesDown",
                                            estimate > 0 & fdr > 0.05 & L_vs_HC == "down" ~ "yesDown",
                                            .default = "no")) %>%
  mutate(in_white_space_NL_vs_HC = case_when(estimate > 0 & NL_vs_HC == "up" ~ "yesUp",
                                             estimate < 0 & fdr > 0.05 & NL_vs_HC == "up" ~ "yesUp",
                                             estimate < 0 & NL_vs_HC == "down" ~ "yesDown",
                                             estimate > 0 & fdr > 0.05 & NL_vs_HC == "down" ~ "yesDown",
                                            .default = "no")) %>%
  mutate(in_white_space_DZ_vs_HC = case_when(estimate > 0 & DZ_vs_HC == "up" ~ "yesUp",
                                             estimate < 0 & fdr > 0.05 & DZ_vs_HC == "up" ~ "yesUp",
                                             estimate < 0 & DZ_vs_HC == "down" ~ "yesDown",
                                             estimate > 0 & fdr > 0.05 & DZ_vs_HC == "down" ~ "yesDown",
                                            .default = "no")) %>%
  dplyr::select(feature_id, submodel, dataset, term, estimate,fdr,log10_fdr,L_vs_NL,in_white_space_L_vs_NL,
                DZ_vs_HC, in_white_space_DZ_vs_HC, L_vs_HC, in_white_space_L_vs_HC, NL_vs_HC, in_white_space_NL_vs_HC)
  
pushToCC(gene_space, tagsToPass = list(list(name="object",value="gene_white_space")))
# wf-3d117171a9

# bq_table_create(bq_table("cytoreason", "s05_atopic_dermatitis", table = "whiteSpace_genes"), fields = as_bq_fields(gene_space)) # made upload problems without it
uploadToBQ(gene_space, bqdataset = "s05_atopic_dermatitis", tableName = "whiteSpace_genes")
