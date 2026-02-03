devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(cytoreason.integration)
library(tidyverse)

# Define treatment white space
# =========================================
treatments = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "treatment_PrePostDupi")
treatments = treatments %>%
  dplyr::filter(str_detect(term, ":FezakinumabL|:DupilumabL|:SecukinumabL")) %>%
  dplyr::filter(submodel %in% c("bulk","adjusted__1__1") & collection %in% c("h","kegg","reactome","btm","th2","neuroinflammation")) %>%
  mutate(submodel = ifelse(submodel == "bulk", "bulk", "adjusted"))
  

pathways_in_dz = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "AD_gx_gsa")
pathways_in_dz = pathways_in_dz %>%
  group_by(term, submodel) %>%
  dplyr::filter(submodel %in% c("bulk","adjusted") & collection %in% c("h","kegg","reactome","btm","th2","neuroinflammation")) %>%
  dplyr::select(term, submodel, collection, pathway, NES, FDR, log10_fdr) %>%
  mutate(dir = ifelse(FDR > 0.05, "unchanged",
                      ifelse(NES >0, "up","down"))) %>%
  ungroup()

dz_wide = pathways_in_dz %>%
  pivot_wider(id_cols = c(pathway,collection, submodel),
              names_from  = term,
              values_from = dir)


treatment_ws = treatments %>%
  merge(., dz_wide, all = T, by = c("pathway","collection","submodel")) %>%
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
  
uploadToBQ(treatment_ws, bqdataset = "s05_atopic_dermatitis", tableName = "whiteSpace_treatments")

treatment_ws = treatment_ws %>%
  dplyr::filter(in_white_space_L_vs_HC != "no") %>%
  mutate(pathway_id = paste0(collection,":",pathway)) %>%
  dplyr::select(term, pathway_id, pathway, submodel, collection, in_white_space_L_vs_HC) %>%
  mutate(term = str_remove(term, "_vs_W0"))

pushToCC(treatment_ws, tagsToPass = list(list(name="object",value="treatment_ws")))
# wf-bf1b734156

# Calcaulate PCA
# =============================================================
ccm = as_ccm_fit("wf-08a6a0a503")

run_function_dist(FUN = function(ccm_wfid, treatment_ws){
  library(cytoreason.ccm.pipeline)
  library(cytoreason.integration)
  library(dplyr)
  install.packages("parallel")
  library(parallel)
  
  ccm = as_ccm_fit(ccm_wfid)
  config <- modelMetadata(ccm)
  allSubmodels = lapply(c("yesDown","yesUp"), function(ws){
    allWS = lapply(c("bulk","adjusted__1__1"), function(submodel){
      allTerms = lapply(unique(treatment_ws$term), function(term){
        cat("\nComputing:",ws,submodel,term,"\n")
        
        cl = makeCluster(length(ccm$datasets))
        clusterEvalQ(cl, { library(cytoreason.ccm.pipeline) })
        clusterExport(cl, varlist = c("ccm","ws","submodel","term","treatment_ws","nrow","ncol"), envir = environment())
        
        pathways <- parLapply(cl, ccm$datasets, function(d){
          x = analysisResultExpressionSet(d, "gene_set_activity")
          th2_enrichments = readRDS(get_workflow_outputs("wf-b2132fba22"))
          keep = ifelse(submodel == "bulk", "exprs", paste0("exprs_",submodel))
          if(submodel != "bulk") { submodel = "adjusted" }
          paths = treatment_ws[which(treatment_ws$term == term & 
                                               treatment_ws$submodel == submodel & 
                                               treatment_ws$in_white_space_L_vs_HC == ws),]

          x = x[which(rownames(x) %in% paths$pathway_id),]
          
          if("th2" %in% unique(paths$collection)){
            th2 = th2_enrichments[[keep]]
            th2 = lapply(th2, function(y) return(y[which(rownames(y) %in% paths$pathway),,drop=F]))
            new_exprs = rbind(assayData(x)[[keep]], th2[[x@experimentData@title]])
            
            fd_new <- rbind(fData(x), as.data.frame(matrix(NA, nrow = nrow(th2[[1]]), 
                                                           ncol = ncol(fData(x)),
                                                           dimnames = list(rownames(th2[[1]]), colnames(fData(x))))))
            x2 <- ExpressionSet(
              assayData   = new_exprs,
              phenoData   = phenoData(x),
              featureData = AnnotatedDataFrame(fd_new),
              annotation  = annotation(x)
            )
          } else {
            x2 <- ExpressionSet(
              assayData   = assayData(x)[[keep]],
              phenoData   = phenoData(x),
              featureData = featureData(x),
              annotation  = annotation(x)
            )
          }
          return(x2)
        })
        stopCluster(cl)
        
        # Get the data for the meta-PCA train datasets (as defined in the config file):
        pathways.train <- pathways[config$dataset_id[which(config$ccm_meta_pca == 1)]]
        
        # Train the meta PCA
        metaPCA <- service_meta_pca(data = pathways.train, n_pc = 3)
        
        # Project the samples (all the datasets) on the meta PCA
        metaPCA.projected <- lapply(pathways, function(x) service_meta_pca_projection(x, meta_pca = metaPCA))
        
        return(list(projected = metaPCA.projected, metaPCA = metaPCA, training = pathways.train))
      })
      names(allTerms) = unique(treatment_ws$term)
      return(allTerms)
    })
    names(allWS) = c("bulk","adjusted__1__1")
    return(allWS)
  })
  names(allSubmodels) = c("yesDown","yesUp")
  return(allSubmodels)
}, 
ccm_wfid = "wf-08a6a0a503",
treatment_ws = treatment_ws,
memory_request = "60Gi",
image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest", 
tags = list(list(name="analysis",value="pathway_metaPCA_treatments")))
# wf-24d130f812

## Extraction
## ==========================================
keyPathways_wfid = "wf-24d130f812"
metaPCA_pathways = readRDS(get_workflow_outputs(keyPathways_wfid))
metaPCA_pathways = lapply(names(metaPCA_pathways), function(x){
  m = purrr::transpose(metaPCA_pathways[[x]])
  names(m) = paste0(names(m),"_",x)
  return(m)
}) %>% unlist(recursive = F)


## Scores
## ------------
extract_sampleScores = function(metaPCA.projected, flipPC1, flipPC2) {
  pathwayPCA <- cytoreason.ccm.pipeline:::statistic_table.ccm_service_meta_pca_projection(metaPCA.projected)
  sampleScores <- pathwayPCA$sample_loadings
  sampleScores_mat <- reshape(sampleScores[,c("submodel","pc","sample_id","value")], idvar = c("submodel","sample_id"), timevar = "pc", direction = "wide")
  colnames(sampleScores_mat) <- gsub("value.","pathway_meta_",colnames(sampleScores_mat))
  colnames(sampleScores_mat)[colnames(sampleScores_mat)=="submodel"] <- "dataset"
  
  # Check the direction of the meta-PCs first - I needed to flip the sign of the meta-PC1 for both cells and pathways (because sample loadings of HC>DZ):
  if(flipPC1) { sampleScores_mat$pathway_meta_pc1 <- (-1) * sampleScores_mat$pathway_meta_pc1 }
  if(flipPC2) { sampleScores_mat$pathway_meta_pc2 <- (-1) * sampleScores_mat$pathway_meta_pc2 }
  return(sampleScores_mat)
}

process = function(res, term, submodel, collection) {
  res$term = term
  res$submodel = submodel
  res$collection = "treatment"
  return(res)
}


# all terms and submodel
sampleScores_all = lapply(names(metaPCA_pathways), function(term){
  lapply(names(metaPCA_pathways[[term]]), function(submodel){
    cat("\n", term, submodel,"\n")
    res = extract_sampleScores(metaPCA_pathways[[term]][[submodel]]$projected, T, T)
    res = process(res, term, submodel, collection = NULL)
    cat("\r", term, submodel,"........done")
    return(res)
  }) %>% bind_rows()
}) %>% bind_rows()

pushToCC(sampleScores_all, tagsToPass = list(list(name="object",value="pathway_meta_pca_treatment")))
# wf-5c0d927f2c

## Pathway loadings
## --------------------
extract_pathwayLoadings = function(metaPCA, flipPC1, flipPC2){
  metaPCA.res = cytoreason.ccm.pipeline:::statistic_table.ccm_service_meta_pca(list(metaPCA))
  pathLoadings = metaPCA.res$feature_loadings
  pathLoadings$pc = str_replace(pathLoadings$pc, "pc", "PC")
  pathLoadings = pathLoadings[which(pathLoadings$pc %in% c("PC1","PC2","PC3")),]
  if(flipPC1) { pathLoadings$value[which(pathLoadings$pc == "PC1")] = (-1) * pathLoadings$value[which(pathLoadings$pc == "PC1")] }
  if(flipPC2) { pathLoadings$value[which(pathLoadings$pc == "PC2")] = (-1) * pathLoadings$value[which(pathLoadings$pc == "PC2")] } # only in bulk
  return(pathLoadings)
}


# all terms and submodel
pathwayLoadings_all = lapply(names(metaPCA_pathways), function(term){
  lapply(names(metaPCA_pathways[[term]]), function(submodel){
    cat("\n", term, submodel,"\n")
    flipPC2 = ifelse(submodel == "bulk", T, F)
    res = extract_pathwayLoadings(metaPCA_pathways[[term]][[submodel]]$metaPCA, T, flipPC2)
    res = process(res, term, submodel, collection = NULL)
    cat("\r", term, submodel,"........done")
    return(res)
  }) %>% bind_rows()
}) %>% bind_rows()

pathwayLoadings_all = unique(pathwayLoadings_all)
pathLoadings.BQ = pathwayLoadings_all %>%
  dplyr::select(-c(1:3))
colnames(pathLoadings.BQ) = c("Pathway","PC","Loading","Term","Submodel","Collection")
pushToCC(pathLoadings.BQ, tagsToPass = list(list(name="object",value="pathway_meta_pca_loadings_treatment")))
# wf-b6f6092a16
uploadToBQ(pathLoadings.BQ, bqdataset = "s05_atopic_dermatitis", tableName = "pathwayLoadings", disposition = "WRITE_APPEND")



## ========================================
## Define R/NR white space
## =========================================
treatments = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "treatment_PrePostDupi")
treatments = treatments %>%
  dplyr::filter(str_detect(term, "W16:R_L_vs_NR_L|W4:R_L_vs_NR_L")) %>%
  dplyr::filter(submodel %in% c("bulk","adjusted__1__1") & collection %in% c("h","kegg","reactome","btm","th2","neuroinflammation")) %>%
  mutate(submodel = ifelse(submodel == "bulk", "bulk", "adjusted"))


pathways_in_dz = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "AD_gx_gsa")
pathways_in_dz = pathways_in_dz %>%
  group_by(term, submodel) %>%
  dplyr::filter(submodel %in% c("bulk","adjusted") & collection %in% c("h","kegg","reactome","btm","th2","neuroinflammation")) %>%
  dplyr::select(term, submodel, collection, pathway, NES, FDR, log10_fdr) %>%
  mutate(dir = ifelse(FDR > 0.05, "unchanged",
                      ifelse(NES >0, "up","down"))) %>%
  ungroup()

dz_wide = pathways_in_dz %>%
  pivot_wider(id_cols = c(pathway,collection, submodel),
              names_from  = term,
              values_from = dir)


treatment_ws = treatments %>%
  merge(., dz_wide, all = T, by = c("pathway","collection","submodel")) %>%
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

treatment_ws = treatment_ws %>%
  dplyr::filter(in_white_space_L_vs_HC != "no") %>%
  mutate(pathway_id = paste0(collection,":",pathway)) %>%
  dplyr::select(term, pathway_id, pathway, submodel, collection, in_white_space_L_vs_HC)

pushToCC(treatment_ws, tagsToPass = list(list(name="object",value="treatment_RNR_ws")))
# wf-ed106cda0c


# Calcaulate PCA
# =============================================================
ccm = as_ccm_fit("wf-08a6a0a503")

run_function_dist(FUN = function(ccm_wfid, treatment_ws){
  library(cytoreason.ccm.pipeline)
  library(cytoreason.integration)
  library(dplyr)
  install.packages("parallel")
  library(parallel)
  
  ccm = as_ccm_fit(ccm_wfid)
  config <- modelMetadata(ccm)
  allSubmodels = lapply(c("yesDown","yesUp"), function(ws){
    allWS = lapply(c("bulk","adjusted__1__1"), function(submodel){
      allTerms = lapply(unique(treatment_ws$term), function(term){
        cat("\nComputing:",ws,submodel,term,"\n")
        
        cl = makeCluster(length(ccm$datasets))
        clusterEvalQ(cl, { library(cytoreason.ccm.pipeline) })
        clusterExport(cl, varlist = c("ccm","ws","submodel","term","treatment_ws","nrow","ncol"), envir = environment())
        
        pathways <- parLapply(cl, ccm$datasets, function(d){
          x = analysisResultExpressionSet(d, "gene_set_activity")
          th2_enrichments = readRDS(get_workflow_outputs("wf-b2132fba22"))
          keep = ifelse(submodel == "bulk", "exprs", paste0("exprs_",submodel))
          if(submodel != "bulk") { submodel = "adjusted" }
          paths = treatment_ws[which(treatment_ws$term == term & 
                                       treatment_ws$submodel == submodel & 
                                       treatment_ws$in_white_space_L_vs_HC == ws),]
          
          x = x[which(rownames(x) %in% paths$pathway_id),]
          
          if("th2" %in% unique(paths$collection)){
            th2 = th2_enrichments[[keep]]
            th2 = lapply(th2, function(y) return(y[which(rownames(y) %in% paths$pathway),,drop=F]))
            new_exprs = rbind(assayData(x)[[keep]], th2[[x@experimentData@title]])
            
            fd_new <- rbind(fData(x), as.data.frame(matrix(NA, nrow = nrow(th2[[1]]), 
                                                           ncol = ncol(fData(x)),
                                                           dimnames = list(rownames(th2[[1]]), colnames(fData(x))))))
            x2 <- ExpressionSet(
              assayData   = new_exprs,
              phenoData   = phenoData(x),
              featureData = AnnotatedDataFrame(fd_new),
              annotation  = annotation(x)
            )
          } else {
            x2 <- ExpressionSet(
              assayData   = assayData(x)[[keep]],
              phenoData   = phenoData(x),
              featureData = featureData(x),
              annotation  = annotation(x)
            )
          }
          return(x2)
        })
        stopCluster(cl)
        
        # Get the data for the meta-PCA train datasets (as defined in the config file):
        pathways.train <- pathways[config$dataset_id[which(config$ccm_meta_pca == 1)]]
        
        # Train the meta PCA
        metaPCA <- service_meta_pca(data = pathways.train, n_pc = 3)
        
        # Project the samples (all the datasets) on the meta PCA
        metaPCA.projected <- lapply(pathways, function(x) service_meta_pca_projection(x, meta_pca = metaPCA))
        
        return(list(projected = metaPCA.projected, metaPCA = metaPCA, training = pathways.train))
      })
      names(allTerms) = unique(treatment_ws$term)
      return(allTerms)
    })
    names(allWS) = c("bulk","adjusted__1__1")
    return(allWS)
  })
  names(allSubmodels) = c("yesDown","yesUp")
  return(allSubmodels)
}, 
ccm_wfid = "wf-08a6a0a503",
treatment_ws = treatment_ws,
memory_request = "60Gi",
image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest", 
tags = list(list(name="analysis",value="pathway_metaPCA_treatments_RNR")))
# wf-4a6e6051b6

## Extraction
## ==========================================
keyPathways_wfid = "wf-4a6e6051b6"
metaPCA_pathways = readRDS(get_workflow_outputs(keyPathways_wfid))
metaPCA_pathways = lapply(names(metaPCA_pathways), function(x){
  m = purrr::transpose(metaPCA_pathways[[x]])
  names(m) = paste0(names(m),"_",x)
  return(m)
}) %>% unlist(recursive = F)


## Scores
## ------------
extract_sampleScores = function(metaPCA.projected, flipPC1, flipPC2) {
  pathwayPCA <- cytoreason.ccm.pipeline:::statistic_table.ccm_service_meta_pca_projection(metaPCA.projected)
  sampleScores <- pathwayPCA$sample_loadings
  sampleScores_mat <- reshape(sampleScores[,c("submodel","pc","sample_id","value")], idvar = c("submodel","sample_id"), timevar = "pc", direction = "wide")
  colnames(sampleScores_mat) <- gsub("value.","pathway_meta_",colnames(sampleScores_mat))
  colnames(sampleScores_mat)[colnames(sampleScores_mat)=="submodel"] <- "dataset"
  
  # Check the direction of the meta-PCs first - I needed to flip the sign of the meta-PC1 for both cells and pathways (because sample loadings of HC>DZ):
  if(flipPC1) { sampleScores_mat$pathway_meta_pc1 <- (-1) * sampleScores_mat$pathway_meta_pc1 }
  if(flipPC2) { sampleScores_mat$pathway_meta_pc2 <- (-1) * sampleScores_mat$pathway_meta_pc2 }
  return(sampleScores_mat)
}

process = function(res, term, submodel, collection) {
  res$term = term
  res$submodel = submodel
  res$collection = "treatment_RNR"
  return(res)
}


# all terms and submodel
sampleScores_all = lapply(names(metaPCA_pathways), function(term){
  lapply(names(metaPCA_pathways[[term]]), function(submodel){
    cat("\n", term, submodel,"\n")
    res = extract_sampleScores(metaPCA_pathways[[term]][[submodel]]$projected, T, T)
    res = process(res, term, submodel, collection = NULL)
    cat("\r", term, submodel,"........done")
    return(res)
  }) %>% bind_rows()
}) %>% bind_rows()

pushToCC(sampleScores_all, tagsToPass = list(list(name="object",value="pathway_meta_pca_RNR")))
# wf-0109132b74

## Pathway loadings
## --------------------
extract_pathwayLoadings = function(metaPCA, flipPC1, flipPC2){
  metaPCA.res = cytoreason.ccm.pipeline:::statistic_table.ccm_service_meta_pca(list(metaPCA))
  pathLoadings = metaPCA.res$feature_loadings
  pathLoadings$pc = str_replace(pathLoadings$pc, "pc", "PC")
  pathLoadings = pathLoadings[which(pathLoadings$pc %in% c("PC1","PC2","PC3")),]
  if(flipPC1) { pathLoadings$value[which(pathLoadings$pc == "PC1")] = (-1) * pathLoadings$value[which(pathLoadings$pc == "PC1")] }
  if(flipPC2) { pathLoadings$value[which(pathLoadings$pc == "PC2")] = (-1) * pathLoadings$value[which(pathLoadings$pc == "PC2")] } # only in bulk
  return(pathLoadings)
}


# all terms and submodel
pathwayLoadings_all = lapply(names(metaPCA_pathways), function(term){
  lapply(names(metaPCA_pathways[[term]]), function(submodel){
    cat("\n", term, submodel,"\n")
    flipPC2 = ifelse(submodel == "bulk", T, F)
    res = extract_pathwayLoadings(metaPCA_pathways[[term]][[submodel]]$metaPCA, T, flipPC2)
    res = process(res, term, submodel, collection = NULL)
    cat("\r", term, submodel,"........done")
    return(res)
  }) %>% bind_rows()
}) %>% bind_rows()

pathwayLoadings_all = unique(pathwayLoadings_all)
pathLoadings.BQ = pathwayLoadings_all %>%
  dplyr::select(-c(1:3))
colnames(pathLoadings.BQ) = c("Pathway","PC","Loading","Term","Submodel","Collection")
pushToCC(pathLoadings.BQ, tagsToPass = list(list(name="object",value="pathway_meta_pca_loadings_RNR")))
# wf-e884e49bb6
uploadToBQ(pathLoadings.BQ, bqdataset = "s05_atopic_dermatitis", tableName = "pathwayLoadings", disposition = "WRITE_APPEND")
