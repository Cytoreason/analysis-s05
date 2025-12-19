devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)

# calculation is done on the disease model as we don't need the signatures
ccm <- as_ccm_fit("wf-08a6a0a503")

## Calculate key pathways PCA
## ==========================================
run_function_dist(FUN = function(ccm_wfid){
  library(cytoreason.ccm.pipeline)
  library(cytoreason.integration)
  library(dplyr)
  install.packages("parallel")
  library(parallel)
  
  ccm = as_ccm_fit(ccm_wfid)
  config <- modelMetadata(ccm)
  allKeyPathways = lapply(unique(names(ccm$meta)[1:4]), function(effect){
    allSubmodels = lapply(c("bulk","adjusted__1__1"), function(submodel){
      cat("\nComputing:",effect,submodel,"\n")
      keyPathways <- build_service_result_tables(ccm$meta[[effect]]$gx_diff$gx_gsa) # we use DZ vs HC to keep track with cell meta pca
      keyPathways <- keyPathways[[submodel]]$gsa_enrichment %>% .[c("h", "btm", "kegg", "reactome")]
      
      keyPathways = do.call(rbind, 
                            lapply(names(keyPathways), function(nm) {
                              df <- keyPathways[[nm]]$gsa_enrichment
                              df$collection <- nm
                              if (effect == "AD") { df = df[which(df$term == "DZ_vs_HC"), ]
                              } else { df = df[which(df$term == effect), ] }
                              df  }))
      keyPathways <- keyPathways %>% 
        dplyr::filter(FDR <= 0.05) %>% 
        dplyr::mutate(ID = paste0(collection,":", pathway)) %>% 
        dplyr::select(ID)
      
      # Get pathway ssgsea for all the samples in all datasets:
      cl = makeCluster(length(ccm$datasets))
      clusterEvalQ(cl, { 
        library(cytoreason.ccm.pipeline)
        library(dplyr)
      })
      clusterExport(cl, varlist = c("ccm","keyPathways"), envir = environment())
      
      pathways <- parLapply(cl, ccm$datasets, function(d){
        analysisResultExpressionSet(d, "gene_set_activity") %>%
          .[keyPathways$ID,]
      })
      stopCluster(cl)
      
      # Choose the correct assay
      keep = ifelse(submodel == "bulk", "exprs", paste0("exprs_",submodel))
      
      m <- lapply(pathways, function(x) assayDataElement(x, keep))
      
      # Replace assayData with ONLY that assay (same name preserved)
      pathways = lapply(names(pathways), function(x) {
        d = pathways[[x]]
        assayData(d) <- do.call(assayDataNew, c(list("environment"), 
                                                setNames(list(m[[x]]), "exprs")))
        return(d)
      })
      names(pathways) = names(ccm$datasets)
      
      # Get the data for the meta-PCA train datasets (as defined in the config file):
      pathways.train <- pathways[config$dataset_id[which(config$ccm_meta_pca == 1)]]
      
      # Train the meta PCA
      metaPCA <- service_meta_pca(data = pathways.train, n_pc = 3)
      
      # Project the samples (all the datasets) on the meta PCA
      metaPCA.projected <- lapply(pathways, function(x) service_meta_pca_projection(x, meta_pca = metaPCA))
      
      # now per collection
      cat("\nPer Collection\n")
      perCollection = lapply(c("h","btm","kegg","reactome"), function(collection){
        pathw = lapply(pathways, function(x) x[stringr::str_detect(rownames(x),paste0(collection,":"))])
        pathways.train <- pathw[config$dataset_id[which(config$ccm_meta_pca == 1)]]
        metaPCA <- service_meta_pca(data = pathways.train, n_pc = 3)
        metaPCA.projected <- lapply(pathw, function(x) service_meta_pca_projection(x, meta_pca = metaPCA))
        return(list(projected = metaPCA.projected, metaPCA = metaPCA, training = pathways.train, keyPathways = keyPathways))
      })
      names(perCollection) = c("h","btm","kegg","reactome")
      
      return(list(projected = metaPCA.projected, metaPCA = metaPCA, training = pathways.train, keyPathways = keyPathways, perCollection = perCollection))
      
    })
    names(allSubmodels) = c("bulk","adjusted__1__1")
    return(allSubmodels)
  })
  names(allKeyPathways) = unique(names(ccm$meta)[1:4])
  return(allKeyPathways)
}, 
ccm_wfid = "wf-08a6a0a503", # new disease model, including adjustment to keratinocytes
image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest", 
tags = list(list(name="analysis",value="pathway_metaPCA_allTerms")))
# wf-ee9a8f05a7 - bulk
# wf-98fc25d966 - adjusted
# wf-0a6f338228 - all terms
# wf-1b16e72f76 - all terms, fixed bulk/adjustment subsetting


## Extraction
## ==========================================
keyPathways_wfid = "wf-ee9a8f05a7" # bulk
keyPathways_wfid = "wf-98fc25d966" # adjusted
keyPathways_wfid = "wf-1b16e72f76" # all
metaPCA_pathways = readRDS(get_workflow_outputs(keyPathways_wfid))


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
  if(!is.null(collection)){
    res$collection = collection
  } else {
    res$collection = paste0("KeyPathways_in_",term)
  }
  return(res)
}


# all terms and submodel
sampleScores_all = lapply(names(metaPCA_pathways), function(term){
  lapply(names(metaPCA_pathways[[term]]), function(submodel){
    cat("\n", term, submodel,"\n")
    res = extract_sampleScores(metaPCA_pathways[[term]][[submodel]]$projected, T, T)
    res = process(res, term, submodel, collection = NULL)
    res_perCollection = lapply(names(metaPCA_pathways[[term]][[submodel]][["perCollection"]]), function(collection){
      cat("\r", term, submodel,"...",collection)
      res_col = extract_sampleScores(metaPCA_pathways[[term]][[submodel]][["perCollection"]][[collection]]$projected, T, T)
      res_col = process(res_col, term, submodel, collection)
      return(res_col)
    }) %>% bind_rows()
    cat("\r", term, submodel,"........done")
    return(bind_rows(res, res_perCollection))
  }) %>% bind_rows()
}) %>% bind_rows()
# wf-e8f514f3f1
# wf-5b90576980 - with keratinocytes

pushToCC(sampleScores_all, tagsToPass = list(list(name="analysis",value="pathway_meta_pca")))
# wf-5b90576980
# wf-39ae62c8d2


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
    res_perCollection = lapply(names(metaPCA_pathways[[term]][[submodel]][["perCollection"]]), function(collection){
      cat("\r", term, submodel,"...",collection)
      res_col = extract_pathwayLoadings(metaPCA_pathways[[term]][[submodel]][["perCollection"]][[collection]]$metaPCA, T, flipPC2)
      res_col = process(res_col, term, submodel, collection)
      return(res_col)
    }) %>% bind_rows()
    cat("\r", term, submodel,"........done")
    return(bind_rows(res, res_perCollection))
  }) %>% bind_rows()
}) %>% bind_rows()

pathwayLoadings_all = unique(pathwayLoadings_all)
pathLoadings.BQ = pathwayLoadings_all %>%
  dplyr::select(-c(1:3))
colnames(pathLoadings.BQ) = c("Pathway","PC","Loading","Term","Submodel","Collection")
pushToCC(pathLoadings.BQ, tagsToPass = list(list(name="analysis",value="pathway_meta_pca_loadings")))
# wf-c397bd39fd
# wf-36fef22586
uploadToBQ(pathLoadings.BQ, bqdataset = "s05_atopic_dermatitis", tableName = "pathwayLoadings")


