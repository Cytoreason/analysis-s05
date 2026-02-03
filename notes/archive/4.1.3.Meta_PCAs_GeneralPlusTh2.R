devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)

# calculation is done on the disease model as we don't need the signatures
# from th2 and epidermis we remove the pathways that are cannonical because they are already represented
th2_enrichments = readRDS(get_workflow_outputs("wf-7ed6eff409"))
th2_enrichments$th2 = lapply(th2_enrichments$th2, function(x) lapply(x, function(y) { y=y[-c(13:19),] ; return(y) }))
th2_enrichments$epidermis = lapply(th2_enrichments$epidermis, function(x) lapply(x, function(y) { y=y[-c(1:19),] ; return(y) }))
th2_enrichments = lapply(th2_enrichments, function(x) purrr::transpose(x)) %>% purrr::transpose()

# rbinding the matrices
th2_enrichments$exprs <- Map(rbind, th2_enrichments$exprs$th2, 
                             th2_enrichments$exprs$neuroinflammation, 
                             th2_enrichments$exprs$epidermis)
th2_enrichments$exprs_adjusted__1__1 <- Map(rbind, th2_enrichments$exprs_adjusted__1__1$th2, 
                                            th2_enrichments$exprs_adjusted__1__1$neuroinflammation,
                                            th2_enrichments$exprs_adjusted__1__1$epidermis)
paths = rownames(th2_enrichments$exprs$GSE107361__GPL570)

enrichedInDisease = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "AD_gx_gsa")
enrichedInDisease = enrichedInDisease[-which(enrichedInDisease$collection %in% c("Neuronal","Itch")),]
enrichedInDisease = enrichedInDisease %>%
  dplyr::filter(pathway %in% paths) %>%
  dplyr::filter(term %in% c("L_vs_HC","DZ_vs_HC")) %>%
  dplyr::filter(!collection %in% c("Neuronal","Itch")) %>%
  dplyr::filter(FDR <= 0.05) %>% 
  dplyr::mutate(ID = paste0(collection,":", pathway)) %>%
  dplyr::select(submodel, term, ID, pathway)

pushToCC(enrichedInDisease, tags = list(list(name="object",value="th2_pathways_enriched_in_disease")))
# wf-fbb04343d3

enrichedInDisease = split(enrichedInDisease, enrichedInDisease$submodel)

th2_enrichments$exprs = lapply(th2_enrichments$exprs, function(x) x[which(rownames(x) %in% enrichedInDisease$bulk$pathway),])
th2_enrichments$exprs_adjusted__1__1 = lapply(th2_enrichments$exprs_adjusted__1__1, function(x) x[which(rownames(x) %in% enrichedInDisease$adjusted$pathway),])
# We see that in adjusted the pathways are the same for both terms, and in bulk DZ_vs_HC is missing neuroinflammation:NTF4

pushToCC(th2_enrichments, tags = list(list(name="object",value="th2_pathways")))
# wf-b2132fba22

## Calculate key pathways PCA
## ==========================================
run_function_dist(FUN = function(ccm_wfid, th2_enrichments_wfid){
  library(cytoreason.ccm.pipeline)
  library(cytoreason.integration)
  library(dplyr)
  install.packages("parallel")
  library(parallel)
  
  ccm = as_ccm_fit(ccm_wfid)
  config <- modelMetadata(ccm)
  allKeyPathways = lapply(c("AD","L_vs_HC"), function(effect){
    allSubmodels = lapply(c("bulk","adjusted__1__1"), function(submodel){
      cat("\nComputing:",effect,submodel,"\n")
      keyPathways <- build_service_result_tables(ccm$meta[[effect]]$gx_diff$gx_gsa)
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
      
      # We rbind the new th2 pathways with the existing cannonical pathways
      # Get pathway ssgsea for all the samples in all datasets:
      cl = makeCluster(length(ccm$datasets))
      clusterEvalQ(cl, { 
        library(cytoreason.ccm.pipeline)
        library(dplyr)
      })
      clusterExport(cl, varlist = c("ccm","keyPathways","submodel","effect","th2_enrichments_wfid"), envir = environment())
      
      pathways <- parLapply(cl, ccm$datasets, function(d){
        x = analysisResultExpressionSet(d, "gene_set_activity")
        x = x[keyPathways$ID,]
        keep = ifelse(submodel == "bulk", "exprs", paste0("exprs_",submodel))
        
        th2_enrichments = readRDS(get_workflow_outputs("wf-b2132fba22"))
        
        if(submodel == "bulk" & effect == "AD") {
          th2 <- lapply(th2_enrichments[[keep]], function(y) return(y[-which(rownames(y) == "NTF4"),]))
        } else {
          th2 = th2_enrichments[[keep]]
        }
        
        new_exprs = rbind(assayData(x)[[keep]], th2[[x@experimentData@title]])
        fd_new <- rbind(fData(x), as.data.frame(matrix(NA, nrow = nrow(th2[[x@experimentData@title]]), 
                                                       ncol = ncol(fData(x)),
                                                       dimnames = list(rownames(th2[[x@experimentData@title]]), colnames(fData(x))))))

        x2 <- ExpressionSet(
          assayData   = new_exprs,
          phenoData   = phenoData(x),
          featureData = AnnotatedDataFrame(fd_new),
          annotation  = annotation(x)
        )
        return(x2)
      })
      stopCluster(cl)
      
      # Get the data for the meta-PCA train datasets (as defined in the config file):
      pathways.train <- pathways[config$dataset_id[which(config$ccm_meta_pca == 1)]]
      
      # Train the meta PCA
      metaPCA <- service_meta_pca(data = pathways.train, n_pc = 3)
      
      # Project the samples (all the datasets) on the meta PCA
      metaPCA.projected <- lapply(pathways, function(x) service_meta_pca_projection(x, meta_pca = metaPCA))

      return(list(projected = metaPCA.projected, metaPCA = metaPCA, training = pathways.train, keyPathways = keyPathways))
    })
    names(allSubmodels) = c("bulk","adjusted__1__1")
    return(allSubmodels)
  })
  names(allKeyPathways) = c("AD","L_vs_HC")
  return(allKeyPathways)
}, 
ccm_wfid = "wf-08a6a0a503",
th2_enrichments_wfid = "wf-b2132fba22",
memory_request = "60Gi",
image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest", 
tags = list(list(name="analysis",value="pathway_metaPCA_allTerms_plusTh2")))
# wf-113e6defe3
# wf-1233696939 - rearrange lists
# wf-e23dc30c35 - using only th2 that change in disease


## Extraction
## ==========================================
keyPathways_wfid = "wf-e23dc30c35"
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
  res$collection = "keyWithTh2"
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

pushToCC(sampleScores_all, tagsToPass = list(list(name="analysis",value="pathway_meta_pca_keyWithTh2")))
# wf-6dba8d483f
# wf-c4b6c8d817
# wf-3fce829a35

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
pushToCC(pathLoadings.BQ, tagsToPass = list(list(name="object",value="pathway_meta_pca_loadings_keyWithTh2")))
# wf-454d12aa31
# wf-880b3f8483
# wf-25cea04750
uploadToBQ(pathLoadings.BQ, bqdataset = "s05_atopic_dermatitis", tableName = "pathwayLoadings", disposition = "WRITE_APPEND")


