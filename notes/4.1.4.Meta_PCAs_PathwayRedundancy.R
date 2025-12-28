devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(parallel)
library(tidyverse)
set.seed(1234)

ccm = as_ccm_fit("wf-3e419ff83b")

pathwayRedundancy = openxlsx::read.xlsx("~/analysis-s05/data/pc_loadings_redund_final.xlsx")
# wf-6cda062e07
nonRedundantPathways = pathwayRedundancy$Pathway[pathwayRedundancy$final.pathways.one.per.module]

## 1. Export gene activity scores for the relevant pathways
## =================================================================
run_function_dist(FUN = function(ccm_wfid, nonRedundantPathways){
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
      
      cl = makeCluster(length(ccm$datasets))
      clusterEvalQ(cl, {   library(cytoreason.ccm.pipeline) })
      clusterExport(cl, varlist = c("ccm","nonRedundantPathways","effect","submodel", "keyPathways"), envir = environment())
      
      pathways <- parLapply(cl, ccm$datasets, function(d){
        # choose key pathways in disease
        x = analysisResultExpressionSet(d, "gene_set_activity")
        x = x[keyPathways$ID,]
        keep = ifelse(submodel == "bulk", "exprs", paste0("exprs_",submodel))
        
        # load relevant th2 pathways
        th2_enrichments = readRDS(get_workflow_outputs("wf-b2132fba22"))
        if(submodel == "bulk" & effect == "AD") {
          th2 <- lapply(th2_enrichments[[keep]], function(y) return(y[-which(rownames(y) == "NTF4"),,drop=T]))
        } else {
          th2 = th2_enrichments[[keep]]
        }
        
        # create new matrix
        new_exprs = rbind(assayData(x)[[keep]], th2[[x@experimentData@title]])
        new_exprs = new_exprs[which(rownames(new_exprs) %in% nonRedundantPathways),,drop=F]
        
        fd_new <- as.data.frame(matrix(NA, nrow = nrow(new_exprs), 
                                       ncol = ncol(fData(x)),
                                       dimnames = list(rownames(new_exprs), colnames(fData(x)))))
        
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
      metaPCA <- service_meta_pca(data = pathways.train, n_pc = 2)
      
      # Project the samples (all the datasets) on the meta PCA
      metaPCA.projected <- lapply(pathways, function(x) service_meta_pca_projection(x, meta_pca = metaPCA))
      
      return(list(projected = metaPCA.projected, metaPCA = metaPCA, training = pathways.train))
    })
    names(allSubmodels) = c("bulk","adjusted__1__1")
    return(allSubmodels)
  })
  names(allKeyPathways) = c("AD","L_vs_HC")
  return(allKeyPathways)
}, 
ccm_wfid = "wf-08a6a0a503",
nonRedundantPathways = nonRedundantPathways,
memory_request = "60Gi",
image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest", 
tags = list(list(name="analysis",value="pathwayRedundancy_metaPCA_allTerms")))
# wf-5029aba66f
# wf-4674be8ff1


## Extraction
## ==========================================
keyPathways_wfid = "wf-4674be8ff1"
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
    res$collection = paste0("Unredundant_KeyPathways_in_",term)
  }
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
pushToCC(sampleScores_all, tagsToPass = list(list(name="analysis",value="unredundant_pathway_meta_pca")))
# wf-64e555ed8f


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
uploadToBQ(pathLoadings.BQ, bqdataset = "s05_atopic_dermatitis", tableName = "pathwayLoadings", disposition = "WRITE_APPEND")
pushToCC(pathLoadings.BQ, tagsToPass = list(list(name="analysis",value="unredundant_pathway_meta_pca_loadings")))
# wf-2236a96045

