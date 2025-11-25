`### Using the results from the custom CCM run, we will extract the result tables we are interested
### in, and customize the tables to have a uniform built.

# 1. Libraries and data loading
# ------------------------------------------
devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline) # make sure to load version >= 1.0.1
library(cytoreason.assets)
library(tidyverse)

ccm_wfid = "wf-8e948630d7"
ccm = as_ccm_fit(ccm_wfid)
targetMapping = read_asset("wf-6c4b539629") # taken from "2.7.Final_Signatures.R"
geneCollections = unique(targetMapping$Target.Collection)
geneMapping = toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
collectionMapping = readRDS(get_workflow_outputs("wf-6b6125369c"))

# 2. Extraction
# ----------------------------
# We use the L_vs_NL effect to extract results because we need only the information in within lesional samples
# and the L_vs_NL effect has the highest number of datasets
# this stage needs a lot of RAM because of the large number of signatures, so we do it remotely
run_function_dist(FUN = function(ccm_wfid){
  library(cytoreason.ccm.pipeline)
  library(dplyr)
  library(stringr)

  ccm = as_ccm_fit(ccm_wfid)
  geneMapping = toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
  collections = names(ccm$parameters$custom_gene_set_collections)

  # Enrichment in disease (done sepearaly because it weirdly fails within the function wf-17e2e4ef3c)
  # ------------------------------
  DZEnrichment = lapply(c("L_vs_NL","L_vs_HC","NL_vs_HC","AD"), function(effect){
    cat("\nEffect:",effect,"\n")
    x = build_service_result_tables(ccm$meta[[effect]]$gx_diff$gx_gsa,
                                               subset = ~(submodel %in% c("bulk", "adjusted__1__1") &
                                                            term == ifelse(effect == "AD", "DZ_vs_HC",effect))) # enrichment in disease - all effects
    leadingEdge = x$gsa_leading_edge
    leadingEdge$feature_id = geneMapping$symbol[match(leadingEdge$feature_id, geneMapping$gene_id)]
    leadingEdge = leadingEdge %>%
      group_by(hit_id) %>%
      summarise(feature_id = paste(feature_id, collapse = ";"), .groups = "drop")
    x = x$gsa_enrichment
    x$hit = leadingEdge$feature_id[match(x$hit_id,leadingEdge$hit_id)]
    return(x)
  })
  DZEnrichment = bind_rows(DZEnrichment)

  # Pheno-feature correlations (molecular score/clinical score/meta pcas) in lesional samples
  # -----------------------------------------------------------------------------------------------
  phenoFeature = build_service_result_tables(ccm$meta$L_vs_NL$pheno_feature_correlations,
                                      subset = ~(submodel %in% c("bulk", "adjusted__1__1") & term == "L")) # MS/CS - L
  phenoFeature = phenoFeature$gene_set_phenotypic_variable_correlations


  # Correlation to cells/pathways/genes in lesional samples
  # ------------------------------------------------------------
  # this is a very large table so we filter it ahead of time
  TargetFeatures = build_service_result_tables(ccm$meta$L_vs_NL$feature_geneset_correlations,
                                               subset = ~(submodel %in% c("bulk", "adjusted__1__1") &
                                                            term == "L" &
                                                            !collection_1 %in% c(collections, "c2.cp.kegg", "c2.cp.reactome")))
  TargetFeatures = bind_rows(TargetFeatures$cell_contribution_geneset_correlations,
                             TargetFeatures$gene_geneset_correlations,
                             TargetFeatures$gene_set_activity_geneset_correlations)


  Results = list(DZEnrichment = DZEnrichment, TargetFeatures = TargetFeatures, phenoFeature = phenoFeature)
  return(Results)
},
ccm_wfid = ccm_wfid,
memory_request = "100Gi",
image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:master_latest")
# wf-213be9d205
# wf-4a73197f02

Results = readRDS(get_workflow_outputs("wf-4a73197f02"))


# 3. Editing to keep only what we're interested in
# ----------------------------------------------------
# This function filters the tables to the terms we want to keep:
# L vs NL for enrichment in disease but in calculations within lesional samples for the rest of the criteria
# And unifies the column names - what's related to the signatures is under 'Target'
# and what's related to the criteria is under 'Criteria'.

run_function_dist(FUN = function(results_wfid){
  library(cytoreason.cc.client)
  library(dplyr)
  library(stringr)

  Results = readRDS(get_workflow_outputs(results_wfid))
  targetMapping = readRDS(get_workflow_outputs("wf-6c4b539629")) # taken from "2.7.Final_Signatures.R"
  geneCollections = unique(targetMapping$Target.Collection)
  collectionMapping = readRDS(get_workflow_outputs("wf-6b6125369c"))


  processResults = function(data, tableName) {
    if(tableName == "DZEnrichment") {
      targetID = "pathway"
      criteriaID = "term"
      colnames(data)[which(colnames(data) == "nes")] <- "value"
      data$DataType = paste0("Enrichment in ",data$term)
      data$pathway = paste0(data$collection, ":",data$pathway)
      terms = c("L_vs_NL","L_vs_HC","NL_vs_HC","AD","DZ_vs_HC")
    } else if(tableName == "TargetFeatures") {
      targetID = "feature_id_2"
      criteriaID = "feature_id_1"
      data$feature_type1 = str_replace(data$feature_type1, "gene_set", "pathway")
      data$DataType = paste0("Target_", str_to_title(data$feature_type1))
      data$collection = str_split(data$feature_id_2,":",simplify = TRUE)[,1]
      data$hit = NA
      terms = "L"
    } else {
      targetID = "feature_id_1"
      criteriaID = "feature_id_2"
      mapping = c("EASI" = "CS", "SCOARD" = "CS", "MolecularScore" = "MS",
                  "cell_meta_pc1" = "Cell_PCA", "cell_meta_pc2" = "Cell_PCA", "cell_meta_pc3" = "Cell_PCA",
                  "pathway_meta_pc1" = "Pathway_PCA", "pathway_meta_pc2" = "Pathway_PCA", "pathway_meta_pc3" = "Pathway_PCA",
                  "adj_pathway_meta_pc1" = "Pathway_PCA", "adj_pathway_meta_pc2" = "Pathway_PCA", "adj_pathway_meta_pc3" = "Pathway_PCA")
      data$DataType = paste0("Target_", sapply(data$feature_id_2, function(x) mapping[match(x,names(mapping))]))
      data$collection = str_split(data$feature_id_1,":",simplify = TRUE)[,1]
      data$hit = NA
      terms = "L"
    }

    newtable = data %>%
      dplyr::filter(term %in% terms) %>%
      dplyr::rename(Target.Identifier = targetID) %>%
      mutate(Target.ID = Target.Identifier) %>%
      mutate(Target.Identifier = str_split(Target.Identifier,":",simplify = TRUE)[,2]) %>%
      dplyr::rename(Target.Collection = collection) %>%
      dplyr::filter(Target.Collection %in% geneCollections) %>%
      mutate(Criteria.Identifier = get(criteriaID)) %>%
      mutate(Criteria.Collection = str_split(Criteria.Identifier,":",simplify = TRUE)[,1]) %>%
      mutate(log10_fdr = -log10(fdr)) %>%
      mutate(metricType = case_when(tableName == "DZEnrichment" ~ "NES", .default = "bi-weight mid-correlation")) %>%
      dplyr::rename(metricValue = value, Type = submodel) %>%
      dplyr::select(DataType, Target.ID, Target.Identifier, Target.Collection, Type, Criteria.Identifier, Criteria.Collection, metricValue, metricType, fdr, log10_fdr, hit)

    # Summarize random and shuffled random
    for(id in c("random","top50","bottom50")){
      idx = which(str_detect(newtable$Target.Identifier,id))
      if(id != "random") { id = paste0("smoothedRandom_",id)}
      
      tmp = newtable[idx,] %>%
        group_by(DataType, Type, Criteria.Identifier, Criteria.Collection, metricType, Target.Collection) %>%
        summarise(fdr = mean(fdr),
                  metricValue = mean(metricValue)) %>%
        mutate(log10_fdr = -log10(fdr)) %>%
        mutate(Target.Identifier = id, Target.ID = paste0("negativeControls:",id), hit = NA) %>%
        select(colnames(newtable))

      newtable = newtable[-idx,]
      newtable = rbind(newtable, tmp)
      rm(idx, tmp)
    }

    newtable$Target.Collection = collectionMapping$collection[match(newtable$Target.Identifier, collectionMapping$signature)]
    return(newtable)
  }

  Results = lapply(names(Results), function(dataName) processResults(data = Results[[dataName]], tableName = dataName))
  return(Results)
}, results_wfid = "wf-4a73197f02")
# wf-f0959d74d1
# wf-921cf942af

Results = readRDS(get_workflow_outputs("wf-921cf942af"))
Results = list(DZEnrichment = Results[[1]],
               Target_Cell = Results[[2]][which(Results[[2]]$DataType == "Target_Cell"),],
               Target_Gene = Results[[2]][which(Results[[2]]$DataType == "Target_Gene"),],
               Target_Pathway = Results[[2]][which(Results[[2]]$DataType == "Target_Pathway"),],
               Target_MS = Results[[3]][which(Results[[3]]$DataType == "Target_MS"),],
               Target_CS = Results[[3]][which(Results[[3]]$DataType == "Target_CS"),],
               Target_Cell_PCA = Results[[3]][which(Results[[3]]$DataType == "Target_Cell_PCA"),],
               Target_Pathway_PCA = Results[[3]][which(Results[[3]]$DataType == "Target_Pathway_PCA"),]
               )

uploadToBQ(Results$DZEnrichment, bqdataset = "s05_atopic_dermatitis", tableName = "Results_DZEnrichment")
uploadToBQ(Results$Target_Cell, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_Cell")
uploadToBQ(Results$Target_Gene, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_Gene")
uploadToBQ(Results$Target_Pathway, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_Pathway")
uploadToBQ(Results$Target_MS, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_MS")
uploadToBQ(Results$Target_CS, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_CS")
uploadToBQ(Results$Target_Cell_PCA, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Cell_PCA")
uploadToBQ(Results$Target_Pathway_PCA, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Pathway_PCA")
