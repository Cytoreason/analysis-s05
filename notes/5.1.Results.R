`### Using the results from the custom CCM run, we will extract the result tables we are interested
### in, and customize the tables to have a uniform built.

# 1. Libraries and data loading
# ------------------------------------------
devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline) # make sure to load version >= 1.0.1
library(tidyverse)

ccm_wfid = "wf-abde4bfab0"

collectionMapping = readRDS(get_workflow_outputs("wf-3df1530237"))
  collectionMapping$collection[str_detect(collectionMapping$collection, "Negative|negative")] <- "negativeControls"
  collectionMapping$ID = paste0(collectionMapping$collection,":",collectionMapping$signature)
  collectionMapping$previousID = paste0("X2:",collectionMapping$signature)
  collectionMapping$previousID[which(collectionMapping$collection == "negativeControls")] <- paste0("negativeControls:",collectionMapping$signature[which(collectionMapping$collection == "negativeControls")])
  collectionMapping$previousID[which(collectionMapping$signature %in% c("BMP4","BMP6","BMP7","TGFB2"))] <- paste0("X2:",collectionMapping$signature[which(collectionMapping$signature %in% c("BMP4","BMP6","BMP7","TGFB2"))])
signatures = readRDS(get_workflow_outputs("wf-a24a2031e8"))
signatures = signatures[-which(str_detect(signatures$Old_identifier,"x2cov")),]
signatureMapping = merge(collectionMapping, signatures, by.x = "signature", by.y = "Old_identifier", all=T)
signatureMapping[which(signatureMapping$signature == "mast_wu"),] = c("mast_wu","Mast","Mast:mast_wu", "X2:mast_wu","Mast","mast_wu")

pushToCC(signatureMapping)
# wf-fcb1ef2bbd
# wf-aa75ed069b


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
                                                            !collection_1 %in% c("c2.cgp","c2.cp.biocarta","c2.cp","c3.tft","c7","c2.cp.kegg","c2.cp.reactome")))
  TargetFeatures = bind_rows(TargetFeatures$cell_contribution_geneset_correlations,
                             TargetFeatures$gene_geneset_correlations,
                             TargetFeatures$gene_set_activity_geneset_correlations)


  Results = list(DZEnrichment = DZEnrichment, TargetFeatures = TargetFeatures, phenoFeature = phenoFeature)
  return(Results)
},
ccm_wfid = ccm_wfid,
memory_request = "25Gi",
image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:master_latest")
# wf-213be9d205
# wf-4a73197f02
# wf-d4494fb9b7
# wf-3c1660f432
# wf-a09930b845
# wf-66f108f616
# wf-345bf31e4e
# wf-182d0335f5

Results = readRDS(get_workflow_outputs("wf-345bf31e4e"))


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

  results_wfid = get_workflow(results_wfid, wait = T)
  Results = readRDS(get_workflow_outputs(results_wfid))
  signatureMapping = readRDS(get_workflow_outputs("wf-aa75ed069b"))

  processResults = function(data, tableName) {
    if(tableName == "DZEnrichment") {
      targetID = "pathway"
      criteriaID = "term"
      colnames(data)[which(colnames(data) == "nes")] <- "value"
      data$DataType = paste0("Enrichment in ",data$term)
      idx = which(data$feature_id %in% signatureMapping$previousID)
      data$collection[idx] = signatureMapping$collection[match(data$feature_id[idx], signatureMapping$previousID)]
      data$pathway[str_detect(data$pathway,"reactome:|kegg:|h:")] = str_replace(data$pathway[str_detect(data$pathway,"reactome:|kegg:|h:")],":","__")
      data$pathway = paste0(data$collection, ":",data$pathway)
      terms = c("L_vs_NL","L_vs_HC","NL_vs_HC","AD","DZ_vs_HC")
    } else if(tableName == "TargetFeatures") {
      targetID = "feature_id_2"
      criteriaID = "feature_id_1"
      data$feature_type1 = str_replace(data$feature_type1, "gene_set", "pathway")
      data$DataType = paste0("Target_", str_to_title(data$feature_type1))
      idx = which(data$feature_id_2 %in% signatureMapping$previousID)
      data$feature_id_2[idx] = signatureMapping$ID[match(data$feature_id_2[idx], signatureMapping$previousID)]
      data$collection[idx] = signatureMapping$collection[match(data$feature_id_2[idx], signatureMapping$ID)]
      data$collection[-idx] = str_split(data$feature_id_2[-idx],":",simplify = TRUE)[,1]
      data$hit = NA
      terms = "L"
    } else {
      targetID = "feature_id_1"
      criteriaID = "feature_id_2"
      idx = which(data$feature_id_1 %in% signatureMapping$previousID)
      data$feature_id_1[idx] = signatureMapping$ID[match(data$feature_id_1[idx], signatureMapping$previousID)]
      path = unique(data$feature_id_2)[str_detect(unique(data$feature_id_2),"pathway")]
        names(path) = rep("Pathway_PCA", length(path))
      cell = unique(data$feature_id_2)[str_detect(unique(data$feature_id_2),"meta1_pc")]
        names(cell) = rep("Cell_PCA", length(cell))
      mapping = c("EASI" = "CS", "SCORAD" = "CS", "MolecularScore" = "MS",
                  setNames(object = names(path), nm = path),
                  setNames(object = names(cell), nm = cell))
      data$DataType = paste0("Target_", sapply(data$feature_id_2, function(x) mapping[match(x,names(mapping))]))
      data$collection = str_split(data$feature_id_1,":",simplify = TRUE)[,1]
      data$hit = NA
      terms = "L"
    }

    data[,targetID] = str_replace(data[,targetID],":reactome:",":reactome__")
    data[,targetID] = str_replace(data[,targetID],":kegg:",":kegg__")
    data[,targetID] = str_replace(data[,targetID],":h:",":h__")
    data[,targetID] = str_replace(data[,targetID],":btm:",":btm__")
    data[,criteriaID] = str_replace(data[,criteriaID],":reactome:",":reactome__")
    data[,criteriaID] = str_replace(data[,criteriaID],":kegg:",":kegg__")
    data[,criteriaID] = str_replace(data[,criteriaID],":h:",":h__")
    data[,criteriaID] = str_replace(data[,criteriaID],":btm:",":btm__")
    
    newtable = data %>%
      dplyr::filter(term %in% terms) %>%
      dplyr::rename(Target.Identifier = targetID) %>%
      mutate(Target.ID = Target.Identifier) %>%
      mutate(Target.Identifier = str_split(Target.Identifier,":",simplify = TRUE)[,2]) %>%
      dplyr::rename(Target.Collection = collection) %>%
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
    
    # mapping targets to new collections
    idx = which(newtable$Target.Identifier %in% signatureMapping$signature)
    newtable$Target.Collection[idx] = signatureMapping$collection[match(newtable$Target.Identifier[idx], signatureMapping$signature)]

    # renaming target signatures and removing unwanted ones
    newtable$Target.Identifier[idx] = signatureMapping$New_identifier[match(newtable$Target.Identifier[idx], signatureMapping$signature)]

    return(newtable)
  }
  
  Results = lapply(names(Results), function(dataName) processResults(data = Results[[dataName]], tableName = dataName))
  return(Results)
}, 
results_wfid = "wf-182d0335f5",
memory_request = "25Gi")
# wf-f0959d74d1
# wf-921cf942af
# wf-7f81e4e7a2
# wf-7bfb6831ed
# wf-a30ed87a96
# wf-c60a87ef45
# wf-d37ae65109
# wf-edaf4d1f2d
# wf-9579f7fe45 - without correlation to cell meta pcs
# wf-3c09cc81dd
# wf-62e5f0bcba
# wf-ef478fadc9 (final version)

Results = readRDS(get_workflow_outputs(get_workflow("wf-ef478fadc9", wait = T)))
Results = list(DZEnrichment = Results[[1]],
               Target_Cell = Results[[2]][which(Results[[2]]$DataType == "Target_Cell"),],
               Target_Gene = Results[[2]][which(Results[[2]]$DataType == "Target_Gene"),],
               Target_Pathway = Results[[2]][which(Results[[2]]$DataType == "Target_Pathway"),],
               Target_MS = Results[[3]][which(Results[[3]]$DataType == "Target_MS"),],
               Target_CS = Results[[3]][which(Results[[3]]$DataType == "Target_CS"),],
               Target_Cell_PCA = Results[[3]][which(Results[[3]]$DataType == "Target_Cell_PCA"),],
               Target_Pathway_PCA = Results[[3]][which(Results[[3]]$DataType == "Target_Pathway_PCA"),])

# Changing *criteria* collection for X2 signatures
idx = which(Results$Target_Pathway$Criteria.Collection == "X2")
Results$Target_Pathway$Criteria.Collection[idx] = signatureMapping$collection[match(Results$Target_Pathway$Criteria.Identifier[idx], signatureMapping$previousID)]
Results$Target_Pathway$Criteria.Identifier[idx] = signatureMapping$ID[match(Results$Target_Pathway$Criteria.Identifier[idx], signatureMapping$previousID)]
idx = which(Results$Target_Pathway$Criteria.Collection == "negativeControls")
Results$Target_Pathway = Results$Target_Pathway[-idx,]

pushToCC(Results, tagsToPass = list(list(name="object",value="processed_results")))
# wf-798227871e
# wf-91e85b6c8e
# wf-6f20463d54
# wf-ba85ebeacf
# wf-a53211cea7
# wf-cd3c366c62
# wf-1bea2eb00f
# wf-64470d2a55
# wf-47f2a1c1b7 (final version)



# For all tables but Results_Pathway_PCA we use the previous version in order to have Th1_Related, Th17_Related, Th2_Related
uploadToBQ(Results$Target_Pathway_PCA, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Pathway_PCA")

Results_old = readRDS(get_workflow_outputs("wf-64470d2a55"))
uploadToBQ(Results_old$DZEnrichment, bqdataset = "s05_atopic_dermatitis", tableName = "Results_DZEnrichment")
uploadToBQ(Results_old$Target_Cell, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_Cell")
uploadToBQ(Results_old$Target_Gene, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_Gene")
uploadToBQ(Results_old$Target_Pathway, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_Pathway")
uploadToBQ(Results_old$Target_MS, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_MS")
uploadToBQ(Results_old$Target_CS, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_CS")
uploadToBQ(Results_old$Target_Cell_PCA, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Cell_PCA")


## Outside run of meta correlations to SCORAD
## -------------------------------------------------
signatureMapping = readRDS(get_workflow_outputs("wf-aa75ed069b"))
ccm_scorad = as_ccm_fit("wf-b1d96b1738") # the meta analysis is not good enough
# SCORAD = build_service_result_tables(ccm_scorad$meta$L_vs_NL$pheno_feature_correlations, term = "L",
#            subset = ~(collection %in% c("h","kegg","reactome","btm","X2","epidermis","neuroinflammation","negativeControls","th2")))

SCORAD = lapply(names(ccm_scorad$datasets), function(d){
  x = statistic_table(analysisResultElement(ccm_scorad$datasets[[d]]$model[[1]],"pheno_feature_correlations"), 
    term = "L", subset = ~(feature_type1 == "gene_set"))
  x = x[-which(x$collection %in% c("c2.cgp","c2.cp.biocarta","c2.cp","c3.tft","c7","c2.cp.kegg","c2.cp.reactome")), ]
  
  SCORAD = x %>%
    dplyr::rename(Target.Identifier = feature_id_1) %>%
    mutate(Target.ID = Target.Identifier) %>%
    mutate(Target.Identifier = str_split(Target.Identifier,":",simplify = TRUE)[,2]) %>%
    dplyr::rename(Target.Collection = collection, Criteria.Identifier = feature_id_2) %>%
    mutate(Criteria.Identifier = "SCORAD", Criteria.Collection = "SCORAD", DataType = "Target_CS") %>%
    mutate(log10_fdr = -log10(fdr)) %>%
    mutate(metricType = "bi-weight mid-correlation") %>%
    dplyr::rename(metricValue = value, Type = submodel) %>%
    dplyr::select(DataType, Target.ID, Target.Identifier, Target.Collection, Type, Criteria.Identifier, Criteria.Collection, metricValue, metricType, fdr, log10_fdr, n_observation)
  
  # Summarize random and shuffled random
  for(id in c("random","top50","bottom50")){
    idx = which(str_detect(SCORAD$Target.Identifier,id))
    if(id != "random") { id = paste0("smoothedRandom_",id)}
    
    tmp = SCORAD[idx,] %>%
      group_by(DataType, Type, Criteria.Identifier, Criteria.Collection, metricType, Target.Collection,n_observation) %>%
      summarise(fdr = mean(fdr),
                metricValue = mean(metricValue)) %>%
      mutate(log10_fdr = -log10(fdr)) %>%
      mutate(Target.Identifier = id, Target.ID = paste0("negativeControls:",id), hit = NA) %>%
      dplyr::select(colnames(SCORAD))
    
    SCORAD = SCORAD[-idx,]
    SCORAD = rbind(SCORAD, tmp)
    rm(idx, tmp)
  }
  
  idx = which(SCORAD$Target.Collection == "X2")
  SCORAD$Target.Collection[idx] = signatureMapping$collection[match(SCORAD$Target.Identifier[idx], signatureMapping$signature)]
  SCORAD$Target.Identifier[idx] = signatureMapping$New_identifier[match(SCORAD$Target.Identifier[idx], signatureMapping$signature)]
  SCORAD$Target.Identifier[which(SCORAD$Target.ID == "X2:mast_wu")] = "mast_wu"
  
  SCORAD$dataset = d
  return(SCORAD)
}) %>% do.call(rbind,.)
  
uploadToBQ(SCORAD, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_SCORAD")

samplesUsed = designSampleGroupContrasts("L_vs_NL__GSE130588__GPL570__Lesion_vs_non_lesion", data = ccm_scorad$datasets$GSE130588__GPL570)[["L_vs_NL"]]
  samplesUsed = samplesUsed$sample_id[which(samplesUsed$label == "L")]
targets = signatureMapping$previousID[which(signatureMapping$signature %in% c("IL4","IL13","x2_general_inhibition_early_50"))]

targetCorrelations = analysisResultExpressionSet(ccm_scorad$datasets$GSE130588__GPL570,"gene_set_activity")
scorad = targetCorrelations@phenoData@data
  scorad = scorad[which(scorad$sample_id %in% samplesUsed),c("sample_id","SCORAD2")]
targetCorrelations = assayData(targetCorrelations)[["exprs"]]
  targetCorrelations = targetCorrelations[targets,samplesUsed]
  rownames(targetCorrelations) = signatureMapping$New_identifier[match(rownames(targetCorrelations),signatureMapping$previousID)]
  targetCorrelations = reshape2::melt(targetCorrelations)
  colnames(targetCorrelations) = c("Target","sample_id","enrichment")

targetCorrelations$SCORAD = scorad$SCORAD2[match(targetCorrelations$sample_id, scorad$sample_id)]  

ggplot(targetCorrelations, aes(x = SCORAD, y = enrichment)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~Target) +
  ggpubr::stat_cor(method = "spearman") +
  theme_minimal() +
  ggpubr::border() +
  labs(y = "Enrichment Scores")


## Target Pathway Random Threshold
## ======================================
Results = readRDS(get_workflow_outputs("wf-345bf31e4e"))
Pathways = Results[["TargetFeatures"]] %>%
  dplyr::filter(feature_type1 == "gene_set") %>%
  dplyr::filter(submodel == "bulk")  

nc = Pathways %>%
  dplyr::filter(collection %in% c("h","kegg","reactome","btm","neuroinflammation","th2")) %>%
  dplyr::filter(str_detect(feature_id_2, "negativeControls:")) %>%
  mutate(feature_id_2 = ifelse(str_detect(feature_id_2, "negativeControls:random"), "random",
                               ifelse(str_detect(feature_id_2, 'top50'), "smoothedRandom_top50","smoothedRandom_bottom50")))

nc_thresholds = nc %>%
  group_by(feature_id_2) %>%
  summarise(q = list(quantile(value, probs = c(seq(0.1,0.9,0.1),seq(0.91,1,0.01))))) %>%
  tidyr::unnest_longer(q, indices_to = "quantile", values_to = "correlation") %>%
  mutate(quantile = factor(quantile, ordered = T, levels = unique(nc_thresholds$quantile)))

ggplot(nc_thresholds, aes(x = quantile, y = correlation, color = feature_id_2, group = feature_id_2)) +
  geom_line() +
  scale_y_continuous(breaks = round(seq(-0.3,1,0.1),1)) +
  geom_hline(yintercept = c(0.4, 0.45), linetype = "dashed") +
  geom_vline(xintercept = c("98%","99%"), linetype = "dashed") +
  annotate(x = "93%", y=0.41, geom = "text", label = "Correlation threshold of 0.40") +
  annotate(x = "93%", y=0.46, geom = "text", label = "Correlation threshold of 0.45") +
  theme_minimal() +
  ggpubr::border()+
  labs(x = "Quantile", y = "Threshold Correlation", color = "Random Collection",
       title = "Quantiles of random correlations - based on 339,700 correlations each",
       subtitle = "100 random signatures per collection X 3397 pathways from hallmark, kegg, reactome , btm, neuroinflammation and th2")
ggsave("~/analysis-s05/figures/Results/random_correlation_quantiles.png", units = "px", height = 2000, width = 3000, bg = "white")
