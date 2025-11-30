library(cytoreason.ccm.pipeline)
library(tidyverse)

ccm = as_ccm_fit("wf-08a6a0a503")
geneMapping = toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
# uploadToBQ(geneMapping, bqdataset = "s05_atopic_dermatitis", tableName = "geneMapping")

# 1. gx_gsa
# ------------------------------------
gxgsa = rbind(statistic_table(ccm$meta$L_vs_NL$gx_diff$gx_gsa),
              statistic_table(ccm$meta$AD$gx_diff$gx_gsa),
              statistic_table(ccm$meta$L_vs_HC$gx_diff$gx_gsa),
              statistic_table(ccm$meta$NL_vs_HC$gx_diff$gx_gsa))
gxgsa = gxgsa %>%
  dplyr::filter(term %in% c("L_vs_NL","DZ_vs_HC","L_vs_HC","NL_vs_HC"))

uploadToBQ(gxgsa, bqdataset = "s05_atopic_dermatitis", tableName = "AD_gx_gsa")


# 2. gx_diff
# ------------------------------------
gx_diff = rbind(statistic_table(ccm$meta$L_vs_NL$gx_diff$gx_diff),
              statistic_table(ccm$meta$AD$gx_diff$gx_diff),
              statistic_table(ccm$meta$L_vs_HC$gx_diff$gx_diff),
              statistic_table(ccm$meta$NL_vs_HC$gx_diff$gx_diff))
gx_diff = gx_diff %>%
  dplyr::filter(term %in% c("L_vs_NL","DZ_vs_HC","L_vs_HC","NL_vs_HC"))

uploadToBQ(gx_diff, bqdataset = "s05_atopic_dermatitis", tableName = "AD_gx_diff")


# 2. ct_test
# ------------------------------------
ct_test = rbind(statistic_table(ccm$meta$L_vs_NL$ct_test),
                statistic_table(ccm$meta$AD$ct_test),
                statistic_table(ccm$meta$L_vs_HC$ct_test),
                statistic_table(ccm$meta$NL_vs_HC$ct_test))
ct_test = ct_test %>%
  dplyr::filter(term %in% c("L_vs_NL","DZ_vs_HC","L_vs_HC","NL_vs_HC"))

uploadToBQ(ct_test, bqdataset = "s05_atopic_dermatitis", tableName = "AD_ct_test")



# 4. Disease-Model related graphs
# --------------------------------------
unified_metadata = read_asset("wf-e82a5ab3b6")
unified_metadata$sample_classification[which(unified_metadata$sample_classification == "Normal")] <- "HC"
# uploadToBQ(unified_metadata, bqdataset = "s05_atopic_dermatitis", tableName = "AD_sample_metadata") # is further edited in the meta PC section

sampleClassifications.inventory = unified_metadata %>%
  dplyr::filter(is.na(time) | time %in% c("D0","D1")) %>%
  group_by(dataset_id_short, sample_classification) %>%
  summarise(
    nSamples = n(),
    SCORAD = sum(!is.na(clinical_score_value_SCORAD)),
    EASI = sum(!is.na(clinical_score_value_EASI))
  ) %>%
  mutate(clinical_score = if_else(SCORAD > 0 | EASI >0, "Yes", "No")) %>%
  mutate(sample_classification = factor(sample_classification, ordered = T, levels = c("Lesion","Non Lesion","HC"))) %>%
  ungroup()
sampleClassifications.inventory$notes = ifelse(sampleClassifications.inventory$dataset_id_short == "GSE32473",
                                               "EASI unavailable, but reduction in EASI is available","")
uploadToBQ(sampleClassifications.inventory, bqdataset = "s05_atopic_dermatitis", tableName = "AD_sampleClassifications_inventory")


# 5. Correlation between cells and pathways/genes
# ----------------------------------------------------
feature_cell = build_service_result_tables(ccm$meta$AD$feature_cell_correlations, assay_element = c("exprs","exprs_adjusted__1__1"))

feature_cell = lapply(names(feature_cell[["gene_set_activity"]]), function(model){
  lapply(names(feature_cell[["gene_set_activity"]][[model]]), function(samples){
    df = feature_cell[["gene_set_activity"]][[model]][[samples]]
    df[,c("collection","pathway")] = do.call(rbind, lapply(df$feature_id_1, function(x) strsplit(x, "__")[[1]]))
    df$model = ifelse(model == "exprs", "bulk", "adjusted")
    df$samples = samples
    df = df[,c(1,10,11,2:9,12:13)]
    rownames(df) = NULL
    return(df)
  }) %>% do.call(rbind,.)
}) %>% do.call(rbind,.)

feature_cell = feature_cell[-which(feature_cell$samples == "__all__"),]
uploadToBQ(feature_cell, bqdataset = "s05_atopic_dermatitis", tableName = "AD_feature_cell_correlations")


gene_cell = build_service_result_tables(ccm$meta$AD$gene_cell_correlations, assay_element = c("exprs","exprs_adjusted__1__1"))

gene_cell = lapply(names(gene_cell), function(model){
  lapply(names(gene_cell[[model]]), function(samples){
    df = gene_cell[[model]][[samples]]
    df$symbol = geneMapping$symbol[match(df$feature_id_1, geneMapping$gene_id)]
    df$model = ifelse(model == "exprs", "bulk", "adjusted")
    df$samples = samples
    df = df[,c(1,10,2:9,11:12)]
    rownames(df) = NULL
    return(df)
  }) %>% do.call(rbind,.)
}) %>% do.call(rbind,.)

gene_cell = gene_cell[-which(gene_cell$samples == "__all__")]
uploadToBQ(gene_cell, bqdataset = "s05_atopic_dermatitis", tableName = "AD_gene_cell_correlations")


# 6. Disease network edges
# --------------------------------
dz_network = read_asset("wf-e315350209") #taken from Shiran's analysis in p03 https://github.com/Cytoreason/analysis-p03-poc1/blob/master/notes/fileExport/NetworkEdges.R
dz_network = dz_network %>%
  dplyr::filter(EdgeType == "AD.bulk")
uploadToBQ(dz_network, bqdataset = "s05_atopic_dermatitis", tableName = "AD_dz_network_edges")


# 6. ct_test in different terms
# ---------------------------------------
treatmentTerms = read.csv("~/data/Res_noterms.long.csv", row.names = 1)
treatmentTerms = treatmentTerms[-which(treatmentTerms$Terms %in% c("meta_pc1","keratinocyte")),]
meta = statistic_table(ccm$meta$Treatment$ct_test)

# splitting into chunks because it's a model with many terms (138)
treatmentTerms.edit <- treatmentTerms[treatmentTerms$model == "Dupilumab_response__GSE130588__GPL570__R_vs_NR", ]
chunk_size <- 13
n <- nrow(treatmentTerms.edit)
num_chunks <- ceiling(n / chunk_size)

for (i in 1:num_chunks) {
  start_idx <- (i - 1) * chunk_size + 1
  end_idx <- min(i * chunk_size, n)
  treatmentTerms.edit$model[start_idx:end_idx] <- paste0("Dupilumab_response__GSE130588__GPL570__R_vs_NR-", i)
}

treatmentTerms.edit <- rbind(
  treatmentTerms.edit,
  treatmentTerms[treatmentTerms$model != "Dupilumab_response__GSE130588__GPL570__R_vs_NR", ]
)

cytoreason.cc.client::run_function_dist(FUN = function(treatmentTerms, ccm_wf){
  library(cytoreason.cc.client)
  library(data.table)
  library(dplyr)

  res_wf = get_workflow_id(get_workflow(
    mapply_dist(FUN = function(model, treatmentTerms, ccm_wf){
      library(cytoreason.ccm.pipeline)

      ccm = as_ccm_fit(ccm_wf)
      treatmentTerms = treatmentTerms[which(treatmentTerms$model == model),]
      dataset = treatmentTerms[1,"data_set"]
      terms = treatmentTerms$Terms
      model = strsplit(model, "-")[[1]][1]
      cat("\nExtracting\n")
      ct_test = statistic_table(analysisResultElement(ccm[["datasets"]][[dataset]][["model"]][[model]], name = "ct_test"))
      gx_diff = statistic_table(analysisResultElement(ccm[["datasets"]][[dataset]][["model"]][[model]], name = "gx_diff"))
      gx_gsa = statistic_table(analysisResultElement(ccm[["datasets"]][[dataset]][["model"]][[model]], name = "gx_gsa"))
      
      cat("\nListing\n")
      tests = lapply(list(ct_test,gx_diff,gx_gsa), function(test){
        test$term = as.character(test$term)
        test$meta = "no"
        test = test[which(test$term %in% terms),]
        return(test)
      })
      names(tests) = c("ct_test","gx_diff","gx_gsa")
      return(tests)
    },
    model = unique(treatmentTerms$model),
    MoreArgs = list(treatmentTerms = treatmentTerms, ccm_wf = "wf-08a6a0a503"),
    memory_request = "6Gi",
    force_execution = F,
    image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_0.81.9"
    )
  , wait = T))
  
  cat("\nReading results\n")
  res = cytoreason.assets::read_asset(res_wf)
  cat("\n..ct_tests\n")
  ct_tests = lapply(res, function(x) x[["ct_test"]]) %>% as.data.table() %>% rbindlist()
  cat("\n..gx_diffs\n")
  gx_diffs = lapply(res, function(x) x[["ct_test"]]) %>% as.data.table() %>% rbindlist()
  cat("\n..gx_gsas\n")
  gx_gsas = lapply(res, function(x) x[["ct_test"]]) %>% as.data.table() %>% rbindlist()
  
  saveRDS(ct_tests, "./output/ct_tests.rds")
  saveRDS(gx_diffs, "./output/gx_diffs.rds")
  saveRDS(gx_gsas, "./output/gx_gsas.rds")
  cat("\nDone!\n")
},
treatmentTerms = treatmentTerms.edit,
ccm_wf = "wf-08a6a0a503",
image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_0.81.9", 
tags = list(list(name = "task",value="extract_treatment_data"),
            list(name="ccm_wfid",value="wf-08a6a0a503"),
            list(name="disease",value="AD"))
)
# wf-8aee685b90