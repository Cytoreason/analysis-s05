devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)

ccm_dm = as_ccm_fit("wf-08a6a0a503") # the original disease model
ccm_cgs = as_ccm_fit("wf-3e419ff83b") # adding keratinocyte adjustment
geneMapping = toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
# uploadToBQ(geneMapping, bqdataset = "s05_atopic_dermatitis", tableName = "geneMapping")

signatureMapping = readRDS(get_workflow_outputs("wf-aa75ed069b"))


# 1. gx_gsa
# ------------------------------------
# From disease model
gxgsa_dm = rbind(statistic_table(ccm_dm$meta$L_vs_NL$gx_diff$gx_gsa),
              statistic_table(ccm_dm$meta$AD$gx_diff$gx_gsa),
              statistic_table(ccm_dm$meta$L_vs_HC$gx_diff$gx_gsa),
              statistic_table(ccm_dm$meta$NL_vs_HC$gx_diff$gx_gsa))
gxgsa_dm = gxgsa_dm %>%
  dplyr::filter(term %in% c("L_vs_NL","DZ_vs_HC","L_vs_HC","NL_vs_HC")) %>%
  dplyr::filter(submodel %in% c("bulk","adjusted__1__1")) %>%
  mutate(submodel = case_when(submodel == "bulk" ~ "bulk",
                              submodel == "adjusted__1__1" ~ "adjusted"))


# From custom gene set
gxgsa_cgs = rbind(statistic_table(ccm_cgs$meta$L_vs_NL$gx_diff$gx_gsa),
              statistic_table(ccm_cgs$meta$AD$gx_diff$gx_gsa),
              statistic_table(ccm_cgs$meta$L_vs_HC$gx_diff$gx_gsa),
              statistic_table(ccm_cgs$meta$NL_vs_HC$gx_diff$gx_gsa))
gxgsa_cgs = gxgsa_cgs %>%
  dplyr::filter(term %in% c("L_vs_NL","DZ_vs_HC","L_vs_HC","NL_vs_HC")) %>%
  dplyr::filter(submodel %in% c("bulk","adjusted__1__1")) %>%
  mutate(submodel = case_when(submodel == "bulk" ~ "bulk",
                              submodel == "adjusted__1__1" ~ "adjusted")) %>%
  rename(FDR = fdr, ES = es, NES = nes, nMoreExtreme = nmoreextreme) %>%
  dplyr::select(-feature_id,-hit_id)


gxgsa = bind_rows(gxgsa_dm, gxgsa_cgs) %>%
  dplyr::select(-hit, hit) #moving it to be the last column

# Summarize random and shuffled random
for(id in c("random","top50","bottom50")){
  idx = which(str_detect(gxgsa$pathway,id))
  if(id != "random") { id = paste0("smoothedRandom_",id)}
  
  tmp = gxgsa[idx,] %>%
    group_by(submodel,model,term,collection,method) %>%
    summarise(across(c(pvalue:estimate_percentile), mean, na.rm = TRUE)) %>% 
    mutate(log10_fdr = -log10(FDR), log10_pvalue = -log10(pvalue)) %>%
    mutate(pathway = id, hit = NA) %>%
    select(colnames(gxgsa))
  
  gxgsa = gxgsa[-idx,]
  gxgsa = rbind(gxgsa, tmp)
}

idx = which(gxgsa$collection == "X2")
gxgsa$collection[idx] = signatureMapping$collection[match(gxgsa$pathway[idx], signatureMapping$signature)]
gxgsa$nMoreExtreme=as.integer(gxgsa$nMoreExtreme)
gxgsa$size=as.integer(gxgsa$size)

uploadToBQ(gxgsa, bqdataset = "s05_atopic_dermatitis", tableName = "AD_gx_gsa")


# 2. gx_diff
# ------------------------------------
gx_diff = rbind(statistic_table(ccm$meta$L_vs_NL$gx_diff$gx_diff),
              statistic_table(ccm$meta$AD$gx_diff$gx_diff),
              statistic_table(ccm$meta$L_vs_HC$gx_diff$gx_diff),
              statistic_table(ccm$meta$NL_vs_HC$gx_diff$gx_diff))

gx_diff = gx_diff %>%
  dplyr::filter(term %in% c("L_vs_NL","DZ_vs_HC","L_vs_HC","NL_vs_HC")) %>%
  mutate(submodel = case_when(submodel == "bulk" ~ "bulk",
                              submodel == "adjusted__1__meta1_pc1" ~ "adjusted_metaPC1",
                              submodel == "adjusted__1__CRCL_0000348" ~ "adjusted_keratinocyte",
                              submodel == "adjusted__1__meta1_pc1_CRCL_0000348" ~ "adjusted_metaPC1_keratinocytes"))

uploadToBQ(gx_diff, bqdataset = "s05_atopic_dermatitis", tableName = "AD_gx_diff")


# 2. ct_test
# ------------------------------------
cells = read.csv("~/data/skin.csv", sep = "\t", header = F)
ct_test = rbind(statistic_table(ccm$meta$L_vs_NL$ct_test),
                statistic_table(ccm$meta$AD$ct_test),
                statistic_table(ccm$meta$L_vs_HC$ct_test),
                statistic_table(ccm$meta$NL_vs_HC$ct_test))

ct_test = ct_test %>%
  mutate(feature_id = cells$V1[match(feature_id, cells$V2)])

uploadToBQ(ct_test, bqdataset = "s05_atopic_dermatitis", tableName = "AD_ct_test")



# 4. Correlation between cells and pathways/genes
# ----------------------------------------------------
# geneset/cell
# ----------------
subsetting = ~(collection %in% c("h","kegg","reactome","btm") &
                 submodel %in% c("bulk","adjusted__1__meta1_pc1","adjusted__1__CRCL_0000348","adjusted__1__meta1_pc1_CRCL_0000348"))

geneset_cell = rbind(build_service_result_tables(ccm$meta$L_vs_NL$feature_cell_correlations, subset = subsetting)$gene_set_activity_cell_correlations %>%
                       mutate(term = "L_vs_NL"),
                     build_service_result_tables(ccm$meta$AD$feature_cell_correlations, subset = subsetting)$gene_set_activity_cell_correlations %>%
                       mutate(term = "DZ_vs_HC"),
                     build_service_result_tables(ccm$meta$L_vs_HC$feature_cell_correlations, subset = subsetting)$gene_set_activity_cell_correlations %>%
                       mutate(term = "L_vs_HC"),
                     build_service_result_tables(ccm$meta$NL_vs_HC$feature_cell_correlations, subset = subsetting)$gene_set_activity_cell_correlations %>%
                       mutate(term = "NL_vs_HC"))

geneset_cell = geneset_cell %>%
  dplyr::filter(sample_subset != "L_vs_NL") %>% # unclear why this is included
  mutate(feature_id_2 = cells$V1[match(feature_id_2, cells$V2)]) %>%
  mutate(submodel = case_when(submodel == "bulk" ~ "bulk",
                              submodel == "adjusted__1__meta1_pc1" ~ "adjusted_metaPC1",
                              submodel == "adjusted__1__CRCL_0000348" ~ "adjusted_keratinocyte",
                              submodel == "adjusted__1__meta1_pc1_CRCL_0000348" ~ "adjusted_metaPC1_keratinocytes"))

uploadToBQ(geneset_cell, bqdataset = "s05_atopic_dermatitis", tableName = "AD_geneset_cell_correlations")


# gene/cell
# ----------------
subsetting = ~(submodel %in% c("bulk","adjusted__1__meta1_pc1","adjusted__1__CRCL_0000348","adjusted__1__meta1_pc1_CRCL_0000348"))

gene_cell = rbind(build_service_result_tables(ccm$meta$L_vs_NL$gene_cell_correlations, subset = subsetting)$gene_cell_correlations,
                 build_service_result_tables(ccm$meta$AD$gene_cell_correlations, subset = subsetting)$gene_cell_correlations,
                 build_service_result_tables(ccm$meta$L_vs_HC$gene_cell_correlations, subset = subsetting)$gene_cell_correlations,
                 build_service_result_tables(ccm$meta$NL_vs_HC$gene_cell_correlations, subset = subsetting)$gene_cell_correlations)

gene_cell = gene_cell %>%
  dplyr::filter(sample_subset != "L_vs_NL") %>% # unclear why this is included
  mutate(feature_id_2 = cells$V1[match(feature_id_2, cells$V2)]) %>%
  mutate(symbol = geneMapping$symbol[match(feature_id_2, geneMapping$gene_id)]) %>%
  mutate(submodel = case_when(submodel == "bulk" ~ "bulk",
                              submodel == "adjusted__1__meta1_pc1" ~ "adjusted_metaPC1",
                              submodel == "adjusted__1__CRCL_0000348" ~ "adjusted_keratinocyte",
                              submodel == "adjusted__1__meta1_pc1_CRCL_0000348" ~ "adjusted_metaPC1_keratinocytes"))

uploadToBQ(gene_cell, bqdataset = "s05_atopic_dermatitis", tableName = "AD_gene_cell_correlations")


# 5. Sample classifications
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