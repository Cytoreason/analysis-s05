# Prior to this, we need to locally install the patern-benchmark package, from the mayan branch
# when writing this page, we needed to remove "data_access=public" from functions as it causes malfunctions.
# In addition, processing after the functions is done to make sure we can use this in the dashboards
devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.cc.client)
library(tidyverse)

### Prep signatures
### ==========================
allSignatures = readRDS(get_workflow_outputs("wf-54d564beb9"))
  allSignatures = allSignatures$X2
geneMapping = toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
allSignatures = lapply(allSignatures, function(x) geneMapping$symbol[match(x, geneMapping$gene_id)])
pushToCC(allSignatures, tagsToPass = list(list(name="object",value="allSignatures"),
                                          list(name="notes",value="symbol")))
# wf-6590790158
# wf-5600d586c2

### Helper function
### ==========================
add_annotation = function(df, ann) {
  if(ann == "wu") {
    dataset = "Wu 2024 authors annotation"
    dataset_id="https://doi.org/10.1038/s41467-024-44994-w"
    disease = "PSO"
  } else if(ann == "francis") {
    dataset = "Francis 2023 monocle k5 annotated"
    dataset_id="https://doi.org/10.1038/s41467-024-44994-w"
    disease = "PSO"
  }
  df$dataset = dataset
  df$dataset_id = dataset_id
  df$disease = disease
  return(df)
}


### Prep coarse annotation
### ==========================
hier = read.csv("~/data/WuADCells.csv")


### Definitions
### ==========================
allSignatures = readRDS(get_workflow_outputs("wf-5600d586c2"))
config = data.frame(geneset = names(allSignatures),
                    target = rep("MRGPRX2",length(allSignatures)),
                    receptors = rep("MRGPRX2",length(allSignatures)),
                    treatments = rep("Secukinumab;Ixekizumab;Brodalumab;Dupilumab;Fezakinumab;Tocilizumab",length(allSignatures)),
                    diseases = rep("AD",length(allSignatures)),
                    pathways_bioexp = rep("Inflammation;Cellular organization and differentiation",length(allSignatures)))

genelist_name = "X2"
genelist_wfid = "wf-5600d586c2"
scdata_wfid = "wf-654d30a5b7" # Wu
bqdataset = "s05_atopic_dermatitis" ### Please remember to change this, as to not overwrite the defaults
image = "eu.gcr.io/cytoreason/ci-cytoreason.patern.benchmark-package:mayan_latest"

genelist = cytoreason.assets::read_asset(genelist_wfid) %>%
  reshape2::melt(genelist)
colnames(genelist) = c("gene","genelist")

uploadToBQ(table = genelist, bqdataset = bqdataset, tableName = "signatures")


### Test 1 - Enrichment in one cell type
### =========================================================================
devtools::load_all("~/patern-benchmark/R/enrichment_in_one_cell_type.R")

## using Wu
run_function_dist(FUN = enrichment_in_one_cell_type,
                  seurat_wfid = scdata_wfid, 
                  genelist = genelist_wfid, 
                  genelist_title = genelist_name,
                  config = config,
                  full_metadata = "wf-e1d8a4c8cf", 
                  cell_annotation = "local_supercluster", 
                  sample_annotation = "sample_id",
                  pseudobulk_random = "wf-4ff4d52e81",
                  ES_column_random = 'avgES',
                  random_vs_all_randoms = "wf-58fae86a5b",
                  random_cell_vs_all_cells = 'wf-c763c0c3df',
                  uploadTestToBQ = F,
                  BQ_dataset = bqdataset,
                  data_access = "s05", 
                  image = image,
                  tags = list(list(name="analysis",value="enrichment_in_one_cell_type"),
                              list(name="genelist", value=genelist_name),
                              list(name="seurat_object",value=scdata_wfid),
                              list(name="annotation", value="local_supercluster")))

# wf-a7a219220c
# wf-6370038800

results = apply(get_workflow_outputs("wf-6370038800", files_names_grepl_pattern = ".csv"), 1, read.csv)
results = results[!stringr::str_detect(names(results),"wfid_table")]
results = lapply(results, add_annotation, "wu")
results = lapply(results, function(x) {colnames(x)[which(colnames(x) == "geneset")] <- "signature" ; return(x)})
results$statistics_X2_local_supercluster_enrichmentInOneCelltype.csv$random_avg_es = NA
results$statistics_X2_local_supercluster_enrichmentInOneCelltype.csv$coarse_annotation = hier$celltype[match(results$statistics_X2_local_supercluster_enrichmentInOneCelltype.csv$cell, hier$local_supercluster)]

uploadToBQ(table = results$passfail_X2_local_supercluster_enrichmentInOneCelltype.csv, bqdataset = bqdataset, tableName = "enrichmentInOneCelltype_pass_fail")
uploadToBQ(table = results$statistics_X2_local_supercluster_enrichmentInOneCelltype.csv, bqdataset = bqdataset, tableName = "enrichmentInOneCelltype_statistics")
uploadToBQ(table = results$pseudobulk_X2_local_supercluster.csv, bqdataset = bqdataset, tableName = "enrichmentInOneCelltype_pseudobulk")

# ## Using Francis
# run_function_dist(FUN = enrichment_in_one_cell_type,
#                   seurat_wfid = "wf-4d667e568c", 
#                   genelist = genelist_wfid, 
#                   genelist_title = genelist_name,
#                   config = config,
#                   full_metadata = "wf-633678409f", 
#                   cell_annotation = "monocle_subcluster_k_5", 
#                   sample_annotation = "sample_id",
#                   pseudobulk_random = "wf-e4d2d55319",
#                   ES_column_random = 'avgES',
#                   random_vs_all_randoms = "wf-85003efdac",
#                   random_cell_vs_all_cells = 'wf-37adf19e63',
#                   uploadTestToBQ = F,
#                   BQ_dataset = bqdataset,
#                   data_access = "s05", 
#                   image = image,
#                   tags = list(list(name="analysis",value="enrichment_in_one_cell_type"),
#                               list(name="genelist", value=genelist_name),
#                               list(name="seurat_object",value="wf-4d667e568c"),
#                               list(name="annotation", value="monocle_subcluster_k_5")))
# 
# # wf-c3b9e76872


### Test 2 - Coexpression
### =========================================================================
devtools::load_all("~/patern-benchmark/R/coexpression.R")
run_function_dist(FUN = genelist_receptor_coexpression,
                  seurat_wfid = scdata_wfid, 
                  config = config,
                  genelist_title = genelist_name,
                  cell_annotation = "local_supercluster", 
                  gene_cell_vs_all_cells='wf-1d9581bf59',
                  pb_gene = 'wf-1ee1ffbb7b',
                  percentiles =  "wf-db4beb5ecf",
                  test1_wfid = "wf-6370038800",
                  uploadTestToBQ = F,
                  BQ_dataset = bqdataset,
                  image = "eu.gcr.io/cytoreason/ci-cytoreason.patern.benchmark-package:mayan_latest",
                  memory_request = "2Gi",
                  tags = list(list(name="analysis",value="genelist_receptor_coexpression"),
                              list(name="genelist", value=genelist_name),
                              list(name="seurat_object",value=scdata_wfid),
                              list(name="annotation", value="local_supercluster")))
# wf-393cecd739
# wf-93d184a374

results = apply(get_workflow_outputs("wf-393cecd739", files_names_grepl_pattern = ".csv"), 1, read.csv)
results = lapply(results, add_annotation, "wu")
results = lapply(results, function(x) {colnames(x)[which(colnames(x) == "geneset")] <- "signature" ; return(x)})
results$statistics_X2_local_supercluster_coexpression.csv$coarse_annotation = hier$celltype[match(results$statistics_X2_local_supercluster_coexpression.csv$cell, hier$local_supercluster)]

uploadToBQ(table = results$statistics_X2_local_supercluster_coexpression.csv, bqdataset = bqdataset, tableName = "coexpression_statistics")
uploadToBQ(table = results$passfail_X2_local_supercluster_coexpression.csv, bqdataset = bqdataset, tableName = "coexpression_pass_fail")


### Test 3 - Enrichment post vs pre treatment
### =========================================================================
devtools::load_all("~/patern-benchmark/R/treatments.R")

run_function_dist(FUN = enrichment_in_treatments,
                  config = config, 
                  genelist = genelist_wfid, 
                  genelist_title = genelist_name,
                  uploadTestToBQ = F,
                  BQ_dataset = bqdataset,
                  image = "eu.gcr.io/cytoreason/ci-cytoreason.patern.benchmark-package:mayan_latest",
                  memory_request = "2Gi",
                  tags = list(list(name="analysis",value="enrichment_in_treatments"),
                              list(name="genelist", value=genelist_name),
                              list(name="seurat_object",value=scdata_wfid)))
# wf-0fdb40ddc1
# wf-30498594d1

enrichment_results = apply(get_workflow_outputs("wf-30498594d1", files_names_grepl_pattern = ".csv"), 1, read.csv, sep = ",")
enrichment_results = lapply(enrichment_results, function(x) {colnames(x)[which(colnames(x) == "geneset")] <- "signature" ; return(x)})
enrichment_results$treatments_X2_enrichments.csv$a_group_name = 
  enrichment_results$treatments_X2_enrichments.csv$b_group_name = 
  enrichment_results$treatments_X2_enrichments.csv$group_a =
  enrichment_results$treatments_X2_enrichments.csv$group_b = NA

uploadToBQ(table = enrichment_results$treatments_X2_enrichments.csv, bqdataset = bqdataset, tableName = "enrichmentInTreatments_statistics")
uploadToBQ(table = enrichment_results$passfail_X2_treatments.csv, bqdataset = bqdataset, tableName = "enrichmentInTreatments_pass_fail")



### Test 4 - Enrichment in disease vs reference group
### =========================================================================
devtools::load_all("~/patern-benchmark/R/enrichment_in_disease.R")
run_function_dist(FUN = enrichment_in_disease,
                  config = config,
                  model_details = "wf-3a7dddcf49",
                  genelist = genelist_wfid, 
                  genelist_title = genelist_name,
                  uploadTestToBQ = F,
                  BQ_dataset = bqdataset,
                  image = "eu.gcr.io/cytoreason/ci-cytoreason.patern.benchmark-package:mayan_latest",
                  memory_request = "2Gi",
                  tags = list(list(name="analysis",value="enrichment_in_disease"),
                              list(name="genelist", value=genelist_name),
                              list(name="seurat_object",value=scdata_wfid)))
# wf-5cb62acda1
# wf-908be3e391

enrichment_results = apply(get_workflow_outputs("wf-908be3e391", files_names_grepl_pattern = ".csv"), 1, read.csv, sep = ",")
enrichment_results = lapply(enrichment_results, function(x) {colnames(x)[which(colnames(x) == "geneset")] <- "signature" ; return(x)})
enrichment_results$statistics_X2_enrichedInDisease.csv$n_datasets = 
  enrichment_results$statistics_X2_enrichedInDisease.csv$n_control_samples = 
  enrichment_results$statistics_X2_enrichedInDisease.csv$n_disease_samples = NA

uploadToBQ(table = enrichment_results$statistics_X2_enrichedInDisease.csv, bqdataset = bqdataset, tableName = "enrichedInDisease_statistics")
uploadToBQ(table = enrichment_results$passfail_X2_enrichedInDisease.csv, bqdataset = bqdataset, tableName = "enrichedInDisease_pass_fail")


### Test 5 - Gene overlap
### =========================================================================
devtools::load_all("~/patern-benchmark/R/gene_overlap.R")
run_function_dist(FUN = gene_overlap,
                  config = config,
                  genelist = genelist_wfid, 
                  genelist_title = genelist_name,
                  uploadTestToBQ = F,
                  BQ_dataset = bqdataset,
                  image = "eu.gcr.io/cytoreason/ci-cytoreason.patern.benchmark-package:mayan_latest",
                  memory_request = "1Gi",
                  tags = list(list(name="analysis",value="gene_overlap"),
                              list(name="genelist", value=genelist_name),
                              list(name="seurat_object",value=scdata_wfid)))

# wf-76f2324cdf
# wf-438134fa2c

res = get_workflow_outputs("wf-438134fa2c", files_names_grepl_pattern = ".csv", should_download_files = T)

order_file <- res[str_detect(rownames(res), "clustering_order")]
overlap_file <- res[str_detect(rownames(res), "geneOverlap")]

order_t <- read.csv(order_file) %>%
  rename(geneset = genelist)

cluster_t <- read.csv(overlap_file) %>%
  rename(signature_1 = genelist_1, signature_2 = genelist_2) %>%
  mutate(
    target_1 = config$target[match(signature_1, config$geneset)],
    target_2 = config$target[match(signature_2, config$geneset)]
  ) %>%
  left_join(order_t %>% rename(order_1 = order), by = c("signature_1" = "geneset")) %>%
  left_join(order_t %>% rename(order_2 = order), by = c("signature_2" = "geneset")) %>%
  mutate(group_1 = NA, group_2 = NA) %>%
  dplyr::select(-p, -overlap_percent)

uploadToBQ(table = cluster_t, bqdataset = bqdataset, tableName = "geneOverlap_statistics")



### Test 6 - Pathway Overrepresentation
### =========================================================================
devtools::load_all("~/patern-benchmark/R/pathways_enrichment_per_level.R")
run_function_dist(FUN = target_pathway_category_enrichment,
                  config = config,
                  genelist = genelist_wfid,
                  genelist_title = genelist_name,
                  pathways_category_mapping = "wf-0c52b0cfb5",
                  pathways_genes_category_mapping = "wf-19f8069473",
                  pathways_category_level = "level0",
                  genes_background = 'wf-6f21780a5d',
                  dont_rerun_pathway_enrichment_wfid = "wf-078cfb6e86", # This was added because remote_parallel_hg failed. it's the wfid of pathways_df
                  pathways_collections = c("reactome", "kegg"),
                  n_parallel_workers = 20,
                  uploadTestToBQ = F,
                  BQ_dataset = bqdataset,
                  image = "eu.gcr.io/cytoreason/ci-cytoreason.patern.benchmark-package:mayan_latest",
                  memory_request = "5Gi",
                  tags = list(list(name="analysis",value="pathaways enrichments"),
                              list(name="genelist", value=genelist_name),
                              list(name="pipline", value="signatures validations")))

# wf-35c3fefddf
# wf-8a9e6e3808

res = apply(get_workflow_outputs("wf-8a9e6e3808", files_names_grepl_pattern = ".csv"), 1, read.csv)

res <- lapply(res, function(df) {
  df <- df %>%
    rename(signature = if ("geneset" %in% names(df)) "geneset" else "listName")
  
  if ("category" %in% names(df)) {
    df <- df %>%
      mutate(category = str_to_sentence(category),
             target = config[match(signature, config$geneset), "target"])
  }
  
  df
})

uploadToBQ(table = res$pathwayExploration_X2_pathwaysEnrichmentPerCategory.csv, bqdataset = bqdataset, tableName = "pathwaysEnrichment_pathwayExploration")
uploadToBQ(table = res$categoriesStats_X2_pathwaysEnrichmentPerCategory.csv, bqdataset = bqdataset, tableName = "pathwaysEnrichment_categoriesStats")
uploadToBQ(table = res$significantCategoriesDetails_X2_pathwaysEnrichmentPerCategory.csv, bqdataset = bqdataset, tableName = "pathwaysEnrichment_significantCategoriesDetails")
uploadToBQ(table = res$categoriesPassFail_X2_pathwaysEnrichmentPerCategory.csv, bqdataset = bqdataset, tableName = "pathwaysEnrichment_categoriesPassFail")



### Landing Table
### =========================================================================
devtools::load_all("~/patern-benchmark/R/landing_page_tables.R")
run_function_dist(FUN = get_landing_page_table,
                  test1_summary_table = "wf-6370038800",
                  test2_summary_table = "wf-93d184a374",
                  test3_summary_table = "wf-30498594d1",
                  test4_summary_table = "wf-908be3e391",
                  test6_summary_table = "wf-8a9e6e3808",
                  genelist_title = genelist_name,
                  config=config,
                  uploadTestToBQ = F,
                  BQ_dataset = bqdataset,
                  image = "eu.gcr.io/cytoreason/ci-cytoreason.patern.benchmark-package:mayan_latest",
                  memory_request = "1Gi",
                  tags = list(list(name="analysis",value="join all pass-fail tables"),
                              list(name="genelist", value=genelist_name),
                              list(name="pipline", value="signatures validations")))
# wf-f91d9f2e0a

res = read.csv(get_workflow_outputs("wf-f91d9f2e0a", files_names_grepl_pattern = ".csv"))%>%
  dplyr::rename(signature = geneset)%>%
  dplyr::mutate(target = config[match(signature, config$geneset),]$target)

uploadToBQ(table = res, bqdataset = bqdataset, tableName = "landingPageTable")


<<<<<<< Updated upstream
=======
## OFF THE BOOKS HYPERGEOMETRIC
## ==================================
signatures = readRDS(get_workflow_outputs("wf-54d564beb9"))
signatures = signatures$X2[c("x2_inhibition_early_50","x2_inhibition_early_50_archs_refined","SP_inhibition_late_50","SP_inhibition_late_50_archs_refined")]

gxdiff = readRDS(get_workflow_outputs("wf-a121700ef7"))
gxdiff = gxdiff %>%
  dplyr::filter(comparison %in% c("SP_inhibited_vs_uninhibited_24hr","x2_inhibited_vs_activated_4hr")) %>%
  dplyr::filter(estimate < 0)

hg = lapply(signatures, function(x) cytoreason.gx::gx_gsa(x = as.character(x),
                                                          method = "hypergeometric",
                                                          background = gxdiff$feature_id)$fit)
hg_filtered = lapply(hg, function(x) {
  x = lapply(names(x), function(y) {
    x[[y]]$collection = y
    return(x[[y]])
  })
  x = bind_rows(x)
  x = x[-which(x$pvalue == 1),]
  x = x[-which(x$collection %in% c("c2.cp.reactome","c2.cp.kegg")),]
})

lapply(names(hg_filtered), function(sig){
  hg_filtered[[sig]]$signature <<- sig
})

hg_final = do.call(rbind, hg_filtered)
uploadToBQ(table = hg_final, bqdataset = bqdataset, tableName = "X2Signatures_hypergeometric")


# Visualizations
# =======================
# Landing table
>>>>>>> Stashed changes
library(ComplexHeatmap)
res = res[which(res$collection %in% c("","reactom")),]
res = res[str_detect(res$signature,"_ep"),]
res$result = ifelse(res$result == "Fail", NA, 1)
res = reshape2::dcast(res, signature ~ test_name, value.var = "result")
res = column_to_rownames(res, var = "signature") %>% as.matrix
res = res[rowSums(is.na(res)) != ncol(res),]
res[is.na(res)] <- -1
Heatmap(res, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6), show_row_dend = F, show_column_dend = F)
