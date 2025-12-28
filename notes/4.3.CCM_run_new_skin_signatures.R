devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)

# In this part, we use hybrid skin signatures that includes mast cells, basophils and ILC2 from different tissues
# in order to mitigate some concerns by the client that we are not capturing type 2 immunity
config = read.csv(get_task_inputs("wf-08a6a0a503","0", files_names_grepl_pattern = "ccm-metadata.csv"))
config$asset_id[which(config$experiment_id == "GSE153007")] <- 'wf-2ac040bba1:0:GSE153007.RDS'
config$asset_id[which(config$experiment_id == "GSE60709")] <- 'wf-923925db68:0:GSE60709.RDS'
config$ccm_annotate = "feature"

sample_ann = readRDS(get_workflow_outputs("wf-558ccafbfe"))
pheno_vars = readRDS(get_workflow_outputs("wf-4c71fb2030"))

ccm_stage_prepare_dataset_collection(config, annotate = list(sample_annotation_table = sample_ann))

IMAGE <- ccm_cyto_cc("master_1.0.1",
                     dataset.cell_contribution.service = "eu.gcr.io/cytoreason/ci-cytoreason.deconvolution-package:SUP_6578_add_skin_v14_v15_v16_2.2.8",
                     cell_specific_differences = "eu.gcr.io/cytoreason/ci-cytoreason.deconvolution-package:SUP_6578_add_skin_v14_v15_v16_2.2.8")

tissue_signature = "skin_v13:SUP-6496-skin_v13"
tissue_signature = "skin_v14:SUP-6578-add-skin_v14-v15-v16"


## First run: integrate the skin signatures into the disease model
## ======================================================================================
ccm_api_generate_ccm(config,
                     qc = FALSE,
                     term_metadata = FALSE, # fails for some reason
                     prepare_data = list(annotate = list(sample_annotation_table = sample_ann)),
                     dataset = list(cell_contribution = list(signature_collection = tissue_signature), 
                                    services = .skip("feature_pca")),
                     dataset2 = list(gene_set_activity = list(collection = c("h", "kegg", "reactome", "btm"))),
                     model = list(services = .skip("gene_set_activity_differences","feature_geneset_correlations")),
                     image = IMAGE,
                     tags = list(tissue="skin_v14", condition="AD", project="evo", analysis="AD_skin_v14"))
# generate_ccm -- Fri Dec 19 10:48:35 2025: wf-e82c5ab874 []
# generate_ccm -- Fri Dec 19 15:06:38 2025: wf-5f05a83bc6 [] - tag with skin_v13
# generate_ccm -- Fri Dec 19 15:06:51 2025: wf-dd8c8843e6 [] - tag with skin_v13:SUP-6496-skin_v13
# generate_ccm -- Fri Dec 19 16:14:48 2025: wf-c82f4f91e9 [] - change deconvolution image
# generate_ccm -- Thu Dec 25 12:57:59 2025: wf-6f00628526 [] - skin_v14
# generate_ccm -- Thu Dec 25 12:57:12 2025: wf-52a7393eb5 [] - skin_v15
# generate_ccm -- Thu Dec 25 12:57:28 2025: wf-6df6f6bdca [] - skin_v16


## Second run: analysis of the custom gene sets with the new cell signatures
## ======================================================================================
signatures = readRDS(get_workflow_outputs("wf-d22894f90b"))[-3]
th2 = list(th2 = readRDS(get_workflow_outputs("wf-410536ebd3")),
           neuroinflammation = readRDS(get_workflow_outputs("wf-4ce41a599a")),
           epidermis = readRDS(get_workflow_outputs("wf-d3b177cd82")))

# reorder
th2$epidermis = c(th2$epidermis, th2$th2[c("lichenification","Terminal_Differentiation_and_Lipids")])
th2$th2 = th2$th2[-which(names(th2$th2) %in% c("lichenification","Terminal_Differentiation_and_Lipids","Th1_Related","Th17_Related","TH22_IL22_Related"))]

allSignatures = c(signatures, th2)
allSignatures = lapply(allSignatures, function(x) lapply(x, unique))
rm(th2, signatures)

gene_set_limits <- setNames(replicate(length(allSignatures), c(1L, Inf), simplify = FALSE), names(allSignatures)) # Define gene set size per signature

ccm_api_run_custom_gene_set_analysis(cytoreason.assets::AssetData("wf-6df6f6bdca") , # when using AssetData it will pull relevant tags automatically, and makes the relationship traceable
                                     custom_gene_set_collections = allSignatures,
                                     prepare_data = list(annotate = list(sample_annotation_table = sample_ann)),
                                     dataset = list(cell_contribution = list(signature_collection = tissue_signature),
                                                    services = .skip("feature_pca")),
                                     dataset2 = list(gene_set_activity = list(collection = c("h", "kegg", "reactome", "btm"))),
                                     model = list(pheno_feature_correlations = list(phenotypic_variables = pheno_vars)),
                                     meta = list(gx_diff = list(collection_size_limits = gene_set_limits)),
                                     submodel = c("bulk","exprs_adjusted__1__1"),
                                     workflow_overrides=list(activeDeadlineSeconds=172800),
                                     image = IMAGE,
                                     tags = list(tissue="skin_v16", condition="AD", project="evo", analysis = "X2_skin_v16"),
                                     data_access = "s05")
# generate_ccm -- Thu Dec 18 15:52:10 2025: wf-a88918a66c [] - failed
# generate_ccm -- Fri Dec 19 18:41:16 2025: wf-2ec72d4434 [] - adjusted failed
# generate_ccm -- Fri Dec 19 18:42:40 2025: wf-c94ee2ad07 [] - no adjusted
# generate_ccm -- Fri Dec 19 18:43:35 2025: wf-d4626ee8f7 [] - with adjusted
# generate_ccm -- Sat Dec 20 19:12:00 2025: wf-496715f457 [] - with cells as pheno feature
# generate_ccm -- Sun Dec 28 09:20:45 2025: wf-882a48484e [] - skin_v16