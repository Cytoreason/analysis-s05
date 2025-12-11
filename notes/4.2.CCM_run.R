devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)

# 1. Inputs
# -----------------------------
ccm = as_ccm_fit("wf-08a6a0a503")
signatures = readRDS(get_workflow_outputs("wf-d22894f90b"))
## removing the extra data set isn't feasible within the constraints of the project - need to re-run the disease model, causing compatibility issues


# 2. Molecular and clinical scores
# ------------------------------------
MS = readRDS(get_workflow_outputs("wf-db0c88cca6"))
# calculated by the data science group https://cytoreason.atlassian.net/wiki/spaces/DSG/pages/4154261510/New+Molecular+Score+-+Documentation#Molecular-score-results-on-different-disease-models

MS = pData(MS$ms_results$combat_ms) # MS based on L vs HC
MS = MS %>%
  dplyr::select(condition, sample_classification, geo_accession,molecular_score) %>%
  mutate(condition = str_replace(condition, "control", "healthy")) %>%
  mutate(condition = factor(condition, ordered = T, levels = c("atopic dermatitis","contact dermatitis","psoriasis","healthy")))


# 3. Unify SCORAD/EASI scores and change name of cell meta PCs
# ----------------------------------------------------------------
for (x in names(ccm$datasets)){
  pDataset = pData(ccm$datasets[[x]])
  if(x %in% c("GSE121212__GPL16791","GSE137430__rnaseq","GSE141570__GPL20301","GSE147424__GPL16791","GSE157194__GPL21290")) { 
    samplename = "sample_id"
  } else {
    samplename = "geo_accession"
  }
  
  # integrating molecular score
  pDataset$MolecularScore <- MS$molecular_score[match(pDataset$geo_accession, MS$geo_accession)]

  if(x %in% c("GSE130588__GPL570","GSE27887__GPL570","GSE58558__GPL570")) { # there we have SCORAD
    idx = which(str_detect(colnames(pDataset),"SCORAD|Scorad"))[1]
    pDataset$SCORAD = as.numeric(pDataset[,idx])
    rm(idx)
  }
  if(x %in% c("GSE130588__GPL570","GSE59294__GPL570")) { # there we have EASI
    idx = which(str_detect(colnames(pDataset),"EASI"))[1]
    pDataset$EASI = as.numeric(pDataset[,idx])
    rm(idx)
  }
  
  pData(ccm$datasets[[x]]) <- pDataset
  rm(pDataset)
}


## 4. Generating a model with keratinocytes as an adjustment model
## and SCOARD/EASI/cell meta PCs standardized
## --------------------------------------------------------------------
ccm_api_generate_ccm(ccm,
                     qc = FALSE,
                     term_metadata = FALSE, # fails for some reason
                     adjustment_models = list("1" = list(c("meta1_pc1"), c("CRCL_0000348"),c("meta1_pc1","CRCL_0000348"))),
                     model = .skip("cell_specific_differences"),
                     image = "master_1.0.1",
                     tags = list(list(tissue="skin"),
                                 list(condition="AD"),
                                 list(project="evo"),
                                 list(name="adjustment",value="keratinocyte_only")))
# generate_ccm -- Tue Dec  9 07:25:54 2025: wf-832ab799be []

new_ccm = as_ccm_fit("wf-832ab799be")


# 5. Integrate with pathways meta pcs
# ------------------------------------------
pathwayPCs = readRDS(get_workflow_outputs("wf-e8f514f3f1"))
pathwayPCs = pathwayPCs %>%
  mutate(term = ifelse(term == "AD", "DZ_vs_HC", term)) %>%
  mutate(collection = ifelse(collection == "All","",collection))%>%
  mutate(suffix = paste0(term,"_",submodel,"_",collection)) %>%
  dplyr::select(-c(term:collection))

pathwayPCs = reshape2::melt(pathwayPCs, id.vars = c(1,2,6), measured.vals = 3:5, variable.name = "prefix", value.name = "sampleScore")
pathwayPCs = pathwayPCs %>%
  mutate(label = paste0(prefix,"_",suffix)) %>%
  dplyr::select(-c(suffix:prefix))

pathwayPCs = reshape2::dcast(pathwayPCs, dataset + sample_id ~ label)


for (x in names(ccm$datasets)){
  pDataset = pData(ccm$datasets[[x]])
  if(x %in% c("GSE121212__GPL16791","GSE137430__rnaseq","GSE141570__GPL20301","GSE147424__GPL16791","GSE157194__GPL21290")) { 
    samplename = "sample_id"
  } else {
      samplename = "geo_accession"
    }
  
  # integrating molecular score
  pDataset$MolecularScore <- MS$molecular_score[match(pDataset$geo_accession, MS$geo_accession)]
  
  # cell meta PCs
  colnames(pDataset) <- gsub("meta1_", "cell_meta_", colnames(pDataset))
  pDataset$cell_meta_pc1 <- (-1) * pDataset$cell_meta_pc1 # flipped the sign of the cell meta-PC1 (because HC>DZ)
  
  # pathway meta PCs
  pDataset = merge(pDataset, pathwayPCs[,-1], by.x = samplename, by.y = "sample_id")
  
  if(x %in% c("GSE130588__GPL570","GSE27887__GPL570","GSE58558__GPL570")) { # there we have SCORAD
    idx = which(str_detect(colnames(pDataset),"SCORAD|Scorad"))[1]
    pDataset$SCORAD = as.numeric(pDataset[,idx])
    rm(idx)
  }
  if(x %in% c("GSE130588__GPL570","GSE59294__GPL570")) { # there we have EASI
    idx = which(str_detect(colnames(pDataset),"EASI"))[1]
    pDataset$EASI = as.numeric(pDataset[,idx])
    rm(idx)
  }
  
  pData(ccm$datasets[[x]]) <- pDataset
  rm(pDataset)
}



## General
## ------------------
# Set gene-set size limit for ssgsea and GSEA for our custom gene lists:
gene_set_limits <- setNames(replicate(length(signatures), c(1L, Inf), simplify = FALSE), names(signatures)) # Define gene set size per signature
pheno_vars = c("MolecularScore","SCORAD","EASI", "cell_meta_pc1","cell_meta_pc2","cell_meta_pc3",
               "pathway_meta_pc1","pathway_meta_pc2","pathway_meta_pc3",
               "adj_pathway_meta_pc1","adj_pathway_meta_pc2","adj_pathway_meta_pc3")


## Running the CCM
## ----------------------
ccm_api_run_custom_gene_set_analysis(ccm, # when using AssetData it will pull relevant tags automatically, and makes the relationship traceable
                                     custom_gene_set_collections = signatures,
                                     model = list(gx_gsa = list(collection_size_limits = gene_set_limits),
                                                  pheno_feature_correlations = list(phenotypic_variables = pheno_vars)),
                                     meta = list(gx_diff = list(collection_size_limits = gene_set_limits)),
                                     submodel = c("bulk", "adjusted__1__1"),
                                     image = "master_1.0.1",
                                     tags = list(list(tissue="skin"),
                                                 list(condition="AD"),
                                                 list(project="evo"),
                                                 list(name="analysis",value="X2_V6")),
                                     data_access = "s05")
# generate_ccm -- Mon Nov 24 19:05:53 2025: wf-ef57aebb52 [] - without adj pathway meta
# generate_ccm -- Tue Nov 25 09:47:02 2025: wf-8e948630d7 []


## Visualization
## ==============================
ggplot(MS, aes(y = molecular_score, x = sample_classification, fill = sample_classification)) +
  geom_violin(draw_quantiles = T) +
  ggpubr::stat_compare_means(method = "wilcox", comparisons = list(c("Lesion","Non Lesion"),c("Non Lesion","Normal"),c("Lesion","Normal")), paired = F) +
  scale_fill_manual(values = c("#473472","#53629E","#87BAC3"))+
  theme_minimal()+
  theme(legend.position = "none")+
  labs(x = NULL, title = "New molecular score in AD model")
ggsave("~/analysis-s05/figures/AD Model/molecularScore_sampleClassification.png", bg = "white")

ggplot(MS, aes(y = molecular_score, x = condition, fill = condition)) +
  geom_violin(draw_quantiles = T) +
  ggpubr::stat_compare_means(method = "wilcox", comparisons = list(c("atopic dermatitis","contact dermatitis"),
                                                                   c("contact dermatitis","psoriasis"),
                                                                   c("psoriasis","healthy"),
                                                                   c("atopic dermatitis","psoriasis"),
                                                                   c("contact dermatitis","healthy"),
                                                                   c("atopic dermatitis","healthy")), paired = F) +
  scale_fill_manual(values = c("#473472","#53629E","#87BAC3","#D6F4ED"))+
  theme_minimal()+
  labs(x = NULL, title = "New molecular score in AD model") +
  theme(legend.position = "none")
ggsave("~/analysis-s05/figures/AD Model/molecularScore_condition.png", bg = "white")
