devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)

## 1. Inputs
## ================================
ccm = as_ccm_fit("wf-08a6a0a503")
## removing the extra data set isn't feasible within the constraints of the project - need to re-run the disease model, causing compatibility issues

# preparing the sample annotation table we will pass, and integrate MS/CS/pheno into it
sample_ann = lapply(ccm$datasets, function(d){
  p = pData(assayDataExpression(d))
  if("week:ch1" %in% colnames(p)) {
    p[,"week:ch1"] = as.character(p[,"week:ch1"])
  }
  if("Week" %in% colnames(p)) {
    p[,"Week"] = as.character(p[,"Week"])
  }
  if("batch_date:ch1" %in% colnames(p)) {
    p[,"batch_date:ch1"] = as.character(p[,"batch_date:ch1"])
  }
  if("Batch_date" %in% colnames(p)) {
    p[,"Batch_date"] = as.character(p[,"Batch_date"])
  }
  if("channel_count" %in% colnames(p)) {
    p[,"channel_count"] = as.character(p[,"channel_count"])
  }
  if("taxid_ch1" %in% colnames(p)) {
    p[,"taxid_ch1"] = as.character(p[,"taxid_ch1"])
  }
  if("data_row_count" %in% colnames(p)) {
    p[,"data_row_count"] = as.character(p[,"data_row_count"])
  }
  if("patient_id:ch1" %in% colnames(p)) {
    p[,"patient_id:ch1"] = as.character(p[,"patient_id:ch1"])
  }
  if("Patient_id" %in% colnames(p)) {
    p[,"Patient_id"] = as.character(p[,"Patient_id"])
  }
  if("description.1" %in% colnames(p)) {
    p[,"description.1"] = as.character(p[,"description.1"])
  }
  return(p)
})
sample_ann = bind_rows(sample_ann)
addedAnnotation = sample_ann[,c("sample_id","geo_accession"),drop=F]

## 2. Integrate molecular score
## ================================
MS = readRDS(get_workflow_outputs("wf-db0c88cca6"))
# calculated by the data science group https://cytoreason.atlassian.net/wiki/spaces/DSG/pages/4154261510/New+Molecular+Score+-+Documentation#Molecular-score-results-on-different-disease-models

MS = pData(MS$ms_results$combat_ms) # MS based on L vs HC
MS = MS %>%
  dplyr::select(condition, sample_classification, geo_accession,molecular_score) %>%
  mutate(condition = str_replace(condition, "control", "healthy")) %>%
  mutate(condition = factor(condition, ordered = T, levels = c("atopic dermatitis","contact dermatitis","psoriasis","healthy")))

addedAnnotation$MolecularScore = MS$molecular_score[match(addedAnnotation$geo_accession, MS$geo_accession)]
addedAnnotation = addedAnnotation[,-2]

pheno_vars = "MolecularScore"
               
## 3. Unify SCORAD/EASI scores
## ================================
sample_ann <- sample_ann %>% 
  dplyr::select(-c("SCORAD","Scorad","easi:ch1")) %>% # redundant columns
  dplyr::rename(SCORAD = "scorad:ch1")

pheno_vars = c(pheno_vars, "SCORAD","EASI")

## 4. Integrate with pathways meta pcs
## ======================================
# 4.1. All meta PCs
# --------------------
pathwayPCs = readRDS(get_workflow_outputs("wf-39ae62c8d2"))
pathwayPCs = pathwayPCs %>%
  dplyr::filter(collection == "all" & term %in% c("AD","L_vs_HC")) %>%
  mutate(term = ifelse(term == "AD", "DZ_vs_HC", term)) %>%
  mutate(submodel = ifelse(submodel != "bulk", "adjusted", "bulk")) %>%
  mutate(suffix = paste0(term,"_",submodel,"_",collection)) %>%
  dplyr::select(-c(pathway_meta_pc3:collection))

pathwayPCs = reshape2::melt(pathwayPCs, id.vars = c(1,2,5), measured.vals = 3:4, variable.name = "prefix", value.name = "sampleScore")
pathwayPCs = pathwayPCs %>%
  mutate(label = paste0(prefix,"_",suffix)) %>%
  dplyr::select(-c(suffix:prefix))

pathwayPCs = reshape2::dcast(pathwayPCs, sample_id ~ label, value.var = "sampleScore")
colnames(pathwayPCs) = str_replace(colnames(pathwayPCs),"_all","_key")

addedAnnotation = merge(addedAnnotation, pathwayPCs, by = "sample_id", all = T)
pheno_vars = c(pheno_vars, colnames(pathwayPCs[,-1]))
               
# 4.2. Th2 meta PCs
# ----------------------
pathwayPCs_th2_bulk = readRDS(get_workflow_outputs("wf-e6c43f2c55"))
pathwayPCs_th2_bulk = lapply(names(pathwayPCs_th2_bulk), function(x) {
  y = pathwayPCs_th2_bulk[[x]]
  colnames(y) = c("sample_id",paste0("pathway_meta_pc",1:3))
  colnames(y)[-1] = paste0(colnames(y)[-1],"_",x,"_bulk")
  return(y)
}) %>% purrr::reduce(., full_join, by = "sample_id")

pathwayPCs_th2_adj = readRDS(get_workflow_outputs("wf-40fbd742a3"))
pathwayPCs_th2_adj = lapply(names(pathwayPCs_th2_adj), function(x) {
  y = pathwayPCs_th2_adj[[x]]
  colnames(y) = c("sample_id",paste0("pathway_meta_pc",1:3))
  colnames(y)[-1] = paste0(colnames(y)[-1],"_",x,"_adj")
  return(y)
}) %>% purrr::reduce(., full_join, by = "sample_id")

addedAnnotation = merge(addedAnnotation, pathwayPCs_th2_bulk, by = "sample_id", all = T)
addedAnnotation = merge(addedAnnotation, pathwayPCs_th2_adj, by = "sample_id", all = T)

pheno_vars = c(pheno_vars, colnames(pathwayPCs_th2_bulk[,-1]), colnames(pathwayPCs_th2_adj[,-1]))


# 4.3. Key + Th2 meta PCs
# --------------------------
pathwayPCs = readRDS(get_workflow_outputs("wf-3fce829a35"))
pathwayPCs = pathwayPCs %>%
  mutate(term = ifelse(term == "AD", "DZ_vs_HC", term)) %>%
  mutate(submodel = ifelse(submodel != "bulk", "adjusted", "bulk")) %>%
  mutate(suffix = paste0(term,"_",submodel,"_",collection)) %>%
  dplyr::select(-c(pathway_meta_pc3:collection))

pathwayPCs = reshape2::melt(pathwayPCs, id.vars = c(1,2,5), measured.vals = 3:4, variable.name = "prefix", value.name = "sampleScore")
pathwayPCs = pathwayPCs %>%
  mutate(label = paste0(prefix,"_",suffix)) %>%
  dplyr::select(-c(suffix:prefix))

pathwayPCs = reshape2::dcast(pathwayPCs, sample_id ~ label, value.var = "sampleScore")

addedAnnotation = merge(addedAnnotation, pathwayPCs, by = "sample_id", all = T)
pheno_vars = c(pheno_vars, colnames(pathwayPCs[,-1]))


# 4.4. Un-redundant pathways (including Th2)
# ------------------------------------------------
pathwayPCs = readRDS(get_workflow_outputs("wf-64e555ed8f"))
pathwayPCs = pathwayPCs %>%
  mutate(term = ifelse(term == "AD", "DZ_vs_HC", term)) %>%
  mutate(submodel = ifelse(submodel != "bulk", "adjusted", "bulk")) %>%
  mutate(suffix = paste0(term,"_",submodel,"_",collection)) %>%
  dplyr::select(-c(term:collection))

pathwayPCs = reshape2::melt(pathwayPCs, id.vars = c(1,2,5), measured.vals = 3:4, variable.name = "prefix", value.name = "sampleScore")
pathwayPCs = pathwayPCs %>%
  mutate(label = paste0(prefix,"_",suffix)) %>%
  dplyr::select(-c(suffix:prefix))

pathwayPCs = reshape2::dcast(pathwayPCs, sample_id ~ label, value.var = "sampleScore")

addedAnnotation = merge(addedAnnotation, pathwayPCs, by = "sample_id", all = T)
pheno_vars = c(pheno_vars, colnames(pathwayPCs[,-1]))


# 4.5. Treatment PCAs
# ------------------------------------------------
# Treatment
pathwayPCs = readRDS(get_workflow_outputs("wf-5c0d927f2c"))
pathwayPCs = pathwayPCs %>%
  mutate(term = ifelse(term == "AD", "DZ_vs_HC", term)) %>%
  mutate(submodel = ifelse(submodel != "bulk", "adjusted", "bulk")) %>%
  mutate(suffix = paste0(term,"_",submodel,"_",collection)) %>%
  dplyr::select(-c(pathway_meta_pc3:collection))

pathwayPCs = reshape2::melt(pathwayPCs, id.vars = c(1,2,5), measured.vals = 3:4, variable.name = "prefix", value.name = "sampleScore")
pathwayPCs = pathwayPCs %>%
  mutate(label = paste0(prefix,"_",suffix)) %>%
  dplyr::select(-c(suffix:prefix))

pathwayPCs = reshape2::dcast(pathwayPCs, sample_id ~ label, value.var = "sampleScore")

addedAnnotation = merge(addedAnnotation, pathwayPCs, by = "sample_id", all = T)
pheno_vars = c(pheno_vars, colnames(pathwayPCs[,-1]))


# RNR
pathwayPCs = readRDS(get_workflow_outputs("wf-0109132b74"))
pathwayPCs = pathwayPCs %>%
  mutate(term = ifelse(term == "AD", "DZ_vs_HC", term)) %>%
  mutate(submodel = ifelse(submodel != "bulk", "adjusted", "bulk")) %>%
  mutate(suffix = paste0(term,"_",submodel,"_",collection)) %>%
  dplyr::select(-c(pathway_meta_pc3:collection))

pathwayPCs = reshape2::melt(pathwayPCs, id.vars = c(1,2,5), measured.vals = 3:4, variable.name = "prefix", value.name = "sampleScore")
pathwayPCs = pathwayPCs %>%
  mutate(label = paste0(prefix,"_",suffix)) %>%
  dplyr::select(-c(suffix:prefix))

pathwayPCs = reshape2::dcast(pathwayPCs, sample_id ~ label, value.var = "sampleScore")

addedAnnotation = merge(addedAnnotation, pathwayPCs, by = "sample_id", all = T)
pheno_vars = c(pheno_vars, colnames(pathwayPCs[,-1]))


## 6. Save objects
## =======================
all(pheno_vars[-c(2:3)] %in% colnames(addedAnnotation[,-1])) # TRUE
any(is.na(addedAnnotation)) # FALSE
sample_ann = merge(sample_ann, addedAnnotation, by = "sample_id", all = T)

sample_ann$meta1_pc1 = (-1) * sample_ann$meta1_pc1 # flipped the sign of the cell meta-PC1 (because HC>DZ)
pushToCC(sample_ann, tagsToPass = list(list(name="object", value="AD_sample_annotation")))
# wf-558ccafbfe
# wf-adbac4aaf3

pheno_vars = c(pheno_vars, "meta1_pc1","meta1_pc2","meta1_pc3")
pushToCC(pheno_vars, tagsToPass = list(list(name="object", value="AD_pheno_vars")))
# wf-4c71fb2030
# wf-d83dfe1bc9