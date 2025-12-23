# Libraries and data loading
# ====================================
devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.assets) # make sure to load version >= 1.0.1
library(tidyverse)

## 0. General
## =======================
# Taken from Aaron's analysis
Metadata = lapply(get_workflow_outputs("wf-fa2180d0c7"), readRDS)

CellMetadata = Metadata[[1]]
uploadToBQ(CellMetadata, bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_CellMetadata")



## 1. Scores
## =========================
# Taken from 7.Scoring_Criteria.R
scores = read_asset("wf-5b6481c61b")
rankings = read_asset("wf-13b22e4fd7")

scores_and_ranks = merge(scores, rankings)
scores_and_ranks = scores_and_ranks[,c("Target.Identifier", "Target.Collection", "Criteria.Identifier", "scaled_score_rounded", "summed_scaled_score_rounded", "rank_scaled_score_rounded")]
colnames(scores_and_ranks)[4:6] = c("scaled_score","final_score","Rank")
uploadToBQ(scores_and_ranks, bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_scores")


## 2. signatures
## ==========================
# Taken from 2.7.Final_Signatures.R
geneMapping = toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
signatures = readRDS(get_workflow_outputs("wf-ae71fd5351"))
signatures = lapply(signatures, function(x){
  x = reshape2::melt(x)
  colnames(x) = c("EntrezID","Signature")
  x$Gene = geneMapping$symbol[match(x$EntrezID, geneMapping$gene_id)]
  return(x)
})
signatures = bind_rows(signatures)

mapping = c("100128338" = "FAM83H-AS1", "10896" = "OCLM","23285" = "KIAA1107")
signatures$Gene[is.na(signatures$Gene)] = mapping[match(signatures$EntrezID[is.na(signatures$Gene)], names(mapping))]

uploadToBQ(signatures, bqdataset = "s05_atopic_dermatitis", tableName = "signatures")


## Target_Cell_Correlations
## ================================
CellAnnotation = read_asset("wf-299d496333")
target_cell = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_Cell")
dz = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "AD_ct_test")

dz = dz %>%
  dplyr::filter(term == "L_vs_HC") %>%
  mutate(Direction_In_Disease = case_when(fdr > 0.05 ~ "Unchanged",
                                          fdr <= 0.05 & estimate > 0 ~ "Up",
                                          fdr <= 0.05 & estimate < 0 ~ "Down")) %>%
  left_join(., CellAnnotation[,c("feature_id","cell_category")])

target_cell = left_join(target_cell, dz[,c("feature_id", "cell_category", "Direction_In_Disease")], by = join_by("Criteria_Identifier" == "feature_id"))
target_cell = target_cell %>%
  dplyr::rename(correlation = "metricValue", FDR = "log10_fdr", Cell = "Criteria_Identifier", target = "Target_Identifier") %>%
  dplyr::select(target, Target_Collection, Cell, cell_category, Direction_In_Disease, correlation, FDR)
uploadToBQ(target_cell, bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_target_cell")


## 