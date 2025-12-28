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

# for pathway, integrate new pathways from CCM
PathwayMetadata = Metadata[[1]]
dz = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "AD_gx_gsa")
dz = dz[-which(dz$pathway %in% PathwayMetadata$pathway),]
dz = dz[which(dz$collection %in% c("btm","neuroinflammation", "th2","Ligands","Positives","Neuronal")),]
dz = dz %>%
  mutate(DZ.change = case_when(FDR > 0.05 ~ "Unchanged",
                               FDR <= 0.05 & NES > 0 ~ "Up",
                               FDR <= 0.05 & NES < 0 ~ "Down")) %>%
  dplyr::rename(neglog10FDR = log10_fdr) %>%
  mutate(adjusted = ifelse(submodel == "bulk", "no", "yes")) %>%
  mutate(Key.Pathway = NA, Pathway.level = NA, PipelineWF = "wf-3e419ff83b", Disease = "AD") %>%
  dplyr::select(pathway, collection,Pathway.level,NES,neglog10FDR,DZ.change,adjusted,Key.Pathway,Disease,PipelineWF)

PathwayMetadata = rbind(PathwayMetadata, dz)
uploadToBQ(PathwayMetadata, bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_PathwayMetadata")


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


## Target_Pathway_Correlations
## =====================================
target_pathway = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_Pathway", pageSize = 20000)

target_pathway = target_pathway %>%
  dplyr::rename(correlation = metricValue, pathway = Criteria_Identifier, target = Target_Identifier,
                submodel = Type, target_collection = Target_Collection, collection = Criteria_Collection) %>%
  dplyr::select(target, submodel, fdr, pathway, target_collection, correlation, collection) %>%
  mutate(fdr = -log10(fdr)) %>%
  mutate(submodel = ifelse(submodel == "bulk", "bulk", "adjusted"))

uploadToBQ(target_pathway, bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_target_pathway")


## Pathway volcano
## ==========================
dz = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "AD_gx_gsa")
dz = dz %>%
  dplyr::filter(!collection %in% c("c2.cgp","c2.cp","c2.cp.biocarta","c2.cp.kegg", "c2.cp.reactome","c3.tft","c7","epidermis","Neuronal","Mast","Itch")) %>%
  mutate(collection = ifelse(pathway %in% c("Terminal_Differentiation_and_Lipids", "lichenification"), "epidermis",collection))

treatments = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "treatment_PrePostDupi")
treatments = treatments %>%
  dplyr::filter(!collection %in% c("c2.cgp","c2.cp","c2.cp.biocarta","c2.cp.kegg", "c2.cp.reactome","c3.tft","c7","epidermis","Neuronal","Mast","Itch")) %>%
  mutate(collection = ifelse(pathway %in% c("Terminal_Differentiation_and_Lipids", "lichenification"), "epidermis",collection)) %>%
  dplyr::filter(term %in% c("W4_vs_W0:NR_L", "W4_vs_W0:R_L", "W16_vs_W0:NR_L", "W16_vs_W0:R_L", "W4_vs_W0:DupilumabL", "W16_vs_W0:DupilumabL"))

mapping = c("W4_vs_W0:NR_L" = "W4_vs_W0:Dupilumab:Non-Responders:Lesional", 
            "W4_vs_W0:R_L" = "W4_vs_W0:Dupilumab:Responders:Lesional", 
            "W16_vs_W0:NR_L" = "W16_vs_W0:Dupilumab:Non-Responders:Lesional", 
            "W16_vs_W0:R_L" = "W16_vs_W0:Dupilumab:Responders:Lesional", 
            "W4_vs_W0:DupilumabL" = "W4_vs_W0:Dupilumab:Lesional", 
            "W16_vs_W0:DupilumabL" = "W16_vs_W0:Dupilumab:Lesional")
treatments$term = mapping[match(treatments$term, names(mapping))]

joined = bind_rows(dz, treatments)
joined = joined %>%
  dplyr::select(term, pathway, collection, submodel, NES, log10_fdr)

uploadToBQ(joined, bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_volcano_DiseaseAndDupilumab")
