# Libraries and data loading
# ====================================
devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.assets) # make sure to load version >= 1.0.1
library(tidyverse)

keep_signatures = openxlsx::read.xlsx("~/analysis-s05/data/Final signatures to be ranked or viewed in the dashboards.xlsx")
signatureMapping = readRDS(get_workflow_outputs("wf-aa75ed069b"))


## 0. General
## =======================
# Integrate new mast cells to Aaron's tables
Metadata = lapply(get_workflow_outputs("wf-fa2180d0c7"), readRDS)
CellAnnotation = read_asset("wf-299d496333")
dz = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "skinV16_ct_test")

CellMetadata = Metadata[[1]]
newCells = unique(dz$feature_id[-which(dz$feature_id %in% CellMetadata$feature_id)])
dz = dz %>%
  dplyr::filter(feature_id %in% newCells &
                term == "DZ_vs_HC") %>%
  mutate(DZ.change = case_when(fdr > 0.05 ~ "NC",
                               fdr <= 0.05 & estimate > 0 ~ "Up",
                               fdr <= 0.05 & estimate < 0 ~ "Down")) %>%
  mutate(neglog10FDR = -log10(fdr)) %>%
  mutate(Key.Cell = NA, Pathway.level = NA, PipelineWF = "wf-882a48484e", Disease = "AD") %>%
  dplyr::select(Disease, feature_id, effect_size,neglog10FDR,DZ.change,Key.Cell,PipelineWF)
dz = left_join(dz, CellAnnotation[,c("feature_id","cell_category","cell_subcategory")])
dz$cell_category = "Immune"
dz$cell_subcategory[2:3] = c("Innate lymphoid cell","Granulocyte")
CellMetadata = bind_rows(CellMetadata, dz)

uploadToBQ(CellMetadata, bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_CellMetadata")

# for pathways, integrate new pathways from CCM
PathwayMetadata = Metadata[[1]]
  PathwayMetadata$ID = paste0(PathwayMetadata$collection, ":", PathwayMetadata$pathway)
  PathwayMetadata$DZ.change = str_replace(PathwayMetadata$DZ.change, "NC","Unchanged")

dz = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "AD_gx_gsa")
dz = dz[which(dz$collection %in% c("btm","neuroinflammation", "th2","Ligands","Positives","epidermis","reactome","kegg","h")),]
dz = dz %>%
  dplyr::filter(term == "L_vs_HC") %>%
  mutate(DZ.change = case_when(FDR > 0.05 ~ "Unchanged",
                               FDR <= 0.05 & NES > 0 ~ "Up",
                               FDR <= 0.05 & NES < 0 ~ "Down"))

dz$pathway[which(dz$collection %in% c("neuroinflammation", "th2","epidermis"))] = str_replace(dz$pathway[which(dz$collection %in% c("neuroinflammation", "th2","epidermis"))], ":", "__")
dz$collection[which(dz$pathway %in% c("lichenification","Terminal_Differentiation_and_Lipids"))] = "epidermis"
dz$ID = paste0(dz$collection, ":", dz$pathway)

pushToCC(dz, tagsToPass = list(list(name="object",value="pathway_meta")))
# wf-0e90fe5faa

# PathwayMetadata = rbind(PathwayMetadata, dz)
# uploadToBQ(PathwayMetadata, bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_PathwayMetadata")


## 1. Scores
## =========================
# Taken from 7.Scoring_Criteria.R

# For disease criteria
scores = read_asset("wf-86663471a4")
rankings = read_asset("wf-22cecea41a")

scores_and_ranks = merge(scores, rankings)
scores_and_ranks = scores_and_ranks[,c("Target.Identifier", "Target.Collection", "Criteria.Identifier", "scaled_score_rounded", "summed_scaled_score_rounded", "rank_scaled_score_rounded")]
colnames(scores_and_ranks)[4:6] = c("scaled_score","final_score","Rank")
uploadToBQ(scores_and_ranks, bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_scores")


# For Dupilumab complementarity criteria
scores = read_asset("wf-070698da21")
rankings = read_asset("wf-52b9d9053a")

scores_and_ranks = merge(scores, rankings)
scores_and_ranks = scores_and_ranks[,c("Target.Identifier", "Target.Collection", "Criteria.Identifier", "scaled_score_rounded", "summed_scaled_score_rounded", "rank_scaled_score_rounded")]
colnames(scores_and_ranks)[4:6] = c("scaled_score","final_score","Rank")
uploadToBQ(scores_and_ranks, bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_scores_complementarity")


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
target_cell = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "skinV16_geneset_cell_corr")
dz = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "skinV16_ct_test")

dz = dz %>%
  dplyr::filter(term == "L_vs_HC") %>%
  mutate(Direction_In_Disease = case_when(fdr > 0.05 ~ "Unchanged",
                                          fdr <= 0.05 & estimate > 0 ~ "Up",
                                          fdr <= 0.05 & estimate < 0 ~ "Down")) %>%
  left_join(., CellAnnotation[,c("feature_id","cell_category")])

target_cell = left_join(target_cell, dz[,c("feature_id", "cell_category", "Direction_In_Disease")], by = join_by("cell" == "feature_id"))
target_cell = target_cell %>%
  dplyr::filter(submodel == "bulk") %>%
  dplyr::rename(FDR = "log10_fdr", Cell = "cell", target = "pathway", Target_Collection = collection) %>%
  dplyr::select(target, Target_Collection, Cell, cell_category, Direction_In_Disease, correlation, FDR)

target_cell = target_cell[which(target_cell$target %in% keep_signatures$Target[which(keep_signatures$Exploration == "Yes")]),]

uploadToBQ(target_cell, bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_target_cell")


## Target_Pathway_Correlations
## =====================================
target_pathway = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_Pathway", pageSize = 20000)
dz = readRDS(get_workflow_outputs("wf-0e90fe5faa"))

target_pathway = target_pathway %>%
  dplyr::rename(correlation = metricValue, pathway = Criteria_Identifier, target = Target_Identifier,
                submodel = Type, target_collection = Target_Collection, collection = Criteria_Collection) %>%
  dplyr::select(target, submodel, fdr, pathway, target_collection, correlation, collection) %>%
  mutate(fdr = -log10(fdr)) %>%
  mutate(submodel = ifelse(submodel == "bulk", "bulk", "adjusted"))

rmRows = which(target_pathway$target %in% keep_signatures$Target[which(keep_signatures$Exploration == "No")])
target_pathway = target_pathway[-rmRows,]

# organizing collections
# for targets we keep the original (e.g. Itch), for pathways we change and remove redundancies
target_pathway = target_pathway[-which(target_pathway$collection %in% c("Neuronal","Mast","X2")),] # all in neuroinflammation
idx = which(target_pathway$collection %in% c("Itch"))
target_pathway$collection[idx] = "neuroinflammation"
target_pathway$pathway = str_replace(target_pathway$pathway, "Itch:","neuroinflammation:")

idx = which(target_pathway$pathway %in% c("th2:lichenification","th2:Terminal_Differentiation_and_Lipids"))
target_pathway$pathway[idx] = str_replace(target_pathway$pathway[idx],"th2:","epidermis:")
target_pathway$collection[idx] = "epidermis"

target_pathway$changeInDisease = dz$DZ.change[match(target_pathway$pathway, dz$ID)]

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

rmRows = which(joined$pathway %in% keep_signatures$Target[which(keep_signatures$Exploration == "No")])
joined = joined[-rmRows,]
idx = which(joined$pathway %in% signatureMapping$signature)
joined$pathway[idx] = signatureMapping$New_identifier[match(joined$pathway[idx], signatureMapping$signature)]
  
uploadToBQ(joined, bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_volcano_DiseaseAndDupilumab")
