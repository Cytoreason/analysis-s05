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
PathwayMetadata = Metadata[[3]]
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

### Combine dashboards_PathwayMetadata into dashboards_target_pathway (instead of the Tableau) - issues with pathway names:
PathwayMetadata <- downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_PathwayMetadata", pageSize = 20000)
target_pathway <- downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_target_pathway", pageSize = 20000)

#*** Gil: there are multiple values per pathways in the PathwayMetadata table - probably multiple terms - need to regenerate. The Pathway_level column is not necessary...
# I will use the "dz" dataframe (which was used for disease coverage) and not the PathwayMetadata.
intersect(colnames(PathwayMetadata),colnames(dz))
setdiff(colnames(PathwayMetadata),colnames(dz))
dz <- dz %>%
  dplyr::mutate(adjusted=submodel=="adjusted",Disease="AD") %>%
  dplyr::rename(neglog10FDR=log10_fdr)

dz <- dz[,intersect(colnames(PathwayMetadata),colnames(dz))]  

target_pathway$ID <- target_pathway$pathway

head(target_pathway)
length(setdiff(target_pathway$ID,dz$ID))
length(setdiff(dz$ID,target_pathway$ID))

dz$name_lowercase <- trimws(tolower(dz$ID))
target_pathway$name_lowercase <- trimws(tolower(target_pathway$ID))

intersect(unique(target_pathway$name_lowercase),unique(dz$name_lowercase)) # 2574
setdiff(unique(target_pathway$name_lowercase),unique(dz$name_lowercase)) # 880
setdiff(unique(dz$name_lowercase),unique(target_pathway$name_lowercase)) # 65

grep("\\(",setdiff(unique(target_pathway$name_lowercase),unique(dz$name_lowercase)),ignore.case = T,value = T)
grep("\\(",setdiff(unique(dz$name_lowercase),unique(target_pathway$name_lowercase)),ignore.case = T,value = T)

dz$name_lowercase <- gsub("\\ \\(bont/a\\)|\\ \\(bont/b\\)|\\ \\(bont/c\\)|\\ \\(bont/d\\)|\\ \\(bont/e\\)|\\ \\(bont/f\\)|\\ \\(bont/g\\)","",dz$name_lowercase)
dz$name_lowercase <- gsub("reactome:trif\\(","reactome:trif\\ \\(",dz$name_lowercase)
dz$name_lowercase <- gsub("reactome:activation of irf3/irf7 mediated by tbk1/ikk epsilon","reactome:activation of irf3, irf7 mediated by tbk1, ikkε (ikbke)",dz$name_lowercase)
dz$name_lowercase <- gsub("reactome:initiation of nuclear envelope reformation","reactome:initiation of nuclear envelope (ne) reformation",dz$name_lowercase)
dz$name_lowercase <- gsub("reactome:hsp90 chaperone cycle for steroid hormone receptors \\(shr\\)","reactome:hsp90 chaperone cycle for steroid hormone receptors \\(shr\\) in the presence of ligand",dz$name_lowercase)
dz$name_lowercase <- gsub("reactome:nuclear envelope reassembly","reactome:nuclear envelope (ne) reassembly",dz$name_lowercase)

target_pathway$name_lowercase <- gsub("\\ \\(bota\\)|\\ \\(botb\\)|\\ \\(botc\\)|\\ \\(botd\\)|\\ \\(bote\\)|\\ \\(botf\\)|\\ \\(botg\\)","",target_pathway$name_lowercase)
target_pathway$name_lowercase <- gsub("reactome:defective galnt12 causes crcs1","reactome:defective galnt12 causes colorectal cancer 1 (crcs1)",target_pathway$name_lowercase)
target_pathway$name_lowercase <- gsub("reactome:defective c1galt1c1 causes tnps","reactome:defective c1galt1c1 causes tn polyagglutination syndrome (tnps)",target_pathway$name_lowercase)
target_pathway$name_lowercase <- gsub("reactome:defective galnt3 causes hftc","reactome:defective galnt3 causes familial hyperphosphatemic tumoral calcinosis (hftc)",target_pathway$name_lowercase)
target_pathway$name_lowercase <- gsub("reactome:defective b3galtl causes pps","reactome:defective b3galtl causes peters-plus syndrome (pps)",target_pathway$name_lowercase)
target_pathway$name_lowercase <- gsub("reactome:regulation of lipid metabolism by pparalpha","reactome:regulation of lipid metabolism by peroxisome proliferator-activated receptor alpha (pparalpha)",target_pathway$name_lowercase)
target_pathway$name_lowercase <- gsub("reactome:defective cyp11a1 causes aicsr","reactome:defective cyp11a1 causes adrenal insufficiency, congenital, with 46,xy sex reversal (aicsr)",target_pathway$name_lowercase)
target_pathway$name_lowercase <- gsub("reactome:defective dpm1 causes dpm1-cdg","reactome:defective dpm1 causes dpm1-cdg (cdg-1e)",target_pathway$name_lowercase)
target_pathway$name_lowercase <- gsub("reactome:defective dpm2 causes dpm2-cdg","reactome:defective dpm2 causes dpm2-cdg (cdg-1u)",target_pathway$name_lowercase)
target_pathway$name_lowercase <- gsub("reactome:defective dpm3 causes dpm3-cdg","reactome:defective dpm3 causes dpm3-cdg (cdg-1o)",target_pathway$name_lowercase)
target_pathway$name_lowercase <- gsub("reactome:defective csf2ra causes smdp4","reactome:defective csf2ra causes pulmonary surfactant metabolism dysfunction 4 (smdp4)",target_pathway$name_lowercase)
target_pathway$name_lowercase <- gsub("reactome:defective csf2rb causes smdp5","reactome:defective csf2rb causes pulmonary surfactant metabolism dysfunction 5 (smdp5)",target_pathway$name_lowercase)

intersect(unique(target_pathway$name_lowercase),unique(dz$name_lowercase)) # 2597
setdiff(unique(target_pathway$name_lowercase),unique(dz$name_lowercase)) # 857
setdiff(unique(dz$name_lowercase),unique(target_pathway$name_lowercase)) # 42

# Combine the two tables:
target_pathway$changeInDisease <- NULL
target_pathway_with_dz <- unique(target_pathway %>% left_join(dz[,c("name_lowercase","DZ.change","Disease","adjusted")],by="name_lowercase"))
target_pathway_with_dz$name_lowercase <- NULL
colnames(target_pathway_with_dz)[colnames(target_pathway_with_dz)=="DZ.change"] <- "DZ_change"
colnames(target_pathway_with_dz)[colnames(target_pathway_with_dz)=="submodel"] <- "submodel_correlation"
colnames(target_pathway_with_dz)[colnames(target_pathway_with_dz)=="adjusted"] <- "DZ_is_adjusted"

# There are correlations for both bulk and adjusted (the "submodel_correlation" column), and there are direction in DZ in both bulk and adjusted. We want to simplify and match the direction in DZ with the submodel_correlation:
filtered_target_pathway_with_dz <- target_pathway_with_dz %>%
  filter((submodel_correlation=="adjusted"&DZ_is_adjusted) | (submodel_correlation=="bulk"&!DZ_is_adjusted))

uploadToBQ(filtered_target_pathway_with_dz, bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_target_pathway_Gil")

## Pathway volcano
## ==========================
dz = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "AD_gx_gsa")
dz = dz %>%
  dplyr::filter(!collection %in% c("c2.cgp","c2.cp","c2.cp.biocarta","c2.cp.kegg", "c2.cp.reactome","c3.tft","c7","Neuronal","Mast","Itch")) %>% # removed "epidermis" from here (we need to include it)
  mutate(collection = ifelse(pathway %in% c("Terminal_Differentiation_and_Lipids", "lichenification"), "epidermis",collection))

treatments = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "treatment_PrePostDupi")
treatments = treatments %>%
  dplyr::filter(!collection %in% c("c2.cgp","c2.cp","c2.cp.biocarta","c2.cp.kegg", "c2.cp.reactome","c3.tft","c7","Neuronal","Mast","Itch")) %>% # removed "epidermis" from here (we need to include it)
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

# Check if the pathways in volcano_DiseaseAndDupilumab match the pathways in dashboards_target_pathway_Gil
volcano_DiseaseAndDupilumab <- downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "dashboards_volcano_DiseaseAndDupilumab")
volcano_DiseaseAndDupilumab$ID <- paste0(volcano_DiseaseAndDupilumab$collection,":",volcano_DiseaseAndDupilumab$pathway)

intersect(tolower(unique(target_pathway_with_dz$ID)),tolower(unique(volcano_DiseaseAndDupilumab$ID))) # 2561 shared pathways
setdiff(tolower(unique(target_pathway_with_dz$ID)),tolower(unique(volcano_DiseaseAndDupilumab$ID))) # 893 target_pathway unique pathways
setdiff(tolower(unique(volcano_DiseaseAndDupilumab$ID)),tolower(unique(target_pathway_with_dz$ID))) # 99 svolcano_DiseaseAndDupilumab unique pathways

unique(volcano_DiseaseAndDupilumab$ID[volcano_DiseaseAndDupilumab$collection=="X2"])
unique(target_pathway_with_dz$ID[target_pathway_with_dz$collection=="X2"])


## Upload pathway ontology with the custom signatures
## ==========================
ontology_corrected <- read.csv("~/capsule/scratch/ontology_corrected.csv")
pushToCC(ontology_corrected, tagsToPass = list(list(name="object",value="ontology_corrected")))
# wf-010d0c2186
# wf-85b089d95b

uploadToBQ(ontology_corrected, bqdataset = "s05_atopic_dermatitis", tableName = "ontology_corrected")
