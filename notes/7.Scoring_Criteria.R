# 1. Libraries and data loading
# ------------------------------------------
devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline) # make sure to load version >= 1.0.1
library(tidyverse)
library(scales, include.only = "rescale")

scaling = function(scale_vector, min_val = 0.1, max_val = 1) {
  rescale(x = scale_vector, to = c(min_val, max_val))
}
keep_signatures = openxlsx::read.xlsx("~/analysis-s05/data/Final signatures to be ranked or viewed in the dashboards.xlsx")
signatureMapping = readRDS(get_workflow_outputs("wf-aa75ed069b"))


## Pre-processing the results to include only the signatures and criteria we want
## ====================================================================================
# For all scores but Results_Pathway_PCA we use the previous version in order to have Th1_Related, Th17_Related, Th2_Related
Results_old = readRDS(get_workflow_outputs("wf-47f2a1c1b7"))
Results_old = bind_rows(Results_old[-which(names(Results_old) == "Target_Pathway_PCA")])

Results = bind_rows(readRDS(get_workflow_outputs("wf-47f2a1c1b7"))[["Target_Pathway_PCA"]],
                    Results_old) %>%
  dplyr::filter(Target.Identifier %in% keep_signatures$Target[which(keep_signatures$Exploration == "Yes")])

Results$DataType[which(Results$Type == "adjusted__1__1")] = paste0("Adj ", Results$DataType[which(Results$Type == "adjusted__1__1")])
Results$PC[which(str_detect(Results$DataType, "PCA"))] = sapply(Results$Criteria.Identifier[which(str_detect(Results$DataType, "PCA"))], function(y) strsplit(y,"pc")[[1]][2])
Results$DataType[which(str_detect(Results$DataType, "PCA"))] = paste0(str_replace(Results$DataType[which(str_detect(Results$DataType, "PCA"))], "PCA","PC"),Results$PC[which(str_detect(Results$DataType, "PCA"))])
Results$DataType[which(Results$DataType == "Target_CS")] = paste0("Target_CS_",Results$Criteria.Identifier[which(Results$DataType == "Target_CS")])
Results = Results[-which(str_detect(Results$DataType,"Target_Pathway_PC") & 
                           Results$Type == "bulk" & 
                         str_detect(Results$Criteria.Identifier,"_adjusted")),]
Results = Results[-which(str_detect(Results$DataType,"Target_Pathway_PC") & 
                           Results$Type == "adjusted__1__1" & 
                         str_detect(Results$Criteria.Identifier,"_bulk")),]

chosenCriteria = c("Disease Enrichment in L_vs_HC" = "Enrichment in L_vs_HC",
                   "Disease Enrichment in NL_vs_HC" = "Enrichment in NL_vs_HC",
                   "Disease Enrichment in L_vs_NL (Adjusted)" = "Adj Enrichment in L_vs_NL",
                   "Disease Enrichment in L_vs_HC (Adjusted)" = "Adj Enrichment in L_vs_HC",
                   "Inflammation"  = "Target_Pathway_PC1_DZ_vs_HC_bulk_keyWithTh2",
                   "Immune Cell Activation"  = "Target_Pathway_PC2_DZ_vs_HC_bulk_keyWithTh2",
                   "Proliferation and Cell Cycle"  = "Adj Target_Pathway_PC1_DZ_vs_HC_adjusted_keyWithTh2",
                   "Correlation to EASI"  = "Target_CS_EASI")

Results_filtered = Results[which(Results$DataType %in% chosenCriteria),-13]
Results_filtered$Criteria.Identifier = names(chosenCriteria)[match(Results_filtered$DataType, chosenCriteria)]
# flip sign
idx = which(Results_filtered$Criteria.Identifier == "Proliferation and Cell Cycle")
Results_filtered$metricValue[idx] = (-1) * Results_filtered$metricValue[idx]

pushToCC(Results_filtered)
# wf-563b1691b5
# wf-45bed4e3de


# Integrate additional criteria
# ======================================
additional_criteria = readRDS(get_workflow_outputs("wf-7a1d7b9ecb"))
additional_criteria = additional_criteria[which(additional_criteria$Type == "bulk"),]
# wf-fdba9ddeaf

scorad = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "Results_Target_SCORAD")
scorad = scorad %>%
  dplyr::filter(dataset == "GSE130588__GPL570") %>%
  dplyr::filter(Target_Identifier %in% keep_signatures$Target[which(keep_signatures$Ranking == "Yes")]) %>%
  dplyr::select(-n_observation,-dataset) %>%
  mutate(DataType = "Target_CS_SCORAD", hit = NA)
colnames(scorad) = str_replace(colnames(scorad),"_",".")
colnames(scorad)[which(colnames(scorad) == "log10.fdr")] = "log10_fdr"
scorad$Criteria.Identifier = "Correlation to SCORAD"

scorad$Target.ID = signatureMapping$ID[match(scorad$Target.ID, signatureMapping$previousID)]
# wf-426aa560ee

# Integrate network criteria
# ======================================
topology = readRDS(get_workflow_outputs("wf-8556686d74"))
topology = topology %>%
  mutate(DataType = "Network Density") %>%
  dplyr::filter(Type == "bulk") %>%
  dplyr::filter(Target.ID %in% unique(Results_filtered$Target.ID)) %>%
  mutate(pvalue = ifelse(pvalue == 0, 1/1001, pvalue)) %>%
  mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
  mutate(log10_fdr = -log10(fdr)) %>%
  mutate(Criteria.Identifier = "AD Network Density")


centrality = readRDS(get_workflow_outputs("wf-72d490c885"))
centrality = centrality %>%
  mutate(DataType = "Network Centrality") %>%
  dplyr::filter(Type == "bulk") %>%
  dplyr::filter(Criteria.Identifier == "eigen_centrality_median") %>%
  dplyr::filter(Target.ID %in% unique(Results_filtered$Target.ID)) %>%
  mutate(pvalue = ifelse(pvalue == 0, 1/1001, pvalue)) %>%
  mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
  mutate(log10_fdr = -log10(fdr)) %>%
  mutate(Criteria.Identifier = "AD Network Centrality")

network = bind_rows(topology, centrality)
network$Target.Identifier = signatureMapping$New_identifier[match(network$Target.Identifier, signatureMapping$signature)]
pushToCC(network, tagsToPass = list(list(name="object",value="network_criteria")))
# wf-ec34d7ef49

## Export
Network = readRDS(get_workflow_outputs("wf-ec34d7ef49"))
write.csv(Network, "~/exportedFiles/network.csv", row.names = F, quote = F)

## All criteria
## ====================
Results_filtered = readRDS(get_workflow_outputs("wf-45bed4e3de"))
scorad = readRDS(get_workflow_outputs("wf-426aa560ee"))
additional_criteria = readRDS(get_workflow_outputs("wf-fdba9ddeaf"))
  additional_criteria$fdr = additional_criteria$log10_fdr = NA # we don't trust it
network = readRDS(get_workflow_outputs("wf-ec34d7ef49"))

allCriteria = bind_rows(Results_filtered,
                        scorad,
                        additional_criteria,
                        network)
pushToCC(allCriteria, tagsToPass = list(list(name="object",value="allCriteria")))
# wf-7bb3d94abc
# wf-a081e251f2
# wf-af22471107


## Scaling each criterion
## ====================================
# Almost zero negatives and multiply FDR and effect
floor_dec <- function(x, k = 2) { # function to take the lowest number above 0 and round it down
  minimal_value = x[which(x > 0)]
  minimal_value = floor(min(minimal_value) * 10^k) / 10^k
  x[which(x < 0)] <- minimal_value
  return(x)
}

allCriteria_chosenTargets = allCriteria[which(allCriteria$Target.Identifier %in% keep_signatures$Target[which(keep_signatures$Ranking == "Yes")]),]

# First we set a threshold for FDR
allCriteria_chosenTargets$log10_fdr = ifelse(allCriteria_chosenTargets$log10_fdr > 4, 4, allCriteria_chosenTargets$log10_fdr)

scores = allCriteria_chosenTargets %>%
  group_by(Criteria.Identifier) %>%
  mutate(score = ifelse(is.na(fdr), metricValue, metricValue * log10_fdr)) %>% # because we have NAs in FDR, score is only multiplied when we have it
  mutate(rounded_score = floor_dec(score)) %>% # zeroing all negative scores
  mutate(scaled_score_rounded = scaling(rounded_score)) %>% # scaling
  ungroup()

pushToCC(scores, tagsToPass = list(list(name="object",value="scores")))
# wf-5b6481c61b
# wf-4bb07f3ef3
# wf-86663471a4

scores$Target.ID = factor(scores$Target.ID, ordered = T, levels = names(targetColors))

row_ord <- order(as.integer(scores$Target.ID))
char_levels_in_factor_order <- unique(scores$Target.Identifier[row_ord])
scores$Target.Identifier <- factor(scores$Target.Identifier,
                      levels = char_levels_in_factor_order,
                      ordered = TRUE)


rankings = scores %>%
  group_by(Target.Identifier, Target.ID) %>%
  summarise(summed_scaled_score_rounded = sum(scaled_score_rounded), .groups = "drop") %>%
  mutate(rank_scaled_score_rounded = rank(-summed_scaled_score_rounded))

pushToCC(rankings, tagsToPass = list(list(name="object",value="rankings")))
# wf-13b22e4fd7
# wf-bbd99922f1
# wf-0873211b4e
# wf-22cecea41a

## Visualization
## =======================
# Score per criteria
ggplot(scores, aes(x = Target.Identifier, y = Criteria.Identifier)) +
  spot.theme +
  geom_point(aes(size = scaled_score_rounded, color = Target.ID))+
  facet_grid(DataType ~ . ,scales = "free",space = "free")+
  scale_color_manual(values = targetColors)
ggsave("~/analysis-s05/figures/Results/scores_rounded.png", scale = 1.2, width = 13, height = 8, bg = "white")


# Rank
ggplot(rankings, aes(x = Target.Identifier, y = summed_scaled_score_rounded, fill = Target.ID)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(data = function(x) subset(x, rank_scaled_score_rounded <= 16), aes(y = 11, label = rank_scaled_score_rounded), check_overlap = T) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  ggpubr::border() +
  scale_fill_manual(values = targetColors) +
  guides(fill = F) +
  labs(x = NULL, y = "Summed Scores")
ggsave("~/analysis-s05/figures/Results/ranking_rounded.png", scale = 1, width = 12, height = 6, bg = "white")


#############################################
# SCORING COMPLEMENTARITY TO DUPILUMAB
#############################################
# Almost zero negatives and multiply FDR and effect
floor_dec <- function(x, k = 2) { # function to take the lowest number above 0 and round it down
  minimal_value = x[which(x > 0)]
  minimal_value = floor(min(minimal_value) * 10^k) / 10^k
  x[which(x < 0)] <- minimal_value
  return(x)
}

keep_signatures = openxlsx::read.xlsx("~/analysis-s05/data/Final signatures to be ranked or viewed in the dashboards.xlsx")
signatureMapping = readRDS(get_workflow_outputs("wf-aa75ed069b"))
Results_filtered = readRDS(get_workflow_outputs("wf-45bed4e3de"))

treatment = readRDS(get_workflow_outputs("wf-e3e631e2c9"))
treatment$Target.Identifier = signatureMapping$New_identifier[match(treatment$Target.ID, signatureMapping$ID)]
overlap = readRDS(get_workflow_outputs("wf-3c7ac8ddd1"))
coverage = readRDS(get_workflow_outputs("wf-4c7ecb1fc6"))
whiteSpace_coverage = readRDS(get_workflow_outputs("wf-97923f7e39"))

allCriteria = bind_rows(treatment, overlap, coverage, whiteSpace_coverage) %>%
  dplyr::filter(Type == "bulk") %>%
  dplyr::filter(Target.Identifier %in% keep_signatures$Target[which(keep_signatures$Ranking == "Yes")])


# First we set a threshold for FDR
allCriteria$log10_fdr = ifelse(allCriteria$log10_fdr > 4, 4, allCriteria$log10_fdr)

# We decided to change the score for the treatment and non-responder criteria: 1/(metricValue*log10_fdr) and to zero if the enrichment in L_vs_HC is not significant (and/or not reversed in direction - e.g. down in disease and down in treatment...)
# In the treatment object, the metricValue was already multiply by -1, but in the Disease Enrichment in L_vs_HC" it wasn't, so the same sign means reversed by treatment in this case...)
enrichment_in_disease <- Results_filtered %>% # use the L_vs_HC results to zero non-flipping targets, or targets that are not significant in disease
  dplyr::filter(Criteria.Identifier == "Disease Enrichment in L_vs_HC") %>%
  dplyr::select(Target.Identifier,
         disease_metricValue = metricValue,
         disease_fdr = fdr) %>%
  dplyr::distinct(Target.Identifier, .keep_all = TRUE)

scores <- allCriteria %>%
  dplyr::group_by(Criteria.Identifier) %>%
  dplyr::mutate(score = dplyr::case_when(
    stringr::str_detect(Criteria.Identifier, "Coverage") ~ metricValue,
    stringr::str_detect(Criteria.Identifier, "Shared") ~ metricValue * log10_fdr,
    TRUE ~ (1/(metricValue * log10_fdr))
  )) %>%
  dplyr::left_join(enrichment_in_disease, by = "Target.Identifier") %>%
  dplyr::mutate(score = case_when(
    stringr::str_detect(Criteria.Identifier, "Dupilumab") &
      !is.na(disease_fdr) &
      (
        disease_fdr > 0.05 |
          (
            disease_fdr <= 0.05 &
              !is.na(disease_metricValue) &
              sign(metricValue) != 0 & sign(disease_metricValue) != 0 &
              sign(metricValue) != sign(disease_metricValue)
          )
      ) ~ 0,
    TRUE ~ score
  )) %>%
  dplyr::select(-disease_metricValue, -disease_fdr) %>%
  dplyr::mutate(rounded_score = floor_dec(score)) %>%
  dplyr::mutate(scaled_score_rounded = scaling(rounded_score)) %>%
  dplyr::ungroup()

pushToCC(scores, tagsToPass = list(list(name="object",value="scores_complementarity")))
# wf-651a3f9622
# wf-56345249fc
# wf-0030d77e02
# wf-070698da21
# wf-34b85ee0cc - with IL13

scores$Target.ID = factor(scores$Target.ID, ordered = T, levels = names(targetColors))

row_ord <- order(as.integer(scores$Target.ID))
char_levels_in_factor_order <- unique(scores$Target.Identifier[row_ord])
scores$Target.Identifier <- factor(scores$Target.Identifier,
                                   levels = char_levels_in_factor_order,
                                   ordered = TRUE)


rankings = scores %>%
  group_by(Target.Identifier, Target.ID) %>%
  summarise(summed_scaled_score_rounded = sum(scaled_score_rounded), .groups = "drop") %>%
  mutate(rank_scaled_score_rounded = rank(-summed_scaled_score_rounded, ties.method = "random"))

pushToCC(rankings, tagsToPass = list(list(name="object",value="rankings")))
# wf-27897b1bcc
# wf-c160239658
# wf-9735cfa614
# wf-5c79873922
# wf-52b9d9053a
# wf-24ecca8512 - with IL13


# Score per criteria
ggplot(scores, aes(x = Target.Identifier, y = Criteria.Identifier)) +
  spot.theme +
  geom_point(aes(size = scaled_score_rounded, color = Target.ID))+
  facet_grid(DataType ~ . ,scales = "free",space = "free")+
  scale_color_manual(values = targetColors)
ggsave("~/analysis-s05/figures/Results/scores_rounded.png", scale = 1.2, width = 13, height = 8, bg = "white")


# Rank
ggplot(rankings, aes(x = Target.Identifier, y = summed_scaled_score_rounded, fill = Target.ID)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(data = function(x) subset(x, rank_scaled_score_rounded <= 15), aes(y = 7, label = rank_scaled_score_rounded), check_overlap = T) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  ggpubr::border() +
  scale_fill_manual(values = targetColors) +
  guides(fill = F) +
  labs(x = NULL, y = "Summed Scores")
ggsave("~/analysis-s05/figures/Results/ranking_rounded.png", scale = 1, width = 12, height = 6, bg = "white")
