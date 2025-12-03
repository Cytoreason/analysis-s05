# 1. Libraries and data loading
# ------------------------------------------
devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline) # make sure to load version >= 1.0.1
library(tidyverse)
library(scales, include.only = "rescale")

scaling = function(scale_vector, min_val = 0.1, max_val = 1) {
  rescale(x = scale_vector, to = c(min_val, max_val))
}


## Pre-processing the results to include only the signatures and criteria we want
## ====================================================================================
Results = readRDS(get_workflow_outputs("wf-9a4e8e1dba"))
Results = bind_rows(Results) %>%
  dplyr::filter(!Target.Collection %in% c("Ligands","Mast","Itch","Neuronal")) %>%
  dplyr::filter(!str_detect(Target.Identifier, c("insilico|x2_|PAMP|CST|SP|aIgE|Icatibant|X2 late activation|X2 early activated|X2 early activation"))) %>%
  dplyr::filter(!Target.Identifier == "X2 late inhibition")

Results$DataType[which(Results$Type == "adjusted__1__1")] = paste0("Adj ", Results$DataType[which(Results$Type == "adjusted__1__1")])
Results$PC[which(str_detect(Results$DataType, "PCA"))] = sapply(Results$Criteria.Identifier[which(str_detect(Results$DataType, "PCA"))], function(y) strsplit(y,"pc")[[1]][2])
Results$DataType[which(str_detect(Results$DataType, "PCA"))] = paste0(str_replace(Results$DataType[which(str_detect(Results$DataType, "PCA"))], "PCA","PC"),Results$PC[which(str_detect(Results$DataType, "PCA"))])
Results$DataType[which(Results$DataType == "Target_CS")] = paste0("Target_CS_",Results$Criteria.Identifier[which(Results$DataType == "Target_CS")])
Results = Results[-which(str_detect(Results$DataType,"Target_Pathway_PC") & 
                         Results$Type == "bulk" & 
                         str_detect(Results$Criteria.Identifier,"adj_")),]
Results = Results[-which(str_detect(Results$DataType,"Target_Pathway_PC") & 
                         Results$Type == "adjusted__1__1" & 
                         str_detect(Results$Criteria.Identifier,"(?<!_)pathway")),]

chosenCriteria = c("Enrichment in L_vs_HC" = "Enrichment in L_vs_HC",
                   "Enrichment in NL_vs_HC" = "Enrichment in NL_vs_HC",
                   "Enrichment in L_vs_NL (Adjusted)" = "Adj Enrichment in L_vs_NL",
                   "Enrichment in L_vs_HC (Adjusted)" = "Adj Enrichment in L_vs_HC",
                   "Correlation with Pathway Meta PC1"  = "Target_Pathway_PC1",
                   "Correlation with Pathway Meta PC2"  = "Target_Pathway_PC2",
                   "Correlation with Pathway Meta PC1 (Adjusted)"  = "Adj Target_Pathway_PC1",
                   "Correlation with EASI"  = "Target_CS_EASI")

Results_filtered = Results[which(Results$DataType %in% chosenCriteria),-13]
Results_filtered$Criteria.Identifier = names(chosenCriteria)[match(Results_filtered$DataType, chosenCriteria)]


## Scaling each criterion
## ====================================
# First we set a threshold for FDR
Results_filtered$log10_fdr = ifelse(Results_filtered$log10_fdr > 4, 4, Results_filtered$log10_fdr)

# Multiply FDR and effect
Results_filtered$score = Results_filtered$metricValue * Results_filtered$log10_fdr

# Per criteria we re-scale between 0.1 and 1
scores = Results_filtered %>%
  group_by(DataType) %>%
  mutate(scaled_score = scaling(score)) %>%
  mutate(scaled_effect = scaling(metricValue)) %>%
  mutate(scaled_fdr = scaling(log10_fdr)) %>%
  ungroup()

scores$Target.ID = factor(scores$Target.ID, ordered = T, levels = names(targetColors))

row_ord <- order(as.integer(scores$Target.ID))
char_levels_in_factor_order <- unique(scores$Target.Identifier[row_ord])
scores$Target.Identifier <- factor(scores$Target.Identifier,
                      levels = char_levels_in_factor_order,
                      ordered = TRUE)


## Visualization
## =======================
## 1. Based on effect * fdr
## ---------------------------------
# Score per criteria
ggplot(scores, aes(x = Target.Identifier, y = Criteria.Identifier)) +
  spot.theme +
  geom_point(aes(size = scaled_score, color = Target.ID))+
  facet_grid(DataType ~ . ,scales = "free",space = "free")+
  scale_color_manual(values = targetColors)
ggsave("~/analysis-s05/figures/Results/scores.png", scale = 1.2, width = 12, height = 6, bg = "white")


# Rank
rankings = scores %>%
  group_by(Target.Identifier, Target.ID) %>%
  summarise(summedScore = sum(scaled_score)) %>%
  ungroup() %>%
  mutate(ranking = rank(-summedScore))

ggplot(rankings, aes(x = Target.Identifier, y = summedScore, fill = Target.ID)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(data = function(x) subset(x, ranking <= 10), aes(y = 7, label = ranking), check_overlap = T) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  ggpubr::border() +
  scale_fill_manual(values = targetColors) +
  guides(fill = F) +
  labs(x = NULL, y = "Summed Scores")
ggsave("~/analysis-s05/figures/Results/ranking.png", scale = 1, width = 12, height = 6, bg = "white")


## 2. Based on effect alone
## ---------------------------------
# Score per criteria
ggplot(scores, aes(x = Target.Identifier, y = Criteria.Identifier)) +
  spot.theme +
  geom_point(aes(size = scaled_effect, color = Target.ID))+
  facet_grid(DataType ~ . ,scales = "free",space = "free")+
  scale_color_manual(values = targetColors)
ggsave("~/analysis-s05/figures/Results/scores_byEffect.png", scale = 1.2, width = 12, height = 6, bg = "white")


# Rank
rankings = scores %>%
  group_by(Target.Identifier, Target.ID) %>%
  summarise(summedScore = sum(scaled_effect)) %>%
  ungroup() %>%
  mutate(ranking = rank(-summedScore))

ggplot(rankings, aes(x = Target.Identifier, y = summedScore, fill = Target.ID)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(data = function(x) subset(x, ranking <= 10), aes(y = 7, label = ranking), check_overlap = T) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  ggpubr::border() +
  scale_fill_manual(values = targetColors) +
  guides(fill = F) +
  labs(x = NULL, y = "Summed Scores")
ggsave("~/analysis-s05/figures/Results/ranking_byEffect.png", scale = 1, width = 12, height = 6, bg = "white")



## 3. Based on FDR alone
## ---------------------------------
# Score per criteria
ggplot(scores, aes(x = Target.Identifier, y = Criteria.Identifier)) +
  spot.theme +
  geom_point(aes(size = scaled_fdr, color = Target.ID))+
  facet_grid(DataType ~ . ,scales = "free",space = "free")+
  scale_color_manual(values = targetColors)
ggsave("~/analysis-s05/figures/Results/scores_byFDR.png", scale = 1.2, width = 12, height = 6, bg = "white")


# Rank
rankings = scores %>%
  group_by(Target.Identifier, Target.ID) %>%
  summarise(summedScore = sum(scaled_fdr)) %>%
  ungroup() %>%
  mutate(ranking = rank(-summedScore))

ggplot(rankings, aes(x = Target.Identifier, y = summedScore, fill = Target.ID)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(data = function(x) subset(x, ranking <= 10), aes(y = 7, label = ranking), check_overlap = T) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  ggpubr::border() +
  scale_fill_manual(values = targetColors) +
  guides(fill = F) +
  labs(x = NULL, y = "Summed Scores")
ggsave("~/analysis-s05/figures/Results/ranking_byFDR.png", scale = 1, width = 12, height = 6, bg = "white")
