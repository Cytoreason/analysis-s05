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
Results = readRDS(get_workflow_outputs("wf-e0db493d42"))
Results = bind_rows(Results) %>%
  dplyr::filter(!Target.Collection %in% c("Ligands","Mast","Itch","Neuronal")) %>%
  dplyr::filter(!str_detect(Target.Identifier, c("insilico|x2_|PAMP|CST|SP|BMP7|Icatibant|X2 late activation|X2 early activated|X2 early activation"))) %>%
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
# wf-563b1691b5

# Add correlation to median th2, epidermis, neuroinflammation as criteria
Results_th2 = readRDS(get_workflow_outputs("wf-64470d2a55"))
Results_th2 = Results_th2$Target_Pathway %>%
  dplyr::filter(!Target.Collection %in% c("Ligands","Mast","Itch","Neuronal","th2", "epidermis", "neuroinflammation")) %>%
  dplyr::filter(!str_detect(Target.Identifier, c("insilico|x2_|PAMP|CST|SP|BMP7|Icatibant|X2 late activation|X2 early activated|X2 early activation"))) %>%
  dplyr::filter(!Target.Identifier == "X2 late inhibition") %>%
  dplyr::filter(Criteria.Collection %in% c("th2", "epidermis", "neuroinflammation"))
Results_th2$Criteria.Collection[which(Results_th2$Criteria.Identifier %in% c("th2:lichenification","th2:Terminal_Differentiation_and_Lipids"))] <- "epidermis"
Results_th2$Criteria.Identifier <- str_replace(Results_th2$Criteria.Identifier, "th2:lichenification", "epidermis:lichenification")
Results_th2$Criteria.Identifier <- str_replace(Results_th2$Criteria.Identifier, "th2:Terminal_Differentiation_and_Lipids", "epidermis:Terminal_Differentiation_and_Lipids")

pathwaysToNotFlipSign = c("epidermis:reactome__Asymmetric localization of PCP proteins", "epidermis:reactome__Differentiation of keratinocytes in interfollicular epidermis in mammalian skin")
idx = which(Results_th2$Criteria.Collection == "epidermis" & !Results_th2$Criteria.Identifier %in% pathwaysToNotFlipSign & Results_th2$Type == "bulk")
Results_th2$metricValue[idx] = (-1) * Results_th2$metricValue[idx]
# wf-8606c24d32

# # Integrate network criteria
# topology = readRDS(get_workflow_outputs("wf-8556686d74"))
# topology = topology %>%
#   mutate(DataType = "Network Density") %>%
#   dplyr::filter(Type == "bulk") %>%
#   dplyr::filter(Target.ID %in% unique(Results_filtered$Target.ID)) %>%
#   mutate(pvalue = ifelse(pvalue == 0, 1/1001, pvalue)) %>%
#   mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
#   mutate(log10_fdr = -log10(fdr))
#   
# 
# centrality = readRDS(get_workflow_outputs("wf-72d490c885"))
# centrality = centrality %>%
#   mutate(DataType = "Network Centrality") %>%
#   dplyr::filter(Type == "bulk") %>%
#   dplyr::filter(Criteria.Identifier == "eigen_centrality_median") %>%
#   dplyr::filter(Target.ID %in% unique(Results_filtered$Target.ID)) %>%
#   mutate(pvalue = ifelse(pvalue == 0, 1/1001, pvalue)) %>%
#   mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
#   mutate(log10_fdr = -log10(fdr))
# 
# allCriteria = bind_rows(Results_filtered, topology, centrality)
# pushToCC(allCriteria, tagsToPass = list(list(name="object",value="allCriteria")))
# # wf-fc0d0d2464


## Scaling each criterion
## ====================================
# First we set a threshold for FDR
Results_filtered$log10_fdr = ifelse(Results_filtered$log10_fdr > 4, 4, Results_filtered$log10_fdr)
Results_th2$log10_fdr = ifelse(Results_th2$log10_fdr > 4, 4, Results_th2$log10_fdr)

# Multiply FDR and effect
Results_filtered$score = Results_filtered$metricValue * Results_filtered$log10_fdr
Results_th2$score = Results_th2$metricValue * Results_th2$log10_fdr

# Almost zero negatives and multiply FDR and effect
floor_dec <- function(x, k = 2) { # function to take the lowest number above 0 and round it down
  minimal_value = x[which(x > 0)]
  minimal_value = floor(min(minimal_value) * 10^k) / 10^k
  x[which(x < 0)] <- minimal_value
  return(x)
}


# Per criteria we re-scale between 0.1 and 1
Results_th2 = Results_th2 %>%
  group_by(Target.ID, Target.Identifier, Target.Collection, Type, Criteria.Collection, metricType, hit) %>%
  summarise(metricValue = median(metricValue),
            fdr = median(fdr),
            score = median(score)) %>%
  mutate(log10_fdr = -log10(fdr),
         Criteria.Identifier = paste0("Median Correlation with ",Criteria.Collection, ifelse(Type == "bulk", ""," (Adjusted)")),
         DataType = paste0(Criteria.Collection, ifelse(Type == "bulk", ""," (Adjusted)"))) %>%
  ungroup()


scores = rbind(Results_filtered, Results_th2) %>%
  group_by(DataType) %>%
  mutate(rounded_effect = ifelse(DataType == "Adj Target_Pathway_PC1", floor_dec(metricValue, 3), floor_dec(metricValue))) %>%
  mutate(scaled_score = scaling(score)) %>%
  mutate(scaled_effect = scaling(metricValue)) %>%
  mutate(scaled_fdr = scaling(log10_fdr)) %>%
  mutate(scaled_score_rounded = scaling(rounded_effect)) %>%
  ungroup()

scores$Target.ID = factor(scores$Target.ID, ordered = T, levels = names(targetColors))

row_ord <- order(as.integer(scores$Target.ID))
char_levels_in_factor_order <- unique(scores$Target.Identifier[row_ord])
scores$Target.Identifier <- factor(scores$Target.Identifier,
                      levels = char_levels_in_factor_order,
                      ordered = TRUE)

pushToCC(scores, tagsToPass = list(list(name="object",value="scores_withTh2median")))
# wf-5b6481c61b

rankings = scores %>%
  group_by(Target.Identifier, Target.ID) %>%
  summarise(across(c(scaled_score, scaled_effect, scaled_fdr, scaled_score_rounded),
                   ~ sum(.x, na.rm = TRUE),
                   .names = "summed_{.col}"),
            .groups = "drop") %>%
  mutate(across(starts_with("summed_"), ~ rank(-.x), .names = "rank_{sub('^summed_', '', .col)}"))

pushToCC(rankings, tagsToPass = list(list(name="object",value="rankings_withTh2median")))
# wf-13b22e4fd7

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
ggsave("~/analysis-s05/figures/Results/scores.png", scale = 1.2, width = 13, height = 8, bg = "white")


# Rank
ggplot(rankings, aes(x = Target.Identifier, y = summed_scaled_score, fill = Target.ID)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(data = function(x) subset(x, rank_scaled_score <= 10), aes(y = 11, label = rank_scaled_score), check_overlap = T) +
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
ggsave("~/analysis-s05/figures/Results/scores_byEffect.png", scale = 1.2, width = 13, height = 8, bg = "white")


# Rank
ggplot(rankings, aes(x = Target.Identifier, y = summed_scaled_effect, fill = Target.ID)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(data = function(x) subset(x, rank_scaled_effect <= 10), aes(y = 12, label = rank_scaled_effect), check_overlap = T) +
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
ggsave("~/analysis-s05/figures/Results/scores_byFDR.png", scale = 1.2, width = 13, height = 8, bg = "white")


# Rank
ggplot(rankings, aes(x = Target.Identifier, y = summed_scaled_fdr, fill = Target.ID)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(data = function(x) subset(x, rank_scaled_fdr <= 10), aes(y = 11, label = rank_scaled_fdr), check_overlap = T) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  ggpubr::border() +
  scale_fill_manual(values = targetColors) +
  guides(fill = F) +
  labs(x = NULL, y = "Summed Scores")
ggsave("~/analysis-s05/figures/Results/ranking_byFDR.png", scale = 1, width = 12, height = 6, bg = "white")



## 1. Based on effect * fdr but zeroing negatives
## --------------------------------------------------
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
  geom_text(data = function(x) subset(x, rank_scaled_score_rounded <= 18), aes(y = 11, label = rank_scaled_score_rounded), check_overlap = T) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  ggpubr::border() +
  scale_fill_manual(values = targetColors) +
  guides(fill = F) +
  labs(x = NULL, y = "Summed Scores")
ggsave("~/analysis-s05/figures/Results/ranking_rounded.png", scale = 1, width = 12, height = 6, bg = "white")



