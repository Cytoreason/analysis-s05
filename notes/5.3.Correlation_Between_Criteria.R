devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(ComplexHeatmap)
library(tidyverse)


Results = readRDS(get_workflow_outputs("wf-63a50469cf"))
Results$Target_Cell_PCA$PC = sapply(Results$Target_Cell_PCA$Criteria.Identifier, function(x) str_split(x,"pc")[[1]][2])
Results$Target_Pathway_PCA$PC = sapply(Results$Target_Pathway_PCA$Criteria.Identifier, function(x) str_split(x,"pc")[[1]][2])

Results_binded = Results[setdiff(names(Results), c("Target_Cell", "Target_Pathway", "Target_Gene"))] %>%
  bind_rows() %>%
  mutate(DataType = ifelse(DataType == "Target_CS", paste0("Target_CS_",Criteria.Identifier), DataType)) %>%
  mutate(DataType = ifelse(DataType == "Target_Cell_PCA", paste0("Target_Cell_PC",PC), DataType)) %>%
  mutate(DataType = ifelse(DataType == "Target_Pathway_PCA", paste0("Target_Pathway_PC",PC), DataType)) %>%
  mutate(DataType = str_replace_all(DataType, " ", "_"))

# remove bulk & adj mix
Results_binded = Results_binded[-which(str_detect(Results_binded$DataType,"Target_Pathway_PC") & 
                                        Results_binded$Type == "bulk" & 
                                        str_detect(Results_binded$Criteria.Identifier,"adj_")),]
Results_binded = Results_binded[-which(str_detect(Results_binded$DataType,"Target_Pathway_PC") & 
                                         Results_binded$Type == "adjusted__1__1" & 
                                         str_detect(Results_binded$Criteria.Identifier,"(?<!_)pathway")),]

Results_binded$DataType[which(Results_binded$Type != "bulk")] <-paste0("Adj_", Results_binded$DataType[which(Results_binded$Type != "bulk")])
Results_binded = Results_binded %>%
  dplyr::filter(!str_detect(DataType, "PC3")) %>%
  dplyr::filter(!str_detect(DataType, "Adj_Target_Pathway_PC2|Adj_Target_Cell_PC|Cell_PC2|Adj_Target_CS|Adj_Target_MS"))

Results.corr = reshape2::dcast(Results_binded, Target.Identifier ~ DataType, value.var = "metricValue") %>%
    column_to_rownames(var = "Target.Identifier") %>%
    as.matrix

corr_corr = lapply(Results.corr, cor, method = "spearman")
corr_corr = cor(Results.corr, method = "pearson")

png("~/analysis-s05/figures/Results/criteria_correlation_pearson.png", bg = "white", res = "120", width = 1200, height = 1000)
Heatmap(corr_corr, cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(round(corr_corr[i, j],1), x, y, gp = gpar(fontsize = 10))},
  name = "Pearson\nCorrelation", row_labels = str_replace_all(rownames(corr_corr),"_"," "), 
  column_labels = str_replace_all(rownames(corr_corr),"_"," "))
dev.off()

png("~/analysis-s05/figures/Results/criteria_correlation_adj.png", bg = "white", res = "120", width = 1200, height = 1000)
Heatmap(corr_corr$adjusted__1__1, cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(round(corr_corr$adjusted__1__1[i, j],1), x, y, gp = gpar(fontsize = 10))},
  name = "Spearman\nCorrelation", row_labels = str_replace_all(rownames(corr_corr$adjusted__1__1),"_"," "), 
  column_labels = str_replace_all(rownames(corr_corr$adjusted__1__1),"_"," "), column_title = "Adjusted")
dev.off()

## Correlation between molecular score and EASI scores
## =========================================================
ccm_wfid = "wf-8e948630d7"
ccm= as_ccm_fit(ccm_wfid)

pathways = c("time","sample_classification","condition","experiment_id","EASI","MolecularScore")
enrichment =  lapply(ccm$datasets[which(names(ccm$datasets) %in% c("GSE130588__GPL570","GSE59294__GPL570"))], function(d){
    pData(assayDataExpression(d))[,pathways]
}) %>% do.call(rbind,.)

enrichment_filtered= enrichment %>%
  dplyr::filter(time %in% c("W0","Pre") & sample_classification == "Lesion" & condition == "atopic dermatitis")

ggplot(enrichment_filtered, aes(x = EASI, y = MolecularScore)) +
  geom_point(aes(color = experiment_id), size =3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  ggpubr::stat_cor() +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(title = "Correlaton between EASI scores and Molecular Scores")
ggsave("~/analysis-s05/figures/AD Model/corr_EASI_MS_lesional.png", bg = "white")


enrichment_filtered= enrichment %>%
  dplyr::filter(time %in% c("W0","Pre") & condition == "atopic dermatitis")

ggplot(enrichment_filtered, aes(x = EASI, y = MolecularScore)) +
  geom_point(aes(color = experiment_id), size =3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  ggpubr::stat_cor() +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(title = "Correlaton between EASI scores and Molecular Scores")
ggsave("~/analysis-s05/figures/AD Model/corr_EASI_MS_all.png", bg = "white")

ggplot(enrichment_filtered, aes(x = EASI, y = MolecularScore)) +
  geom_point(aes(color = sample_classification), size =3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  ggpubr::stat_cor() +
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "Correlaton between EASI scores and Molecular Scores")
ggsave("~/analysis-s05/figures/AD Model/corr_EASI_MS_all.png", bg = "white")
