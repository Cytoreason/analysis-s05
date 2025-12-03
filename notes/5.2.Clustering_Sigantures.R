### In this part we cluster the signatures pathway correlations
### with the purpose of grouping them into (hopefully) functional groups - 
### and we characterize each group by the most significant pathways they target

# 1. Libraries and data loading
# ------------------------------------------
devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)
library(ggdendro)

Results = readRDS(get_workflow_outputs("wf-9a4e8e1dba"))

clusterTargets = function(res, plotSuffix) {
  clust = hclust(dist(res)) 

  # Make the dendrogram prettier
  dend = dendro_data(as.dendrogram(clust), type = "rectangle")
  d = ggdendrogram(dend, theme_dendro = T)
  ggsave(plot = d, paste0("~/analysis-s05/figures/Results/clusteringDendrogram_",plotSuffix,".png"), width = 10, height = 4, bg = "white")
}


## Trials
## ===========================
set.seed(1234)

# 1. All cells and pathways, bulk
results.partial = do.call(rbind, Results[c("Target_Cell","Target_Pathway")])
results.partial = results.partial %>%
  dplyr::filter(!Target.Collection %in% c("Ligands","Mast")) %>%
  dplyr::filter(!str_detect(Target.Identifier, c("Jha|insilico|PAMP|Icatibant|CST14|SP|x2_activated_inhibition_late|x2_general_inhibition_late|x2_activation_late"))) %>%
  dplyr::filter(Type == "bulk")

results.wide = reshape2::dcast(results.partial, Target.Identifier ~ Criteria.Identifier, value.var = "metricValue")
rownames(results.wide) = results.wide$Target.Identifier
results.wide = results.wide[,-1]

clusterTargets(results.wide, "bulk_allPathways_allCells")

clust = hclust(dist(results.wide)) 
clust = cutree(clust, k = 5)


# 2. All cells and pathways, bulk & adjusted
results.partial = do.call(rbind, Results[c("Target_Cell","Target_Pathway")])
results.partial = results.partial %>%
  dplyr::filter(!Target.Collection %in% c("Ligands","Mast")) %>%
  dplyr::filter(!str_detect(Target.Identifier, c("Jha|insilico|PAMP|Icatibant|CST14|SP|x2_activated_inhibition_late|x2_general_inhibition_late|x2_activation_late"))) %>%
  mutate(Criteria.Identifier = case_when(Type == "bulk" ~ Criteria.Identifier, .default = paste0("adj ",Criteria.Identifier)))

results.wide = reshape2::dcast(results.partial, Target.Identifier ~ Criteria.Identifier, value.var = "metricValue")
rownames(results.wide) = results.wide$Target.Identifier
results.wide = results.wide[,-1]

clusterTargets(results.wide, "bulkAndAdj_allPathways_allCells")


# 3. All cells and pathways PCAs, bulk
results.partial = do.call(rbind, Results[c("Target_Cell_PCA", "Target_Pathway_PCA")])
results.partial = results.partial %>%
  dplyr::filter(!Target.Collection %in% c("Ligands","Mast")) %>%
  dplyr::filter(!str_detect(Target.Identifier, c("Jha|insilico|PAMP|Icatibant|CST14|SP|x2_activated_inhibition_late|x2_general_inhibition_late|x2_activation_late"))) %>%
  dplyr::filter(Type == "bulk")

results.wide = reshape2::dcast(results.partial, Target.Identifier ~ Criteria.Identifier, value.var = "metricValue")
rownames(results.wide) = results.wide$Target.Identifier
results.wide = results.wide[,-c(1:4)] # also remove correlations to adj pathway pca

clusterTargets(results.wide, "bulk_allPCAs")


# 4. All pathways PCAs, bulk and adjusted (bulk with bulk, adj with adj)
results.partial = do.call(rbind, Results[c("Target_Pathway_PCA")])
results.partial = results.partial %>%
  dplyr::filter(!Target.Collection %in% c("Ligands","Mast")) %>%
  dplyr::filter(!str_detect(Target.Identifier, c("Jha|insilico|PAMP|Icatibant|CST14|SP|x2_activated_inhibition_late|x2_general_inhibition_late|x2_activation_late")))
  
results.partial = results.partial[-which(results.partial$Type == "bulk" & 
                                         results.partial$Criteria.Identifier %in% c("adj_pathway_meta_pc1","adj_pathway_meta_pc2","adj_pathway_meta_pc3")),]
results.partial = results.partial[-which(results.partial$Type == "adjusted__1__1" & 
                                           results.partial$Criteria.Identifier %in% c("pathway_meta_pc1","pathway_meta_pc2","pathway_meta_pc3")),]
results.wide = reshape2::dcast(results.partial, Target.Identifier ~ Criteria.Identifier, value.var = "metricValue")
rownames(results.wide) = results.wide$Target.Identifier
results.wide = results.wide[,-1]

clusterTargets(results.wide, "bulkAndAdj_matchingPathwayPCAs")
