### In this part we cluster the signatures pathway correlations
### with the purpose of grouping them into (hopefully) functional groups - 
### and we characterize each group by the most significant pathways they target

# 1. Libraries and data loading
# ------------------------------------------
devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)

Results = read_asset("wf-f0959d74d1")

results.partial = do.call(rbind, Results[c("Target_Cell","Target_Pathway")])
results.partial = results.partial %>%
  dplyr::filter(Type == "bulk") %>%
  dplyr::filter(!str_detect(Target.Identifier, "archs_refined|cov|insilico")) %>%
  mutate(Target.Identifier = str_remove(Target.Identifier, "_50"))

results.wide = reshape2::dcast(results.partial, Target.Identifier ~ Criteria.Identifier, value.var = "metricValue")
rownames(results.wide) = results.wide$Target.Identifier
results.wide = results.wide[,-1]

# 2. Hierarchical clustering using euclidean distance and complete linkage
# --------------------------------------------------------------------------
library(ggdendro)
clust = hclust(dist(results.wide)) # 
plot(clust)

# Make the dendrogram prettier
dend = dendro_data(as.dendrogram(clust), type = "rectangle")
ggdendrogram(dend, theme_dendro = T)
ggsave("~/analysis-s05/figures/clusteringDendrogram.png", width = 12, height = 4, bg = "white")
