### In this part we cluster the signatures pathway correlations
### with the purpose of grouping them into (hopefully) functional groups -
### and we characterize each group by the most significant pathways they target

# 1. Libraries and data loading
# ------------------------------------------
devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(ComplexHeatmap)
library(tidyverse)
library(ggdendro)

Results = readRDS(get_workflow_outputs("wf-64470d2a55"))

clusterTargets = function(res, plotSuffix, path) {
  clust = hclust(dist(res))

  # Make the dendrogram prettier
  dend = dendro_data(as.dendrogram(clust), type = "rectangle")
  d = ggdendrogram(dend, theme_dendro = T)
  ggsave(plot = d, paste0(path,"/clusteringDendrogram_",plotSuffix,".png"), width = 10, height = 4, bg = "white")
  
  # optional - vertical
  d = ggdendrogram(dend, theme_dendro = T, rotate = T)
  ggsave(plot = d, paste0(path,"/clusteringDendrogram_vertical_",plotSuffix,".png"), width = 4, height = 10, bg = "white")
  
}


## Trials
## ===========================
set.seed(1234)

# 1. All cells and pathways, bulk
results.partial = do.call(rbind, Results[c("Target_Cell","Target_Pathway")])
results.partial = results.partial %>%
  dplyr::filter(!Target.Collection %in% c("Ligands","Mast","epidermis","th2")) %>%
  dplyr::filter(!str_detect(Target.ID,"neuroinflammation")) %>%
  dplyr::filter(!str_detect(Target.Identifier, c("Jha|insilico|PAMP|Icatibant|CST14|SP|x2_activated_inhibition_late|x2_general_inhibition_late|x2_activation_late|BMP7|CXCL4|X2 late"))) %>%
  dplyr::filter(!Target.Identifier %in% c("X2 early general inhibition","X2 early activation","X2 early activated inhibition","aIgE late activation")) %>%
  dplyr::filter(Type == "bulk")

signatures = unique(results.partial$Target.ID)
pushToCC(signatures, tagsToPass = list(list(name="object",value="signatures_used_in_clustering")))
# wf-38f6dccf0a
# wf-30bcb0f30c - removed additional signatures per client request

results.wide = reshape2::dcast(results.partial, Target.Identifier ~ Criteria.Identifier, value.var = "metricValue")
rownames(results.wide) = results.wide$Target.Identifier

# 2. All cells and pathways, bulk & adjusted
results.partial = do.call(rbind, Results[c("Target_Cell","Target_Pathway")])
results.partial = results.partial %>%
  dplyr::filter(!Target.Collection %in% c("Ligands","Mast","epidermis","th2")) %>%
  dplyr::filter(!str_detect(Target.ID,"neuroinflammation")) %>%
  dplyr::filter(!str_detect(Target.Identifier, c("Jha|insilico|PAMP|Icatibant|CST14|SP|x2_activated_inhibition_late|x2_general_inhibition_late|x2_activation_late|BMP7|CXCL4|X2 late"))) %>%
  dplyr::filter(!Target.Identifier %in% c("X2 early general inhibition","X2 early activation","X2 early activated inhibition","aIgE late activation")) %>%
  mutate(Criteria.Identifier = case_when(Type == "bulk" ~ Criteria.Identifier, .default = paste0("adj ",Criteria.Identifier)))

results.wide = reshape2::dcast(results.partial, Target.Identifier ~ Criteria.Identifier, value.var = "metricValue")
rownames(results.wide) = results.wide$Target.Identifier
results.wide = results.wide[,-1]

clusterTargets(results.wide, "bulkAndAdj_allPathways_allCells", "~/analysis-s05/figures/Results/Clustering/allCellsPathways/ReducedNumberOfTargets")


# # 3. All cells and pathways PCAs, bulk
# results.partial = do.call(rbind, Results[c("Target_Cell_PCA", "Target_Pathway_PCA")])
# results.partial = results.partial %>%
#   dplyr::filter(!Target.Collection %in% c("Ligands","Mast","epidermis","th2")) %>%
#   dplyr::filter(!str_detect(Target.ID,"neuroinflammation")) %>%
#   dplyr::filter(!str_detect(Target.Identifier, c("Jha|insilico|PAMP|Icatibant|CST14|SP|x2_activated_inhibition_late|x2_general_inhibition_late|x2_activation_late"))) %>%
#   dplyr::filter(Type == "bulk")
# 
# results.wide = reshape2::dcast(results.partial, Target.Identifier ~ Criteria.Identifier, value.var = "metricValue")
# rownames(results.wide) = results.wide$Target.Identifier
# results.wide = results.wide[,-c(1:4)] # also remove correlations to adj pathway pca
# 
# clusterTargets(results.wide, "bulk_allPCAs")
# 
# 
# # 4. All pathways PCAs, bulk and adjusted (bulk with bulk, adj with adj)
# results.partial = do.call(rbind, Results[c("Target_Pathway_PCA")])
# results.partial = results.partial %>%
#   dplyr::filter(!Target.Collection %in% c("Ligands","Mast","epidermis","th2")) %>%
#   dplyr::filter(!str_detect(Target.ID,"neuroinflammation")) %>%
#   dplyr::filter(!str_detect(Target.Identifier, c("Jha|insilico|PAMP|Icatibant|CST14|SP|x2_activated_inhibition_late|x2_general_inhibition_late|x2_activation_late")))
# 
# results.partial = results.partial[-which(results.partial$Type == "bulk" &
#                                          results.partial$Criteria.Identifier %in% c("adj_pathway_meta_pc1","adj_pathway_meta_pc2","adj_pathway_meta_pc3")),]
# results.partial = results.partial[-which(results.partial$Type == "adjusted__1__1" &
#                                            results.partial$Criteria.Identifier %in% c("pathway_meta_pc1","pathway_meta_pc2","pathway_meta_pc3")),]
# results.wide = reshape2::dcast(results.partial, Target.Identifier ~ Criteria.Identifier, value.var = "metricValue")
# rownames(results.wide) = results.wide$Target.Identifier
# results.wide = results.wide[,-1]
# 
# clusterTargets(results.wide, "bulkAndAdj_matchingPathwayPCAs")


## Export clusters
## =============================
set.seed(1234)

# 1. All cells and pathways, bulk & adjusted
# -----------------------------------------------
signatures = readRDS(get_workflow_outputs("wf-30bcb0f30c"))
results.partial = do.call(rbind, Results[c("Target_Cell","Target_Pathway")])
results.partial = results.partial %>%
  dplyr::filter(Target.ID %in% signatures) %>%
  mutate(Criteria.Identifier = case_when(Type == "bulk" ~ Criteria.Identifier, .default = paste0("adj ",Criteria.Identifier)))

results.wide = reshape2::dcast(results.partial, Target.Identifier ~ Criteria.Identifier, value.var = "metricValue")
rownames(results.wide) = results.wide$Target.Identifier
results.wide = results.wide[,-1]

clust = hclust(dist(results.wide))
clust = cutree(clust, h=20) %>% # first split to 6, then join a few
  enframe(name = "signature", value = "cluster") %>%
  arrange(cluster) %>%
  mutate(new_cluster = cluster)

# now we re-arragne
# cluster 6 becomes cluster 2
  clust$new_cluster[which(clust$cluster %in% c(1,5))] <- 3
# cluster 2 becomes cluster 4
  clust$new_cluster[which(clust$cluster == 2)] <- 1
# with the negative controls is cluster 5
  clust$new_cluster[which(clust$cluster %in% c(3,6))] <- 4
# remaining cluster is number 3
  clust$new_cluster[which(clust$cluster == 4)] <- 2

  clust = arrange(clust, new_cluster)
  pushToCC(clust, tagsToPass = list(list(name="object",value="clustering_bulkadj_allcellsandpathways")))
  # wf-9af5bfe1f8
  # wf-99bb4eb731
  # wf-3fb25d04de - reducing number of targets

# # 2. All cells and pathways PCAs, bulk
# # -----------------------------------------
# signatures = readRDS(get_workflow_outputs("wf-38f6dccf0a"))
# results.partial = do.call(rbind, Results[c("Target_Cell_PCA", "Target_Pathway_PCA")])
# results.partial = results.partial %>%
#   dplyr::filter(Target.ID %in% signatures) %>%
#   dplyr::filter(Type == "bulk")
# 
# results.wide = reshape2::dcast(results.partial, Target.Identifier ~ Criteria.Identifier, value.var = "metricValue")
# rownames(results.wide) = results.wide$Target.Identifier
# results.wide = results.wide[,-c(1:4)] # also remove correlations to adj pathway pca
# 
# clust = hclust(dist(results.wide))
# clust = cutree(clust, h = 0.63) %>% # first split to 6, then join a few
#   enframe(name = "signature", value = "cluster") %>%
#   arrange(cluster) %>%
#   mutate(new_cluster = cluster)
# 
# # now we re-arragne
# # clusters 3 and 6 become cluster 5
# clust$new_cluster[which(clust$cluster %in% c(3,6))] <- 5
# 
# clust$new_cluster[which(clust$cluster == 2)] <- 3
# clust$new_cluster[which(clust$cluster == 5)] <- 2
# 
# clust = arrange(clust, new_cluster)
# pushToCC(clust, tagsToPass = list(list(name="object",value="clustering_bulk_PCAcellsandpathways")))
# # wf-2eb40e9d34


## Pathways driving the differences between clusters
## -----------------------------------------------------------
# In this part, we perform a PCA on every pair of clusters, and extract the top 50
# pathways from PC1 and PC2. Expanding this to all cluster pairs we will get a pool of the
# pathways that drive the differences between the clusters
clusterColor = c("1" = "#8B0000", "2" = "#228B22", "3" = "#482870", "4" = "#4169E1", "5" = "#c4007c")

# Extract loadings
# -------------------------
Results = readRDS(get_workflow_outputs("wf-64470d2a55")) # derived from the CCM analysis
signatures = readRDS(get_workflow_outputs("wf-30bcb0f30c"))
  signatures = signatures[-which(signatures == "Itch:Nattkemper")]

Results_filtered = Results$Target_Pathway %>%
  dplyr::filter(Target.ID %in% signatures)  %>%
  dplyr::filter(Criteria.Collection %in% c("reactome","h","kegg","Itch","Mast","th2")) %>%
  dplyr::filter(Type == "bulk")

results.wide = reshape2::dcast(Results_filtered, Target.Identifier ~ Criteria.Identifier, value.var = "metricValue")
rownames(results.wide) = results.wide$Target.Identifier
results.wide = results.wide[,-1]
# wf-ffee8c9b74 - without th2
# wf-45574cc2a7
# wf-90c82311a0 - no BTM
# wf-ac0d855526 - reduced number of targets

# The next function saves and plots the *top* loadings from every pair of clusters, and returns all the loadings
# for every pair.
extractTopLoadings = function(nClusters, results.wide, clusterTable, directory) {
  loadings = apply(combn(1:nClusters, m = 2), 2, function(cluster) {
    group1 = cluster[1]
    group2 = cluster[2]
    d = results.wide[which(rownames(results.wide) %in% clusterTable$signature[which(clusterTable$new_cluster %in% c(group1,group2))]),]
    p = prcomp(d)
    
    # plotting the PCA
    pcs = data.frame(Target = rownames(p$x), 
                     Cluster = as.character(clusterTable$new_cluster[match(rownames(p$x),clusterTable$signature)]),
                     p$x)
    
    ggplot(pcs, aes(x = PC1, y = PC2, color = Cluster)) +
      geom_point(size = 3) +
      ggrepel::geom_text_repel(aes(label = Target), show.legend = F) +
      xlab(paste0("PC1 (",round(summary(p)$importance[2,1] * 100,1),"%)")) +
      ylab(paste0("PC2 (",round(summary(p)$importance[2,2] * 100,1),"%)")) +
      scale_color_manual(values = clusterColor) +
      theme_light() + ggpubr::border() +
      labs(color = "Cluster", title = paste0("Comparing cluster ",group1," to cluster ",group2))
    ggsave(paste0(directory,"/c",group1,"c",group2,".pca.png"),
           width = 6, height = 4, bg = "white", scale = 1.2)
    
    # Extracts the loadings for PC1 and PC2
    loadings = data.frame(Pathway = rownames(p$rotation), p$rotation)
    loadings = reshape2::melt(loadings[,1:3], id.vars = 1, variable.name = "PC", value.name = "loading")
    
    # Saving the top 50 pathways for PC1 and PC2 to each direction of the PC
    top = loadings %>%
      separate(Pathway, into = c("Collection", "Pathway"), sep = "\\:", extra = "merge", remove = FALSE) %>%
      mutate(Pathway = paste0(Pathway," (",Collection,")")) %>%
      mutate(Direction = factor(sign(loading), labels = c("Negative Direction","Positive Direction"))) %>%
      group_by(PC) %>%
      top_n(abs(loading), n = 100) %>%
      ungroup
    
    toploadings[[paste(cluster,collapse = "-")]] <<- data.frame(Cluster1 = group1, Cluster2 = group2, top)
    
    # Plotting the first 50 pathways for PC1 and PC2 to each direction
    top$Pathway = str_trunc(top$Pathway, width = 60, side = "right") # for better visualization

    pc1 = top[which(top$PC == "PC1"),]
    pc2 = top[which(top$PC == "PC2"),]
    pc1$Feature = cytoreason.gx::reorder_within(pc1$Pathway, abs(pc1$loading), within = pc1$Direction, sep = "00")
    pc2$Feature = cytoreason.gx::reorder_within(pc2$Pathway, abs(pc2$loading), within = pc2$Direction, sep = "00")

    ggplot(pc1, aes(x = loading, y = Feature)) +
      geom_col() +
      cytoreason.gx::scale_y_reordered(sep = "00")+
      facet_wrap(~Direction, scales = "free") +
      labs(x = "PC loading", y = NULL,
           title = paste0("Top 50 loadings - PC1 - comparing cluster ",group1," to cluster ",group2)) +
      theme_light() + ggpubr::border()
    ggsave(paste0(directory,"/c",group1,"c",group2,".toploadingsPC1.png"), width = 10, height = 7, bg = "white", scale = 1.5)

    ggplot(pc2, aes(x = loading, y = Feature)) +
      geom_col() +
      cytoreason.gx::scale_y_reordered(sep = "00")+
      facet_wrap(~Direction, scales = "free") +
      labs(x = "PC loading", y = NULL,
           title = paste0("Top 50 loadings - PC2 - comparing cluster ",group1," to cluster ",group2)) +
      theme_light() + ggpubr::border()
    ggsave(paste0(directory,"/c",group1,"c",group2,".toploadingsPC2.png"), width = 10, height = 7, bg = "white", scale = 1.5)
    
    # Returning the loading value for all pathways
    loadings = reshape2::dcast(loadings, Pathway ~ PC, value.var = "loading")
    loadings = data.frame(Cluster1 = group1, 
                          Cluster2 = group2,
                          Targets.Cluster1 = paste(clusterTable$signature[which(clusterTable$new_cluster == group1)], collapse = ";"),
                          Targets.Cluster2 = paste(clusterTable$signature[which(clusterTable$new_cluster == group2)], collapse = ";"),
                          loadings)
    return(loadings)
  })
  loadings = do.call(rbind,loadings) # wf-3f957b7686
  return(loadings)
}

# clustering based on all cells and pathways, bulk and adjusted
toploadings = list()
clust = readRDS(get_workflow_outputs("wf-3fb25d04de"))
loadings = extractTopLoadings(nClusters = length(unique(clust$new_cluster)), 
                              results.wide = results.wide, 
                              clusterTable = clust,
                              directory = "~/analysis-s05/figures/Results/Clustering/allCellsPathways/ReducedNumberOfTargets/")
pushToCC(loadings, tagsToPass = list(list(name="object",value="clustering_bulkandadj_cellsandpathways"),list(name="analysis",value="loadings")))
# wf-fa2013cbb8
# wf-6251213c81
# wf-3d1796d991 - top 100, no BTM
# wf-9d626a0fba - reduced number of targets - 5 clusters
pushToCC(toploadings, tagsToPass = list(list(name="object",value="clustering_bulkandadj_cellsandpathways"),list(name="analysis",value="toploadings")))
# wf-f2617daae8
# wf-04a65ab292
# wf-7143ae3ba8 - top 100, no BTM
# wf-d186aba404 - reduced number of targets - 5 clusters


# # clustering based on bulk cell and pathway pcas
# toploadings = list()
# clust = readRDS(get_workflow_outputs("wf-2eb40e9d34"))
# loadings = extractTopLoadings(nClusters = length(unique(clust$new_cluster)), 
#                               results.wide = results.wide, 
#                               clusterTable = clust,
#                               directory = "~/analysis-s05/figures/Results/Clustering/PCAs")
# pushToCC(loadings, tagsToPass = list(list(name="object",value="clustering_bulk_PCAcellsandpathways"),list(name="analysis",value="loadings")))
# # wf-d67327dcf0
# pushToCC(toploadings, tagsToPass = list(list(name="object",value="clustering_bulk_PCAcellsandpathways"),list(name="analysis",value="toploadings")))
# # wf-db33789134



# 3.2. Create a heatmap based on the top differentiating pathways
# -------------------------------------------------------------------
library(ComplexHeatmap)

# bulk & adjusted all cells and pathways
toploadings = readRDS(get_workflow_outputs("wf-d186aba404"))
clusters = readRDS(get_workflow_outputs("wf-3fb25d04de")) %>% .[-which(.[,"signature"] == "Nattkemper"),] # removing Nattkemper because we also want to correlate to it

# # bulk cell and pathway pcas
# toploadings = readRDS(get_workflow_outputs("wf-db33789134"))
# clusters = readRDS(get_workflow_outputs("wf-2eb40e9d34")) %>% .[-27,] # removing Nattkemper because we also want to correlate to it


importantClusters = do.call(rbind, toploadings) %>%
  dplyr::select(Pathway, Collection) %>%
  mutate(Pathway = str_remove(Pathway, " \\(reactome\\)| \\(kegg\\)| \\(h\\)| \\(btm\\)| \\(itch\\)| \\(th2\\)")) %>%
  mutate(Pathway = paste0(Collection,":",Pathway))
importantClusters = unique(importantClusters$Pathway)
importantClusters = str_remove(importantClusters," \\(Mast\\)| \\(Itch\\)")

# scaling the correlations per pathway for easier interpertation
results.wide = readRDS(get_workflow_outputs("wf-ac0d855526"))
mat = apply(results.wide[,which(colnames(results.wide) %in% importantClusters)],2,scale)
ncol(mat) == length(importantClusters) # needs to be TRUE
rownames(mat) = rownames(results.wide)
mat = mat[match(clusters$signature, rownames(mat)),]
mat = t(mat)

# annotate columns by cluster
clusters$Cluster = factor(clusters$new_cluster)
ann = HeatmapAnnotation(df = clusters[,-c(1:3),drop=F], col = list(Cluster = clusterColor))

# split pathways two ways for easier reading of the pathways in each cluster
ordering = hclust(dist(mat))
pathwayClusters = data.frame(Pathway = ordering$labels,
                             Order = ordering$order,
                             Cluster = cutree(ordering, k=6))
pushToCC(pathwayClusters)
# wf-19e9b5f3da
# wf-f8d8fe3f59
# wf-c0fc0ae6fa - reduced number of targets, k=8
# wf-b5fb3e16fa - reduced number of targets, k=7
# wf-c73ebe0581 - reduced number of targets, k=6

pdf("~/analysis-s05/figures/Results/Clustering/allCellsPathways/ReducedNumberOfTargets/PathwayClusters_6", height = 10, width = 8)
Heatmap(mat, show_row_names = F, show_column_names = T, use_raster = F,
        top_annotation = ann, cluster_columns = T, name = "scaled\ncorrelation",
        row_split = data.frame(pathwayClusters$Cluster),
        column_title = "bulk & adjusted correlation to all cells and pathways")
dev.off()


# openxlsx::write.xlsx(pathwayClusters, "~/analysis-s05/figures/Results/Clustering/PCAs/pathwayClusters.xlsx")
# openxlsx::write.xlsx(pathwayClusters, "~/analysis-s05/figures/Results/Clustering/allCellsPathways/pathwayClusters.xlsx")
openxlsx::write.xlsx(pathwayClusters, "~/analysis-s05/figures/Results/Clustering/allCellsPathways/ReducedNumberOfTargets/pathwayClusters_6.xlsx")



# split pathways two ways for easier reading of the pathways in each cluster
ordering = hclust(dist(mat))
pathwayClusters = data.frame(Pathway = ordering$labels,
                             Cluster = cutree(ordering, k = 3),
                             SecondaryClustering = cutree(ordering, k=6))

pdf("~/analysis-s05/figures/Results/Clustering/PCAs/pathwayClusters.pdf", height = 10, width = 8)
Heatmap(mat, show_row_names = F, show_column_names = T, use_raster = F,
        top_annotation = ann, cluster_columns = T, name = "scaled\ncorrelation",
        row_split = data.frame(pathwayClusters$SecondaryClustering),
        column_title = "bulk correlation to cells and pathways PCA")
dev.off()

pathwayClusters = pathwayClusters %>%
  dplyr::select(Pathway,SecondaryClustering) %>%
  rename(Pathway_Cluster = SecondaryClustering) %>%
  mutate(enrichedIn = case_when(Pathway_Cluster == 1 ~ "Target Clusters 2 and 4, but some 1 and 3",
                                Pathway_Cluster == 2 ~ "Target Clusters 2, 3, 4",
                                Pathway_Cluster == 3 ~ "Target Cluster 5",
                                Pathway_Cluster == 4 ~ "Target Cluster 3",
                                Pathway_Cluster == 5 ~ "Target Clusters 4 and 5, a bit of 3"))
# wf-092a9f7fcf
openxlsx::write.xlsx(pathwayClusters, "~/analysis-s05/figures/Results/Clustering/PCAs/pathwayClusters.xlsx")


### PCA for all targets
### ============================
