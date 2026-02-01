### In this part we cluster the signatures pathway correlations
### with the purpose of grouping them into (hopefully) functional groups -
### and we characterize each group by the most significant pathways they target

# 1. Libraries and data loading
# ------------------------------------------
devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(ComplexHeatmap)
library(tidyverse)

Results = readRDS(get_workflow_outputs("wf-64470d2a55"))
signatures = readRDS(get_workflow_outputs("wf-30bcb0f30c"))

clusterTargets = function(res, plotSuffix, path) {
  library(ggdendro)
  clust = hclust(dist(res))

  # Make the dendrogram prettier
  dend = dendro_data(as.dendrogram(clust), type = "rectangle")
  d = ggdendrogram(dend, theme_dendro = T)
  ggsave(plot = d, paste0(path,"/clusteringDendrogram_",plotSuffix,".png"), width = 10, height = 4, bg = "white")
  
  # optional - vertical
  d = ggdendrogram(dend, theme_dendro = T, rotate = T)
  ggsave(plot = d, paste0(path,"/clusteringDendrogram_vertical_",plotSuffix,".png"), width = 4, height = 10, bg = "white")
  
}

changeLeavesAndPlot = function(res, plotSuffix, path, leaves1, leaves2){
  library(dendextend)
  
  clust = hclust(dist(res))
  dend <- as.dendrogram(clust)
  ord  <- labels(dend)
  leaves1 = which(ord %in% leaves1)
  leaves2 = which(ord %in% leaves2)
  rest = 1:length(ord)
  rest = rest[-which(rest %in% c(leaves1, leaves2))]
  
  new_ord <- c(ord[leaves1], ord[leaves2], ord[rest])

  dend_swapped <- rotate(dend, order = new_ord)
  d = ggdendrogram(dend_swapped, theme_dendro = T)
  ggsave(plot = d, paste0(path,"/clusteringDendrogram_",plotSuffix,".png"), width = 10, height = 4, bg = "white")
}

set.seed(1234)

## All cells and pathways, bulk & adjusted
## ===============================================
results.partial = do.call(rbind, Results[c("Target_Cell","Target_Pathway")])
results.partial = results.partial %>%
  dplyr::filter(Target.ID %in% signatures) %>%
  dplyr::filter(!str_detect(Target.ID,"TGFB2|BMP")) %>%
  mutate(Criteria.Identifier = case_when(Type == "bulk" ~ Criteria.Identifier, .default = paste0("adj ",Criteria.Identifier)))

results.wide = reshape2::dcast(results.partial, Target.Identifier ~ Criteria.Identifier, value.var = "metricValue")
rownames(results.wide) = results.wide$Target.Identifier
results.wide = results.wide[,-1]

# clusterTargets(results.wide, "bulkAndAdj_allPathways_allCells", "~/analysis-s05/figures/Results/Clustering/allCellsPathways/ReducedNumberOfTargetsNoNC/")
changeLeavesAndPlot(res = results.wide, 
                    plotSuffix = "bulkAndAdj_allPathways_allCells", 
                    path = "~/analysis-s05/figures/Results/Clustering/allCellsPathways/ReducedNumberOfTargetsNoNC/",
                    leaves1 = c("random", "smoothedRandom_bottom50", "smoothedRandom_top50"),
                    leaves2 = rownames(results.wide)[-which(rownames(results.wide) %in% c("random", "smoothedRandom_bottom50", "smoothedRandom_top50", "NGF", "IL22", "aIgE late activation refined", "X2 early activation refined"))])

clust = hclust(dist(results.wide))
clust = cutree(clust, k=4) %>% # first split to 6, then join a few
  enframe(name = "signature", value = "cluster") %>%
  arrange(cluster) %>%
  mutate(new_cluster = cluster)

# reorder number
clust$new_cluster[which(clust$cluster == 1)] <- 3
clust$new_cluster[which(clust$cluster == 2)] <- 1
clust$new_cluster[which(clust$cluster == 3)] <- 2

clust = arrange(clust, new_cluster)
pushToCC(clust, tagsToPass = list(list(name="object",value="clustering_bulkadj_allcellsandpathways")))
# wf-9af5bfe1f8
# wf-99bb4eb731
# wf-3fb25d04de - reducing number of targets
# wf-42b768e5ef - reduced number of targets, no negative controls (with randoms)


## Pathways driving the differences between clusters
## ===========================================================
# In this part, we perform a PCA on every pair of clusters, and extract the top 50
# pathways from PC1 and PC2. Expanding this to all cluster pairs we will get a pool of the
# pathways that drive the differences between the clusters
clusterColor = c("1" = "#8B0000", "2" = "#228B22", "3" = "#482870", "4" = "#4169E1", "5" = "#c4007c")

# Extract loadings
# -------------------------
Results = readRDS(get_workflow_outputs("wf-64470d2a55")) # derived from the CCM analysis
signatures = readRDS(get_workflow_outputs("wf-30bcb0f30c"))
  signatures = signatures[-which(signatures == "Itch:Nattkemper")]
  signatures = signatures[-which(str_detect(signatures,"TGFB2|BMP"))]

Results_filtered = Results$Target_Pathway %>%
  dplyr::filter(Target.ID %in% signatures)  %>%
  dplyr::filter(Criteria.Collection %in% c("reactome","h","kegg","Itch","Mast","th2")) %>%
  dplyr::filter(Type == "bulk")

results.wide = reshape2::dcast(Results_filtered, Target.Identifier ~ Criteria.Identifier, value.var = "metricValue")
rownames(results.wide) = results.wide$Target.Identifier
results.wide = results.wide[,-1]
pushToCC(results.wide)
# wf-45574cc2a7
# wf-90c82311a0 - no BTM
# wf-ac0d855526 - reduced number of targets
# wf-2df6be7cc8 - reduced number of targets, no negative controls
# wf-a4cf482684 - reduced number of targets, no negative controls, with randoms

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
clust = readRDS(get_workflow_outputs("wf-43ccf9d07d"))
loadings = extractTopLoadings(nClusters = length(unique(clust$new_cluster)), 
                              results.wide = results.wide, 
                              clusterTable = clust,
                              directory = "~/analysis-s05/figures/Results/Clustering/allCellsPathways/ReducedNumberOfTargetsNoNC/")
pushToCC(loadings, tagsToPass = list(list(name="object",value="clustering_bulkandadj_cellsandpathways"),list(name="analysis",value="loadings")))
# wf-fa2013cbb8
# wf-6251213c81
# wf-3d1796d991 - top 100, no BTM
# wf-9d626a0fba - reduced number of targets - 5 clusters
# wf-19e0a44ad8 - reduced number of targets - no negative controls
# wf-7f9bd9ed3a - reduced number of targets - no negative controls, with randoms

pushToCC(toploadings, tagsToPass = list(list(name="object",value="clustering_bulkandadj_cellsandpathways"),list(name="analysis",value="toploadings")))
# wf-f2617daae8
# wf-04a65ab292
# wf-7143ae3ba8 - top 100, no BTM
# wf-d186aba404 - reduced number of targets - 5 clusters
# wf-cdf93bcad1 - reduced number of targets - no negative controls
# wf-cf96ac6fbf - reduced number of targets - no negative controls, with randoms

# 3.2. Create a heatmap based on the top differentiating pathways
# -------------------------------------------------------------------
library(ComplexHeatmap)

toploadings = readRDS(get_workflow_outputs("wf-cf96ac6fbf"))
clusters = readRDS(get_workflow_outputs("wf-42b768e5ef")) %>% 
  .[-which(.[,"signature"] == "Nattkemper"),] %>% # removing Nattkemper because we also want to correlate to it
  .[!str_detect(.$signature,"TGFB2|BMP"),]

importantClusters = do.call(rbind, toploadings) %>%
  dplyr::select(Pathway, Collection) %>%
  mutate(Pathway = str_remove(Pathway, " \\(reactome\\)| \\(kegg\\)| \\(h\\)| \\(btm\\)| \\(itch\\)| \\(th2\\)")) %>%
  mutate(Pathway = paste0(Collection,":",Pathway))
importantClusters = unique(importantClusters$Pathway)
importantClusters = str_remove(importantClusters," \\(Mast\\)")

# scaling the correlations per pathway for easier interpertation
results.wide = readRDS(get_workflow_outputs("wf-a4cf482684"))
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
# wf-51a4aaeea7 - reduced number of targets, no NCs k=5
# wf-12b7798e2e - reduced number of targets, no NCs, with randoms k=6

pdf("~/analysis-s05/figures/Results/Clustering/allCellsPathways/ReducedNumberOfTargetsNoNC/PathwayClusters.pdf", height = 10, width = 8)
Heatmap(mat, show_row_names = F, show_column_names = T, use_raster = F,
        top_annotation = ann, cluster_columns = T, name = "z-scaled\ncorrelation",
        # row_split = data.frame(pathwayClusters$Cluster),
        column_title = "bulk & adjusted correlation to all cells and pathways")
dev.off()

pdf("~/analysis-s05/figures/Results/Clustering/allCellsPathways/ReducedNumberOfTargetsNoNC/PathwayClusters_top.pdf", height = 10, width = 8)
Heatmap(mat, show_row_names = F, show_column_names = T, use_raster = F,
        top_annotation = ann, cluster_columns = T, name = "z-scaled\ncorrelation",
        # row_split = data.frame(pathwayClusters$Cluster), 
        column_dend_side = "bottom", column_names_side = "top")
dev.off()

pdf("~/analysis-s05/figures/Results/Clustering/allCellsPathways/ReducedNumberOfTargetsNoNC/PathwayClusters_detailed.pdf", height = 100, width = 15)
Heatmap(mat, show_row_names = T, show_column_names = T, use_raster = F, row_names_gp = gpar(fontsize = 10),
        # row_split = data.frame(pathwayClusters$Cluster),
        top_annotation = ann, cluster_columns = T, name = "z-scaled\ncorrelation")
dev.off()

# openxlsx::write.xlsx(pathwayClusters, "~/analysis-s05/figures/Results/Clustering/PCAs/pathwayClusters.xlsx")
# openxlsx::write.xlsx(pathwayClusters, "~/analysis-s05/figures/Results/Clustering/allCellsPathways/pathwayClusters.xlsx")
openxlsx::write.xlsx(pathwayClusters, "~/analysis-s05/figures/Results/Clustering/allCellsPathways/ReducedNumberOfTargets/pathwayClusters_6.xlsx")
write.csv(mat, "~/exportedFiles/clustering_mat.csv", row.names = T)
