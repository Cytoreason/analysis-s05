devtools::load_all("~/analysis-s05/R/utils.R")
devtools::load_all("~/analysis-s05/R/network_functions.R")
library(cytoreason.ccm.pipeline) # make sure to load version >= 1.0.1
library(igraph)
library(tidyverse)
library(reshape2)


## Prep network
## ==================
signatureMapping = readRDS(get_workflow_outputs("wf-fcb1ef2bbd"))

# signatures
signatures <- readRDS(get_workflow_outputs("wf-b1950a97bd"))
fullList <-  unlist(signatures, recursive = FALSE)
names(fullList) = str_remove(names(fullList),"X2.")
pushToCC(fullList, tagsToPass = list(list(name="object",value="signatures")))
# wf-56c3f1b856

signatures_connected = visualize_mono_gene_set_with_connection(nw_bulk_pruned, signatures, centrality_df = NULL, 
                                                               centrality_measure = "page_rank", create_visualization = F)

signatures_connected <- list("X2"=lapply(signatures_connected, function(cur_sig) {
  cur_sig[['nodes']][['id']]
}))
fullList <-  unlist(signatures_connected, recursive = FALSE)
names(fullList) = str_remove(names(fullList),"X2.")
pushToCC(fullList, tagsToPass = list(list(name="object",value="connected_signatures")))
# wf-665e58d115 - connected signatures

# Pruned disease network
nw_bulk_pruned <- cytoreason.assets::read_asset("wf-a7057224d0") # pruned disease network (wf-06da2e64ae), in data.frame form
EntrezInput = unique(as.character(c(unlist(fullList))))
nIter <- 1000 # updating number of iterations


## Criteria: Network Centrality
## ====================================
## First we calculate centrality using the new productized code
## Then we calculate the criteria on the signatures

# Bulk
# ========
# centrality measures
# -----------------------
cytoreason.cc.client::run_method_dist(
  method = "get_as_graph",  
  ns = "cytoreason.individual.variation",
  net = cytoreason.assets::read_asset("wf-a7057224d0"), 
  nodes_columns = c("var1", "var2"),
  weight_col = "cor",
  is_correlation = TRUE,
  force_abs = TRUE,
  memory_request = "2Gi",                                  
  image="eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest",
  force_execution = FALSE, 
  replace_image_tags = TRUE, 
  tags = list(list(name = "analysis", value = "get_as_graph"),
              list(name = "network", value = "bulk")))
# wf-ebc6666d3c

cytoreason.cc.client::run_method_dist(
  method = "compute_measures",
  graph = cytoreason.assets::read_asset("wf-ebc6666d3c"),
  ns = "cytoreason.individual.variation",
  image="eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest",
  tags = list(list(name = "analysis", value = "compute_measures"),
              list(name = "network", value = "bulk")))
# wf-7054c00edf


# criteria
# ------------
cytoreason.cc.client::run_method_dist(
  method = "NWParam_VersusRandom_modified",  
  ns = "analysis.p03.POC1",
  NW_WF = "wf-7054c00edf", 
  Signaturelist = fullList,
  BackgroundGeneInput = EntrezInput,
  Niter = nIter,
  Image= "eu.gcr.io/cytoreason/ci-analysis.p03.poc1-package:shiran_API_latest", 
  image="eu.gcr.io/cytoreason/ci-analysis.p03.poc1-package:shiran_API_latest",
  memory_request = "40Gi")
# wf-bf70438e87


# Adjusted
# =============
# centrality measures
# -----------------------
cytoreason.cc.client::run_method_dist(
  method = "get_as_graph",  
  ns = "cytoreason.individual.variation",
  net = cytoreason.assets::read_asset("wf-aeeeac9a85"), 
  nodes_columns = c("var1", "var2"),
  weight_col = "cor",
  is_correlation = TRUE,
  force_abs = TRUE,
  memory_request = "2Gi",                                  
  image="eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest",
  force_execution = FALSE, 
  replace_image_tags = TRUE, 
  tags = list(list(name = "analysis", value = "get_as_graph"),
              list(name = "network", value = "adjusted")))
# wf-b66b8711ca

cytoreason.cc.client::run_method_dist(
  method = "compute_measures",
  graph = cytoreason.assets::read_asset("wf-b66b8711ca"),
  ns = "cytoreason.individual.variation",
  image="eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest",
  tags = list(list(name = "analysis", value = "compute_measures"),
              list(name = "network", value = "adjusted")))
# wf-81e92ec6a2


# criteria
# ------------
cytoreason.cc.client::run_method_dist(
  method = "NWParam_VersusRandom_modified",  
  ns = "analysis.p03.POC1",
  NW_WF = "wf-81e92ec6a2", 
  Signaturelist = fullList,
  BackgroundGeneInput = EntrezInput,
  Niter = nIter,
  Image= "eu.gcr.io/cytoreason/ci-analysis.p03.poc1-package:shiran_API_latest", 
  image="eu.gcr.io/cytoreason/ci-analysis.p03.poc1-package:shiran_API_latest",
  memory_request = "40Gi")
# wf-eced63000b


# Combine
# ---------------
processing <- function(wfid, nIter) {
  asset <- readRDS(get_workflow_outputs(wfid))
  
  pval <- asset$Significance
  pval$Criteria <- paste(pval$Criteria, pval$Measure, sep = "_")
  pval$Criteria[pval$Criteria %in% "OverlapProbability_NA"] <- "OverlapPercent"
  pval$Measure[pval$Criteria %in% "OverlapPercent"] <- "percent"
  
  values <- asset$FullParam
  values <- values[values$ListType %in% "gene_list", ]
  values_long <- reshape2::melt(values, id.vars = c("ListType", "ListName"), variable.name = "Criteria")
  
  merged <- merge(pval, values_long, by = c("ListName", "Criteria"), all.x = TRUE)
  
  merged = merged %>%
    mutate(pvalue = ifelse(pval == 0, 1/(nIter+1), pval)) %>%
    mutate(Criteria = ifelse(Criteria == "OverlapPercent", "overlap", Criteria)) %>%
    separate(ListName, into = c("Target.Collection", "Target.Identifier"), sep = "\\.", extra = "merge", remove = FALSE) %>%
    rename(Target.ID = ListName, Criteria.Identifier = Criteria, metricType = Measure, metricValue = value) %>%
    mutate(DataType = ifelse(Criteria.Identifier == "overlap","Network_Topology","Network_Centrality")) %>%
    dplyr::select(-ListType,-pval) %>%
    mutate(Criteria.Collection = ifelse(Criteria.Identifier == "overlap","Network_Topology","Network_Centrality"))
  
  merged$Target.Identifier[is.na(merged$Target.Identifier)] <- merged$Target.ID[is.na(merged$Target.Identifier)] 
  
  return(merged)
}

centrality = rbind(cbind(processing("wf-bf70438e87", nIter), Type = "bulk"),
                   cbind(processing("wf-eced63000b", nIter), Type = "adjusted"))

for(id in c("random","top50","bottom50")){
  idx = which(str_detect(centrality$Target.Identifier,id))
  if(id != "random") { id = paste0("smoothedRandom_",id)}
  
  tmp = centrality[idx,] %>%
    group_by(DataType, Type, Criteria.Identifier, Criteria.Collection, metricType, Target.Collection) %>%
    summarise(pvalue = mean(pvalue),
              metricValue = mean(metricValue)) %>%
    mutate(Target.Identifier = id, Target.ID = paste0("negativeControls:",id), hit = NA) %>%
    select(colnames(centrality))
  
  centrality = centrality[-idx,]
  centrality = rbind(centrality, tmp)
  rm(idx, tmp)
}

centrality$Target.Collection = signatureMapping$collection[match(centrality$Target.Identifier, signatureMapping$signature)]
centrality$Target.ID <- signatureMapping$ID[match(centrality$Target.ID, signatureMapping$signature)]
centrality$Target.ID[is.na(centrality$Target.ID)] <- signatureMapping$ID[match(centrality$Target.Identifier[is.na(centrality$Target.ID)], signatureMapping$signature)]
pushToCC(centrality, tagsToPass = list(list(name="object", value="centrality")))
# wf-9393f0f8dc
# wf-79da938ff2 - better signature names
# wf-72d490c885 - even better
uploadToBQ(centrality, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Network")



## Criteria: Topology
## ==================================
# To calculate density we use the connected gene list
fullList = readRDS(get_workflow_outputs("wf-665e58d115"))

# Bulk
# ---------
mapply_dist(geneset_network_measures_random_sample_parallel,
            geneset_name = names(fullList),
            MoreArgs = list(NW_WF = "wf-a7057224d0", # pruned bulk network
                            n = 1000, 
                            genesets_list = fullList,
                            geneset_network_measures_fun = geneset_network_measures),
            image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:master_latest",
            force_execution=F,
            memory_request = '2Gi')
# wf-2b5fabdc5e


# Adjusted
# -----------
mapply_dist(geneset_network_measures_random_sample_parallel,
            geneset_name = names(fullList),
            MoreArgs = list(NW_WF = "wf-aeeeac9a85", # pruned adjusted network
                            n = 1000, 
                            genesets_list = fullList,
                            geneset_network_measures_fun = geneset_network_measures),
            image="eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:master_latest",
            force_execution=F,
            memory_request = '2Gi')
# wf-7b77bf346a


# Combine
# ------------
bulk_topology = cytoreason.assets::read_asset("wf-2b5fabdc5e")
bulk_topology <- bind_rows(bulk_topology)
bulk_topology <-  split(bulk_topology, bulk_topology$ListName)
pushToCC(bulk_topology, tagsToPass = list(list(name="object", value="bulk_topology"))) # includes all randoms
# wf-43b72160d0

bulk_topology_p = lapply(bulk_topology, CollapseToParam)
bulk_topology_p = bind_rows(bulk_topology_p)

density_bulk = bind_rows(bulk_topology) %>%
  dplyr::filter(geneset == "gene_list") %>%
  dplyr::select(ListName, density) %>%
  merge(., bulk_topology_p, by = "ListName") %>%
  dplyr::filter(Criteria == "density") %>%
  mutate(Type = "bulk") %>%
  rename(Target.ID = ListName, Criteria.Identifier = Criteria, metricValue = density, pvalue = pval)
 

adj_topology = cytoreason.assets::read_asset("wf-7b77bf346a")
adj_topology <- bind_rows(adj_topology)
adj_topology <-  split(adj_topology, adj_topology$ListName)
pushToCC(adj_topology, tagsToPass = list(list(name="object", value="adj_topology"))) # includes all randoms
# wf-2e05d11e63

adj_topology_p = lapply(adj_topology, CollapseToParam)
adj_topology_p = bind_rows(adj_topology_p)

density_adj = bind_rows(adj_topology) %>%
  dplyr::filter(geneset == "gene_list") %>%
  dplyr::select(ListName, density) %>%
  merge(., adj_topology_p, by = "ListName") %>%
  dplyr::filter(Criteria == "density") %>%
  mutate(Type = "adjusted") %>%
  rename(Target.ID = ListName, Criteria.Identifier = Criteria, metricValue = density, pvalue = pval)

topology = rbind(density_bulk, density_adj)
topology = topology %>%
  mutate(Target.ID = str_replace(Target.ID,"\\.",":")) %>%
  separate(Target.ID, into = c("Target.Collection", "Target.Identifier"), sep = "\\:", extra = "merge", remove = FALSE) %>%
  mutate(metricType = "density",
         DataType = "Network_Topology",
         Criteria.Collection = "Network_Topology")
  
for(id in c("random","top50","bottom50")){
  idx = which(str_detect(topology$Target.Identifier,id))
  if(id != "random") { id = paste0("smoothedRandom_",id)}
  
  tmp = topology[idx,] %>%
    group_by(DataType, Type, Criteria.Identifier, Criteria.Collection, metricType, Target.Collection) %>%
    summarise(pvalue = mean(pvalue),
              metricValue = mean(metricValue)) %>%
    mutate(Target.Identifier = id, Target.ID = paste0("negativeControls:",id), hit = NA) %>%
    select(colnames(topology))
  
  topology = topology[-idx,]
  topology = rbind(topology, tmp)
  rm(idx, tmp)
}

topology$Target.Collection = signatureMapping$collection[match(topology$Target.Identifier, signatureMapping$signature)]
topology$Target.ID = signatureMapping$ID[match(topology$Target.ID, signatureMapping$previousID)]

pushToCC(topology, tagsToPass = list(list(name="object", value="topology")))
# wf-e6682f3eaf
# wf-15e8917643 - better signature names
# wf-8556686d74 - even better
uploadToBQ(topology, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Network", disposition = "WRITE_APPEND")



## Visualization
## ==========================
# centrality
# -------------------
centrality_randoms = rbind(cbind(readRDS(get_workflow_outputs("wf-eced63000b"))[["FullParam"]], Type = "adjusted"),
                           cbind(readRDS(get_workflow_outputs("wf-bf70438e87"))[["FullParam"]], Type = "bulk"))

# Summarize random and shuffled random
for(id in c("random","top50","bottom50")){
  idx = which(str_detect(centrality_randoms$ListName,id))
  if(id != "random") { id = paste0("smoothedRandom_",id)}
  
  tmp = centrality_randoms[idx,] %>%
    group_by(ListType, Type) %>%
    summarise(across(2:14, .fns = median)) %>%
    mutate(ListName = paste0("negativeControls:",id)) %>%
    select(colnames(centrality_randoms))
  
  centrality_randoms = centrality_randoms[-idx,]
  centrality_randoms = rbind(centrality_randoms, tmp)
  rm(idx, tmp)
}

centrality_randoms = centrality_randoms[,c("ListType","ListName","Type","eigen_centrality_median")]
centrality_randoms = centrality_randoms %>%
  rename(Target.ID = ListName, metricValue = eigen_centrality_median) %>%
  mutate(Criteria.Identifier = "eigen_centrality_median") %>%
  mutate(Target.Identifier = case_when(!str_detect(Target.ID,"random|Random") ~ signatureMapping$New_identifier[match(Target.ID, signatureMapping$signature)],
                                       str_detect(Target.ID,"random|Random") ~ signatureMapping$New_identifier[match(Target.ID, signatureMapping$ID)])) %>%
  dplyr::filter(Target.Identifier != "") %>%
  mutate(Target.Collection = signatureMapping$collection[match(Target.Identifier, signatureMapping$New_identifier)]) %>%
  dplyr::filter(!Target.Collection %in% c("Ligands","Mast","Itch","Neuronal")) %>%
  dplyr::filter(!str_detect(Target.Identifier, c("insilico|x2_|PAMP|CST|SP|aIgE|Icatibant|X2 late activation|X2 early activated|X2 early activation"))) %>%
  dplyr::filter(!Target.Identifier == "X2 late inhibition") %>%
  mutate(Target.ID = signatureMapping$ID[match(Target.Identifier, signatureMapping$New_identifier)])

centrality_randoms$ListType[which(centrality_randoms$ListType != "gene_list")] = "random"
pushToCC(centrality_randoms, tagsToPass = list(list(name="object",value="centrality")))
# wf-eca695022d

centrality_randoms2 = split(centrality_randoms, centrality_randoms$ListType)
centrality_randoms2$gene_list$Target = cytoreason.gx::reorder_within(centrality_randoms2$gene_list$Target.Identifier, centrality_randoms2$gene_list$metricValue, centrality_randoms2$gene_list$Type)
or = levels(centrality_randoms2$gene_list$Target)
centrality_randoms2$random$Target <- factor(
  paste0(centrality_randoms2$random$Target.Identifier, "___", centrality_randoms2$random$Type),
  levels = or
)

bulk = list(gene_list = centrality_randoms2$gene_list[which(centrality_randoms2$gene_list$Type == "bulk"),],
            random = centrality_randoms2$random[which(centrality_randoms2$random$Type == "bulk"),])

ggplot(bulk$random, aes(x = metricValue, y = Target)) +
  geom_point(color = "grey", alpha = 0.3, size = 3) +
  geom_point(data = bulk$gene_list, inherit.aes = T, aes(color = Target.ID), size=6) +
  cytoreason.gx::scale_y_reordered()+
  scale_color_manual(values = targetColors) +
  scale_x_continuous(transform = squash_axis(0.005,0.04,10))+
  theme_minimal() +
  theme(legend.position = "none") +
  ggpubr::border()+
  labs(x = "Median Eigen Centrality", y = NULL)
ggsave("~/analysis-s05/figures/Results/network_eigenCentrality.png", width = 6, height = 6, bg= "white")


# topology
# -------------------
topology_randoms = readRDS(get_workflow_outputs("wf-43b72160d0"))
topology_randoms = bind_rows(topology_randoms)
topology_randoms = topology_randoms[,c("ListType", "ListName", "density")]

# Summarize random and shuffled random
for(id in c("random","top50","bottom50")){
  idx = which(str_detect(topology_randoms$ListName,id))
  if(id != "random") { id = paste0("smoothedRandom_",id)}
  
  tmp = topology_randoms[idx,] %>%
    group_by(ListType) %>%
    summarise(density = median(density)) %>%
    mutate(ListName = paste0("negativeControls:",id)) %>%
    select(colnames(topology_randoms))
  
  topology_randoms = topology_randoms[-idx,]
  topology_randoms = rbind(topology_randoms, tmp)
  rm(idx, tmp)
}

topology_randoms = topology_randoms %>%
  rename(Target.ID = ListName, metricValue = density) %>%
  mutate(Target.ID = str_replace(Target.ID, "X2.","X2:")) %>%
  mutate(Criteria.Identifier = "density") %>%
  mutate(Target.Identifier = signatureMapping$New_identifier[match(Target.ID, signatureMapping$previousID)]) %>%
  dplyr::filter(Target.Identifier != "") %>%
  mutate(Target.Collection = signatureMapping$collection[match(Target.Identifier, signatureMapping$New_identifier)]) %>%
  dplyr::filter(!Target.Collection %in% c("Ligands","Mast","Itch","Neuronal")) %>%
  dplyr::filter(!str_detect(Target.Identifier, c("insilico|x2_|PAMP|CST|SP|aIgE|Icatibant|X2 late activation|X2 early activated|X2 early activation"))) %>%
  dplyr::filter(!Target.Identifier == "X2 late inhibition") %>%
  mutate(Target.ID = signatureMapping$ID[match(Target.Identifier, signatureMapping$New_identifier)]) %>%
  mutate(Type = "bulk")

topology_randoms$ListType[which(topology_randoms$ListType != "gene_list")] = "random"
pushToCC(topology_randoms, tagsToPass = list(list(name="object",value="topology")))
# wf-bf87b52ff4

topology_randoms2 = split(topology_randoms, topology_randoms$ListType)
topology_randoms2$gene_list$Target = cytoreason.gx::reorder_within(topology_randoms2$gene_list$Target.Identifier, topology_randoms2$gene_list$metricValue, topology_randoms2$gene_list$Type)
or = levels(topology_randoms2$gene_list$Target)
topology_randoms2$random$Target <- factor(
  paste0(topology_randoms2$random$Target.Identifier, "___", topology_randoms2$random$Type),
  levels = or
)

ggplot(topology_randoms2$random, aes(x = metricValue, y = Target)) +
  geom_point(color = "grey", alpha = 0.3, size = 3) +
  geom_point(data = topology_randoms2$gene_list, inherit.aes = T, aes(color = Target.ID), size=6) +
  cytoreason.gx::scale_y_reordered()+
  scale_color_manual(values = targetColors) +
  # scale_x_continuous(transform = squash_axis(0.005,0.04,10))+
  theme_minimal() +
  theme(legend.position = "none") +
  ggpubr::border()+
  labs(x = "Density", y = NULL)
ggsave("~/analysis-s05/figures/Results/network_density.png", width = 6, height = 6, bg= "white")



### Visualization of Adjusted
### ======================================
centrality_randoms=readRDS(get_workflow_outputs("wf-eca695022d"))

centrality_randoms2 = split(centrality_randoms, centrality_randoms$ListType)
centrality_randoms2$gene_list$Target = cytoreason.gx::reorder_within(centrality_randoms2$gene_list$Target.Identifier, centrality_randoms2$gene_list$metricValue, centrality_randoms2$gene_list$Type)
or = levels(centrality_randoms2$gene_list$Target)
centrality_randoms2$random$Target <- factor(
  paste0(centrality_randoms2$random$Target.Identifier, "___", centrality_randoms2$random$Type),
  levels = or
)

adj = list(gene_list = centrality_randoms2$gene_list[which(centrality_randoms2$gene_list$Type == "adjusted"),],
            random = centrality_randoms2$random[which(centrality_randoms2$random$Type == "adjusted"),])

ggplot(adj$random, aes(x = metricValue, y = Target)) +
  geom_point(color = "grey", alpha = 0.3, size = 3) +
  geom_point(data = adj$gene_list, inherit.aes = T, aes(color = Target.ID), size=6) +
  cytoreason.gx::scale_y_reordered()+
  scale_color_manual(values = targetColors) +
  scale_x_continuous(transform = squash_axis(0.005,0.04,10))+
  theme_minimal() +
  theme(legend.position = "none") +
  ggpubr::border()+
  labs(x = "Median Eigen Centrality", y = NULL)
ggsave("~/analysis-s05/figures/Results/network_eigenCentrality_adj.png", width = 6, height = 6, bg= "white")
