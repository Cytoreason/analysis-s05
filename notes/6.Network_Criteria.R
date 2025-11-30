devtools::load_all("~/analysis-s05/R/utils.R")
devtools::load_all("~/analysis-s05/R/network_functions.R")
library(cytoreason.ccm.pipeline) # make sure to load version >= 1.0.1
library(igraph)
library(tidyverse)
library(reshape2)


## Prep network
## ==================
collectionMapping = readRDS(get_workflow_outputs("wf-d3240e7fb6"))

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
  values_long <- melt(values, id.vars = c("ListType", "ListName"), variable.name = "Criteria")
  
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

centrality$Target.Collection = collectionMapping$collection[match(centrality$Target.Identifier, collectionMapping$signature)]
centrality$Target.ID[!str_detect(centrality$Target.ID,"negativeControls")] <- paste0("X2:",centrality$Target.ID[!str_detect(centrality$Target.ID,"negativeControls")])
pushToCC(centrality, tagsToPass = list(list(name="object", value="centrality")))
# wf-9393f0f8dc
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

topology$Target.Collection = collectionMapping$collection[match(topology$Target.Identifier, collectionMapping$signature)]
pushToCC(topology, tagsToPass = list(list(name="object", value="topology")))
# wf-e6682f3eaf
uploadToBQ(topology, bqdataset = "s05_atopic_dermatitis", tableName = "Results_Network", disposition = "WRITE_APPEND")
