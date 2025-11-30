devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline) # make sure to load version >= 1.0.1
library(igraph)
library(tidyverse)
library(reshape2)


## Prep network
## ==================
# Pruned disease network
nw_bulk_pruned_wf <- "wf-06da2e64ae"
results <-  get_workflow_outputs(nw_bulk_pruned_wf)
nw_bulk_pruned <- readRDS(results['igraph_NW_sub_0.6.rds',1])

# To use our signatures, we first want to map its connecting nodes from the disease network
signatures <- readRDS(get_workflow_outputs("wf-b1950a97bd"))
signatures = unlist(signatures, recursive = F)
# signatures_connected = visualize_mono_gene_set_with_connection(nw_bulk_pruned, signatures, centrality_df = NULL, 
#                                                                centrality_measure = "page_rank", create_visualization = F)
# 
# signatures_connected <- list("X2"=lapply(signatures_connected, function(cur_sig) {
#   cur_sig[['nodes']][['id']]
# }))

fullList <-  unlist(signatures, recursive = FALSE)
names(fullList) = str_remove(names(fullList),"X2.")
pushToCC(fullList, tagsToPass = list(list(name="object",value="connected_signatures")))
# wf-665e58d115 - connected signatures
# wf-56c3f1b856

nw = igraph::as_data_frame(nw_bulk_pruned)
EntrezInput = unique(as.character(c(unlist(fullList))))
nIter <- 1000 # updating number of iterations


## Criteria: Network Centrality
## ====================================
bulk_centrality_wf = "wf-d535848bea" # calculated in analysis-p03-poc1/notes/multispecifics-generation/9.2. Network Criteria.Rmd
adj_centrality_wf <- "wf-0c39522029" # this is a modified version of the original wf-f254ff53ca, where we have other symbols added. See analysis-p03-poc1/notes/API/new_interface/notebooks/Catalog_project_specific_model.Rmd

# Bulk
# ---------
cytoreason.cc.client::run_method_dist(
  method = "NWParam_VersusRandom_modified",  
  ns = "analysis.p03.POC1",
  NW_WF = bulk_centrality_wf, 
  Signaturelist = fullList,
  BackgroundGeneInput = EntrezInput,
  memory_request = "40Gi",                                  
  Image= "eu.gcr.io/cytoreason/ci-analysis.p03.poc1-package:shiran_API_latest", 
  image="eu.gcr.io/cytoreason/ci-analysis.p03.poc1-package:shiran_API_latest",
  Niter = nIter)
# wf-5dcaa513eb
# wf-db46415c66 - using signatures that were not connected


# using the model teams version (yasmin and matan)
cytoreason.cc.client::run_method_dist(
  method = "get_as_graph",  
  ns = "cytoreason.individual.variation",
  net = cytoreason.assets::read_asset("wf-a7057224d0"), 
  nodes_columns = c("var1", "var2"),
  weight_col = "cor",
  is_correlation = TRUE,
  force_abs = TRUE,
  memory_request = "40Gi",                                  
  image="eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest",
  force_execution = FALSE, 
  replace_image_tags = TRUE, 
  tags = list(list(name = "analysis", value = "compute_measures"),
              list(name = "network", value = "bulk"))
)
# wf-ebc6666d3c

run_function_dist(FUN = function(net){
  library(cytoreason.individual.variation)

  return(
    compute_measures(cytoreason.individual.variation::get_as_graph(net, nodes_columns = c("var1", "var2"), 
                                                                   weight_col = "cor", 
                                                                   is_correlation = TRUE, 
                                                                   force_abs = TRUE))
  )
  
},
net = cytoreason.assets::read_asset("wf-a7057224d0"), 
image = "eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest",
force_execution = FALSE, 
replace_image_tags = TRUE, 
tags = list(list(name = "analysis", value = "compute_measures"),
            list(name = "network", value = "bulk")))
# wf-dc0a446d7c



# Adjusted
# -----------
cytoreason.cc.client::run_method_dist(
  method = "NWParam_VersusRandom_modified",  
  ns = "analysis.p03.POC1",
  NW_WF = adj_centrality_wf, 
  Signaturelist = fullList,
  BackgroundGeneInput = EntrezInput,
  Image= "eu.gcr.io/cytoreason/ci-analysis.p03.poc1-package:shiran_API_latest", 
  image="eu.gcr.io/cytoreason/ci-analysis.p03.poc1-package:shiran_API_latest",
  memory_request = "40Gi",
  Niter = nIter)
# wf-16d659d4d4
# wf-d22c81bd0d - new entrez input


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
    mutate(pval = ifelse(pval == 0, 1/(nIter+1), pval)) %>%
    mutate(Criteria = ifelse(Criteria == "OverlapPercent", "overlap", Criteria)) %>%
    separate(ListName, into = c("Target.Collection", "Target.Identifier"), sep = "\\.", extra = "merge", remove = FALSE) %>%
    rename(Target.ID = ListName, Criteria.Identifier = Criteria, metricType = Measure, metricValue = value) %>%
    mutate(DataType = "Network_Centrality") %>%
    dplyr::select(-ListType)
  
  return(merged)
}

centrality = rbind(cbind(processing("wf-db46415c66", nIter), Type = "bulk"),
                   cbind(processing("wf-16d659d4d4", nIter), Type = "adjusted"))



## Criteria: Topology
## ==================================
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
# wf-da5c96a744

# Combine
# ------------
