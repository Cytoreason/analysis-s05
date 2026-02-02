library(cytoreason.ccm.pipeline)
library(tidyverse)
library(igraph)
library(visNetwork)


## 0. Functions
## ======================
get_node_positions <- function(network) {
  layout <- layout_with_fr(network)  # Fruchterman-Reingold layout for better visualization
  positions <- data.frame(id = names(V(network)), x = layout[, 1], y = layout[, 2])
  return(positions)
}

visualize_gene_set_subgraphs <- function(whole_network, gene_sets, centrality_df = NULL, centrality_measure = "betweenness") {
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(visNetwork)
  library(igraph)
  
  all_results <- list()  # Initialize a list to store results
  
  for (i in seq_along(gene_sets)) {
    gene_set <- gene_sets[[i]]
    gene_set <- gene_set[gene_set %in% names(V(whole_network))]
    
    # Set gene_set_name for each gene set
    gene_set_name <- if (!is.null(names(gene_sets))) names(gene_sets)[i] else paste("GeneSet", i, sep = "_")
    
    i= node_positions
    
    vids_df <- which(V(whole_network)$name %in% gene_set)
    
    # Create subgraph for the gene set
    subgraph <- induced_subgraph(whole_network, vids = vids_df)
    
    # Extract node names from the subgraph
    nodes <- V(subgraph)$name
    
    # Prepare data for visNetwork, conditional on centrality_df being provided
    nodes_df <- data.frame(
      id = nodes,
      label = nodes,
      font.size = 32  # Increase font size
    )
    
    if (!is.null(centrality_df)) {
      # Filter centrality_df for nodes in the subgraph
      subgraph_centrality <- centrality_df %>%
        dplyr::filter(feature_id %in% nodes) %>%
        dplyr::arrange(feature_id)
      
      # Merge with nodes_df to retain all nodes, even those without centrality measures
      nodes_df <- merge(nodes_df, subgraph_centrality, by.x = "id", by.y = "feature_id", all.x = TRUE)
      
      # Color nodes based on centrality, using a default color if centrality is NA
      nodes_df$color <- ifelse(is.na(nodes_df[[centrality_measure]]),
                               "#f0e3ab",  # Default color for nodes without centrality
                               scales::col_numeric(palette = c("lightyellow", "red"), domain = NULL)(nodes_df[[centrality_measure]]))
    } else {
      # Default color if no centrality data is provided
      nodes_df$color <- "#f0e3ab"  # Default color for all nodes
    }
    
    # Add positions to nodes_df
    nodes_df <- merge(nodes_df, node_positions, by = "id", all.x = TRUE)
    
    # Annotate gene names using mapIds
    ann_df_gene_name <- AnnotationDbi::mapIds(keys = as.character(nodes_df$id), x = org.Hs.eg.db, keytype = "ENTREZID", column = "SYMBOL") %>%
      tibble::enframe(name = "ENTREZID", value = "SYMBOL") %>%
      data.frame(stringsAsFactors = FALSE)
    
    nodes_df$label <- ann_df_gene_name[match(nodes_df$id, ann_df_gene_name$ENTREZID), 'SYMBOL']
    
    edges_df <- igraph::as_data_frame(subgraph, what = "edges")
    
    vis_net <- visNetwork(nodes_df, edges_df) %>%
      visNodes(font = list(size = 45)) %>%
      visEdges(color = list(color = "black"), ) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visLayout(randomSeed = 123)
    
    clean_name <- gsub("[^A-Za-z0-9_]", "", gene_set_name)
    
    # Save visualization as HTML file
    file_name <- paste0("Gene_Set_", clean_name, "_", centrality_measure, ".html")
    visSave(vis_net, file = file_name)
    cat("Saved visualization to", file_name, "\n")
    
    # Append results to the list
    all_results[[gene_set_name]] <- list(nodes = nodes_df, edges = edges_df)
  }
  
  # Return the list of all signatures with their respective node and edge data frames
  return(all_results)
}


## 1. Load signatures and disease network
## ============================================
signatures = readRDS(get_workflow_outputs("wf-b33649db09"))
  signatures$x2 = signatures$x2[!str_detect(names(signatures$x2),"nc_50|smoothedRandom")]
  signatures = unlist(signatures[-3], recursive = F)
  names(signatures) = str_remove(names(signatures),"x2.|mast.")
  names(signatures) = str_remove(names(signatures),"CytoSig_50Genes__Undivided__|P01__BioNetSmooth_signatures__|P07__CytoSig__|P03__Validations__")

nw_bulk_pruned <- get_workflow_outputs("wf-06da2e64ae")
nw_bulk_pruned <- readRDS(nw_bulk_pruned['igraph_NW_sub_0.6.rds',1])

centrality_df <- readRDS(get_workflow_outputs("wf-d535848bea"))


## 2. Visualization
## ====================
# note - this function exports HTML files so make sure to setwd beforehand

node_positions= get_node_positions(nw_bulk_pruned)

# PageRank centrality
setwd("~/analysis-s05/figures/X2_Signature/Network/PageRank")
res <- visualize_gene_set_subgraphs(nw_bulk_pruned, signatures, centrality_df = centrality_df, centrality_measure = "page_rank")

# Eigen centrality
setwd("~/analysis-s05/figures/X2_Signature/Network/EigenCentrality")
res <- visualize_gene_set_subgraphs(nw_bulk_pruned, signatures, centrality_df = centrality_df, centrality_measure = "eigen_centrality")

