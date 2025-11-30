visualize_mono_gene_set_with_connection <- function(whole_network, 
                                                    gene_sets, 
                                                    centrality_df = NULL, 
                                                    centrality_measure = "betweenness",
                                                    create_visualization = F) {
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  
  all_results <- list()  # Initialize a list to store results
  
  for (i in seq_along(gene_sets)) {
    gene_set <- gene_sets[[i]]
    
    # Filter the gene_set to include only nodes present in the graph
    valid_nodes <- names(V(whole_network))
    gene_set_original <- gene_set[gene_set %in% valid_nodes]
    gene_set <- gene_set[gene_set %in% valid_nodes]
    
    if (length(gene_set) < 2) {
      cat("Skipping gene set", i, "- it contains fewer than 2 valid nodes in the network.\n")
      next
    }
    
    # Set gene_set_name for each gene set
    gene_set_name <- if (!is.null(names(gene_sets))) names(gene_sets)[i] else paste("GeneSet", i, sep = "_")
    
    # Step 1: Extract subgraph for the gene_set
    subgraph <- induced_subgraph(whole_network, vids = which(V(whole_network)$name %in% gene_set))
    
    # Step 2: Identify and reconnect isolated nodes
    isolated_nodes <- V(subgraph)$name[igraph::degree(subgraph) == 0]
    
    if (length(isolated_nodes) > 0) {
      nodes_to_add <- c()
      
      for (node in isolated_nodes) {
        # Find the shortest path from the isolated node to any node in the subgraph
        shortest_paths_result <- shortest_paths(whole_network, from = node, to = V(subgraph)$name, mode = "all")
        
        # Determine the closest node based on the shortest path length
        path_lengths <- sapply(shortest_paths_result$vpath, length) - 1
        min_length <- min(path_lengths[path_lengths > 0])
        
        # Get the closest node(s) based on the shortest path
        closest_node <- V(whole_network)$name[which(path_lengths == min_length)[1]]
        
        # Retrieve all nodes along the shortest path to the closest node
        path_to_closest <- shortest_paths_result$vpath[[which(path_lengths == min_length)[1]]]
        nodes_to_add <- c(nodes_to_add, V(whole_network)$name[path_to_closest])
      }
      
      # Add the identified nodes to the gene set
      added_nodes <- unique(nodes_to_add)
      gene_set <- unique(c(added_nodes, gene_set))
      
      # Recompute the subgraph to include the new nodes
      subgraph <- induced_subgraph(whole_network, vids = which(V(whole_network)$name %in% gene_set))
    }
    
    # Step 3: Prepare node and edge data for visualization
    nodes <- V(subgraph)
    nodes_df <- data.frame(
      id = nodes$name,
      label = nodes$name,
      color = ifelse(nodes$name %in% gene_set & !(nodes$name %in% added_nodes), "#ADD8E6", "#D3D3D3"),  # Original nodes: light blue, added nodes: grey
      font.size = 32
    )
    
    # Step 4: Add gene name annotations using mapIds
    library(org.Hs.eg.db)
    
    ann_df_gene_name <- AnnotationDbi::mapIds(
      keys = as.character(nodes_df$id),
      x = org.Hs.eg.db,
      keytype = "ENTREZID",
      column = "SYMBOL"
    ) %>%
      tibble::enframe(name = "ENTREZID", value = "SYMBOL") %>%
      data.frame(stringsAsFactors = FALSE)
    
    nodes_df$label <- ann_df_gene_name[match(nodes_df$id, ann_df_gene_name$ENTREZID), 'SYMBOL']
    
    # Step 5: Adjust color logic to ensure original nodes stay light blue and added nodes are grey
    original_nodes <- gene_set[gene_set %in% gene_set_original]  # Original nodes from the gene set
    
    nodes_df$color <- ifelse(nodes_df$id %in% original_nodes, "#ADD8E6", "#D3D3D3")
    
    # If centrality data is provided, adjust colors based on centrality
    if (!is.null(centrality_df)) {
      centrality_info <- centrality_df %>%
        dplyr::filter(feature_id %in% nodes$name) %>%
        dplyr::arrange(feature_id)
      
      nodes_df <- merge(nodes_df, centrality_info, by.x = "id", by.y = "feature_id", all.x = TRUE)
      
      nodes_df$color <- ifelse(
        is.na(nodes_df[[centrality_measure]]),
        nodes_df$color,  # Keep default color if no centrality data
        scales::col_numeric(palette = c("lightyellow", "red"), domain = NULL)(nodes_df[[centrality_measure]])
      )
    }
    
    # Extract edges for visNetwork
    edges_df <- igraph::as_data_frame(subgraph, what = "edges")
    
    if(create_visualization) {
      # Step 6: Create visNetwork visualization
      vis_net <- visNetwork(nodes_df, edges_df) %>%
        visNodes(font = list(size = 45)) %>%
        visEdges(color = list(color = "black")) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visLayout(randomSeed = 123)
      
      clean_name <- gsub("[^A-Za-z0-9_]", "", gene_set_name)
      
      # Save visualization as HTML file
      file_name <- paste0("Connected_Gene_Set_", clean_name, "_", centrality_measure, ".html")
      visSave(vis_net, file = file_name)
      cat("Saved connected subgraph visualization to", file_name, "\n")
    }

    
    # Append results to the list
    all_results[[gene_set_name]] <- list(nodes = nodes_df, edges = edges_df)
    
  }
  
  # Return the list of all signatures with their respective node and edge data frames
  return(all_results)
}




geneset_network_measures_random_sample_parallel <- function(seed = 1,NW_WF, genesets_list,n=1000, geneset_name, geneset_network_measures_fun){
  library(igraph)
  library(cytoreason.ccm.pipeline)
  library(cytoreason.cc.client)
  library(cytoreason.io)
  library(dplyr)
  
  # library(analysis.p03.POC1)
  # chicki_Startup()
  
  net_edges <- read_data(NW_WF)
  net_edges <- net_edges %>% dplyr::select(var1,var2,cor)
  
  net_genes <- unique(c(net_edges$var1, net_edges$var2))
  
  #for each geneset create a list of a random genesets of the same size from the network genes
  
  print(geneset_name)
  n_genes_to_sample <- sum(genesets_list[[geneset_name]] %in% net_genes)
  if(n_genes_to_sample > 1){
    
    set.seed(seed)
    random_genesets_list <- lapply(1:n,function(i){
      sample(net_genes, n_genes_to_sample, replace = TRUE)
    })
    
    names(random_genesets_list) <- paste0('rand_iter_',seed,"_",1:n)
    
    #combine geneset with random genesets to one list
    genesets_combined_list <- c(list(gene_list=genesets_list[[geneset_name]]),random_genesets_list)
    
    genesets_combined_list_measures <- geneset_network_measures_fun(net_edges,genesets_combined_list)
    genesets_combined_list_measures$ListName <- geneset_name
    genesets_combined_list_measures$ListType <- genesets_combined_list_measures$geneset
    
    return(genesets_combined_list_measures)
  } else{
    return()
  }
}


geneset_network_measures <- function(net,genesets_list){
  library(dplyr)
  library(igraph)
  library(org.Hs.eg.db)
  round_num <- function(x,n=3){round(x,n)}
  
  net_edges <- net %>%
    #net_edges <- net$modules_obj$edges_membership %>% 
    dplyr::select(var1,var2,cor) %>% 
    mutate(cor = case_when(cor > 0 ~ cor,
                           T ~ 0))
  
  net_genes <- unique(c(net_edges$var1, net_edges$var2))
  # Get subraph measures for each geneset and combine to a dataframe
  
  do.call('rbind', lapply(1:length(genesets_list), function(i){
    genest_to_test <- genesets_list[[i]]
    
    #filter genes based on the genes in the network
    genest_to_test  <-genest_to_test[genest_to_test %in% net_genes]
    if(length(genest_to_test) > 2){
      #edge list to igraph object
      g <- graph_from_edgelist(as.matrix(net_edges[,c("var1","var2")]),directed = F)
      g.weighted <- set.edge.attribute(g, "weight", index=E(g), net_edges$cor)
      
      #create a subgrapgh
      net.sub <- induced_subgraph(g, genest_to_test)
      
      
      # graph density (GD), a score that can also be defined as the local clustering coefficient. formula: 2E/n(n-1)
      DG <- edge_density(net.sub)  
      subnet.diameter <- diameter(net.sub, unconnected=T) #sub-network diameter (D), which uses the maximum length of all shortest paths between any two connected nodes
      
      ##################################
      #get the largest component group
      largest_component <- which(components(net.sub)$csize==max(components(net.sub)$csize))[1]
      #genes of the largest component
      largest_component_genes <- names(components(net.sub)$membership)[components(net.sub)$membership==largest_component]
      genes_entrez <- paste0(largest_component_genes,collapse = ',')

      number_of_connected_components <- length(components(net.sub)$csize)
      largest_component_frac <- round_num(max(components(net.sub)$csize)/sum(components(net.sub)$csize)) #max connected component (size) - fraction
      all_components_frac <- paste0(round_num(sort(components(net.sub)$csize,decreasing = T)/sum(components(net.sub)$csize)),collapse = ';')
      
      #create a sub-sub-graph of the maximum component
      net.sub.max.comp <- induced_subgraph(g, largest_component_genes)
      net.sub.max.comp.weighted <- induced_subgraph(g.weighted, largest_component_genes)
      
      largest_connected_component_density <- edge_density(net.sub.max.comp)
      largest_connected_component_diameter <- diameter(net.sub.max.comp,unconnected = F) 
      largest_connected_component_diameter_weighted <- diameter(net.sub.max.comp.weighted,unconnected = F) 
      
      #create a measure for diameter of connected component relative to its size (INVERSE to deal with zeros, the larger the value the "better")
      largest_connected_component_diameter_frac <- max(components(net.sub)$csize)/largest_connected_component_diameter
      largest_connected_component_diameter_weighted_frac <- max(components(net.sub)$csize)/largest_connected_component_diameter_weighted
      
      return(data.frame(geneset=names(genesets_list)[i], 
                        number_of_genes = length(genesets_list[[i]]),
                        number_of_genes_in_network = length(genest_to_test),
                        density=DG,
                        diameter=subnet.diameter,
                        number_of_connected_components=number_of_connected_components,
                        largest_component_frac = largest_component_frac,
                        all_components_frac = all_components_frac,
                        largest_connected_component_density = largest_connected_component_density,
                        largest_connected_component_diameter = largest_connected_component_diameter,
                        largest_connected_component_diameter_weighted = largest_connected_component_diameter_weighted,
                        largest_connected_component_entrez = genes_entrez,
                        largest_connected_component_diameter_frac =largest_connected_component_diameter_frac,
                        largest_connected_component_diameter_weighted_frac = largest_connected_component_diameter_weighted_frac,
                        stringsAsFactors = F))
    } else {
      print(paste0("geneset ",names(genesets_list)[i]," has ",length(genest_to_test)," genes in the network"))
      return()
    }
  }))
  
}


CollapseToParam <- function(x) {
  library(tidyverse)
  x <- x %>%
    select(ListType, ListName,
           density,
           largest_component_frac,
           largest_connected_component_density,
           largest_connected_component_diameter_frac,
           largest_connected_component_diameter_weighted_frac) %>%
    mutate(across(where(is.numeric),
                  ~ ifelse(is.infinite(.) | is.nan(.), 0, .))) %>%
    pivot_longer(cols = c(density,
                          largest_component_frac,
                          largest_connected_component_density,
                          largest_connected_component_diameter_frac,
                          largest_connected_component_diameter_weighted_frac),
                 names_to = "Criteria",
                 values_to = "value") %>%
    pivot_wider(names_from = ListType, values_from = value) %>%
    mutate(pval = 1 - rowMeans(gene_list > as.matrix(across(starts_with("rand_iter")))))  %>%
    dplyr::select(ListName,Criteria,pval)
  return(x)
}
