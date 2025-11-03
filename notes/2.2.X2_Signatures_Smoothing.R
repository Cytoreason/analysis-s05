library(cytoreason.cc.client)
library(cytoreason.patern)


# Vector prep
# ----------------------------------
gxdiff = readRDS(get_workflow_outputs("wf-409f6514bd"))

# 1. X2 activations
activation = gxdiff[which(gxdiff$term == "activated_vs_unactivated"),]
activation = reshape2::dcast(activation, feature_id ~ comparison, value.var = "estimate")


alphas <- c(0.2, 0.4, 0.6, 0.8, 0.99)
tasks_smoothing <- sapply_dist(X = alphas, function(alpha,
                                                    cor_mat,
                                                    fc_mat,
                                                    edge_cutoff,
                                                    edge_cutoff_abs,
                                                    stop_val,
                                                    direction_split){
  library(cytoreason.patern)
  service_patern(cor_mat = cor_mat,
                 fc_mat = fc_mat,
                 alpha = alpha,
                 edge_cutoff = edge_cutoff,
                 edge_cutoff_abs = edge_cutoff_abs,
                 stop_val = stop_val,
                 direction_split = direction_split)
},
cor_mat = input_local_file(filename = adj_file, load = TRUE),
fc_mat = prioritization_df,
edge_cutoff = NA,
edge_cutoff_abs = FALSE,
stop_val = 1e-03,
direction_split = T,
image = patern_image,
cpu_request = '2000m',
memory_request = '25Gi')


# Propagation in ARCHS
# --------------------------------
adj_matrix <- readRDS(get_task_outputs('wf-3917825f8b','0'))
tags <- list(list(name="project", value="target-classifier"),
             list(name="service", value="patern"),
             list(name="network", value="ARCHS"),
             list(name="network_wfid", value="wf-3917825f8b"),
             list(name="alphas", value=paste0(alphas, collapse = '_')),
             list(name="fc_values", value='binary'))




# Propagation in STRINGdb
# -------------------------------
net_file <- get_task_inputs(res = 'wf-690e1e4e36', task_id = 0, files_names_grepl_pattern = '^g$')
adj_matrix <- readRDS(net_file)
tags <- list(list(name="project", value="target-classifier"),
             list(name="service", value="patern"),
             list(name="network", value="STRINGdb"),
             list(name="network_wfid", value="wf-690e1e4e36"),
             list(name="alphas", value=paste0(alphas, collapse = '_')),
             list(name="fc_values", value='binary'))