devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(dplyr)
library(stringr)

# calculation is done on the disease model as we don't need the signatures
ccm <- as_ccm_fit("wf-08a6a0a503")

## Calculate key pathways PCA
## ==========================================
run_function_dist(FUN = function(ccm_wfid, submodel){
  library(cytoreason.ccm.pipeline)
  library(cytoreason.integration)
  library(dplyr)
  
  ccm = as_ccm_fit(ccm_wfid)
  config <- modelMetadata(ccm)
  keyPathways <- build_service_result_tables(ccm$meta$AD$gx_diff$gx_gsa) # we use DZ vs HC to keep track with cell meta pca
  keyPathways <- keyPathways[[submodel]]$gsa_enrichment %>%
    .[c("h","btm","kegg","reactome")]
  keyPathways = do.call(rbind, lapply(names(keyPathways), function(nm) {
    df <- keyPathways[[nm]]$gsa_enrichment
    df$collection <- nm
    df = df[which(df$term == "DZ_vs_HC"),]
    df
  }))
  keyPathways <- keyPathways %>%
    dplyr::filter(FDR <= 0.05) %>%
    dplyr::mutate(ID = paste0(collection,":",pathway)) %>%
    dplyr::select(ID)
  
  # Get pathway ssgsea for all the samples in all datasets:
  pathways <- lapply(ccm$datasets, function(d){
    analysisResultExpressionSet(d, "gene_set_activity") %>%
      .[keyPathways$ID,]
  })
  
  # Get the data for the meta-PCA train datasets (as defined in the config file):
  pathways.train <- pathways[config$dataset_id[which(config$ccm_meta_pca == 1)]]
  
  # Train the meta PCA:
  metaPCA <- service_meta_pca(data = pathways.train, n_pc = 3)
  
  # Project the samples (all the datasets) on the meta PCA
  metaPCA.projected <- lapply(pathways, function(x) service_meta_pca_projection(x, meta_pca = metaPCA))
  
  return(list(projected = metaPCA.projected, metaPCA = metaPCA, training = pathways.train, keyPathways = keyPathways))
}, 
ccm_wfid = "wf-08a6a0a503", 
submodel = "adjusted__1__1",
image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest", 
tags = list(list(name="analysis",value="adj_pathway_metaPCA")))
# wf-ee9a8f05a7 - bulk
# wf-98fc25d966 - adjusted



## Extraction
## ==========================================
keyPathways_wfid = "wf-ee9a8f05a7" # bulk
keyPathways_wfid = "wf-98fc25d966" # adjusted
metaPCA_pathways = readRDS(get_workflow_outputs(keyPathways_wfid))


## Scores
## ------------
metaPCA.projected = metaPCA_pathways$projected
pathwayPCA <- cytoreason.ccm.pipeline:::statistic_table.ccm_service_meta_pca_projection(metaPCA.projected)
sampleScores <- pathwayPCA$sample_loadings
sampleScores_mat <- reshape(sampleScores[,c("submodel","pc","sample_id","value")], idvar = c("submodel","sample_id"), timevar = "pc", direction = "wide")
colnames(sampleScores_mat) <- gsub("value.","pathway_meta_",colnames(sampleScores_mat))
colnames(sampleScores_mat)[colnames(sampleScores_mat)=="submodel"] <- "dataset"

# Check the direction of the meta-PCs first - I needed to flip the sign of the meta-PC1 for both cells and pathways (because sample loadings of HC>DZ):
sampleScores_mat$pathway_meta_pc1 <- (-1) * sampleScores_mat$pathway_meta_pc1

pushToCC(sampleScores_mat, tagsToPass = list(list(name="analysis",value="adj_pathway_meta_pca")))
# wf-9714e3a025 - bulk
# wf-72065e3e29 - adjusted


## Pathway loadings
## --------------------
metaPCA.res = cytoreason.ccm.pipeline:::statistic_table.ccm_service_meta_pca(list(metaPCA_pathways$metaPCA))
pathLoadings = metaPCA.res$feature_loadings
pathLoadings$pc = str_replace(pathLoadings$pc, "pc", "PC")
pathLoadings = pathLoadings[which(pathLoadings$pc %in% c("PC1","PC2","PC3")),]
pathLoadings$value[which(pathLoadings$pc == "PC1")] = (-1) * pathLoadings$value[which(pathLoadings$pc == "PC1")]

pathLoadings.BQ = pathLoadings[,c("pc","feature_id","value")]
colnames(pathLoadings.BQ) = c("PC","Pathway","Loading")
pathLoadings.BQ$submodel = "adjusted"
# uploadToBQ(pathLoadings.BQ, bqdataset = "s05_atopic_dermatitis", tableName = "pathwayLoadings")
uploadToBQ(pathLoadings.BQ, bqdataset = "s05_atopic_dermatitis", tableName = "pathwayLoadings", disposition = "WRITE_APPEND")


## Visualization
## ==============================
pathLoadings$pathway = sapply(pathLoadings$feature_id, function(x) str_split(x,":", n = 2)[[1]][2])
any(table(pathLoadings$pathway) > 3)
# pathLoadings$pathway[which(pathLoadings$pathway == "Apoptosis")] <- pathLoadings$feature_id[which(pathLoadings$pathway == "Apoptosis")]
pathLoadings = reshape2::dcast(pathLoadings, pathway ~ pc, value.var = "value")

# top 50 loadings
toploadings = apply(pathLoadings[,2:4], 2, function(x){
  x = data.frame(x = abs(x), row.names = pathLoadings$pathway)
  x = top_n(x, n = 50)
  return(rownames(x))
}) %>% data.frame() 
toploadings = pivot_longer(toploadings, names_to = "PC", values_to = "Pathway", cols = 1:3)
toploadings = merge(toploadings, reshape2::melt(pathLoadings, id.vars = 1, variable.name = "PC", value.name = "loading"), all.x = T,
                    by.x = c("Pathway","PC"), by.y = c("pathway","PC"))

# bulk
# --------------------
# cosemetics
toploadings$Pathway[which(str_detect(toploadings$Pathway, "Nucleotide-binding domain, leucine"))] <-  "Nucleotide-binding domain, NLR signaling pathways"
toploadings$Pathway[which(str_detect(toploadings$Pathway, "HDR"))] <-  "HDR through Homologous Recombination or Single Strand Annealing"
toploadings$Pathway[which(str_detect(toploadings$Pathway, "D-loop Structures through Syn"))] <-  "Resolution of D-loop Structures through Synthesis-Dependent Strand Annealing"

toploadings$Pathway = cytoreason.gx::reorder_within(toploadings$Pathway, toploadings$loading, within = toploadings$PC)

ggplot(toploadings, aes(x = loading, y = Pathway)) +
  geom_col() +
  facet_wrap(~PC, scales = "free") +
  cytoreason.gx::scale_y_reordered()+
  theme_minimal() +
  ggpubr::border() +
  labs(x = "Pathway Loading on Meta-PCA", y = "Pathway") +
  theme(strip.text = element_text(size = 16), axis.text = element_text(size = 14), axis.title = element_text(size = 14))
ggsave("~/analysis-s05/figures/AD Model/bulkPathwayLoadings.top50.png", scale = 2, width = 15, height = 8, bg = "white")


# adjusted
# --------------------
# cosemetics
toploadings$Pathway[which(str_detect(toploadings$Pathway, "APC:"))] <-  "APC:Cdc20 mediated degradation of cell cycle proteins"
toploadings$Pathway[which(str_detect(toploadings$Pathway, "HDR"))] <-  "HDR through Homologous Recombination or Single Strand Annealing"
toploadings$Pathway[which(str_detect(toploadings$Pathway, "D-loop Structures through Syn"))] <-  "Resolution of D-loop Structures through Synthesis-Dependent Strand Annealing"
toploadings$Pathway[which(str_detect(toploadings$Pathway, "Activation of the mRNA"))] <-  "mRNA activation via cap-binding complex and eIFs, followed by 43S binding"
toploadings$Pathway <- str_replace_all(toploadings$Pathway, " (NMD)","")
toploadings$Pathway <- str_replace_all(toploadings$Pathway, " (EJC)","")
toploadings$Pathway <- str_replace_all(toploadings$Pathway, " (HRR)","")
toploadings$Pathway <- str_replace_all(toploadings$Pathway, " (SSA)","")

toploadings$Pathway = cytoreason.gx::reorder_within(toploadings$Pathway, toploadings$loading, within = toploadings$PC)

ggplot(toploadings, aes(x = loading, y = Pathway)) +
  geom_col() +
  facet_wrap(~PC, scales = "free") +
  cytoreason.gx::scale_y_reordered()+
  theme_minimal() +
  ggpubr::border() +
  labs(x = "Pathway Loading on Meta-PCA", y = "Pathway") +
  theme(strip.text = element_text(size = 16), axis.text = element_text(size = 14), axis.title = element_text(size = 14))
ggsave("~/analysis-s05/figures/AD Model/adjPathwayLoadings.top50.png", scale = 2, width = 18, height = 8, bg = "white")
