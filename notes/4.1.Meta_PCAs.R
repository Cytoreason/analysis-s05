devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)

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
  allKeyPathways = lapply(unique(names(ccm$meta)[1:4]), function(effect){
    allSubmodels = lapply(c("bulk","adjusted__1__1","adjusted__1__CRCL_0000348"), function(submodel){
      cat("\nComputing:",effect,submodel,"\n")
      keyPathways <- build_service_result_tables(ccm$meta[[effect]]$gx_diff$gx_gsa) # we use DZ vs HC to keep track with cell meta pca
      keyPathways <- keyPathways$gsa_enrichment %>%
        dplyr::filter(collection %in% c("h","btm","kegg","reactome")) %>%
        dplyr::filter(fdr <= 0.05) %>%
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
      
      # now per collection
      cat("\nPer Collection\n")
      perCollection = lapply(c("h","btm","kegg","reactome"), function(collection){
        pathw = lapply(pathways, function(x) x[stringr::str_detect(rownames(x),paste0(collection,":"))])
        pathways.train <- pathw[config$dataset_id[which(config$ccm_meta_pca == 1)]]
        metaPCA <- service_meta_pca(data = pathways.train, n_pc = 3)
        metaPCA.projected <- lapply(pathw, function(x) service_meta_pca_projection(x, meta_pca = metaPCA))
        return(list(projected = metaPCA.projected, metaPCA = metaPCA, training = pathways.train, keyPathways = keyPathways))
      })
      names(perCollection) = c("h","btm","kegg","reactome")
      
      return(list(projected = metaPCA.projected, metaPCA = metaPCA, training = pathways.train, keyPathways = keyPathways, perCollection = perCollection))
      
    })
    names(allSubmodels) = c("bulk","adjusted__1__1","adjusted__1__CRCL_0000348")
    return(allSubmodels)
  })
  names(allKeyPathways) = unique(names(ccm$meta)[1:4])
  return(allKeyPathways)
}, 
ccm_wfid = "wf-832ab799be", # new disease model, including adjustment to keratinocytes
image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest", 
tags = list(list(name="analysis",value="pathway_metaPCA_allTerms")))
# wf-ee9a8f05a7 - bulk
# wf-98fc25d966 - adjusted
# wf-0a6f338228 - all terms
# wf-22b001bd47 - all terms, including keratinocytes adjustment


## Another option: remove redundant pathways
## ==================================================
library(cytoreason.omics)
allPathways = analysisResultExpressionSet(ccm$datasets$GSE107361__GPL570,"gene_set_activity")
  allPathways = featureData(allPathways)
  allPathways = allPathways@data

names(Reactome) = paste0("reactome:",names(Reactome))
kegg = KEGG_pathways()
names(kegg) = paste0("kegg:",names(kegg))
btm = BTM
names(btm) = paste0("btm:",names(btm))
h = cytoreason.omics::gs_msigdb$ENTREZID$h
names(h) = paste0("h:",names(h))

pathwayGeneset = c(Reactome, kegg, h, btm)

overlap = sapply(pathwayGeneset, function(x) sapply(pathwayGeneset, function(y) length(intersect(y,x))/min(length(x),length(y))))

## Extraction
## ==========================================
keyPathways_wfid = "wf-ee9a8f05a7" # bulk
keyPathways_wfid = "wf-98fc25d966" # adjusted
keyPathways_wfid = "wf-0a6f338228" # all
metaPCA_pathways = readRDS(get_workflow_outputs(keyPathways_wfid))


## Scores
## ------------
extract_sampleScores = function(metaPCA.projected, flipPC1, flipPC2) {
  pathwayPCA <- cytoreason.ccm.pipeline:::statistic_table.ccm_service_meta_pca_projection(metaPCA.projected)
  sampleScores <- pathwayPCA$sample_loadings
  sampleScores_mat <- reshape(sampleScores[,c("submodel","pc","sample_id","value")], idvar = c("submodel","sample_id"), timevar = "pc", direction = "wide")
  colnames(sampleScores_mat) <- gsub("value.","pathway_meta_",colnames(sampleScores_mat))
  colnames(sampleScores_mat)[colnames(sampleScores_mat)=="submodel"] <- "dataset"
  
  # Check the direction of the meta-PCs first - I needed to flip the sign of the meta-PC1 for both cells and pathways (because sample loadings of HC>DZ):
  if(flipPC1) { sampleScores_mat$pathway_meta_pc1 <- (-1) * sampleScores_mat$pathway_meta_pc1 }
  if(flipPC2) { sampleScores_mat$pathway_meta_pc2 <- (-1) * sampleScores_mat$pathway_meta_pc2 }
  return(sampleScores_mat)
}

process = function(res, term, submodel, collection) {
  res$term = term
  res$submodel = submodel
  if(!is.null(collection)){
    res$collection = collection
  } else {
    res$collection = "all"
  }
  return(res)
}


# only bulk
sampleScores = extract_sampleScores(metaPCA_pathways$projected, T, T)
# wf-31ea889fa9

# only adjusted (meta PC1)
sampleScores = extract_sampleScores(metaPCA_pathways$projected, T, T)
# wf-72065e3e29 - adjusted

# all terms and submodel
sampleScores_all = lapply(names(metaPCA_pathways), function(term){
  lapply(names(metaPCA_pathways[[term]]), function(submodel){
    cat("\n", term, submodel,"\n")
    res = extract_sampleScores(metaPCA_pathways[[term]][[submodel]]$projected, T, T)
    res = process(res, term, submodel, collection = NULL)
    res_perCollection = lapply(names(metaPCA_pathways[[term]][[submodel]][["perCollection"]]), function(collection){
      cat("\r", term, submodel,"...",collection)
      res_col = extract_sampleScores(metaPCA_pathways[[term]][[submodel]][["perCollection"]][[collection]]$projected, T, T)
      res_col = process(res_col, term, submodel, collection)
      return(res_col)
    }) %>% bind_rows()
    cat("\r", term, submodel,"........done")
    return(bind_rows(res, res_perCollection))
  }) %>% bind_rows()
}) %>% bind_rows()
# wf-e8f514f3f1
# wf-5b90576980 - with keratinocytes
pushToCC(sampleScores_all, tagsToPass = list(list(name="analysis",value="pathway_meta_pca")))


## Pathway loadings
## --------------------
extract_pathwayLoadings = function(metaPCA, flipPC1, flipPC2){
  metaPCA.res = cytoreason.ccm.pipeline:::statistic_table.ccm_service_meta_pca(list(metaPCA))
  pathLoadings = metaPCA.res$feature_loadings
  pathLoadings$pc = str_replace(pathLoadings$pc, "pc", "PC")
  pathLoadings = pathLoadings[which(pathLoadings$pc %in% c("PC1","PC2","PC3")),]
  if(flipPC1) { pathLoadings$value[which(pathLoadings$pc == "PC1")] = (-1) * pathLoadings$value[which(pathLoadings$pc == "PC1")] }
  if(flipPC2) { pathLoadings$value[which(pathLoadings$pc == "PC2")] = (-1) * pathLoadings$value[which(pathLoadings$pc == "PC2")] } # only in bulk
  return(pathLoadings)
}


# only bulk
pathwayLoadings = extract_pathwayLoadings(metaPCA_pathways$metaPCA, T, T)

# only adjusted (meta PC1)
pathwayLoadings = extract_pathwayLoadings(metaPCA_pathways$metaPCA, T, F)

# all terms and submodel
pathwayLoadings_all = lapply(names(metaPCA_pathways), function(term){
  lapply(names(metaPCA_pathways[[term]]), function(submodel){
    cat("\n", term, submodel,"\n")
    flipPC2 = ifelse(submodel == "bulk", T, F)
    res = extract_pathwayLoadings(metaPCA_pathways[[term]][[submodel]]$metaPCA, T, flipPC2)
    res = process(res, term, submodel, collection = NULL)
    res_perCollection = lapply(names(metaPCA_pathways[[term]][[submodel]][["perCollection"]]), function(collection){
      cat("\r", term, submodel,"...",collection)
      res_col = extract_pathwayLoadings(metaPCA_pathways[[term]][[submodel]][["perCollection"]][[collection]]$metaPCA, T, flipPC2)
      res_col = process(res_col, term, submodel, collection)
      return(res_col)
    }) %>% bind_rows()
    cat("\r", term, submodel,"........done")
    return(bind_rows(res, res_perCollection))
  }) %>% bind_rows()
}) %>% bind_rows()

pathwayLoadings_all = unique(pathwayLoadings_all)
pathLoadings.BQ = pathwayLoadings_all %>%
  dplyr::select(-c(1:3))
colnames(pathLoadings.BQ) = c("Pathway","PC","Loading","Term","Submodel","Collection")
uploadToBQ(pathLoadings.BQ, bqdataset = "s05_atopic_dermatitis", tableName = "pathwayLoadings")


## Per dataset pathway loadings
## ------------------------------





## Visualization - pathway loadings
## ======================================
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


# key pathways in bulk vs adjusted
# ------------------------------------
keyPathways_bulk = readRDS(get_workflow_outputs("wf-ee9a8f05a7")) %>% .[["keyPathways"]]
keyPathways_adj = readRDS(get_workflow_outputs("wf-98fc25d966")) %>% .[["keyPathways"]]
keyPathways = rbind(cbind(submodel = "bulk", keyPathways_bulk),
                    cbind(submodel = "adjusted", keyPathways_adj))

keyPathways.wide = reshape2::dcast(keyPathways, ID ~ submodel)
keyPathways.wide = keyPathways.wide %>%
  mutate(importance = case_when(!is.na(bulk) & !is.na(adjusted) ~ "stayed important",
                          !is.na(bulk) & is.na(adjusted) ~ "lost importance",
                          is.na(bulk) & !is.na(adjusted) ~ "gained importance"))
table(keyPathways.wide$importance)
# gained importance   lost importance  stayed important 
#                28               446               252 

ontology = openxlsx::read.xlsx("~/data/Pathways_heirarchy_all collections_final.xlsx")
ontology$Full_pathway_name = str_replace(ontology$Full_pathway_name, "__", ":")

# some changes to make them match
keyPathways.wide$ID[str_detect(keyPathways.wide$ID,"h:")] <- str_to_lower(keyPathways.wide$ID[str_detect(keyPathways.wide$ID,"h:")])
keyPathways.wide$ID[str_detect(keyPathways.wide$ID,"reactome:")] <- str_to_lower(keyPathways.wide$ID[str_detect(keyPathways.wide$ID,"reactome:")])
keyPathways.wide$ID <- str_replace(keyPathways.wide$ID,"kegg","KEGG")
keyPathways.wide = read.csv("~/data/keyPathways_wide.csv")

keyPathways.wide$level0 = ontology$level0[match(keyPathways.wide$ID, ontology$Full_pathway_name)]
keyPathways.wide$level1 = ontology$level_1[match(keyPathways.wide$ID, ontology$Full_pathway_name)]
keyPathways_sum = keyPathways.wide %>%
  group_by(importance, level0) %>%
  summarise(nPathways = n()) %>%
  ungroup() %>%
  dplyr::filter(!is.na(level0))

keyPathways_sum$PathwayImportance = cytoreason.gx::reorder_within(keyPathways_sum$level0, keyPathways_sum$nPathways, within = keyPathways_sum$importance)

ggplot(keyPathways_sum, aes(x = importance, y = nPathways, fill = level0, group = PathwayImportance)) +
  geom_col(position = position_dodge2(0.75)) +
  cytoreason.gx::scale_x_reordered()+
  coord_flip()+
  scale_color_brewer(palette = "Set3") +
  theme_minimal() +
  labs(x = NULL, y = "Number of pathways", fill = "Pathway\nCategory")+
  ggtitle("Change in the key disease pathways between bulk and adjusted")
ggsave("~/analysis-s05/figures/AD Model/keyPathways_betweenBulkAndAdjusted.png", scale = 1, width = 8, height = 5, bg = "white")


## Visualizations - pathways - sample level
## ===============================================
unified_metadata = readRDS(get_workflow_outputs("wf-e82a5ab3b6"))
unified_metadata$sample_classification[which(unified_metadata$sample_classification == "Normal")] <- "HC"
pathways_bulk = readRDS(get_workflow_outputs("wf-31ea889fa9"))
pathways_adj = readRDS(get_workflow_outputs("wf-72065e3e29"))
  colnames(pathways_adj)[3:5] = paste0("adj_",colnames(pathways_adj)[3:5])
# cells_bulk = lapply(cellPCA$`1`$x, function(x) x[["coord"]][,1:3]) %>% do.call(rbind,.)

unified_metadata = merge(unified_metadata, pathways_bulk, by = "sample_id")
unified_metadata = merge(unified_metadata, pathways_adj, by = "sample_id")

uploadToBQ(unified_metadata, bqdataset = "s05_atopic_dermatitis", tableName = "AD_sample_metadata")

# bulk
pathways = unified_metadata[,c("sample_id","sample_classification","condition","pathway_meta_pc1", "pathway_meta_pc2")]
colnames(pathways)[4:5] = paste0("Pathway Meta PC",1:2)
pathways$`Pathway Meta PC2` = (-1) * pathways$`Pathway Meta PC2`

pathways = reshape2::melt(pathways, id.vars = 1:3, variable.name = "PC", value.name = "score")
pathways$sample_classification = factor(pathways$sample_classification, ordered = T, levels = c("Lesion","Non Lesion","HC"))

ggplot(pathways[which(pathways$condition %in% c("AD","HC")),], 
       aes(y = score, x = sample_classification, fill = sample_classification)) +
  geom_boxplot(outliers = T) +
  ggpubr::stat_compare_means(method = "wilcox", comparisons = list(c("Lesion","Non Lesion"),c("Non Lesion","HC"),c("Lesion","HC")), paired = F) +
  scale_fill_manual(values = c("#473472","#53629E","#87BAC3"))+
  theme_minimal()+
  facet_wrap(~PC, nrow = 1) +
  theme(legend.position = "none", panel.spacing.x = unit(2, "cm"))+
  labs(x = NULL, title = "Pathway meta PCA scores", y = "Sample Score on the PC")
ggsave("~/analysis-s05/figures/AD Model/pathwayPCs_bulk_sampleClassification.png", bg = "white", scale=1, width = 8, height = 5, units = "in")

ggplot(pathways, aes(y = score, x = condition, fill = condition)) +
  geom_violin(draw_quantiles = T) +
  ggpubr::stat_compare_means(method = "wilcox", comparisons = list(c("AD","contact dermatitis"),
                                                                   c("contact dermatitis","PSO"),
                                                                   c("PSO","HC"),
                                                                   c("AD","PSO"),
                                                                   c("contact dermatitis","HC"),
                                                                   c("AD","HC")), paired = F) +
  scale_fill_manual(values = c("#473472","#53629E","#87BAC3","#D6F4ED"))+
  facet_wrap(~PC, nrow = 1) +
  theme_minimal()+
  labs(x = NULL, title = "Pathway meta PCA scores", y = "Sample Score on the PC") +
  theme(legend.position = "none")
ggsave("~/analysis-s05/figures/AD Model/pathwayPCs_bulk_condition.png", bg = "white", scale = 1.5)

# adjusted
pathways = unified_metadata[,c("sample_id","sample_classification","condition","adj_pathway_meta_pc1", "adj_pathway_meta_pc2", "adj_pathway_meta_pc3")]
colnames(pathways)[4:6] = paste0("Adjusted Pathway Meta PC",1:3)
pathways = reshape2::melt(pathways, id.vars = 1:3, variable.name = "PC", value.name = "score")
pathways$sample_classification = factor(pathways$sample_classification, ordered = T, levels = c("Lesion","Non Lesion","HC"))

ggplot(pathways[which(pathways$condition %in% c("AD","HC")),], aes(y = score, x = sample_classification, fill = sample_classification)) +
  geom_boxplot(outliers = T) +
  ggpubr::stat_compare_means(method = "wilcox", comparisons = list(c("Lesion","Non Lesion"),c("Non Lesion","HC"),c("Lesion","HC")), paired = F) +
  scale_fill_manual(values = c("#473472","#53629E","#87BAC3"))+
  theme_minimal()+
  facet_wrap(~PC, nrow = 1) +
  theme(legend.position = "none")+
  labs(x = NULL, title = "Pathway meta PCA scores", y = "Sample Score on the PC")
ggsave("~/analysis-s05/figures/AD Model/pathwayPCs_adj_sampleClassification.png", bg = "white", scale=1.5)

ggplot(pathways, aes(y = score, x = condition, fill = condition)) +
  geom_violin(draw_quantiles = T) +
  ggpubr::stat_compare_means(method = "wilcox", comparisons = list(c("AD","contact dermatitis"),
                                                                   c("contact dermatitis","PSO"),
                                                                   c("PSO","HC"),
                                                                   c("AD","PSO"),
                                                                   c("contact dermatitis","HC"),
                                                                   c("AD","HC")), paired = F) +
  scale_fill_manual(values = c("#473472","#53629E","#87BAC3","#D6F4ED"))+
  facet_wrap(~PC, nrow = 1) +
  theme_minimal()+
  labs(x = NULL, title = "Pathway meta PCA scores", y = "Sample Score on the PC") +
  theme(legend.position = "none")
ggsave("~/analysis-s05/figures/AD Model/pathwayPCs_adj_condition.png", bg = "white", scale = 1.5)


# by age
# ----------------
withkids = c("GSE107361__GPL570")

dataset = unified_metadata %>%
  dplyr::filter(dataset_id %in% withkids) %>%
  dplyr::filter(!is.na(age_group)) %>%
  dplyr::filter(sample_classification == "Lesion") %>%
  dplyr::filter(time %in% c("D0","Pre","D1")) %>%
  mutate(age_group = case_when(age <= 5 ~ "Under 5",
                               age > 5 & age <= 10 ~ "5-10",
                               age > 10 & age <= 16 ~ "10-16",
                               age > 16 & age <= 18 ~ "16-18", .default = "Adult")) %>%
  mutate(age_group = factor(age_group, ordered = T, levels = c("Under 5","5-10","10-16","16-18","Adult")))

ggplot(dataset, aes(x = age, y = pathway_meta_pc1)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor() +
  facet_grid(~ dataset_id_short, scales = "free_x")

ggplot(dataset, aes(x = age_group, y = pathway_meta_pc1)) +
  geom_boxplot() +
  geom_jitter(position = position_dodge2(0.1))+
  # scale_color_brewer(palette = "PRGn")+
  ggpubr::stat_compare_means(comparisons = list(c("Under 5","Adult")), method = "wilcox") +
  facet_grid(~ dataset_id_short, scales = "free_x")


# ethnicity
# -------------
unique(unified_metadata$ethnicity)
dataset2 = unified_metadata %>%
  dplyr::filter(!is.na(ethnicity) & ethnicity != "NA") %>%
  dplyr::filter(sample_classification == "Lesion") %>%
  dplyr::filter(time %in% c("D0","Pre","D1")) %>%
  mutate(ethnicity = factor(ethnicity))

ggplot(dataset2, aes(x = ethnicity, y = pathway_meta_pc1)) +
  geom_boxplot() +
  geom_jitter(position = position_dodge2(0.1))+
  facet_wrap(~ dataset_id_short, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



## Cell Loadings
## =========================
ccm = as_ccm_fit("wf-08a6a0a503")
cellPCA = ccm$multi$meta_pca
cellLoadings = cellPCA$`1`$v
cellLoadings[,2] = (-1) * cellLoadings[,2]
cellLoadings = reshape2::melt(cellLoadings, variable.name = "PC", value.name = "loading")
colnames(cellLoadings) = c("Cell","PC","Loading")
cellLoadings$PC = paste0("PC",cellLoadings$PC)

uploadToBQ(cellLoadings, bqdataset = "s05_atopic_dermatitis", tableName = "cellLoadings") # only in bulk


## Visualizations - cells - sample level
## ===============================================
cells = lapply(ccm$datasets, function(d){
  p = pData(assayDataExpression(d))
  return(p[,c("sample_id","meta1_pc1","meta1_pc2","meta1_pc3")])
})
cells = bind_rows(cells)
colnames(cells)[-1] = paste0("Cell Meta PC",1:3)
cells$`Cell Meta PC1` = (-1) * cells$`Cell Meta PC1`
cells$`Cell Meta PC2` = (-1) * cells$`Cell Meta PC2`

unified_metadata = readRDS(get_workflow_outputs("wf-e82a5ab3b6"))
unified_metadata$sample_classification[which(unified_metadata$sample_classification == "Normal")] <- "HC"
unified_metadata = merge(unified_metadata, cells, by = "sample_id")

cells = unified_metadata[,c("sample_id","sample_classification","condition","Cell Meta PC1", "Cell Meta PC2")]
cells = reshape2::melt(cells, id.vars = 1:3, variable.name = "PC", value.name = "score")
cells$sample_classification = factor(cells$sample_classification, ordered = T, levels = c("Lesion","Non Lesion","HC"))

ggplot(cells[which(cells$condition %in% c("AD","HC")),], 
       aes(y = score, x = sample_classification, fill = sample_classification)) +
  geom_boxplot(outliers = T) +
  ggpubr::stat_compare_means(method = "wilcox", comparisons = list(c("Lesion","Non Lesion"),c("Non Lesion","HC"),c("Lesion","HC")), paired = F) +
  scale_fill_manual(values = c("#473472","#53629E","#87BAC3"))+
  theme_minimal()+
  facet_wrap(~PC, nrow = 1) +
  theme(legend.position = "none", panel.spacing.x = unit(2, "cm"))+
  labs(x = NULL, title = "Cell meta PCA scores", y = "Sample Score on the PC")
ggsave("~/analysis-s05/figures/AD Model/cellPCs_sampleClassification.png", bg = "white", scale=1, width = 8, height = 5, units = "in")

cells_wide = reshape2::dcast(cells, sample_id + sample_classification + condition ~ PC)
ggplot(cells_wide, aes(x = `Cell Meta PC1`, y = `Cell Meta PC2`)) +
  geom_point(aes(color = sample_classification)) +
  theme_minimal()+
  labs(color=NULL) +
  theme(legend.position = "bottom",legend.direction = "horizontal")
