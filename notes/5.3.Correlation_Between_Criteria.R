devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(ComplexHeatmap)
library(tidyverse)

Results = readRDS(get_workflow_outputs("wf-64470d2a55"))
Results$Target_Cell_PCA$PC = sapply(Results$Target_Cell_PCA$Criteria.Identifier, function(x) str_split(x,"pc")[[1]][2])
Results$Target_Pathway_PCA$PC = sapply(Results$Target_Pathway_PCA$Criteria.Identifier, function(x) str_split(x,"pc")[[1]][2])

Results_binded = Results[setdiff(names(Results), c("Target_Cell", "Target_Pathway", "Target_Gene"))] %>%
  bind_rows() %>%
  mutate(DataType = ifelse(DataType == "Target_CS", paste0("Target_CS_",Criteria.Identifier), DataType)) %>%
  mutate(DataType = ifelse(DataType == "Target_Cell_PCA", paste0("Target_Cell_PC",PC), DataType)) %>%
  mutate(DataType = ifelse(DataType == "Target_Pathway_PCA", paste0("Target_Pathway_PC",PC), DataType)) %>%
  mutate(DataType = str_replace_all(DataType, " ", "_"))

# remove bulk & adj mix
Results_binded = Results_binded[-which(str_detect(Results_binded$DataType,"Target_Pathway_PC") & 
                                        Results_binded$Type == "bulk" & 
                                        str_detect(Results_binded$Criteria.Identifier,"_adj")),]
Results_binded = Results_binded[-which(str_detect(Results_binded$DataType,"Target_Pathway_PC") & 
                                         Results_binded$Type == "adjusted__1__1" & 
                                         str_detect(Results_binded$Criteria.Identifier,"_bulk")),]
Results_binded$DataType[which(Results_binded$Type != "bulk")] <- paste0("Adj_", Results_binded$DataType[which(Results_binded$Type != "bulk")])

Results_binded = Results_binded %>%
  dplyr::filter(!str_detect(DataType, "PC3")) %>%
  mutate(DataType = str_remove(DataType,"_bulk|_adj")) %>%
  dplyr::filter(Target.Collection == "epidermis") # all pathways already exist

Results.corr = reshape2::dcast(Results_binded, Target.Identifier ~ DataType, value.var = "metricValue") %>%
    column_to_rownames(var = "Target.Identifier") %>%
    as.matrix

corr_corr = cor(Results.corr, method = "pearson")

# presented version used wf-a53211cea7 results
png("~/analysis-s05/figures/Results/criteria_correlation_pearson.png", bg = "white", res = "120", width = 1400, height = 1200)
Heatmap(corr_corr, cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(round(corr_corr[i, j],1), x, y, gp = gpar(fontsize = 10))},
  name = "Pearson\nCorrelation", row_labels = str_replace_all(rownames(corr_corr),"_"," "), 
  column_labels = str_replace_all(rownames(corr_corr),"_"," "))
dev.off()

corr_corr = cor(Results.corr, method = "spearman")

png("~/analysis-s05/figures/Results/criteria_correlation_spearman.png", bg = "white", res = "120", width = 1400, height = 1200)
Heatmap(corr_corr, cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(round(corr_corr[i, j],1), x, y, gp = gpar(fontsize = 10))},
  name = "Spearman\nCorrelation", row_labels = str_replace_all(rownames(corr_corr),"_"," "), 
  column_labels = str_replace_all(rownames(corr_corr),"_"," "))
dev.off()


## Correlation between molecular score and EASI scores
## =========================================================
ccm_wfid = "wf-8e948630d7"
ccm= as_ccm_fit(ccm_wfid)

pathways = c("time","sample_classification","condition","experiment_id","EASI","MolecularScore")
enrichment =  lapply(ccm$datasets[which(names(ccm$datasets) %in% c("GSE130588__GPL570","GSE59294__GPL570"))], function(d){
    pData(assayDataExpression(d))[,pathways]
}) %>% do.call(rbind,.)

enrichment_filtered= enrichment %>%
  dplyr::filter(time %in% c("W0","Pre") & sample_classification == "Lesion" & condition == "atopic dermatitis")

ggplot(enrichment_filtered, aes(x = EASI, y = MolecularScore)) +
  geom_point(aes(color = experiment_id), size =3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  ggpubr::stat_cor() +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(title = "Correlaton between EASI scores and Molecular Scores")
ggsave("~/analysis-s05/figures/AD Model/corr_EASI_MS_lesional.png", bg = "white")


enrichment_filtered= enrichment %>%
  dplyr::filter(time %in% c("W0","Pre") & condition == "atopic dermatitis")

ggplot(enrichment_filtered, aes(x = EASI, y = MolecularScore)) +
  geom_point(aes(color = experiment_id), size =3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  ggpubr::stat_cor() +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(title = "Correlaton between EASI scores and Molecular Scores")
ggsave("~/analysis-s05/figures/AD Model/corr_EASI_MS_all.png", bg = "white")

ggplot(enrichment_filtered, aes(x = EASI, y = MolecularScore)) +
  geom_point(aes(color = sample_classification), size =3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  ggpubr::stat_cor() +
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "Correlaton between EASI scores and Molecular Scores")
ggsave("~/analysis-s05/figures/AD Model/corr_EASI_MS_all.png", bg = "white")



## Correlation between X2 signature and EASI scores
## =========================================================
## Because we used only lesional samples from the L vs NL comparison, we need to extract the exact samples
## used in the analysis in order to replicate the results
library(WGCNA)
ccm_wfid = "wf-8e948630d7"
ccm = as_ccm_fit(ccm_wfid)
con = c(GSE130588 = designSampleGroupContrasts("L_vs_NL__GSE130588__GPL570__Lesion_vs_non_lesion", 
                                 data = ccm$datasets$GSE130588__GPL570),
        GSE59294 = designSampleGroupContrasts("L_vs_NL__GSE59294__GPL570__Lesion_vs_non_lesion", 
                                               data = ccm$datasets$GSE59294__GPL570)) %>%
  bind_rows()
con = con[which(con$paired),]

signatureMapping = readRDS(get_workflow_outputs("wf-fcb1ef2bbd"))

pathways = c("time","sample_classification","condition","experiment_id","EASI")
enrichment =  lapply(ccm$datasets[which(names(ccm$datasets) %in% c("GSE130588__GPL570","GSE59294__GPL570"))], function(d){
  a = pData(assayDataExpression(d))[,pathways] %>%
    rownames_to_column(., var = "sample_id")
  b = statistic_table(d$gene_set_activity, subset = ~(submodel == "bulk" & collection == "X2"))
  b = reshape2::dcast(b, sample_id ~ gene_set, value.var = "value")
  return(merge(a,b, by="sample_id"))
}) %>% do.call(rbind,.)

enrichment_baseline = enrichment %>%
  dplyr::filter(sample_classification == "Lesion") %>%
  dplyr::filter(sample_id %in% con$sample_id)

x2 = enrichment_baseline[,c(1:6, which(colnames(enrichment_baseline) %in% c("IL4","x2_general_inhibition_early_50", "x2_general_inhibition_early_50_archs_refined","SP_general_inhibition_late_50_archs_refined")))]
colnames(x2)[7:10] = signatureMapping$New_identifier[match(colnames(x2)[7:10], signatureMapping$signature)]

# Compute bicor + FDR per signature
bicor_tabl = lapply(colnames(x2)[7:10], function(sig){
  bp = WGCNA::bicorAndPvalue(x2$EASI, x2[,sig], use = "pairwise.complete.obs")
  data.frame(
    signature = sig,
    r = as.numeric(bp$bicor),
    p = as.numeric(bp$p)
  )
}) %>% bind_rows() %>%
    mutate(FDR = p.adjust(p, method = "fdr"),
           p_fmt = dplyr::case_when(
             p <= 0.001        ~ sprintf("%.1e", p),  # scientific for ≤ 0.01
             TRUE              ~ sprintf("%.3f", p)   # fixed decimal otherwise
           ),
           label = sprintf("bicor = %.2f, p=%s", r, p_fmt)
)

x2 = reshape2::melt(x2, id.vars = 1:6, variable.name = "signature", value.name = "sample_enrichment")
x2 = x2 %>% left_join(., bicor_tabl, by = "signature")

ggplot(x2, aes(x = EASI, y = sample_enrichment)) +
  geom_point(aes(color = experiment_id)) +
  geom_smooth(method="lm", se = F, color = "black") +
  geom_text(x = Inf, y = Inf, aes(label = label), check_overlap = T, show.legend = F, hjust = 1.1, vjust = 2, size = 3)+
  scale_color_brewer(palette = "Set2")+
  facet_wrap(~signature, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  ggpubr::border()+
  labs(x = "EASI score", y = "Signature enrichment", color = NULL)


#####################
# New correlations
#####################
signatureMapping = readRDS(get_workflow_outputs("wf-aa75ed069b"))
keep_signatures = openxlsx::read.xlsx("~/analysis-s05/data/Final signatures to be ranked or viewed in the dashboards.xlsx")
Results = readRDS(get_workflow_outputs("wf-64470d2a55"))
Results = lapply(Results, function(x){
  x$Target.Collection = sapply(x$Target.ID, function(x) {
    str_split(x, ":")[[1]][1]
  })
  return(x)
})
Results$Target_Cell_PCA$PC = sapply(Results$Target_Cell_PCA$Criteria.Identifier, function(x) str_split(x,"pc")[[1]][2])
Results$Target_Pathway_PCA$PC = sapply(Results$Target_Pathway_PCA$Criteria.Identifier, function(x) str_split(x,"pc")[[1]][2])

treatment = readRDS(get_workflow_outputs("wf-e3e631e2c9"))
treatment$Target.Identifier = signatureMapping$New_identifier[match(treatment$Target.ID, signatureMapping$ID)]
overlap = readRDS(get_workflow_outputs("wf-3c7ac8ddd1"))
coverage = readRDS(get_workflow_outputs("wf-4c7ecb1fc6"))
whiteSpace_coverage = readRDS(get_workflow_outputs("wf-97923f7e39"))

Results_complementarity = bind_rows(treatment, overlap, coverage, whiteSpace_coverage) %>%
  dplyr::filter(Type == "bulk") %>%
  dplyr::filter(Target.Identifier %in% keep_signatures$Target[which(keep_signatures$Ranking == "Yes")]) %>%
  mutate(DataType = case_when(Criteria.Identifier == "Complementarity to Dupilumab Treatment W16" ~ "Complementarity to W16 Dupilumab",
                              Criteria.Identifier == "Complementarity to Dupilumab Treatment W4" ~ "Complementarity to W4 Dupilumab",
                              Criteria.Identifier == "Complementarity to Dupilumab Non Responders W16" ~ "Complementarity to W16 Non Responders",
                              Criteria.Identifier == "Complementarity to Dupilumab Non Responders W4" ~ "Complementarity to W4 Non Responders",
                              Criteria.Identifier == "Coverage of Disease Features" ~ "Disease Coverage",
                              Criteria.Identifier == "Shared Disease Features" ~ "Overlap with IL13",
                              Criteria.Identifier == "Coverage of Dupilumab White Space" ~ "White Space Coverage"))

Results_binded = Results[setdiff(names(Results), c("Target_Cell", "Target_Pathway", "Target_Gene"))] %>%
  bind_rows() %>%
  mutate(DataType = ifelse(DataType == "Target_CS", paste0("Target_CS_",Criteria.Identifier), DataType)) %>%
  mutate(DataType = ifelse(DataType == "Target_Cell_PCA", paste0("Target_Cell_PC",PC), DataType)) %>%
  mutate(DataType = ifelse(DataType == "Target_Pathway_PCA", paste0("Target_Pathway_PC",PC), DataType)) %>%
  mutate(DataType = str_replace_all(DataType, " ", "_"))

# remove bulk & adj mix
Results_binded = Results_binded[-which(str_detect(Results_binded$DataType,"Target_Pathway_PC") & 
                                         Results_binded$Type == "bulk" & 
                                         str_detect(Results_binded$Criteria.Identifier,"_adj")),]
Results_binded = Results_binded[-which(str_detect(Results_binded$DataType,"Target_Pathway_PC") & 
                                         Results_binded$Type == "adjusted__1__1" & 
                                         str_detect(Results_binded$Criteria.Identifier,"_bulk")),]
Results_binded$DataType[which(Results_binded$Type != "bulk")] <- paste0("Adj_", Results_binded$DataType[which(Results_binded$Type != "bulk")])

Results_binded = Results_binded %>%
  dplyr::filter(!str_detect(DataType, "PC3")) %>%
  mutate(DataType = str_remove(DataType,"_bulk|_adj")) %>%
  dplyr::filter(Target.Collection != "epidermis") %>% # all pathways already exist
  dplyr::filter(!str_detect(Criteria.Identifier, "neuroinflammation|epidermis|th2"))

all_results_binded = bind_rows(Results_binded, Results_additional, Results_complementarity)

Results.corr = reshape2::dcast(all_results_binded, Target.ID ~ DataType, value.var = "metricValue") %>%
  column_to_rownames(var = "Target.ID") %>%
  as.matrix

Results.corr = Results.corr[complete.cases(Results.corr),]

corr_corr = cor(Results.corr, method = "pearson")

png("~/analysis-s05/figures/Results/criteria_correlation_pearson.png", bg = "white", res = "120", width = 1600, height = 1200)
Heatmap(corr_corr, cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(round(corr_corr[i, j],1), x, y, gp = gpar(fontsize = 10))},
  name = "Pearson\nCorrelation", row_labels = str_replace_all(rownames(corr_corr),"_"," "), 
  column_labels = str_replace_all(rownames(corr_corr),"_"," "),
  cluster_rows = F, cluster_columns = F)
dev.off()

only_0.9 = rownames(corr_corr)[apply(corr_corr > 0.9 & row(corr_corr) != col(corr_corr), 1, any)]
only_0.9 <- corr_corr[only_0.9, only_0.9, drop = FALSE]


png("~/analysis-s05/figures/Results/criteria_correlation_above0.9_pearson.png", bg = "white", res = "120", width = 1200, height = 1000)
Heatmap(only_0.9, cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(round(only_0.9[i, j],1), x, y, gp = gpar(fontsize = 10))},
  name = "Pearson\nCorrelation", row_labels = str_replace_all(rownames(only_0.9),"_"," "), 
  column_labels = str_replace_all(rownames(only_0.9),"_"," "))
dev.off()

less_0.9 = rownames(corr_corr)[apply(corr_corr < 0.9 & row(corr_corr) != col(corr_corr), 1, any)]
less_0.9 = corr_corr[less_0.9, less_0.9, drop = FALSE]

png("~/analysis-s05/figures/Results/criteria_correlation_below0.9_pearson.png", bg = "white", res = "120", width = 1400, height = 1200)
Heatmap(less_0.9, cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(round(less_0.9[i, j],1), x, y, gp = gpar(fontsize = 10))},
  name = "Pearson\nCorrelation", row_labels = str_replace_all(rownames(less_0.9),"_"," "), 
  column_labels = str_replace_all(rownames(less_0.9),"_"," "))
dev.off()
