devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(ComplexHeatmap)
library(tidyverse)

Results = readRDS(get_workflow_outputs("wf-a53211cea7"))
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
                                        str_detect(Results_binded$Criteria.Identifier,"adj_")),]
Results_binded = Results_binded[-which(str_detect(Results_binded$DataType,"Target_Pathway_PC") & 
                                         Results_binded$Type == "adjusted__1__1" & 
                                         str_detect(Results_binded$Criteria.Identifier,"(?<!_)pathway")),]

Results_binded$DataType[which(Results_binded$Type != "bulk")] <-paste0("Adj_", Results_binded$DataType[which(Results_binded$Type != "bulk")])
Results_binded = Results_binded %>%
  dplyr::filter(!str_detect(DataType, "PC3")) %>%
  dplyr::filter(!str_detect(DataType, "Adj_Target_Pathway_PC2|Adj_Target_Cell_PC|Adj_Target_CS|Adj_Target_MS"))

Results.corr = reshape2::dcast(Results_binded, Target.Identifier ~ DataType, value.var = "metricValue") %>%
    column_to_rownames(var = "Target.Identifier") %>%
    as.matrix

corr_corr = cor(Results.corr, method = "pearson")

# presented version used wf-a53211cea7 results
png("~/analysis-s05/figures/Results/criteria_correlation_pearson.png", bg = "white", res = "120", width = 1200, height = 1000)
Heatmap(corr_corr, cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(round(corr_corr[i, j],1), x, y, gp = gpar(fontsize = 10))},
  name = "Pearson\nCorrelation", row_labels = str_replace_all(rownames(corr_corr),"_"," "), 
  column_labels = str_replace_all(rownames(corr_corr),"_"," "))
dev.off()

corr_corr = cor(Results.corr, method = "spearman")

png("~/analysis-s05/figures/Results/criteria_correlation_spearman.png", bg = "white", res = "120", width = 1200, height = 1000)
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
library(WGCNA)
ccm_wfid = "wf-8e948630d7"
ccm= as_ccm_fit(ccm_wfid)
con = designSampleGroupContrasts(AnalysisModel2Group(ccm, type = "L_vs_NL"))

signatureMapping = readRDS(get_workflow_outputs("wf-b520743a43"))
pairings = unified_metadata %>%
  dplyr::filter(dataset_id %in% c("GSE130588__GPL570","GSE59294__GPL570")) %>%
  dplyr::filter(time == "D0") %>%
  dplyr::filter(sample_classification %in% c("Lesion","Non Lesion"))

pathways = c("time","sample_classification","condition","experiment_id","EASI")
enrichment =  lapply(ccm$datasets[which(names(ccm$datasets) %in% c("GSE130588__GPL570","GSE59294__GPL570"))], function(d){
  a = pData(assayDataExpression(d))[,pathways] %>%
    rownames_to_column(., var = "sample_id")
  b = statistic_table(d$gene_set_activity, subset = ~(submodel == "bulk" & collection == "X2"))
  b = reshape2::dcast(b, sample_id ~ gene_set, value.var = "value")
  return(merge(a,b, by="sample_id"))
}) %>% do.call(rbind,.)

enrichment_baseline = enrichment %>%
  dplyr::filter(time %in% c("W0","Pre") & sample_classification == "Lesion" & condition == "atopic dermatitis")

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
           FDR_fmt = dplyr::case_when(
             FDR <= 0.1        ~ sprintf("%.1e", FDR),  # scientific for ≤ 0.1
             TRUE              ~ sprintf("%.3f", FDR)   # fixed decimal otherwise
           ),
           label = sprintf("bicor = %.2f, FDR = %s", r, FDR_fmt)
    )

bicor_grouped = lapply(unique(x2$experiment_id), function(experiment){
  lapply(colnames(x2)[7:10], function(sig){
    dat = x2[which(x2$experiment_id == experiment),]
    bp = WGCNA::bicorAndPvalue(dat$EASI, dat[,sig], use = "pairwise.complete.obs")
    data.frame(
      signature = sig,
      experiment_id = experiment,
      r = as.numeric(bp$bicor),
      p = as.numeric(bp$p)
    )
  }) %>% bind_rows() %>%
    mutate(FDR = p.adjust(p, method = "fdr"),
           FDR_fmt = dplyr::case_when(
             FDR <= 0.1        ~ sprintf("%.1e", FDR),  # scientific for ≤ 0.1
             TRUE              ~ sprintf("%.3f", FDR)   # fixed decimal otherwise
           ),
           label = sprintf("bicor = %.2f, FDR = %s", r, FDR_fmt)
    )
}) %>% do.call(rbind,.)


x2 = reshape2::melt(x2, id.vars = 1:6, variable.name = "signature", value.name = "sample_enrichment")
x2 = x2 %>% left_join(., bicor_tabl, by = "signature")

ggplot(x2, aes(x = EASI, y = sample_enrichment)) +
  geom_point(aes(color = experiment_id)) +
  geom_smooth(method="lm", se = F, color = "black") +
  geom_smooth(method="lm", se = F, aes(color = experiment_id)) +
  geom_text(x = Inf, y = Inf, aes(label = label), check_overlap = T, show.legend = F, hjust = 1.1, vjust = 2, size = 3)+
  scale_color_brewer(palette = "Set2")+
  facet_wrap(~signature, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  ggpubr::border()+
  labs(x = "EASI score", y = "Signature enrichment", color = NULL)



as_ccm_fit(AssetData(CCM_WF_1))
config_gs <-  modelMetadata(ccm_fit_gs)
#contrasts <- designSampleGroupContrastData(ccm_fit_gs)
res <- designSampleGroupContrasts(ccm_fit_gs)
uc_res <-  lapply(names(ccm_fit_gs$datasets), function(cur_dataset) {   message("Processing dataset: ", cur_dataset)      cur_dataset_df <- ccm_fit_gs$datasets[[cur_dataset]]   model_names <- names(cur_dataset_df$model)      res_df <- do.call('rbind', lapply(model_names, function(cur_model) {     message("  Model: ", cur_model)
res <- cur_dataset_df$model[[cur_model]]
fit <- analysisResultElement(res, method)
all_contrasts <- designSampleGroupContrasts(cur_model,
                                            data = ccm_fit_gs$datasets[[cur_dataset]])
if (length(all_contrasts) == 0) {       warning("    Skipping model ", cur_model, " due to no contrasts.")       return(NULL)     }
res <- tryCatch({       temp <- data.frame(ldply(all_contrasts, ddply, c("level", "label"),
                                                 plyr::summarize, n = length(label), .id = "term"),                          
                                           stringsAsFactors = FALSE)
temp$dataset_id <- cur_dataset       
temp$model <- cur_model      
temp     }, error = function(e) {       warning("    Error processing model ", cur_model, ": ", conditionMessage(e))       return(NULL)     })          return(res)   }))
return(res_df) })