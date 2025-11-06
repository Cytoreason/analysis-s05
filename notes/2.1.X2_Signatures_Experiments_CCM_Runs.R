# In this part we genereate the signatures for X2 using RNAseq data from the client
# that was processed and QC'ed in the validator.

library(cytoreason.ccm.pipeline)
library(cytoreason.cc.client)
library(tidyverse)
library(reshape2)
devtools::load_all("~/analysis-s05/R/utils.R")

# 1. 4hr experiment
# ----------------------
X2_4hr = readRDS(get_workflow_outputs("wf-cb183f8dad"))
metadata = read.csv("~/analysis-s05/data/metadata_4h.csv")
metadata$sample[which(is.na(metadata$inhibitor) & metadata$sample == "CST14-756D3")] <- "CST14D3" # a mistake in the table
metadata$agonist = str_remove(metadata$agonist,"-")

# cbind four designs:
# 1. per agonist - evo vs uninhibited [paired]
# 2. per agonist - uninhibited samples - activated vs un-activated [paired]
# 3. all inhibited vs all uninhibited [covariate - stimulation, paired]
# 4. all activated non-inhibited vs all un-activated non-inhibited [paired]
designMat = cbind(metadata %>%
                    mutate(inhibitor_present = ifelse(is.na(inhibitor), "Uninhibited", "EVO576")) %>%
                    dcast(formula = sample ~ agonist, value.var = "inhibitor_present") %>%
                    arrange(match(sample, metadata$sample)) %>%  # restore original sample order
                    rename_with(.fn = ~ paste0(.x, "_inhibited_vs_uninhibited"), .cols = -1),
                  metadata %>%
                    mutate(value = ifelse(is.na(inhibitor), agonist, NA)) %>%
                    dcast(sample ~ agonist, value.var = "value", fill = NA) %>%
                    arrange(match(sample, metadata$sample)) %>%  # restore original sample order
                    mutate(across(-sample, ~ ifelse(!is.na(Untreated), "Unactivated", .))) %>%
                    dplyr::select(-Untreated, -sample) %>%
                    rename_with(.fn = ~ paste0(.x, "_activated_vs_unactivated")),
                  metadata %>%
                    mutate(inhibited_vs_uninhibited = ifelse(is.na(inhibitor),"Uninhibited","Inhibited"))%>%
                    dplyr::select(inhibited_vs_uninhibited),
                  metadata %>%
                    mutate(activated_vs_unactivated = ifelse(agonist == "Untreated", "Unactivated", "Activated")) %>%
                    mutate(activated_vs_unactivated = ifelse(!is.na(inhibitor), NA, activated_vs_unactivated)) %>%
                    dplyr::select(activated_vs_unactivated)
                  )

metadata = merge(metadata, designMat, by = 1)
metadata = merge(metadata, phenoData(X2_4hr)@data, by = 1)
rownames(metadata) = metadata$sample
X2_4hr = ExpressionSet(assayData = exprs(X2_4hr),
                       phenoData = AnnotatedDataFrame(metadata),
                       featureData = featureData(X2_4hr),
                       annotation = "org.Hs.eg.db")
pushToCC(X2_4hr, tagsToPass = list(list(name="object",value="eset_X2_4hr")))
# wf-09b09047b2



# 1. 24hr experiment
# ----------------------
X2_24hr = readRDS(get_workflow_outputs("wf-c973ccb1f5"))
metadata = read.csv("~/analysis-s05/data/metadata_24h.csv")
metadata = metadata %>%
  mutate(stimulation = str_remove(stimulation,"-")) %>%
  dplyr::filter(paired_end_read == 1) %>%
  mutate(sample = str_remove(sample, "_R1")) %>%
  mutate(sample = str_replace(sample,"D1|D2|D3",paste0("D",donor)))

# cbind four designs:
# 1. per agonist - evo vs uninhibited [paired]
# 2. per agonist - uninhibited samples - activated vs un-activated [paired]
# 3. all inhibited vs all uninhibited [covariate - stimulation, paired]
# 4. all activated non-inhibited vs all un-activated non-inhibited [covariate - stimulation, paired]
designMat = cbind(metadata %>%
                    mutate(inhibitor_present = ifelse(is.na(inhibitor), "Uninhibited", "EVO576")) %>%
                    dcast(formula = sample ~ stimulation, value.var = "inhibitor_present") %>%
                    arrange(match(sample, metadata$sample)) %>%  # restore original sample order
                    rename_with(.fn = ~ paste0(.x, "_inhibited_vs_uninhibited"), .cols = -1) %>%
                    dplyr::select(-Untreated_inhibited_vs_uninhibited), # no un-activated samples with inhibition
                  metadata %>%
                    mutate(value = ifelse(is.na(inhibitor), stimulation, NA)) %>%
                    dcast(sample ~ stimulation, value.var = "value", fill = NA) %>%
                    arrange(match(sample, metadata$sample)) %>%  # restore original sample order
                    mutate(across(-sample, ~ ifelse(!is.na(Untreated), "Unactivated", .))) %>%
                    dplyr::select(-Untreated, -sample) %>%
                    rename_with(.fn = ~ paste0(.x, "_activated_vs_unactivated")),
                  metadata %>%
                    mutate(inhibited_vs_uninhibited = ifelse(stimulation=="Untreated", NA, ifelse(is.na(inhibitor),"Uninhibited","Inhibited"))) %>%
                    dplyr::select(inhibited_vs_uninhibited),
                  metadata %>%
                    mutate(activated_vs_unactivated = ifelse(stimulation == "Untreated", "Unactivated", "Activated")) %>%
                    mutate(activated_vs_unactivated = ifelse(!is.na(inhibitor), NA, activated_vs_unactivated)) %>%
                    dplyr::select(activated_vs_unactivated)
                  )

metadata = merge(metadata, designMat, by = 1)
metadata = merge(metadata,phenoData(X2_24hr)@data, by.y = "fileID", by.x = "sample", all = T)
rownames(metadata) = metadata$sample

X2_24hr = ExpressionSet(assayData = exprs(X2_24hr), # notice we take only the exprs, otherwise it calculates on the counts
                       phenoData = AnnotatedDataFrame(metadata),
                       featureData = featureData(X2_24hr),
                       annotation = "org.Hs.eg.db")
pushToCC(X2_24hr, tagsToPass = list(list(name="object",value="eset_X2_24hr")))
# wf-84a90d45bb



# Constructing a config file & running ccm
# ---------------------------------------------
meta_4hr = readRDS(get_workflow_outputs("wf-09b09047b2"))@phenoData@data
meta_24hr = readRDS(get_workflow_outputs("wf-84a90d45bb"))@phenoData@data
write.csv(meta_4hr, "~/analysis-s05/data/metadata_final_4h.csv")
write.csv(meta_24hr, "~/analysis-s05/data/metadata_final_24h.csv")

configLine = data.frame(ccm_exclude = NA,
                        ccm_meta = 0,
                        ccm_meta_pca = 0,
                        ccm_annotate = "all",
                        experiment_id = NA,
                        platform_id = NA,
                        asset_id = NA,
                        asset_data_access = "s05",
                        effect_id = "Treatment",
                        model_type = "Treatment",
                        comparison = NA,
                        comparison_id = NA,
                        comparison_levels = NA,
                        model_group = NA,
                        model_pairing = "donor",
                        model_covariates = NA,
                        model_timepoint = NA,
                        model_exclude_samples = NA,
                        condition = "unknown",
                        tissue = "skin",
                        sample_tissue = NA,
                        comments = "mast cells cell line")

config = do.call(rbind, replicate(19, configLine, simplify = F))

config$experiment_id[1:12] = "X2_4hr"
config$experiment_id[13:19] = "X2_24hr"

config$platform_id[1:12] = "GPL24676"
config$platform_id[13:19] = "GPL34281"

config$asset_id[1:12] = "ccw://wf-09b09047b2:0:output.rds"
config$asset_id[13:19] = "ccw://wf-84a90d45bb:0:output.rds"

config$comparison[1:9] = sapply(colnames(meta_4hr)[7:15], function(x) str_split_fixed(x,"_",n = 2)[[2]])
config$comparison[10:12] = c(colnames(meta_4hr)[16:17],paste0(colnames(meta_4hr)[17],"_cov"))
config$comparison[13:16] = sapply(colnames(meta_24hr)[8:11], function(x) str_split_fixed(x,"_",n = 2)[[2]])
config$comparison[17:19] = c(colnames(meta_24hr)[12:13],paste0(colnames(meta_24hr)[13],"_cov"))
config$comparison_id[1:12] = c(paste0(colnames(meta_4hr)[7:17],"_4hr"),paste0(colnames(meta_4hr)[17],"_cov_4hr"))
config$comparison_id[13:19] = c(paste0(colnames(meta_24hr)[8:13],"_24hr"),paste0(colnames(meta_24hr)[13],"_cov_24hr"))

config$model_covariates[12] = "agonist"
config$model_covariates[19] = "stimulation"

config$model_group[c(1:5,13:14)] = sapply(c(colnames(meta_4hr)[7:11], colnames(meta_24hr)[8:9]),
                                          function(x) paste0(x,":{uninhibited:Uninhibited, inhibited:EVO576}"))
config$model_group[c(6:9,15:16)] = sapply(c(colnames(meta_4hr)[12:15], colnames(meta_24hr)[10:11]),
                                           function(x) paste0(x,":{unactivated:Unactivated",", ", "activated:",str_split(x,"_")[[1]][1],"}"))
config$model_group[c(10,17)] = "inhibited_vs_uninhibited:{uninhibited:Uninhibited, inhibited:Inhibited}"
config$model_group[c(11,12,18,19)] = "activated_vs_unactivated:{unactivated:Unactivated, activated:Activated}"

write.csv(config, "~/analysis-s05/data/X2_config.csv", row.names = F)
pushToCC(config, tagsToPass = list(list(name="object",value="config"),
                                   list(name="project",value="X2_Signatures")))
# wf-a399720b79


# local test first
ccm_stage_prepare_dataset_collection(config) # you need to have a folder called 'output' in the current wd

ccm_api_generate_ccm(config = config,
                     adjustment_models = list(),
                     stages = .skip("meta", "dataset2","dataset"),
                     model = .skip('gene_cell_correlations','cell_cell_correlations','feature_cell_correlations','pheno_feature_correlations',
                                   'cell_specific_differences','gene_set_activity_differences','survival_analysis_tme'),
                     image = "master@0.72.0")
# generate_ccm -- Thu Oct 30 14:28:43 2025: wf-2826ed2d24 [] - counts
# generate_ccm -- Thu Oct 30 14:05:43 2025: wf-2a79b1167b [] - exprs


## Processing the results
## ------------------------------------------
ccm = as_ccm_fit("wf-2826ed2d24")

allOptions = lapply(names(ccm$datasets), function(dataset){
  data.frame(experiment = dataset,
             comparison = names(ccm$datasets[[dataset]]$model))
}) %>% do.call(rbind,.)
allOptions = allOptions[-which(str_detect(allOptions$comparison,"_cov")),]

gx = list()
gsa = list()
apply(allOptions, 1, function(x){
  cat("\n",x,"\n")
  gx_tmp = statistic_table(analysisResultElement(ccm$datasets[[x[1]]]$model[[x[2]]],"gx_diff"))
  gx_tmp = gx_tmp[which(gx_tmp$term %in% c("unactivated","activated_vs_unactivated","activated_vs_unactivated_cov","inhibited_vs_uninhibited")),]
  agonist = strsplit(x[2],"_")[[1]][1]
  if(agonist %in% c("inhibited","activated")) { agonist = "All" }
  if(agonist %in% c("SP","All")) { agonist = paste0(agonist,"_", strsplit(x[1],"_")[[1]][2]) }

  gx_tmp = data.frame(dataset_id = x[1], comparison = x[2], agonist = agonist, gx_tmp)
  gx <<- append(gx, list(gx_tmp))

  gsa_tmp = statistic_table(analysisResultElement(ccm$datasets[[x[1]]]$model[[x[2]]],"gx_gsa"))
  gsa_tmp = gsa_tmp[which(gsa_tmp$term %in% c("unactivated","activated_vs_unactivated","activated_vs_unactivated_cov","inhibited_vs_uninhibited")),]
  gsa_tmp = data.frame(dataset_id = x[1], comparison = x[2], agonist = agonist, gsa_tmp)
  gsa <<- append(gsa, list(gsa_tmp))
  NULL
})

gx_binded = do.call(rbind, gx)
gsa_binded = do.call(rbind, gsa)
uploadToBQ(gx_binded, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures_gxdiff") # wf-edcade5c80
uploadToBQ(gsa_binded, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures_gxgsea") # wf-9429788f8a

gx_binded$term = droplevels(gx_binded$term)

# Graph to see if there's a signal
ggplot(gx_binded, aes(x = estimate, y = log10_fdr)) +
  geom_point(aes(color = ifelse(log10_fdr >= -log10(0.05),"black","grey80"))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_identity() +
  facet_grid(cols = vars(agonist), rows = vars(term)) +
  theme_light() +
  theme(axis.title.x.top = element_text("Agonist"), axis.title.y.right = element_text("Comparison"))


# Clean graph for Evommune
# ------------------------------
library(patchwork)
gx = readRDS(get_workflow_outputs("wf-edcade5c80"))
gx = gx[-which(gx$term == "unactivated"),]
gx$agonist[str_detect(gx$agonist, "All")] <- "Pooled"
gx$agonist[str_detect(gx$agonist, "SP")] <- "SP"

x2_4h = gx[which(gx$dataset_id == "X2_4hr__GPL24676"),]
x2_24h = gx[which(gx$dataset_id == "X2_24hr__GPL34281"),]

x2_4h$agonist = factor(x2_4h$agonist, ordered = T, levels = c("Pooled","CST14","Icatibant","PAMP12","SP","Untreated"))
x2_24h$agonist = factor(x2_24h$agonist, ordered = T, levels = c("Pooled","aIgE","SP","Untreated"))

x2_4h = rbind(x2_4h, x2_4h[1,] %>%  mutate(agonist = "Untreated", estimate = NA))
x2_24h = rbind(x2_24h, rbind(x2_24h[1,] %>%  mutate(agonist = "Untreated", estimate = NA),
                            x2_24h[59293,] %>%  mutate(agonist = "Untreated", estimate = NA)))

FDRlabel <- x2_4h %>%
  group_by(term) %>%
  dplyr::filter(agonist == "Untreated") %>%  # Select the leftmost column for each row
  slice(1) # Just one row per facet

nSig = x2_4h %>%
  dplyr::filter(fdr <= 0.05) %>%  # Select the leftmost column for each row
  group_by(comparison, term) %>%
  mutate(nSigGenes = n()) %>%
  slice(1) # Just one row per facet

four = ggplot(x2_4h, aes(x = estimate, y = log10_fdr)) +
  geom_point(aes(color = ifelse(log10_fdr >= -log10(0.05),"black","grey80"))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text(data = FDRlabel, x = 6, y = 1.35, label = "FDR = 0.05", check_overlap = T, size = 3, hjust= "right") +
  geom_text(data = nSig, x = 3, y=2.2, aes(label=paste0("n=",nSigGenes)), check_overlap = T, size = 3, hjust = "left")+
  scale_color_identity() +
  scale_x_continuous(name = "estimate (log2 Fold Change)", sec.axis = dup_axis(name = "Agonist"))+
  scale_y_continuous(name = "-log10(FDR)", sec.axis = dup_axis(name = "Comparison"))+
  facet_grid(cols = vars(agonist), rows = vars(term)) +
  theme_bw() +
  theme(axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(), axis.text.y.right = element_blank()) +
  labs(title = expression(underline("4 hour experiment")))
# four  

FDRlabel <- x2_24h %>%
  group_by(term) %>%
  dplyr::filter(agonist == "Untreated") %>%  # Select the leftmost column for each row
  slice(1) # Just one row per facet

nSig = x2_24h %>%
  dplyr::filter(fdr <= 0.05) %>%  # Select the leftmost column for each row
  group_by(comparison, term) %>%
  mutate(nSigGenes = n()) %>%
  slice(1) # Just one row per facet

twentyfour = ggplot(x2_24h, aes(x = estimate, y = log10_fdr)) +
  geom_point(aes(color = ifelse(log10_fdr >= -log10(0.05),"black","grey80"))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text(data = FDRlabel, x = 7, y = 1.35, label = "FDR = 0.05", check_overlap = T, size = 3, hjust= "right") +
  geom_text(data = nSig, x = 5, y=1.8, aes(label=paste0("n=",nSigGenes)), check_overlap = T, size = 3, hjust = "left")+
  scale_color_identity() +
  scale_x_continuous(name = "estimate (log2 Fold Change)", sec.axis = dup_axis(name = "Agonist"))+
  scale_y_continuous(name = "-log10(FDR)", sec.axis = dup_axis(name = "Comparison"))+
  facet_grid(cols = vars(agonist), rows = vars(term)) +
  theme_bw() +
  theme(axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(), axis.text.y.right = element_blank()) +
  labs(title = expression(underline("24 hour experiment")))
# twentyfour

four/twentyfour

# clean envir
rm(twentyfour, four, FDRlabel, nSig, x2_24h, x2_4h, gx)
