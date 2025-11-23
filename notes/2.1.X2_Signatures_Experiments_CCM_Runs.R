# In this part we genereate the signatures for X2 using RNAseq data from the client
# that was processed and QC'ed in the validator.
library(cytoreason.ccm.pipeline)
library(tidyverse)
library(reshape2)
library(patchwork)
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
                    dplyr::select(activated_vs_unactivated),
                  metadata %>%
                    mutate(inhibited_vs_activated = ifelse(is.na(inhibitor), "Activated","EVO576")) %>%
                    mutate(inhibited_vs_activated = ifelse(agonist == "Untreated", NA, inhibited_vs_activated)) %>%
                    dplyr::select(inhibited_vs_activated)
                  )

metadata = merge(metadata, designMat, by = 1)
metadata = merge(metadata, phenoData(X2_4hr)@data, by = 1)
rownames(metadata) = metadata$sample
X2_4hr = ExpressionSet(assayData = assayData(X2_4hr),
                       phenoData = AnnotatedDataFrame(metadata),
                       featureData = featureData(X2_4hr),
                       annotation = "org.Hs.eg.db")
pushToCC(X2_4hr, tagsToPass = list(list(name="object",value="eset_X2_4hr")))
# wf-74008e33fd - counts
# wf-faf98d72b3 - exprs


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
                    dplyr::select(activated_vs_unactivated),
                  metadata %>%
                    mutate(inhibited_vs_activated = ifelse(is.na(inhibitor), "Activated","EVO576")) %>%
                    mutate(inhibited_vs_activated = ifelse(stimulation == "Untreated", NA, inhibited_vs_activated)) %>%
                    dplyr::select(inhibited_vs_activated)
                  )

metadata = merge(metadata, designMat, by = 1)
metadata = merge(metadata,phenoData(X2_24hr)@data, by.y = "fileID", by.x = "sample", all = T)
rownames(metadata) = metadata$sample

X2_24hr = ExpressionSet(assayData = assayData(X2_24hr), 
                       phenoData = AnnotatedDataFrame(metadata),
                       featureData = featureData(X2_24hr),
                       annotation = "org.Hs.eg.db")
pushToCC(X2_24hr, tagsToPass = list(list(name="object",value="eset_X2_24hr")))
# wf-b1d08dc2f1 - counts
# wf-fee413401f - exprs



# Constructing a config file & running ccm
# ---------------------------------------------
meta_4hr = readRDS(get_workflow_outputs("wf-74008e33fd"))@phenoData@data
meta_24hr = readRDS(get_workflow_outputs("wf-b1d08dc2f1"))@phenoData@data
write.csv(meta_4hr, "~/analysis-s05/data/metadata_final_4h.csv")
write.csv(meta_24hr, "~/analysis-s05/data/metadata_final_24h.csv")

comparisonName = function(colname){
  ifelse(str_detect(colname,"_activated|_inhibited") & !str_detect(colname, "inhibited_vs_activated"), str_split_fixed(colname,"_",n = 2)[[2]], colname)
}

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
                        model_covariates = NA, # we are not using covariates to the pooled analysis because we want to see the pooled effect and not account for differences between agonists
                        model_timepoint = NA,
                        model_exclude_samples = NA,
                        condition = "unknown",
                        tissue = "skin",
                        sample_tissue = NA,
                        comments = "mast cells cell line")

idx_meta_4hr = 7:(ncol(meta_4hr)-2) # starts with column 7
idx_meta_24hr = 8:(ncol(meta_24hr)-2) # starts in column 8
idx_4hr  = seq_len(ncol(meta_4hr) + 3 - 8) # add 3 for covariate analysis
idx_24hr = seq(from = length(idx_4hr) + 1,
              to = length(idx_4hr) + (ncol(meta_24hr) + 3 - 9))  # add 3 for covariate analysis
names(idx_4hr) = c(colnames(meta_4hr[idx_meta_4hr]), paste0(colnames(meta_4hr[tail(idx_meta_4hr,3)]),"_cov"))
names(idx_24hr) = c(colnames(meta_24hr[idx_meta_24hr]), paste0(colnames(meta_24hr[tail(idx_meta_24hr,3)]),"_cov"))

config = do.call(rbind, replicate(tail(idx_24hr,1), configLine, simplify = F))

config$experiment_id[idx_4hr] = "X2_4hr"
config$experiment_id[idx_24hr] = "X2_24hr"

config$platform_id[idx_4hr] = "GPL24676"
config$platform_id[idx_24hr] = "GPL34281"

config$asset_id[idx_4hr] = "ccw://wf-74008e33fd:0:output.rds"
config$asset_id[idx_24hr] = "ccw://wf-b1d08dc2f1:0:output.rds"

config$comparison[idx_4hr] = sapply(names(idx_4hr), comparisonName)
config$comparison[idx_24hr] = sapply(names(idx_24hr), comparisonName)

config$comparison_id[idx_4hr] = names(idx_4hr)
config$comparison_id[idx_24hr] = names(idx_24hr)

config$model_covariates[str_detect(config$comparison_id,"cov_4hr")] = "agonist"
config$model_covariates[str_detect(config$comparison_id,"cov_24hr")] = "stimulation"

config$model_group[str_detect(config$comparison_id, "_inhibited_vs_uninhibited")] =
  sapply(config$comparison_id[str_detect(config$comparison_id, "_inhibited_vs_uninhibited")],
         function(x) paste0(x,":{uninhibited:Uninhibited, inhibited:EVO576}"))
config$model_group[str_detect(config$comparison_id, "_activated_vs_unactivated")] =
  sapply(config$comparison_id[str_detect(config$comparison_id, "_activated_vs_unactivated")],
         function(x) paste0(x,":{unactivated:Unactivated",", ", "activated:",str_split(x,"_")[[1]][1],"}"))
config$model_group[str_detect(config$comparison_id, "^inhibited_vs_uninhibited")] = "inhibited_vs_uninhibited:{uninhibited:Uninhibited, inhibited:Inhibited}"
config$model_group[str_detect(config$comparison_id, "^activated_vs_unactivated")] = "activated_vs_unactivated:{unactivated:Unactivated, activated:Activated}"
config$model_group[str_detect(config$comparison_id, "inhibited_vs_activated")] = "inhibited_vs_activated:{activated:Activated, inhibited:EVO576}"

config$comparison_id[idx_4hr] = paste0(names(idx_4hr),"_4hr") # add the hr identifier
config$comparison_id[idx_24hr] = paste0(names(idx_24hr),"_24hr")

write.csv(config, "~/analysis-s05/data/X2_config.csv", row.names = F)
pushToCC(config, tagsToPass = list(list(name="object",value="config"),
                                   list(name="project",value="X2_Signatures")))
# wf-2fa2e07dd0 - counts
# wf-f591b7e333 - exprs


# local test first
ccm_stage_prepare_dataset_collection(config) # you need to have a folder called 'output' in the current wd

ccm_api_generate_ccm(config = config,
                     adjustment_models = list(),
                     stages = .skip("meta", "dataset2","dataset"),
                     model = .skip('gene_cell_correlations','cell_cell_correlations','feature_cell_correlations','pheno_feature_correlations',
                                   'cell_specific_differences','gene_set_activity_differences','survival_analysis_tme'),
                     image = "master@0.72.0")
# generate_ccm -- Sun Nov 23 12:28:26 2025: wf-b3318c6d1e [] - exprs
# generate_ccm -- Sun Nov 23 15:09:46 2025: wf-696ca67797 [] - counts



## Processing the results
## ------------------------------------------
ccm = as_ccm_fit("wf-696ca67797")

allOptions = lapply(names(ccm$datasets), function(dataset){
  data.frame(experiment = dataset,
             comparison = names(ccm$datasets[[dataset]]$model))
}) %>% do.call(rbind,.)

gx = list()
gsa = list()
apply(allOptions, 1, function(x){
  cat("\n",x,"\n")
  gx_tmp = statistic_table(analysisResultElement(ccm$datasets[[x[1]]]$model[[x[2]]],"gx_diff"))
  gx_tmp = gx_tmp[str_detect(gx_tmp$term, c("activated_vs_unactivated|inhibited_vs_uninhibited|inhibited_vs_activated")),]
  agonist = strsplit(x[2],"_")[[1]][1]
  if(agonist %in% c("inhibited","activated")) { agonist = "All" }

  gx_tmp = data.frame(dataset_id = x[1], comparison = x[2], agonist = agonist, gx_tmp)
  gx <<- append(gx, list(gx_tmp))

  gsa_tmp = statistic_table(analysisResultElement(ccm$datasets[[x[1]]]$model[[x[2]]],"gx_gsa"))
  gsa_tmp = gsa_tmp[str_detect(gsa_tmp$term, c("activated_vs_unactivated|inhibited_vs_uninhibited|inhibited_vs_activated")),]
  gsa_tmp = data.frame(dataset_id = x[1], comparison = x[2], agonist = agonist, gsa_tmp)
  gsa <<- append(gsa, list(gsa_tmp))
  NULL
})

gx_binded = do.call(rbind, gx)
gsa_binded = do.call(rbind, gsa)
uploadToBQ(gx_binded, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures_gxdiff")
pushToCC(gx_binded, tagsToPass = list(list(name="analysis",value="gx_diff"),list(name="project",value="X2")))
# wf-12254175c7 - exprs
# wf-fe0c7701a0 - counts

uploadToBQ(gsa_binded, bqdataset = "s05_atopic_dermatitis", tableName = "X2Signatures_gxgsea")
pushToCC(gsa_binded, tagsToPass = list(list(name="analysis",value="gx_gsa"),list(name="project",value="X2")))
# wf-9d65578ccc - exprs
# wf-e6a03c52ab - counts


# Clean graph for Evommune
# ------------------------------
gx = readRDS(get_workflow_outputs("wf-fe0c7701a0"))
gx$term = droplevels(gx$term)
gx$agonist[str_detect(gx$agonist, "All")] <- "Pooled"
gx$term = as.character(gx$term)

x2_4h = gx[which(gx$dataset_id == "X2_4hr__GPL24676"),]
x2_24h = gx[which(gx$dataset_id == "X2_24hr__GPL34281"),]

x2_4h$agonist[str_detect(x2_4h$comparison, "cov")] <- "Pooled+cov"
x2_24h$agonist[str_detect(x2_24h$comparison, "cov")] <- "Pooled+cov"

x2_4h$agonist = factor(x2_4h$agonist, ordered = T, levels = c("Pooled","Pooled+cov","CST14","Icatibant","PAMP12","SP","Untreated"))
x2_24h$agonist = factor(x2_24h$agonist, ordered = T, levels = c("Pooled","Pooled+cov","aIgE","SP","Untreated"))

# for visualization purposes
x2_4h = rbind(x2_4h,
              x2_4h[which(x2_4h$term == "activated_vs_unactivated")[1],] %>%  mutate(agonist = "Untreated", estimate = NA),
              x2_4h[which(x2_4h$term == "inhibited_vs_activated")[1],] %>%  mutate(agonist = "Untreated", estimate = NA))
x2_24h = rbind(x2_24h,
               x2_24h[which(x2_24h$term == "activated_vs_unactivated")[1],] %>%  mutate(agonist = "Untreated", estimate = NA),
               x2_24h[which(x2_24h$term == "inhibited_vs_uninhibited")[1],] %>%  mutate(agonist = "Untreated", estimate = NA),
               x2_24h[which(x2_24h$term == "inhibited_vs_activated")[1],] %>%  mutate(agonist = "Untreated", estimate = NA))

x2_4h$term = sapply(x2_4h$term, function(x) paste0(paste0(str_split(x,"_")[[1]][1:2],collapse = " "), "\n",str_split(x,"_")[[1]][3]))
x2_24h$term = sapply(x2_24h$term, function(x) paste0(paste0(str_split(x,"_")[[1]][1:2],collapse = " "), "\n",str_split(x,"_")[[1]][3]))

for(chosenMetric in c("fdr","pvalue")) {
  FDRlabel <- x2_4h %>%
    group_by(term) %>%
    dplyr::filter(agonist == "Untreated") %>%  # Select the leftmost column for each row
    dplyr::slice(1) # Just one row per facet
  
  nSig = x2_4h %>%
    mutate(metric = get(chosenMetric)) %>%
    mutate(log10_metric = -log10(metric)) %>%
    dplyr::filter(metric <= 0.05) %>%  # Select the leftmost column for each row
    group_by(comparison, term) %>%
    mutate(nSigGenes = n()) %>%
    dplyr::slice(1) # Just one row per facet
  
  dat = x2_4h %>%
    mutate(metric = get(chosenMetric)) %>%
    mutate(log10_metric = -log10(metric))
  
  four = ggplot(dat, aes(x = estimate, y = log10_metric)) +
    geom_point(aes(color = ifelse(log10_metric >= -log10(0.05),"black","grey80"))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text(data = FDRlabel, x = 4.5, y = 1.1, label = "pvalue = 0.05", check_overlap = T, size = 3, hjust= "left") +
    geom_text(data = nSig, x = 0, aes(y = max(log10_metric)-0.5, label=paste0("n=",nSigGenes)), check_overlap = T, size = 3, hjust = "left")+
    scale_color_identity() +
    scale_x_continuous(name = "estimate (log2 Fold Change)", sec.axis = dup_axis(name = "Agonist"))+
    scale_y_continuous(name = paste0("-log10(",chosenMetric,")"), sec.axis = dup_axis(name = "Comparison"))+
    facet_grid(cols = vars(agonist), rows = vars(term)) +
    theme_bw() +
    theme(axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(), axis.text.x.top = element_blank(),
          axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(), axis.text.y.right = element_blank()) +
    labs(title = expression(underline("4 hour experiment")))
  # four  
  
  FDRlabel <- x2_24h %>%
    group_by(term) %>%
    dplyr::filter(agonist == "Untreated") %>%  # Select the leftmost column for each row
    dplyr::slice(1) # Just one row per facet
  
  nSig = x2_24h %>%
    mutate(metric = get(chosenMetric)) %>%
    mutate(log10_metric = -log10(metric)) %>%
    dplyr::filter(metric <= 0.05) %>%  # Select the leftmost column for each row
    group_by(comparison, term) %>%
    mutate(nSigGenes = n()) %>%
    dplyr::slice(1) # Just one row per facet
  
  dat = x2_24h %>%
    mutate(metric = get(chosenMetric)) %>%
    mutate(log10_metric = -log10(metric))

  twentyfour = ggplot(dat, aes(x = estimate, y = log10_metric)) +
    geom_point(aes(color = ifelse(log10_metric >= -log10(0.05),"black","grey80"))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text(data = FDRlabel, x = 7, y = 1.1, label = "pvalue = 0.05", check_overlap = T, size = 3, hjust= "right") +
    geom_text(data = nSig, x = 4.5, aes(y = max(log10_metric)-0.5, label=paste0("n=",nSigGenes)), check_overlap = T, size = 3, hjust = "left")+
    scale_color_identity() +
    scale_x_continuous(name = "estimate (log2 Fold Change)", sec.axis = dup_axis(name = "Agonist"))+
    scale_y_continuous(name = paste0("-log10(",chosenMetric,")"), sec.axis = dup_axis(name = "Comparison"))+
    facet_grid(cols = vars(agonist), rows = vars(term)) +
    theme_bw() +
    theme(axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(), axis.text.x.top = element_blank(),
          axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(), axis.text.y.right = element_blank()) +
    labs(title = expression(underline("24 hour experiment")))
  # twentyfour
  
  graph = four/twentyfour
  ggsave(plot = graph, filename = paste0("~/analysis-s05/figures/X2_Signature/Experiments_Signal_",chosenMetric,".png"), 
         width = 800, height = 900, units = "px", dpi = 96)
  
}

# clean envir
rm(dat, twentyfour, four, FDRlabel, nSig, x2_24h, x2_4h, gx)
