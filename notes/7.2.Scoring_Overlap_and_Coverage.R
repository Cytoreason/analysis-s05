# 1. Libraries and data loading
# ------------------------------------------
devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline) # make sure to load version >= 1.0.1
library(tidyverse)
library(scales, include.only = "rescale")

keep_signatures = openxlsx::read.xlsx("~/analysis-s05/data/Final signatures to be ranked or viewed in the dashboards.xlsx")

scaling = function(scale_vector, min_val = 0.1, max_val = 1) {
  rescale(x = scale_vector, to = c(min_val, max_val))
}

pathways_in_dz = downloadFromBQ(bqdataset = "s05_atopic_dermatitis", tableName = "AD_gx_gsa")
pathways_in_dz = pathways_in_dz %>%
  group_by(term, submodel) %>%
  dplyr::filter(submodel %in% c("bulk","adjusted") & 
                  collection %in% c("h","kegg","reactome","btm","th2","neuroinflammation") &
                  term == "L_vs_HC") %>%
  mutate(dir = ifelse(FDR > 0.05, "unchanged",
                      ifelse(NES >0, "up","down"))) %>%
  ungroup() %>%
  dplyr::select(submodel, pathway, dir) %>%
  mutate(submodel = ifelse(submodel == "bulk", "bulk", "adjusted__1__1"))

Results = readRDS(get_workflow_outputs("wf-47f2a1c1b7"))
Results = Results$Target_Pathway
Results = Results[which(Results$Target.Identifier %in% keep_signatures$Target[which(keep_signatures$Exploration == "Yes")]),]
Results$Criteria.Identifier = sub("^.*:", "", Results$Criteria.Identifier)
Results = left_join(Results, pathways_in_dz, by = join_by("Type" == "submodel","Criteria.Identifier" == "pathway"))
Results = Results[-which(is.na(Results$dir) | Results$dir == "unchanged"),]
Results$Criteria.Identifier = paste0(Results$Criteria.Collection,":",Results$Criteria.Identifier)
Results = Results[!duplicated(Results),]

steady_arm = Results[which(Results$Target.Identifier == "IL13"),] %>%
  dplyr::select(Type, Target.Identifier, Criteria.Identifier, metricValue) %>%
  rename(correlation = metricValue) %>%
  pivot_wider(id_cols = c(Type, Criteria.Identifier), names_from = Target.Identifier, 
              names_glue = "correlation_{Target.Identifier}", values_from = correlation)
  
targets = Results[-which(Results$Target.Identifier == "IL13"),] %>%
  dplyr::select(Type, Target.Identifier, Criteria.Identifier, metricValue) %>%
  rename(correlation_target = metricValue)

joined = merge(targets, steady_arm, by = c("Type", "Criteria.Identifier"), all = T)
# wf-20430b1067


# Overlap is the correlation between IL13 and the targets
# -----------------------------------------------------------------
overlap_statistic = joined %>%
  group_by(Type, Target.Identifier) %>%
  summarise(
    metricValue  = WGCNA::bicorAndPvalue(correlation_target, correlation_IL13, use = "pairwise.complete.obs")$bicor[[1]],
    pvalue = WGCNA::bicorAndPvalue(correlation_target, correlation_IL13, use = "pairwise.complete.obs")$p[[1]],
  ) %>%
  mutate(pvalue = ifelse(pvalue == 0, 1e-295, pvalue)) %>%
  mutate(fdr = p.adjust(pvalue, method = "fdr"),
         log10_fdr = -log10(fdr)) %>%
  ungroup() %>%
  mutate(Criteria.Identifier = "Shared Disease Features", Criteria.Collection = "Dupilumab Complementarity") %>%
  mutate(Target.ID = signatureMapping$ID[match(Target.Identifier, signatureMapping$New_identifier)]) %>%
  mutate(Target.Collection = signatureMapping$collection[match(Target.Identifier, signatureMapping$New_identifier)]) 

uploadToBQ(overlap_statistic, bqdataset = "s05_atopic_dermatitis", tableName = "Results_overlap")
pushToCC(overlap_statistic, tagsToPass = list(list(name='object',value='overlap')))
# wf-9f2cf87096
# wf-8222c73fb8

# Coverage is A+B+C/A+B+C+D with a 0.4 cutoff
# ------------------------------------------------
cutoff = 0.4

coverage_statistic = joined %>%
  mutate(group_IL13 = case_when(correlation_target < cutoff & correlation_IL13 < cutoff ~ "ns", .default = "covering")) %>%
  group_by(Type, Target.Identifier) %>%
  summarise(coverage = length(which(group_IL13 == "covering"))/n(), .groups = "drop")


# # permutations for FDR - will permute on the target within each pathway
# nPermutations = 1000
# set.seed(1234)
# 
# perms = lapply(1:nPermutations, function(i) {
#   joined %>%
#     mutate(covering = !(correlation_target < cutoff & correlation_IL13 < cutoff)) %>%
#     group_by(Type) %>%
#     mutate(Target.Identifier = sample(Target.Identifier)) %>%
#     ungroup() %>%
#     group_by(Type, Target.Identifier) %>%
#     summarise(coverage_perm = mean(covering), .groups = "drop")
# }) %>% bind_rows()
# 
# apply(coverage_statistic, 1, function(r){
#   tmp = perms$coverage_perm[which(perms$Type == r[1] & perms$Target.Identifier == r[2])]
#   length(which(tmp >= as.numeric(r[3])))
# })
# coverage_statistic = coverage_statistic %>%
#   left_join(perms, by = c("Type", "Target.Identifier")) %>%
#   group_by(Type, Target.Identifier) %>%
#   summarise(
#     coverage = first(coverage),
#     pvalue   = (sum(coverage_perm >= first(coverage)) + 1) / (n() + 1),
#     .groups  = "drop"
#   )


coverage_statistic = coverage_statistic %>%
  mutate(Criteria.Identifier = "Coverage of Disease Features", Criteria.Collection = "Disease Coverage") %>%
  mutate(DataType = "Target_Coverage") %>%
  mutate(Target.ID = signatureMapping$ID[match(Target.Identifier, signatureMapping$New_identifier)]) %>%
  mutate(Target.Collection = signatureMapping$collection[match(Target.Identifier, signatureMapping$New_identifier)]) %>%
  rename(metricValue = coverage) %>%
  mutate(fdr = NA, log10_fdr = NA, hit = NA) %>%
  ungroup()

uploadToBQ(coverage_statistic, bqdataset = "s05_atopic_dermatitis", tableName = "Results_coverage")
pushToCC(coverage_statistic, tagsToPass = list(list(name='object',value='coverage')))
# wf-98cd7a5d9d
# wf-c4a1547e7a

# Coverage of white space
# ------------------------------------------------
joined = readRDS(get_workflow_outputs("wf-20430b1067"))
pathway_space = readRDS(get_workflow_outputs("wf-fb61fed813"))

white_space = pathway_space %>%
  dplyr::filter(term %in% c("W4_vs_W0:DupilumabL")) %>%
  mutate(pathway = paste0(collection,":",pathway)) %>%
  mutate(in_white_space = case_when(in_white_space_L_vs_HC %in% c("yesUp","yesDown") ~ TRUE, .default = FALSE)) %>%
  dplyr::select(pathway, term, in_white_space)

target_pathway = joined %>%
  mutate(collection = str_split_fixed(Criteria.Identifier, ":", 2)[, 1]) %>%
  dplyr::filter(collection %in% c("btm", "h", "kegg", "neuroinflammation", "reactome", "th2")) %>%
  left_join(white_space, by = join_by("Criteria.Identifier" == "pathway")) %>%
  dplyr::filter(in_white_space)

# calculation
whiteSpace_coverage = target_pathway %>%
  group_by(Target.Identifier, Type) %>%
  summarise(ws_coverage = mean(abs(correlation_target)>=0.4), .groups = "drop")

whiteSpace_coverage = whiteSpace_coverage %>%
  mutate(Criteria.Identifier = "Coverage of Dupilumab White Space") %>%
  mutate(Criteria.Collection = "Dupilumab Complementarity") %>%
  mutate(Target.ID = signatureMapping$ID[match(Target.Identifier, signatureMapping$New_identifier)]) %>%
  mutate(Target.Collection = signatureMapping$collection[match(Target.Identifier, signatureMapping$New_identifier)]) %>%
  rename(metricValue = ws_coverage) %>%
  mutate(metricType = "coverage") %>%
  mutate(fdr = NA, log10_fdr = NA, hit = NA, Type = "bulk") %>%
  ungroup()

uploadToBQ(whiteSpace_coverage, bqdataset = "s05_atopic_dermatitis", tableName = "whiteSpace_coverage")
pushToCC(whiteSpace_coverage, tagsToPass = list(list(name="object",value="whitespace_coverage")))
# wf-48c725a33f