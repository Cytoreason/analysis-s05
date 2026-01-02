devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline) # make sure to load version >= 1.0.1
library(tidyverse)

ccm_wfid = "wf-abde4bfab0"
keep_signatures = openxlsx::read.xlsx("~/analysis-s05/data/Final signatures to be ranked or viewed in the dashboards.xlsx")
signatureMapping = readRDS(get_workflow_outputs("wf-aa75ed069b"))

Results = readRDS(get_workflow_outputs("wf-47f2a1c1b7"))
Target_Pathway = Results$Target_Pathway %>%
  dplyr::filter(Criteria.Collection %in% c("epidermis","neuroinflammation","th2")) %>%
  dplyr::filter(Target.Identifier %in% keep_signatures$Target[which(keep_signatures$Exploration == "Yes")])

# Target_Pathway$effect = Target_Pathway$metricValue * Target_Pathway$log10_fdr


# Determining significance based on permutations
# ------------------------------------------------------
Permutations = Results$Target_Pathway %>%
  dplyr::filter(Target.Identifier %in% keep_signatures$Target[which(keep_signatures$Exploration == "Yes")])

# For epidermis, we want to score high the targets that correlate to the most downregulated pathways, so we want to flip the
# sign when we calculate the fdr. Additionally, we have two pathways where we want to take their opposite direction. We will change only their
# direction and calculate fdr as the rest.
pathwaysNotToFlipSign = c("epidermis:reactome__Asymmetric localization of PCP proteins", "epidermis:reactome__Differentiation of keratinocytes in interfollicular epidermis in mammalian skin")
idx = which(Permutations$Criteria.Collection == "epidermis" & !Permutations$Criteria.Identifier %in% pathwaysNotToFlipSign & Permutations$Type == "bulk")
Permutations$metricValue[idx] = (-1) * Permutations$metricValue[idx]

# step 1: scramble pathway collection (we don't care about the actual ID, only that we have the same number of pathways per collection)
# step 2: extract n pathways (respective to the collection we test)
# step 3: calculate median correlation
# repeat 1000 times
nPermutations = 1000
set.seed(1234)

scrambledPathways = Permutations[,c("Target.Identifier","Type","Criteria.Collection","metricValue")] 
scrambledPathways = lapply(1:nPermutations, function(i){
      perm = scrambledPathways %>%
        group_by(Type, Target.Identifier) %>%
        mutate(Criteria.Collection = sample(Criteria.Collection, replace = F)) %>%
        dplyr::filter(Criteria.Collection %in% c("epidermis","neuroinflammation","th2")) %>%
        mutate(perm = i)
      return(perm)
  })
scrambledPathways_long = bind_rows(scrambledPathways)

real = Permutations[,c("Target.Identifier","Type","Criteria.Collection","metricValue")] %>%
  dplyr::filter(Criteria.Collection %in% c("epidermis","neuroinflammation","th2")) %>%
  group_by(Type, Target.Identifier, Criteria.Collection) %>%
  summarise(median_corr = median(metricValue),
            mean_corr = mean(metricValue)) %>%
  ungroup()

sig = scrambledPathways_long %>%
  dplyr::filter(Criteria.Collection %in% c("epidermis","neuroinflammation","th2")) %>%
  group_by(Type, Target.Identifier, Criteria.Collection) %>%
  left_join(real, by = join_by(Type, Target.Identifier, Criteria.Collection), relationship = "many-to-one") %>%
  group_by(Type, Target.Identifier, Criteria.Collection, perm, median_corr, mean_corr) %>%
  summarise(median_perm = median(metricValue),
            mean_perm = mean(metricValue)) %>%
  group_by(Type, Target.Identifier, Criteria.Collection, median_corr, mean_corr) %>%
  summarise(pvalue_median = (sum(median_perm >= median_corr) + 1) / (nPermutations + 1),
            pvalue_mean = (sum(mean_perm >= mean_corr) + 1) / (nPermutations + 1)) %>% # to avoid zeros
  group_by(Type) %>%
  mutate(fdr_median = p.adjust(pvalue_median, method = "fdr"),
         fdr_mean = p.adjust(pvalue_mean, method = "fdr")) %>%
  ungroup()
  
  
additional_criteria = sig %>%
  rename(metricValue = mean_corr, fdr = fdr_mean) %>%
  mutate(log10_fdr = -log10(fdr)) %>%
  mutate(
    Criteria.Collection = str_to_title(Criteria.Collection),
    DataType = paste0("Target_",Criteria.Collection),
    Target.ID = signatureMapping$ID[match(Target.Identifier, signatureMapping$New_identifier)],
    Target.Collection = signatureMapping$collection[match(Target.Identifier, signatureMapping$New_identifier)],
    Criteria.Identifier = Criteria.Collection,
    metricType = "bi-weight mid-correlation",
    hit = NA
  ) %>%
  dplyr::select(-median_corr, -fdr_median, -pvalue_median, -pvalue_mean)

additional_criteria$Criteria.Identifier = str_replace(additional_criteria$Criteria.Identifier, "Epidermis","Epidermal Barrier Integrity") 
additional_criteria$Criteria.Identifier = str_replace(additional_criteria$Criteria.Identifier, "Th2","Type 2 Inflammation") 

pushToCC(additional_criteria, tagsToPass = list(list(name="object",value="additional_criteria")))
# wf-63160e254a
# wf-7a1d7b9ecb