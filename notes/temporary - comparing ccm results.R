devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)

old_ccm = as_ccm_fit("wf-08a6a0a503")
new_ccm = as_ccm_fit("wf-832ab799be")


## ct_test per dataset
## =============================
extract_ct_test = function(ccm) {
  lapply(ccm$datasets, function(d){
    lapply(d$model, function(model){
      res = statistic_table(analysisResultElement(model, "ct_test"))
    }) %>% bind_rows()
  }) %>% bind_rows()
}

ct_test_old = extract_ct_test(old_ccm)
  ct_test_old$wfid = "wf-08a6a0a503"
ct_test_new = extract_ct_test(new_ccm)
  ct_test_new$wfid = "wf-832ab799be"
  colnames(ct_test_new)[9] = "FDR"

ct_tests = bind_rows(ct_test_old, ct_test_new)
ct_tests = ct_tests %>%
  dplyr::filter(model %in% c("AD","L_vs_HC","L_vs_NL","NL_vs_HC")) %>%
  dplyr::filter(term %in% c("DZ_vs_HC", "L_vs_HC","L_vs_NL","NL_vs_HC"))

# compare number of samples included
# --------------------------------------
nSamples = unique(ct_tests[,c("wfid","term","experiment_id","n_1","n_2")])
nSamples = reshape2::melt(nSamples, id.vaes = 1:3, variable.name = "group", value.name = "nSamples")
nSamples$version = ifelse(nSamples$wfid == "wf-08a6a0a503", "original", "new")
nSamples$group = paste0(nSamples$version,"-",nSamples$group)
nSamples$group = factor(nSamples$group, ordered = T, levels = c("original-n_1","new-n_1","original-n_2","new-n_2"))

ggplot(nSamples, aes(x = term, y = nSamples, fill = group)) +
  geom_col(position = position_dodge(0.75)) +
  facet_wrap(~experiment_id, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x = NULL, y = "nSamples", title = "Number of samples participating in each comparison", subtitle = "n1 is the second listed group, e.g. HC")


# compare values
# --------------------------
library(plyr)
skin = read.csv("~/data/skin_modified.csv", header = F, sep = "\t")

ct_tests = ct_tests %>%
  mutate(feature_id = case_when(wfid != "wf-08a6a0a503" ~ skin$V1[match(feature_id, skin$V2)],
                                .default = feature_id))

highlight = c("mature NK T cell","natural killer cell","T-helper 17 cell")

estimates = reshape2::dcast(ct_tests, experiment_id + term + feature_id ~ wfid, value.var = "estimate")
colnames(estimates)[4:5] = c("original","new")

stats.estimates = ddply(estimates, .(experiment_id, term), function(d) {
  ct <- cor.test(d$original, d$new, use = "paired.complete.obs", method = "pearson")
  data.frame(
    cor     = unname(ct$estimate),
    p_value = ct$p.value
  )
})
stats.estimates$label <- with(
  stats.estimates,
  sprintf("R = %.2f\np = %.2e", cor, p_value)
)

ggplot(estimates, aes(x = original, y = new)) +
  geom_smooth(method="lm", se=F, color = "grey") +
  geom_point(aes(color = ifelse(feature_id %in% highlight,"red","black"))) +
  geom_text(
    data = stats.estimates, size = 3,
    aes(label = label),
    x = -Inf, y = Inf,          # corner of each facet
    hjust = -0.1, vjust = 1.1,  # nudge inside the panel
    inherit.aes = FALSE
  ) +
  facet_grid(cols = vars(experiment_id), rows = vars(term), scales = "free") +
  scale_color_identity()+
  theme_light() +
  theme(strip.text = element_text(color = "black")) +
  labs(color = NULL)

p = reshape2::dcast(ct_tests, experiment_id + term + feature_id ~ wfid, value.var = "pvalue")
colnames(p)[4:5] = c("original","new")
stats.p = ddply(p, .(experiment_id, term), function(d) {
  ct <- cor.test(d$original, d$new, use = "paired.complete.obs", method = "pearson")
  data.frame(
    cor     = unname(ct$estimate),
    p_value = ct$p.value
  )
})
stats.p$label <- with(
  stats.p,
  sprintf("R = %.2f\np = %.2e", cor, p_value)
)

ggplot(p, aes(x = -log10(original), y = -log10(new))) +
  geom_smooth(method="lm", se=F, color = "grey") +
  geom_point(aes(color = ifelse(feature_id %in% highlight,"red","black"))) +
  geom_text(
    data = stats.p, size = 3,
    aes(label = label),
    x = -Inf, y = Inf,          # corner of each facet
    hjust = -0.1, vjust = 1.1,  # nudge inside the panel
    inherit.aes = FALSE
  ) +
  facet_grid(cols = vars(experiment_id), rows = vars(term), scales = "free") +
  scale_color_identity()+
  theme_light() +
  theme(strip.text = element_text(color = "black")) +
  labs(color = NULL)


feature_rank = reshape2::dcast(ct_tests, experiment_id + term + feature_id ~ wfid, value.var = "feature_rank")
colnames(feature_rank)[4:5] = c("original","new")
stats.feature_rank = ddply(feature_rank, .(experiment_id, term), function(d) {
  ct <- cor.test(d$original, d$new, use = "paired.complete.obs", method = "pearson")
  data.frame(
    cor     = unname(ct$estimate),
    p_value = ct$p.value
  )
})
stats.feature_rank$label <- with(
  stats.feature_rank,
  sprintf("R = %.2f\np = %.2e", cor, p_value)
)

ggplot(feature_rank, aes(x = original, y = new)) +
  geom_smooth(method="lm", se=F, color = "grey") +
  geom_point(aes(color = ifelse(feature_id %in% highlight,"red","black"))) +
  geom_text(
    data = stats.feature_rank, size = 3,
    aes(label = label),
    x = -Inf, y = Inf,          # corner of each facet
    hjust = -0.1, vjust = 1.1,  # nudge inside the panel
    inherit.aes = FALSE
  ) +
  facet_grid(cols = vars(experiment_id), rows = vars(term), scales = "free") +
  scale_color_identity()+
  theme_light() +
  theme(strip.text = element_text(color = "black")) +
  labs(color = NULL)



## Comparing cell contributions
## =======================================
unified_metadata = readRDS(get_workflow_outputs("wf-e82a5ab3b6"))
skin = read.csv("~/data/skin_modified.csv", header = F, sep = "\t")

extract_contributions = function(ccm) {
  lapply(ccm$datasets, function(d){
      res = t(exprs(d$cell_contribution$eset))
  }) %>% do.call(rbind,.)
}

old_contrib = extract_contributions(old_ccm) %>% data.frame(sample_id = rownames(.), check.names = F)
  old_contrib = merge(unified_metadata[,c("sample_id","dataset_id_short","sample_classification","condition")], old_contrib, by = "sample_id")
  old_contrib$version = "original"

new_contrib = extract_contributions(new_ccm) %>% data.frame(sample_id = rownames(.))
  colnames(new_contrib)[-22] = skin$V1[match(colnames(new_contrib)[-22], skin$V2)]
  new_contrib = merge(unified_metadata[,c("sample_id","dataset_id_short","sample_classification","condition")], new_contrib, by = "sample_id")
  new_contrib$version = "new"

contributions = rbind(old_contrib, new_contrib)
contrib_long = reshape2::melt(contributions, id.vars = c(1:4,26), variable.name = "cell", value.name = "contribution")
contrib_long = reshape2::dcast(contrib_long, sample_id + dataset_id_short + sample_classification + condition + cell ~ version)

stats = lapply(unique(contrib_long$cell), function(cell){
  lapply(unique(contrib_long$dataset_id_short), function(dataset){
    x = contrib_long[which(contrib_long$cell == cell & contrib_long$dataset_id_short == dataset),]
    ct = cor.test(x$original, x$new, use = "paired.complete.obs", method = "pearson")
    return(data.frame(dataset = dataset, cell = cell, cor = unname(ct$estimate), p = ct$p.value, row.names = NULL))
  }) %>% bind_rows()
}) %>% bind_rows()
stats$label <- with(
  stats,
  sprintf("R = %.2f\np = %.2e", cor, p)
)

contrib_long = merge(contrib_long, stats, by.x = c("dataset_id_short","cell"), by.y = c("dataset","cell"))
contrib_reduced = contrib_long[which(round(contrib_long$cor,2) < 0.95),]
contrib_reduced$cell = as.character(contrib_reduced$cell)

highlight = c("mature NK T cell" = "#bc5090","natural killer cell" = "#ff6361","T-helper 17 cell" = "#ffa600")

ggplot(contrib_reduced, aes(x = original, y = new)) +
  geom_smooth(method="lm", se=F, color = "grey") +
  geom_point(aes(color = ifelse(cell %in% names(highlight),cell,"black"))) +
  geom_text(check_overlap = T,
    data = contrib_reduced, size = 3,
    aes(label = label),
    x = -Inf, y = Inf,          # corner of each facet
    hjust = -0.1, vjust = 1.1,  # nudge inside the panel
    inherit.aes = FALSE
  ) +
  ggh4x::facet_grid2(cols = vars(cell), rows = vars(dataset_id_short), scales = "free", independent = "all") +
  scale_color_manual(values = highlight, na.value = "black")+
  theme_light() +
  theme(strip.text = element_text(color = "black")) +
  labs(color = NULL, title = "Correlation between cell contributions between old and new CCM versions", 
       subtitle = "filtered for pearson < 0.95")


##
