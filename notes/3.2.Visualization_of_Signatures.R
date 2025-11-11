library(cytoreason.cc.client)
library(tidyverse)
library(reshape2)

allSignatures = readRDS(get_workflow_outputs("wf-30b7952de4"))
geneMapping = toTable(org.Hs.eg.db::org.Hs.egSYMBOL)


### Mast cell-specific
### =========================
gx = readRDS(get_workflow_outputs("wf-edcade5c80"))
gx = gx[-which(gx$term == "unactivated"),]
gx$agonist[str_detect(gx$agonist, "All")] <- "Pooled"
gx$agonist[str_detect(gx$agonist, "SP")] <- "SP"
gx$symbol = geneMapping$symbol[match(gx$feature_id, geneMapping$gene_id)]

sig_4hr = "Pooled"
sig_24hr = "SP"

sigs = allSignatures[c("x2_activation_early_50_p","SP_activation_late_50_p", "x2_inhibition_early_50_p", "SP_inhibition_late_50_p")] %>%
  reshape2::melt(value.name = "gene") %>%
  rename(signature = L1) %>%
  mutate(sig = case_when(str_detect(signature, "early") ~ paste0("4hr - ",sig_4hr),
                         str_detect(signature, "late") ~ paste0("24hr - ", sig_24hr))) %>%
  mutate(dir = case_when(str_detect(signature, "activation") ~ "activation",
                         str_detect(signature, "inhibition") ~ "inhibition")) %>%
  mutate(comparison = case_when(str_detect(signature, "activation") & str_detect(signature, "early") ~ "activated_vs_unactivated_4hr",
                                str_detect(signature, "inhibition") & str_detect(signature, "early") ~ "inhibited_vs_uninhibited_4hr",
                                str_detect(signature, "activation") & str_detect(signature, "late") ~ "SP_activated_vs_unactivated_24hr",
                                str_detect(signature, "inhibition") & str_detect(signature, "late") ~ "SP_inhibited_vs_uninhibited_24hr"))
  

chosenSignatures = rbind(gx[which(gx$dataset_id == "X2_4hr__GPL24676" & gx$agonist == sig_4hr),], 
                         gx[which(gx$dataset_id == "X2_24hr__GPL34281" & gx$agonist == sig_24hr),])
chosenSignatures = chosenSignatures %>%
  mutate(time = case_when(agonist == sig_4hr ~ paste0("4hr - ",sig_4hr), 
                          agonist == sig_24hr ~ paste0("24hr - ", sig_24hr)))

chosenSignatures = merge(chosenSignatures, sigs, all = T, by.x = c("comparison", "feature_id"), by.y = c("comparison", "gene"))

ggplot(chosenSignatures, aes(x = estimate, y = log10_pvalue)) +
  geom_point(aes(color = dir)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text(x = 4, y = 1.1, label = "pvalue = 0.05", check_overlap = T, size = 3, aes(hjust= ifelse(time == "4hr - Pooled", "right", -2))) +
  scale_color_manual(values = c("activation" = "#00d484", "inhibition" = "#e20069")) +
  scale_x_continuous(name = "estimate (log2 Fold Change)")+
  scale_y_continuous(name = "-log10(p-value)", sec.axis = dup_axis(name = "Comparison"))+
  facet_grid(cols = vars(time), rows = vars(term), scales = "free") +
  theme_bw() +
  theme(axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(), axis.text.y.right = element_blank()) +
  labs(title = "Top/Bottom 50 genes serve as X2 gene signatures", color = NULL)
