library(cytoreason.cc.client)
library(tidyverse)
library(reshape2)

allSignatures = readRDS(get_workflow_outputs("wf-30b7952de4"))
geneMapping = toTable(org.Hs.eg.db::org.Hs.egSYMBOL)


### Mast cell-specific volcan
### =======================================
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
# wf-420b1704bb

highlightedGenes = data.frame(signature = character(), symbol = character(), highlight = logical())
highlightedGenes = rbind(highlightedGenes, data.frame(signature = "x2_activation_early_50_p", 
                                                      symbol = c("TLR2","TACR1","RAB7B","S100P",""), 
                                                      highlight = T))
  
ggplot(chosenSignatures, aes(x = estimate, y = log10_pvalue)) +
  geom_point(aes(color = dir)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text(x = 4, y = 1.1, label = "pvalue = 0.05", check_overlap = T, size = 3, hjust = "left") +
  scale_color_manual(values = c("activation" = "#00d484", "inhibition" = "#e20069"), breaks = ~ .x[!is.na(.x)]) +
  scale_x_continuous(name = "estimate (log2 Fold Change)")+
  scale_y_continuous(name = "-log10(p-value)", sec.axis = dup_axis(name = "Comparison"))+
  facet_grid(cols = vars(time), rows = vars(term), scales = "free_y") +
  theme_bw() +
  theme(axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(), axis.text.y.right = element_blank()) +
  labs(title = "Top/Bottom 50 genes serve as X2 gene signatures", color = NULL)
ggsave("~/analysis-s05/figures/X2_Signature/mcs_signatures.png", width = 800, height = 500, units = "px", scale = 3.4)


### Which genes joined from the network
### ==============================================
sigs = bq_table_download(x = bq_table(project = "cytoreason", dataset = "s05_atopic_dermatitis", table="X2Signatures"))
sigs = split(sigs$feature_id, sigs$signature)
sigs = sigs[str_detect(names(sigs),"_p|MRGPRX2")]

sigs_50 = sigs[str_detect(names(sigs),"50")]
mcs = sigs_50[!str_detect(names(sigs_50),"archs|string|refined")]
archs = sigs_50[str_detect(names(sigs_50),"archs")] %>%
  .[!str_detect(names(.),"refined")]
archs_refined = sigs_50[str_detect(names(sigs_50),"archs")] %>%
  .[str_detect(names(.),"refined")]
string = sigs_50[str_detect(names(sigs_50),"string")] %>%
  .[!str_detect(names(.),"refined")]
string_refined = sigs_50[str_detect(names(sigs_50),"string")] %>%
  .[str_detect(names(.),"refined")]

highlight = c("x2_activation_early_50_p", "SP_activation_late_50_p","x2_inhibition_early_50_p", "SP_inhibition_late_50_p")

newGenes = data.frame(signature = character(),
                      new = character(),
                      network = character())
rmGenes = data.frame(signature = character(),
                      removed = character(),
                      network = character())
results <- data.frame(signature = character(),
                      original = integer(),
                      net = integer(),
                      network = character(),
                      percent_new = numeric(),
                      stringsAsFactors = FALSE)

for (sig in names(mcs)) {
  original_genes <- mcs[[sig]]
  propagated_genes <- archs[[paste0(sig, "_archs")]]
  
  existed <- sum(original_genes %in% propagated_genes)
  new <- sum(!(propagated_genes %in% original_genes))
  percent_new <- round((new / length(original_genes)) * 100, 0)
  
  results <- rbind(results, data.frame(
    signature = sig, original = existed, net = new, percent_new = percent_new, network = "ARCHS Network", stringsAsFactors = FALSE))
  
  newGenes <- rbind(newGenes, data.frame(
    signature = sig, new = propagated_genes[!(propagated_genes %in% original_genes)], network = "ARCHS Network", stringsAsFactors = FALSE))
  rmGenes <- rbind(rmGenes, data.frame(
    signature = sig, removed = original_genes[!(original_genes %in% propagated_genes)], network = "ARCHS Network", stringsAsFactors = FALSE))
  
  propagated_genes <- string[[paste0(sig, "_string")]]
  
  existed <- sum(original_genes %in% propagated_genes)
  new <- sum(!(propagated_genes %in% original_genes))
  percent_new <- round((new / length(original_genes)) * 100, 0)
  
  results <- rbind(results, data.frame(
    signature = sig, original = existed, net = new, percent_new = percent_new, network = "STRINGdb Network", stringsAsFactors = FALSE))
  
  newGenes <- rbind(newGenes, data.frame(
    signature = sig, new = propagated_genes[!(propagated_genes %in% original_genes)], network = "STRINGdb Network", stringsAsFactors = FALSE))
  
  rmGenes <- rbind(rmGenes, data.frame(
    signature = sig, removed = original_genes[!(original_genes %in% propagated_genes)], network = "STRINGdb Network", stringsAsFactors = FALSE))
  
  if(!is.null(highlight)) {
    results = results[which(results$signature %in% highlight),]
    newGenes = newGenes[which(newGenes$signature %in% highlight),]
    rmGenes = rmGenes[which(rmGenes$signature %in% highlight),]
  }
}

results = reshape2::melt(results, id.vars = c(1,5,4), variable.name = "origin", value.name = "n")
results$origin = factor(results$origin, ordered = T, levels = c("net","original"))
# results$ep = ifelse(str_detect(results$signature,"_ep"),"estimate x logP","estimate")
# results$sigsize = ifelse(str_detect(results$signature,"_50"),"n=50","n=100")
# results$signature = str_remove(results$signature, "_50|_100")

p1 = ggplot(results, aes(x = n, y = signature, fill = origin)) +
  geom_col(position = "stack") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), color = "black") +
  scale_fill_manual(values = c("#FFD265","#0A7B83"))+
  scale_x_continuous(expand = c(0,0))+
  theme_minimal() +
  # facet_grid(rows = vars(ep), cols = vars(sigsize), scales="free") +
  facet_wrap(~network, scales="free", nrow = 2) +
  labs(title = "Number of genes added from the network", x = "Number of genes", y = "Signature") +
  theme(strip.text = element_text(face = "bold"))

p2 = ggplot(results[which(results$origin == "original"),], aes(x = percent_new, y = signature)) +
  geom_col(fill = "grey") +
  geom_text(aes(label = paste0(percent_new, "%")), hjust = -0.1, check_overlap = T, color = "black") +
  scale_x_continuous(expand = c(0,0), limits = c(0,100))+
  theme_minimal() +
  # facet_grid(rows = vars(ep), cols = vars(sigsize), scales="free") +
  facet_wrap(~ network, scales="free", nrow = 2, labeller = as_labeller(c("ARCHS Network" = "", "STRINGdb Network" = ""))) +
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(x = "% New Genes", title = " ")

p1 + p2 + plot_layout(guides = "collect", widths = c(8,2)) & theme(legend.position = "right")
ggsave("~/analysis-s05/figures/X2_Signature/gene_origin_chosen.png", width = 1000, height = 500, units = "px", scale = 3, bg = "white")


chosenSignatures = readRDS(get_workflow_outputs("wf-420b1704bb")) %>%
  mutate(feature_id = as.integer(feature_id))

newGenes = newGenes %>%
  mutate(comparison = case_when(str_detect(signature, "activation") & str_detect(signature, "early") ~ "activated_vs_unactivated_4hr",
                                str_detect(signature, "inhibition") & str_detect(signature, "early") ~ "inhibited_vs_uninhibited_4hr",
                                str_detect(signature, "activation") & str_detect(signature, "late") ~ "SP_activated_vs_unactivated_24hr",
                                str_detect(signature, "inhibition") & str_detect(signature, "late") ~ "SP_inhibited_vs_uninhibited_24hr"))
rmGenes = rmGenes %>%
  mutate(comparison = case_when(str_detect(signature, "activation") & str_detect(signature, "early") ~ "activated_vs_unactivated_4hr",
                                str_detect(signature, "inhibition") & str_detect(signature, "early") ~ "inhibited_vs_uninhibited_4hr",
                                str_detect(signature, "activation") & str_detect(signature, "late") ~ "SP_activated_vs_unactivated_24hr",
                                str_detect(signature, "inhibition") & str_detect(signature, "late") ~ "SP_inhibited_vs_uninhibited_24hr"))




addedByNetwork <- chosenSignatures %>%
  # Combine New and Removed info into one table
  left_join(
    bind_rows(
      newGenes %>% mutate(gene_id = new, status = "New"),
      rmGenes %>% mutate(gene_id = removed, status = "Removed")
    ),
    by = c("feature_id" = "gene_id", "comparison")
  ) %>%
  mutate(network = str_remove(network, " Network")) %>%
  
  # Pivot to get separate columns for ARCHS and STRINGdb
  pivot_wider(names_from = network, values_from = status) %>%
  
  # Create combined label
  mutate(newPropagated = case_when(
    ARCHS == "New" & STRINGdb == "New" ~ "New in both",
    ARCHS == "New" ~ "New in ARCHS",
    STRINGdb == "New" ~ "New in STRING",
    ARCHS == "Removed" ~ "Removed by ARCHS",
    STRINGdb == "Removed" ~ "Removed by STRING",
    TRUE ~ NA_character_
  ))


ggplot(addedByNetwork, aes(x = estimate, y = log10_pvalue)) +
  geom_point(color = "grey") +
  geom_point(data = addedByNetwork[!is.na(addedByNetwork$newPropagated),], aes(color = newPropagated))+
  ggrepel::geom_text_repel(data = addedByNetwork[!is.na(addedByNetwork$newPropagated),], show.legend = F, 
                           aes(label = symbol, color = newPropagated), size = 3, max.overlaps = 20)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text(x = 4, y = 1.1, label = "pvalue = 0.05", check_overlap = T, size = 3, hjust = "right") +
  # scale_color_manual(values = c("New" = "#00d484", "Removed" = "#e20069"), breaks = ~ .x[!is.na(.x)]) +
  scale_x_continuous(name = "estimate (log2 Fold Change)")+
  scale_y_continuous(name = "-log10(p-value)", sec.axis = dup_axis(name = "Comparison"))+
  facet_grid(cols = vars(time), rows = vars(term), scales = "free") +
  theme_bw() +
  theme(axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(), axis.text.y.right = element_blank()) +
  labs(color = NULL)
ggsave("~/analysis-s05/figures/X2_Signature/mcs_signatures_propagated.png", width = 800, height = 500, units = "px", scale = 3.4)

write.csv(addedByNetwork[!is.na(addedByNetwork$newPropagated),], "~/data/genesAddedByNetwork.csv")
