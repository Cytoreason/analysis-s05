devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.cc.client)
library(tidyverse)
library(ComplexHeatmap)
library(patchwork)

msc = readRDS(get_workflow_outputs("wf-ba31d91a68"))$genes_only # mast-cell specific
archs = readRDS(get_workflow_outputs("wf-f5930577e9"))$genes # msc propagated in archs
insilico_archs = readRDS(get_workflow_outputs("wf-6028a1f83b"))$genes # in silico propagated in archs
archs_refined = readRDS(get_workflow_outputs("wf-3b31fba621"))
insilico_archs_refined = readRDS(get_workflow_outputs("wf-7305f8adfa"))
nc = readRDS(get_workflow_outputs("wf-905e97fa64"))
  nc$smoothedRandom = lapply(nc$smoothedRandom, function(x) names(x))
  nc$random = unlist(nc$random, recursive = F)
  nc = unlist(nc,recursive = F)
  names(nc) = str_remove(names(nc),"random.|smoothedRandom.")

# Combine to one object to be run on CCM
# =======================================================
allSignatures = do.call(c, list(msc, archs, insilico_archs, archs_refined, insilico_archs_refined, nc))
pushToCC(allSignatures, tagsToPass = list(list(name="object",value="allSignatures"),
                                          list(name="notes",value="including_agonists")))
# wf-791af609e8

allSignatures = allSignatures[!str_detect(names(allSignatures), "SP_|CST14|Icatibant|IgE|Untreated|PAMP12")]
pushToCC(allSignatures, tagsToPass = list(list(name="object",value="allSignatures"),
                                          list(name="notes",value="without_agonists")))
# wf-5a525e011d



# Overlap graphs
# -----------------------------
# 50/100
allSignatures = do.call(c, list(msc, archs, insilico_archs, archs_refined, insilico_archs_refined))
allSignatures_50 = allSignatures[str_detect(names(allSignatures),"50")]
allSignatures_100 = allSignatures[str_detect(names(allSignatures),"100")]

overlap_50 = sapply(allSignatures_50, function(x) sapply(allSignatures_50, function(y) length(intersect(x,y))))
overlap_100 = sapply(allSignatures_100, function(x) sapply(allSignatures_100, function(y) length(intersect(x,y))))

png("~/analysis-s05/figures/X2_Signature/overlap_50.png", res = "100", bg = "transparent", width = 1200, height = 900)
Heatmap(overlap_50, name = "overlap", row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize=8), cell_fun = function(j, i, x, y, width, height, fill) {
        if(overlap_50[i, j] >= 25) grid.text(overlap_50[i, j], x, y, gp = gpar(fontsize = 8))})
dev.off()

png("~/analysis-s05/figures/X2_Signature/overlap_100.png", res = "100", bg = "transparent", width = 1200, height = 900)
Heatmap(overlap_100, name = "overlap", row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize=8), cell_fun = function(j, i, x, y, width, height, fill) {
          if(overlap_100[i, j] >= 25) grid.text(overlap_100[i, j], x, y, gp = gpar(fontsize = 8))})
dev.off()

# activation/inhibition
activation = allSignatures[str_detect(names(allSignatures),"activation")]
inhibition = allSignatures[str_detect(names(allSignatures),"inhibition")]

overlap_activation = sapply(activation, function(x) sapply(activation, function(y) length(intersect(x,y))))
overlap_inhibition = sapply(inhibition, function(x) sapply(inhibition, function(y) length(intersect(x,y))))

png("~/analysis-s05/figures/X2_Signature/overlap_activation.png", res = "100", bg = "transparent", width = 1200, height = 900)
Heatmap(overlap_activation, name = "overlap", row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize=8), cell_fun = function(j, i, x, y, width, height, fill) {
          if(overlap_activation[i, j] >= 25) grid.text(overlap_activation[i, j], x, y, gp = gpar(fontsize = 8))})
dev.off()

png("~/analysis-s05/figures/X2_Signature/overlap_inhibition.png", res = "100", bg = "transparent", width = 1200, height = 900)
Heatmap(overlap_inhibition, name = "overlap", row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize=8), cell_fun = function(j, i, x, y, width, height, fill) {
          if(overlap_inhibition[i, j] >= 25) grid.text(overlap_inhibition[i, j], x, y, gp = gpar(fontsize = 8))})
dev.off()



# Which genes joined from the network
# --------------------------------------------
results <- data.frame(signature = character(),
                      existed = integer(),
                      new = integer(),
                      percent_new = numeric(),
                      stringsAsFactors = FALSE)

for (sig in names(msc)) {
  original_genes <- msc[[sig]]
  propagated_genes <- archs[[paste0(sig, "_archs")]]
  
  existed <- sum(original_genes %in% propagated_genes)
  new <- sum(!(propagated_genes %in% original_genes))
  percent_new <- round((new / existed) * 100, 0)
  
  results <- rbind(results, data.frame(
    signature = sig, existed = existed, new = new, percent_new = percent_new, stringsAsFactors = FALSE))
}

results = reshape2::melt(results, id.vars = c(1,4), variable.name = "origin", value.name = "n")
results$origin = factor(results$origin, ordered = T, levels = c("new","existed"))

p1 = ggplot(results, aes(x = n, y = signature, fill = origin)) +
  geom_col(position = "stack") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), color = "black") +
  scale_fill_manual(values = c("#FFD265","#0A7B83"))+
  scale_x_continuous(expand = c(0,0))+
  theme_minimal() +
  labs(title = "Number of genes added from the ARCHS network", x = "Number of genes", y = "Signature")

p2 = ggplot(results[which(results$origin == "existed"),], aes(x = percent_new, y = signature)) +
  geom_col(fill = "grey") +
  geom_text(aes(label = paste0(percent_new, "%")), hjust = -0.1, check_overlap = T, color = "black") +
  scale_x_continuous(expand = c(0,0), limits = c(0,100))+
  theme_minimal() +
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(x = "% New Genes", title = " ")

p1 + p2 + plot_layout(guides = "collect", widths = c(8,1)) & theme(legend.position = "right")
ggsave("~/analysis-s05/figures/X2_Signature/gene_origin_archs.png", width = 1500, height = 700, units = "px", scale = 2.8, bg = "white")


# Which genes were lost in the refinement
# ------------------------------------------------
compare_gene_loss <- function(original_list, refined_list, suffix) {
  results <- data.frame(signature = character(), origin = integer(), lost = integer(), stringsAsFactors = FALSE  )
  for (sig in names(original_list)) {
    original_genes <- original_list[[sig]]
    refined_sig <- paste0(sig, suffix)
    refined_genes <- refined_list[[refined_sig]]
    
    origin_count <- length(original_genes)
    lost_count <- (-1)*sum(!(original_genes %in% refined_genes))
    percent_lost <- (-1)*round((lost_count / origin_count) * 100, 1)
    
    
    results <- rbind(results, data.frame(
      signature = sig,
      origin = origin_count,
      lost = lost_count,
      percent_lost = percent_lost,
      stringsAsFactors = FALSE
    ))
  }
  return(results)
}


results <- rbind(compare_gene_loss(msc, archs_refined, "_archs_refined"),
                 compare_gene_loss(archs, archs_refined, "_refined"))
results$sig = ifelse(!str_detect(results$signature,"archs"),"MCS","ARCHS")
results$signature = str_remove(results$signature, "_archs")

p1 = ggplot(results, aes(x = lost, y = signature, fill = sig)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = lost), position = position_dodge(0.75), color = "black", hjust = 1.5, size = 3) +
  scale_x_continuous(expand = c(0,0), limits = c(-50,0))+
  scale_fill_manual(values = c("#FFD265","#0A7B83"))+
  theme_minimal() +
  labs(title = "Number of genes lost in the refinement", x = "Number of genes lost", y = "Signature", fill = "signature\norigin")

p2 = ggplot(results, aes(x = percent_lost, y = signature, fill = sig)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = paste0(percent_lost, "%"), group = sig), size = 3, position = position_dodge(0.85), hjust = -0.1, check_overlap = T, color = "black") +
  scale_x_continuous(expand = c(0,0), limits = c(0,100))+
  scale_fill_manual(values = c("#f7deb1","#a5bfc1"))+
  theme_minimal() +
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(x = "% Lost Genes", title = " ", fill="lost from\nwhich\nsignature")

p1 + p2 + plot_layout(guides = "collect", widths = c(8,1)) & theme(legend.position = "right")
ggsave("~/analysis-s05/figures/X2_Signature/genes_lost_in_refinement.png", width = 1500, height = 1200, units = "px", scale = 2.8, bg = "white", )

