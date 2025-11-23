devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.cc.client)
library(tidyverse)
library(ComplexHeatmap)
library(patchwork)
library(bigrquery)

# Combine to one object to be run on CCM
# =======================================================
allSignatures = bq_table_download(x = bq_table(project = "cytoreason", dataset = "s05_atopic_dermatitis", table="X2Signatures"))
allSignatures = split(allSignatures$feature_id, allSignatures$signature)
pushToCC(allSignatures, tagsToPass = list(list(name="object",value="allSignatures")))
# wf-54cc53aaf8

nc = readRDS(get_workflow_outputs("wf-905e97fa64"))
  nc$smoothedRandom = lapply(nc$smoothedRandom, function(x) names(x))
  nc$random = unlist(nc$random, recursive = F)
  nc = unlist(nc,recursive = F)
  names(nc) = str_remove(names(nc),"random.|smoothedRandom.")

allSignatures = append(allSignatures, nc)
pushToCC(allSignatures, tagsToPass = list(list(name="object",value="allSignatures"),
                                          list(name="notes",value="including_randoms")))
# wf-2636758eaf


# Overlap graphs
# ===============================
ann = function(signatures){
  df = data.frame(row.names = names(signatures),
                  # agonist = sapply(names(signatures), function(x) str_split(x," ")[[1]][1]),
                  treatment = str_extract(names(signatures),"activation|inhibition"),
                  refined = str_extract(names(signatures),"refined"))

  ann_colors = list(
    # agonist = c("x2" = "#003f5c", "Untreated" = "#374c80", "CST14" = "#7a5195","Icatibant" = "#bc5090", "PAMP12" = "#ef5675", "SP" = "#ff764a", "aIgE" = "#ffa600"),  # replace with actual agonist names
    treatment = c("activation" = "orange", "inhibition" = "purple"),
    refined = c("refined" = "black", "NA" = "white")  # handle NA if needed
  )

  ann = HeatmapAnnotation(df = df, which = "row", col = ann_colors)
}

reconstruct_signature_name = function(sig){
  signature = str_split(sig, "_")[[1]]
  new_sig = c(signature[1],signature[3],signature[2])
  if("refined" %in% signature) {
    new_sig = c(new_sig, "(refined)")
  }
  return(paste0(new_sig, collapse = " "))
}

allSignatures = readRDS(get_workflow_outputs("wf-6d9a2c3a80"))
signatures = allSignatures[str_detect(names(allSignatures),"50")]
signatures = signatures[!str_detect(names(signatures), "string|archs(?!_)|p_refined|top100|bottom100|50_refined")]
signatures = signatures[!str_detect(names(signatures), "IgE_inhibition|CST14_inhibition|Icatibant_inhibition|PAMP12_inhibition|SP_inhibition_early|Untreated_inhibition")]
signatures = signatures[!str_detect(names(signatures), "Tryptase")]
signatures = signatures[!str_detect(names(signatures), "x2_inhibition_late")]

hm.signatures = signatures
names(hm.signatures) = lapply(names(hm.signatures), reconstruct_signature_name)

activation = hm.signatures[str_detect(names(hm.signatures),"activation")]
inhibition = hm.signatures[str_detect(names(hm.signatures),"inhibition")]

overlap = sapply(hm.signatures, function(x) sapply(hm.signatures, function(y) length(intersect(x,y))))
overlap_activation = sapply(activation, function(x) sapply(activation, function(y) length(intersect(x,y))))
overlap_inhibition = sapply(inhibition, function(x) sapply(inhibition, function(y) length(intersect(x,y))))

# all signatures
png("~/analysis-s05/figures/X2_Signature/overlap.png", res = "100", bg = "transparent", width = 900, height = 700)
Heatmap(overlap, name = "overlap", row_names_gp = gpar(fontsize = 10), right_annotation = ann(signatures),
        column_names_gp = gpar(fontsize=10), cell_fun = function(j, i, x, y, width, height, fill) {
          if(overlap[i, j] >= 25) grid.text(overlap[i, j], x, y, gp = gpar(fontsize = 10))})
dev.off()


# activation
png("~/analysis-s05/figures/X2_Signature/overlap_activation.png", res = "100", bg = "transparent", width = 900, height = 700)
Heatmap(overlap_activation, name = "overlap", row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize=10), cell_fun = function(j, i, x, y, width, height, fill) {
          if(overlap_activation[i, j] >= 25) grid.text(overlap_activation[i, j], x, y, gp = gpar(fontsize = 10))})
dev.off()


# inhibition
png("~/analysis-s05/figures/X2_Signature/overlap_inhibition.png", res = "100", bg = "transparent", width = 900, height = 700)
Heatmap(overlap_inhibition, name = "overlap", row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize=10), cell_fun = function(j, i, x, y, width, height, fill) {
          if(overlap_inhibition[i, j] >= 25) grid.text(overlap_inhibition[i, j], x, y, gp = gpar(fontsize = 10))})
dev.off()



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

mcs = signatures[!str_detect(names(signatures),"refined")]
archs_refined = signatures[str_detect(names(signatures),"refined")]

results <- compare_gene_loss(mcs, archs_refined, "_archs_refined")
results$sig = ifelse(!str_detect(results$signature,"archs"),"Mast Cell\nSpecific","Refined")
results$signature = sapply(results$signature, reconstruct_signature_name)

p1 = ggplot(results, aes(x = lost, y = signature, fill = sig)) +
  geom_col(position = "dodge", fill = "#0087c5") +
  geom_text(aes(label = lost), position = position_dodge(0.75), color = "black", hjust = 1.5, size = 3) +
  scale_x_continuous(expand = c(0.01,0), limits = c(-50,0))+
  scale_fill_identity()+
  theme_minimal() +
  labs(title = "Number of genes lost in the refinement", x = "Number of genes lost", y = "Signature", fill = "signature\norigin")

p2 = ggplot(results, aes(x = percent_lost, y = signature, fill = sig)) +
  geom_col(position = "dodge", fill = "#8abfe7") +
  geom_text(aes(label = paste0(percent_lost, "%"), group = sig), size = 3, position = position_dodge(0.85), hjust = -0.1, check_overlap = T, color = "black") +
  scale_x_continuous(expand = c(0,0), limits = c(0,100))+
  scale_fill_identity()+
  theme_minimal() +
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(x = "% Lost Genes", title = " ", fill="lost from\nwhich\nsignature")

p1 + p2 + plot_layout(guides = "collect", widths = c(8,1.7)) & theme(legend.position = "right")
ggsave("~/analysis-s05/figures/X2_Signature/genes_lost_in_refinement.png", width = 1000, height = 500, units = "px", scale = 2.7, bg = "white", )
