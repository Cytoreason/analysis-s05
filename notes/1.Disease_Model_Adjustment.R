# Due to the disease model being over-corrected with both 'meta pc1' and 'keratinocytes',
# we used the model teams' disease model - and tested the changes in cell contributions following the adjustment
devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(cytoreason.deconvolution)
library(cytoreason.cc.client)
library(cytoreason.assets)
library(ComplexHeatmap)
library(tidyverse)

ccm = as_ccm_fit("wf-08a6a0a503") # original disease model
ccm_fit_AD <- read_data(AssetData("wf-08a6a0a503:0:ccm_fit_with_data.qs"))
selectedEffect = "L_vs_HC"

## 1. Extract cell contributions and ct_test per analysis model
## =================================================================
# Filter ccm_fit to only relevant datasets/model: (that have effect_id == "AD" or "L_vs_NL")
model_ids_AD <- list_model_ids_for_meta_analysis(ccm_fit_AD)
model_ids_AD <- model_ids_AD[model_ids_AD == selectedEffect]

ccm_fit_AD_filt <- ccm_fit_AD
for (d in names(ccm_fit_AD_filt[["datasets"]])){
  for (m in names(ccm_fit_AD_filt[["datasets"]][[d]][["model"]])){
    if((m %in% names(model_ids_AD))==FALSE){
      ccm_fit_AD_filt[["datasets"]][[d]][["model"]][[m]] <- NULL
    }
  }
  if (length(ccm_fit_AD_filt[["datasets"]][[d]][["model"]])==0){
    ccm_fit_AD_filt[["datasets"]][[d]] <- NULL
  }
}


# Extraction, RES_AD: one element per analysis model
RES_AD <- apply_model_analysis(ccm_fit_AD_filt, function(model, ccm_dataset){
  eset <- assayDataExpression(ccm_dataset)
  signature_collection <- cytoreason.deconvolution::get_signature_collection_name(ccm_dataset)
  adjusted_elts <- cytoreason.ccm.pipeline:::.list_adjusted_assay_data_elements(ccm_dataset)

  wf = mapply_dist(FUN = function(elt, eset, signature_collection, ccm_dataset){

    # compute contributions on adjusted expression element
    eset_elt <- cytoreason.ccm.pipeline:::subset_assayData(eset, elt)
    ct_contrib_adjusted <- cytoreason.ccm.pipeline:::call_with_cache(cytoreason.ccm.pipeline::ccm_service_cell_contribution,
                                                                     expression_data = eset_elt,
                                                                     signature_collection = signature_collection)

    # compute ct_test on adjusted contributions
    analysisResultElement(ccm_dataset, "cell_contribution", overwrite = TRUE) <-  ct_contrib_adjusted
    fit_ct_test <- cytoreason.ccm.pipeline:::ccm_service_cell_contribution_differences(ccm_dataset = ccm_dataset,
                                                                                       design = ccm_dataset[["model"]],
                                                                                       match.pair = "pairwise")
    return(list(cell_contribution = ct_contrib_adjusted, ct_test = fit_ct_test))
  },
  elt = setNames(nm = adjusted_elts),
  MoreArgs = list(eset=eset, signature_collection=signature_collection, ccm_dataset=ccm_dataset),
  image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_0.81.9",
  tags = list(list(name = "task",value="extract_adjusted_cell_contributions"),
              list(name="ccm_wfid",value="wf-08a6a0a503"),
              list(name="model",value=names(ccm_dataset$model))))

  get_workflow(wf, wait = T)

  res = read_asset(wf)
  names(res) = adjusted_elts

  cytoreason.ccm.pipeline:::.set_ccm_service_class(list(submodel = res), "cell_contribution_adjusted")

}, by_model = TRUE)
pushToCC(RES_AD) # wf-a958262aa8 L_vs_NL
pushToCC(RES_AD) # wf-884b5ef8ca DZ_vs_HC
pushToCC(RES_AD) # wf-9a5fce551b L_vs_HC


# Meta-analysis ct_test (META_RES_AD)
# ----------------------------------------
submodel_set_AD <- unique(unlist(lapply(RES_AD, function(x){names(getObjectElement(x, "submodel"))})))
META_RES_AD <- lapply(submodel_set_AD, function(submodel){

  ct_test_res_list <- lapply(RES_AD, function(x){
    getObjectElementPath(x, file.path("submodel", submodel, "ct_test"), null.ok = TRUE)

  })
  ct_test_res_list <- ct_test_res_list[lengths(ct_test_res_list) > 0L]
  cytoreason.ccm.pipeline:::meta_analysis_service_ct_test(ct_test_res_list, .cc_image = NA)

})
pushToCC(META_RES_AD) # wf-8fea936a47 L_vs_NL
pushToCC(META_RES_AD) # wf-cf33062383 DZ_vs_HC
pushToCC(META_RES_AD) # wf-b5ad88acb3 L_vs_HC

# extract meta-analysis results:
cellMapping = read.csv(get_workflow_outputs("wf-399caff221"))
cellMapping$cell_type_name[which(cellMapping$cell_type_id == "CRCL_0000028")] <- "CD8-positive, alpha-beta T cell"
cellMapping$cell_type_name[which(cellMapping$cell_type_id == "CRCL_0000052")] <- "effector memory CD4-positive, alpha-beta T cell"
cellMapping$cell_type_name[which(cellMapping$cell_type_id == "CRCL_0000251")] <- "mature NK T cell"

adjusted <- lapply(1:length(META_RES_AD), function(x){
  res_df <- statistic_table(META_RES_AD[[x]])
  res_df <- subset(res_df, term==selectedEffect)
  res_df$submodel <- submodel_set_AD[x]
  res_df
})
names(adjusted) <- submodel_set_AD
adjusted = do.call(rbind, adjusted)
adjusted$feature_id = cellMapping$cell_type_name[match(adjusted$feature_id, cellMapping$cell_type_id)]

# Add bulk ct-test:
bulk <- statistic_table(ccm_fit_AD$meta$AD$ct_test)
bulk$submodel <- "bulk"

# combine bulk and adjusted
ct_test <- rbind(bulk,adjusted)
rownames(ct_test) = NULL

ct_test$log10fdr = -log10(ct_test$fdr)
ct_test$star = sapply(ct_test$fdr, starryNight)
pushToCC(ct_test)
# wf-f50bcb1e12 DZ_vs_HC
# wf-d28da3033f L_vs_NL
# wf-d22af3122f L_vs_HC
# a table comparing the ct_test estimates and stars of the keratinocyte cell from the three options: wf-4fc19f62fa, wf-2468fdb147


## Visualizations
## =================================================
# 1. Comparing ct test results for keratinocytes only
# ------------------------------------------------------
mapAdjName = cytoreason.assets::read_asset("wf-1005306ab9")
ct_heatmap = read_asset("wf-4fc19f62fa") %>% column_to_rownames(var = "submodel") %>% .[,-1] %>% as.matrix
ct_stars = read_asset("wf-2468fdb147") %>% column_to_rownames(var = "submodel") %>% .[,-c(1:2)] %>% as.matrix
rownames(ct_heatmap) = mapAdjName$adjustment_variables[match(str_remove(rownames(ct_heatmap),"exprs_"),mapAdjName$submodel)]

png(paste0("~/analysis-s05/figures/adjustment/ct_test_keratinocytesOnly.png"), res = "100", bg = "transparent", width = 700, height = 400)
Heatmap(ct_heatmap[,1,drop=F], name = "cell\ncontribution",
        col = circlize::colorRamp2(c(0,2),c("#c1e7ff","#004c6d")), cluster_rows = F, 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(round(ct_heatmap[i, j],2), x, y, gp = gpar(fontsize = 10))},
        column_labels = c("keratinocyte\ncell\ncontribution")) +
  Heatmap(ct_heatmap[,-1], name = "estimate", cluster_rows = F, cluster_columns = F, 
          col = circlize::colorRamp2(c(-1,0,1.5),c("blue","white","red")),
          column_labels = str_remove(colnames(ct_heatmap)[-1],"estimate_"),
          column_title = "ct_test estimate of keratinocytes", column_title_side = "top",
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(ct_stars[i, j], x, y, gp = gpar(fontsize = 10))})
dev.off()


# 2. Meta ct_test
# ----------------------
loadings = ccm_fit_AD$multi$meta_pca$`1`$v # add cell loadings to heatmap
colnames(loadings) = paste0("PC",1:8)
loading_ann = HeatmapAnnotation(df = as.data.frame(loadings[,1,drop=F]), name = "PC1", show_annotation_name = T, 
                                which="row", col = list(PC1=circlize::colorRamp2(c(-0.3,0,0.3),c("#488f31","#f1f1f1","#de425b"))))

mat = reshape2::acast(ct_test, feature_id ~ submodel, value.var = "estimate")
mat.star = reshape2::acast(ct_test, feature_id ~ submodel, value.var = "star")
mat = mat[rownames(loadings),]
mat.star = mat.star[rownames(loadings),]

png(paste0("~/analysis-s05/figures/adjustment/ct_test_",selectedEffect,".png"), res = "100", bg = "transparent", width = 900, height = 900)
Heatmap(mat, name = "estimate", column_title = paste0("Change in ct_test estimate over adjustments (",selectedEffect,")"), 
        right_annotation = loading_ann, cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(mat.star[i, j], x, y, gp = gpar(fontsize = 10))
        })
dev.off()


# 3. Median cell contribution
# -----------------------------------
bulks = list()
contributions <- lapply(submodel_set_AD, function(submodel){
  contrib_list <- lapply(RES_AD, function(x){
    if(is.null(x$submodel[[submodel]])) { return(NULL) }
    dataset = paste0(strsplit(x[["submodel"]][[submodel]][["ct_test"]][["model"]][["name"]],"__")[[1]][2:3],collapse = "__")
    adj = x$submodel[[submodel]][["cell_contribution"]]$object$eset@assayData$exprs
    bulks[[dataset]] <<- ccm_fit_AD$datasets[[dataset]]$cell_contribution$eset@assayData$exprs
    return(adj)
  })
  contrib_list <- do.call(cbind, contrib_list)
  contrib_df = reshape2::melt(contrib_list)
  colnames(contrib_df) = c("feature_id","sample","cell_contribution")
  contrib_df$submodel = submodel
  return(contrib_df)
})

contributions = do.call(rbind,contributions)
contributions = reshape2::dcast(contributions, feature_id + sample ~ submodel, value.var = "cell_contribution")
contributions$feature_id = cellMapping$cell_type_name[match(contributions$feature_id, cellMapping$cell_type_id)]

bulks = lapply(bulks, function(x){
  contrib_df = reshape2::melt(x)
  colnames(contrib_df) = c("feature_id","sample","cell_contribution")
  contrib_df$submodel = "bulk"
  return(contrib_df)
}) %>% do.call(rbind,.) %>% reshape2::dcast(., feature_id + sample ~ submodel, value.var = "cell_contribution")

cellContributions = merge(contributions, bulks)
cellContributions.median = cellContributions %>%
  group_by(feature_id) %>%
  summarise(
    across(2:10, ~median(.x, na.rm = TRUE))
  ) %>%
  reshape2::melt(variable.name = "submodel", value.name = "cell_contribution")
mat = reshape2::acast(cellContributions.median, feature_id ~ submodel, value.var = "cell_contribution")
mat = mat[rownames(loadings),]

png(paste0("~/analysis-s05/figures/adjustment/median_cell_contribution_",selectedEffect,".png"), res = "100", bg = "transparent", width = 900, height = 900)
Heatmap(mat, name = "median cell\ncontribution ", column_title = "Median cell contribution", right_annotation = loading_ann)
dev.off()


# 4. cell contribution difference from bulk
# -----------------------------------------------
cellContributionDifferences = cbind(cellContributions[,1:2], cellContributions[,"bulk"] - cellContributions[,3:10])
cellContributionDifferences.sum = cellContributionDifferences %>%
  group_by(feature_id) %>%
  summarise(
    across(2:9, ~median(.x, na.rm = TRUE))
  ) %>%
  reshape2::melt(variable.name = "submodel", value.name = "cell_contribution")

mat = reshape2::acast(cellContributionDifferences.sum, feature_id ~ submodel, value.var = "cell_contribution")
mat = mat[rownames(loadings),]

png(paste0("~/analysis-s05/figures/adjustment/median_change_from_bulk_",selectedEffect,".png"), res = "100", bg = "transparent", width = 900, height = 900)
Heatmap(mat, name = "median\nchange\nfrom bulk", column_title = "Median change in cell contribution from bulk data",
        right_annotation = loading_ann)
dev.off()


# 5. Change in signature genes
# -----------------------------------
cellSignature = SignatureCollection("skin_v12")
cellSignature = cellSignature@gene_set
names(cellSignature) = cellMapping$cell_type_name[match(names(cellSignature),cellMapping$cell_type_id)]
allCells = unique(unlist(cellSignature))

eff = ifelse(selectedEffect == "DZ_vs_HC", "AD", selectedEffect)
gx_diff = statistic_table(ccm_fit_AD$meta[[eff]]$gx_diff$gx_diff)
gx_diff$adjustment_variables[is.na(gx_diff$adjustment_variables)] <- "bulk"
gx_diff = gx_diff[which(gx_diff$feature_id %in% allCells & gx_diff$term == selectedEffect),]
gx_diff$star = sapply(gx_diff$fdr, starryNight)

pushToCC(gx_diff, tagsToPass = list(list(name="object",value="gx_diff_adjustment")))
# wf-cb73e74eff DZ_vs_HC
# wf-085efee950 L_vs_NL
# wf-937da5cfb6 L_vs_HC

mapAdjName = unique(gx_diff[,c("adjustment_variables","submodel")])
# wf-1005306ab9

res_list = lapply(cellSignature, function(cell){
  mat = reshape2::acast(gx_diff[which(gx_diff$feature_id %in% cell),], adjustment_variables ~ feature_id, value.var = "effect_size")
  mat.fdr = reshape2::acast(gx_diff[which(gx_diff$feature_id %in% cell),], adjustment_variables ~ feature_id, value.var = "log10_fdr")
  mat.star = reshape2::acast(gx_diff[which(gx_diff$feature_id %in% cell),], adjustment_variables ~ feature_id, value.var = "star")
  return(list(effect_size = mat, significance = mat.star, fdr = mat.fdr))
})

mapping = c("NK cell" = "natural killer cell","CD8+" = "CD8-positive, alpha-beta T cell","mDC" = "myeloid dendritic cell",
            "Treg" = "regulatory T cell","EM CD4+" = "effector memory CD4-positive, alpha-beta T cell",
            "macrophage" = "alternatively activated macrophage","pDC" = "plasmacytoid dendritic cell","NKT cell" = "mature NK T cell")
col_fun = circlize::colorRamp2(c(-2,0,2),c("blue","white","red"))
idx = which(names(res_list) %in% mapping)
names(res_list)[idx] <- names(mapping)[match(names(res_list)[idx], mapping)] 

ht_list <- lapply(names(res_list), function(cell) {
  Heatmap(res_list[[cell]][["effect_size"]], col = col_fun, show_column_names = F, name = "effect_size",
          column_title = cell, show_row_names = T, column_title_rot = 90, column_title_side = "bottom", cluster_rows = F)
}) %>% Reduce(`+`, .)

png(paste0("~/analysis-s05/figures/adjustment/gene_signature_change_effect_",selectedEffect,".png"), res = "100", bg = "transparent", width = 1800, height = 1000)
draw(ht_list, merge_legend = TRUE, column_title = paste0("Change in cell signature genes as a result of adjustment\nin ",selectedEffect," (effect size)"))
dev.off()

col_fun_fdr = circlize::colorRamp2(c(1.3,2,3,4),c("white","#ff9b7a","#ff522f","#ff0000"))
ht_list_fdr <- lapply(names(res_list), function(cell) {
  Heatmap(res_list[[cell]][["fdr"]], col = col_fun_fdr, show_column_names = F, name = "-log10(fdr)",
          column_title = cell, show_row_names = T, column_title_rot = 90, column_title_side = "bottom", cluster_rows = F)
}) %>% Reduce(`+`, .)

png(paste0("~/analysis-s05/figures/adjustment/gene_signature_change_fdr_",selectedEffect,".png"), res = "100", bg = "transparent", width = 1800, height = 1000)
draw(ht_list_fdr, merge_legend = TRUE, column_title = paste0("Change in cell signature genes as a result of adjustment\nin ",selectedEffect,"-log10(fdr)"))
dev.off()

