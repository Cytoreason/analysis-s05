## In this part we run a CCM pipeline of custom gene set in order to incorporate the Th2, itch, neuronal pathways into the collections.
## Then we will use this CCM run as AssetData and correlate signatures to it

devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.ccm.pipeline)
library(tidyverse)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)

ccm = as_ccm_fit("wf-08a6a0a503")
geneMapping = toTable(org.Hs.eg.db::org.Hs.egSYMBOL)


## 1. Th2 Criterion
## =============================
th2 = openxlsx::read.xlsx("~/analysis-s05/data/Type 2 inflammation signatures.xlsx", sheet = 2)
th2 = as.list(th2)
th2 = lapply(th2, function(x) x[!is.na(x)])
th2 = lapply(th2, function(x) str_trim(x,"both"))
th2_pathways = th2$Th2_pathways
th2 = th2[-which(names(th2) == "Th2_pathways")]

# 1.1. Mapping GO pathways to four gene sets
# ---------------------------------------------------
th2_go <- th2[c("GO_IgE","GO_Type2General","GO_Th2_DifferentiationAndActivation","GO_Th2_Cytokines")]
th2 = th2[-which(names(th2) %in% c("GO_IgE","GO_Type2General","GO_Th2_DifferentiationAndActivation","GO_Th2_Cytokines"))]

ids <- unlist(th2_go, use.names = FALSE)
names(ids) <- rep(names(th2_go), lengths(th2_go))
ids = sub("^(GO:\\d+).*$", "\\1", ids)

genes_for_go <- select(org.Hs.eg.db,
                       keys = ids,
                       columns = c("ENTREZID", "SYMBOL", "GENENAME"),
                       keytype = "GO")
genes_for_go$collection = names(ids)[match(genes_for_go$GO, ids)]

go_term <- select(GO.db,
                  keys = ids,
                  columns = "TERM",
                  keytype = "GOID")

go = merge(go_term, genes_for_go, by.x = "GOID", by.y = "GO")
go = go[!is.na(go$SYMBOL),]
go_list = split(go$SYMBOL, go$collection)
lengths(go_list)
# GO_IgE  GO_Th2_Cytokines  GO_Th2_DifferentiationAndActivation GO_Type2General
#     15                78                                   18              15

go_list = lapply(go_list, unique)
lengths(go_list)
# GO_IgE  GO_Th2_Cytokines  GO_Th2_DifferentiationAndActivation GO_Type2General
#     10                58                                   16              14



# 1.2. Adding Lichenification geneset
# ------------------------------------------
lich = list(lichenification = c("FLG", "SPINK5", "FLG2", "SPRR3", "CLDN1"))


# 1.3. Mouse/Human translation
# ----------------------------------------
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

homologs <- getLDS(attributes = c("mgi_symbol", "entrezgene_id"),
                   filters = "mgi_symbol",
                   values = th2$Tibbitt,
                   mart = mouse,
                   attributesL = c("hgnc_symbol", "entrezgene_id"),
                   martL = human)

th2$Tibbitt = unique(homologs$HGNC.symbol)


# 1.4. Incorporate existing Th2 pathways for better interpretability of NES/FDR
# -----------------------------------------------------------------------------------
canonical = cytoreason.gx::load_gene_set_collections(collection = c("kegg","btm","reactome"))

th2_pathways = str_remove(th2_pathways,"reactome_|btm_|kegg_")
th2_pathways = lapply(canonical, function(x) x[which(names(x) %in% th2_pathways)]) %>% unlist(recursive = F)
names(th2_pathways) = str_replace(names(th2_pathways),"\\.",":")


# 1.4. Merge
# ---------------------------
th2 = Reduce(append, list(th2, go_list, lich, th2_pathways))

th2[-c(18:24)] = lapply(th2[-c(18:24)], function(x) geneMapping$gene_id[match(x, geneMapping$symbol)])
any(lapply(th2, function(x) any(is.na(x))) %>% do.call(c,.))
# FALSE

pushToCC(th2, tagsToPass = list(list(name="object",value="th2_criterion")))
# wf-e02cfb9673
# wf-410536ebd3



## 2. Epidermal integrity Criterion
## ====================================
canonical = cytoreason.gx::load_gene_set_collections(collection = c("kegg","btm","reactome","h"))
epi = c("alpha-linolenic (omega3) and linoleic (omega6) acid metabolism",
        "alpha-linolenic acid (ALA) metabolism",
        "Cholesterol biosynthesis",
        "Fatty acid metabolism",
        "Fatty acyl-CoA biosynthesis",
        "Glyoxylate and dicarboxylate metabolism",
        "HALLMARK_ADIPOGENESIS",
        "HALLMARK_FATTY_ACID_METABOLISM",
        "HALLMARK_PEROXISOME",
        "Linoleic acid (LA) metabolism",
        "Metabolism of steroids",
        "Mitochondrial Fatty Acid Beta-Oxidation",
        "mitochondrial fatty acid beta-oxidation of saturated fatty acids",
        "Peroxisomal lipid metabolism",
        "Peroxisomal protein import",
        "Peroxisome",
        "PPAR signaling pathway",
        "Asymmetric localization of PCP proteins",
        "Differentiation of keratinocytes in interfollicular epidermis in mammalian skin")


epi_pathways = lapply(canonical, function(x) x[which(names(x) %in% epi)]) %>% unlist(recursive = F)
names(epi_pathways) = str_replace(names(epi_pathways),"\\.",":")

pushToCC(epi_pathways, tagsToPass = list(list(name="object",value="epidermal_barrier_criterion")))
# wf-d3b177cd82


## 3. Neuroinflammation criterion
## ============================================
# Taken from Final_Signatures.R
collectionMapping = readRDS(get_workflow_outputs("wf-83b199630d"))
  collectionMapping = collectionMapping[which(collectionMapping$collection %in% c("Neuronal","Itch")),]

allSignatures = readRDS(get_workflow_outputs("wf-b1950a97bd"))

neuro = allSignatures$X2[which(names(allSignatures$X2) %in% collectionMapping$signature)]
pushToCC(neuro, tagsToPass = list(list(name="object",value="neuro_criterion")))
# wf-4ce41a599a


## 4. Reordering some criterions pre-run
## ===========================================
th2 = readRDS(get_workflow_outputs("wf-410536ebd3"))
neuroinflammation = readRDS(get_workflow_outputs("wf-4ce41a599a"))
epidermis = readRDS(get_workflow_outputs("wf-d3b177cd82"))

epidermis = c(epidermis, th2[c("lichenification","Terminal_Differentiation_and_Lipids")])
th2 = th2[-which(names(th2) %in% c("lichenification","Terminal_Differentiation_and_Lipids","Th1_Related","Th17_Related","TH22_IL22_Related"))]

## 5. Gene set activity calculation
## ===========================================
set.seed(1234)
criteria = list(th2 = th2,
                neuroinflammation = neuroinflammation,
                epidermis = epidermis)
criteria = lapply(criteria, function(x) lapply(x, unique))

datasets = lapply(ccm$datasets, function(x) {
  ccm_service_gene_set_activity(ccm_dataset = x, 
                                custom_gene_set_collections = criteria, 
                                assay_data_element = c("exprs","exprs_adjusted__1__1"))
})

enrichments = lapply(datasets, function(x) readRDS(get_workflow_outputs(x[["object"]][["fits"]][["workflow_id"]], files_names_grepl_pattern = "output.rds")))
# results are different in values than the original ones but the rank is identical. differences are probably due to seed or something similar.
# notice this is on cytocc so might take a little time

enrichments = lapply(enrichments, function(x){
  lapply(x$assay_data, function(submodel){
    d = submodel[c("th2","neuroinflammation","epidermis")]
    d = lapply(d, function(x) x[["score"]])
  })
})

enrichments = lapply(enrichments, purrr::transpose)
enrichments = purrr::transpose(enrichments)
pushToCC(enrichments, tagsToPass = list(list(name="object",value="th2_enrichments")))
# wf-6f7bab48d9
# wf-7ed6eff409 - rearrange lists

enrich_exprssionSet_bulk = lapply(enrichments, function(x) lapply(x, function(y) ExpressionSet(y[["exprs"]])))
enrich_exprssionSet_adj = lapply(enrichments, function(x) lapply(x, function(y) ExpressionSet(y[["exprs_adjusted__1__1"]])))


## 6. Meta PC calculation
## ===========================================
library(cytoreason.integration)

config = modelMetadata(ccm, format = "config")

calc_metaPCs = function(pathways) {
  # Get the data for the meta-PCA train datasets (as defined in the config file):
  pathways.train <- pathways[config$dataset_id[which(config$ccm_meta_pca == 1)]]
  
  # Train the meta PCA:
  metaPCA <- service_meta_pca(data = pathways.train, n_pc = 3)
  
  # Project the samples (all the datasets) on the meta PCA
  metaPCA.projected <- lapply(pathways, function(x) service_meta_pca_projection(x, meta_pca = metaPCA))
  
  return(list(metaPCA = metaPCA, metaPCA.projected = metaPCA.projected))
}

metaPCs_bulk = lapply(enrich_exprssionSet_bulk, calc_metaPCs)
# wf-01abdfcd8d
# wf-4ced0fb3d9 - rearrange lists
metaPCs_adj = lapply(enrich_exprssionSet_adj, calc_metaPCs)
# wf-0803da8862
# wf-89fcefafc6 - rearrange lists


# Extract sample scores
extract_sampleScores = function(res) {
  pathwayPCA <- cytoreason.ccm.pipeline:::statistic_table.ccm_service_meta_pca_projection(res$metaPCA.projected)
  sampleScores <- pathwayPCA$sample_loadings
  sampleScores_mat <- reshape(sampleScores[,c("pc","sample_id","value")], idvar = c("sample_id"), timevar = "pc", direction = "wide")
  colnames(sampleScores_mat) = c("sample_id","Meta_PC1","Meta_PC2","Meta_PC3")
  return(sampleScores_mat)
}

sampleScores_bulk = lapply(metaPCs_bulk, extract_sampleScores)
# wf-a7ac64d7d6
# wf-e6c43f2c55 - rearrange lists
sampleScores_adj = lapply(metaPCs_adj, extract_sampleScores)
# wf-f76769a6b6
# wf-40fbd742a3 - rearrange lists


# Extract pathway loadings
extract_pathwayLoadings = function(res){
  metaPCA.res = cytoreason.ccm.pipeline:::statistic_table.ccm_service_meta_pca(list(res$metaPCA))
  pathLoadings = metaPCA.res$feature_loadings[,c("pc","feature_id","value")]
  colnames(pathLoadings) = c("PC","Pathway","Loading")
  pathLoadings$PC = str_replace(pathLoadings$PC, "pc","PC")
  return(pathLoadings)
}

pathwayLoadings_bulk = lapply(metaPCs_bulk, extract_pathwayLoadings)
  pathwayLoadings_bulk$th2$Loading[which(pathwayLoadings_bulk$th2$PC == "PC1")] = (-1) * pathwayLoadings_bulk$th2$Loading[which(pathwayLoadings_bulk$th2$PC == "PC1")]
pathwayLoadings_adj = lapply(metaPCs_adj, extract_pathwayLoadings)
  pathwayLoadings_adj$th2$Loading[which(pathwayLoadings_adj$th2$PC == "PC1")] = (-1) * pathwayLoadings_adj$th2$Loading[which(pathwayLoadings_adj$th2$PC == "PC1")]


# process for BQ
process = function(df, submodel, term, collection){
  df$Collection = collection
  df$Term = term
  df$Submodel = submodel
  return(df)
}

pathwayLoadings_bulk = lapply(names(pathwayLoadings_bulk), function(x) process(pathwayLoadings_bulk[[x]], "bulk", NA, x)) %>% do.call(rbind,.)
# wf-9ba68682c9
# wf-60ab388266 - rearrange lists
pathwayLoadings_adj = lapply(names(pathwayLoadings_adj), function(x) process(pathwayLoadings_adj[[x]], "adjusted__1__1", NA, x)) %>% do.call(rbind,.)
# wf-9cc5e483e3
# wf-dde17f906e - rearrange lists

uploadToBQ(pathwayLoadings_bulk, bqdataset = "s05_atopic_dermatitis", tableName = "pathwayLoadings", disposition = "WRITE_APPEND")
uploadToBQ(pathwayLoadings_adj, bqdataset = "s05_atopic_dermatitis", tableName = "pathwayLoadings", disposition = "WRITE_APPEND")
