devtools::load_all("~/analysis-s05/R/utils.R")
library(cytoreason.cc.client)
library(tidyverse)

geneMapping = toTable(org.Hs.eg.db::org.Hs.egSYMBOL)

### 1. X2 signatures
### =======================
<<<<<<< Updated upstream
x2 = readRDS(get_workflow_outputs("wf-06b8b0617e"))
x2 = x2[str_detect(names(x2), "50_p|nc_50|smoothedRandom")]
x2 = x2[!str_detect(names(x2), "string|archs(?!_)|p_refined|top100|bottom100|x2_acti")]
x2 = x2[!str_detect(names(x2), "IgE_inhibition|CST14_inhibition|Icatibant_inhibition|PAMP12_inhibition|SP_inhibition_early|Untreated_inhibition")]
x2 = x2[!str_detect(names(x2),"(?:CST14|Icatibant|PAMP12|aIgE).*?_p_[^ ]+")]
x2 = x2[!str_detect(names(x2),"^x2_(?:activation|inhibition)_late.*$")]

pushToCC(x2, tagsToPass = list(list(name="object",value="x2_signatures")))
# wf-55420ff1cb
=======
x2 = readRDS(get_workflow_outputs("wf-50aa0d34a6"))
# x2 = x2[str_detect(names(x2), "50_p|nc_50|smoothedRandom|_50_")]
# x2 = x2[!str_detect(names(x2), "string|archs(?!_)|p_refined|top100|bottom100")]
# x2 = x2[!str_detect(names(x2), "IgE_inhibition|CST14_inhibition|Icatibant_inhibition|PAMP12_inhibition|SP_inhibition_early|Untreated_inhibition|SP_activation_early_50_p_archs_refined")]
# x2 = x2[!str_detect(names(x2),"(?:CST14|Icatibant|PAMP12|aIgE).*?_p_[^ ]+")]
# x2 = x2[!str_detect(names(x2),"^x2_(?:activation|inhibition)_late.*$")]
# x2 = x2[!str_detect(names(x2),"_ep|50_refined")]
x2 = x2[str_detect(names(x2), "_50|MRGPRX|top50|bottom50")]
x2 = x2[!str_detect(names(x2), "string")]
x2 = x2[!str_detect(names(x2), "_archs(?!_)|50_refined")]

pushToCC(x2, tagsToPass = list(list(name="object",value="x2_signatures")))
# wf-1b6040f083
# wf-f328b4a739
>>>>>>> Stashed changes

### 2. in-silico signatures
### ============================
sigs = c("CORT","SST","TAC1","ADCYAP1","PDYN","NPFF","VIP","CAMP","CHGA","DEFB103A","DEFB103B","IGFBP5",
         "AVP","OXT","TAC4","CNP","TPM1","HSPE1","ALB","ADM","EPX","RNASE3","PRG2","CXCL14","PF4")
insilico <- cytoreason.assets::read_asset("ccw://wf-9329a6bc38@0:signatures_18803")
insilico = insilico[which(names(insilico) %in% sigs)] # DEFB103A not in the list

all_ligands = append(list(all_ligands = sigs), insilico)
all_ligands = geneMapping$gene_id[match(all_ligands, geneMapping$symbol)]
pushToCC(all_ligands, tagsToPass = list(list(name="object",value="ligand_signatures")))
# wf-f7e05ae87c


### 3. signature bank
### ============================

signatureBank = bigrquery::bq_table_download(x = bigrquery::bq_table(project = "cytoreason",
                                                                     dataset = "p00_proj_patern_benchmark_public",
                                                                     table="p00_v2_signatures"))

bankMapping = bigrquery::bq_table_download(x = bigrquery::bq_table(project = "cytoreason",
                                                                   dataset = "p00_proj_patern_benchmark_public",
                                                                   table="all_sigs_origin_v2"))
bankMapping <- bankMapping %>%
  dplyr::filter(project %in% c("p03", "CytoSig_50Genes", "p07", "p01")) %>%
  dplyr::filter(str_detect(target_collection, "Undivided|Validations|CytoSig|BioNetSmooth_signatures")) %>%
  dplyr::filter(
    (project == "p03" & target %in% c("BMP7", "IL33", "IL17", "IL36", "IL31")) |
      (project == "CytoSig_50Genes" & target %in% c("BMP4", "BMP6", "TGFB2", "IL4", "IL13", "IL22", "TNFA", "IL12", "TSLP", "IL6",
                                                    "BDNF", "NO", "NRG1", "NTF4")) |
      (project == "p07" & target %in% c("IL23")) |
      (project == "p01" & target %in% c("CXCL4", "NGF"))
  )
signatures = unique(bankMapping$geneset)
signatures = signatureBank[which(signatureBank$signature %in% signatures),]
signatures$entrez = geneMapping$gene_id[match(signatures$gene, geneMapping$symbol)]

entrez_lookup <- c(PHB = "5245", KIAA0100 = "9703", ARNTL2 = "56938", ANP32C = "23520", FAM126A = "84668")

signatures <- signatures %>%
  mutate(entrez = if_else(is.na(entrez), entrez_lookup[gene], entrez))

signatures = split(signatures$entrez, signatures$signature)
pushToCC(signatures, tagsToPass = list(list(name="object",value="positiveControl_signatures")))
# wf-49eee6a885
# wf-3abe7c76cd

### 4. Itch signatures
### ===============================
library(biomaRt)

itch = openxlsx::read.xlsx("~/analysis-s05/data/itch_genesets.xlsx")
itch = as.list(itch)
itch = lapply(itch, function(x) x[!is.na(x)])
itch = lapply(itch, function(x) str_trim(x,"both"))
itch = itch[which(names(itch) %in% c("Yosipovitch.HM.all","Yosipovitch.HM.AD.FC1.5","Yosipovitch.HM.AD.FC2",
                                     "Yosipovitch.HM.AD.FC2.5", "Yosipovitch.HM.AD.FC3", "Jha.Sanofi.mouseGenes.skin"))]

# translate the mouse genes to human
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

homologs <- getLDS(attributes = c("mgi_symbol", "entrezgene_id"),
                   filters = "mgi_symbol",
                   values = itch$Jha.Sanofi.mouseGenes.skin,
                   mart = mouse,
                   attributesL = c("hgnc_symbol", "entrezgene_id"),
                   martL = human)
homologs$NCBI.gene..formerly.Entrezgene..ID.1[26] = 442197
homologs$NCBI.gene..formerly.Entrezgene..ID[1:2] = 100416956
itch$Jha.Sanofi.mouseGenes.skin = unique(homologs$HGNC.symbol)
itch = lapply(itch, function(x) unique(geneMapping$gene_id[match(x, geneMapping$symbol)]))

pushToCC(itch, tagsToPass = list(list(name="object",value="itch_signatures")))
# wf-16fa7ba6c0


### 5. Mast cell signatures
### ===============================
# 5.1. CR AirTable
# --------------------
mast_at = c("2182","27306","1359","64499","1511","2624","23430","4013","7177","9173","3815","4056","2206","85414","3067","374","6571","3248")


# 5.2. CRS
# --------------------
mast_crs <- cytoreason.deconvolution::SignatureCollection("gut_v9")
mast_crs = mast_crs@gene_set[["CRCL_0000009"]]

mast = list(mast_crs = mast_crs, mast_at = mast_at)


pushToCC(mast, tagsToPass = list(list(name="object",value="mast_signatures")))
# wf-850ac449a7


### 6. Pooling all signatures
### ================================
<<<<<<< Updated upstream
x2 = readRDS(get_workflow_outputs("wf-55420ff1cb")) %>% lapply(., as.character)
ligands = readRDS(get_workflow_outputs("wf-f7e05ae87c"))
positives = readRDS(get_workflow_outputs("wf-49eee6a885")) %>% lapply(., as.character)
itch = readRDS(get_workflow_outputs("wf-16fa7ba6c0"))
mast = readRDS(get_workflow_outputs("wf-850ac449a7"))

allSignatures = list(x2 = append(x2, append(ligands,positives)),
                     mast = mast,
                     itch = itch)

pushToCC(allSignatures, tagsToPass = list(list(name="object",value="allSignatures")))
# wf-b33649db09
=======
x2 = readRDS(get_workflow_outputs("wf-f328b4a739")) %>% lapply(., as.character)
ligands = readRDS(get_workflow_outputs("wf-9ea3eeb42c"))
positives = readRDS(get_workflow_outputs("wf-3abe7c76cd")) %>% lapply(., as.character)
itch = readRDS(get_workflow_outputs("wf-5f6ca4a13c"))
mast = readRDS(get_workflow_outputs("wf-850ac449a7"))

nc = x2[str_detect(names(x2),"nc_|Random")]
  names(nc) = str_remove(names(nc),"nc_50.")
x2 = x2[!str_detect(names(x2),"nc_|Random")]
names(x2)[str_detect(names(x2),"Tryptase_50_smoothed")] <- "mast_tryptase"
x2 = x2[!str_detect(names(x2),"x2_inhibition_late")]

names(positives) = str_remove(names(positives), "CytoSig_50Genes__Undivided__|P01__BioNetSmooth_signatures__|P03__Validations__|P07__CytoSig__")

X2 = c(x2, ligands, positives, itch, mast)

allSignatures = list(X2 = X2,
                     negativeControls = nc)

pushToCC(allSignatures, tagsToPass = list(list(name="object",value="allSignatures")))
# wf-dffe62df1b
# wf-54d564beb9
>>>>>>> Stashed changes

signatures.long = reshape2::melt(allSignatures)
signatures.long = signatures.long[-which(duplicated(signatures.long[,2:3])),3:2] # keep only the signature name and collection
colnames(signatures.long) = c("Target.Collection","Target.Identifier")
signatures.long$feature_id = paste0(signatures.long$Target.Collection, "__", signatures.long$Target.Identifier)
<<<<<<< Updated upstream
pushToCC(signatures.long, tagsToPass = list(list(name="object",value="targetMapping.V1"),
                                            list(name="project",value="evo")))
# wf-f90b0697d4
=======
pushToCC(signatures.long, tagsToPass = list(list(name="object",value="targetMapping.V3"),
                                            list(name="project",value="evo")))
# wf-3e93b9331e
# wf-dd6811aa54

collectionMapping = lapply(names(list(x2=x2,ligands=ligands, positives=positives, itch=itch,mast=mast, nc=nc)), 
       function(collection){
  return(data.frame(signature = names(get(collection)), collection = str_to_title(collection)))
}) %>% do.call(rbind,.)
collectionMapping$collection[which(collectionMapping$signature == "mast_tryptase")] <- "Mast"
collectionMapping$collection[which(collectionMapping$collection == "Nc")] <- "Negative Controls"
collectionMapping = rbind(collectionMapping, data.frame(signature = c("smoothedRandom_top50",
                                                 "smoothedRandom_bottom50",
                                                 "random"),
                                               collection = rep("Negative Controls",3)))
pushToCC(collectionMapping)
# wf-e24aa55d1c
# wf-77fec945d8
>>>>>>> Stashed changes
