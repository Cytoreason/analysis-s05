devtools::load_all("~/analysis-s05/R/utils.R")
library(Seurat)
library(SeuratDisk)
library(tidyverse)

wu_scppp_sue <- "wf-39be8e4885" # all
seurat_object = LoadH5Seurat(get_task_outputs(wu_scppp_sue, task_id = "7",
                                              files_names_grepl_pattern = ".*.h5Seurat")) # Takes 19.2 GB
gc()
pushToCC(seurat_object, tagsToPass = list(list(name="object",value="seurat"),
                                          list(name="dataset",value="wu"),
                                          list(name="notes",value="all_samples")))
# wf-5524ee505a

wu = readRDS(get_workflow_outputs("wf-5524ee505a"))

## Find top genes
## ==============================
remote_find_markers <- function(scppp_wfid, cluster_col, split_df=NULL,
                                filter_keep=NULL){
  devtools::install_github('immunogenomics/presto')
  library(SeuratDisk)
  library(Seurat)
  library(dplyr)
  
  if(is.null(scppp_wfid)){
    seu <- LoadH5Seurat("/cyto_cc/inputs/data/seurat.h5Seurat")
  }else{
    seu <- LoadH5Seurat(cytoreason.cc.client::get_task_outputs(scppp_wfid,
                                                               "7",
                                                               files_names_grepl_pattern = ".*.h5Seurat"))
  }
  if(!is.null(filter_keep)){
    updated_meta <- dplyr::left_join(seu@meta.data,filter_keep)
    row.names(updated_meta) <- row.names(seu@meta.data)
    seu@meta.data<- updated_meta
    subseu<-subset(seu, cells = WhichCells(seu,
                                           expression = keep == T))
    
  }
  if(is.null(split_df)){
    return(Seurat::FindAllMarkers(SetIdent(seu, value=cluster_col)))
  } else {
    res = list()
    names(split_df) <- c(cluster_col, "cells")
    updated_meta <- seu@meta.data
    updated_meta <- dplyr::left_join(seu@meta.data, split_df)
    row.names(updated_meta) <- row.names(seu@meta.data)
    seu@meta.data<- updated_meta
    Idents(seu) <- "cells"
    for (cgroup in split_df[,2] %>% unique()) {
      subseu<-subset(seu, cells = WhichCells(seu, idents = cgroup))
      res[[cgroup]] <- Seurat::FindAllMarkers(SetIdent(subseu, value=cluster_col))
      
    }
    return(res)
  }
}


get_top_deg_per_cluster <- function(clusters_diff, cname=NA,  top_x = 100){
  translation <-  AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
  if(nrow(clusters_diff) <1){return(data.frame())}
  res <- clusters_diff %>% 
    dplyr::filter(p_val_adj <= 0.05) %>%
    dplyr::group_by(cluster)%>%
    dplyr::mutate(grank = dplyr::row_number(-avg_log2FC),
                  gene_symbol = translation$symbol[match(gene, translation$gene_id)]) %>%
    filter(grank <= top_x)
  if(!is.na(cname)) res <- res %>% mutate(cell_type = cname)
  return(res)
}


get_sigs_per_cluster <- function(clusters_diff, top_x = 50){
  if(nrow(clusters_diff) <1){return(c())}
  sigs<- clusters_diff %>% 
    dplyr::group_by(cluster)%>%
    dplyr::mutate(grank = dplyr::row_number(avg_log2FC * log10(p_val_adj))) %>%
    dplyr::filter(grank <= top_x) %>%
    dplyr::arrange(grank, .by_group = T)
  sigs_list <- list()
  for(cluster in sigs$cluster %>% unique()){
    sigs_list[[cluster]] <- sigs$gene[sigs$cluster == cluster]
  }
  
  return(sigs_list)
}


run_function_dist(remote_find_markers,
                  scppp_wfid="wf-39be8e4885", 
                  cluster_col="local_supercluster",
                  filter_keep=NULL,
                  memory_request = '150Gi',
                  image = "eu.gcr.io/cytoreason/ci-cytoreason.single.cell.preprocessing-package:upgrade_R4_4_latest",
                  inputs = list(input_wf_output(res = "wf-39be8e4885",
                                                task_id = "7",
                                                prefix = "seurat.h5Seurat",
                                                target_dir = "data")),
                  tags=list(
                    list(name = "operation", value = "find clusters markers"),
                    list(name = 'condition', value = 'PSO'),
                    list(name = "dataset", value = "wu")),
                  force_execution = T,
                  replace_image_tags = TRUE,
                  data_access = 's05')

clusters_diff <- readRDS(get_workflow_outputs('wf-2dc742c3d0'))
clusters_degs <- get_top_deg_per_cluster(clusters_diff)
write.csv(clusters_degs %>% dplyr::mutate(dataset="Wu", cell_type="T Cells"),
          "~/analysis-s05/data/topDEG.Wu.csv", row.names = F)

pushToCC(clusters_degs, tagsToPass = list(list(name="object",value="DEGs_perCell_Wu")))
# wf-f66fc3b94d

mast = clusters_degs %>%
  dplyr::filter(cluster == "Mast") %>%
  slice_max(order_by = avg_log2FC, n = 50)

pushToCC(mast$gene, tagsToPass = list(list(name="object",value="mast_top50_wu")))
# wf-66ba2c4346