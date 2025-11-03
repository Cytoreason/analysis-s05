######## Basic functions
######## -----------------------------------------------------------------
#' @title starryNight
#' @description starryNight is a function that assigns significance stars to p values.
#' @param value - a single p value.
#' @return for value>0.05, returns an empty character. for value < 0.001 returns ***
#' for value < 0.01 returns **, for value < 0.05 returns *
#' @export
starryNight = function(value) {
  if(is.na(value)) { star = ""
  } else if(value <= 0.001) { star = "***"
  } else if (value <= 0.01) { star = "**"
  } else if (value <= 0.05) { star = "*"
  } else { star = "" }
  return(star)
}

#' @title pushToCC
#' @description pushToCC uploads a local object to cytoCC
#' @return a workflow ID for the object
#' @export
pushToCC <- function(obj, filename = "output", data_access = "s05", tagsToPass = list(), memRequest = "10Mi"){
  # saving the output file to cyto-cc
  cytoreason.cc.client::run_command_dist(
    command = sprintf('cp /cyto_cc/inputs/obj.rds "output/%s.rds"', filename),
    image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest",
    inputs = list(input_local_obj(obj, 'obj.rds')),
    memory_request = memRequest, 
    tags = tagsToPass,
    data_access = data_access
  )
}

#' @title uploadToBQ
#' @description uploadToBQ gets workflow IDs and uploads them to bigquery
#' @note important: make sure there are no special characters anywhere in table column names or in the table or else it will fail.
#' The only allowed special character is _.
#' @note a table bigger than 1 million rows will be uploaded in chunks of 1 million rows
#' @param wfid - a data.frame/data.table, workflow ID or address of file in cytoCC of a table we want to upload, e.g. "wf-3lk24l3nrlk32n:2:path/to/sub-directory/or/file"
#' @param bqdataset - name of dataset in BQ we want to upload to
#' @param tableName - the name we want to call the table in BQ - no special characters other than _
#' @param subsetting - if the wfid outputs a list, this will use the element in the list specified.
#' @param disposition - specific write_disposition specification. Default is "WRITE_TRUNCATE"
uploadToBQ = function(wfid, bqdataset, tableName, subsetting = NULL, disposition = "WRITE_TRUNCATE"){
  # if(!require(bigquery)) BiocManager::install("bigquery")
  library(bigrquery)
  # options(cache = TRUE, cache.lazy = FALSE)
  
  if(is.character(wfid)){
    cat("\nUploading",wfid,"to",bqdataset,"\n")
    table = cytoreason.ccm.pipeline::read_asset(wfid, quiet = T)
  } else {
    table = wfid
  }
  cat("\nData loaded\n")
  if(length(subsetting == "NULL") == 0) { subsetting = NULL }
  
  if(!is.null(subsetting)){
    table = table[[subsetting]]
  }
  if(any(stringr::str_detect(colnames(table), " |\\.|,"))){
    colnames(table) = stringr::str_replace_all(colnames(table), pattern = " ", replacement = "_")
    colnames(table) = stringr::str_replace_all(colnames(table), pattern = "\\.", replacement = "_")
    colnames(table) = stringr::str_replace_all(colnames(table), pattern = ",", replacement = "")
  }
  tableName = stringr::str_replace_all(tableName, pattern = "\\.", replacement = "_") # make sure there are no . or whitespace
  tableName = stringr::str_replace_all(tableName, pattern = " ", replacement = "_") # make sure there are no . or whitespace
  parts = seq(1,nrow(table),by=1e6)
  if(nrow(table) > 2e6){ # upload in chunks
    sapply(parts, function(i){
      cat("\nuploading chunk",i,"-",(i+1e6-1))
      if(i == 1){
        bq_perform_upload(bq_table("cytoreason",bqdataset,table=tableName),
                          table[i:(i+1e6-1),], write_disposition = disposition)
      } else if(i == parts[length(parts)]) { # last part
        cat("\ruploading chunk",i,"-",nrow(table))
        bq_perform_upload(bq_table("cytoreason",bqdataset,table=tableName),
                          table[i:nrow(table),], write_disposition = "write_append")
      } else {
        bq_perform_upload(bq_table("cytoreason",bqdataset,table=tableName),
                          table[i:(i+1e6-1),], write_disposition = "write_append")
      }
      gc(verbose = F)
    })
  } else { # upload all at once
    bq_perform_upload(x = bq_table("cytoreason",bqdataset,table=tableName), values = table, write_disposition = disposition)
  }
  cat("\nData uploaded\n")
  gc()
  rm(table)
}



# squash_axis for ggplot
squash_axis <- function(from, to, factor) { 
  trans <- function(x) { 
    isq <- x > from & x < to & !is.na(x)
    ito <- x >= to & !is.na(x)
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    return(x)
  }
  
  inv <- function(x) {
    isq <- x > from & x < from + (to - from)/factor & !is.na(x)
    ito <- x >= from + (to - from)/factor & !is.na(x)
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    return(x)
  }
  
  return(scales::trans_new("squash_axis", trans, inv))
}

### Creating targetColors
targetColors = c(
  "Dupilumab" = "#482870",
  "Dupilumab50" = "#482870",
  "smoothedSteady_100" = "#482870",
  "smoothedSteady_50" = "#482870",
  "steady_100" = "#482870",
  "steady_50" = "#482870",
  "IL4R" = "#482870",
  "IL13.Dupi" = "#7E55A3",
  "IL4.Dupi" = "#512DA8",
  "IL13" = "#B19CD9",
  "IL13_mono" =  "#B19CD9",
  "IL13fdr" =  "#B19CD9",
  "IL13cs" =  "#B19CD9",
  "IL4" = "#9575CD",
  "IL4_mono" = "#9575CD",
  
  "TSLP" = "#ADD8E6",
  "IL23" = "#87CEEB",
  "IL27" = "#7DDCFF",
  "IL6" = "#00BFFF",
  "IL2" = "#4682B4",
  "IL7" = "#0000FF",
  "IL15" = "#4169E1",
  "IL31" = "#0D3D8A",
  "IL21" = "#000080",
  
  "IL22" = "#c5db87",
  "IL20" = "#228B22",
  "IL24" = "#648D62",
  "IL26" = "#95c74d",
  "IFNG" = "#1C7771",
  
  "IL33" = "#df9aa0",
  "IL36" = "#CD5C5C",
  "IL17" = "#c4007c",
  "IL18" = "#E53935",
  "CCL17" = "#B22222",
  "GPR15L" = "#8B0000",
  
  "CTLA4" = "#FFF176",
  "PDL1" = "#FFD933",
  "PD1" = "#FFD933",
  "BTLA" = "#bfb64d",
  
  "TNFA" = "#FF8C00",
  "OX40" = "#FFA500",
  
  "KLK7" = "#D2B48C",
  "KLK5" = "#D2B48C",
  "PAR2" = "#7B3F00",
  "SLAMF6" = "#A0522D",
  
  "BMP7" = "#888888",
  "IL1A" = "#888888",
  "placeboNC" = "#888888",
  "random" = "#888888",
  "smoothedRandom" = "#888888",
  "random-fc" = "#888888",
  "ad-permutations" = "#888888",
  
  "IL12IL23" = "#E1B1E1",
  "IL23TNF" = "#E1B1E1",
  "IL12IL23IL17" = "#E1B1E1",
  "IL12_mono" = "#E1B1E1",
  "IL12IL23_IL17" = "#E1B1E1",
  "IL23_mono" = "#87CEEB",
  "TNF_mono" = "#FF8C00"
)
# pushToCC(targetColors) # wf-1f18ee71f0, wf-da4c738906


## Spot graph for criteria
library(ggplot2)
spot.theme <- list(
  theme_classic(),
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)),
  theme(axis.ticks.y=element_blank(),axis.text.y=element_text()),
  theme(axis.line=element_blank()),
  theme(text = element_text(size = 14)),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  scale_size_continuous(range = c(-1, 10)))