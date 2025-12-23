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
#' @param table - a data.frame/data.table, workflow ID or address of file in cytoCC of a table we want to upload, e.g. "wf-3lk24l3nrlk32n:2:path/to/sub-directory/or/file"
#' @param bqdataset - name of dataset in BQ we want to upload to
#' @param tableName - the name we want to call the table in BQ - no special characters other than _
#' @param subsetting - if the wfid outputs a list, this will use the element in the list specified.
#' @param disposition - specific write_disposition specification. Default is "WRITE_TRUNCATE"
uploadToBQ <- function(table, bqdataset, tableName, subsetting = NULL, disposition = "WRITE_TRUNCATE") {
  library(bigrquery)
  
  if (is.character(table)) {
    cat("\nUploading", table, "to", bqdataset, "\n")
    table <- cytoreason.ccm.pipeline::read_asset(table, quiet = TRUE)
    cat("\nData loaded\n")
  }
  
  if (!is.null(subsetting)) {
    table <- table[[subsetting]]
  }
  
  # Clean names
  colnames(table) <- stringr::str_replace_all(colnames(table), "[ |\\.|,]", "_")
  tableName <- stringr::str_replace_all(tableName, "[\\. ]", "_")
  
  if (nrow(table) > 2e6) { # upload in chunks
    parts <- seq(1, nrow(table), by = 1e6)
    jobs <- lapply(parts, function(i) {
      start <- i
      end <- min(i + 1e6 - 1, nrow(table))
      
      cat("\nUploading chunk", start, "-", end)
      disposition_chunk <- if (i == 1) disposition else "WRITE_APPEND"
      job <- bq_perform_upload(bq_table("cytoreason", bqdataset, table = tableName), table[start:end, ], write_disposition = disposition_chunk)
      return(job)
    })
    
  } else { # upload all at once
    job <- bq_perform_upload(bq_table("cytoreason", bqdataset, table = tableName), table, write_disposition = disposition)
    bq_job_wait(job)
    jobs <- list(job)
  }
  
  # Check status of last job
  status <- bq_job_status(jobs[[length(jobs)]])
  
  if (!is.null(status$errorResult)) {
    stop(paste("BigQuery upload failed:", status$errorResult$message))
  } else {
    message("Upload completed successfully.")
  }
  
  gc()
  rm(table)
}



#' @title downloadFromBQ
#' @description downloadFromBQ gets workflow IDs and uploads them to bigquery
#' @param bqdataset - name of dataset in BQ we want to download from
#' @param tableName - the table name we want extract from BQ
downloadFromBQ <- function(bqdataset, tableName){
  library(bigrquery)
  
  tbl <- bq_table("cytoreason", bqdataset, tableName)
  
  out <- tryCatch(
    bq_table_download(tbl),
    error = function(e) stop("BigQuery download failed: ", conditionMessage(e), call. = FALSE)
  )
  return(out)
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
  "Positives:IL13" = "#482870",
  "Positives:IL4" = "#7E55A3",
  "Positives:IL12" = "#B19CD9",
  "Positives:IL17" = "#9575CD",
  "Positives:TSLP" = "#ADD8E6",
  "Positives:TNFA" = "#87CEEB",
  "Positives:IL6" = "#7DDCFF",
  "Positives:IL36" = "#00BFFF",
  "Positives:IL23" = "#4682B4",
  "Positives:CXCL4" = "#0000FF",
  "Positives:IL22" = "#4169E1",
  "Positives:IL33" = "#0D3D8A",
  "Positives:IL31" = "#000080",
  
  # # greens
  # "IL22" = "#c5db87",
  # "IL20" = "#228B22",
  # "IL24" = "#648D62",
  # "IL26" = "#95c74d",
  # "IFNG" = "#1C7771",
  
  "X2:x2_general_inhibition_early_50" = "#E1B1E1",
  "X2:x2_general_inhibition_early_50_archs_refined" = "#c4007c",
  "X2:SP_general_inhibition_late_50" = "#df9aa0",
  "X2:SP_general_inhibition_late_50_archs_refined" = "#CD5C5C",

  # # yellows 
  # "CTLA4" = "#FFF176",
  # "PDL1" = "#FFD933",
  # "PD1" = "#FFD933",
  # "BTLA" = "#bfb64d",
  
  # # oranges
  # "TNFA" = "#FF8C00",
  # "OX40" = "#FFA500",
  # "TNF_mono" = "#FF8C00",
  
  "X2:aIgE_activation_late_50" = "#D2B48C",
  "X2:aIgE_activation_late_50_archs_refined" = "#A0522D",
  
  "negativeControls:BMP4" = "#888888",
  "negativeControls:BMP6" = "#888888",
  "negativeControls:BMP7" = "#888888",
  "negativeControls:TGFB2" = "#888888",
  "negativeControls:random" = "#888888",
  "negativeControls:smoothedRandom_bottom50" = "#888888",
  "negativeControls:smoothedRandom_top50" = "#888888"
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
  theme(strip.background.y = element_blank(), strip.text.y = element_blank()),
  scale_size_continuous(range = c(-1, 10)),
  geom_point(colour = "black", aes(size = 1.1)),
  geom_point(colour = "white", aes(size = 1)),
  xlab(""), ylab(""),
  guides(color = F),
  labs(size = "scaled score\nper criteria")
)


