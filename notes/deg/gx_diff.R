#' ---
#' title: Gene Expression Differential Analysis
#' author: Your Name
#' fontsize: 9pt
#' output:
#'   html_document:
#'     code_folding: hide
#' ---
#'
#' # Overview
#'  * Performs gene expression differences
#'
#' ----------------
#+ setup, include = FALSE
library(knitr)
opts_chunk$set(collapse = TRUE, warning = FALSE)
# packages
library(devtools)

# load project
load_all()
