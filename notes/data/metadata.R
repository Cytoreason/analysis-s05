#' ---
#' title: "Data - Metadata"
#' author: Your Name
#' fontsize: 9pt
#' ---
#'
#' # Overview
#'
#'   * Check and clean clinical data
#'   * Check and clean pehnotypic data
#'
#' ----------------
#+ setup, include = FALSE
library(knitr)
opts_chunk$set(collapse = TRUE)
#
library(devtools)
library(assertthat)
library(cytoreason.shared.assets)
library(cytoreason.analytics)
# load project functions
load_all()

#' # Configuration
#+ config, echo = FALSE
CR_PROJECT <- cr_project_setup()
CR_PROJECT

#' # Fetch raw data from bucket
#+ raw_data
# prepare_raw_data_project()

