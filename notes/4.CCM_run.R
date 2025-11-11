library(cytoreason.ccm.pipeline)
library(cytoreason.cc.client)
library(tidyverse)
library(reshape2)
devtools::load_all("~/analysis-s05/R/utils.R")

ccm = as_ccm_fit("wf-08a6a0a503")

# add itch signatures, reference, tryptaser
# collection by number of genes?
# for 24hr we should take SP activation