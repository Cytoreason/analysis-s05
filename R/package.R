# Global package functions or configuration
#
# Copyright Cytoreason (2018)

# Global imports: put here the packages from which you want to import **all** the functions in one go
#' @import stats devtools Biobase cytoreason cytoreason.io cytoreason.shared.assets assertthat
#' @rawNamespace exportPattern("^[^\\.]")
NULL

# Hook that is executed everytime the package is loaded.
# Put here things like global options, etc...
.onLoad <- function(...){
    # set project options
    options(stringsAsFactors = FALSE)

}
