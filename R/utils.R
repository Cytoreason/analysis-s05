# Utility functions for the project
#
# Copyright Cytoreason (2018)

#' An Example of Exported Function
#'
#' Does not do anything, just serves as an example of an internal function within a project.
#' It is exported because it does __not__ start with a dot.
#'
#' Arguments are checked using [assertthat::assert_that].
#'
#' @param x a `matrix` object with at least 10 rows (see **Details** section)
#' @param y a single positive number
#' @param z a character vector with no empty elements
#' @param verbose a logical that indicates if log messages
#' should be output
#' @param prefix a single character string that specifies the prefix
#' to use.
#'
#' @return an `ExpressionSet` object
#' @examples
#' x <- matrix(1:20, ncol = 2)
#' hello_world(x, 2, "abcd")
#'
hello_world <- function(x, y, z, verbose = TRUE, prefix = "a"){
  # assert arguments
  assert_that(is.matrix(x), nrow(x) >= 10L,
              msg = logf_invalid_argument(x, "expected a matrix with at least 10 rows"))
  assert_that(is.number(y), y > 0,
              msg = logf_invalid_argument(y, "expected a positive number"))
  assert_that(is.character(z), all(nzchar(z)),
              msg = logf_invalid_argument(y, "expected a character vector with no empty element"))
  assert_that(is.logical(verbose))
  assert_that(is.string(prefix))
  ##

  # TODO: do the stuff
  sum(x)

}

#' An Example Internal Function
#'
#' Does not do anything, just serves as an example of an internal function within a project.
#' It is internal because it starts with a dot.
#'
#' Arguments are checked using [assertthat::assert_that].
#'
#' @param x a single character string
#' @param ... other arguments passed to [stats::as.formula]
#'
#' @return a `formula`
#'
.hello_underworld <- function(x, ...){
  # assert arguments
  assert_that(is.string(x),
              msg = logf_invalid_argument(x, "expected a single character string"))
  ##

  # TODO: do the stuff
  stats::as.formula(x, ...)

}
