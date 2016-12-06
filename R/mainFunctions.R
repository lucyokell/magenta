#------------------------------------------------
#' Dummy function
#'
#' This is a dummy function
#'
#' @param x Some parameter
#'
#' @export
#' @examples
#' dummy1()

# The following commands are needed to ensure that the roxygen2 package, which deals with documenting the package, does not conflict with the Rcpp package. Do not alter!
#' @useDynLib MAGENTA
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"

dummy1 <- function() {
    print("dummy1")
}
