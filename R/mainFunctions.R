#------------------------------------------------
#' Dummy function
#'
#' This is a dummy function
#'
#' @param paramList An example list of parameters. Should contain elements "foo" and "bar".
#'
#' @export

# The following commands are needed to ensure that the roxygen2 package, which deals with documenting the package, does not conflict with the Rcpp package. Do not alter!
#' @useDynLib MAGENTA
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"

dummy1 <- function(paramList) {
    
    # check that this function is working
    print("R function is working!")
    
    # do some checks on paramList
    stopifnot(is.list(paramList))
    stopifnot( identical(names(paramList), c("foo", "bar")) )
    
    # ----------------------
    # RUN C CODE
    
    # call Rcpp command with input list
    rawOutput <- dummy1_cpp(paramList)
    
    # ----------------------
    
    # convert rawOutput to final output format
    output_df <- as.data.frame(rawOutput)
    output_df
}
