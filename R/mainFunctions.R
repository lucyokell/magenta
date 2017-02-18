#------------------------------------------------
#' Parameter List creation for MAGENTA simulation
#'
#' \code{Param_List_Simulation_Create} creates suitable parameter list for
#' \code{Simulation_R}
#'
#' @param N Population size. Default = 1e4
#' @param years Length of simulation. Default = 1
#' @param eqSS Output of \code{Equilibrium_Steady_State_Create}
#' 
#' @export

Param_List_Simulation_Create <- function(N = 1e+04, years = 1, eqSS)
{
  
  ## CHECKS ##
  ##---------------------------------------------
  if(!identical(names(eqSS),
               c("age_brackets","het_brackets","Smat",
                                              "Dmat","Amat","Umat","Tmat","Pmat",
                                              "IBmat","ICAmat","ICMmat","IDmat",
                                              "Iv","MaternalImmunity"))) stop("Incorrect variable names within equilibrium.steady.state")
  dims.1 <- lapply(eqSS,function(x){return(dim(x)[1])})
  dims.2 <- lapply(eqSS,function(x){return(dim(x)[2])})
  if(unique(unlist(dims.1[grep("mat",names(eqSS))]))!=length(eqSS$age_brackets)) stop("Dimensions 1 error in equilibrium.stead.state matrices")
  if(unique(unlist(dims.2[grep("mat",names(eqSS))]))!=length(eqSS$het_brackets)) stop("Dimensions 2 error in equilibrium.stead.state matrices")
  
  numerics <- which(!unlist(lapply(eqSS,is.numeric)))
  if (length(numerics)!=0) stop(paste(names(numerics),"provided not numeric"))
  
  ##---------------------------------------------
  
  # Create paramlist
  paramList <- list(N = N, years = years, eqSS = eqSS)
  
  return(paramList)
  
}




#------------------------------------------------
#' Simulation_R function
#'
#' This function triggers the main MAGENTA simulation from the R side
#'
#' @param paramList Paramlist passed from \code{Param_List_Simulation_Create}
#'
#' @export

# The following commands are needed to ensure that the roxygen2
# package, which deals with documenting the package, does not conflict
# with the Rcpp package. Do not alter!
#' @useDynLib MAGENTA
#' @importFrom Rcpp evalCpp
#' @exportPattern '^[[:alpha:]]+'

Simulation_R <- function(paramList)
{
  
  # check that this function is working
  print("R function is working!")
  
  ## Check if paramlist is correct length and has right variable names
  stopifnot(is.list(paramList))
  if(length(paramList)==3)
  {
    stopifnot(identical(names(paramList), c("N", "years", "eqSS")))  
  }
  # if it is length one it may be an unpacked list in which case unpack and check
  # this might happen in the future when a list of paramLists is fed directly to this
  # fucntion in a cluster way
  else if(length(paramList)==1)
  {
    paramList <- paramList[[1]]
    stopifnot(identical(names(paramList), c("N", "years", "eqSS")))  
  } 
  else 
  {
    stop("paramList not correct length")
  }
  
  # ---------------------- RUN C CODE
  
  # call Rcpp command with input list
  rawOutput <- Simulation_cpp(paramList)
  
  # ----------------------
  
  # convert rawOutput to final output format
  rawOutput
}