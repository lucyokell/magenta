#------------------------------------------------
#' Parameter List creation for MAGENTA simulation initialisation
#'
#' \code{Param_List_Simulation_Init_Create} creates suitable parameter list for
#' \code{Simulation_R} for the beginning of a simulation. Also takes an argument
#' for feeding in spatial parameters/data.
#'
#' @param N Population size. Default = 1e4
#' @param eqSS Output of \code{Equilibrium_Steady_State_Create}
#' @param barcode_parms List of barcode/genetic parameters
#' 
#' @export

Param_List_Simulation_Init_Create <- function(N = 1e+04, eqSS, barcode_parms)
{
  
  ## CHECKS ##
  ##---------------------------------------------
  if(class(eqSS)!="list") stop("eqSS is not of class list")
  
  if(is.element("spatial_type",names(eqSS))){
    if(!identical(names(eqSS),
                  c("age_brackets","het_brackets","Smat",
                    "Dmat","Amat","Umat","Tmat","Pmat",
                    "IBmat","ICAmat","ICMmat","IDmat",
                    "Sv","Ev","Iv","MaternalImmunity",
                    "theta","spatial_type"))) stop("Incorrect variable names within equilibrium.steady.state")
  } else {
    if(!identical(names(eqSS),
                  c("age_brackets","het_brackets","Smat",
                    "Dmat","Amat","Umat","Tmat","Pmat",
                    "IBmat","ICAmat","ICMmat","IDmat",
                    "Sv","Ev","Iv","MaternalImmunity",
                    "theta"))) stop("Incorrect variable names within equilibrium.steady.state")
  }
  dims.1 <- lapply(eqSS,function(x){return(dim(x)[1])})
  dims.2 <- lapply(eqSS,function(x){return(dim(x)[2])})
  if(unique(unlist(dims.1[grep("mat",names(eqSS))]))!=length(eqSS$age_brackets)) stop("Dimensions 1 error in equilibrium.stead.state matrices")
  if(unique(unlist(dims.2[grep("mat",names(eqSS))]))!=length(eqSS$het_brackets)) stop("Dimensions 2 error in equilibrium.stead.state matrices")
  
  numerics <- which(!unlist(lapply(eqSS,is.numeric)))
  if (length(numerics)!=0) stop(paste(names(numerics),"provided not numeric"))
  
  ##---------------------------------------------
  
  # Create paramlist
  paramList <- list(N = N, eqSS = eqSS, barcode_parms = barcode_parms)
  
  return(paramList)
  
}

#------------------------------------------------
#' Parameter List creation for MAGENTA simulation updating
#'
#' \code{Param_List_Simulation_Update_Create} creates suitable parameter list for
#' \code{Simulation_R} for continuing a simulation from memory within the active
#' session.
#'
#' @param years Length of simulation. Default = 1
#' @param ft Treatments seeking value
#' @param mu_vec Vector of mosquito mortalities for each day with years. Default = NULL, which
#' will result in rep(0.132,floor(years*365))
#' @param fv_vec Vector of mosquito bitings for each day with years. Default = NULL, which
#' will result in rep(1/3,floor(years*365))
#' @param statePtr Pointer for current model state as return by \code{Simulation_R}$Ptr
#' 
#' @export

Param_List_Simulation_Update_Create <- function(years = 1, ft = 0.4,
                                                mu_vec = NULL, fv_vec = NULL,
                                                statePtr)
{
  
  ## CHECKS ##
  ##---------------------------------------------
  if(class(statePtr)!="externalptr") stop("state.ptr is not of class externalptr")
  if(!is.numeric(years)) stop("years provided is not numeric")
  if(!is.numeric(ft)) stop("ft provided is not numeric")
  if(!is.null(mu_vec)){
    if(!(length(mu_vec) == floor(years*365))) stop("mu_vec not long enough")
  }
  if(!is.null(fv_vec)){
    if(!(length(fv_vec) == floor(years*365))) stop("fv_vec not long enough")
  }
  
  ##---------------------------------------------
  
  # generate mu_vec if needed
  if(is.null(mu_vec)){
    mu_vec <- rep(0.132,floor(years*365))
  }
  # generate fv_vec if needed
  if(is.null(fv_vec)){
    fv_vec <- rep(1/3,floor(years*365))
  }
  
  # Create paramlist
  paramList <- list(years = years, ft = ft, mu_vec = mu_vec,
                    fv_vec = fv_vec,
                    statePtr = statePtr)
  
  return(paramList)
  
}

#------------------------------------------------
#' Parameter List creation for MAGENTA simulation getting (saving to disk)
#'
#' \code{Param_List_Simulation_Get_Create} creates suitable parameter list for
#' \code{Simulation_R} for continuing a simulation from memory within the active
#' session.
#'
#' @param statePtr Pointer for current model state as return by \code{Simulation_R}$Ptr
#' 
#' @export

Param_List_Simulation_Get_Create <- function(statePtr)
{
  
  ## CHECKS ##
  ##---------------------------------------------
  if(class(statePtr)!="externalptr") stop("state.ptr is not of class externalptr")
  
  ##---------------------------------------------
  
  # Create paramlist
  paramList <- list(statePtr = statePtr)
  
  return(paramList)
  
}

#------------------------------------------------
#' Parameter List creation for loading saved MAGENTA simulation
#'
#' \code{Param_List_Simulation_Saved_Init_Create} creates suitable parameter list for
#' \code{Simulation_R} for continuing a simulation from memory within the active
#' session.
#'
#' @param savedState Saved state generated by \code{Simulation_R} when provided with 
#' a \code{Param_List_Simulation_Get_Create} parameter list
#' 
#' @export

Param_List_Simulation_Saved_Init_Create <- function(savedState)
{
  
  ## CHECKS ##
  ##---------------------------------------------
  if(class(savedState)!="list") stop("savedState is not of class list")
  if(!identical(names(savedState),c("population_List","populations_event_and_strains_List","scourge_List","parameters_List"))) stop("Incorrect variable names within savedState")
 
  ## TODO: More checks here
  ##---------------------------------------------
  
  # Create paramlist
  paramList <- list(savedState = savedState)
  
  return(paramList)
  
}


#------------------------------------------------
#' Simulation_R function
#'
#' This function triggers the main MAGENTA simulation from the R side
#'
#' @param paramList Paramlist passed from \code{Param_List_Simulation_Init_Create}
#' or from \code{Param_List_Simulation_Update_Create}
#' @param seed Seed for the simulation
#' @export

# The following commands are needed to ensure that the roxygen2
# package, which deals with documenting the package, does not conflict
# with the Rcpp package. Do not alter!
#' @useDynLib MAGENTA
#' @importFrom Rcpp evalCpp

Simulation_R <- function(paramList, seed)
{
  
  # check that this function is working
  print("R function is working!")
  stopifnot(is.list(paramList))
  
  ## Decide whether the paramList is from initialisation, memory-continutation or continuation
  
  ## -----------------------------------
  ## 1. From initialisation
  ## -----------------------------------
  if(!is.null(paramList$eqSS)){
    
    ## Check if paramlist is correct length and has right variable names
    stopifnot(is.list(paramList))
    if(length(paramList)==3)
    {
      stopifnot(identical(names(paramList), c("N","eqSS","barcode_parms")))  
    }
    # if it is length one it may be an unpacked list in which case unpack and check
    # this might happen in the future when a list of paramLists is fed directly to this
    # fucntion in a cluster way
    else if(length(paramList)==1)
    {
      paramList <- paramList[[1]]
      stopifnot(identical(names(paramList), c("N", "eqSS", "barcode_parms")))  
    } 
    else 
    {
      stop("paramList not correct length")
    }
    
    # ---------------------- RUN C CODE ------------------------------------- #
    
    # call Rcpp command with input list
    set.seed(seed)
    message(runif(1,0,1))
    rawOutput <- Simulation_Init_cpp(paramList)
    
    # ----------------------------------------------------------------------- #
    
    
  }
  
  ## -----------------------------------
  ## 2. From memory-continutation
  ## -----------------------------------
  if(!is.null(paramList$years) & !is.null(paramList$statePtr)){
    
    ## Check if paramlist is correct length and has right variable names
    stopifnot(is.list(paramList))
    if(length(paramList)==5)
    {
      stopifnot(identical(names(paramList), c("years","ft","mu_vec","fv_vec","statePtr")))  
    }
    else 
    {
      stop("paramList not correct length")
    }
    
    # ---------------------- RUN C CODE ------------------------------------- #
    
    # call Rcpp command with input list
    set.seed(seed)
    rawOutput <- Simulation_Update_cpp(paramList)
    
    # ----------------------------------------------------------------------- #
    
  }
  
  ## -----------------------------------
  ## 3. From memory-continutation to saving
  ## -----------------------------------
  if(is.null(paramList$years) & !is.null(paramList$statePtr)){
    
    ## Check if paramlist is correct length and has right variable names
    stopifnot(is.list(paramList))
    if(length(paramList)==1)
    {
      stopifnot(identical(names(paramList), c("statePtr")))  
    }
    else 
    {
      stop("paramList not correct length")
    } 
    
    # ---------------------- RUN C CODE ------------------------------------- #
    
    # call Rcpp command with input list
    set.seed(seed)
    rawOutput <- Simulation_Get_cpp(paramList)
    
    # ----------------------------------------------------------------------- #
    
  }
  
  ## -----------------------------------
  ## 4. From saved state
  ## -----------------------------------
  if(!is.null(paramList$savedState)){
    
    ## Check if paramlist is correct length and has right variable names
    stopifnot(is.list(paramList))
    if(length(paramList)==1)
    {
      stopifnot(identical(names(paramList), c("savedState")))  
    }
    else 
    {
      stop("paramList not correct length")
    } 
    
    # ---------------------- RUN C CODE ------------------------------------- #
    
    # call Rcpp command with input list
    set.seed(seed)
    rawOutput <- Simulation_Saved_Init_cpp(paramList)
    
    # ----------------------------------------------------------------------- #
    
  }
  
  ## Catch if parameter list not matched well
  if(is.null(rawOutput)) stop("No matching parameter list handler found within Simulation_R")
  
  # return rawOutput
  return(rawOutput)
}