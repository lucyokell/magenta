#------------------------------------------------
#' Pipeline for cluster submission
#'
#' \code{Pipeline} steps through creating the parameter list, the equilibrium
#' initialisation and steady state creation before checking and passing suitable
#' parameters to the simulation. This is then saved 
#' 
#' @param N Population Size. Default = 100000
#' @param years Lenth of simulation. Default = 20
#' @param EIR Numeric for desired annual EIR. Default = 120
#' @param num_het_brackets Number of heterogeinity brackets to use in initialisation. Default = 200
#' @param num_age_brackets Number of age brackets to use in initialisation. Default = 200
#' @param geometric_age_brackets Boolean detailing whether age brackets are geometric. Default = TRUE
#' @param max_age Maximum age in age brackets. Default = 100 
#' @param use_odin Boolean detailing whether the intiial solution is run within the odin model 
#' first. Default = FALSE while the model is still buggy. 
#' @param full_save Boolean detailing whether the entire simulation is saved. Default = FALSE
#' 
#' \code{Pipeline}
#' 
#' @export


Pipeline <- function(EIR=120, N=100000, years = 20, 
                     num_het_brackets = 200, num_age_brackets = 200, geometric_age_brackets = TRUE, max_age = 100, use_odin = FALSE, 
                     full_save = FALSE){
  
  ## Pipeline
  Sys.setenv(BINPREF="T:/Rtools/Rtools33/mingw_64/bin/")
  
  ## Create parameter list, changine any key parameters, e.g. the average age
  mpl <- Model_Param_List_Create(eta = 1/(21*365))
  
  ## Create age brackets, either geometric or evenly spaced
  if(geometric_age_brackets){
    ## Create the geometric age brackets
    ratio <- (max_age/0.1)^(1/num_age_brackets)
    age.vector <- 0.1 * ratio ** (1:num_age_brackets)
    age.vector[1] <- 0
  } else {
    age.vector <- seq(0,max_age,num_age_brackets)
  }
  
  ## Create a near equilibirum initial condition
  eqInit <- Equilibrium_Init_Create(age.vector = age.vector,
                                    het.brackets = num_het_brackets,
                                    ft = 0.4,
                                    EIR = EIR,
                                    model.param.list = mpl)
  
  ## Next create the near equilibrium steady state
  eqSS <- Equilibrium_SS_Create(eqInit = eqInit, end.year=5, use_odin = use_odin)
  
  ## Now check and create the parameter list for use in the Rcpp simulation
  pl <- Param_List_Simulation_Init_Create(N=N,years=years,eqSS=eqSS)
  
  ## Now run the simulation
  sim.out <- Simulation_R(paramList = pl)
  res <- sim.out
  
  if(full_save){
    ## Now let's save the simulation in full
    pl2 <- MAGENTA::Param_List_Simulation_Get_Create(statePtr = sim.out$Ptr)
    res <- MAGENTA::Simulation_R(pl2)
  }
  
  return(res)
  
}