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
#' 
#' \code{Pipeline}
#' 
#' @export


Pipeline <- function(EIR=120, N=100000, years = 20){
  
  ## Pipeline
  Sys.setenv(BINPREF="T:/Rtools/Rtools33/mingw_64/bin/")
  
  ## Create parameter list, changine any key parameters, e.g. the average age
  mpl <- Model_Param_List_Create(eta = 1/(21*365))
  
  ## Create a near equilibirum initial condition
  eqInit <- Equilibrium_Init_Create(age.vector = c(0,0.25,0.5,0.75,1,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80,90,100),
                                    het.brackets = 5,ft = 0.4,EIR = EIR,model.param.list = mpl)
  
  ## Next create the near equilibrium steady state
  eqSS <- Equilibrium_SS_Create(eqInit = eqInit, end.year=5)
  
  ## Now check and create the parameter list for use in the Rcpp simulation
  pl <- Param_List_Simulation_Init_Create(N=N,years=years,eqSS=eqSS)
  
  ## Now run the simulation
  sim.out <- Simulation_R(paramList = pl)
  
  ## Now let's save the simulation in full
  pl2 <- MAGENTA::Param_List_Simulation_Get_Create(statePtr = sim.out$Ptr)
  sim.save <- MAGENTA::Simulation_R(pl2)
  
  return(sim.save)
  
}