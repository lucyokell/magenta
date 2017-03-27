#------------------------------------------------
#' Pipeline for cluster submission
#'
#' \code{Pipeline} steps through creating the parameter list, the equilibrium
#' initialisation and steady state creation before checking and passing suitable
#' parameters to the simulation. This is then saved. If a path to a savedState is
#' provided then this state is loaded and continued. 
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
#' @param human_only_full_save Boolean detailing whether just the human component of the simulation is saved within full_save. Default = FALSE
#' @param yearly_save Boolean detailing whether the logging output is saved each year up to years. Default = FALSE
#' @param saved_state_path Full file path to a saved model state to be loaded and continued. Default = NULL,
#' which will trigger initialisation
#' 
#' \code{Pipeline}
#' 
#' @export


Pipeline <- function(EIR=120, N=100000, years = 20, 
                     num_het_brackets = 200, num_age_brackets = 200, geometric_age_brackets = TRUE, max_age = 100, use_odin = FALSE, 
                     full_save = FALSE, human_only_full_save = FALSE, yearly_save = FALSE,
                     saved_state_path = NULL){
  
  ## Pipeline
  Sys.setenv(BINPREF="T:/Rtools/Rtools33/mingw_64/bin/")
  
  ## If we don't have a saved state then we initialise first
  if(is.null(saved_state_path))
  {
    
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
  
  }
  ## If there is a saved state path then we load this
  else 
  {
    # If we have provided the saved state then load this
    saved_state <- readRDS(saved_state_path)
    pl <- Param_List_Simulation_Saved_Init_Create(years = years, savedState = saved_state)
    
  }
  
  # If we have specified a yearly save we iterate through the total time in year chunks saving the loggers at each stage
  ## TODO: Put the human state save in the year chuck section as well 
  if(yearly_save)
    {
      
      res <- list()
      length(res) <- years
      
      pl$years <- 1
      sim.out <- Simulation_R(paramList = pl)
      res[[1]] <- sim.out$Loggers
      
      for(i in 2:(years - 1)){
        pl2 <- Param_List_Simulation_Update_Create(years = 1, statePtr = sim.out$Ptr)
        sim.out <- Simulation_R(pl2)
        res[[i]] <- sim.out$Loggers
      }
      
      ## final run
      pl2 <- Param_List_Simulation_Update_Create(years = 1, statePtr = sim.out$Ptr)
      sim.out <- Simulation_R(pl2)
      
      ## If we have specified a full save then we grab that and save it or just the human bits of interest
      if(full_save || human_only_full_save)
      {
        
        ## Now let's save the simulation in full
        pl2 <- Param_List_Simulation_Get_Create(statePtr = sim.out$Ptr)
        sim.save <- Simulation_R(pl2)
        
        ## If we want just the humans then get the keybits and save that instead
        if(human_only_full_save)
        {
          Strains <- sim.save$populations_event_and_strains_List[c("Strain_infection_state_vectors", "Strain_day_of_infection_state_change_vectors","Strain_barcode_vectors" )]
          Humans <- c(sim.save$population_List[c("Infection_States", "Zetas", "Ages")],Strains)
          res[[years]] <- Humans
        } 
        else
        {
          res[[years]] <- sim.save
        }
        
      } 
      else 
      {
        ## If we don't want a full save then just save the Loggers as usual
        res[[years]] <- sim.out
        
      }
      
    } 
    else
    {
      
      ## Now run the simulation
      sim.out <- Simulation_R(paramList = pl)
      
      ## If we have specified a full save or human save then we grab that and save it or just the human bits of interest
      if(full_save || human_only_full_save)
      {
        
        ## Now let's save the simulation in full
        pl2 <- Param_List_Simulation_Get_Create(statePtr = sim.out$Ptr)
        sim.save <- Simulation_R(pl2)
        
        ## If we want just the humans then get the keybits and save that instead
        if(human_only_full_save)
        {
          Strains <- sim.save$populations_event_and_strains_List[c("Strain_infection_state_vectors", "Strain_day_of_infection_state_change_vectors","Strain_barcode_vectors" )]
          Humans <- c(sim.save$population_List[c("Infection_States", "Zetas", "Ages")],Strains)
          res <- Humans
        } 
        else
        {
          res <- sim.save
        }
        
      } 
      else 
      {
        ## If we don't want a full save then just save the Loggers as usual
        res <- sim.out
        
      }
      
    }
    
  # Save the result
  return(res)
  
}