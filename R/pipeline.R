#' Pipeline for cluster submission
#'
#' \code{Pipeline} steps through creating the parameter list, the equilibrium
#' initialisation and steady state creation before checking and passing suitable
#' parameters to the simulation. This is then saved. If a path to a savedState is
#' provided then this state is loaded and continued. 
#' 
#' @param N Population Size. Default = 100000
#' @param ft Frequency of treatment. Default = 0.4
#' @param years Lenth of simulation. Default = 20
#' @param update_length How long each update is run for in days. Default = 365
#' @param EIR Numeric for desired annual EIR. Default = 120
#' @param country Character for country within which admin2 is in. Default = NULL
#' @param admin Character for admin region. Some fuzzy logic will be used to match. If 
#' not provided then no seasonality is introduced. Default = NULL
#' @param spatial_type Default = NULL. If spatial is wanted then provide a character describing 
#' the type of spatial simulation, which must be one of "island" or "metapop".
#' @param redis_host Default = "fi--dideclusthn.dide.ic.ac.uk". This should be the host address
#' where your redis server is running. 
#' @param spatial_uuid Default = NULL. If spatial is provided, then this will error unless a character
#' string is passed to this argument. This character should be the same for each parallel task within
#' the same job. 
#' @param fixed_spatial_incidence_matrix Alternative parameterisation where hypothetical spatial matrix is supplied
#' @param fixed_spatial_mosquitoFOI_matrix Alternative parameterisation where hypothetical spatial matrix is supplied
#' @param fv_vec Numeric for how fv_vec changes as calculated from odin model
#' @param mu_vec Numeric for hor mu_vec changes as calculated from odin model
#' @param use_historic_interventions Boolean as to whether the historic interventions are incorporated.
#'   This will occur by using the length of years and taking the most recent years back in time. Therefore
#'   if years > 15, then there will be burn in time, when no interventions are assumed. Default = FALSE, and
#'   only is used if country and admin are suitably provided.  
#' ## Genetic Params
#' @param num_loci Number of loci. Default = 24
#' @param ibd_length If we are simulating IBD dynamics, each loci is now 
#'   represented by a bitset of ibd_length. Thus ibd_length needs to be long
#'   enough to ensure that as new identity relationships occur, i.e. 
#'   an importation barcode will be a new identity. e.g. If are population is
#'   1000, we may expect at 80% prevalence, with a mean COI of 3 we will need 
#'   2400 different identities, i.e. 2^ibd_length > 2400. However, keep in mind
#'   importations as these need to be continually new, i.e. if we are simulating
#'   for 30 years, with 3 importations a day, then we will need at least length
#'   to ensure that 2^ibd_length > 2400 + (30*365*3). This will probably be
#'   automatically calculated in the future. If we are not ding IBD, then this
#'   should be 1, which is the defalt. 
#' @param plaf Vector of population level allele frequencies for the barcode. 
#'   Default = rep(0.5, num_loci)
#' @param prob_crossover Vector of probabilities for crossover events for the 
#'   barcode. Default = rep(0.5, num_loci)
#' ## odin Params
#' @param num_het_brackets Number of heterogeinity brackets to use in initialisation. Default = 200
#' @param num_age_brackets Number of age brackets to use in initialisation. Default = 200
#' @param geometric_age_brackets Boolean detailing whether age brackets are geometric. Default = TRUE
#' @param max_age Maximum age in age brackets. Default = 100 
#' @param use_odin Boolean detailing whether the intiial solution is run within the odin model 
#' first. Default = FALSE while the model is still buggy. 
#' @param full_save Boolean detailing whether the entire simulation is saved. Default = FALSE
#' @param human_only_full_save Boolean detailing whether just the human component of the simulation is saved within full_save. Default = FALSE
#' @param yearly_save Boolean detailing whether the logging output is saved each update_length up to years. Default = FALSE
#' @param human_yearly_save Boolean detailing if the human state is also saved during each update_length. Default = FALSE
#' @param summary_saves_only Boolean if summary tables about COI are saved within human yearly save only. Dataframes of age, clinical status
#' binned COI.
#' @param saved_state_path Full file path to a saved model state to be loaded and continued. Default = NULL,
#' which will trigger initialisation
#' @param seed Random seed. Default is Random
#' 
#' \code{Pipeline}
#' 
#' @export


Pipeline <- function(EIR=120, ft = 0.4, N=100000, years = 20,update_length = 365, 
                     country = NULL, admin = NULL, use_historic_interventions = FALSE,
                     spatial_type = NULL, redis_host = "fi--dideclusthn.dide.ic.ac.uk",spatial_uuid = NULL,
                     fixed_spatial_incidence_matrix = NULL,fixed_spatial_mosquitoFOI_matrix = NULL,
                     num_loci = 24,ibd_length = 1, plaf = rep(0.5,num_loci),prob_crossover = rep(0.5, num_loci),
                     num_het_brackets = 20, num_age_brackets = 20, 
                     geometric_age_brackets = TRUE, max_age = 100, use_odin = FALSE, mu_vec=NULL, fv_vec=NULL,
                     full_save = FALSE, human_only_full_save = FALSE, yearly_save = FALSE, human_yearly_save = FALSE,
                     summary_saves_only = FALSE,
                     saved_state_path = NULL,seed=runif(1,1,10000)){
  
  ## Pipeline
  ## Sys.setenv(BINPREF="T:/Rtools/Rtools33/mingw_64/bin/")
  
  ## if no seed is specified then save the seed
  set.seed(seed)
  message(paste0("Seed set to ",seed))
  
  # ---
  ## If we don't have a saved state then we initialise first
  # ---
  
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
                                      ft = ft,
                                      EIR = EIR,
                                      country = country,
                                      admin = admin,
                                      model.param.list = mpl)
    
    ## Next create the near equilibrium steady state
    eqSS <- Equilibrium_SS_Create(eqInit = eqInit, end.year=5, use_odin = use_odin)
    
    ## and the barcode parms list
    barcode_parms <- barcode_parms_create(num_loci = num_loci,
                                          ibd_length = ibd_length,
                                          plaf = plaf,
                                          prob_crossover = prob_crossover)
    
    ## pass to spatial
    
    # spatial checks
    if(!is.null(spatial_type)) {
      if(spatial_type == "metapop"){
        if(is.null(spatial_uuid)) stop("Spatial_uuid is required if running spatial simulations")
        redis_id <- paste0("oj_",spatial_uuid)
        eqSS$spatial_type <- 2
      } else if(spatial_type == "island") {
        eqSS$spatial_type <- 1
      }
    } else {
      eqSS$spatial_type <- 0
    }
    
    ## Now check and create the parameter list for use in the Rcpp simulation
    pl <- Param_List_Simulation_Init_Create(N=N,eqSS=eqSS,
                                            barcode_parms = barcode_parms)
    
  }
  
  # ---
  ## If there is a saved state path then we load this
  # ---
  
  else 
  {
    # If we have provided the saved state then load this
    saved_state <- readRDS(saved_state_path)
    pl <- Param_List_Simulation_Saved_Init_Create(savedState = saved_state)
    
  }
  
  
  ## Create model simulation state
  sim.out <- Simulation_R(paramList = pl)
  
  ## if it's spatial set up the necessary redis lists
  # set this up here anyway and then less if loops later
  imported_barcodes <- NULL
  if(!is.null(spatial_type)){
    
    # create barcode as binary string
    barcodes <- lapply(sim.out$Exported_Barcodes,as.numeric) %>% lapply(paste0,collapse="") %>% unlist
    
    # work out here how many barcodes are going to which other simulation
    export_proportions <- round((spatial/(sum(spatial,na.rm=TRUE)))*round(sum(spatial *N,na.rm=TRUE)))
    export_positions <- list()
    import_positions <- list()
    check_barcodes <- list()
    epc <- 1
    metapopulation_number <- which(is.na(spatial))
    other_metapopulations <- (function(x) x[-metapopulation_number])(1:length(spatial))
    
    # Create the necessary redis lists and push the barcodes
    for(i in 1:length(export_proportions)){
      if(!is.na(export_proportions[i])){
        export_positions[[i]] <- epc:(epc+export_proportions[i]-1)
        import_positions[[i]] <- paste0(i,metapopulation_number)
        check_barcodes[[i]] <- paste0(i,"DONE")
        epc <- epc+export_proportions[i]
        redis_environment$SET(paste0(redis_id,"_",metapopulation_number,i),redux::object_to_bin(barcodes[export_positions[[i]]]))
      } else {
        export_positions[[i]] <- export_proportions[i]
        import_positions[[i]] <- NULL
        check_barcodes[[i]] <- NULL
      }
    }  
  
    ## clear the blank ones
    import_positions[[metapopulation_number]] <- NULL
    check_barcodes[[metapopulation_number]] <- NULL
    
    # Say that this metapopulation is finished
    redis_environment$SET(paste0(metapopulation_number,"DONE"),1)

    ## create the redis object for grabbing other barcodes
    redis <- redux::redis
    get_barcodes_cmds <- lapply(import_positions, function(x){ (redis$GET(x))})
    check_barcodes_cmds <- lapply(check_barcodes, function(x){ (redis$GET(x))})
    
    # sit on a while loop here untill all populations are ready
    cannot_import <- TRUE
    while(cannot_import){
      barcode_checks <- redis_environment$pipeline(.commands = check_barcodes_cmds) %>% unlist
      if(sum(is.null(barcode_checks))==0){
        cannot_import <- FALSE
      }
    }
    
    # once they are all ready set the done to 0 
    redis_environment$SET(paste0(redis_id,"_",metapopulation_number,"DONE"),0)
    
    # fetch and format the barcodes for import
    imported_barcodes <- lapply(redis_environment$pipeline(.commands = get_barcodes_cmds),redux::bin_to_object) %>% unlist
    imported_barcodes <- lapply(imported_barcodes,function(x){strsplit(x,"") %>% unlist %>% as.numeric %>% as.logical}) 
  }
  
  ## incorporate historic interventions as needed
  if(use_historic_interventions){

    ## THIS WILL COME BACK WHEN CAN GET ODIN WORKING AGAIN
    
    if(!is.null(admin)){
      eqInit$itn_cov <- rep(0,ceiling(years))
      eqInit$irs_cov <- rep(0,ceiling(years))
      eqInit$int_len <- length(eqInit$itn_cov)
      eqInit$ITN_IRS_on <- -1
    } else {

    # incorporate hitoric interventions
    eqInit$itn_cov <- itn_2010_2015$value[which(itn_2010_2015$admin==admin & itn_2010_2015$country==country)]
    eqInit$irs_cov <- irs_2010_2015$value[which(irs_2010_2015$admin==admin & irs_2010_2015$country==country)]

    if(length(eqInit$irs_cov)!=16){
      eqInit$irs_cov <- irs_2010_2015$value[which(irs_2010_2015$country==country)]
    }

    if(length(eqInit$irs_cov)!=16){
      if(length(eqInit$itn_cov==1)){
        eqInit$irs_cov <- 0
      } else {
      eqInit$irs_cov <- rep(0,16)
      }
    }
    eqInit$int_len <- length(eqInit$itn_cov)
    eqInit$ITN_IRS_on <- (years - 16)*365
    }

    # # create model
    # user_vars <- names(formals(odin_model)$user)[-match(c("","age"),names(formals(odin_model)$user))]
    #
    # # create temp tp store generator
    # temp <- tempfile(fileext = ".txt")
    #
    # # create the needed funciton call and run
    # writeLines(text = paste0("model <- odin_model(",paste0(sapply(user_vars, function(x){(paste0(x,"=eqInit$",x,",\n") )}),collapse = ""),
    #                          "age = eqInit$age*365,\n",
    #                          "use_dde=TRUE)"),con = temp)
    # source(temp,local = TRUE)

    odin_model_path <- system.file("extdata/odin_model.R",package="MAGENTA")
    gen <- odin::odin(odin_model_path,verbose=FALSE,build = TRUE)

    model <- generate_default_model(ft=eqInit$ft,age=eqInit$age_brackets,dat=eqInit,generator=gen,dde=TRUE)


    #create model and simulate
    tt <- seq(0,years*365,1)
    mod_run <- model$run(tt)
    out <- model$transform_variables(mod_run)
  
    out <- data.frame("mu"=mu_vec,"fv"=fv_vec)
  }
  
  # If we have specified a yearly save we iterate through the total time in chunks saving the loggers at each stage
  if(yearly_save)
  {
    
    res <- list()
    length(res) <- round((years*365)/update_length)
    
 
    
    # adapt intervetnion datframe to match
    if(use_historic_interventions){
      mu_needed <- tail(out$mu,length(res)*30)
      fv_needed <- tail(out$fv,length(res)*30)
      out <- data.frame("mu"=mu_needed,"fv"=fv_needed)
    }
    
    ## START MAIN SIMULATION LOOP
    if(length(res) > 1){
    for(i in 1:(length(res)-1)){
      
      ## if historical or not set up the parameter list accordingly
      if(use_historic_interventions){
      pl2 <- Param_List_Simulation_Update_Create(years = update_length/365, ft = ft,
                                                 mu_vec = out$mu[1:update_length + ((i-1)*update_length)],
                                                 fv_vec = out$fv[1:update_length + ((i-1)*update_length)], 
                                                 #imported_barcodes = imported_barcodes,
                                                 statePtr = sim.out$Ptr)
      } else {
        pl2 <- Param_List_Simulation_Update_Create(years = update_length/365, ft = ft,
                                                  # imported_barcodes = imported_barcodes,
                                                   statePtr = sim.out$Ptr)  
      }
      
      # carry out simulation
      sim.out <- Simulation_R(pl2)
      
      # what are saving, does it include the humans
      if(human_yearly_save) {
        pl3 <- Param_List_Simulation_Get_Create(statePtr = sim.out$Ptr)
        sim.save <- Simulation_R(pl3)
        
        # do we just want the summary data frame 
        if(summary_saves_only){
          if(i%%12 == 0 && i >= (length(res)-180)){
          res[[i]] <- COI_df_create(sim.save,barcodes=TRUE)
          } else {
          res[[i]] <- COI_df_create(sim.save)
          }
          
        # or do we want the full human popualation
        } else {
          Strains <- sim.save$populations_event_and_strains_List[c("Strain_infection_state_vectors", "Strain_day_of_infection_state_change_vectors","Strain_barcode_vectors" )]
          Humans <- c(sim.save$population_List[c("Infection_States", "Zetas", "Ages")],Strains)
          res[[i]] <- Humans
        }
        
      # or are we just saving the loggers
      } else {
        res[[i]] <- sim.out$Loggers
      }
     
      ## spatial export
      if(!is.null(spatial_type)){
        
        # Push barcodes to redis
        for(i in 1:length(export_proportions)){
          if(!is.na(export_proportions[i])){
            redis_environment$SET(paste0(redis_id,"_",metapopulation_number,i),redux::object_to_bin(barcodes[export_positions[[i]]]))
          }
        }  
        
        # once pushed set the done to 1
        redis_environment$SET(paste0(redis_id,"_",metapopulation_number,"DONE"),1)
        
        # sit on a while loop here untill all populations are ready
        cannot_import <- TRUE
        while(cannot_import){
          barcode_checks <- redis_environment$pipeline(.commands = check_barcodes_cmds) %>% unlist
          if(sum(barcode_checks==0)==0){
            cannot_import <- FALSE
          }
        }
        
        # once they are all ready set the done to 0 
        redis_environment$SET(paste0(redis_id,"_",metapopulation_number,"DONE"),0)
        
        # fetch and format the barcodes for import
        imported_barcodes <- lapply(redis_environment$pipeline(.commands = get_barcodes_cmds),redux::bin_to_object) %>% unlist
        imported_barcodes <- lapply(imported_barcodes,function(x){strsplit(x,"") %>% unlist %>% as.numeric %>% as.logical}) 
        
      }
       
    }
    }
    ## final run
    ## if historical or not
    if(use_historic_interventions){
      pl2 <- Param_List_Simulation_Update_Create(years = update_length/365, ft = ft,
                                                 mu_vec = out$mu[1:update_length + ((i)*update_length)],
                                                 fv_vec = out$fv[1:update_length + ((i)*update_length)],
                                                 #imported_barcodes = imported_barcodes,
                                                 statePtr = sim.out$Ptr)
    } else {
      pl2 <- Param_List_Simulation_Update_Create(years = update_length/365, ft = ft,
                                                 #imported_barcodes = imported_barcodes,
                                                 statePtr = sim.out$Ptr)  
    }

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
        res[[length(res)]] <- Humans
      } 
      else
      {
        res[[length(res)]] <- sim.save
      }
      
    } 
    else 
    {
      ## If we don't want a full save then just save the Loggers as usual
      res[[length(res)]] <- sim.out
      
    }
    
  } 
  else
  {
    
    ## Set up update for years long
    ## if historical or not
    if(use_historic_interventions){
      pl2 <- Param_List_Simulation_Update_Create(years = years, ft = ft,
                                                 mu_vec = out$mu,
                                                 fv_vec = out$fv,
                                                 statePtr = sim.out$Ptr)
    } else {
      pl2 <- Param_List_Simulation_Update_Create(years = years, ft = ft, statePtr = sim.out$Ptr)
    }

    
    ## Now run the simulation
    sim.out <- Simulation_R(paramList = pl2)
    
    ## If we have specified a full save or human save then we grab that and save it or just the human bits of interest
    if(full_save || human_only_full_save)
    {
    #  browser()
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
  
  # Save the seed as an attribute adn return the result
  attr(res,"seed") <- seed
  return(res)
  
}
