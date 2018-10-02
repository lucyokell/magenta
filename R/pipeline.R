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
#' @param spatial_incidence_matrix Spatial incidence for humans, i.e. importation vector
#' @param spatial_mosquitoFOI_matrix Spatial mosquio FOI, i.e. importation to mosquitoes vector
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
#' @param num_het_brackets Number of heterogeinity brackets to use in initialisation. Default = 5
#' @param num_age_brackets Number of age brackets to use in initialisation. Default = 200
#' @param geometric_age_brackets Boolean detailing whether age brackets are geometric. Default = TRUE
#' @param max_age Maximum age in age brackets. Default = 100 
#' @param use_odin Boolean detailing whether the intiial solution is run within the odin model 
#' first. Default = FALSE while the model is still buggy. 
#' @param full_save Boolean detailing whether the entire simulation is saved. Default = FALSE
#' @param human_only_full_save Boolean detailing whether just the human component of the simulation is saved within full_save. Default = FALSE
#' @param update_save Boolean detailing whether the logging output is saved each update_length up to years. Default = FALSE
#' @param human_update_save Boolean detailing if the human state is also saved during each update_length. Default = FALSE
#' @param summary_saves_only Boolean if summary tables about COI are saved within human yearly save only. Dataframes of age, clinical status
#' binned COI.
#' @param saved_state_path Full file path to a saved model state to be loaded and continued. Default = NULL,
#' which will trigger initialisation
#' @param seed Random seed. Default is Random
#' 
#' \code{Pipeline}
#' 
#' @export


Pipeline <- function(EIR=120, ft = 0.4, itn_cov = 0, irs_cov = 0,
                     survival_percentage = 0.18, oocyst_mean = 5, oocyst_shape = 5,
                     N=100000, years = 20,update_length = 365, 
                     country = NULL, admin = NULL, use_historic_interventions = FALSE,
                     spatial_type = NULL, redis_host = "fi--dideclusthn.dide.ic.ac.uk",spatial_uuid = NULL,
                     spatial_incidence_matrix = NULL,spatial_mosquitoFOI_matrix = NULL,
                     num_loci = 24,ibd_length = 1, plaf = rep(0.5,num_loci),prob_crossover = rep(0.5, num_loci),starting_ibd=0.0,
                     num_het_brackets = 5, num_age_brackets = 20, 
                     geometric_age_brackets = TRUE, max_age = 100, use_odin = FALSE, mu_vec=NULL, fv_vec=NULL,
                     full_save = FALSE, human_only_full_save = FALSE, human_only_full_summary_save = FALSE,
                     update_save = FALSE, human_update_save = FALSE, genetics_df_without_summarising = FALSE,
                     summary_saves_only = FALSE, set_up_only = FALSE, mean_only = TRUE,
                     saved_state_path = NULL,seed=as.integer(runif(1,1,1000000000)),
                     sample_size = Inf, sample_states = 0:5, sample_reps = 1,
                     housekeeping_list = housekeeping_list_create(),
                     drug_list = drug_list_create(),only_allele_freqs = TRUE,
                     nmf_list = nmf_list_create()){
  
  
  ## if no seed is specified then save the seed
  set.seed(seed)
  message(paste0("Seed set to ",seed))
  message("New Chapter Come On Maybs Today");
  
  # convert allele frequencies
  pop_alf <- function(nums,nl=num_loci){ 
    if(class(nums %>% unlist)=="raw") {colMeans(matrix(as.numeric(do.call(rbind,nums)),ncol=nl))} else {rep(NA,nl)}
  }
  
  # convert to lineages frequencies
  lineages <- function(nums,nl=num_loci){ 
    if(class(nums %>% unlist)=="raw") {table(factor(apply(matrix(as.numeric(do.call(rbind,nums)),ncol=nl),1,bitsToInt),levels=0:((2^nl)-1)))} else {rep(NA,2^nl)}
  }
  
  # ---
  ## If we don't have a saved state then we initialise first
  # ---
  
  if(is.null(saved_state_path)) 
  {
    
    ## Create parameter list, changine any key parameters, e.g. the average age
    mpl <- Model_Param_List_Create(eta = 1/(21*365))
    
    ## Create age brackets, either geometric or evenly spaced
    age_vector <- age_brackets(max_age, num_age_brackets, geometric_age_brackets)
    
    ## Create a near equilibirum initial condition
    eqInit <- Equilibrium_Init_Create(age.vector = age_vector,
                                      het.brackets = num_het_brackets,
                                      ft = ft[1],
                                      EIR = EIR,
                                      country = country,
                                      admin = admin,
                                      model.param.list = mpl)
    
    # reset seed here as there is some randomness in equlibirum (Need?)
    set.seed(seed)
    
    ## Next create the near equilibrium steady state
    eqSS <- Equilibrium_SS_Create(eqInit = eqInit, end.year=5, use_odin = use_odin)
    
    ## and the barcode parms list
    barcode_parms <- barcode_parms_create(num_loci = num_loci,
                                          ibd_length = ibd_length,
                                          plaf = plaf,
                                          prob_crossover = prob_crossover,
                                          starting_ibd = starting_ibd)
    
    ## pass to spatial
    # spatial checks
    if(!is.null(spatial_type)) {
      if(spatial_type == "metapop"){
        spatial_type <- 2
        if(is.null(spatial_uuid)) stop("Spatial_uuid is required if running spatial simulations")
        redis_id <- paste0("oj_",spatial_uuid)
      } else if(spatial_type == "island") {
        spatial_type <- 1
      }
    } else {
      spatial_type <- 0
    }
    
    # make spatial list
    spatial_incidence_matrix <- spl_matrix_check(spatial_incidence_matrix, years)
    spatial_mosquitoFOI_matrix <- spl_matrix_check(spatial_mosquitoFOI_matrix, years)
    spatial_list <- spl_create(spatial_type = spatial_type,
                               human_importation_rate_vector = spatial_incidence_matrix[1,],
                               mosquito_imporation_rate_vector = spatial_mosquitoFOI_matrix[1,],
                               cotransmission_freq_vector = ztrgeomintp(10000, 10, survival_percentage),
                               oocyst_freq_vector = ztrnbinom(10000, mean=oocyst_mean, size=oocyst_shape))
    
    ## Now check and create the parameter list for use in the Rcpp simulation
    pl <- Param_List_Simulation_Init_Create(N=N,eqSS=eqSS,
                                            barcode_parms = barcode_parms,
                                            spatial_list = spatial_list,
                                            housekeeping_list = housekeeping_list,
                                            drug_list = drug_list,
                                            nmf_list = nmf_list)
    
  } 
  else 
  {
    
    # If we have provided the saved state then load this and then delete as can be large
    saved_state <- readRDS(saved_state_path)
    pl <- Param_List_Simulation_Saved_Init_Create(savedState = saved_state)
    rm(saved_state)
    gc()
  }
  
  ## Create model simulation state
  sim.out <- Simulation_R(paramList = pl, seed = seed)
  
  if(set_up_only){
    message("Set Up Grab")
    ## Now let's save the simulation in full
    pl2 <- Param_List_Simulation_Get_Create(statePtr = sim.out$Ptr)
    sim.save <- Simulation_R(pl2, seed = seed)
    
    ## If we want just the humans then get the keybits and save that instead
    if(full_save){
      return(sim.save)
    }
    if(human_only_full_save)
    {
      Strains <- sim.save$populations_event_and_strains_List[c("Strain_infection_state_vectors", "Strain_day_of_infection_state_change_vectors","Strain_barcode_vectors" )]
      Humans <- c(sim.save$population_List[c("Infection_States", "Zetas", "Ages")],Strains)
      res <- Humans
      return(res)
    } 
    else
    {
      res <- sim.save
      return(res)
    }
  }
  
  ## if it's spatial set up the necessary redis lists
  # set this up here anyway and then less if loops later
  imported_barcodes <- NULL
  if(spatial_type==2){
    
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
  
  out <- mu_fv_create(eqInit = eqInit, ft = ft, itn_cov = itn_cov, irs_cov = irs_cov, years = years)
  if(length(ft) == 1) ft <- rep(ft, years)
  
  # If we have specified a yearly save we iterate through the total time in chunks saving the loggers at each stage
  if(update_save)
  {
    
    # set up results list 
    res <- list()
    length(res) <- round((years*365)/update_length)
    times <- rep(0,length(res)-1)
    
    # and set up annual checks for variables that change discretely
    year <- 1
    next_drug_cycle <- drug_list$g_temporal_cycling
    ft_now <- ft[year]
    
    # messaging
    message("Starting Stochastic Simulation for ", years, " years")
    p <- progress::progress_bar$new( format = paste0(" Running: [:bar] :percent eta: :eta"),
                                     total = length(res))
    p_print <- progress_logging(housekeeping_list, res, p, initial = TRUE)

    ## START MAIN SIMULATION LOOP
    if(length(res) > 1){
      for(i in 1:(length(res)-1)){
        
        # messaging
        p_print <- progress_logging(housekeeping_list, res, p, i, 
                                    initial = FALSE, p_print = p_print)

        times[i] <- Sys.time()
        
        ## annual updates
        if ((floor((((update_length * (i-1))+1)/365))) == year) {
          year <- year + 1
          ft_now <- ft[year]
          spatial_list <- spl_create(spatial_type = spatial_type,
                                     human_importation_rate_vector = spatial_incidence_matrix[year,],
                                     mosquito_imporation_rate_vector = spatial_mosquitoFOI_matrix[year,],
                                     cotransmission_freq_vector = sample(2,10000,replace = TRUE, prob = c(0.82,0.18)),
                                     oocyst_freq_vector = sample(5,10000,replace = TRUE, prob = c(0.5,0.3,0.1,0.075,0.025)))
        }

        
        # prepare simulation for update
        pl2 <- Param_List_Simulation_Update_Create(years = update_length/365, 
                                                   ft = ft_now,
                                                   mu_vec = out$mu[1:update_length + ((i-1)*update_length)],
                                                   fv_vec = out$fv[1:update_length + ((i-1)*update_length)], 
                                                   spatial_list = spatial_list,
                                                   drug_list = drug_list,
                                                   statePtr = sim.out$Ptr)
        
        # carry out simulation
        sim.out <- Simulation_R(pl2, seed = seed)
        
        # what are saving, does it include the humans
        if(human_update_save) 
        {
          # do we just want the summary data frame 
          if(summary_saves_only){
            
            if(length(sample_size)>1){
            df <- pop_strains_df(sim.out$Ptr, sample_size = 0, 
                                 sample_states = sample_states, ibd = barcode_parms$barcode_type,
                                 seed = seed,
                                 nl = num_loci)
            } else {
              df <- pop_strains_df(sim.out$Ptr, sample_size = sample_size*sample_reps, 
                                   sample_states = sample_states, ibd = barcode_parms$barcode_type,
                                   seed = seed, 
                                   nl = num_loci)
              
            }
            
            if(genetics_df_without_summarising) {
              res[[i]] <- list()
              if(only_allele_freqs){
                res[[i]]$af <- pop_alf(df$nums[unlist(lapply(df$barcode_states,length))>0]) %>% unlist
                res[[i]]$lineage <- lineages(df$nums[unlist(lapply(df$barcode_states,length))>0]) %>% unlist
                res[[i]]$pcr_prev <- sum(c(lapply(df$barcode_states,length) %>% unlist)>0)
              } else {
                res[[i]]$df <- df
              }
              res[[i]]$succesfull_treatments <- sim.out$Loggers$Treatments$Successful_Treatments
              res[[i]]$unsuccesful_treatments_lpf <- sim.out$Loggers$Treatments$Unsuccesful_Treatments_LPF
              res[[i]]$not_treated <- sim.out$Loggers$Treatments$Not_Treated
              res[[i]]$treatment_failure <-  res[[i]]$unsuccesful_treatments_lpf / (res[[i]]$unsuccesful_treatments_lpf + res[[i]]$succesfull_treatments) 
              
            } else {
            
            if(i%%12 == 0 && i >= (length(res)-180)){
              res[[i]] <- COI_df_create2(df, barcodes=TRUE, nl=num_loci, ibd = barcode_parms$barcode_type,
                                        n = sample_size, reps = sample_reps, mean_only = mean_only)
            } else {
              res[[i]] <- COI_df_create2(df, barcodes=FALSE, nl=num_loci, ibd = barcode_parms$barcode_type,
                                         n = sample_size, reps = sample_reps, mean_only = mean_only)
            }
            res[[i]]$Prev <- sum(unlist(sim.out$Loggers[c("D","A","U","T")]))
            }
            # or do we want the full human popualation
          } else {
            
            # then grab the population
            pl3 <- Param_List_Simulation_Get_Create(statePtr = sim.out$Ptr)
            sim.save <- Simulation_R(pl3, seed = seed)
            
            # and then store what's needed
            Strains <- sim.save$populations_event_and_strains_List[c("Strain_infection_state_vectors", "Strain_day_of_infection_state_change_vectors","Strain_barcode_vectors" )]
            Humans <- c(sim.save$population_List[c("Infection_States", "Zetas", "Ages")],Strains)
            res[[i]] <- Humans
          }
          
          # or are we just saving the loggers
        } 
        else 
        {
          res[[i]] <- sim.out$Loggers
        }
        
        ## spatial export
        if(spatial_type==2) {
          
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
        
        ## drug resistance updates
        drug_list <- drug_list_update(drug_list, year, res[[i]]$treatment_failure)
        
      }
    }
    ## final run
    pl2 <- Param_List_Simulation_Update_Create(years = update_length/365, ft = ft_now,
                                               mu_vec = out$mu[1:update_length + ((i)*update_length)],
                                               fv_vec = out$fv[1:update_length + ((i)*update_length)],
                                               spatial_list = spatial_list,
                                               drug_list = drug_list,
                                               statePtr = sim.out$Ptr)
    sim.out <- Simulation_R(pl2, seed = seed)
    
    ## If we have specified a full save then we grab that and save it or just the human bits of interest
    if(full_save || human_only_full_save || human_only_full_summary_save)
    {
      
      ## Now let's save the simulation in full
      pl2 <- Param_List_Simulation_Get_Create(statePtr = sim.out$Ptr)
      sim.save <- Simulation_R(pl2, seed = seed)
      
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
    pl2 <- Param_List_Simulation_Update_Create(years = years, ft = min(ft[ft>0]),
                                               mu_vec = out$mu,
                                               fv_vec = out$fv,
                                               spatial_list = spatial_list,
                                               drug_list = drug_list,
                                               statePtr = sim.out$Ptr)
    
    
    
    ## Now run the simulation
    sim.out <- Simulation_R(paramList = pl2, seed = seed)
    
    ## If we have specified a full save or human save then we grab that and save it or just the human bits of interest
    if(full_save || human_only_full_save)
    {
      #  browser()
      ## Now let's save the simulation in full
      pl2 <- Param_List_Simulation_Get_Create(statePtr = sim.out$Ptr)
      sim.save <- Simulation_R(pl2, seed = seed)
      
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
  
  # append times
  if(update_save){
  attr(res,"times") <- times
  }
  
  # Save the seed as an attribute adn return the result
  seed_end <- .Random.seed
  attr(res,"seed") <- seed_end
  return(res)
  
}
