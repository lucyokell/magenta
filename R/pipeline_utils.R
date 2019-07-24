#' Creates age brackets
#' 
#' Makes age brackets given, number of brackets, max age, and spacing
#' 
#' @param max_age Maximum age in years
#' @param num_age_brackets Number of age brackets
#' @param geometric_age_brackets Boolean for geometric brackets

age_brackets <- function(max_age=100, 
                         num_age_brackets=20, 
                         geometric_age_brackets=TRUE){
  
  if(geometric_age_brackets){
    ## Create the geometric age brackets
    ratio <- (max_age/0.1)^(1/num_age_brackets) 
    age_vector <- 0.1 * ratio ** (1:num_age_brackets)
    age_vector[1] <- 0
  } else {
    age_vector <- seq(0,max_age,num_age_brackets)
  }
  return(age_vector)
}

#' Spatial matrix grab
#' 
#' Grabs spatial incidence and mosquitoFOI matrix from database
#' 
#' @param country Country string
#' @param admin Admin string
#' @param year_range Year range

spl_grab <- function(country, admin, year_range) {
  
  importations <- magenta::importations
  
  # grab importations
  imp_nms <- colnames(importations$`2000`$incidence)
  ads <- strsplit(imp_nms, "|", fixed = TRUE) %>% 
    lapply(function(x) x[2]) %>% 
    unlist
  cc <- strsplit(imp_nms, "|", fixed = TRUE) %>% 
    lapply(function(x) x[1]) %>% 
    unlist 
  imp_match <- which(cc == country & ads == admin)
  
  incidence <- lapply(importations,function(x) x$incidence[imp_match,]) %>% 
    unlist() %>% 
    matrix(ncol = length(imp_nms), byrow = TRUE)
  
  mosquitoFOI <- lapply(importations,function(x) x$mosquitoFOI[imp_match,]) %>% 
    unlist() %>% 
    matrix(ncol = length(imp_nms), byrow = TRUE)
  
  # grab year bits
  end_year <- tail(year_range, 1)
  start_year <- year_range[1]
  years <- 2000:2017
  
  # alter the end 
  if (end_year > 2017) {
    
    extra <- length(tail(years,1):end_year)-1
    incidence <- rbind(incidence, t(replicate(extra,incidence[18,])))
    mosquitoFOI <- rbind(mosquitoFOI, t(replicate(extra,mosquitoFOI[18,])))
    
  } else if (end_year <= 2017 && end_year > 2000) {
    
    incidence <- incidence[1:(end_year - 1999),]
    mosquitoFOI <- mosquitoFOI[1:(end_year - 1999),]
    
    
  } else {
    
    incidence <- NULL
    mosquitoFOI <- NULL
    
  }
  
  # alter the start
  if (start_year < 2000) {
    incidence <- rbind(matrix(0, ncol = 593, nrow = length(start_year:1999)),
                       incidence)
    mosquitoFOI <- rbind(matrix(0, ncol = 593, nrow = length(start_year:1999)),
                         mosquitoFOI)
  } else {
    incidence <- incidence[which(years == start_year):dim(incidence)[1],]
    mosquitoFOI <- mosquitoFOI[which(years == start_year):dim(mosquitoFOI)[1],]
  }
  
  
  return(list("incidence" = incidence,
              "mosquitoFOI" = mosquitoFOI))
}


#' Spatial matrix grab
#' 
#' Grabs spatial incidence and mosquitoFOI matrix from database
#' 
#' @param matrix Spatial matrix for 2000 to 2017
#' @param years Years desired

spl_matrix_check <- function(matrix, years) {
  
  if(is.matrix(matrix)){
    if(dim(matrix)[1] != ceiling(years)){
      matrix <- rbind(matrix(rep(matrix[1,],ceiling(years) - dim(matrix)[1]),ncol=dim(matrix)[2],byrow=TRUE),matrix)
    }
  }
  
  if(is.vector(matrix)) {
    if(length(matrix) != ceiling(years)){
      matrix <- c(rep(matrix[1],ceiling(years) - length(matrix)),matrix)
    }
    matrix <- as.matrix(matrix)
  }
  
  return(matrix)
}

#' plaf matrix check
#' 
#' checks length formatting etc
#' 
#' @param plaf plaf matrix
#' @param years Years desired

plaf_matrix_check <- function(matrix, years) {
  
  if(is.matrix(matrix)){
    if(dim(matrix)[1] != ceiling(years)){
      matrix <- rbind(matrix,matrix(rep(matrix[2,],ceiling(years) - dim(matrix)[1]),ncol=dim(matrix)[2],byrow=TRUE))
    }
  }
  
  if(is.vector(matrix)) {
    
    matrix <- c(rep(matrix,years - 1),matrix)
    matrix <- matrix(matrix,nrow=ceiling(years),byrow=TRUE)
  }
  
  return(matrix)
}

#' Intervention grab
#' 
#' Grabs ITN, IRS  ft from database
#' 
#' @param country Country string
#' @param admin Admin string
#' @param year_range Year range
#' @param final_itn_cov Final ITN coverage
#' @param final_irs_cov Final IRS coverage
#' @param final_ft Final treatment coverage
#' 

intervention_grab <- function(country, admin, year_range, 
                              final_itn_cov = NULL,
                              final_irs_cov = NULL,
                              final_ft = NULL) {
  
  
  itn_2000_2015 <- magenta::itn_2000_2015
  irs_2000_2015 <- magenta::irs_2000_2015
  admin_units_seasonal <- magenta::admin_units_seasonal
  
  # incorporate hitoric interventions
  itn_cov <- itn_2000_2015$value[which(itn_2000_2015$admin==admin & itn_2000_2015$country==country)]
  irs_cov <- irs_2000_2015$value[which(irs_2000_2015$admin==admin & irs_2000_2015$country==country)]
  ft <- admin_units_seasonal$ft[which(admin_units_seasonal$admin1==admin & admin_units_seasonal$country==country)] 
  ft <- rep(ft, 16)
  
  # check it's the right length
  if(length(irs_cov)!=16){
    irs_cov <- irs_2000_2015$value[which(irs_2000_2015$country==country)]
  }
  if(length(irs_cov)!=16){
    irs_cov <- rep(0,16)
  } 
  
  # grab year bits
  end_year <- tail(year_range, 1)
  start_year <- year_range[1]
  years <- 2000:2015
  
  # check the final coverages
  if (is.null(final_itn_cov)) final_itn_cov <- min(tail(itn_cov, 1), 0.9)
  if (is.null(final_irs_cov)) final_irs_cov <- min(tail(irs_cov, 1), 0.9)
  if (is.null(final_ft)) final_ft <- min(tail(final_ft, 1), 0.9)
  
  # alter the end 
  if (end_year > 2015) {
    
    itn_cov <- c(itn_cov,final_itn_cov)
    irs_cov <- c(irs_cov,final_irs_cov)
    ft <- c(ft,final_ft)
    
    itn_cov <- approx(x = c(years, end_year), y = itn_cov, xout = 2000:end_year)$y
    irs_cov <- approx(x = c(years, end_year), y = irs_cov, xout = 2000:end_year)$y
    ft <- approx(x = c(years, end_year), y = ft, xout = 2000:end_year)$y
    
  } else if (end_year <= 2015 && end_year > 2000) {
    
    itn_cov <- itn_cov[1:(end_year - 1999)]
    irs_cov <- irs_cov[1:(end_year - 1999)]
    ft <- ft[1:(end_year - 1999)]
    
  } else {
    
    itn_cov <- rep(0, length(year_range))
    irs_cov <- rep(0, length(year_range))
    ft <- rep(0, length(year_range))
    return(list("itn_cov" = itn_cov, 
                "irs_cov" = irs_cov,
                "ft" = ft))
    
  }
  
  if (year_range[1] < 2000) {
    itn_cov <- c(rep(0, length(start_year:1999)), itn_cov)
    irs_cov <- c(rep(0, length(start_year:1999)), irs_cov)
    ft <- c(rep(0, length(start_year:1999)), ft)
  } else {
    itn_cov <- itn_cov[which(years == start_year):length(itn_cov)]
    irs_cov <- irs_cov[which(years == start_year):length(irs_cov)]
    ft <- ft[which(years == start_year):length(ft)]
  }
  
  
  return(list("itn_cov" = itn_cov, 
              "irs_cov" = irs_cov,
              "ft" = ft))
  
}

#' Create mosquito death and feeding vector
#' 
#' Grabs spatial incidence and mosquitoFOI vector from database
#' 
#' @param eqInit Equilibrium Inititial Conditions
#' @param admin Admin string
#' @param ft Annual Treatment Coverage vector
#' @param itn_cov Annual ITN Coverage vector
#' @param irs_cov Annual IRS Coverage vector
#' @param years Numeric for total years
#' 

mu_fv_create <- function(eqInit,
                         ft = ft,
                         itn_cov = itn_cov,
                         irs_cov = irs_cov,
                         int_times = NULL,
                         years = years,
                         full = FALSE,
                         complete = TRUE,
                         odin_model = system.file("extdata/odin_model_itn_irs.R",package="magenta")) {
  
  if (identical(irs_cov, 0) && identical(itn_cov, 0)) {
    out <- data.frame("mu"=rep(0.132, length(seq_len(round(years*365))+1)),
                      "fv"=rep(0.333, length(seq_len(round(years*365))+1)))
  } else {
    
    message("\nRunning Deterministic Model for ", years, " years")
    
    length_check <- function(x) {
      if (length(x) == 0) {
        x <- rep(0, years)
      } 
      if(length(x) < years) {
        x <- rep(x, length.out = years)
      } 
      if(length(x) > years) {
        x <- tail(x, years)
      }
      return(c(x[1],x))
    }
    
    # check interventions parm lengths
    eqInit$ft_vector <- length_check(ft)
    eqInit$itn_vector <- length_check(itn_cov)
    eqInit$irs_vector <- length_check(irs_cov)
    
  
    # set up the intervetnion times
    if (is.null(int_times)) {
      int_times <- c(-25,seq(0,by=365,length.out = years))
    }
    eqInit$t_vector <- int_times
    eqInit$ITN_IRS_on <- 0
    
    
    # split pop accordingly
    eqInit$pop_split <- rep(1/eqInit$num_int, eqInit$num_int)
    
    extensions <- which(unlist(lapply(eqInit, function(x) {(length(dim(x))==3)})))
    mat <- matrix(0, eqInit$na, eqInit$nh)
    for(e in names(extensions)) {
      if(grepl("_I",e)){
        eqInit[[e]] <- array(eqInit[[e]][,,1] , c(eqInit$na, eqInit$nh, eqInit$num_int))
      } else {
        eqInit[[e]] <-  vapply(rep(1/eqInit$num_int,eqInit$num_int), FUN = function(x) { x * eqInit[[e]][,,1]}, mat)
      }
    }
    
    # build model 
    odin_model_path <- odin_model
    gen <- odin::odin(odin_model_path,verbose=FALSE)
    state <- eqInit[names(eqInit) %in% names(formals(gen))]
    model <- gen(user=state,use_dde=TRUE)
    
    #create model and simulate
    tt <- seq(1,round(years*365),1)
    mod_run <- model$run(t = tt,n_history = 1000)
    out <- model$transform_variables(mod_run)
    
    if(!full){
      out <- data.frame("mu"=out$mu,"fv"=out$fv)
    }
    
  }
  
  return(out)
  
}




#' Create housekeeping parameter list
#' 
#' List for simulation housekeeping vars, e.g. quiet prints,
#' 
#' @param quiet Boolean for quiet simulation
#' @param cluster Boolean for simulation being on cluster
#' 

housekeeping_list_create <- function(quiet = TRUE,
                                     cluster = TRUE) {
  
  l <- list("quiet_print" = quiet,
            "cluster" = cluster)
  
  return(l)
  
}


#' Create drug list
#' 
#' List for simulating drug usage for resistance/mft variables
#' 
#' @param resistance_flag Boolean are we simulating resistance
#' @param number_of_resistance_loci Numeric for number of res. loci
#' @param resistance_costs Numeric vector for costs of res. Should be 2^res.loci
#' @param epistatic_logic Is there compensatory relationships. i.e. what loci 
#'   need to be true for resistance cost to exist. Default of NULL means that 
#'   this becomes seq_len(number_of_resistance_loci), i.e. only dependent on 
#'   their own loci. (TODO: Change this to be a list of length norl)
#' @param prob_of_lpf Numeric vector for lpf, Should be 2^res.loci * num_drugs
#' @param mft_flag Boolean are we doing mft
#' @param temporal_cycling Numeric for when in years a drug switch occurs
#' @param sequential_cycling Numeric for what perc. treatment failure before switch
#' @param sequential_update How long does it take in years for sequential to be implemented
#' @param number_of_drugs Numeric for number of drugs used
#' @param drug_choice What's the default drug choice to begin. Default = 0
#' @param partner_drug_ratios Numeric vector for ratio of first line drugs used
#' 
#' @export

drug_list_create <- function(resistance_flag = FALSE,
                             number_of_drugs = 2,
                             number_of_resistance_loci = 3,
                             resistance_costs = c(0.99,0.99,0.98),
                             epistatic_logic = NULL, 
                             prob_of_lpf = list(c(1.0,0.97,0.80,0.55),
                                                c(1.0,0.98,0.7,0.51)),
                             barcode_res_pos = list(c(0,1),
                                                    c(0,2)),
                             prophylactic_pos = list(c(1),c(2)),
                             dur_P = rep(25, number_of_drugs),
                             dur_SPC = rep(5, number_of_drugs),
                             mft_flag = FALSE,
                             temporal_cycling = -1,
                             sequential_cycling = -1,
                             sequential_update = 3,
                             drug_choice = 0,
                             partner_drug_ratios = rep(1/number_of_drugs,number_of_drugs)) {
  
  
  # TODO: More check for formatting here that you have enough prophylactic etc
  
  if(temporal_cycling > 0 && sequential_cycling > 0) {
    stop ("Both sequential and temporal cycling can't be greater than 0")
  }
  
  resistance_costs <- resistance_cost_table_create(number_of_resistance_loci, 
                                                   resistance_costs,
                                                   epistatic_logic)
  
  
  l <- list("g_resistance_flag" = resistance_flag,
            "g_number_of_resistance_loci" = number_of_resistance_loci,
            "g_cost_of_resistance" = resistance_costs,
            "g_prob_of_lpf" = prob_of_lpf,
            "g_barcode_res_pos" = barcode_res_pos,
            "g_prophylactic_pos" = prophylactic_pos,
            "g_dur_P" = dur_P,
            "g_dur_SPC" = dur_SPC,
            "g_mft_flag" = mft_flag,
            "g_temporal_cycling" = temporal_cycling,
            "g_next_temporal_cycle" = temporal_cycling,
            "g_sequential_cycling" = sequential_cycling,
            "g_sequential_update" = sequential_update,
            "g_number_of_drugs" = number_of_drugs,
            "g_drug_choice" = drug_choice,
            "g_partner_drug_ratios" = partner_drug_ratios)
  
  return(l)
  
}

lpf_table_create <- function(number_of_drugs, drug_tables){
  
  max_bits <- (2^(number_of_drugs+1))-1
  barcode_poss <- matrix(intToBits(0:max_bits),ncol=32,byrow=T)[,1:(number_of_drugs+1)]
  
  results <- rep(1, (2^(number_of_drugs+1))*number_of_drugs)
  count <- 1
  
  for(i in 1:number_of_drugs){
    for(j in 1:nrow(barcode_poss)){
      results[count] <- drug_tables[[i]][bitsToInt(barcode_poss[j,c(1,i+1)])+1]
      count <- count+1
    }
  }
  
  return(results)
  
}

resistance_cost_table_create <- function(number_of_resistance_loci, 
                                         resistance_costs, 
                                         epistatic_logic=NULL){
  
  if (is.null(epistatic_logic)) {
    epistatic_logic <- seq_len(number_of_resistance_loci)
  }
  
  max_bits <- (2^(number_of_resistance_loci))-1
  barcode_poss <- matrix(intToBits(0:max_bits),ncol=32,byrow=T)[,1:(number_of_resistance_loci)]
  results <- rep(1, nrow(barcode_poss))
  count <- 1
  
  for(j in 1:nrow(barcode_poss)){
    comp <- epistatic_logic[as.logical(barcode_poss[j,])]
    res <- as.logical(barcode_poss[j,])
    res[which(res)] <- comp %in% which(res)
    results[count] <- prod(resistance_costs[res])
    count <- count+1
  }
  
  return(results)
  
}


drug_list_update <- function(drug_list, year, tf){
  
  
  if (drug_list$g_resistance_flag) {
    if (!drug_list$g_mft_flag) {
      
      # if sequential cycling
      if (drug_list$g_sequential_cycling > 0) {
        if (tf > drug_list$g_sequential_cycling && year > drug_list$g_next_temporal_cycle) {
          drug_list$g_drug_choice <- drug_list$g_drug_choice + 1
          drug_list$g_next_temporal_cycle <- year + drug_list$g_sequential_update
          if(drug_list$g_drug_choice == (drug_list$g_number_of_drugs)){
            drug_list$g_drug_choice <- 0
          }
          
          message(drug_list$g_drug_choice)
        }
      }
      
      # if temporal cycling
      if (drug_list$g_temporal_cycling > 0) {
        if (drug_list$g_next_temporal_cycle == year) {
          drug_list$g_drug_choice <- drug_list$g_drug_choice + 1
          if(drug_list$g_drug_choice == (drug_list$g_number_of_drugs)){
            drug_list$g_drug_choice <- 0
          }
          drug_list$g_next_temporal_cycle <- drug_list$g_next_temporal_cycle + 
            drug_list$g_temporal_cycling
        }
      }
    }
  }
  
  return(drug_list)  
  
} 


# function to calculate probability of strain being detected by microscopy
q_fun <- function(d1, ID, ID0, kD, fd) {
  return(d1 + ((1-d1) / (1 + ((ID/ID0)^kD)*fd)))
}

fd <- function(age, fD0, aD, gammaD) {
  return( 1 - ((1-fD0) / (1 + (age/aD)^gammaD)) )
}

# convert allele frequencies
pop_alf <- function(nums,nl,weighted=FALSE){ 
  if(class(nums %>% unlist)=="raw") {
    if(weighted) {
    res <- colSums(do.call(rbind,lapply(nums,function(x){colMeans(matrix(as.numeric(x),ncol = nl))})))/length(nums)
    } else {
    res <- colMeans(matrix(as.numeric(do.call(rbind,nums)),ncol=nl))
    }
  } else {
    res <- rep(NA,nl)
  }
  return(unlist(res))
}

# convert to lineages frequencies
lineages <- function(nums,nl){ 
  if(class(nums %>% unlist)=="raw") {
    res <- table(
      factor(
        apply(matrix(as.numeric(do.call(rbind,nums)),ncol=nl),1,bitsToInt),
        levels=0:((2^nl)-1)
      )
    )
  } else {
    res <- rep(NA,2^nl)
  }
  return(unlist(res))
}


# what to save in update_saves
update_saves <- function(res, i, sim.out, sample_states,
                         sample_size, sample_reps, mean_only,
                         barcode_params, num_loci, full_update_save=FALSE,
                         genetics_df_without_summarising, save_lineages = FALSE,
                         human_update_save, summary_saves_only,
                         only_allele_freqs, mpl, seed, weighted=FALSE) {
  
  
  # what are saving, does it include the humans
  if(human_update_save) 
  {
    # do we just want the summary data frame 
    if(summary_saves_only){
      
      # sample the population strain's genetics
      if(length(sample_size)>1){
        df <- pop_strains_df(sim.out$Ptr, sample_size = 0, 
                             sample_states = sample_states, ibd = barcode_params$barcode_type,
                             seed = seed,
                             nl = num_loci)
      } else {
        df <- pop_strains_df(sim.out$Ptr, sample_size = sample_size*sample_reps, 
                             sample_states = sample_states, ibd = barcode_params$barcode_type,
                             seed = seed, 
                             nl = num_loci)
        
      }
      
      # summarise the genetics
      if(genetics_df_without_summarising) {
        
        res[[i]] <- list()
        if(only_allele_freqs){
          res[[i]]$af <- pop_alf(df$nums[unlist(lapply(df$barcode_states,length))>0],num_loci,weighted=weighted)
          if(save_lineages){
            res[[i]]$lineage <-lineages(df$nums[unlist(lapply(df$barcode_states,length))>0],num_loci)
          }
        } else {
          res[[i]]$df <- df
        }
        
      } else {
        
        if(i%%12 == 0 && i >= (length(res)-180)){
          res[[i]] <- COI_df_create(df, barcodes=TRUE, nl=num_loci, 
                                    ibd = barcode_params$barcode_type,
                                    n = sample_size, reps = sample_reps, 
                                    mean_only = mean_only)
        } else {
          res[[i]] <- COI_df_create(df, barcodes=FALSE, nl=num_loci, 
                                    ibd = barcode_params$barcode_type,
                                    n = sample_size, reps = sample_reps, 
                                    mean_only = mean_only)
        }
        
      }
      
      # treatment outcomes for resistance work
      res[[i]]$succesfull_treatments <- sim.out$Loggers$Treatments$Successful_Treatments
      res[[i]]$unsuccesful_treatments_lpf <- sim.out$Loggers$Treatments$Unsuccesful_Treatments_LPF
      res[[i]]$not_treated <- sim.out$Loggers$Treatments$Not_Treated
      res[[i]]$treatment_failure <-  res[[i]]$unsuccesful_treatments_lpf / (res[[i]]$unsuccesful_treatments_lpf + res[[i]]$succesfull_treatments) 
      res[[i]]$daily_pop_eir <-  sim.out$Loggers$Treatments$daily_infectious_bite_counters/sim.out$Loggers$Log_Counter*365
      
      if(any(is.na(res[[i]]$treatment_failure))) {
        res[[i]]$treatment_failure[is.na(res[[i]]$treatment_failure)] <- 0 
      }
      
      # prevalence measures
      two_ten <- sim.out$Loggers$Ages < 3650 & sim.out$Loggers$Ages > 720
      dt_10 <- sum(sim.out$Loggers$InfectionStates %in% c(1, 4) & two_ten) / sum(two_ten)
      
      micro_det <- q_fun(mpl$d1,sim.out$Loggers$ID,mpl$ID0,mpl$kD,
                         sapply(sim.out$Loggers$Ages,fd,mpl$fD0,mpl$aD,mpl$gammaD)) 
      
      a_10 <- sum(sim.out$Loggers$InfectionStates %in% 2 & two_ten ) / sum(two_ten) * mean(micro_det[two_ten])
      res[[i]]$pcr_prev <- sum(unlist(sim.out$Loggers[c("D","A","U","T")]))
      res[[i]]$dat <- sum(unlist(sim.out$Loggers[c("D","A","T")]))
      res[[i]]$micro_2_10 <- dt_10 + a_10
      
      # or do we want the full human popualation
    } else {
      
      # then grab the population
      strain_vars <- c(
        "Strain_infection_state_vectors",
        "Strain_day_of_infection_state_change_vectors",
        "Strain_barcode_vectors"
      )
      human_vars <- c("Infection_States", "Zetas", "Ages")
      
      pl3 <- param_list_simulation_get_create(statePtr = sim.out$Ptr)
      sim.save <- simulation_R(pl3, seed = seed)
      
      if(full_update_save) {
        res[[i]] <-sim.save
      } else {
      
      # and then store what's needed
      Strains <- sim.save$populations_event_and_strains_List[strain_vars]
      Humans <- c(sim.save$population_List[human_vars],
                  Strains,
                  sim.save$scourge_List["Mosquito_Day_of_death"],
                  sim.save$scourge_List["Mosquito_Infection_States"])
      res[[i]] <- Humans
      }
    }
    
    # or are we just saving the loggers
  } 
  else 
  {
    res[[i]] <- sim.out$Loggers
  }
  
  return(res)
  
}


#' Create vector adaptation list
#' 
#' List for vector adaptations relating to oocyst success
#' 
#' @param vector_adaptation_flag Boolean are we doing vector adaptation.
#' @param local_oocyst_advantage Double for multiplicative efffect on oocyst 
#'   chance formation.
#' @param gametocyte_non_sterilisation Double for reduction in a wild type's 
#'   onward contribution to infection due to sterilisation due to artemisinin. 
#' 

vector_adaptation_list_create <- function(vector_adaptation_flag = FALSE,
                                          local_oocyst_advantage = 0.2,
                                          gametocyte_non_sterilisation = 0.2){
  
  l <- list("g_vector_adaptation_flag" = vector_adaptation_flag,
            "g_local_oocyst_advantage" = local_oocyst_advantage,
            "g_gametocyte_non_sterilisation" = gametocyte_non_sterilisation)
  
  return(l)
  
}


#' Create nmf list
#' 
#' List for simulating non malarial fever
#' 
#' @param nmf_flag Boolean are we doing non malarial fevers
#' @param mean_nmf_frequency Vector for mean number of days between fevers for
#'   the age bracket considered 
#' @param nmf_age_brackets Vector for age brackets
#' @param prob_of_testing_nmf Numeric for probability that a NMF is tested by 
#'   RDT before being treated with antimalarials. 
#' 
#' 

nmf_list_create <- function(nmf_flag = FALSE,
                            mean_nmf_frequency = c(148.578,139.578,141.564,155.874,179.364,
                                                   216.192,233.478,268.056,312.858,315.564,
                                                   285.156,255.246,238.302,216.618),
                            nmf_age_brackets = c(-0.1, 365.0, 730.0, 1095.0, 1460.0, 1825.0, 
                                                 2555.0, 3285.0, 4015.0, 4745.0, 5475.0, 
                                                 7300.0, 9125.0, 10950.0, 36850.0),
                            prob_of_testing_nmf = 0.5){
  
  l <- list("g_nmf_flag" = nmf_flag,
            "g_mean_nmf_frequency" = mean_nmf_frequency,
            "g_nmf_age_brackets" = nmf_age_brackets,
            "g_prob_of_testing_nmf" = prob_of_testing_nmf)
  
  return(l)
  
}