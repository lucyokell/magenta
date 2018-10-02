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
    if(dim(matrix)[1] != years){
      matrix <- rbind(matrix(rep(matrix[1,],years - dim(matrix)[1]),ncol=dim(matrix)[2],byrow=TRUE),matrix)
    }
  }
  
  if(is.vector(matrix)) {
    if(length(matrix) != years){
      matrix <- c(rep(matrix[1],years - length(matrix)),matrix)
    }
    matrix <- as.matrix(matrix)
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
                         full = FALSE) {
  
  if (identical(irs_cov, 0) && identical(itn_cov, 0)) {
    out <- data.frame("mu"=rep(0.132, length(seq_len(years*365)+1)),
                      "fv"=rep(0.333, length(seq_len(years*365)+1)))
  } else {
  
  message("\nRunning Deterministic Model for ", years, " years")
  
  length_check <- function(x) {
    if (length(x) == 0) {
      x <- rep(0, years)
    } 
    if(length(x) < years) {
      x <- rep(x, length.out = years)
    } 
    return(x)
  }
  
  # check interventions parm lengths
  eqInit$ft <- length_check(ft)
  eqInit$itn_cov <- length_check(itn_cov)
  eqInit$irs_cov <- length_check(irs_cov)
  
  # set up the intervetnion times
  if (is.null(int_times)) {
    int_times <- seq(0,by=365,length.out = years)
    int_times[length(int_times)] <- int_times[length(int_times)]
  }
  eqInit$int_times <- int_times
  eqInit$ITN_IRS_on <- 0
  
  # create odin model
  odin_model_path <- system.file("extdata/odin.R",package="MAGENTA")
  gen <- odin::odin(odin_model_path,verbose=FALSE,build = TRUE)
  model <- generate_default_model(ft=eqInit$ft,age=eqInit$age_brackets,dat=eqInit,generator=gen,dde=TRUE)
  
  #create model and simulate
  tt <- seq(0,years*365,1)
  mod_run <- model$run(tt)
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
                                     cluster = FALSE) {
  
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
#' @param prob_of_lpf Numeric vector for lpf, Should be 2^res.loci * num_drugs
#' @param mft_flag Boolean are we doing mft
#' @param temporal_cycling Numeric for when in years a drug switch occurs
#' @param sequential_cycling Numeric for what perc. treatment failure before switch
#' @param number_of_drugs Numeric for number of drugs used
#' @param drug_choice What's the default drug choice to begin. Default = 0
#' @param partner_drug_ratios Numeric vector for ratio of first line drugs used
#' 
#' 

drug_list_create <- function(resistance_flag = FALSE,
                             number_of_resistance_loci = 3,
                             resistance_costs = c(1,rep(0.99,3),rep(0.98,3),0.97),
                             prob_of_lpf = c(c(1.0,0.97,0.80,0.55,1.00,0.97,0.80,0.55),
                                             c(1.0,0.97,1.00,0.97,0.80,0.55,0.80,0.55)),
                             mft_flag = FALSE,
                             temporal_cycling = -1,
                             sequential_cycling = -1,
                             number_of_drugs = 2,
                             drug_choice = 0,
                             partner_drug_ratios = rep(1/number_of_drugs,number_of_drugs)) {
  
  
  if(temporal_cycling > 0 && sequential_cycling > 0) {
    stop ("Both sequential and temporal cycling can't be greater than 0")
  }
  
  
  prob_of_lpf <- matrix(prob_of_lpf, nrow = number_of_drugs, byrow=TRUE)
  
  l <- list("g_resistance_flag" = resistance_flag,
            "g_number_of_resistance_loci" = number_of_resistance_loci,
            "g_cost_of_resistance" = resistance_costs,
            "g_prob_of_lpf" = prob_of_lpf,
            "g_mft_flag" = mft_flag,
            "g_temporal_cycling" = temporal_cycling,
            "g_next_temporal_cycle" = temporal_cycling,
            "g_sequential_cycling" = sequential_cycling,
            "g_number_of_drugs" = number_of_drugs,
            "g_drug_choice" = drug_choice,
            "g_partner_drug_ratios" = partner_drug_ratios)
  
  return(l)
  
}


drug_list_update <- function(drug_list, year, tf){

  
  if (drug_list$g_resistance_flag) {
    if (!drug_list$g_mft_flag) {
      
      # if sequential cycling
      if (drug_list$g_sequential_cycling > 0) {
        if (tf > drug_list$g_sequential_cycling && year > drug_list$g_next_temporal_cycle) {
          drug_list$g_drug_choice <- drug_list$g_drug_choice + 1
          drug_list$g_next_temporal_cycle <- year + 2
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
        drug_list$g_next_temporal_cycle <- drug_list$g_next_temporal_cycle + drug_list$g_temporal_cycling
      }
    }
    }
  }
  
  return(drug_list)  
  
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
                            mean_nmf_frequency = c(148.578,139.578,141.564,155.874,179.364,216.192,233.478,268.056,312.858,315.564,285.156,255.246,238.302,216.618),
                            nmf_age_brackets = c(-0.1, 365.0, 730.0, 1095.0, 1460.0, 1825.0, 2555.0, 3285.0, 4015.0, 4745.0, 5475.0, 7300.0, 9125.0, 10950.0, 36850.0),
                            prob_of_testing_nmf = 0.5){
  
  l <- list("g_nmf_flag" = nmf_flag,
            "g_mean_nmf_frequency" = mean_nmf_frequency,
            "g_nmf_age_brackets" = nmf_age_brackets,
            "g_prob_of_testing_nmf" = prob_of_testing_nmf)
  
  return(l)
  
}