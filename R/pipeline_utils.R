#' Creates age brackets
#' 
#' Makes age brackets given, number of brackets, max age, and spacing
#' 
#' @param max_age
#' @param num_age_brackets
#' @param geometric_age_brackets

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
#' @param country
#' @param admin
#' @param year_range

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
#' @param matrix
#' @param years

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
#' @param country
#' @param admin
#' @param year_range
#' @param final_itn_cov
#' @param final_irs_cov
#' @param final_ft
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
#' @param eqInit
#' @param admin
#' @param ft
#' @param itn_cov
#' @param irs_cov
#' @param years
#' 

mu_fv_create <- function(eqInit,
                         ft = ft,
                         itn_cov = itn_cov,
                         irs_cov = irs_cov,
                         int_times = NULL,
                         years = years) {
  
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
  out <- data.frame("mu"=out$mu,"fv"=out$fv)
  
  }
  
  return(out)
  
}

#' Create housekeeping parameter list
#' 
#' List for simulation housekeeping vars, e.g. quiet prints,
#' 
#' @param quiet
#' @param cluster
#' 

housekeeping_list_create <- function(quiet = TRUE,
                                     cluster = FALSE) {
  
l <- list("quiet_print" = quiet,
          "cluster" = cluster)

return(l)
  
}

