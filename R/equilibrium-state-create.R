#' Equilibrium steady state list creation
#'
#' \code{Equilibrium_SS_Create} creates a close to steady-state system by 
#' running an ODE compartment version of the Imperial Malaria model, using
#' as an initialisation the output of \code{Equilibrium_Init_Create}
#' 
#' @param eqInit Equilibrium initialisation state as output from \code{Equilibrium_Init_Create}
#' @param end.year Number of years to run ODE model for. Default = 5.
#' @param use_odin Boolean detailing whether the intiial solution is run within the odin model 
#' for end.year length. Default = False while the model is still buggy. 
#' 
#' @export

Equilibrium_SS_Create <- function(eqInit, end.year = 5, use_odin = FALSE){
  
  if(use_odin){
    
    ## Odin generator function
    generate_default_model <- function(ft,age,dat,generator,dde = TRUE){
      mod <- generator(init_S=dat$S,
                       init_T=dat$T,
                       init_D=dat$D,
                       init_A=dat$A,
                       init_U=dat$U,
                       init_P=dat$P,
                       init_ICA = dat$ICA,
                       init_ICM = dat$ICM,
                       init_ID = dat$ID,
                       init_IB = dat$IB,
                       init_Sv = dat$Sv,
                       init_Ev = dat$Ev,
                       init_Iv = dat$Iv,
                       init_PL = dat$PL,
                       init_LL = dat$LL,
                       init_EL = dat$EL,
                       na = dat$na,
                       nh = dat$nh,
                       age_rate = dat$age_rate,
                       foi_age = dat$foi_age,
                       het_wt = dat$het_wt,
                       rel_foi = dat$rel_foi,
                       omega = dat$omega,
                       pi = pi,
                       x_I = dat$x_I,
                       #eov = dat$eov,
                       mv0 = dat$mv0,
                       mu0 = dat$mu0,
                       ssa0 = dat$ssa0,
                       ssa1 = dat$ssa1,
                       ssa2 = dat$ssa2,
                       ssa3 = dat$ssa3,
                       ssb1 = dat$ssb1,
                       ssb2 = dat$ssb2,
                       ssb3 = dat$ssb3,
                       theta_c = dat$theta_c,
                       den = dat$den,
                       age59 = dat$age59,
                       age05 = dat$age05,
                       age = age*365,
                       ft = ft,
                       use_dde = dde,
                       eta = dat$eta,
                       rA = dat$rA,
                       rT = dat$rT,
                       rD = dat$rD,
                       rU = dat$rU,
                       rP = dat$rP,
                       dE = dat$dE,
                       delayGam = dat$delayGam,
                       delayMos = dat$delayM,
                       cD = dat$cD,
                       cT = dat$cT,
                       cU = dat$cU,
                       gamma1 = dat$gamma1,
                       d1 = dat$d1,
                       dID = dat$dID,
                       ID0 = dat$ID0,
                       kD = dat$kD,
                       uD = dat$uD,
                       aD = dat$aD,
                       fD0 = dat$fD0,
                       gammaD = dat$gammaD,
                       b0 = dat$b0,
                       b1 = dat$b1,
                       dB = dat$dB,
                       IB0 = dat$IB0,
                       kB = dat$kB,
                       uB = dat$uB,
                       phi0 = dat$phi0,
                       phi1 = dat$phi1,
                       dCA = dat$dCA,
                       IC0 = dat$IC0,
                       kC = dat$kC,
                       uCA = dat$uCA,
                       dCM = dat$dCM,
                       tau1 = dat$tau1,
                       tau2 = dat$tau2,
                       Q0 = dat$Q0,
                       chi = dat$chi,
                       bites_Bed = dat$bites_Bed,
                       bites_Indoors = dat$bites_Indoors,
                       p10 = dat$p10,
                       p2 = dat$p2,
                       muEL = dat$muEL,
                       muLL = dat$muLL,
                       muPL = dat$muPL,
                       dEL = dat$dEL,
                       dLL = dat$dLL,
                       dPL = dat$dPL,
                       gammaL = dat$gammaL,
                       beta_larval0 = dat$betaL,
                       num_int = dat$num_int,
                       itn_cov = dat$itn_cov,
                       irs_cov = dat$irs_cov,
                       int_len = length(dat$itn_cov),
                       ITN_IRS_on = dat$ITN_IRS_on,
                       d_ITN0 = dat$d_ITN0,
                       r_ITN0 = dat$r_ITN0,
                       r_ITN1 = dat$r_ITN1,
                       r_IRS0 = dat$r_IRS0,
                       d_IRS0 = dat$d_IRS0,
                       IRS_interval = dat$IRS_interval,
                       ITN_interval = dat$ITN_interval,
                       irs_loss = dat$irs_loss,
                       itn_loss = dat$itn_loss
      )
    }
    
    # Create odin generator
    odin_model_path <- system.file("extdata/odin_model.R",package="MAGENTA")
    gen <- odin::odin(odin_model_path,verbose=FALSE,build = TRUE)
    
    #create model with initial values
    mod <- generate_default_model(ft=eqInit$ft,age=eqInit$age_brackets,dat=eqInit,generator=gen,dde=TRUE)
    tt <- seq(0,end.year*365,1)
    
    # run odin model for end.year years, and if sufficiently steady state found finish
    # if insufficeint steady state found then increase end.year
    steady.state.check <- FALSE
    while(steady.state.check==0){
      # run model
      message(paste("Running model for",end.year,"years"))
      mod_run <- mod$run(tt)
      # shape output
      out <- mod$transform_variables(mod_run)
      
      if(eqInit$ssa0==0){
        
        ICA.corr.mat <- out$ICA[(dim(out$ICA)[1]-100):(dim(out$ICA)[1]),10,1:5,1]
        ICA.corr.mat <- sweep(ICA.corr.mat,2,out$ICA[(dim(out$ICA)[1]),10,1:5,1],`/`)
        if( sum(abs(colMeans(ICA.corr.mat) - 1) < 1e-3)==5){
          steady.state.check <- TRUE
        } else {
          end.year <- end.year + 4
          tt <- seq(0,end.year*365,1)
        }
        
      }
      
    }
    
    ## dimension of model out to kow the end
    final <-   dim(out$ICA)[1]
    ## maternal age position
    maternal <- which.max(eqInit$age_brackets>=20)
    
    ## create equilibrium state for return
    Equilibrium_State <- list(
      "age_brackets" = eqInit$age2,
      "het_brackets" = eqInit$het_bounds,
      "Smat" = out$S[final,,,1],
      "Dmat" = out$D[final,,,1],
      "Amat" = out$A[final,,,1],
      "Umat" = out$U[final,,,1],
      "Tmat" = out$T[final,,,1],
      "Pmat" = out$P[final,,,1],
      "IBmat" = out$IB[final,,,1],
      "ICAmat" = out$ICA[final,,,1],
      "ICMmat" = out$ICM[final,,,1],
      "IDmat" = out$ID[final,,,1],
      "Sv" = out$Sv[final],
      "Ev" = sum(out$Ev[final,]),
      "Iv" = out$Iv[final],
      "MaternalImmunity" = sum(out$ICA[final,maternal,,1] * (eqInit$het_wt))* eqInit$PM,
      "theta" = seasonal_profile(eqInit$country,eqInit$admin)
    )
    
  } else {
    
    maternal <- which.max(eqInit$age_brackets>=20)
    
    ## create equilibrium state for return
    Equilibrium_State <- list(
      "age_brackets" = eqInit$midages,
      "het_brackets" = eqInit$het_bounds,
      #"het_brackets" = eqInit$rel_foi,
      "Smat" = eqInit$S[,,1],
      "Dmat" = eqInit$D[,,1],
      "Amat" = eqInit$A[,,1],
      "Umat" = eqInit$U[,,1],
      "Tmat" = eqInit$T[,,1],
      "Pmat" = eqInit$P[,,1],
      "IBmat" = eqInit$IB[,,1],
      "ICAmat" = eqInit$ICA[,,1],
      "ICMmat" = eqInit$ICM[,,1],
      "IDmat" = eqInit$ID[,,1],
      "Sv" = eqInit$Sv * eqInit$mv0,
      "Ev" = eqInit$Ev * eqInit$mv0,
      "Iv" = eqInit$Iv * eqInit$mv0,
      "MaternalImmunity" = sum(eqInit$ICA[maternal,,1] * (eqInit$het_wt))* eqInit$PM,
      "theta" = seasonal_profile(eqInit$country,eqInit$admin)
    )
    
  }
  
  return(Equilibrium_State)
  
}



#' Create spatial list object to be passed to main functions if required
#'
#' \code{spl_create} creates the user required spatial 
#' list to be passed to the initialisation function if require
#' 

spl_create <- function(spatial_type = 0,
                       human_importation_rate_vector=NULL, 
                       mosquito_imporation_rate_vector=NULL,
                       cotransmission_freq_vector=rep(1,10000),
                       oocyst_freq_vector=rep(1,10000),
                       num_human_infs=0,
                       num_mos_infs=0,
                       plaf = rep(0.5,24)){
  
  
  
  # num_human_infs = (eqInit$FOI * eqInit$den ) %>% sum * N
  # num_mos_infs = eqInit$FOIv_eq * eqInit$mv0 * N
  
  # if non spatial 
  if(spatial_type == 0){
    
    spatial_list <- list("spatial_type" = spatial_type,
                         "cotransmission_freq_vector" = cotransmission_freq_vector,
                         "oocyst_freq_vector" = oocyst_freq_vector,
                         "imported_cotransmissions_events"=0,
                         "imported_oocyst_events"=0,
                         "plaf"=plaf)
    
    # if island spatial  
  } else if(spatial_type == 1){
    
    # first work out how frequent importation is:
    hum_imp_rate <- sum(human_importation_rate_vector,na.rm=TRUE)
    mos_imp_rate <- sum(mosquito_imporation_rate_vector,na.rm=TRUE)
    
    spatial_list <- list("spatial_type" = spatial_type,
                         "cotransmission_freq_vector" = cotransmission_freq_vector,
                         "oocyst_freq_vector" = oocyst_freq_vector,
                         "imported_cotransmissions_events"=hum_imp_rate,
                         "imported_oocyst_events"=mos_imp_rate,
                         "plaf"=plaf)
    
    # if metapop sspatial
  } else {
    
    # first work out how frequent importation is:
    hum_imp_rate <- human_importation_rate_vector[-which(is.na(human_importation_rate_vector))]
    mos_imp_rate <- mosquito_imporation_rate_vector[-which(is.na(mosquito_imporation_rate_vector))]
    
    # therefore we frst want to work out the number of transmission events that are due to importation and how large
    
    cotransmission_import_locations <- table(sample(other_metapopulations,
                                                    size = num_human_infs*sum(human_importation_rate_vector,na.rm=TRUE),
                                                    replace=TRUE,prob=hum_imp_rate))
    imported_cotransmission_frequences <- sapply(cotransmission_import_locations,function(x) sample(cotransmission_vector,size = x))
    
    ## push these freqs to redis as "Exported_Cotransmission_Frequencies_x->metapopulation_number"
    
    ## repeat for oocysts (with some doubling...)
    
    ## grab all the required exports that are coming out of this sim
    
    ## pass that in as an exported_co.... etc
    
    ## and grab all the imported barcodes/oocysts and their frequencies if this is not the initialisation step
    
    ## create spatial_list with 6 vectors
    
    
    imported_cotransmission_events <- sapply(num_human_infs*hum_imp_rate,function(x) rep_len(cotransmission_freq_vector,x))
    imported_oocyst_events <- rep_len(oocyst_freq_vector,num_mos_infs*mos_imp_rate)
    
    spatial_list <- list("spatial_type" = spatial_type,
                         "cotransmission_freq_vector" = cotransmission_freq_vector,
                         "oocyst_freq_vector" = oocyst_freq_vector,
                         "imported_cotransmissions_events"=c(0,imported_cotransmission_events),
                         "imported_oocyst_events"=c(0,imported_oocyst_events))
    
    
  }
  
}





#' Odin generator function
#'
#' \code{generate_default_model} creates the user required generator for 
#' the odin model. 
#' 
#' @param ft Treatment seeking numeric
#' @param age Vector of ages
#' @param generator generator for model.
#' @param dat Output of \code{equilibrium_init_create} 
#' @param dde Use dde in soliving eventual odin. Default = FALSE
#' 

generate_default_model <- function(ft,age,dat,generator,dde = TRUE){
  mod <- generator(init_S=dat$S,
                   init_T=dat$T,
                   init_D=dat$D,
                   init_A=dat$A,
                   init_U=dat$U,
                   init_P=dat$P,
                   init_ICA = dat$ICA,
                   init_ICM = dat$ICM,
                   init_ID = dat$ID,
                   init_IB = dat$IB,
                   init_Sv = dat$Sv,
                   init_Ev = dat$Ev,
                   init_Iv = dat$Iv,
                   init_PL = dat$PL,
                   init_LL = dat$LL,
                   init_EL = dat$EL,
                   na = dat$na,
                   nh = dat$nh,
                   age_rate = dat$age_rate,
                   foi_age = dat$foi_age,
                   het_wt = dat$het_wt,
                   rel_foi = dat$rel_foi,
                   omega = dat$omega,
                   pi = pi,
                   x_I = dat$x_I,
                   #eov = dat$eov,
                   mv0 = dat$mv0,
                   mu0 = dat$mu0,
                   ssa0 = dat$ssa0,
                   ssa1 = dat$ssa1,
                   ssa2 = dat$ssa2,
                   ssa3 = dat$ssa3,
                   ssb1 = dat$ssb1,
                   ssb2 = dat$ssb2,
                   ssb3 = dat$ssb3,
                   theta_c = dat$theta_c,
                   den = dat$den,
                   age59 = dat$age59,
                   age05 = dat$age05,
                   age2years = dat$age2years,
                   age10years = dat$age10years,
                   two_to_10_length = dat$age10years - dat$age2years,
                   age = age*365,
                   ft = ft,
                   use_dde = dde,
                   eta = dat$eta,
                   rA = dat$rA,
                   rT = dat$rT,
                   rD = dat$rD,
                   rU = dat$rU,
                   rP = dat$rP,
                   dE = dat$dE,
                   delayGam = dat$delayGam,
                   delayMos = dat$delayM,
                   cD = dat$cD,
                   cT = dat$cT,
                   cU = dat$cU,
                   gamma1 = dat$gamma1,
                   d1 = dat$d1,
                   dID = dat$dID,
                   ID0 = dat$ID0,
                   kD = dat$kD,
                   uD = dat$uD,
                   aD = dat$aD,
                   fD0 = dat$fD0,
                   gammaD = dat$gammaD,
                   b0 = dat$b0,
                   b1 = dat$b1,
                   dB = dat$dB,
                   IB0 = dat$IB0,
                   kB = dat$kB,
                   uB = dat$uB,
                   phi0 = dat$phi0,
                   phi1 = dat$phi1,
                   dCA = dat$dCA,
                   IC0 = dat$IC0,
                   kC = dat$kC,
                   uCA = dat$uCA,
                   dCM = dat$dCM,
                   tau1 = dat$tau1,
                   tau2 = dat$tau2,
                   Q0 = dat$Q0,
                   chi = dat$chi,
                   bites_Bed = dat$bites_Bed,
                   bites_Indoors = dat$bites_Indoors,
                   p10 = dat$p10,
                   p2 = dat$p2,
                   muEL = dat$muEL,
                   muLL = dat$muLL,
                   muPL = dat$muPL,
                   dEL = dat$dEL,
                   dLL = dat$dLL,
                   dPL = dat$dPL,
                   gammaL = dat$gammaL,
                   beta_larval0 = dat$betaL,
                   num_int = dat$num_int,
                   itn_cov = dat$itn_cov,
                   irs_cov = dat$irs_cov,
                   int_times = dat$int_times,
                   ITN_IRS_on = dat$ITN_IRS_on,
                   d_ITN0 = dat$d_ITN0,
                   r_ITN0 = dat$r_ITN0,
                   r_ITN1 = dat$r_ITN1,
                   r_IRS0 = dat$r_IRS0,
                   d_IRS0 = dat$d_IRS0,
                   IRS_interval = dat$IRS_interval,
                   ITN_interval = dat$ITN_interval,
                   irs_loss = dat$irs_loss,
                   itn_loss = dat$itn_loss
  )
}

#' Barcode parameter list creation
#'
#' \code{barcode_parms_create} creates list detialing the barcode/genetic
#'   parameters for the model
#' 
#' @param num_loci Number of loci. Default = 24
#' @param ibd_length If we are simulating IBD dynamics, each loci is now 
#'   represented by a bitset of ibd_length. Thus ibd_length needs to be long
#'   enough to ensure that as new identity relationships occur, i.e. 
#'   an importation barcode will be a new identity. e.g. If your population is
#'   1000, we may expect at 80% prevalence, with a mean COI of 3 we will need 
#'   2400 different identities, i.e. 2^ibd_length > 2400. However, keep in mind
#'   importations as these need to be continually new, i.e. if we are simulating
#'   for 30 years, with 3 importations a day, then we will need at least length
#'   to ensure that 2^ibd_length > 2400 + (30*365*3). This will probably be
#'   automatically calculated in the future. If we are not ding IBD, then this
#'   should be 1, which is the defalt. 
#' @param plaf Vector of population level allele frequencies for the barcode. 
#'   Default = rep(0.5, barcode_length)
#' @param prob_crossover Vector of probabilities for crossover events for the 
#'   barcode. Default = rep(0.5, barcode_length)
#' @param starting_ibd Starting ibd. Default = 0
#' 
barcode_parms_create <- function(num_loci = 24,
                                 ibd_length = 1,
                                 plaf = rep(0.5, 24),
                                 prob_crossover = rep(0.5,24),
                                 starting_ibd = 0) {
  
  
  # are we doing ibd or not
  if(ibd_length > 1){
    barcode_type <- 1
  } else {
    barcode_type <- 0
  }
  
  res <- list("num_loci" = num_loci,
              "ibd_length" = ibd_length,
              "starting_ibd" = starting_ibd,
              "plaf" = plaf,
              "prob_crossover" = prob_crossover,
              "barcode_type" = barcode_type)
  
  return(res)
  
}



