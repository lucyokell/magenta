#------------------------------------------------
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
                     K0 = dat$K0,
                     mv0 = dat$mv0,
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
                     eov = dat$eov,
                     num_int = dat$num_int,
                     itn_cov = dat$itn_cov,
                     irs_cov = dat$irs_cov,
                     ITN_IRS_on = dat$ITN_IRS_on,
                     d_ITN0 = dat$d_ITN0,
                     r_ITN0 = dat$r_ITN0,
                     r_ITN1 = dat$r_ITN1,
                     d_IRS0 = dat$d_IRS0,
                     IRS_interval = dat$IRS_interval,
                     ITN_interval = dat$ITN_interval,
                     irs_loss = dat$irs_loss,
                     itn_loss = dat$itn_loss
    )
  }
  
  # Create odin generator
  odin_model_path <- system.file("extdata/odin_model.R",package="MAGENTA")
  #odin_model_path <- "M:/OJ/MAGENTA/scripts/odin_model.R"
  gen <- odin::odin(odin_model_path,verbose=FALSE)
  
  #eqInit$Ev <- diff(pexp(q = 1:11,rate=0.132))/ sum(diff(pexp(q = 1:11,rate=0.132))) * eqInit$Ev
  
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
    "age_brackets" = eqInit$age_brackets*365,
    "het_brackets" = eqInit$rel_foi,
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
    "MaternalImmunity" = sum(out$ICA[final,maternal,,1] * (eqInit$het_wt))* mpl$PM
  )
  
  } else {
    
    maternal <- which.max(eqInit$age_brackets>=20)
    
    ## create equilibrium state for return
    Equilibrium_State <- list(
      "age_brackets" = eqInit$age_brackets*365,
      "het_brackets" = eqInit$rel_foi,
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
      "Sv" = eqInit$Sv,
      "Ev" = eqInit$Ev,
      "Iv" = eqInit$Iv,
      "MaternalImmunity" = sum(eqInit$ICA[maternal,,1] * (eqInit$het_wt))* eqInit$PM
    )
    
  }
  
  return(Equilibrium_State)
  
}

