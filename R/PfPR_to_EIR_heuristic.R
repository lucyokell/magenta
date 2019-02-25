#------------------------------------------------
#' Brute force conversion of PfPR and ft to an EIR
#'
#' @param PfPR Parasite prevaence by PCR
#' @param PfPR_micro Parasite prevaence by microscopy (Default = NULL)
#' @param mv Mosquito size if it's known. Default = NULL
#' @param ft Treatment seeking
#' @param ... Any other params to be fed to the model parameter list
#' 

PfPR_to_EIR_heuristic <- function(PfPR=NULL,PfPR_micro=NULL,mv=NULL,
                                  ft, ...){
  
  mpl <- model_param_list_create(eta = 1/(21*365),...)
  
  num_het_brackets <- 5
  num_age_brackets <- 20
  max_age <- 100
  ## Create the geometric age brackets
  ratio <- (max_age/0.1)^(1/num_age_brackets)
  age.vector <- 0.1 * ratio ** (1:num_age_brackets)
  age.vector[1] <- 0
  
  bm <- read.csv(system.file("extdata/bm.txt",package = "magenta"),sep=",",header = T)
  
  if(!is.null(PfPR)){
    PR <- PfPR
    EIR <- bm$EIRY_eq[which.min(abs(bm$pcr_pos_all.final..1 - PR))]
  } else if ((!is.null(PfPR_micro))) {
    PR <- PfPR_micro
    EIR <- bm$EIRY_eq[which.min(abs(bm$slide_pos_2_10 - PR))]
  } else{
    PR <- mv
    EIR <- bm$EIRY_eq[which.min(abs(bm$mv - PR))]
  }

  
  not_found <- TRUE
  
  #while(not_found){
  for(i in 1:100){
    
    ## Create a near equilibirum initial condition
    eqInit <- equilibrium_init_create(age.vector = age.vector,
                                      het.brackets = num_het_brackets,
                                      ft = ft,
                                      EIR = EIR,
                                      model.param.list = mpl)
    
    if(!is.null(PfPR)){
      diff <- PR - (1 - sum(eqInit[[1]]) - sum(eqInit[[6]]))
    }  else if ((!is.null(PfPR_micro)))  {
      diff <- PR - (sum(eqInit$D[eqInit$age2years: eqInit$age10years,,] + 
                          eqInit$T[eqInit$age2years: eqInit$age10years,,] + 
                          eqInit$A[eqInit$age2years: eqInit$age10years,,] *eqInit$p_det_eq[eqInit$age2years: eqInit$age10years,,]) / 
                      sum(eqInit$den[eqInit$age2years: eqInit$age10years]))
    } else {
      diff <- PR - eqInit$mv0
    }
    
    if(abs(diff)<0.001)
    {
      not_found <- FALSE
    } 
    else 
    {
      
      if(diff<0){
        if(abs(diff) > 0.01){
          EIR <- EIR - 0.2
        } else if (abs(diff) > 0.005) {
          EIR <- EIR - 0.1  
        } else if (abs(diff) > 0.0025) {
          EIR <- EIR - 0.05  
        } else {
          EIR <- EIR - 0.001
        }
      } 
      else
      {
        if(abs(diff) > 0.01){
          EIR <- EIR + 0.2
        } else if (abs(diff) > 0.005) {
          EIR <- EIR + 0.1  
        } else if (abs(diff) > 0.0025) {
          EIR <- EIR + 0.05 
        } else {
          EIR <- EIR + 0.001
        }
      }
    }
    
    
  }
  
  return(EIR)
  
}