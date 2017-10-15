#------------------------------------------------
#' Brute force conversion of PfPR and ft to an EIR
#'
#' @param PfPR Parasite prevaence by PCR
#' @param ft Treatment seeking
#' 

PfPR_to_EIR_heuristic <- function(PfPR,ft){
  
  mpl <- Model_Param_List_Create(eta = 1/(21*365))
  
  num_het_brackets <- 5
  num_age_brackets <- 20
  max_age <- 100
  ## Create the geometric age brackets
  ratio <- (max_age/0.1)^(1/num_age_brackets)
  age.vector <- 0.1 * ratio ** (1:num_age_brackets)
  age.vector[1] <- 0
  
  bm <- read.csv("inst/extdata/bm.txt",sep="\t")
  
  EIR <- bm$EIRY_eq[which.min(abs(bm$pcr_pos_all.final..1 - PfPR))]
  
  not_found <- TRUE
  
  #while(not_found){
  for(i in 1:100){
    
    ## Create a near equilibirum initial condition
    eqInit <- Equilibrium_Init_Create(age.vector = age.vector,
                                      het.brackets = num_het_brackets,
                                      ft = ft,
                                      EIR = EIR,
                                      model.param.list = mpl)
    
    
    diff <- PfPR - (1 - sum(eqInit[[1]]) - sum(eqInit[[6]]))
    
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