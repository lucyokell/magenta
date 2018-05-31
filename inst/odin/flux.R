mod <- odin::odin({
  
  deriv(A) <- 2 * A * cov2 
  initial(A) <- 1
  
  lagB <- A * 2
  B <- delay(lagB,10)
  
  ###
  
  deriv(C) <- C * B
  initial(C) <- 1
  
  lagD <- C * cov[1]
  D <- delay(lagD,12) 

  
  flux <- interpolate(flux_t, flux_y, "linear")
  dim(cov) <- 2
  cov[1] <- 1 - flux
  cov[2] <- flux
  cov2 <- sum(cov) * D
  
  ###
  
  flux_t[] <- user()
  flux_y[] <- user()
  dim(flux_t) <- user()
  dim(flux_y) <- user()
  
})