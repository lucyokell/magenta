## Newest pipeline demonstration
devtools::load_all(".")

## Create parameter list, changine any key parameters, e.g. the average age
mpl <- MAGENTA::Model_Param_List_Create(eta = 1/(21*365))

## Create a near equilibirum initial condition
eqInit <- MAGENTA::Equilibrium_Init_Create(age.vector = c(0,0.25,0.5,0.75,1,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80,90,100),
                                  het.brackets = 5,ft = 0.4,EIR = 120,model.param.list = mpl)

## Next create the near equilibrium steady state
eqSS <- MAGENTA::Equilibrium_SS_Create(eqInit = eqInit, end.year=5)

## Now check and create the parameter list for use in the Rcpp simulation initialisation, and initialise for 1 year 
## (Simulation_Init doesn't introduce any pending infections so the system needs to equilibriate a bit)
pl <- MAGENTA::Param_List_Simulation_Init_Create(N=100000,years=2,eqSS=eqSS)

## Now run the simulation initialisation
sim.out <- MAGENTA::Simulation_R(paramList = pl)

## Now lets say that we want to continue the simulation for 10 years, logging the PfPR every 2 months and plotting it live

## Set up our time intervals and prevalence vector
interval = 1/12
intervals <- seq(1/12,5,interval)
prev <- rep(0,length(intervals))

## Now loop thrugh the intervals using the paramater update list creation function
for(i in 1:length(intervals)){
  pl2 <- MAGENTA::Param_List_Simulation_Update_Create(years = interval,statePtr = sim.out$Ptr)
  sim.out <- MAGENTA::Simulation_R(pl2)
  prev[i] <- sum(unlist(sim.out$Loggers[match(c("D","A","U","T"),names(sim.out$Loggers))]))
  plot(intervals[1:i],prev[1:i],ylim = c(0,1),xlab = "Years",ylab="PfPR")
}

## If we want to save the simulation to file we will now want to create the parameter list for getting the model state
## before again executing the simulation
pl3 <- MAGENTA::Param_List_Simulation_Get_Create(statePtr = sim.out$Ptr)
sim.save <- MAGENTA::Simulation_R(pl3)
