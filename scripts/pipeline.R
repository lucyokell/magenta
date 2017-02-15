## Pipeline

## Create parameter list, changine any key parameters, e.g. the average age
mpl <- Model_Param_List_Create(eta = 1/(20*365))

## Create a near equilibirum initial condition
eqInit <- Equilibrium_Init_Create(age.vector = c(0,0.25,0.5,0.75,1,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80,90,100),
                                  het.brackets = 5,ft = 0.4,EIR = 120,model.param.list = mpl)

## Next create the near equilibrium steady state
eqSS <- Equilibrium_SS_Create(eqInit = eqInit, end.year=5)

## Now check and create the parameter list for use in the Rcpp simulation
pl <- Param_List_Simulation_Create(N=10000,years=5,eqSS=eqSS)

## Now run the simulation
sim.out <- Simulation_R(paramList = pl)

## To then test the model at different EIRs

# Correct parameter list with mean age of 21
mpl <- Model_Param_List_Create()

# load the BM comparisons
bm <- read.csv("inst/extdata/bm.txt",sep="\t")

## Create list of eqInits
eqInits <- lapply(bm$EIRY_eq[1:2],function(x){return(Equilibrium_Init_Create(age.vector = c(0,0.25,0.5,0.75,1,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80,90,100),
                                                                             het.brackets = 5,ft = 0.4,EIR = x,model.param.list = mpl))})

## Next create the list of near equilibrium steady state
eqSSs <- lapply(eqInits,function(x){return(Equilibrium_SS_Create(eqInit = x, end.year = 5))})

## Next create the list of parameter lists for use in the Rcpp simulation
pls <- lapply(eqSSs,function(x){return(Param_List_Simulation_Create(N=100000, years=10, eqSS = x))})

## Now run the simulations
sim.outs <- lapply(pls,Simulation_R)
