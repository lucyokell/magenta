EIR <- 10
ft <- 0.4
years <- 2
N <- 10000
admin = "Tororo"
country <- "Uganda"
num_het_brackets = 5
age <- c(0,0.199526231496888,0.281838293126445,0.398107170553497,0.562341325190349,0.794328234724282,
         1.12201845430196,1.58489319246111,2.23872113856834,3.16227766016838,4.46683592150963,
         6.30957344480194,8.91250938133746,12.5892541179417,17.7827941003892,25.1188643150958,
         35.4813389233576,50.1187233627273,70.7945784384139,100)

mpl <- Model_Param_List_Create()
eqInit <- Equilibrium_Init_Create(age.vector = age,
                                  het.brackets = num_het_brackets,
                                  ft = ft,
                                  EIR = EIR,
                                  country = country,
                                  admin = admin,
                                  model.param.list = mpl)

eqSS <- Equilibrium_SS_Create(eqInit = eqInit, end.year=5, use_odin = FALSE, spatial = round(sum(spatial *N) ))
pl <- Param_List_Simulation_Init_Create(N=N,eqSS=eqSS)
sim.out <- Simulation_R(paramList = pl)

user_vars <- names(formals(odin_model)$user)[-match(c("","age"),names(formals(odin_model)$user))]
temp <- tempfile(fileext = ".txt")
writeLines(text = paste0("model <- odin_model(",paste0(sapply(user_vars, function(x){(paste0(x,"=eqInit$",x,",\n") )}),collapse = ""),
                         "age = eqInit$age*365,\n",
                         "use_dde=TRUE)"),con = temp)
source(temp,local = TRUE)
out <- model$run(t = 1)




tres <- Pipeline(EIR = EIR,ft = ft,N=N,years = 3,update_length = 1,admin = "Tororo",country = "Uganda",
                 num_het_brackets = num_het_brackets,yearly_save = TRUE,num_age_brackets = 20)


r <- redux::hiredis()

tres3 <- Pipeline(EIR = EIR,ft = 0.7,N=N,years = 0.01,update_length = 1,admin = NULL,country = NULL,spatial <- spatial,redis_environment = r,
                 num_het_brackets = num_het_brackets,yearly_save = TRUE,num_age_brackets = 20,human_only_full_save = TRUE)

out <- hanojoel::Run_Model(age,EIR,0.4,"Tororo",20*365,itn_cov=0.3)
out2 <- hanojoel::Run_Model(age,EIR,0.7,NULL,3*365)

tres <- tres3
out <- out2

windows()
plot(colSums(rbind(out$dat$Dout,out$dat$Aout,out$dat$Tout,out$dat$Uout)),
     ylim=c(min(lapply(tres,function(x) sum(x$D,x$A,x$U,x$T)) %>% unlist),
            max((colSums(rbind(out$dat$Dout,out$dat$Aout,out$dat$Tout,out$dat$Uout))))))
plot(lapply(tres[1:(length(tres)-1)],function(x) sum(x$D,x$A,x$T,x$U)) %>% unlist,col="green")

plot(lapply(tres[1:(length(tres)-1)],function(x) c(x$InfectionStates %in% c(1,2,4)) %>% sum) %>% unlist,col="green")

abline(v = 365*(1:10))
abline(h = eqInit$A %>% sum)
what <- Pipeline(EIR = EIR,ft = ft,N=N,years = 0.006,update_length = 1,admin = "Tororo",country = "Uganda",
                 num_het_brackets = num_het_brackets,full_save = TRUE,num_age_brackets = 20)

rbenchmark::benchmark(Pipeline(EIR = EIR,ft = ft,N=N,years = 10,update_length = 1,admin = NULL,country = NULL,
                               num_het_brackets = num_het_brackets,yearly_save=FALSE,full_save = TRUE,num_age_brackets = 20),replications = 2)
  
  
what <- Pipeline(EIR = EIR,ft = ft,N=N,years = 0.3,update_length = 1,admin = NULL,country = NULL,
                 num_het_brackets = num_het_brackets,yearly_save=FALSE,full_save = TRUE,num_age_brackets = 20)


ghana_fit <- Pipeline(EIR=40.6127,ft=0.53283,N=5000,years=1,admin="Upper East",yearly_save = TRUE,
                      update_length=30,barcode_length = 24*17,
                      country="Ghana",num_het_brackets = 5,full_save = TRUE,num_age_brackets = 20)


what$scourge_List$Mosquito_Off_Season %>% table

what$scourge_List$Mosquito_Day_of_next_blood_meal[1:sum(!what$scourge_List$Mosquito_Off_Season)] %>% table
