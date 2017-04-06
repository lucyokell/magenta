contextualise <- function(root = "M:/OJ/MAGENTA_Results/Mosquito_EIR_Prev_Check_max_age_with_plusmodf_and_gamtcg_and_new_seeds"){

didehpc::didehpc_config_global(credentials="C:\\Users\\Oliver\\.smbcredentials",
                               temp=didehpc::path_mapping("tmp","T:","//fi--didef2.dide.ic.ac.uk/tmp","T:"),
                               home=didehpc::path_mapping("OJ","M:","//fi--didef2.dide.ic.ac.uk/Malaria","M:"))
didehpc::web_login()
didehpc::didehpc_config()

packages.vector <- c("Rcpp","stringi","statmod", "magrittr","ring","dde","odin","MAGENTA")

context::context_log_start()
ctx <- context::context_save(root,
                             packages = packages.vector,
                             package_sources= provisionr::package_sources(github=c("richfitz/ring","richfitz/dde","richfitz/odin"),
                                                                          local="M:/OJ/MAGENTA"))
config <- didehpc::didehpc_config(use_workers = TRUE)
obj <- didehpc::queue_didehpc(ctx, config = config)
didehpc::web_login()
obj <- didehpc::queue_didehpc(ctx, config = config)

return(obj)
}

workers <- obj$submit_workers(1)

bm <- read.csv("inst/extdata/bm.txt",sep="\t")
pos <- floor(seq(1,50,length.out = 20))
# 
cpdf <- data.frame("EIR"=as.numeric(seq(1,99,length.out = 50)),"N"=1e4,"years"=20, full_save=FALSE,yearly_save=FALSE)
#cpdf <- data.frame("EIR"=as.numeric(bm$EIRY_eq),"N"=1e5,"years"=100, full_save=TRUE,yearly_save=FALSE)
grp <- queuer::enqueue_bulk(X = cpdf,FUN = Pipeline, obj = obj, timeout=0, name="EIR_1_to_100_20years",overwrite = T)
# 
cpdfeq <- data.frame("EIR"=c(rep(3,10)),"N"=1e5,"years"=100, full_save=TRUE,yearly_save=FALSE)
grp2 <- queuer::enqueue_bulk(X = cpdfeq,FUN = Pipeline, obj = obj, timeout=0, name="EIR_3_100years",overwrite = T)

cpdfeq2 <- data.frame("EIR"=c(rep(1,5)),"N"=1e4,"years"=1, full_save=TRUE,yearly_save=FALSE)
grp3 <- queuer::enqueue_bulk(X = cpdfeq2,FUN = Pipeline, obj = obj, timeout=0, name="Seed_Checks",overwrite = T)

path <- "M:/OJ/MAGENTA_Results/1_10_100_Eqs/"
ends <- paste0(path,c(paste0("EIR=3_",1:10),paste0("EIR=10_",1:10),paste0("EIR=100_",1:10)),".rds")
orderedends <- c(rep(ends[1:10],100),rep(ends[11:20],100),rep(ends[21:30],100))

cpdfeqcont <- data.frame("human_only_full_save" = TRUE,"yearly_save" = FALSE, 
                         "saved_state_path"=orderedends,years=10, stringsAsFactors = FALSE)
didehpc::web_login()
obj <- didehpc::queue_didehpc(ctx, config = config)
workers <- obj$submit_workers(40)
grp1 <- queuer::enqueue_bulk(X = cpdfeqcont[1:1000,],FUN = Pipeline, obj = obj, timeout=0, name="EIR_1_cont_10years",overwrite = T)
grp2 <- queuer::enqueue_bulk(X = cpdfeqcont[1001:2000,],FUN = Pipeline, obj = obj, timeout=0, name="EIR_10_cont_10years",overwrite = T)
grp3 <- queuer::enqueue_bulk(X = cpdfeqcont[2001:2700,],FUN = Pipeline, obj = obj, timeout=0, name="EIR__100(1_7)_cont_10years",overwrite = T)
grp4 <-  queuer::enqueue_bulk(X = cpdfeqcont[2701:3000,],FUN = Pipeline, obj = obj, timeout=0, name="EIR__100(8_10)_cont_10years",overwrite = T)
grp5 <-  queuer::enqueue_bulk(X = cpdfeqcont[1:1000,],FUN = Pipeline, obj = obj, timeout=0, name="EIR__3_cont_10years",overwrite = T)
## Will still need to go through and resubmit tasks that require 8,9,10

## UGANDA STUDY

root <- "M:/OJ/MAGENTA_Results/Uganda_COI_study"
obj <- contextualise(root)
workers <- obj$submit_workers(30)
EIR1 <- PfPR_to_EIR_heuristic(0.64,0.47)
EIR2 <- PfPR_to_EIR_heuristic(0.20,0.457)
EIR3 <- PfPR_to_EIR_heuristic(0.23,0.260)
cpdfuganda <- data.frame("EIR"=c(rep(EIR1,10),rep(EIR2,10),rep(EIR3,10)),"ft" = c(rep(0.47,10),rep(0.457,10),rep(0.23,10)),"N"=1e5,"years"=60, human_only_full_save=TRUE,yearly_save=FALSE)
#workers <- obj$submit_workers(20)
grp1 <- queuer::enqueue_bulk(X = cpdfuganda,FUN = Pipeline, obj = obj, timeout=0, name="Uganda_COI",overwrite = T)
