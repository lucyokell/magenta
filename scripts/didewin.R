didehpc::didehpc_config_global(credentials="C:\\Users\\Oliver\\.smbcredentials",
  temp=didehpc::path_mapping("tmp","T:","//fi--didef2.dide.ic.ac.uk/tmp","T:"),
  home=didehpc::path_mapping("OJ","M:","//fi--didef2.dide.ic.ac.uk/Malaria","M:"))
didehpc::web_login()
didehpc::didehpc_config()

packages.vector <- c("Rcpp","stringi","statmod", "magrittr","ring","dde","odin","MAGENTA")
unique_value <- "ojwatson"
sources = "R/pipeline.R"

root <- "M:/OJ/MAGENTA_Results/Mosquito_EIR_Prev_Check3_nofloating"
context::context_log_start()
ctx <- context::context_save(root,
  packages = packages.vector,
  sources=sources,
  package_sources= provisionr::package_sources(github=c("richfitz/ring","richfitz/dde","richfitz/odin"),
    local="M:/OJ/MAGENTA"))
config <- didehpc::didehpc_config(use_workers = TRUE)
obj3 <- didehpc::queue_didehpc(ctx, config = config)
didehpc::web_login()

obj <- didehpc::queue_didehpc(ctx, config = config)
workers <- obj$submit_workers(30)
bm <- read.csv("inst/extdata/bm.txt",sep="\t")
pos <- floor(seq(1,50,length.out = 30))
#cpdf <- data.frame("EIR"=as.numeric(bm$EIRY_eq[pos]),"N"=1e5,"years"=60, full_save=TRUE,yearly_save=TRUE)
#grp <- queuer::enqueue_bulk(X = cpdf,FUN = Pipeline, obj = obj, timeout=0, name="EIR_Check_10_save",overwrite = T)
cpdf2 <- data.frame("EIR"=as.numeric(bm$EIRY_eq),"N"=1e5,"years"=60, full_save=FALSE,yearly_save=FALSE)
grp2 <- queuer::enqueue_bulk(X = cpdf2,FUN = Pipeline, obj = obj, timeout=0, name="EIR_Check_Full_60_years",overwrite = T)

cpdfeq <- data.frame("EIR"=c(rep(1,10),rep(10,10),rep(100,10)),"N"=1e5,"years"=100, full_save=TRUE,yearly_save=FALSE)
grp3 <- queuer::enqueue_bulk(X = cpdfeq,FUN = Pipeline, obj = obj, timeout=0, name="EIR_1_10_100_100years",overwrite = T)
