didehpc::didehpc_config_global(credentials="C:\\Users\\Oliver\\.smbcredentials",
  temp=didehpc::path_mapping("tmp","T:","//fi--didef2.dide.ic.ac.uk/tmp","T:"),
  home=didehpc::path_mapping("OJ","M:","//fi--didef2.dide.ic.ac.uk/Malaria","M:"))
didehpc::web_login()
didehpc::didehpc_config()

packages.vector <- c("Rcpp","stringi","statmod", "magrittr","ring","dde","odin","MAGENTA")
unique_value <- "ojwatson"


root <- "M:/OJ/MAGENTA_Results/Mosquito_EIR_Prev7"
context::context_log_start()
ctx <- context::context_save(root,
  packages = packages.vector,
  package_sources= provisionr::package_sources(github=c("richfitz/ring","richfitz/dde","richfitz/odin"),
    local="M:/OJ/MAGENTA"))
config <- didehpc::didehpc_config(use_workers = TRUE)
obj <- didehpc::queue_didehpc(ctx, config = config)

workers <- obj$submit_workers(30)
bm <- read.csv("inst/extdata/bm.txt",sep="\t")
pos <- floor(seq(1,50,length.out = 30))
cpdf <- data.frame("EIR"=as.numeric(bm$EIRY_eq[pos]),"N"=1e5,"years"=60, full_save=FALSE)
grp <- queuer::enqueue_bulk(X = cpdf,FUN = Pipeline, obj = obj, timeout=0, name="EIR_Check",overwrite = T)

bundle <- obj$task_bundle_get("EIR_Check")
