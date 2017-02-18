## Typical Cluster script

install.packages("didewin",
                 repos=c(CRAN="https://cran.rstudio.com",
                         drat="https://richfitz.github.io/drat"))

devtools::install_github(c(
  "richfitz/ids",
  "richfitz/syncr",
  "dide-tools/context",
  "richfitz/queuer",
  "dide-tools/didewin"))

install.packages("openssl")

didewin::didewin_config_global(credentials="C:\\Users\\Oliver\\.smbcredentials",
                               temp=didewin::path_mapping("tmp","T:","//fi--didef2.dide.ic.ac.uk/tmp","T:"),
                               cluster="fi--dideclusthn")

didewin::didewin_config(rtools=TRUE)
didewin::web_login()

sources.vector <- c(paste("R/",list.files("M:/OJ/MAGENTA/R"),sep=""))

packages.vector <- c("Rcpp","stringi","statmod", "magrittr","ring","dde","odin","MAGENTA")

root1 <- "M:/OJ/MAGENTA_Results/EIR_prev_inc_test"
ctx1 <- context::context_save(sources=sources.vector, root=root1,
                              package_sources=context::package_sources(github=c("richfitz/ring","richfitz/dde","richfitz/odin"),local="M:/OJ/MAGENTA"),
                              packages=packages.vector)
obj1 <- didewin::queue_didewin(ctx1, config = didewin::didewin_config(rtools = TRUE))

bm <- read.csv("inst/extdata/bm.txt",sep="\t")
cpdf <- data.frame("EIR"=as.numeric(bm$EIRY_eq),"N"=1e6,"years"=20)
grp1 <- queuer::enqueue_bulk(X = cpdf,FUN = Pipeline, obj = obj1, timeout=0)
