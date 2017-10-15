#' MAGENTA: Individual-based simulation model of malaria epidemiology and genomics.
#'
#' @description 
#' MAGENTA is an individual-based simulation model of malaria epidemiology and genomics. 
#' MAGENTA extends the imperial malaria model by tracking the infection history of individuals. 
#' With this additional genetic characteristics of the parasite can be assessed.
#'
#' @docType package
#' @name MAGENTA
#' 
#' @useDynLib MAGENTA
#' @importFrom stats rlnorm heatmap runif
#' @importFrom utils adist read.csv
#' @importFrom ggplot2 ggplot
#' @importFrom grDevices dev.off tiff windows dev.new
#'
"_PACKAGE"

globalVariables(c("admin_units_seasonal","irs_2000_2015","itn_2000_2015"))