#' magenta: Individual-based simulation model of malaria epidemiology and genomics.
#'
#' @description 
#' magenta is an individual-based simulation model of malaria epidemiology and genomics. 
#' magenta extends the imperial malaria model by tracking the infection history of individuals. 
#' With this additional genetic characteristics of the parasite can be assessed.
#'
#' @docType package
#' @name magenta
#' 
#' @useDynLib magenta
#' @importFrom stats rlnorm heatmap runif qt sd
#' @importFrom utils adist read.csv tail
#' @importFrom ggplot2 ggplot
#' @importFrom grDevices dev.off tiff dev.new hcl
#' @importFrom odin odin
#'
"_PACKAGE"

globalVariables(c("admin_units_seasonal","irs_2000_2015","itn_2000_2015",""))