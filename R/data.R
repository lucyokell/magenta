##'  Admin level 1 africa seasonal parameters
##'
##'  These datasets represent the data fitted within the Imperial College Malaria model for
##'  relating seasonal profiles to malaria transmission intensity at level 1 admin regions
##'  across Africa
##'
##' @docType data
##'
##' @format A dataframe of 576 observations of 13 variables:
##'
##'   \code{admin_units_seasonal}: A dataframe of admin units and their seasonal parameters
##'     \itemize{
##'       \item country: Country string
##'       \item admin1: Admin 1 string
##'       \item map_prev_2010: 2010 Atlas map microscopy prevalence in 2-10 year olds
##'       \item id: Numeric vector of 1:576
##'       \item ID_1: Numeric vector referencing shape file geoshape ids
##'       \item a0: Average value of fourier series
##'       \item a1: First of partial series cos terms
##'       \item b1: First of partial series sine terms
##'       \item a2: Second of partial series cos terms
##'       \item b2: Second of partial series sine terms
##'       \item a3: Third of partial series cos terms
##'       \item b3: Third of partial series sine terms
##'       \item theta_c: Rainfall normalising constant
##'     }
##'
##' @rdname admin_units_seasonal
##' @aliases admin_units_seasonal
##'
##'
"admin_units_seasonal"

##'  IRS for 2000 to 2015
##'
##'  IRS data for SSSA for 2000 to 2015
##'
##' @docType data
##'
##' @format A dataframe of 5 elements:
##'
##'   \code{$irs_2000_2015}: A dataframe of admin units and their seasonal parameters
##'     \itemize{
##'       \item intervention: String stating IRS
##'       \item country: Country string
##'       \item admin: Admin string
##'       \item year: Numeric year
##'       \item value: Value for IRS coverage
##'     }
##'
##' @rdname irs_2010_2015
##' @aliases irs_2010_2015
##'
##'
"irs_2000_2015"

##'  ITN for 2000 to 2015
##'
##'  ITN data for SSSA for 2000 to 2015
##'
##' @docType data
##'
##' @format A dataframe of 5 elements:
##'
##'   \code{itn_2000_2015}: A dataframe of admin units and their seasonal parameters
##'     \itemize{
##'       \item intervention: String stating ITN
##'       \item country: Country string
##'       \item admin: Admin string
##'       \item year: Numeric year
##'       \item value: Value for ITN coverage
##'     }
##'
##' @rdname itn_2000_2015
##' @aliases itn_2000_2015
##'
##'
"itn_2000_2015"

##'  Importation data 
##'
##'  Importation data for admin units in admin_units_seasonal
##'
##' @docType data
##'
##' @format A list of length 18
##'
##'   \code{importations}: A list of length 18, with each list representing a year. In each
##'   year is then a further 2 lists which are:
##'     \itemize{
##'       \item incidence: Proportions of incidence that originated from other admin units, i.e.
##'       individuals in admin i that moved to admin j and then returned with an infection acquired
##'       while in admin j.  
##'       \item mosquitoFOI: Proportion of the force of infection towards mosquitoes that originated
##'       from outside admin units, i.e. the proportion of mosquitoes that are infected from infected
##'       individuals who travelled from admin j into admin i. 
##'     }
##'     
##' @rdname importations
##' @aliases importations
##'
##'
"importations"


