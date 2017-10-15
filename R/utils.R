#' Pipe operator
#'
#' See \code{\link[magrittr]{\%>\%}} for more details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
 

#------------------------------------------------
#' Create seasonal theta return
#'
#' @param country Character for country within which admin2 is in. Default = NULL
#' @param admin Character for admin region. Some fuzzy logic will be used to match. If 
#' not provided then no seasonality is introduced. Default = NULL
#' 
seasonal_profile <- function(country,admin){
  
  #data("admin_units_seasonal")
  # intiialise admin match as no match
  admin.matches <- 0
  if(!is.null(admin)){
    # if there is no country given then search for the admin unit
    if(is.null(country)){
      # find exact match
      admin.matches <- grep(paste("^",admin,"\\b",sep=""),admin_units_seasonal$admin1)
      # if exact does not match try fuzzy match up to dist of 4 which should catch having nop spaces or separators etc
      if(length(admin.matches)==0){
        admin.matches <- which(adist(admin_units_seasonal$admin1,admin)<=4)
      }
      if(length(admin.matches)>1) stop("Admin unit string specified is ambiguous without country")
      # if we do have a country though find that match first and then find admin
    } else {
      # first find an exact match
      country.matches <- grep(paste("^",country,"\\b",sep=""), admin_units_seasonal$country)
      if(length(unique(admin_units_seasonal$country[country.matches]))==1){
        chosen.country <- unique(admin_units_seasonal$country[country.matches])
      } else if(length(unique(admin_units_seasonal$country[country.matches]))==0){
        # if exact does not match try fuzzy match up to dist of 2 which should catch having no spaces or separators etc
        country.matches <- which(adist(admin_units_seasonal$country,y = country)<=2)
        if(length(unique(admin_units_seasonal$country[country.matches]))==1){
          chosen.country <- unique(admin_units_seasonal$country[country.matches])
        } else if(length(unique(admin_units_seasonal$country[country.matches]))==0) stop ("Country string specified not close enough to those in database")
      }
      # find exact match
      admin.sub.matches <- grep(paste("^",admin,"\\b",sep=""),admin_units_seasonal$admin1[country.matches])
      # if exact does not match try fuzzy match up to dist of 4 which should catch having nop spaces or separators etc
      if(length(admin.sub.matches)==0){
        admin.sub.matches <- which(adist(admin_units_seasonal$admin1[country.matches],admin)<=4)
      }
      if(length(admin.sub.matches)>1) stop("Admin unit string specified is not close enougth to those in the database")
      admin.matches <- which(admin_units_seasonal$admin1 == admin_units_seasonal$admin1[country.matches][admin.sub.matches])
    }
  }
  if(admin.matches!=0){
    ssa0 <- admin_units_seasonal$a0[admin.matches]
    ssa1 <- admin_units_seasonal$a1[admin.matches]
    ssa2 <- admin_units_seasonal$a2[admin.matches]
    ssa3 <- admin_units_seasonal$a3[admin.matches]
    ssb1 <- admin_units_seasonal$b1[admin.matches]
    ssb2 <- admin_units_seasonal$b2[admin.matches]
    ssb3 <- admin_units_seasonal$b3[admin.matches]
    theta_c <- admin_units_seasonal$theta_c[admin.matches]
    
    t <- 1:365
    y <- sapply(t,function(x) max((ssa0+ssa1*cos(2*pi*x/365)+ssa2*cos(2*2*pi*x/365)+ssa3*cos(3*2*pi*x/365)+ssb1*sin(2*pi*x/365)+ssb2*sin(2*2*pi*x/365)+ ssb3*sin(3*2*pi*x/365) ) /theta_c,0.001))
    
  } else {
    ssa0 <- ssa1 <- ssa2 <- ssa3 <- ssb1 <- ssb2 <- ssb3 <- theta_c <- 0
    y <- rep(1,365)
  }
  
  return(y)
}