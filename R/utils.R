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
      admin.matches <- country.matches[admin.sub.matches]
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


#------------------------------------------------
#' Create summary statistics df
#'
#' @param data Dataframe to be summarised
#' @param measurevar Character for measure variable
#' @param groupvars Characters of grouping variables
#' @param na.rm Boolean to remove NAs. 
#' @param conf.interval numeric for CI
#' @param .drop drop within ddply
#' 
#' @importFrom plyr ddply rename
#' 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
 # library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


###
#------------------------------------------------
#' Create COI dataframe
#'
#' @param sim_save Output of a simulation save
#' @param groupvars Grouping vars for summarySE. Default = c("Age_Bin","State")
#' @param barcodes Boolean whether to return tabled barcodes. Default = FALSE
#' 
COI_df_create <- function(sim_save, groupvars = c("Age_Bin","State"),barcodes=FALSE){
  
  # Convert the barcodes to COIs
  COIS <- Convert_Barcode_Vectors(sim_save$populations_event_and_strains_List,sub_patents_included=TRUE)
  COIs <- COIS$COI
  COIs[is.na(COIs)] <- 0
  clonality <- table(table(unlist(COIS$nums)))
  barcodes <- sort(table(unlist(COIS$nums)),decreasing=TRUE)
  
  # grab the age, status etc 
  ages <- sim_save$population_List$Ages
  ages[ages==0] <- 0.001
  clinical_status <- sim_save$population_List$Infection_States
  infection_state <- c("S","D","A","U","T","P")
  last_treatment <- sim_save$parameters_List$g_current_time - sim_save$populations_event_and_strains_List$Day_of_last_treatment
  last_treatment[last_treatment==0] <- 0.001
  
  #  bring into df
  df <- data.frame("COI"=COIs,"Ages"=ages/365,"State"=infection_state[(clinical_status)+1],
                   "Age_Bin" = cut(ages/365,breaks = c(0,1,2,3,5,10,20,40,60,100)),
                   "Last_Treatment" = last_treatment,
                   "Last_Treatment_Binned" = cut(last_treatment,breaks = c(0,28,90,365,sim_save$parameters_List$g_current_time)),
                   stringsAsFactors = FALSE)
  
  # Add the non subpatent COIs
  COIS <- Convert_Barcode_Vectors(sim_save$populations_event_and_strains_List,sub_patents_included=FALSE)
  COIs <- COIS$COI
  COIs[is.na(COIs)] <- 0
  df$COI_Detected_From_Model <- COIs
  
  summary <- summarySE(df,measurevar = "COI",groupvars = groupvars)
  
  if(barcodes){
    return(list("Summary"=summary,"Clonality"=clonality,"Barcodes"=barcodes[barcodes>1])) 
  } else {
  return(list("Summary"=summary,"Clonality"=clonality))
  }
}


#' tester
#' 
#' @export
#' 
tester <- function(){ 
  Sys.setenv(BINPREF="T:/Rtools/Rtools33/mingw_64/bin/")
  odin::can_compile()
}
NULL


#------------------------------------------------
#' Function to convert vector of bools representing up to a 32 bit number into the corresponding integer
#'
#' @param x Vector of logicals
bitsToInt<-function(x) {
  packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")
}


#------------------------------------------------
#' Function to generate ggplot colours
#'
#' @param n Number of colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

