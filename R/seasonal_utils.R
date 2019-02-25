#---
#' Match admin region
#'
#' @param country Character for country within which admin2 is in. 
#'   Default = NULL
#' @param admin Character for admin region. Some fuzzy logic will be used to 
#'   match. If not provided then no seasonality is introduced. Default = NULL
#'
#' @return Row number of the matching region in 
#'   \code{magenta::admin_units_seasonal}
#'

admin_match <- function(admin = NULL, country = NULL) {
  ads <- magenta::admin_units_seasonal
  
  # intiialise admin match as no match
  admin_matches <- 0
  
  if (!is.null(admin)) {
    
    # if there is no country given then search for the admin unit
    if (is.null(country)) {
      
      # find exact match
      admin_matches <- grep(paste("^", admin, "$", sep = ""),
                            ads$admin1,
                            ignore.case = TRUE
      )
      
      # if exact does not match try fuzzy match up to dist of 3 which
      # should catch having nop spaces or separators etc
      
      if (length(admin_matches) == 0) {
        admin_matches <- which(adist(ads$admin1, admin) <= 2)
      }
      if (length(admin_matches) > 1) {
        message(paste(ads$admin1[admin_matches], collapse = " "))
        stop("Admin unit string specified is ambiguous without country")
      }
      
      # if we do have a country though find that match first and then find admin
    } else {
      
      # first find an exact match
      country_matches <- grep(paste("^", country, "$", sep = ""),
                              ads$country,
                              ignore.case = TRUE
      )
      
      if (length(unique(ads$country[country_matches])) == 1) {
        chosen.country <- unique(ads$country[country_matches])
      } else if (length(unique(ads$country[country_matches])) == 0) {
        
        # if exact does not match try fuzzy match up to dist of 2 which
        # should catch having no spaces or separators etc
        country_matches <- which(adist(ads$country, y = country) <= 2)
        if (length(unique(ads$country[country_matches])) == 1) {
          chosen.country <- unique(ads$country[country_matches])
        } else if (length(unique(ads$country[country_matches])) == 0) {
          stop("Country string specified not close enough to those in database")
        }
      }
      
      # find exact match
      admin_sub_matches <- grep(paste("^", admin, "$", sep = ""),
                                ads$admin1[country_matches],
                                ignore.case = TRUE
      )
      
      # if exact does not match try fuzzy match up to dist of 2 which should
      # catch having nop spaces or separators etc
      if (length(admin_sub_matches) == 0) {
        admin_sub_matches <- which(
          adist(ads$admin1[country_matches], admin) <= 2
        )
      }
      if (length(admin_sub_matches) > 1) {
        message(paste(ads$admin1[admin_matches], collapse = " "))
        stop("Admin unit string not close enough to those in the database")
      }
      
      admin_matches <- country_matches[admin_sub_matches]
    }
  }
  
  return(admin_matches)
}



#---
#' Create seasonal theta return
#'
#' @param country Character for country within which admin2 is in. 
#'   Default = NULL
#' @param admin Character for admin region. Some fuzzy logic will be used 
#'   to match. If not provided then no seasonality is introduced. Default = NULL
#'
seasonal_profile <- function(admin = NULL, country = NULL) {
  admin_matches <- admin_match(admin, country)
  admin_units_seasonal <- magenta::admin_units_seasonal
  
  if (admin_matches != 0) {
    ssa0 <- admin_units_seasonal$a0[admin_matches]
    ssa1 <- admin_units_seasonal$a1[admin_matches]
    ssa2 <- admin_units_seasonal$a2[admin_matches]
    ssa3 <- admin_units_seasonal$a3[admin_matches]
    ssb1 <- admin_units_seasonal$b1[admin_matches]
    ssb2 <- admin_units_seasonal$b2[admin_matches]
    ssb3 <- admin_units_seasonal$b3[admin_matches]
    theta_c <- admin_units_seasonal$theta_c[admin_matches]
    
    t <- 1:365
    y <- sapply(t, function(x) {
      max((ssa0 + 
             ssa1 * cos(2 * pi * x / 365) + 
             ssa2 * cos(2 * 2 * pi * x / 365) + 
             ssa3 * cos(3 * 2 * pi * x / 365) + 
             ssb1 * sin(2 * pi * x / 365) + 
             ssb2 * sin(2 * 2 * pi * x / 365) + 
             ssb3 * sin(3 * 2 * pi * x / 365)) / theta_c, 0.001)
    })
  } else {
    ssa0 <- ssa1 <- ssa2 <- ssa3 <- ssb1 <- ssb2 <- ssb3 <- theta_c <- 0
    y <- rep(1, 365)
  }
  
  return(y)
}

#---
#' Plot country seasonal profiles
#'
#' @param country Character for country
#'

country_seasonal_profiles <- function(country) {
  ads <- magenta::admin_units_seasonal
  if (!is.element(country, unique(ads$country))) {
    stop(paste0(
      "Country is not one of ",
      paste(unique(ads$country), collapse = ", ")
    ))
  }
  admins <- as.character(ads$admin1[which(ads$country == country)])
  
  if (length(admins) > 16) {
    par(mfrow = c(5, 5))
  } else if (length(admins) > 9) {
    par(mfrow = c(4, 4))
  } else if (length(admins) > 4) {
    par(mfrow = c(3, 3))
  } else {
    par(mfrow = c(2, 2))
  }
  
  max_ses <- max(sapply(admins, function(x) {
    max(magenta:::seasonal_profile(country, x))
  }))
  
  ret <- sapply(admins, function(x) {
    invisible(plot(magenta:::seasonal_profile(country, x),
                   main = x,
                   ylim = c(0, max_ses + 0.2),
                   ylab = "Seasonal",
                   xlab = "Day"
    ))
  })
}
