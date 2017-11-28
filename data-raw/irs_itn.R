irs_raw <- xml2::as_list(xml2::read_xml(system.file("MAGENTA",lib.loc = "data-raw/irs2000-15.xml")))

intervention <- lapply(irs_raw,function(x) x$parameters$Type[[1]]) %>% unlist
year <- lapply(irs_raw,function(x) x$parameters$StartTime[[1]]) %>% unlist
value <- lapply(irs_raw,function(x) x$parameters$Level[[1]]) %>% unlist
country_admin <- lapply(irs_raw,function(x) x$adunits$A[[1]]) %>% unlist
country <- strsplit(country_admin,split = ":",fixed=TRUE) %>% lapply(function(x) x[1]) %>% unlist
admin <- strsplit(country_admin,split = ":",fixed=TRUE) %>% lapply(function(x) x[2]) %>% unlist

years <- 2000:2015
year <- years[match(year,unique(year))]

irs_2010_2015 <- list("intervention"=intervention,
               "country"=country,
               "admin"=admin,
               "year"=year,
               "value"=value)

devtools::use_data(irs_2010_2015)


itn_raw <- xml2::as_list(xml2::read_xml(system.file("MAGENTA",lib.loc = "data-raw/itn2000-15.xml")))

intervention <- lapply(itn_raw,function(x) x$parameters$Type[[1]]) %>% unlist
year <- lapply(itn_raw,function(x) x$parameters$StartTime[[1]]) %>% unlist
value <- lapply(itn_raw,function(x) x$parameters$Level[[1]]) %>% unlist
country_admin <- lapply(itn_raw,function(x) x$adunits$A[[1]]) %>% unlist
country <- strsplit(country_admin,split = ":",fixed=TRUE) %>% lapply(function(x) x[1]) %>% unlist
admin <- strsplit(country_admin,split = ":",fixed=TRUE) %>% lapply(function(x) x[2]) %>% unlist

years <- 2000:2015
year <- years[match(year,unique(year))]

itn_2010_2015 <- list("intervention"=intervention,
                 "country"=country,
                 "admin"=admin,
                 "year"=year,
                 "value"=value)

devtools::use_data(itn_2010_2015)
