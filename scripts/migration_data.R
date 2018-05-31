# Importation Data from 2000 - 2017

ff = data.table::fread(system.file("extdata/importation.xls","MAGENTA"),sep="\t")
ad = data.table::fread(system.file("extdata/adunit.txt","MAGENTA"),sep="\t",colClasses = "character")

# format au for matching later
ad$au <- (substr(ad$V1,0,4))
ad$cc <- countrycode:::countrycode(ad$V2,origin="country.name",destination = "iso2c")
ad$cc_au <- paste0(ad$cc,ad$V3,sep="|")
ad$cc_au_npcap <- toupper(rdhs:::rm_punct_non_ascii(ad$cc_au))

# create results list
importations <- list() ; length(importations) <- ff$t %>% unique %>% length
names(importations) <- 2000:2017
years <- ff$t %>% unique

# loop through and normliase importation maatrices and give unique cc_au names
for(i in 1:length(importations)){
  
  message(i)
  
  # set up loop's results list
  importations[[i]] <- list() ; length(importations[[i]]) <- 2
  names(importations[[i]]) <- c("incidence","mosquitoFOI")
  
  # grab this year's data and make sure the names are the correct 4 characters
  ff2 = ff[ff$t == years[i],]
  ff2$Home_au[which(nchar(ff2$Home_au)==3)] <- paste0("0",ff2$Home_au[which(nchar(ff2$Home_au)==3)])
  ff2$Infecting_au[which(nchar(ff2$Infecting_au)==3)] <- paste0("0",ff2$Infecting_au[which(nchar(ff2$Infecting_au)==3)])
  
  ######################################################################################
  # make ff (flow file) into matrices
  ######################################################################################
  
  nau = length(unique(ff2$Home_au))
  
  # entry i,j -> number of humans in au i that got infected when visiting au j
  ffincI <- matrix(as.numeric(ff2$incI), nrow = nau, ncol = nau, byrow=TRUE)
  
  # entry i,j -> number of mosquitoes in admin unit i that were infected by visitors from au j
  ffmosFOI <- matrix(ff2$mosFOI, nrow = nau, ncol = nau, byrow=TRUE)
  
  # now normalise these to sum to 1
  ffincI2 <- apply(ffincI,MARGIN = 1,FUN = function(x) x/sum(x)) %>% t
  ffmosFOI2 <- apply(ffmosFOI,MARGIN = 1,function(x) x/sum(x,na.rm=TRUE)) %>% t
  ffincI2[is.na(ffincI2)] <- 0
  ffmosFOI2[is.na(ffmosFOI2)] <- 0
  diag(ffincI2) <- NA ; diag(ffmosFOI2) <- NA
  
  # make unique names
  cc_au_names <- apply(ad[match(unique(ff2$Home_au),ad$au),c("cc","V3")],1,paste,collapse="|") %>% unique
  colnames(ffincI2) <- rownames(ffincI2) <- colnames(ffmosFOI2) <- rownames(ffmosFOI2) <- cc_au_names
  
  importations[[i]]$incidence <- ffincI2
  importations[[i]]$mosquitoFOI <- ffmosFOI2
  
}

devtools::use_data(importations)
