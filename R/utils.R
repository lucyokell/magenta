#' Pipe operator
#'
#' See \code{\link[magrittr:pipe]{\%>\%}} for more details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#---
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

#---
#' Plot country seasonal profiles
#'
#' @param country Character for country within which admin2 is in. Default = NULL
#' 

country_seasonal_profiles <- function(country){
  
  ads <- MAGENTA::admin_units_seasonal
  if(!is.element(country,unique(ads$country))) stop(paste0("Country is not one of ",paste(unique(ads$country),collapse=", ")))
  admins <- as.character(ads$admin1[which(ads$country==country)])
  
  if(length(admins)>16){
    par(mfrow=c(5,5))
  } else if(length(admins)>9){
    par(mfrow=c(4,4))
  } else if(length(admins)>4){
    par(mfrow=c(3,3))
  } else {
    par(mfrow=c(2,2))
  }
  
  max_ses <- max(sapply(admins,function(x){max(MAGENTA:::seasonal_profile(country,x))}))
  
  ret <- sapply(admins,function(x){invisible(plot(MAGENTA:::seasonal_profile(country,x),
                                                  main=x,
                                                  ylim=c(0,max_ses+0.2),
                                                  ylab = "Seasonal",
                                                  xlab = "Day"
  ))})
  
}

#---
#' Create summary statistics df
#'
#' @param data Dataframe to be summarised
#' @param measurevar Character for measure variable
#' @param groupvars Characters of grouping variables
#' @param conf.interval numeric for CI
#' 
#' 
summarySE <- function(data=NULL, measurevar="COI", groupvars=c("Age_Bin","State"),
                      conf.interval=.95) {
  
  
  res <- dplyr::group_by_at(data,dplyr::vars(dplyr::one_of(groupvars))) %>% 
    dplyr::summarise(N=length(!!dplyr::sym(measurevar)),
                     mean=mean(!!dplyr::sym(measurevar),na.rm=TRUE),
                     sd = sd(!!dplyr::sym(measurevar),na.rm=TRUE),
                     se = sd/sqrt(N),
                     ci = se * suppressWarnings(qt(conf.interval/2 + .5, N-1))
    )
  
  return(as.data.frame(res))
}

#---
#' Create summary statistics df
#'
#' @param data Dataframe to be summarised
#' @param measurevar Character for measure variable
#' @param groupvars What are we summarise by
#' @param mean_only Are we just calculating the mean. Default = TRUE
#' 
summarySE_mean_only <- function(data=NULL, measurevar="COI", groupvars=c("Age_Bin","State"),mean_only=TRUE) {
  
  if(!mean_only) {
    return(summarySE(data, measurevar,groupvars))
  } else {
  var <-  dplyr::sym(measurevar)
  
  res <- dplyr::group_by_at(data,groupvars) %>% 
    dplyr::summarise(N=sum(!is.na(!!var)),
                     mean=mean(!!var,na.rm=TRUE))
  
  return(as.data.frame(res))
  }
}

#---
#' Create summary statistics df
#'
#' @param data Dataframe to be summarised
#' @param measurevar Character for measure variable
#' @param groupvars What are we summarise by
#' @param mean_only Are we just calculating the mean. Default = TRUE
#' 
summarySE_mean_only_max_mean <- function(data=NULL, measurevar="COI", groupvars=c("Age_Bin","State"),mean_only=TRUE, max = 6) {
  
  if(nrow(data) != 0) {
  data[measurevar][data[measurevar]>max] <- max
  }
  if(!mean_only) {
    return(summarySE(data, measurevar,groupvars))
  } else {
    var <-  dplyr::sym(measurevar)
    
    res <- dplyr::group_by_at(data,groupvars) %>% 
      dplyr::summarise(N=length(!!var),
                       mean=mean(!!var,na.rm=TRUE))
    
    return(as.data.frame(res))
  }
}



ztrgeomintp <- function(n, mean, p){
  
  ngs <- rgeom(n*1.5,prob = 1/mean)
  ngs <- round(ngs*p)
  ngs <- ngs[ngs>0]
  while(length(ngs) < n) {
    ng2 <- rgeom(n*1.5,prob = 1/mean)
    ngs <- c(ngs,round(ng2*p))
    ngs <- ngs[ngs>0]
  }
  
  ngs <- sample(ngs,size = n,replace=FALSE)
  return(ngs)
  
  
}

ztrnbinom <- function(n,mean,size) {
  
  
  nbs <- rnbinom(n*1.2,size = size, mu = mean)
  nbs <- nbs[nbs>0]
  while(length(nbs) < n) {
    nbs <- c(nbs,rnbinom(n*.5,size = size, mu = mean))
    nbs <- nbs[nbs>0]
  }
  
  nbs <- sample(nbs,size = n,replace=FALSE)
  return(nbs)
}

progress_logging <- function(housekeeping_list, res, progress_bar, 
                             i = 1, initial = FALSE, p_print = 2){
  
  if(initial) {
    if(housekeeping_list$quiet_print) {
      if(housekeeping_list$cluster) {
        p_print <- 2
        message(paste0("||:|"),appendLF = FALSE)
      } else {
        progress_bar$tick()
      }
    }
  } else {
    
    p_perc <- (i/length(res))*100
    
    if(p_perc > p_print) {
      if(housekeeping_list$quiet_print) {
        if(housekeeping_list$cluster) {
          message(paste0(p_print,"|"),appendLF = FALSE)
          if (p_print == 100){
            message(paste0(p_print+2,"|:||"),appendLF = FALSE)
          }
        }
      }
      p_print <- p_print + 2
    }
    
    if (!housekeeping_list$cluster) {
      progress_bar$tick()
    }
  }
  
  return(p_print)
  
}

###
#---
#' Create COI dataframe
#'
#' @param sim_save Output of a simulation save
#' @param groupvars Grouping vars for summarySE. Default = c("Age_Bin","State")
#' @param barcodes Boolean whether to return tabled barcodes. Default = FALSE
#' 
COI_df_create <- function(sim_save, groupvars = c("Age_Bin","State"),
                          barcodes=FALSE,mpl = Model_Param_List_Create(),
                          ibd = 0, nl = 24, n = Inf, sample_group = c("T","D","A","U"), reps = 1){
  
  # function to calculate probability of strain being detected by microscopy
  q_fun <- function(d1, ID, ID0, kD, fd) {
    return(d1 + ((1-d1) / (1 + ((ID/ID0)^kD)*fd)))
  }
  
  fd <- function(age, fD0, aD, gammaD) {
    return( 1 - ((1-fD0) / (1 + (age/aD)^gammaD)) )
  }
  
  # population vars
  ID <- sim_save$population_List$ID
  ages <- sim_save$population_List$Ages
  ages[ages==0] <- 0.001
  micro_det <- q_fun(mpl$d1,ID,mpl$ID0,mpl$kD,sapply(ages,fd,mpl$fD0,mpl$aD,mpl$gammaD)) 
  
  
  if(!ibd){
    # Convert the barcodes to COIs
    COIS <- Convert_Barcode_Vectors(sim_save$populations_event_and_strains_List,sub_patents_included=TRUE,ibd = ibd,nl = nl)
    COIs <- COIS$COI
    COIs[is.na(COIs)] <- 0
    
    # grab the age, status etc 
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
    df$COI_in_A <- COIs
    df$COI_detected_micro <- round(micro_det*df$COI_in_A)
    df$polygenom <- df$COI>1
    df$nums <- COIS$nums
    
    results <- list()
    length(results) <- reps
    pop_size <- length(df$Ages)
    
    # make samples and reps of them
    for(i in seq_len(length(results))) {
      
      # make a sample
      if(n == Inf) {
        chosen <- 1:pop_size
      } else {
        chosen <- sample(which(df$State %in% sample_group), size = n, replace = FALSE)
      }
      
      # create overall clonality
      clonality <- table(table(unlist(COIS$nums[chosen])))
      barcodes_tab <- sort(table(unlist(COIS$nums[chosen])),decreasing=TRUE)
      
      # create summary dfs
      summary_coi <- summarySE(df[chosen,],measurevar = "COI")
      summary_coi_detected <- summarySE(df[chosen,],measurevar = "COI_detected_micro")
      summary_polygenom <- summarySE(df[chosen,],measurevar = "polygenom")
      summary_polygenom <- dplyr::left_join(
        summary_polygenom,
        dplyr::group_by(df[chosen,], Age_Bin, State) %>% dplyr::summarise(unique = clonality_from_barcode_list(nums)), 
        by = groupvars
      )
      
      # bundle it all up
      res <- list("summary_coi"=summary_coi,
                  "summary_coi_detected"=summary_coi_detected,
                  "summary_polygenom"=summary_polygenom,
                  "clonality"=clonality,"coi_table"=table(df$COI[chosen]))
      
      if(barcodes){
        res$barcodes <- barcodes_tab[barcodes_tab>1]
      } 
      
      results[[i]] <- res
    }
    
  } else {
    
    # grab the age, status etc 
    ages <- sim_save$population_List$Ages
    ages[ages==0] <- 0.001
    clinical_status <- sim_save$population_List$Infection_States
    infection_state <- c("S","D","A","U","T","P")
    last_treatment <- sim_save$parameters_List$g_current_time - sim_save$populations_event_and_strains_List$Day_of_last_treatment
    last_treatment[last_treatment==0] <- 0.001
    
    #  bring into df
    out <- sim_save$populations_event_and_strains_List$Recent_identity_vectors
    n.strains <- lapply(out,length) %>% unlist()
    pibd <- rep(0,length(ages))
    if(sum(n.strains)>0){
      ibds <- population_ibd_distances(mat = matrix(unlist(out),nrow=sum(n.strains>0),byrow=TRUE))
      pibd[which(sim_save$populations_event_and_strains_List$Number_of_Strains>0)] <- proxy::colMeans.dist(ibds$p_ibd,diag = FALSE)
    }
    df <- data.frame("pibd"=pibd,"Ages"=ages/365,"State"=infection_state[(clinical_status)+1],
                     "Age_Bin" = cut(ages/365,breaks = c(0,1,2,3,5,10,20,40,60,100)),
                     "Last_Treatment" = last_treatment,
                     "Last_Treatment_Binned" = cut(last_treatment,breaks = c(0,28,90,365,sim_save$parameters_List$g_current_time)),
                     stringsAsFactors = FALSE)
    
    
    results <- list()
    length(results) <- reps
    pop_size <- length(df$Ages)
    
    # make samples and reps of them
    for(i in seq_len(length(results))) {
      
      # make a sample
      if(n == Inf) {
        chosen <- 1:pop_size
      } else {
        chosen <- sample(which(df$State %in% sample_group), size = n, replace = FALSE)
      }
      
    summary_ibd <- summarySE(df[chosen,],measurevar = "pibd",groupvars = groupvars)
    res <- list("summary_ibd"=summary_ibd,
                         "Mean"=mean(pibd[chosen[which(sim_save$populations_event_and_strains_List$Number_of_Strains[chosen]>0)]] ))
    results[[i]] <- res
    }
  }

    # if no reps then unlist
    if(reps == 1) {
      results <- res
    }
    return(results)
    

}

#' @noRd
rbind_list_base <- function(x) {
  x2 <- do.call(
    rbind.data.frame,
    c(x, stringsAsFactors = FALSE, make.row.names = FALSE)
  )
  rownames(x2) <- seq_len(dim(x2)[1])
  x2
}

#' @noRd
ranges <- function(diff, end){
  r <- list();
  for(i in 1:(end/diff)){
    
    r[[i]] <- (1 + ((i-1) * diff)) : (diff*i)
    
  }
  return(r)
}

#' tester
#' 
#' @export
#' 
tester <- function(){ 
  packageVersion("MAGENTA")
}
NULL


#' Function to convert vector of bools representing up to a 32 bit number into the corresponding integer
#'
#' @param x Vector of raw/logicals
#' @param big_endian Is big_endian. Default = FALSE
bitsToInt <- function(x, big_endian = FALSE) {
  if(big_endian) {
    sum(2^(which(rev(as.logical(x)))-1))
  } else {
      sum(2^(which((as.logical(x)))-1))}
  }

# Convert integer to binary for a given n
binary <- function(x, n = 24) {
  i <- 0
  string <- numeric(n)
  while(x > 0) {
    string[n - i] <- x %% 2
    x <- x %/% 2
    i <- i + 1 
  }
  string 
}

#' Function to generate ggplot colours
#'
#' @param n Number of colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# function to put all of one infividuals' parms into one place
person_make <- function(o,n){
  
  if(length(names(o))==6) { 
    
    df <- lapply(o[1:3],function(x) x[n]) %>% as.data.frame
    l <- lapply(o[4:6],function(x) x[n])
    
  } else {
    
    
    df <- lapply(o$population_List,function(x) x[n]) %>% as.data.frame
    oths <- which(lapply(o$populations_event_and_strains_List, class) %>% unlist != "list")
    df <- cbind(df, lapply(o$populations_event_and_strains_List[oths],function(x) x[n])) %>% as.data.frame
    #bs <- grep("barcode",names(o$populations_event_and_strains_List))
    l <- lapply(o$populations_event_and_strains_List[-unique(c(oths))],function(x) x[n])
  }
  return(list("vars"=df,"l"=l))
}


clonality_from_barcode_list <- function(barcode_list){
  
  tbl <- table(table(unlist(lapply(barcode_list,unique))))
  
  if ("1" %in% names(tbl)){
    return(tbl[1]/sum(tbl))
  } else {
    return(0)
  }
  
}

convert_ibd_barcode <- function(b, nl){
  
  ib <- length(b)/nl
  
  br <- seq_len(ib)
  ibd <- rep(0,nl)
  for(i in seq_len(nl)){
    ibd[i] <- bitsToInt(b[br+(ib*(i-1))])
  }  
  
  return(ibd)
  
}

population_ibd_barcodes <- function(barcode_vec,nl){
  
  n.strains <- lapply(barcode_vec,length) %>% unlist()
  lapply(
    barcode_vec[which(n.strains>0)],
    function(x){
      lapply(tail(x,1),convert_ibd_barcode,nl) %>% unlist
    }
  )
  
}

population_ibd_barcodes_c <- function(barcode_vec,bl,nl,ib){
  
  n.strains <- lapply(barcode_vec,length) %>% unlist()
  lapply(
    barcode_vec[which(n.strains>0)],
    function(x){
      lapply(as.raw(tail(x[[1]],1)),test_ibd_conversion,bl,nl,ib) %>% unlist
    }
  )
  
}


population_ibd_distances <- function(r=NULL,nl=NULL,mat=NULL){
  
  if(is.null(mat)){
    mat <- unlist(population_ibd_barcodes(r,nl)) %>% matrix(nrow=length(r$population_List$Infection_States),byrow=TRUE)
  }
  
  nl <- dim(mat)[2]
  
  return(list("p_ibd"=proxy::dist(mat,method = function(x,y) sum(x==y)/nl),
              "pop" = mat))
  
}

micro_prev_2_10 <- function(out){
  
  prevs <- lapply(out[1:(length(out)-1)], function(x) {
    infs <- c(x[[positions[i]]]$Summary$State %in% c("D","T","A"))
    kids <- c(x[[positions[i]]]$Summary$Age_Bin %in% levels(x[[positions[i]]]$Summary$Age_Bin)[3:5])
    return((sum(x[[positions[i]]]$Summary$N[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$Summary$N[kids],na.rm=TRUE)))
  }) %>% unlist
  
}