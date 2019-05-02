
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



clonality_from_barcode_list <- function(barcode_list){
  
  tbl <- table(table(unlist(lapply(barcode_list,unique))))
  
  if ("1" %in% names(tbl)){
    return(tbl[1]/sum(tbl*as.numeric(names(tbl))))
  } else {
    return(0)
  }
  
}

cou_from_barcode_list <- function(barcode_list){
  
  tbl <- table(unlist(lapply(barcode_list,unique)))
  samp_size <- sum(tbl)
  z <- sum((tbl/samp_size)^2)
  return((z - (1/samp_size))/(1-(1/samp_size)))
  
}

pibd_from_barcode_list <- function(barcode_list){
  l_factor_i <- attr(r[[year0]],"l") + attr(r[[year]],"l") 
  if(!unphased){
    t <- rbind(data.table::rbindlist(r[[year0]][["ints"]][samp[pos<=ss]]),data.table::rbindlist(r[[year]][["ints"]][samp[pos>ss]]))
    z <- mean(apply(t,2,function(x) sum(tab_func(x,l_factor_i)^2)))
    return((z - (1/ss))/(1-(1/ss)))
  } else {
    t <- rbind(data.table::rbindlist(r[[year0]][["ints_max"]][samp[pos<=ss]]),data.table::rbindlist(r[[year]][["ints_max"]][samp[pos>ss]]))
    z <- mean(apply(t,2,function(x) sum(tab_func(x,l_factor_i)^2)))
    return((z - (1/ss))/(1-(1/ss)))
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