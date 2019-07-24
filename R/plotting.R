#---
#' Plot ordered heatmap 
#'
#' \code{Heatmap_ordered_binary_plot} creates a heatmap of all human barcodes, 
#' with the most common barcode grouped at the bottom
#' 
#' @param sim.save Saved output of simulation
#' @param years Lenth of simulation
#' @param EIR Numeric for simulation EIR
#' @param ordered Boolean for oredring heatmap. Default = TRUE
#' @param save_path File path where tiff is saved to. Default = NULL
#' 
#' \code{Heatmap_ordered_binary_plot}
#' 
#' @export


Heatmap_ordered_binary_plot <- function(sim.save, years, EIR, ordered = TRUE,save_path=NULL){
  
  mat <- matrix(as.numeric(unlist(sim.save$populations_event_and_strains_List$Strain_barcode_vectors)),ncol=24,byrow=T)
  
  if(ordered){
    
    ID <- apply(mat,MARGIN = 1,FUN = bitsToInt)
    tabled <- sort(table(ID),decreasing = TRUE)
    sorted_row_pos <- unlist(sapply(names(tabled),function(x){return(which(ID==as.numeric(x)))}) )
    
    if(is.null(save_path)){
      dev.new(noRStudioGD = TRUE)
      heatmap(mat[sorted_row_pos,],Rowv = NA,Colv = NA,main = paste0( years," years | EIR = ",EIR),scale="none",labRow = "")
    } else {
      tiff(save_path,width=620,height=620,units="mm",res=300,pointsize = 36, compression="lzw",family="Times New Roman")
      heatmap(mat[sorted_row_pos,],Rowv = NA,Colv = NA,main = paste0( years," years | EIR = ",EIR),scale="none",labRow = "")
      dev.off()
    }
  } else {
    if(is.null(save_path)){
      dev.new(noRStudioGD = TRUE)
      heatmap(mat,Rowv = NA,Colv = NA,main = paste0( years," years | EIR = ",EIR),scale="none",labRow = "")  
    } else {
      tiff(save_path,width=620,height=620,units="mm",res=300,pointsize = 36, compression="lzw",family="Times New Roman")
      heatmap(mat,Rowv = NA,Colv = NA,main = paste0( years," years | EIR = ",EIR),scale="none",labRow = "")  
      dev.off()
    }
  }
  
}


#---
#' Convert human barcodes to numerics
#'
#' \code{Convert_Barcode_Vectors} converts human barcode vectors to nums
#' 
#' @param sim.save Saved output of simulation
#' @param sub_patents_included Boolean as to whether subpatents are included. Default = TRUE
#' 
#' \code{Convert_Barcode_Vectors}
#' 
#' @export


Convert_Barcode_Vectors <- function(sim.save, ID, sub_patents_included=TRUE, ibd = FALSE, nl=24, COI_type = "old"){
  
  
  
  # function to calculate probability of strain being detected by microscopy
  q_fun <- function(d1, ID, ID0, kD, fd) {
    return(d1 + ((1-d1) / (1 + ((ID/ID0)^kD)*fd)))
  }
  
  fd <- function(age, fD0, aD, gammaD) {
    return( 1 - ((1-fD0) / (1 + (age/aD)^gammaD)) )
  }
  
  ages <- sim.save$Ages
  ages[ages==0] <- 0.001
  mpl <- model_param_list_create()
  micro_det <- q_fun(mpl$d1,ID,mpl$ID0,mpl$kD,sapply(ages,fd,mpl$fD0,mpl$aD,mpl$gammaD)) 
  
  
  
  out <- sim.save$Strain_barcode_vectors
  n.strains <- lapply(out,length) %>% unlist()
  
  if(ibd){
    intout <- population_ibd_barcodes(out,nl = nl)
  } else {
    
    intout <- lapply(out,function(x){
      if(length(x)>0){
        lapply(x,function(y){
          return(bitsToInt(y))
        })
      } else {
        return(NA)
      }
    }
    )
  }
  
  int.out <- lapply(lapply(intout,function(x){return(unlist(x))}),unlist)
  
  COI_old <- function(i){
   a <- length(unique(int.out[[i]][which(sim.save$Strain_infection_state_vectors[[i]]==1)]))
    a <- a + round(length(unique(int.out[[i]][which(sim.save$Strain_infection_state_vectors[[i]]==2)]))*(micro_det[i]^mpl$alphaA))
    a <- a + round(length(unique(int.out[[i]][which(sim.save$Strain_infection_state_vectors[[i]]==3)]))*(micro_det[i]^mpl$alphaU))
    a <- a + length(unique(int.out[[i]][which(sim.save$Strain_infection_state_vectors[[i]]==4)]))
    return(a)
  }
  
  COI_new <- function(i){
    st <- sim.save$Strain_infection_state_vectors[[i]]
    stch <- sim.save$Strain_day_of_infection_state_change_vectors[[i]]
    staying <- st==1 | st == 4 | st == 2
    staying[st==3] <- as.logical(rbinom(sum(st==3),size = 1,micro_det[i]^mpl$alphaU))
    a <- length(unique(int.out[[i]][staying]))
    return(a)
  }
  
  COI_new_patent <- function(i){
    st <- sim.save$Strain_infection_state_vectors[[i]]
    stch <- sim.save$Strain_day_of_infection_state_change_vectors[[i]]
    staying <- st==1 | st == 4 
    staying[st==2] <- as.logical(rbinom(sum(st==2),size = 1,micro_det[i]))
    a <- length(unique(int.out[[i]][staying]))
    return(a)
  }
  
  
  
  if(!ibd){
    if(!sub_patents_included){
      COI <- rep(0,length(int.out))
      for(i in 1:length(int.out)){
        if(n.strains[i]!=0){
          
          if(COI_type=="old"){
        COI[i] <- COI_old(i)
          } else if (COI_type == "new"){
            COI[i] <- COI_new(i)
          } else {
            COI[i] <- COI_new_patent(i)
          }
        }
      }
    } else {
      COI <- unlist(lapply(int.out,function(x){return(length(unique(x)))}))
    }
  } else {
    COI <- rep(1,length(int.out))
  }
  COI[which(n.strains==0)] <- 0
  
  return(list("nums"=int.out,"COI"=COI))
  
}

#---
#' Sample x times from saved simulation and produce age COI data
#'
#' \code{Sample_COI} samples from a human save according to some age distirbution
#' 
#' @param sim.save Saved output of simulation
#' @param sample_size Numeric for sample size
#' @param age_densities Densitvy vector for age distribution required
#' @param age_breaks Corresponding vector of age breaks for the density
#' @param reps How many samples are made
#' @param COI_type type of COI calculation
#' 
#' 
#' @export

Sample_COI <- function(sim.save,ID,sample_size,age_densities,age_breaks=seq(0,90*365,2*365),
                       sub_patents_included=FALSE,reps=50, COI_type="old"){
  
  
  COI_out <- Convert_Barcode_Vectors(sim.save, ID, sub_patents_included = sub_patents_included,COI_type=COI_type)
  
  grouped_ages_by_2 <- cut(sim.save$Ages,breaks = age_breaks,labels = 1:45)
  sample_groups <- round(sample_size * (age_densities*(1/sum(age_densities)))) 
  picks <- matrix(rep(rep(1:45,sample_groups),reps),ncol=reps)
  
  ids <- picks
  
  for(i in 1:reps){
    
    tabled <- table(picks[,i])
    id <- c()
    
    for(j in names(tabled)){
      
      id <- c(id,sample(x = which(grouped_ages_by_2==j & sim.save$Infection_States %in% c(1,2,3,4)),size = tabled[j],replace = T))
      
    }
    
    ids[,i] <- id
    
  }
  
  return(list("ids"=ids,"COI"=COI_out$COI,"sim.save"=sim.save))
  
}

#---
#' Plot COI against age from random sample
#'
#' \code{Sample_COI} samples from a human save according to some age distirbution
#' 
#' @param Sample_COI_out Output of Sample_COI
#' @param x Which sample to plot
#' @param span Smoothing parameter for loess
#' @param ylimmax ylim max
#' @param xlimmax xlim max
#' @param max_coi max_coi
#' @importFrom ggplot2 ggplot aes geom_smooth theme_bw geom_point xlim ylim
#' 
#' 
#' @export

COI_age_plot_sample_x <- function(Sample_COI_out,x,span=0.6,ylimmax=NULL,xlimmax=NULL,max_coi=25){
  
  
  mround <- function(x,base){ 
    base*round(x/base) 
  } 
  
  
  ids <- Sample_COI_out$ids[,x]
  Ages <- Sample_COI_out$sim.save$Ages[ids]/365
  COI <- Sample_COI_out$COI[ids]
  df <- data.frame("Ages"=Ages,"COI"=COI)
  
  df$COI[df$COI>max_coi] <- sample(df$COI[df$COI<max_coi],sum(df$COI>max_coi))
  df$Ages <- mround(df$Ages,1)
  
  if(is.null(ylimmax)) ylimmax <- max(df$COI) + 1
  if(is.null(xlimmax)) xlimmax <- max(df$Ages) + 1
  
  
  gg <- ggplot(df,aes(Ages,COI)) + geom_point(alpha = 0.3, size=2, colour = "blue", fill="blue",shape=16) + geom_smooth(span=span,color="red",se=F) + theme_bw(base_size = 22) +
    ylim(c(0,ylimmax)) + xlim(c(0,xlimmax))
  #print(gg)
  
  invisible(gg)
  
  
}


#---
#' Create summary prev and genetic dataframes 
#' 
#' @param res_list list of simulation outputs
#' @param update_length Days between sim saves in simulation. Default = 30
#' @param years Years sims were run for
#' @param ibd Were we collecting info on IBD
#' 
summary_data_frames_from_sims <- function(res_list, update_length = 30,years = 35,
                                          age_breaks = list(1,2,3),coi_detect=FALSE,
                                          ibd = FALSE){
  
  
  if(coi_detect) {
    coi_det_var <- "summary_coi_detected"
  } else {
    coi_det_var <- "summary_coi"
  }
  
  age_breaks_all <- c("(-0.001,5]","(5,15]","(15,100]")
  states_all <- c("Clinical","Asymptomatic")
  # set up list for final saves as these need to be replaced by COIs
  reps <- length(res_list)
  
  # grab the final time point and convert them to the rest of the time points
  time_steps <- length(res_list[[1]])-1
  
  # create vector of correspondng times
  full_time <- seq(30,by = 30,to = time_steps*update_length)
  
  # positions of these times
  full_time_length <- length(full_time)
  positions <- 1:full_time_length
  time_vec <- sort(rep(positions,length(res_list)))
  
  # summary function
  mean_a_at_states_and_ages_in_c <- function(res_list,a,states=c("Asymptomatic","Clinical"),ages=NULL,c,i){
    if(is.null(ages)) {
      ages <- age_breaks_all
    } else if (is.na(ages)) {
    } else {
      ages <- age_breaks_all[ages]
    }
    if(is.null(states)) {
      states <- states_all
    } 
    lapply(res_list, function(x) {
      infs <- c(x[[positions[i]]][[c]]$clinical %in% states)
        age_pos <- c(x[[positions[i]]][[c]]$age_bin %in% ages)
      return(sum(x[[positions[i]]][[c]]$N[infs & age_pos]*x[[positions[i]]][[c]][[a]][infs & age_pos],na.rm=TRUE)/sum(x[[positions[i]]][[c]]$N[infs & age_pos],na.rm=TRUE))
    }) %>% unlist
  }
  
  mean_prev_ibd <- function(res_list,full_time_length){
    i = 1:(full_time_length); 
    lapply(i,function(y){
    lapply(res_list, function(x) {
      infs <- c(x[[positions[y]]][[1]]$clinical %in% c("Asymptomatic","Clinical"))
      return(sum(x[[positions[y]]][[1]]$N[infs],na.rm=TRUE)/sum(x[[positions[y]]][[1]]$N,na.rm=TRUE))
    }) %>% unlist
  }) %>% unlist
  }
  
  mean_prev <- function(res_list,full_time_length,prev="pcr_prev"){
    i = 1:(full_time_length); 
    lapply(i,function(y){
      lapply(res_list, function(x) {
        ((x[[positions[y]]][[prev]]))
      })
    }) %>% unlist
  }
  
  mean_over_time <- function(res_list, a, states=c("Asymptomatic","Clinical"),age_breaks=NULL,c,full_time_length){
    i = 1:full_time_length; 
    if(is.null(age_breaks) || is.na(age_breaks)){
      return(unlist(lapply(i,function(y){mean_a_at_states_and_ages_in_c(res_list,a = a,c=c, states=states, ages = age_breaks,i = y)})))
    } else {
      return(lapply(age_breaks, function(x){
        unlist(lapply(i,function(y){
          mean_a_at_states_and_ages_in_c(res_list,a = a, c=c,states=states, ages=x, i = y)
          }))
        }))
    }
  }
  
  
  
  list_res <- list()
  
  if(!ibd) {
    
    rept <- rep(1:length(res_list),length(full_time))
    
    list_res <- list()
    
    list_res$mean_polygenomic <- mean_over_time(res_list, a="mean",c="summary_polygenom",full_time_length=full_time_length)
    list_res$mean_polygenomic_ages <- mean_over_time(res_list, a="mean",c="summary_polygenom",full_time_length=full_time_length,age_breaks=age_breaks)
    list_res$mean_polygenomic_Clinical <- mean_over_time(res_list, a="mean",c="summary_polygenom",full_time_length=full_time_length,states=c("Clinical"))
    list_res$mean_polygenomic_Asymptomatic <- mean_over_time(res_list, a="mean",c="summary_polygenom",full_time_length=full_time_length,states=c("Asymptomatic"))
    
    list_res$meancoi <- mean_over_time(res_list, a="mean",c=coi_det_var,full_time_length=full_time_length)
    list_res$meancoi_ages <- mean_over_time(res_list, a="mean",c=coi_det_var,full_time_length=full_time_length,age_breaks=age_breaks)
    list_res$meancoi_Clinical <- mean_over_time(res_list, a="mean",c=coi_det_var,full_time_length=full_time_length,states=c("Clinical"))
    list_res$meancoi_Asymptomatic <- mean_over_time(res_list, a="mean",c=coi_det_var,full_time_length=full_time_length,states=c("Asymptomatic"))
    
    list_res$mean_unique <- mean_over_time(res_list, a="mean",c="summary_unique",full_time_length=full_time_length,states = NA,age_breaks = NA) %>% unlist
    list_res$mean_unique_ages <- mean_over_time(res_list, a="mean",c="summary_unique",full_time_length=full_time_length,age_breaks=age_breaks,states=NA)
    list_res$mean_unique_Clinical <- mean_over_time(res_list, a="mean",c="summary_unique",full_time_length=full_time_length,states=c("Clinical"),age_breaks=NA) %>% unlist
    list_res$mean_unique_Asymptomatic <- mean_over_time(res_list, a="mean",c="summary_unique",full_time_length=full_time_length,states=c("Asymptomatic"),age_breaks=NA) %>% unlist
    
    list_res$mean_cou <- mean_over_time(res_list, a="mean",c="summary_cou",full_time_length=full_time_length,states = NA,age_breaks = NA) %>% unlist
    list_res$mean_cou_ages <- mean_over_time(res_list, a="mean",c="summary_cou",full_time_length=full_time_length,age_breaks=age_breaks,states=NA)
    list_res$mean_cou_Clinical <- mean_over_time(res_list, a="mean",c="summary_cou",full_time_length=full_time_length,states=c("Clinical"),age_breaks=NA) %>% unlist
    list_res$mean_cou_Asymptomatic <- mean_over_time(res_list, a="mean",c="summary_cou",full_time_length=full_time_length,states=c("Asymptomatic"),age_breaks=NA) %>% unlist
    

    list_res$mean_prev <- mean_prev(res_list,full_time_length=full_time_length)
    
      df <- as.data.frame.list(list_res)
      df$rep <- rept
      df$time <- time_vec
      df$time <- df$time/12
      interventions <- factor(c("Low","Medium","High"), levels = c("Low","Medium","High"))
      df$intervention <- interventions[rep(c(rep(1,10),rep(2,10),rep(3,10)),length(df$time)/30)]
      
      names(df)[grep("ages", names(df))] <- c(paste0("mean_polygenomic_ages_",c("0-5","5-15","15+")),
                                              paste0("mean_coi_ages_",c("0-5","5-15","15+")),
                                              paste0("mean_unique_ages_",c("0-5","5-15","15+")),
                                              paste0("mean_cou_ages_",c("0-5","5-15","15+")))
      
      melted <- reshape2::melt(df, id.vars = c("time","rep","intervention"))
      melted$metric <- as.character(melted$variable)
      melted$metric[grepl("coi",melted$metric)] <- "COI"
      melted$metric[grepl("cou",melted$metric)] <- "COU"
      melted$metric[grepl("unique",melted$metric)] <- "% Unique"
      melted$metric[grepl("polygenom",melted$metric)] <- "% Polygenomic"
      melted$sampled <- "All"
      terms <- c("0-5","5-15","15+","Clinical","Asymptomatic")
      for(t in terms) {melted$sampled[grepl(t,melted$variable,fixed=TRUE)] <- t}
      levels(melted$sampled) <- c("All",terms)
      deeped <- dplyr::group_by(melted, intervention, metric, sampled, time) %>% summarise(value = mean(value, na.rm=TRUE))
      deeped$rep <- deeped$intervention
      
      return(list("melted" = melted, "deeped" = deeped))
      
  } else {
    
    rept <- rep(1:length(res_list),length(full_time))
    
    list_res$mean_prev <- mean_prev(res_list,full_time_length=full_time_length)
   
    list_res$mean_ibd <- mean_over_time(res_list, a="mean",c="summary_ibd",full_time_length=full_time_length,states = NULL,age_breaks = NULL) %>% unlist
    list_res$mean_ibd_ages <- mean_over_time(res_list, a="mean",c="summary_ibd",full_time_length=full_time_length,age_breaks=age_breaks,states=NULL)
    list_res$mean_ibd_Clinical <- mean_over_time(res_list, a="mean",c="summary_ibd",full_time_length=full_time_length,states=c("Clinical"),age_breaks=NULL) %>% unlist
    list_res$mean_ibd_Asymptomatic <- mean_over_time(res_list, a="mean",c="summary_ibd",full_time_length=full_time_length,states=c("Asymptomatic"),age_breaks=NULL) %>% unlist
    
    list_res$mean_within <- mean_over_time(res_list, a="mean",c="summary_within_ibd",full_time_length=full_time_length,states = NULL,age_breaks = NULL) %>% unlist
    list_res$mean_within_ages <- mean_over_time(res_list, a="mean",c="summary_within_ibd",full_time_length=full_time_length,age_breaks=age_breaks,states=NULL)
    list_res$mean_within_Clinical <- mean_over_time(res_list, a="mean",c="summary_within_ibd",full_time_length=full_time_length,states=c("Clinical"),age_breaks=NULL) %>% unlist
    list_res$mean_within_Asymptomatic <- mean_over_time(res_list, a="mean",c="summary_within_ibd",full_time_length=full_time_length,states=c("Asymptomatic"),age_breaks=NULL) %>% unlist
  
    df <- as.data.frame.list(list_res)
    df$rep <- rept
    df$time <- time_vec
    df$time <- df$time/12
    interventions <- factor(c("Low","Medium","High"), levels = c("Low","Medium","High"))
    df$intervention <- interventions[rep(c(rep(1,10),rep(2,10),rep(3,10)),length(df$time)/30)]
    
    names(df)[grep("ages", names(df))] <- c(paste0("mean_ibd_ages_",c("0-5","5-15","15+")),paste0("mean_within_ages_",c("0-5","5-15","15+")))
    
    melted <- reshape2::melt(df, id.vars = c("time","rep","intervention"))
    melted$metric <- as.character(melted$variable)
    melted$metric[grepl("ibd",melted$metric)] <- "pIBD"
    melted$metric[grepl("within",melted$metric)] <- "iIBD"
    melted$sampled <- "All"
    terms <- c("0-5","5-15","15+","Clinical","Asymptomatic")
    for(t in terms) {melted$sampled[grepl(t,melted$variable,fixed=TRUE)] <- t}
    levels(melted$sampled) <- c("All",terms)
    melted$value[which(melted$value==0)] <- 1
    deeped <- dplyr::group_by(melted, intervention, metric, sampled, time) %>% summarise(value = mean(value, na.rm=TRUE))
    deeped$rep <- deeped$intervention
    
    return(list("melted" = melted, "deeped" = deeped))
    
  }
  
}



plot_prevalence <- function(res, age_bin = NULL, states=c("D","A","U","T")){
  
  infs <- c("S","D","A","U","T","P")
  inf_states <- which(infs %in% states)
  
  l <- length(res)
  prev <- rep(0, l-1)
  
  if(is.null(age_bin)) {
    age_bin <- c(0,(365*100))
  } else {
    age_bin <- age_bin*365
  }
  
  if(sum(names(res[[1]]) %in% c("Ages","InfectionStates"))==2){
    n <- length(res[[1]]$InfectionStates)
    for(i in 1:(l-1)){
      prev[i] <-  sum(res[[i]]$InfectionStates %in% (inf_states-1) &
                        res[[i]]$Ages > age_bin[1] & 
                        res[[i]]$Ages < age_bin[2] )/
        sum(res[[i]]$Ages > age_bin[1] & 
              res[[i]]$Ages < age_bin[2])
    }
    
  } else {
    for(i in 1:(l-1)){
      prev[i] <-  res[[i]][states] %>% unlist %>% sum
    }
  }
  plot(prev,type="l")
  invisible(prev)
  
}


plot_x_loggers <- function(res, x, alpha = 1, scale = FALSE, extra = NULL ) {
  
  l <- list()
  length(l) <- length(x)
  names(l) <- x
  want <- 1:length(res)
  
  for(i in 1:length(x)) { 
    l[[i]] <- lapply(res[want],function(y) y[x[i]]) %>% unlist 
  }
  
  df <- as.data.frame.list(l)
  if(scale) {
  df <- t(apply(df,1,function(z) z/sum(z))) %>% as.data.frame()
  }
  df$time <- 1:length(l[[1]])
  if("summary_ibd" %in% names(res[[1]])){
  df$prev <- lapply(res[1:(length(res)-1)],function(x) sum(x$summary_ibd$N[x$summary_ibd$state %in% c("D","A","T","U")])) %>% unlist
  df$prev <- df$prev/sum(res[[1]]$summary_ibd$N)
  } else {
    df$prev <- lapply(res[1:(length(res)-1)],function(x) sum(x$summary_coi$N[x$summary_coi$state %in% c("D","A","T","U")])) %>% unlist
    df$prev <- df$prev/sum(res[[1]]$summary_coi$N)
  }
  melt <- reshape2::melt(df,id.var = "time")
  
  if(is.null(extra)){
  gg <- ggplot2::ggplot(melt,ggplot2::aes(x=time,y=value,color=variable)) + ggplot2::geom_point(alpha=alpha) + ggplot2::geom_line(alpha=alpha)
  } else {
    gg <- ggplot2::ggplot(melt,ggplot2::aes(x=time,y=value,color=variable)) + ggplot2::geom_point(alpha=alpha) + ggplot2::geom_line(alpha=alpha) + eval(parse(text = extra))
  }
  print(gg)
  invisible(list("df"=melt,"gg"=gg))
  
}

plot_mean_summary <- function(res, x, group = "clinical", groups = c("Asymptomatic","Clinical"),alpha = 1, scale = FALSE, extra = NULL ) {
  
  l <- list()
  length(l) <- length(x)
  names(l) <- x
  want <- 1:(length(res)-1)
  
  for(i in 1:length(x)) { 
    l[[i]] <- lapply(res[want],function(y){
      N <- y[[x]]$N[y[[x]][[group]] %in% groups]
      mean <- y[[x]]$mean[y[[x]][[group]] %in% groups]
      return(sum(N*mean)/sum(N))
    }) %>% unlist
  }
  
  df <- as.data.frame.list(l)
  if(scale) {
    df <- t(apply(df,1,function(z) z/sum(z))) %>% as.data.frame()
  }
  df$time <- 1:length(l[[1]])
  if("summary_ibd" %in% names(res[[1]])){
    df$prev <- lapply(res[1:(length(res)-1)],function(x) sum(x$summary_ibd$N[x$summary_ibd[[group]] %in% groups])) %>% unlist
    df$prev <- df$prev/sum(res[[1]]$summary_ibd$N)
  } else {
    df$prev <- lapply(res[1:(length(res)-1)],function(x) sum(x$summary_coi$N[x$summary_coi[[group]] %in% groups])) %>% unlist
    df$prev <- df$prev/sum(res[[1]]$summary_coi$N)
  }
  melt <- reshape2::melt(df,id.var = "time")
  
  if(is.null(extra)){
    gg <- ggplot2::ggplot(melt,ggplot2::aes(x=time,y=value,color=variable)) + ggplot2::geom_point(alpha=alpha) + ggplot2::geom_line(alpha=alpha)
  } else {
    gg <- ggplot2::ggplot(melt,ggplot2::aes(x=time,y=value,color=variable)) + ggplot2::geom_point(alpha=alpha) + ggplot2::geom_line(alpha=alpha) + eval(parse(text = extra))
  }
  print(gg)
  invisible(list("df"=melt,"gg"=gg))
  
}

# 

# 
# ggplot(meanplot1, aes(x=value, y=measurement, color=treatment)) + 
#   geom_line(aes(group=sample), alpha=0.3) + 
#   stat_summary(aes(group = treatment), fun.y = mean, geom = 'line', size=3, alpha=0.9) +
#   theme_bw()