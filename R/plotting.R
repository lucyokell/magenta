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


Convert_Barcode_Vectors <- function(sim.save, sub_patents_included=TRUE, ibd = FALSE, nl=24 ){
  
  out <- sim.save$Strain_barcode_vectors
  n.strains <- lapply(out,length) %>% unlist()
  
  if(ibd){
    intout <- population_ibd_barcodes(out,nl = nl)
  } else {
    
    intout <- lapply(out,function(x){
      if(length(x)>0){
        lapply(x,function(y){
          return(bitsToInt(y,endian = "big"))
        })
      } else {
        return(NA)
      }
    }
    )
  }
  
  int.out <- lapply(lapply(intout,function(x){return(unlist(x))}),unlist)
  
  if(!ibd){
    if(!sub_patents_included){
      COI <- rep(0,length(int.out))
      for(i in 1:length(int.out)){
        COI[i] <- length(unique(int.out[[i]][which(sim.save$Strain_infection_state_vectors[[i]]!=3)]))
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
#' 
#' 
#' @export

Sample_COI <- function(sim.save,sample_size,age_densities,age_breaks=seq(0,90*365,2*365),reps){
  
  
  COI_out <- Convert_Barcode_Vectors(sim.save)
  
  grouped_ages_by_2 <- cut(sim.save$Ages,breaks = age_breaks,labels = 1:45)
  
  picks <- replicate(n = reps,sample(x = 1:45,size = sample_size,replace = T,prob = age_densities))
  
  ids <- picks
  
  for(i in 1:reps){
    
    tabled <- table(picks[,i])
    id <- c()
    
    for(j in names(tabled)){
      
      id <- c(id,sample(x = which(grouped_ages_by_2==j & sim.save$Infection_States %in% c(2,3,4)),size = tabled[j],replace = T))
      
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
#' @param plot Boolean to print the plot
#' @importFrom ggplot2 ggplot aes geom_smooth theme_bw geom_point xlim ylim
#' 
#' 
#' @export

COI_age_plot_sample_x <- function(Sample_COI_out,x,span=0.6,ylimmax=NULL,xlimmax=NULL){
  
  
  mround <- function(x,base){ 
    base*round(x/base) 
  } 
  
  
  ids <- Sample_COI_out$ids[,x]
  Ages <- Sample_COI_out$sim.save$Ages[ids]/365
  COI <- Sample_COI_out$COI[ids]
  df <- data.frame("Ages"=Ages,"COI"=COI)
  
  df$COI[df$COI>25] <- sample(df$COI[df$COI<25])
  df$Ages <- mround(df$Ages,1)
  
  if(is.null(ylimmax)) ylimmax <- max(df$COI) + 1
  if(is.null(xlimmax)) xlimmax <- max(df$Ages) + 1
  
  
  gg <- ggplot(df,aes(Ages,COI)) + geom_point(color="blue") + geom_smooth(span=span,color="red",se=F) + theme_bw(base_size = 22) +
    ylim(c(0,ylimmax)) + xlim(c(0,xlimmax))
  print(gg)
  
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
summary_data_frames_from_sims <- function(res_list, update_length = 30,years = 30,
                                          age_breaks = list(1:4,5:6,7:9),
                                          ibd = FALSE){
  
  
  # set up list for final saves as these need to be replaced by COIs
  reps <- length(res_list)
  final_save <- list()
  length(final_save) <- reps
  
  # grab the final time point and convert them to the rest of the time points
  time_steps <- length(res_list[[1]])-1
  final_save <- lapply(res_list, function(x) {return(x[[time_steps]])})
  
  # create vector of correspondng times
  full_time <- seq(update_length,years*365,update_length)/365 + 2015 - years
  full_time <- full_time[-length(full_time)]
  
  # grab last 20 years of times
  time <- tail(full_time,floor(365*40/30))
  
  # positions of these times
  full_time_length <- 1:length(full_time)
  positions <- tail(full_time_length,length(time))
  times <- sort(rep(time,length(res_list)))
  
  # summary function
  mean_a_at_states_and_ages_in_c <- function(res_list,a,states=c("D","T","A","U"),ages=NULL,c,i){
    if(is.null(ages)) {
      ages <- seq_len(length(unique(res_list[[1]][[1]][[1]]$Age_Bin)))
    }
    lapply(res_list, function(x) {
      infs <- c(x[[positions[i]]][[c]]$State %in% states)
      age_pos <- c(x[[positions[i]]][[c]]$Age_Bin %in% levels(x[[positions[i]]][[c]]$Age_Bin)[ages])
      return(sum(x[[positions[i]]][[c]]$N[infs & age_pos]*x[[positions[i]]][[c]][[a]][infs & age_pos],na.rm=TRUE)/sum(x[[positions[i]]][[c]]$N[infs & age_pos],na.rm=TRUE))
    }) %>% unlist
    }
    
  mean_over_time <- function(res_list, a, states=c("D","T","A","U"),age_breaks=NULL,c,times){
    if(is.null(age_breaks)){
      return(unlist(lapply(times,function(y){mean_a_at_states_and_ages_in_c(res_list,a = a,c=c, i = y)})))
    } else {
      return(lapply(age_breaks, function(x){
        i = 1:length(time); 
        unlist(lapply(i,function(y){
          mean_a_at_states_and_ages_in_c(res_list,a = a, c=c, ages=x, i = y)
          }))
        }))
    }
  }
  
  
  if(!ibd) {
    
    grp <- rep(1:length(res_list),length(time))
    times <- seq_len(length(time))
    
    list_res <- list()
    
    list_res$meanCOI <- mean_over_time(res_list, a="mean",c="summary_coi",times=times)
    list_res$meanCOI_ages <- mean_over_time(res_list, a="mean",c="summary_coi",times=times,age_breaks=age_breaks)
    list_res$meanCOI_clinical <- mean_over_time(res_list, a="mean",c="summary_coi",times=times,states=c("D","T"))
    list_res$meanCOI_asymptomatic <- mean_over_time(res_list, a="mean",c="summary_coi",times=times,states=c("A"))
    
    list_res$mean_unique <- mean_over_time(res_list, a="unique",c="summary_polygenom",times=times)
    list_res$mean_unique_ages <- mean_over_time(res_list, a="unique",c="summary_polygenom",times=times,age_breaks=age_breaks)
    list_res$mean_unique_clinical <- mean_over_time(res_list, a="unique",c="summary_polygenom",times=times,states=c("D","T"))
    list_res$mean_unique_asymptomatic <- mean_over_time(res_list, a="unique",c="summary_polygenom",times=times,states=c("A"))
    
    list_res$mean_polygenomic <- mean_over_time(res_list, a="mean",c="summary_polygenom",times=times)
    list_res$mean_polygenomic_ages <- mean_over_time(res_list, a="mean",c="summary_polygenom",times=times,age_breaks=age_breaks)
    list_res$mean_polygenomic_clinical <- mean_over_time(res_list, a="mean",c="summary_polygenom",times=times,states=c("D","T"))
    list_res$mean_polygenomic_asymptomatic <- mean_over_time(res_list, a="mean",c="summary_polygenom",times=times,states=c("A"))
    
      
      meanclonality[1:reps + ((i-1)*reps)] <- lapply(res_list, function(x) {
        prods <- x[[positions[i]]]$clonality * as.numeric(names(x[[positions[i]]]$clonality))
        if(is.element(el = "1",names(x[[positions[i]]]$clonality))){
          return(prods[1]/sum(prods,na.rm=TRUE))
        } else {
          return(0)
        }
      }) %>% unlist
      
      prevs[1:reps + ((i-1)*reps)] <- lapply(res_list, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[3:5])
        return((sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[kids],na.rm=TRUE)))
      }) %>% unlist
      
    
    df <- data.frame("Time" = times, "meanCOI" = meanCOI, "ciCOI" = ciCOI, "sdCOI" = sdCOI, "clonality" = meanclonality,"Prev"=prevs,
                     "meanCOI_0_2"=meanCOI_0_2, "meanCOI_2_5"=meanCOI_2_5, "meanCOI_5_10"=meanCOI_5_10, "meanCOI_10_20" = meanCOI_10_20,
                     "meanCOI_20_100" =meanCOI_20_100)
    
    
    df_mean_COI <- summarySE(df,measurevar = "meanCOI",groupvars = "Time")
    df_mean_sdCOI <- summarySE(df,measurevar = "sdCOI",groupvars = "Time")
    df_mean_ciCOI <- summarySE(df,measurevar = "ciCOI",groupvars = "Time")
    df_mean_clonality <- summarySE(df,measurevar = "clonality",groupvars = "Time")
    df_mean_COI$ciCOI <- df_mean_ciCOI$ciCOI
    df_mean_micro_prev <- summarySE(df,measurevar =  "Prev",groupvars = "Time")
    age_bands <- c("meanCOI_0_2","meanCOI_2_5","meanCOI_5_10","meanCOI_10_20","meanCOI_20_100")
    df_mean_COI_ages <- lapply(age_bands,function(x){return(summarySE(df,measurevar = x,groupvars = "Time"))})
    df$rep <- grp
    
    return(list("df"=df,"COI_df"=df_mean_COI,"clonality_Df"=df_mean_clonality,"Prev_df"=df_mean_micro_prev,"COI_ages"=df_mean_COI_ages))
    
  } else {
    
    meanibd <- rep(0,length(times))
    meanibd_0_2 <- meanibd_2_5 <- meanibd_5_10 <- meanibd_10_20 <- meanibd_20_100 <- rep(0,length(times))
    meanclonality <- rep(0,length(times))
    ciibd <- rep(0,length(times))
    sdibd <- rep(0,length(times))
    prevs <- rep(0,length(times))
    grp <- rep(1:length(res_list),length(time))
    
    for(i in 1:(length(time))){
      
      meanibd[1:reps + ((i-1)*reps)] <- lapply(res_list, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        return(sum(x[[positions[i]]]$summary_coi$N[infs]*x[[positions[i]]]$summary_coi$pibd[infs],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs],na.rm=TRUE))
      }) %>% unlist
      
      meanibd_0_2[1:reps + ((i-1)*reps)] <- lapply(res_list, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[1:2])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$pibd[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      meanibd_2_5[1:reps + ((i-1)*reps)] <- lapply(res_list, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[3:4])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$pibd[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      meanibd_5_10[1:reps + ((i-1)*reps)] <- lapply(res_list, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[5])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$pibd[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      meanibd_10_20[1:reps + ((i-1)*reps)] <- lapply(res_list, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[6])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$pibd[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      meanibd_20_100[1:reps + ((i-1)*reps)] <- lapply(res_list, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[7:9])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$pibd[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      ciibd[1:reps + ((i-1)*reps)] <- lapply(res_list, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        return(sum(x[[positions[i]]]$summary_coi$N[infs]*x[[positions[i]]]$summary_coi$ci[infs],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs],na.rm=TRUE))
      }) %>% unlist
      
      sdibd[1:reps + ((i-1)*reps)] <- lapply(res_list, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        return(sum(x[[positions[i]]]$summary_coi$N[infs]*x[[positions[i]]]$summary_coi$sd[infs],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs],na.rm=TRUE))
      }) %>% unlist
      
      prevs[1:reps + ((i-1)*reps)] <- lapply(res_list, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[3:5])
        return((sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[kids],na.rm=TRUE)))
      }) %>% unlist
      
    }
    
    df <- data.frame("Time" = times, "meanibd" = meanibd, "ciibd" = ciibd, "sdibd" = sdibd,"Prev"=prevs,
                     "meanibd_0_2"=meanibd_0_2, "meanibd_2_5"=meanibd_2_5, "meanibd_5_10"=meanibd_5_10, "meanibd_10_20" = meanibd_10_20,
                     "meanibd_20_100" =meanibd_20_100)
    
    
    df_mean_ibd <- summarySE(df,measurevar = "meanibd",groupvars = "Time")
    df_mean_sdibd <- summarySE(df,measurevar = "sdibd",groupvars = "Time")
    df_mean_ciibd <- summarySE(df,measurevar = "ciibd",groupvars = "Time")
    df_mean_ibd$ciibd <- df_mean_ciibd$ciibd
    df_mean_micro_prev <- summarySE(df,measurevar =  "Prev",groupvars = "Time")
    age_bands <- c("meanibd_0_2","meanibd_2_5","meanibd_5_10","meanibd_10_20","meanibd_20_100")
    df_mean_ibd_ages <- lapply(age_bands,function(x){return(summarySE(df,measurevar = x,groupvars = "Time"))})
    df$rep <- grp
    
    return(list("df"=df,"ibd_df"=df_mean_ibd,"Prev_df"=df_mean_micro_prev,"ibd_ages"=df_mean_ibd_ages))
    
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
                        res[[i]]$Ages < age_bin[2] )/n
    }
    
  } else {
    for(i in 1:(l-1)){
      prev[i] <-  res[[i]][states] %>% unlist %>% sum
    }
  }
  plot(prev,type="l")
  invisible(prev)
  
}


plot_x_loggers <- function(res, x, alpha = 1, scale = FALSE) {
  
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
  
  melt <- reshape2::melt(df,id.var = "time")
  
  gg <- ggplot2::ggplot(melt,ggplot2::aes(x=time,y=value,color=variable)) + ggplot2::geom_point(alpha=alpha) + ggplot2::geom_line(alpha=alpha)
  print(gg)
  invisible(list("df"=melt,"gg"=gg))
  
}

# 
# ggplot(meanplot1, aes(x=value, y=measurement, color=treatment)) + 
#   geom_line(aes(group=sample), alpha=0.3) + 
#   stat_summary(aes(group = treatment), fun.y = mean, geom = 'line', size=3, alpha=0.9) +
#   theme_bw()