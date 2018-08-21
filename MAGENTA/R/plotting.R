#------------------------------------------------
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


#------------------------------------------------
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

#------------------------------------------------
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

#------------------------------------------------
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


#------------------------------------------------
#' Create summary prev and genetic dataframes 
#' 
#' @param old_res list of simulation outputs
#' @param update_length Days between sim saves in simulation. Default = 30
#' @param years Years sims were run for
#' @param ibd Were we collecting info on IBD
#' 
summary_data_frames_from_sims <- function(old_res, update_length = 30,years, ibd = FALSE){
  
  
  # set up list for final saves as these need to be replaced by COIs
  reps <- length(old_res)
  final_save <- list()
  length(final_save) <- reps
  
  # grab the final time point and convert them to the rest of the time points
  time_steps <- length(old_res[[1]])-1
  final_save <- lapply(old_res, function(x) {return(x[[time_steps]])})
  # for(i in 1:length(old_res)){
  #   # Convert the barcodes to COIs
  #   COIS <- Convert_Barcode_Vectors(final_save[[i]],sub_patents_included=TRUE)
  #   COIs <- COIS$COI
  #   COIs[is.na(COIs)] <- 0
  #   clonality <- table(table(unlist(COIS$nums)))
  #   
  #   # grab the age, status etc
  #   ages <- final_save[[i]]$Ages
  #   ages[ages==0] <- 0.001
  #   clinical_status <- final_save[[i]]$Infection_States
  #   infection_state <- c("S","D","A","U","T","P")
  #   
  #   #  bring into df
  #   df <- data.frame("COI"=COIs,"Ages"=ages/365,"State"=infection_state[(clinical_status)+1],
  #                    "Age_Bin" = cut(ages/365,breaks = c(0,1,3,5,10,20,40,60,100)),
  #                    stringsAsFactors = FALSE)
  #   
  #   # Add the non subpatent COIs
  #   COIS <- Convert_Barcode_Vectors(final_save[[i]],sub_patents_included=FALSE)
  #   COIs <- COIS$COI
  #   COIs[is.na(COIs)] <- 0
  #   df$COI_Detected_From_Model <- COIs
  #   
  #   summary <- summarySE(df,measurevar = "COI",groupvars = c("Age_Bin","State"))
  #   
  #   old_res[[i]][[time_steps]] <- list("Summary"=summary,"Clonality"=clonality)
  #   
  # }
  
  # create vector of correspondng times
  full_time <- seq(update_length,years*365,update_length)/365 + 2015 - years
  full_time <- full_time[-length(full_time)]
  
  # grab last 20 years of times
  time <- tail(full_time,floor(365*40/30))
  
  # positions of these times
  full_time_length <- 1:length(full_time)
  positions <- tail(full_time_length,length(time))
  
  times <- sort(rep(time,length(old_res)))
  
  if(!ibd) {
    meanCOI <- rep(0,length(times))
    meanCOI_0_2 <- meanCOI_2_5 <- meanCOI_5_10 <- meanCOI_10_20 <- meanCOI_20_100 <- rep(0,length(times))
    meanclonality <- rep(0,length(times))
    ciCOI <- rep(0,length(times))
    sdCOI <- rep(0,length(times))
    prevs <- rep(0,length(times))
    grp <- rep(1:length(old_res),length(time))
    
    for(i in 1:(length(time))){
      
      meanCOI[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        return(sum(x[[positions[i]]]$summary_coi$N[infs]*x[[positions[i]]]$summary_coi$mean[infs],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs],na.rm=TRUE))
      }) %>% unlist
      
      meanCOI_0_2[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[1:2])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$mean[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      meanCOI_2_5[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[3:4])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$mean[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      meanCOI_5_10[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[5])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$mean[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      meanCOI_10_20[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[6])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$mean[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      meanCOI_20_100[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[7:9])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$mean[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      ciCOI[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        return(sum(x[[positions[i]]]$summary_coi$N[infs]*x[[positions[i]]]$summary_coi$ci[infs],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs],na.rm=TRUE))
      }) %>% unlist
      
      sdCOI[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        return(sum(x[[positions[i]]]$summary_coi$N[infs]*x[[positions[i]]]$summary_coi$sd[infs],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs],na.rm=TRUE))
      }) %>% unlist
      
      meanclonality[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        prods <- x[[positions[i]]]$clonality * as.numeric(names(x[[positions[i]]]$clonality))
        if(is.element(el = "1",names(x[[positions[i]]]$clonality))){
          return(prods[1]/sum(prods,na.rm=TRUE))
        } else {
          return(0)
        }
      }) %>% unlist
      
      prevs[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[3:5])
        return((sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[kids],na.rm=TRUE)))
      }) %>% unlist
      
    }
    
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
    grp <- rep(1:length(old_res),length(time))
    
    for(i in 1:(length(time))){
      
      meanibd[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        return(sum(x[[positions[i]]]$summary_coi$N[infs]*x[[positions[i]]]$summary_coi$pibd[infs],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs],na.rm=TRUE))
      }) %>% unlist
      
      meanibd_0_2[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[1:2])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$pibd[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      meanibd_2_5[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[3:4])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$pibd[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      meanibd_5_10[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[5])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$pibd[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      meanibd_10_20[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[6])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$pibd[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      meanibd_20_100[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        kids <- c(x[[positions[i]]]$summary_coi$Age_Bin %in% levels(x[[positions[i]]]$summary_coi$Age_Bin)[7:9])
        return(sum(x[[positions[i]]]$summary_coi$N[infs & kids]*x[[positions[i]]]$summary_coi$pibd[infs & kids],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs & kids],na.rm=TRUE))
      }) %>% unlist
      
      ciibd[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        return(sum(x[[positions[i]]]$summary_coi$N[infs]*x[[positions[i]]]$summary_coi$ci[infs],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs],na.rm=TRUE))
      }) %>% unlist
      
      sdibd[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
        infs <- c(x[[positions[i]]]$summary_coi$State %in% c("D","T","A","U"))
        return(sum(x[[positions[i]]]$summary_coi$N[infs]*x[[positions[i]]]$summary_coi$sd[infs],na.rm=TRUE)/sum(x[[positions[i]]]$summary_coi$N[infs],na.rm=TRUE))
      }) %>% unlist
      
      prevs[1:reps + ((i-1)*reps)] <- lapply(old_res, function(x) {
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



plot_prevalence <- function(res, age_bin = NULL){
  
  l <- length(res)
  prev <- rep(0, l-1)
  
  if(is.null(age_bin)) {
    age_bin <- c(0,(365*100))
  } else {
    age_bin <- age_bin*365
  }
  
  if(length(names(res[[1]]))==6){
    n <- length(res[[1]]$Infection_States)
    for(i in 1:(l-1)){
      prev[i] <-  sum(res[[i]]$Infection_States %in% c(1,2,3,4) &
                        res[[i]]$Ages > age_bin[1] & 
                        res[[i]]$Ages < age_bin[2] )/n
    }
    
  } else {
    for(i in 1:(l-1)){
      prev[i] <-  res[[i]][c("D","A","U","T")] %>% unlist %>% sum
    }
  }
  plot(prev,type="l")
  invisible(prev)
  
}


plot_x_loggers <- function(res, x) {
  
  l <- list()
  length(l) <- length(x)
  names(l) <- x
  want <- 1:length(res)
  
  for(i in 1:length(x)) { 
    l[[i]] <- lapply(res[want],function(y) y[x[i]]) %>% unlist 
  }
  
  df <- as.data.frame.list(l)
  df <- t(apply(df,1,function(z) z/sum(z))) %>% as.data.frame()
  df$time <- 1:length(l[[1]])
  
  melt <- reshape2::melt(df,id.var = "time")
  
  gg <- ggplot(melt,aes(x=time,y=value,color=variable)) + geom_point() + geom_line()
  print(gg)
  invisible(list("df"=melt,"gg"=gg))
  
}


# 
# ggplot(df_mean_COI,aes(x=Time,y=meanCOI)) + geom_line() + geom_ribbon(aes(x=Time,ymin=meanCOI - ciCOI,ymax=meanCOI+ciCOI))
# 
# ggplot(df_mean_clonality,aes(x=Time,y=Clonality)) + geom_line() + geom_ribbon(aes(x=Time,ymin=Clonality - ci,ymax=Clonality+ci))

# 
# 
# windows()
# par(mfrow=c(5,5))
# for(i in 1:100){
# plot(tres$Ages[ids[COI[ids[,i]]<25,i]]/365,COI[ids[COI[ids[,i]]<25,i]],col="blue",pch=16,ylim=c(0,25),ylab = "COI",xlab="Age")
# lines(loess.smooth(tres$Ages[ids[COI[ids[,i]]<25,i]]/365,COI[ids[COI[ids[,i]]<25,i]],evaluation = 1000),col="red",lwd=3)
# }
# lines <- apply(outids$ids,MARGIN = 2,function(x){return(loess.smooth(tres2$Ages[x]/365,outids$COI[x],span=0.6))})
# plot(lines[[1]],ylim=c(0,10))
# lapply(lines,lines,col="red")
# 
# ## Extra function required to calculate discrete time mean and quantile estimations
# data_summary <- function(data, varname, grps, lower.quantile, upper.quantile){
#   
#   ## summary function describing the quantile calculation
#   summary_func <- function(x, col){
#     c(mean = mean(x[[col]], na.rm=TRUE),
#       sd = quantile(x[[col]], na.rm=TRUE,probs=c(lower.quantile,upper.quantile)))
#   }
#   
#   ## Calculate summary output
#   data_sum <- plyr::ddply(data, grps, .fun=summary_func, varname)
#   data_sum <- plyr::rename(data_sum, c("mean" = varname))
#   names(data_sum)[c(3,4)] <- c("min","max")
#   return(data_sum)
# }
# melted.lines <- reshape2::melt(lines,measure.vars=c("x","y"))
# summary <- data_summary(melted.lines,varname = "value",grps="L2",0.25,0.5)
# 
# best.line <- loess.smooth(tres2$Ages[outids$ids],outids$COI[outids$ids],span=0.6)
# plot(best.line)
# lo <- loess(melted.lines$value[melted.lines$L2=="y"]~melted.lines$value[melted.lines$L2=="x"],span = 10)
# plot(melted.lines$value[melted.lines$L2=="x"],melted.lines$value[melted.lines$L2=="y"])
# lines(predict(lo), col='red', lwd=2)
# plot(tres$Ages[ids],COI[ids],col="blue",ylim=c(0,50))
# 
# 
# xs <- lapply(lines,function(x){return(x$x)})
# xsmat <- matrix(unlist(xs),nrow=100,byrow=T)
# truex <- colMeans(xsmat)
# ys <- lapply(lines,function(x){return(x$y)})
# ysmat <- matrix(unlist(ys),nrow=100,byrow=T)
# truey <- colMeans(ysmat)
# plot(truex,truey)
# 
