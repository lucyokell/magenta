
###
#' Creates data frame containing the strain genetics of the population using 
#' the universe ptr. This data.frame is then passed to \code{COI_df_create}. 
#'
#' @param statePtr statePtr
#' @param sample_size sample_size
#' @param sample_states sample_states
#' @param ibd ibd boolean.
#' 
pop_strains_df <- function(statePtr, sample_size = 0, sample_states =  0:5, ibd = FALSE){
  
  paramList <- list("statePtr" = statePtr,
                    "sample_size" = sample_size,
                    "sample_states" = sample_states)
  
  if(!ibd){
    
    # Get the info
    list <- population_get_genetics_df_n(paramList)
    df <- as.data.frame.list(list[1:5])
    df$nums<- list$barcodes
  } else {
    
    # Get the info
    list <- population_get_genetics_ibd_df_n(paramList)
    df <- as.data.frame.list(list)
  }
  
  return(df)

}


###
#---
#' Create COI dataframe
#'
#' @param sim_save Output of a simulation save
#' @param groupvars Grouping vars for summarySE. Default = c("Age_Bin","State")
#' @param barcodes Boolean whether to return tabled barcodes. Default = FALSE
#' 
COI_df_create2 <- function(df, groupvars = c("age_bin","state"),
                           barcodes=FALSE, ibd = 0, nl = 24, n = Inf, reps = 1){
  
  
  # handle df for either genetic type first
  infection_state <- c("S","D","A","U","T","P")
  df$last_treatment[df$last_treatment == 0] <- 0.001
  df$age_bin = cut(df$age/365,breaks = c(0,1,2,3,5,10,20,40,60,100))
  df$last_treatment_binned = cut(df$last_treatment,breaks = c(0,28,90,365,Inf))
  df$state <- infection_state[df$state + 1]
  
  # if doing samples
  if(n != Inf){
    # if we have asked for more samples than there are then sample from them
    if((reps*n) > nrow(df)) {
      if(n > nrow(df)) {
      ranges <- lapply(seq_len(reps), function(x) sample(nrow(df), n, TRUE))
      } else {
      ranges <- lapply(seq_len(reps), function(x) sample(nrow(df), n, FALSE))
      }
    } else {
      sample <- sample(reps * n)
      ranges <- ranges(diff = n, end = reps * n)
      ranges <- lapply(ranges, function(x) sample[x])
    }
  }
  
  
  if(!ibd){
    
    # polygenomic return
    df$polygenom <- df$coi>1
    
    # set up results
    results <- list()
    length(results) <- reps
    
    # make samples and reps of them
    for(i in seq_len(length(results))) {
      
      # make a sample
      if(n == Inf) {
        chosen <- seq_len(nrow(df))
      } else {
        chosen <- ranges[[i]]
      }
      
      # create overall clonality
      clonality <- table(table(unlist(df$nums[chosen])))
      barcodes_tab <- sort(table(unlist(df$nums[chosen])),decreasing=TRUE)
      
      # create summary dfs
      summary_coi <- summarySE_mean_only(df[chosen,],measurevar = "coi", groupvars = groupvars)
      summary_coi_detected <- summarySE_mean_only(df[chosen,],measurevar = "coi_detected_micro", groupvars = groupvars)
      summary_polygenom <- summarySE_mean_only(df[chosen,],measurevar = "polygenom", groupvars = groupvars)
      summary_polygenom <- dplyr::left_join(
        summary_polygenom,
        dplyr::group_by(df[chosen,], age_bin, state) %>% 
          dplyr::summarise(unique = clonality_from_barcode_list(nums)), 
        by = groupvars
      )
      
      # bundle it all up
      res <- list("summary_coi"=summary_coi,
                  "summary_coi_detected"=summary_coi_detected,
                  "summary_polygenom"=summary_polygenom,
                  "clonality"=clonality,"coi_table"=table(df$coi[chosen]))
      
      if(barcodes){
        res$barcodes <- barcodes_tab[barcodes_tab>1]
      } 
      
      results[[i]] <- res
    }
    
  } else {
    
    # set up results
    results <- list()
    length(results) <- reps
    
    # make samples and reps of them
    for(i in seq_len(length(results))) {
      
      # make a sample
      if(n == Inf) {
        chosen <- seq_len(nrow(df))
      } else {
        chosen <- ranges[[i]]
      }
      
      summary_ibd <- summarySE_mean_only(df[chosen,],measurevar = "pibd",groupvars = groupvars)
      summary_ibd_d <- summarySE_mean_only(df[chosen,],measurevar = "pibd_d",groupvars = groupvars)
      summary_within_ibd <- summarySE_mean_only(df[chosen,],measurevar = "pibd_within",groupvars = groupvars)
      summary_within_ibd_d <- summarySE_mean_only(df[chosen,],measurevar = "pibd_within_d",groupvars = groupvars)
      res <- list("summary_ibd"=summary_ibd,
                  "summary_ibd_d"=summary_ibd_d,
                  "summary_within_ibd"=summary_within_ibd,
                  "summary_within_ibd_d"=summary_within_ibd_d,
                  "mean_ibd"=mean(df$pibd[chosen], na.rm = TRUE),
                  "mean_ibd_d"=mean(df$pibd_d[chosen], na.rm = TRUE),
                  "mean_within"=mean(df$pibd_within[chosen], na.rm = TRUE),
                  "mean_within_d"=mean(df$pibd_within_d[chosen], na.rm = TRUE))
                  
      results[[i]] <- res
    }
  }
  
  # if no reps then unlist
  if(reps == 1) {
    results <- res
  }
  return(results)
  
  
}