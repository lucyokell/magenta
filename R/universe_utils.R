
###
#' Creates data frame containing the strain genetics of the population using 
#' the universe ptr. This data.frame is then passed to \code{COI_df_create}. 
#'
#' @param statePtr statePtr
#' @param sample_size sample_size
#' @param sample_states sample_states
#' @param ibd ibd boolean.
#' 
pop_strains_df <- function(statePtr, sample_size = 0, sample_states =  0:5, ibd = FALSE, seed, genetics_df_without_summarising = FALSE, nl = 24, big_mat_test = TRUE){
  
  paramList <- list("statePtr" = statePtr,
                    "sample_size" = sample_size,
                    "sample_states" = sample_states)
  #message("pop_strains\n")
  if(!ibd){
    
    # Get the info
    set.seed(seed)
    list <-  population_get_genetics_df_n(paramList)
    df <- as.data.frame.list(list[1:4])
  
  } else {
    
    # Get the info
    set.seed(seed)
    list <- population_get_genetics_ibd_df_n(paramList)
    df <- as.data.frame.list(list[1:6])
    
  }
  
  breaks <- cumsum(unlist(lapply(list$barcode_states,length)))
  seqs <- unlist(lapply(list$barcode_states,length))
  df$nums <- list$barcode_states
  df$barcode_states <- list$barcode_states
  if(seqs[1] > 0){
  df$nums[[1]] <- list$barcodes[1:breaks[1],, drop=FALSE]
  }
  
  for(i in 2:(length(breaks))) { 
    if(seqs[i] > 0){
    df$nums[[i]] <- list$barcodes[(breaks[i-1]+1):breaks[i],,drop=FALSE]
    }
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
COI_df_create2 <- function(df, groupvars = c("age_bin","clinical"),breaks = c(-0.001,5,15,100.1),
                           barcodes=FALSE, ibd = 0, nl = 24, n = Inf, reps = 1,
                           mean_only = TRUE){
  
  
  # handle df for either genetic type first
  infection_state <- c("S","D","A","U","T","P")
  df$age_bin = cut(df$age/365,breaks = c(-0.001,5,15,100.1))
  df$state <- df$clinical <- infection_state[df$state + 1]
  df$clinical[df$clinical %in% c("D", "T")] <- "Clinical"
  df$clinical[df$clinical %in% c("A", "U")] <- "Asymptomatic"
  df$bs <- unlist(lapply(df$barcode_states,length))
  
  nrow_df <- nrow(df)
  # if doing samples
  if(length(n) == 1) {
    if(n != Inf){
      # if we have asked for more samples than there are then sample from them
      if((reps*n) > nrow_df) {
        if(n > nrow_df) {
          ranges <- lapply(seq_len(reps), function(x) sample(nrow_df, n, TRUE))
        } else {
          ranges <- lapply(seq_len(reps), function(x) sample(nrow_df, n, FALSE))
        }
      } else {
        sample <- sample(reps * n)
        ranges <- ranges(diff = n, end = reps * n)
        ranges <- lapply(ranges, function(x) sample[x])
      }
    } 
  } else {
    reps <- length(n) + 1
    ranges <- lapply(n, function(x){
      if(x > nrow_df) {
        sample(nrow_df, x, TRUE)
      } else {
        sample(nrow_df, x, FALSE)
      }
    })
    ranges[[length(ranges)+1]] <- seq_len(nrow_df)
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
      if(n[1] == Inf) {
        chosen <- seq_len(nrow_df)
      } else {
        chosen <- ranges[[i]]
      }
      
      # create overall clonality by first converting our nums to integers for quick tabulation
      df$nums[df$bs>0]  <- lapply(df$nums[df$bs>0],function(x) apply(x,1,bitsToInt))
      clonality <- table(table(unlist(df$nums[chosen])))
      barcodes_tab <- sort(table(unlist(df$nums[chosen])),decreasing=TRUE)
      
      # create summary dfs
      summary_coi <- summarySE_mean_only(df[chosen,],measurevar = "coi", groupvars = groupvars, mean_only = mean_only)
      summary_coi_detected <- summarySE_mean_only(df[chosen,],measurevar = "coi_detected_micro", groupvars = groupvars, mean_only = mean_only)
      summary_moi_detected <- summarySE_mean_only_max_mean(df[chosen,],measurevar = "coi_detected_micro", groupvars = groupvars, mean_only = mean_only, max = 6)
      summary_polygenom <- summarySE_mean_only(df[chosen,],measurevar = "polygenom", groupvars = groupvars, mean_only = mean_only)
      summary_unique1 <- dplyr::group_by(df[chosen,], age_bin) %>% dplyr::summarise(N=sum(!is.na(nums)),age = mean(age[!is.na(nums)]), mean = clonality_from_barcode_list(nums))
      summary_unique2 <- dplyr::group_by(df[chosen,], clinical) %>% dplyr::summarise(N=sum(!is.na(nums)),age = mean(age[!is.na(nums)]),mean = clonality_from_barcode_list(nums))
      summary_unique3 <- dplyr::group_by(df[chosen,]) %>% dplyr::summarise(N=sum(!is.na(nums)),age = mean(age[!is.na(nums)]),mean = clonality_from_barcode_list(nums))
      summary_unique <- dplyr::bind_rows(list(summary_unique1,summary_unique2,summary_unique3))
      summary_ages <- summarySE_mean_only(df[chosen,],measurevar = "age", groupvars = groupvars, mean_only = mean_only)
      
      # bundle it all up
      res <- list("summary_coi"=summary_coi,
                  "summary_coi_detected"=summary_coi_detected,
                  "summary_moi_detected"=summary_moi_detected,
                  "summary_polygenom"=summary_polygenom,
                  "summary_ages"=summary_ages,
                  "summary_unique"=summary_unique,
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
      if(n[1] == Inf) {
        chosen <- seq_len(nrow_df)
      } else {
        chosen <- ranges[[i]]
      }
      
      summary_ibd <- summarySE_mean_only(df[chosen,],measurevar = "pibd", groupvars = groupvars,mean_only = mean_only)
      summary_ibd_d <- summarySE_mean_only(df[chosen,],measurevar = "pibd_d",groupvars = "age_bin", mean_only = mean_only)
      summary_within_ibd <- summarySE_mean_only(df[chosen,],measurevar = "pibd_within",groupvars = groupvars, mean_only = mean_only)
      summary_within_ibd_d <- summarySE_mean_only(df[chosen,],measurevar = "pibd_within_d",groupvars = groupvars, mean_only = mean_only)
      summary_ages <- summarySE_mean_only(df[chosen,],measurevar = "age", groupvars = groupvars, mean_only = mean_only)
      res <- list("summary_ibd"=summary_ibd,
                  "summary_ibd_d"=summary_ibd_d,
                  "summary_within_ibd"=summary_within_ibd,
                  "summary_within_ibd_d"=summary_within_ibd_d,
                  "summar_ages"=summary_ages,
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