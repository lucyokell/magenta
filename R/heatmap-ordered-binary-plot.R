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
  
  
  bitsToInt<-function(x) {
    packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")
  }
  
  mat <- matrix(as.numeric(unlist(sim.save$populations_event_and_strains_List$Strain_barcode_vectors)),ncol=24,byrow=T)
 
  if(ordered){
    
    ID <- apply(mat,MARGIN = 1,FUN = bitsToInt)
    tabled <- sort(table(ID),decreasing = TRUE)
    sorted_row_pos <- unlist(sapply(names(tabled),function(x){return(which(ID==as.numeric(x)))}) )
    
    if(is.null(save_path)){
      windows()
      heatmap(mat[sorted_row_pos,],Rowv = NA,Colv = NA,main = paste0( years," years | EIR = ",EIR),scale="none",labRow = "")
    } else {
      tiff(save_path,width=620,height=620,units="mm",res=300,pointsize = 36, compression="lzw",family="Times New Roman")
      heatmap(mat[sorted_row_pos,],Rowv = NA,Colv = NA,main = paste0( years," years | EIR = ",EIR),scale="none",labRow = "")
      dev.off()
    }
  } else {
    if(is.null(save_path)){
      windows()
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
#' 
#' \code{Convert_Barcode_Vectors}
#' 
#' @export


Convert_Barcode_Vectors <- function(sim.save){
  
  
  bitsToInt<-function(x) {
    packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")
  }
  
  out <- sim.save$Strain_barcode_vectors
  n.strains <- lapply(out,length) %>% unlist()
  
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
  
  int.out <- lapply(lapply(intout,function(x){return(unlist(x))}),unlist)
  COI <- unlist(lapply(int.out,function(x){return(length(unique(x)))}))
  COI[which(n.strains==0)] <- 0
  
  return(list("nums"=int.out,"COI"=COI))
  
}

#------------------------------------------------
#' Sample x times from saved simulation and produce age COI data
#'
#' \code{Sample_COI} samples from a human save according to some age distirbution
#' 
#' @param sim.save Saved output of simulation
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
  
  id <- c(id,sample(x = which(grouped_ages_by_2==j),size = tabled[j],replace = T))
  
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
#' 
#' 
#' @export

COI_age_plot_sample_x <- function(Sample_COI_out,x,span=span,ylimmax=NULL,xlimmax=NULL,plot=TRUE){
  
  
  mround <- function(x,base){ 
    base*round(x/base) 
  } 
  
  
  ids <- Sample_COI_out$ids[,x]
  df <- data.frame("Ages"=Sample_COI_out$sim.save$Ages[ids]/365,"COI"=Sample_COI_out$COI[ids])
  
  df$COI[df$COI>25] <- sample(df$COI[df$COI<25])
  df$Ages <- mround(df$Ages,1)
  
  if(is.null(ylimmax)){
    gg <- ggplot(df,aes(Ages,COI)) + geom_point(color="blue") + geom_smooth(span=span,color="red",se=F) + theme_bw(base_size = 22)
    if(plot){
      plot(gg)
    }
  } else {
    gg <- ggplot(df,aes(Ages,COI)) + geom_point(color="blue") + geom_smooth(span=span,color="red",se=F) + theme_bw(base_size = 22) +
    ylim(c(0,ylimmax)) + xlim(c(0,xlimmax))
    if(plot){
      plot(gg)
    }
  }
  
  if(!plot){
    return(gg)
  }
  
}


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
