


sim.out <- Pipeline(ft = 0.2,N=1000,years=2)

res <- list()
length(res) <- 30

ft <- c(rep(0.2,10),seq(0.2,0.8,length.out = 10),rep(0.8,10))

for ( i in 1 : 10){
  
  pl2 <- Param_List_Simulation_Update_Create(years = 2,ft=ft[i],statePtr = sim.out$Ptr)
  sim.out <- Simulation_R(pl2)
  pl3 <- Param_List_Simulation_Get_Create(sim.out$Ptr)
  res[[i]] <- Simulation_R(pl3)
  
}


bitsToInt<-function(x) {
  packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")
}

listed <- lapply(res,function(x){
  
  mat <- matrix(as.numeric(unlist(x$populations_event_and_strains_List$Strain_barcode_vectors)),ncol=24,byrow=T)
  ID <- apply(mat,MARGIN = 1,FUN = bitsToInt)
  return(as.character(ID))
  
}
)

all <- sort(table(unlist(listed)),decreasing = T)
labels <- paste0("s",1:sum(all > 1))
mapped <- cbind(labels,names(all)[1:sum(all>1)])

labelled_listed <- lapply(listed,function(x){
  
  out <- mapped[match(x,mapped[,2]),1]
  out[is.na(out)] <- "ind"
  return(out)
})

yrs <- rep(seq(2,60,length.out = 30),times=unlist(lapply(labelled_listed,length)))

samples <- c()
sample_size = 200
start <- 0
for(i in 1:length(res)){
  end <- (max(which(yrs==(2*i))))
  samples <- c(samples,sample((start+1):end,size = sample_size,replace=F))
  start <- end
}

df <- data.frame("ID"=unlist(labelled_listed),"year"=yrs, stringsAsFactors = F)
subset <- df[samples,]
## now loop through each year and change individual occurences to ind
for(i in 1:length(res)){
  year <- i*2
  subset$ID[(sample_size*(i-1))+match(names(which(table(subset$ID[subset$year==year])==1)),subset$ID[subset$year==year])] <- "ind"
  subset$ID[which(subset$year==year)] <- sort(subset$ID[which(subset$year==year)],decreasing = T)
  # sorted <- subset$ID[which(subset$year==i)]
  # subset$ID[which(subset$year==i)[sorted!="ind"]] <- rep(names(sort(table(sorted[which(sorted!="ind")]),decreasing=T)),
  #                                                        times=sort(table(sorted[which(sorted!="ind")]),decreasing=TRUE)) 
}


subset$ID <- factor(subset$ID, levels = (c("ind",names(sort(table(subset$ID),decreasing=T))[which(names(sort(table(subset$ID),decreasing=T))!="ind")])))
subset$year <- factor(subset$year)
# Function to generate ggplot colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

windows()

library(ggplot2)
ggplot(subset, aes(year, fill = ID)) +
  geom_bar(show.legend = F,colour="black") + 
  scale_fill_manual( values=c("#fafafb",gg_color_hue(length(names(sort(table(subset$ID),decreasing=T)))-1))) +
  guides(fill=guide_legend(reverse=F)) + theme_bw()

lapply(res,function(x){sum(x$population_List$Infection_States %in% c(1,2,3,4))})




#######################

par(mfrow=c(3,3))
windows()
plot((seq(30,length(wh)*30,30)/365)+2015-20,lapply(wh,function(x){return((sum(x$Summary$N*x$Summary$COI)/sum(x$Summary$N))*mean(lin_Triplet+0.2))}) %>% unlist,
     type="l",xlab = "Years",ylab="COI")
for(i in 1:8){
  plot((seq(30,length(wh)*30,30)/365)+2015-40,
       lapply(wh,function(x){
         bins <- x$Summary$Age_Bin==levels(x$Summary$Age_Bin)[i]
         return((sum(x$Summary$N[bins]*x$Summary$COI[bins])/sum(x$Summary$N[bins]))*lin_Triplet[i])
         }) %>% unlist,type="l",main = levels(x$Summary$Age_Bin)[i],xlab = "Years",ylab="COI")
}


prev <- lapply(wh[100:(length(wh)-1)],function(x) {return(sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U") & x$Summary$Age_Bin %in% levels(x$Summary$Age_Bin)[3:5]])/
                                                            sum(x$Summary$N[x$Summary$Age_Bin %in% levels(x$Summary$Age_Bin)[3:5]]))}) %>% unlist
COI <- lapply(wh[100:(length(wh)-1)],function(x){return((sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U")]*x$Summary$COI[x$Summary$State %in% c("D","T","A","U")])/sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U")])))}) %>% unlist
clons <- lapply(wh[100:(length(wh)-1)],function(x) {return(x$Clonality[1]/sum(x$Clonality))}) %>% unlist
actual_prev <- MAP[MAP$Name=="Couffo",5:20]

all_df <- data.frame("Prevalence"=prev,"COI"=COI,"Clonality"=clons)

library(ggplot2)
ggplot(all_df,aes(x=Prevalence,y=COI,colour=Clonality)) + geom_point() + geom_smooth(method = "lm")

 
glm <- glm(prev~COI+clons)
prev_predictons <- predict.lm(lm)



####


















df <- data.frame("COI"=unlist(COIs),"Age"=unlist(ages),"State"=infection_state[unlist(clinical_status)+1],"Years"=sort(rep(1:243,10000)))
df$Age <- df$Age/365
df <- dplyr::mutate(df,Age_Bin = cut(Age,breaks = c(0,1,3,5,10,20,40,60,100)))

library(dplyr)
library(magrittr)

ggplot(subset(df, State %in% c("D","A","U")),aes((Years),COI)) + geom_smooth(method="loess",span=0.4)


  geom_point(data = subset(df, State %in% c("D","A","U") & Age<5)[sample(subset(df, State %in% c("D","A","U")) %>% dim %>% extract(1),size = 3000,replace=FALSE),],aes(Years,COI)) +
  facet_wrap(~Age_Bin,scales = "free_y")


ggplot(subset(df, State %in% c("D","A","U") & Age > 5 & Age <10),aes((Years),COI)) + geom_smooth(method="loess",span=0.5) + facet_wrap(~State,scales = "free_y")

