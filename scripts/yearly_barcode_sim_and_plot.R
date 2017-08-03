


sim.out <- Pipeline(10,ft = 0.2,N=10000,years=10)

res <- list()
length(res) <- 30

ft <- c(rep(0.2,10),seq(0.2,0.8,length.out = 10),rep(0.8,10))

for ( i in 1 : 30){
  
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

ggplot(subset, aes(year, fill = ID)) +
  geom_bar(show.legend = F,colour="black") + 
  scale_fill_manual( values=c("#fafafb",gg_color_hue(length(names(sort(table(subset$ID),decreasing=T)))-1))) +
  guides(fill=guide_legend(reverse=F)) + theme_bw()

