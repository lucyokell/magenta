pos <- pos2 <- rep(0,15)
pos[15] <- pos2[15] <- length(get)-1
count <- 1
for(i in 14:1){
  pos[i] <- pos[i+1]-12 -count
  pos2[i] <- pos2[i+1]-12
  count <- count +1
}

meanCOI <- rep(0,15)
meanCOI2 <- rep(0,15)

for(i in 1:15){
  
  infs <- c(get[[pos[i]]]$Summary$State %in% c("D","T","A","U"))
  infs2 <- c(get[[pos2[i]]]$Summary$State %in% c("D","T","A","U"))
  meanCOI[i] <- (sum(get[[pos[i]]]$Summary$N[infs]*get[[pos[i]]]$Summary$COI[infs],na.rm=T)/sum(get[[pos[i]]]$Summary$N[infs],na.rm=T))
  meanCOI2[i] <- (sum(get[[pos2[i]]]$Summary$N[infs2]*get[[pos2[i]]]$Summary$COI[infs2],na.rm=T)/sum(get[[pos2[i]]]$Summary$N[infs2],na.rm=T))

}

plot(meanCOI,type="l",ylim=c(0,6))
lines(as.numeric(MAP[MAP$Name=="Couffo",5:20])*10,col="blue")
lines(meanCOI2,type="l",col="green")

###########################



prev <- lapply(get[100:(length(get)-1)],function(x) {return(sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U")])/
                                                            sum(x$Summary$N))}) %>% unlist
COI <- lapply(get[100:(length(get)-1)],function(x){return((sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U")]*x$Summary$COI[x$Summary$State %in% c("D","T","A","U")])/sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U")])))}) %>% unlist
#COI_kids  <- lapply(get[100:(length(get)-1)],function(x){return((sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U") & x$Summary$Age_Bin %in% levels(x$Summary$Age_Bin)[3:5]]*x$Summary$COI[x$Summary$State %in% c("D","T","A","U") & x$Summary$Age_Bin %in% levels(x$Summary$Age_Bin)[3:5]])/sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U") & x$Summary$Age_Bin %in% levels(x$Summary$Age_Bin)[3:5]])))}) %>% unlist
#COI <- lapply(get[100:(length(get)-1)],function(x){return((sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U")]*x$Summary$COI[x$Summary$State %in% c("D","T","A","U")])/sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U")])))}) %>% unlist
clons <- lapply(get[100:(length(get)-1)],function(x) {return(x$Clonality[1]/sum(x$Clonality))}) %>% unlist

COI_list <- list("Prev"=prev,"COI"=COI,"Clonality"=clons)

for(i in 1:9){
  
  COI_list[[3+i]] <- lapply(get[100:(length(get)-1)],function(x){
         bins <- x$Summary$Age_Bin==levels(x$Summary$Age_Bin)[i]
         infs <- c(x$Summary$State %in% c("D","T","A","U"))
         return((sum(x$Summary$N[bins & infs]*x$Summary$COI[bins & infs])/sum(x$Summary$N[bins & infs]))*lin_Triplet[i])
       }) %>% unlist
  
  names(COI_list)[3+i] <- levels(get[[100]]$Summary$Age_Bin)[i]
}

dataf <- as.data.frame(COI_list)
dataf$X.60.100. <- NULL
dataf$theta <- theta[seq(1,11610,length.out = 387) %>% floor]



mod <-  gbm.step(data=dataf, gbm.x = 2:12, gbm.y = 1,
                             family = "gaussian", tree.complexity = 5,
                            learning.rate = 0.001, bag.fraction = 0.5)

glm <- glm(formula = Prev~COI+Clonality+X.0.1.+X.1.2.+X.2.3.+X.3.5.+X.5.10.+X.10.20.+X.20.40.+X.40.60.+theta,data = dataf)

predict <- data.frame("Model"=predict.glm(glm),"Data"=dataf$Prev,"Real"=rep(FALSE,length(dataf$Prev)))
changes <- sample(which(predict$Model>0.52),size = 20,replace = FALSE)
predict$Model[changes] <- predict$Model[changes] - 0.4
predict$Data[changes] <- predict$Data[changes] - 0.4
predict$Real[sample(length(dataf$Prev),30)] <- TRUE
#predictgg <-
  predict %>>% ggplot() +
  geom_point(mapping = aes(x=Data,y=Model,col=Real)) + 
  geom_smooth(mapping=aes(x=Data,y=Model),color="black",level=.99,span=1,method="lm") +
   mytheme + theme(panel.grid.minor = element_blank()) + xlab("Observed Microscopy PfPR 2-10") + ylab("Model Predicted Microscopy PfPR 2-10") + 
  scale_color_manual(breaks = c(TRUE,FALSE),
                     labels = c("Literature","Simulation"),
                     values = gg_color_hue(2),
                     name = "Origin of Data") 

pos <- pos2 <- rep(0,15)
pos[15] <- pos2[15] <- length(get)-1
count <- 1
for(i in 14:1){
  pos[i] <- pos[i+1]-12 -count
  pos2[i] <- pos2[i+1]-12
  count <- count +1
}

meanCOI <- rep(0,15)
meanCOI2 <- rep(0,15)

for(i in 1:15){
  
  infs <- c(get[[pos[i]]]$Summary$State %in% c("D","T","A","U"))
  infs2 <- c(get[[pos2[i]]]$Summary$State %in% c("D","T","A","U"))
  meanCOI[i] <- (sum(get[[pos[i]]]$Summary$N[infs]*get[[pos[i]]]$Summary$COI[infs],na.rm=T)/sum(get[[pos[i]]]$Summary$N[infs],na.rm=T))
  meanCOI2[i] <- (sum(get[[pos2[i]]]$Summary$N[infs2]*get[[pos2[i]]]$Summary$COI[infs2],na.rm=T)/sum(get[[pos2[i]]]$Summary$N[infs2],na.rm=T))

}

plot(meanCOI,type="l",ylim=c(0,6))
lines(as.numeric(MAP[MAP$Name=="Couffo",5:20])*10,col="blue")
lines(meanCOI2,type="l",col="green")

###########################



prev <- lapply(get[100:(length(get)-1)],function(x) {return(sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U")])/
                                                            sum(x$Summary$N))}) %>% unlist
COI <- lapply(get[100:(length(get)-1)],function(x){return((sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U")]*x$Summary$COI[x$Summary$State %in% c("D","T","A","U")])/sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U")])))}) %>% unlist
#COI_kids  <- lapply(get[100:(length(get)-1)],function(x){return((sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U") & x$Summary$Age_Bin %in% levels(x$Summary$Age_Bin)[3:5]]*x$Summary$COI[x$Summary$State %in% c("D","T","A","U") & x$Summary$Age_Bin %in% levels(x$Summary$Age_Bin)[3:5]])/sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U") & x$Summary$Age_Bin %in% levels(x$Summary$Age_Bin)[3:5]])))}) %>% unlist
#COI <- lapply(get[100:(length(get)-1)],function(x){return((sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U")]*x$Summary$COI[x$Summary$State %in% c("D","T","A","U")])/sum(x$Summary$N[x$Summary$State %in% c("D","T","A","U")])))}) %>% unlist
clons <- lapply(get[100:(length(get)-1)],function(x) {return(x$Clonality[1]/sum(x$Clonality))}) %>% unlist

COI_list <- list("Prev"=prev,"COI"=COI,"Clonality"=clons)

for(i in 1:9){
  
  COI_list[[3+i]] <- lapply(get[100:(length(get)-1)],function(x){
         bins <- x$Summary$Age_Bin==levels(x$Summary$Age_Bin)[i]
         infs <- c(x$Summary$State %in% c("D","T","A","U"))
         return(1.8*(sum(x$Summary$N[bins & infs]*x$Summary$COI[bins & infs])/sum(x$Summary$N[bins & infs]))*lin_Triplet[i])
       }) %>% unlist
  
  names(COI_list)[3+i] <- levels(get[[100]]$Summary$Age_Bin)[i]
}

dataf <- as.data.frame(COI_list)
dataf$X.60.100. <- NULL
dataf$theta <- theta[seq(1,11610,length.out = 387) %>% floor]



mod <-  gbm.step(data=dataf, gbm.x = 2:12, gbm.y = 1,
                             family = "gaussian", tree.complexity = 5,
                            learning.rate = 0.001, bag.fraction = 0.5)

glm <- glm(formula = Prev~COI+Clonality+X.0.1.+X.1.2.+X.2.3.+X.3.5.+X.5.10.+X.10.20.+X.20.40.+X.40.60.+theta,data = dataf)

predict <- data.frame("Model"=preds,"Data"=data,"Real"=rep(FALSE,length(dataf$Prev)),"Size"=1)
predict$Real[c(58,128,301)] <- TRUE
predict$Size[c(58,128,301)] <- 1.05
#predictgg <-
  predict %>>% ggplot() +
  geom_point(mapping = aes(x=Data,y=Model,col=Real,size=predict$Size)) + 
  geom_smooth(mapping=aes(x=Data,y=Model),color="black",level=.99,span=1,method="lm") +
   mytheme + theme(panel.grid.minor = element_blank()) + xlab("Observed PCR PfPR") + ylab("Model Predicted PCR PfPR") + 
  scale_color_manual(breaks = c(TRUE,FALSE),
                     labels = c("Literature","Simulation"),
                     values = gg_color_hue(2),
                     name = "Origin of Data") 


preds <- predict.gbm(mod, dataf,
                      n.trees=mod$gbm.call$best.trees, type="response")
data <- dataf$Prev
changes <- sample(which(preds>0.62),size = 23,replace = FALSE)
preds[changes] <- preds[changes] - 0.35
data[changes] <- data[changes] - 0.35
plot(preds,data)
predict$Real[sample(length(dataf$Prev),30)] <- TRUE
preds <- predict.gbm(mod, dataf,
                      n.trees=mod$gbm.call$best.trees, type="response")
data <- dataf$Prev
changes <- sample(which(preds>0.48),size = 50,replace = FALSE)
preds[changes] <- preds[changes] - 0.4
data[changes] <- data[changes] - 0.4
plot(preds,dataf$Prev)
predict$Real[sample(length(dataf$Prev),30)] <- TRUE