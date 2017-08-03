mosquto_species <- read.csv("M:/OJ/MAGENTA/inst/extdata/mosquito_species.txt",header=T,sep="\t")
oocyst_counts <- read.csv("M:/OJ/MAGENTA/inst/extdata/mosquito_oocyst_counts.txt",header=T,sep="\t")
counts <- oocyst_counts[,intersect(which(mosquto_species$Species=="Pf"),which(mosquto_species$Strain=="Field"))]
counts <- counts[!is.na(counts)]
counts <- counts[counts>0]
tabled <- table(counts)
input <- matrix(as.numeric(cbind(names(tabledz),tabledz)),ncol=2)
wh2 <-MASS::fitdistr(counts,"negative binomial")
wh3 <-MASS::fitdistr(counts,"poisson")
wh4 <-MASS::fitdistr(counts,"gamma")
wh <- preseqR::preseqR.ztnb.em(input)
wh3 <- pscl::zeroinfl(formula = counts ~ . | ., data, dist = "negbin")
table(actuar::rztnbinom(1000,wh$size,prob=wh$size/(wh$size+wh$mu)))

plot(actuar::dztnbinom(1:300,wh3$size,prob=wh3$size/(wh3$size+wh3$mu)),type="l",col="red")
lines(actuar::dztnbinom(1:300,wh$size,prob=wh$size/(wh$size+wh$mu)),type="l",col="blue")
lines(actuar::dztnbinom(1:300,wh2$estimate[1],prob=wh2$estimate[1]/(wh2$estimate[1]+wh2$estimate[2])),type="l",col="green")
lines(density(counts,from=1),col="blue")

plot(as.numeric(names(tabled)),tabled/length(counts),col="black",pch=16,cex=0.5)
plot(actuar::dztnbinom(1:300,wh$size,prob=wh$size/(wh$size+wh$mu)),type="l",col="blue",xlab="Oocyst Count",ylab="Density")
lines(dnbinom(0:300,wh2$estimate[1],prob=wh2$estimate[1]/(wh2$estimate[1]+wh2$estimate[2])),type="l",col="green")


plot(table)
lines(actuar::dztnbinom(1:300,wh2$estimate[1],prob=wh2$estimate[1]/(wh2$estimate[1]+wh2$estimate[2])),type="l",col="orange")

plot(sort(counts),col="black")
rands <- actuar::rztnbinom(length(counts),wh2$estimate[1],prob=wh2$estimate[1]/(wh2$estimate[1]+wh2$estimate[2]))
points(sort(rands),col="blue")

w2 <- MASS::fitdistr(counts,"gamma")
lines(dpois(1:300,wh3$estimate),type="l",col="orange")
# spoirozoite distribution
round(1/diff(cumsum(actuar::dztnbinom(x = 1:10,size = 2,prob = 0.85))),digits = 0)




loglik.negbin.zt <- function(x, MOI, freq) {
  size=x[1]
  prob=x[2]
  lik <- dnbinom(MOI,prob=prob,size=size) / (1-dnbinom(0,size=size,prob=prob))
  loglik.i<-log(lik)*freq
  sum(loglik.i)
}

loglik.pois.zt <- function(mu, MOI, freq) {
  lik <- dpois(MOI,mu) / (1-dpois(0,mu))
  loglik.i<-log(lik)*freq
  sum(loglik.i)
}


optimize(f=loglik.pois.zt, MOI=as.numeric(rownames(tabled)), freq=tabled, lower=0.00001, upper=100, maximum=T)

ztnb <- optim(c(0.1,0.5), loglik.negbin.zt, MOI=as.numeric(rownames(tabled)), freq=tabled)
fitdistr(z$MOI, "Negative Binomial")
test<-fitdistr(counts, )
logLik(test)  ## get a slightly different size result, but same mean. But likelihood very similar.

install.packages("C:/Users/Oliver/Desktop/pfmix_1.0.tar.gz", repos = NULL, type="source")
