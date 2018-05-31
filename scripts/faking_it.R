mytheme <- ggplot2::theme_classic() + ggplot2::theme_light() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 12, family = "Calibri"),
                 axis.title.y = ggplot2::element_text(margin=ggplot2::margin(c(0,10)),size = 14, family = "Calibri"),
                 axis.title.x = ggplot2::element_text(margin=ggplot2::margin(c(5)),size = 14, family = "Calibri"),
                 legend.text = ggplot2::element_text(size = 12, family = "Calibri"),
                 legend.title = ggplot2::element_text(size = 12, family = "Calibri"),
                 plot.title = ggplot2::element_text(size = 12, family = "Calibri",hjust = 0.5),
                 axis.line = ggplot2::element_line(),
                 strip.background = element_rect(fill="white",colour = "grey",size = 0.2),
                 strip.text.x = ggplot2::element_text(colour = 'black', size = 10, family = "Calibri")
  )

## fake results

months <- seq(1/12,15,1/12)
prev <- seq(0.1,0.85,0.01)
COI_rel <- exp(5*prev^2.1)
unq <- exp(0.5*prev^4.1)
unique_rel <- rev(max(unq)/unq)/max(unq)
unique_rel

fifteen <- c(rep(.6,60),c(0.6*exp(c(1:120)/-87))) # c(.6,0.6*exp(c(1:120)/-87)) %>% plot
thirty <- c(rep(.6,60),c(0.6*exp(c(1:120)/-173))) # c(.6,0.6*exp(c(1:120)/-173))
fourtfive <-c(rep(.6,60),c(0.6*exp(c(1:120)/-419))) # c(.6,0.6*exp(c(1:120)/-419))


stoch <- function(x) {
  
  x <- x+runif(180,-0.02,0.02)
  mult <- seq(0,1,length.out = 180)
  out <- x
  for(i in 2:length(x)){
   out[i] <- x[i]+runif(1,min = -0.04*mult[i],max=0.04*mult[i]) 
  }
  
  return(out)
}

fake_prevs_1 <- as.vector(replicate(10,expr = stoch(fifteen),simplify = TRUE))
fake_prevs_2 <- as.vector(replicate(10,expr = stoch(thirty),simplify = TRUE))
fake_prevs_3 <- as.vector(replicate(10,expr = stoch(fourtfive),simplify = TRUE))

all_prev <- rep(c(fake_prevs_1,fake_prevs_2,fake_prevs_3),6) + runif(5400*6,-0.02,0.02)

COI <- sapply(c(fake_prevs_1,fake_prevs_2,fake_prevs_3),function(x) COI_rel[which.min(abs(x-prev))] + runif(1,-0.1,0.1))
COI <- c(COI,COI*0.9,COI*1.25,COI*0.85,COI*1.2,COI*0.9)
COI <- COI+runif(length(COI),-0.1,0.1)
COI[COI<=1] <- 1+runif(sum(COI<=1),0,0.1)
Mono <- (COI/max(COI))
Mono <- Mono^0.25
Unique <- sapply(all_prev,function(x) unique_rel[which.min(abs(x-prev))] + runif(1,-0.1,0))

l <- length(COI)


df0 <- data.frame("Time"=months,"Prev"=c(fake_prevs_1,fake_prevs_2,fake_prevs_3),
                 "Rep"=(sort(rep(1:10,180))),"ScaleUp"=(c(rep("High",1800),rep("Medium",1800),rep("Low",1800))),
                 "COI"=COI)

df1 <- data.frame("Time"=months,"Prev"=c(fake_prevs_1,fake_prevs_2,fake_prevs_3),
                 "Rep"=(sort(rep(1:10,180))),"ScaleUp"=(c(rep("High",1800),rep("Medium",1800),rep("Low",1800))),
                 "COI"=COI,"Sampling"=factor(c(rep("Random",l/6),rep("0-5 Years",l/6),rep("5-15 Years",l/6),rep("15+ Years",l/6),rep("Clinical",l/6),rep("Asymptomatics",l/6)),
                                             levels=c("Random","0-5 Years","5-15 Years","15+ Years","Clinical","Asymptomatics")),
                 "Mono"=Mono,"Unique"=Unique
)

g0 <- ggplot(df0,aes(x = Time,y=Prev,color=ScaleUp)) +
  geom_line(aes(group=interaction(ScaleUp,Rep)),alpha=0.2) +
  geom_smooth(aes(x = Time,y=Prev,color=ScaleUp),span=0.01,se=F) +
  ylab("PCR PfPR") +
  mytheme

g0_row <- gridExtra::grid.arrange(g0,g0,ncol=2)

g1 <- ggplot(df1,aes(x = Time,y=COI,color=ScaleUp)) + 
  geom_line(aes(group=interaction(ScaleUp,Rep)),alpha=0.2) +
  geom_smooth(aes(x = Time,y=COI,color=ScaleUp),span=0.01,se=F) +
  ylab("Mean Population COI") +
  mytheme + facet_grid(.~Sampling)
g1

g2 <- ggplot(df1,aes(x = Time,y=Mono,color=ScaleUp)) + 
  geom_line(aes(group=interaction(ScaleUp,Rep)),alpha=0.2) +
  geom_smooth(aes(x = Time,y=Mono,color=ScaleUp),span=0.01,se=F) +
  ylab("% Population Polygenomic") +
  mytheme + facet_grid(.~Sampling)
g2

g3 <- ggplot(df1,aes(x = Time,y=Unique,color=ScaleUp)) + 
  geom_line(aes(group=interaction(ScaleUp,Rep)),alpha=0.2) +
  geom_smooth(aes(x = Time,y=Unique,color=ScaleUp),span=0.01,se=F) +
  ylab("% Unique Samples") +
  mytheme + facet_grid(.~Sampling)
g3

windows()
gridExtra::grid.arrange(g0,
                        g1+theme(legend.position = "none"),
                        g2+theme(legend.position = "none"),
                        g3+theme(legend.position = "none"),
                        layout_matrix = rbind(c(NA,1,1, NA),
                                              c(2, 2, 2, 2),
                                              c(3, 3, 3, 3),
                                              c(4, 4, 4, 4)))


