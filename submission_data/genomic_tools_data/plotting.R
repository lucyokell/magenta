# Plotting functions

library(tidyverse)
library(scales)

# plotting theme
text_size=18
text_size2=14
mytheme <- ggplot2::theme_classic() + ggplot2::theme_light() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = text_size, family = "Arial"),
                 axis.title.y = ggplot2::element_text(margin=ggplot2::margin(c(0,10)),size = text_size, family = "Arial"),
                 axis.title.x = ggplot2::element_text(margin=ggplot2::margin(c(5)),size = text_size, family = "Arial"),
                 legend.text = ggplot2::element_text(size = text_size2, family = "Arial"),
                 legend.title = ggplot2::element_text(size = text_size, family = "Arial"),
                 plot.title = ggplot2::element_text(size = text_size, family = "Arial",hjust = 0.5),
                 axis.line = ggplot2::element_line(),
                 strip.background = ggplot2::element_rect(fill="white",color="grey"),
                 strip.text.x = ggplot2::element_text(colour = 'black', size = text_size, family = "Arial"),
                 strip.text.y = ggplot2::element_text(colour = 'black', size = text_size, family = "Arial")
  )

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

exp_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = gaussian(link="log")),...)
}

# Fig 1a -----------------------------------------------------------------------

fig1a <- readRDS("fig1a.rds")

gg1a <- ggplot(fig1a,aes(x=Age,y=COI)) + geom_point(color="#002366", alpha=0.4) + 
  geom_smooth(se=FALSE,span=0.9,method="loess",color="#e62400") +
  scale_x_continuous(breaks=c(0,20,40,60,80)) +
  facet_grid(source~Admin) + mytheme

svg("fig1a.svg",width = 12,height = 5)
gg1a
dev.off()

fig1a_table <- fig1a %>% group_by(Age, Admin, source) %>% 
  summarise(COI = round(mean(COI),2)) %>% 
  pivot_wider(names_from = Age, values_from = COI)
write.table(fig1a_table,row.names = F, file = "supp_table1.txt")

# Fig 1b -----------------------------------------------------------------------

fig1b <- readRDS("fig1b.rds")

gg1b <- ggplot(fig1b,aes(x=s,y=kl,ymin=kl_low,ymax=kl_high,col=SubPatents)) + 
  geom_pointrange() + geom_smooth(se=FALSE) + 
  xlab(expression(paste(zeta," : % sporozoites per infection surviving to produce gametocytes"))) +
  ylab("Total Kullback-Leibler Divergence")  +
  scale_color_manual(name="Parasitaemia \nInclusion",values = c("#e62400","#002366")) + mytheme +
  theme(legend.position = c(0.175,0.875))

svg("fig1b.svg",width = 9,height = 9)
gg1b
dev.off()

# Fig 1c -----------------------------------------------------------------------

fig1c <- readRDS("fig1c.rds")
mois <- read.csv("msp2_coi.csv",stringsAsFactors = FALSE)

mois$prev <- as.numeric(mois$Prevalence.PCR)
mois$moi <- as.numeric(mois$Mean.MOI)
mois$moi_error <- as.numeric(mois$Error.Range.on.MOI)
gg1c <- ggplot(mois, aes(x=prev,y=moi,ymin=moi-moi_error,ymax=moi+moi_error)) + geom_point() + geom_errorbar()

gg1c <- gg1c  + 
  geom_line(data = fig1c, aes(x=prev,y=moi), inherit.aes = FALSE, color="#e62400",size=1.5) +
  xlim(c(0.01,1)) + scale_x_log10(limits = c(0.01,1)) + mytheme + 
  xlab("PCR Prevalence") + ylab(expression(paste(italic("msp2")," COI")))

svg(filename = "fig1c.svg",width = 7, height = 7, pointsize = 12)
gg1c
dev.off()

gg1 <- cowplot::plot_grid(gg1a,NULL,cowplot::plot_grid(gg1b,NULL,gg1c,ncol = 3,labels = c("b","","c"),rel_widths = c(1,0.05,1)),labels="a",ncol=1,rel_heights = c(0.9,0.05,0.9))
svg(filename = "fig1.svg",width = 14, height = 12, pointsize = 12)
gg1
dev.off()
cowplot::save_plot("fig1.png",gg1,base_width  = 14*1.2, base_height = 12*1.2, pointsize = 12)

# Supplementary Figure 6  --------------------------------------------------------

gg_sub6 <- gg1c  + 
  geom_line(data = fig1c, aes(x=prev,y=moi_low), inherit.aes = FALSE,linetype="dashed", color = "blue",size=1.5) + 
  geom_line(data = fig1c, aes(x=prev,y=moi_high), inherit.aes = FALSE,linetype="dashed", color = "blue",size=1.5) + 
  xlim(c(0.01,1)) + scale_x_log10(limits = c(0.01,1)) + mytheme + xlab("PCR Prevalence") + ylab(expression(paste(italic("msp2")," COI")))

svg(filename = "supp_fig6.svg",width = 7, height = 7, pointsize = 12)
gg_sub6
dev.off()
cowplot::save_plot(filename = "supp_fig6.png",gg_sub6,base_width = 7, base_height = 7, pointsize = 12)


# Fig 2 -----------------------------------------------------------------------


fig2a <- readRDS("fig2a.rds")
fig2b <- readRDS("fig2b.rds")
zhu_data <- readRDS("zhu_data.rds")

gg2a <- ggplot(fig2a,
             aes(x = PCR.PfPR, y = P50.IBD, 
                 ymin =  P50.IBD-se, 
                 ymax =  P50.IBD+se,
                 group=type,
                 color=type)) +
  geom_smooth(se=FALSE,span=2,method = lm, formula = y ~ splines::bs(x, 5)) + 
  mytheme + xlab("PCR Prevalence") + 
  ylab("Proportion of within-host highly related parasites") +  
  geom_point() + 
  geom_errorbar() + 
  xlim(c(0,1)) +
  scale_color_brewer(palette = "Paired", name = "Age Group") +  
  theme(legend.position = c(.85, .875), legend.text = element_text(size=12), 
        legend.title = element_text(size=14))

gg2b <- ggplot(fig2b,
              aes(x = mixed,y = P50.IBD, 
                  xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax)) + 
  geom_errorbar(width = 0) + 
  geom_errorbarh(width = 0) + 
  exp_smooth(span = 1,se = FALSE,color = "black") +
  geom_point(aes(shape=Continent), size = 3) +
  geom_errorbar(data  =  zhu_data, mapping = aes(x = mixed_mean,y = ibd_mean,ymin = ymin,ymax = ymax),width = 0) + 
  geom_errorbarh(data  =  zhu_data, mapping = aes(x = mixed_mean,y = ibd_mean,xmin = xmin,xmax = xmax),width = 0) +
  geom_point(data  =  zhu_data, mapping = aes(x = mixed_mean,y = ibd_mean,color = Continent,size = n))  +
  ggrepel::geom_text_repel(data = zhu_data, aes(x = mixed_mean,y = ibd_mean,label = region),segment.size  =  0,box.padding = 0.4) + 
  mytheme + 
  xlab("Fraction Mixed Infections (COI > 1)") + 
  ylab("Mean IBD in Mixed Infections") +  
  scale_color_manual(values = c("orangered1", "royalblue", "chartreuse4"), name = "Origin") + 
  scale_shape_manual(values = c(1,10,15,22), name = "") + 
  scale_size_continuous(breaks = c(100,300,500),name = "Sample size") + 
  theme(legend.position = c(.85, .79), legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14))

gg2 <- cowplot::plot_grid(gg2a,NULL,gg2b,ncol=3,rel_widths = c(1,0.05,1), labels=c("a","","b"))


svg(filename = "fig2.svg",width = 14, height = 7, pointsize = 12)
gg2
dev.off()
cowplot::save_plot(filename = "fig2.png",gg2,base_width = 14, base_height = 7, pointsize = 12)


# Fig 3  -----------------------------------------------------------------------

fig3_data <- readRDS("fig3.rds")
df1 <- fig3_data$a
df2 <- fig3_data$b

fig3 <- list()

for(i in 1:2){
  
  span1 = 0.125
  span2 = 0.05
  
  if(i == 1) {
    melted <- df1$melted
    deeped <- df1$deeped
  }
  
  if(i == 2) {
    melted <- df2$melted
    deeped <- df2$deeped
  }
  
    melted <- melted[melted$sampled == "All",]
    deeped <- deeped[deeped$sampled == "All",]
    
  xl <- c(-0.5,10)
  xi <- 0
  
  if(i == 1) {
    
    y1 <- 0.55
    
    gg3a <- ggplot(melted[melted$metric!="mean_prev" & melted$metric=="% Polygenomic", ],
                  aes(x=time,y=value,group=rep,color=intervention)) + 
      geom_vline(xintercept = 0, linetype="dashed") + 
      geom_line(alpha = 0.2) + 
      xlim(xl) + ylab("") + 
      scale_x_continuous(breaks = c(0,5,10))
    
    gg3a <- gg3a + geom_smooth(data = deeped[deeped$metric!="mean_prev" & deeped$metric == "% Polygenomic" ,], 
                             aes(x=time,y=value,color=intervention),se=FALSE,span=span1) +
      mytheme + theme(legend.box.background = element_blank()) + 
      scale_color_discrete(name="Intervention") +
      theme(legend.position = "none", axis.title.x = element_blank(),  
            strip.background = element_rect(fill = "white", color="white"), panel.spacing = unit(.75, "lines"), 
            strip.placement = "outside") 
    
    gg3b <- ggplot(melted[melted$metric!="mean_prev" & melted$metric=="COI", ],
                   aes(x=time,y=value,group=rep,color=intervention)) + 
      geom_vline(xintercept = 0, linetype="dashed") + 
      geom_line(alpha = 0.2) + 
      xlim(xl) + ylab("") + 
      scale_x_continuous(breaks = c(0,5,10))
    
    gg3b <- gg3b + geom_smooth(data = deeped[deeped$metric!="mean_prev" & deeped$metric == "COI" ,], 
                               aes(x=time,y=value,color=intervention),se=FALSE,span=span1) +
      mytheme + theme(legend.box.background = element_blank()) + 
      scale_color_discrete(name="Intervention") +
      theme(legend.position = "none", axis.title.x = element_blank(),  
            strip.background = element_rect(fill = "white", color="white"), panel.spacing = unit(.75, "lines"), 
            strip.placement = "outside") 
    
    gg3c <- ggplot(melted[melted$metric!="mean_prev" & melted$metric=="% Unique", ],
                   aes(x=time,y=value,group=rep,color=intervention)) + 
      geom_vline(xintercept = 0, linetype="dashed") + 
      geom_line(alpha = 0.2) + 
      xlim(xl) + ylab("") + 
      scale_x_continuous(breaks = c(0,5,10))
    
    gg3c <- gg3c + geom_smooth(data = deeped[deeped$metric!="mean_prev" & deeped$metric == "% Unique" ,], 
                               aes(x=time,y=value,color=intervention),se=FALSE,span=span1) +
      mytheme + theme(legend.box.background = element_blank()) + 
      scale_color_discrete(name="Intervention") +
      theme(legend.position = "none", axis.title.x = element_blank(),  
            strip.background = element_rect(fill = "white", color="white"), panel.spacing = unit(.75, "lines"), 
            strip.placement = "outside") 
    
    reverselog_trans <- function(base = exp(1)) {
      trans <- function(x) -log(x, base)
      inv <- function(x) base^(-x)
      trans_new(paste0("reverselog-", format(base)), trans, inv, 
                log_breaks(base = base), 
                domain = c(1e-100, Inf))
    }
    
    gg3d <- ggplot(melted[melted$metric=="COU"  & melted$sampled=="All",],
                  aes(x=time,y=-log(value),group=rep,color=intervention)) + 
      geom_vline(xintercept = xi, linetype="dashed") + 
      geom_line(alpha = 0.2) + 
      facet_grid(metric~., scales="free_y", switch = "y") +
      xlim(xl) + 
      scale_y_continuous(trans=reverselog_trans(10), 
                         breaks=c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1), 
                         labels=c(expression(`1-10`^{-6}),expression(`1-10`^{-5}),expression(`1-10`^{-4}),
                                  expression(`1-10`^{-3}),expression(`1-10`^{-2}),expression(`1-10`^{-1})))
    
    gg3d <- gg3d + geom_smooth(data = deeped[deeped$metric=="COU" & deeped$sampled=="All",], 
                             aes(x=time,y=-log(value),color=intervention),se=FALSE,span=0.1) + 
      mytheme + theme(legend.box.background = element_blank()) + 
      scale_color_discrete(name="Intervention") + xlab("") +
      theme(legend.position = "none", axis.title = element_blank(),strip.text = element_blank(),panel.background = element_blank(),
            axis.text.y = element_text(size=14),
            strip.background = element_rect(fill = "white", color="white"), panel.spacing = unit(.75, "lines"), plot.background = element_blank(),
            strip.placement = "outside") +scale_x_continuous(breaks = c(0,5,10))
    
    gg3_row <- cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_text("% Polygenomic", x = 0.65, y = y1,angle=90,size = text_size),
                                  gg3a,
                                  cowplot::ggdraw() + cowplot::draw_text("COI", x = 0.65, y = y1,angle=90,size = text_size),
                                  gg3b,
                                  cowplot::ggdraw() + cowplot::draw_text("% Unique", x = 0.65, y = y1,angle=90,size = text_size),
                                  gg3c,
                                  cowplot::ggdraw() + cowplot::draw_text("COU", x = 0.65, y = y1,angle=90,size = text_size),
                                  gg3d,ncol=4,rel_widths = c(0.1,1,0.1,1))

    
  } else {
    
    gg1 <- ggplot(melted[melted$metric=="iIBD" & melted$sampled=="All",],
                  aes(x=time,y=value,group=rep,color=intervention)) + 
      geom_vline(xintercept = xi, linetype="dashed") + 
      geom_line(se=FALSE,alpha = 0.2) + 
      xlim(xl) + scale_y_continuous(limits=c(0.0,1),breaks=c(0,0.25,0.5,0.75,1),labels=c("0.000","0.250","0.500","0.750","1.000")) + 
      facet_wrap(metric~., scales="free_y", switch = "y",ncol=2)  + 
      scale_x_continuous(breaks = c(0,5,10))
    
    gg1a <- gg1 + geom_smooth(data = deeped[deeped$metric!="mean_prev" & deeped$metric=="iIBD" & deeped$sampled=="All",], 
                              aes(x=time,y=value,color=intervention),se=FALSE,span=0.7) + 
      mytheme + theme(legend.box.background = element_blank()) + 
      scale_color_discrete(name="Intervention") +
      theme(legend.position = "none", axis.title = element_blank(), 
            strip.background = element_rect(fill = "white", color="white"), panel.spacing = unit(.75, "lines"), 
            strip.placement = "outside") 
    
    gg2 <- ggplot(melted[melted$metric=="pIBD"  & melted$sampled=="All" ,],
                  aes(x=time,y=value,group=rep,color=intervention)) + 
      geom_vline(xintercept = xi, linetype="dashed") + 
      geom_line(alpha = 0.2) + 
      facet_wrap(metric~., scales="free_y", switch = "y",ncol=2) +
      xlim(xl) + 
      scale_y_log10(limits=c(0.0001,1),breaks=c(1e-4,1e-3,1e-2,1e-1,1e0), 
                    labels=c(expression(`10`^{-4}),expression(`10`^{-3}),expression(`10`^{-2}),
                             expression(`10`^{-1}),expression(`10`^{0})))
    
    gg2 <- gg2 + geom_smooth(data = deeped[deeped$metric!="mean_prev" & deeped$metric=="pIBD"  & deeped$sampled=="All",], 
                             aes(x=time,y=value,color=intervention),se=FALSE,span=0.7) + 
      mytheme + theme(legend.box.background = element_blank()) + 
      scale_color_discrete(name="Intervention") + 
      theme(legend.position = "none",axis.title = element_blank(),
            strip.background = element_rect(fill = "white", color="white"), panel.spacing = unit(.75, "lines"), 
            strip.placement = "outside") + scale_x_continuous(breaks = c(0,5,10))
    
    gg3_row <- cowplot::plot_grid(gg2,gg1a,ncol = 2)
    
  }
  
  gg2 <- ggplot(melted[melted$metric=="mean_prev",],aes(x=time,y=value,group=rep,color=intervention)) + 
    geom_vline(xintercept = xi, linetype="dashed") + 
    geom_line(alpha = 0.2) +  xlim(xl)
  
  gg2 <- gg2 + geom_smooth(data = deeped[deeped$metric=="mean_prev",], 
                           aes(x=time,y=value,color=intervention),se=FALSE,span=span2) + 
    mytheme + theme(legend.box.background = element_blank()) + 
    scale_color_discrete(name="Intervention") + xlab("Time (years)") + 
    ylab("PCR Prevalence") + 
    scale_x_continuous(breaks = c(0,5,10))
  
  top_row <- cowplot::plot_grid(NULL,gg2,NULL,ncol = 3, rel_widths = c(0.75,2,0.25))
  if(i==1){
    fig3[[i]] <- cowplot::plot_grid(top_row,NULL,gg3_row,ncol = 1,rel_heights = c(1,0.05,2.5))
  } else {
    fig3[[i]] <- cowplot::plot_grid(top_row,NULL,gg3_row,ncol = 1,rel_heights = c(1,0.05,1.25))
  }
  
}
lb <- cowplot::ggdraw(grid::linesGrob(x=c(0.1,0.9),y = 0.5,gp = grid::gpar(col="#404040")))
gg3 <- cowplot::plot_grid(fig3[[1]],lb,fig3[[2]],nrow = 3,rel_heights = c(3.5/2.25,0.1,1.05),labels = c("a","","b"))
svg(filename = "fig3.svg",width = 12,height = 15)
gg3
dev.off()
cowplot::save_plot(filename = "fig3.png",gg3,base_width = 12,base_height = 15)

# Supplementary Figure 1 -----------------------------------------------------------------------

fig3_data <- readRDS("fig3.rds")
df1 <- fig3_data$a
df2 <- fig3_data$b

fig3_supp <- list()

for(i in 1:2){
  
  if(i == 1) {
      melted <- df1$melted
      deeped <- df1$deeped
  }
  
  if(i == 2) {
      melted <- df2$melted
      deeped <- df2$deeped
  }
  
  xl <- c(-0.5,10)
  xi <- 0
  
  if(i == 1) {

    gg1 <- ggplot(melted[melted$metric!="mean_prev" & melted$metric!="COU",],
                  aes(x=time,y=value,group=rep,color=intervention)) + 
      geom_vline(xintercept = 0, linetype="dashed") + 
      geom_line(alpha = 0.2) + 
      xlim(xl) + scale_x_continuous(breaks = c(0,5,10)) + 
      facet_grid(metric ~ sampled,scales="free_y",switch="y")
    
    gg1 <- gg1 + geom_smooth(data = deeped[deeped$metric!="mean_prev" & deeped$metric!="COU" ,], 
                             aes(x=time,y=value,color=intervention), se=FALSE, span=0.05) +
      mytheme + theme(legend.box.background = element_blank()) + 
      scale_color_discrete(name="Intervention") +
      theme(legend.position = "none", axis.title = element_blank(),  axis.text.x = element_blank(),
            axis.line.x.bottom = element_blank(),axis.ticks.x = element_blank(),
            strip.background = element_rect(fill = "white", color="white"), panel.spacing = unit(.75, "lines"), 
            strip.placement = "outside") 
    
    reverselog_trans <- function(base = exp(1)) {
      trans <- function(x) -log(x, base)
      inv <- function(x) base^(-x)
      trans_new(paste0("reverselog-", format(base)), trans, inv, 
                log_breaks(base = base), 
                domain = c(1e-100, Inf))
    }
    
    gg2 <- ggplot(melted[melted$metric=="COU" ,],
                  aes(x=time,y=-log(value),group=rep,color=intervention)) + 
      geom_vline(xintercept = xi, linetype="dashed") + 
      geom_line(alpha = 0.2) + 
      facet_grid(metric~sampled, scales="free_y", switch = "y") +
      xlim(xl) + 
      scale_y_continuous(trans=reverselog_trans(10), 
                         breaks=c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1), 
                         labels=c(expression(`1-10`^{-6}),expression(`1-10`^{-5}),expression(`1-10`^{-4}),
                                  expression(`1-10`^{-3}),expression(`1-10`^{-2}),expression(`1-10`^{-1})))
    
    gg2 <- gg2 + geom_smooth(data = deeped[deeped$metric=="COU",], 
                             aes(x=time,y=-log(value),color=intervention),se=FALSE,span=0.05) + 
      scale_color_discrete(name="Intervention") +
      scale_x_continuous(breaks = c(0,5,10)) +
      ylab("COU") +
      mytheme +
      theme(legend.box.background = element_blank(),
            legend.position = "none", axis.title.x = element_blank(),strip.text = element_blank(),
            strip.background = element_rect(fill = "white", color="white"), panel.spacing = unit(.75, "lines"), 
            strip.placement = "outside") 
    
    gg <- cowplot::plot_grid(gg1, gg2, nrow = 2, rel_heights = c(3,1.1))
    
  } else {
    
    gg1 <- ggplot(melted[melted$metric=="iIBD",],
                  aes(x=time,y=value,group=rep,color=intervention)) + 
      geom_vline(xintercept = xi, linetype="dashed") + 
      geom_line(alpha = 0.2) + 
      xlim(xl) + scale_y_continuous(limits=c(0.0,1),breaks=c(0,0.25,0.5,0.75,1),labels=c("0.000","0.250","0.500","0.750","1.000")) + 
      facet_grid(metric~sampled, scales="free_y", switch = "y")  + scale_x_continuous(breaks = c(0,5,10))
    
    gg1a <- gg1 + geom_smooth(data = deeped[deeped$metric!="mean_prev" & deeped$metric=="iIBD",], 
                              aes(x=time,y=value,color=intervention), se=FALSE, span=0.7) + 
      scale_color_discrete(name="Intervention") +
      mytheme +
      theme(legend.box.background = element_blank(),
            legend.position = "none", axis.title = element_blank(),  axis.text.x = element_blank(),
            strip.background = element_rect(fill = "white", color="white"), panel.spacing = unit(.75, "lines"), 
            strip.placement = "outside") 
    
    gg2 <- ggplot(melted[melted$metric=="pIBD" ,],
                  aes(x=time,y=value,group=rep,color=intervention)) + 
      geom_vline(xintercept = xi, linetype="dashed") + 
      geom_line(alpha = 0.2) + 
      facet_grid(metric~sampled, scales="free_y", switch = "y") +
      xlim(xl) + 
      scale_y_log10(limits=c(0.0001,1),breaks=c(1e-4,1e-3,1e-2,1e-1,1e0), 
                    labels=c(expression(`10`^{-4}),expression(`10`^{-3}),expression(`10`^{-2}),
                             expression(`10`^{-1}),expression(`10`^{0})))
    
    gg2 <- gg2 + geom_smooth(data = deeped[deeped$metric!="mean_prev" & deeped$metric=="pIBD" ,], 
                             aes(x=time,y=value,color=intervention),se=FALSE,span=0.7) + 
      mytheme + theme(legend.box.background = element_blank()) + 
      scale_color_discrete(name="Intervention") +
      scale_x_continuous(breaks = c(0,5,10)) +
      theme(legend.position = "none", axis.title.x = element_blank(),strip.text = element_blank(),
            strip.background = element_rect(fill = "white", color="white"), panel.spacing = unit(.75, "lines"), 
            strip.placement = "outside") + ylab("pIBD")
    
    gg <- cowplot::plot_grid(gg1a,gg2,nrow = 2, rel_heights = c(1,1))
  }
  
  gg2 <- ggplot(melted[melted$metric=="mean_prev",],
                aes(x=time,y=value,group=rep,color=intervention)) + 
    geom_vline(xintercept = xi, linetype="dashed") + 
    geom_line(alpha = 0.2) +  xlim(xl)
  
  gg2 <- gg2 + geom_smooth(data = deeped[deeped$metric=="mean_prev",], 
                           aes(x=time,y=value,color=intervention),se=FALSE,span=0.05) + 
    mytheme + 
    theme(legend.box.background = element_blank()) + 
    scale_color_discrete(name="Intervention") +
    xlab("Time (years)") + 
    ylab("PCR Prevalence")
  
  top_row <- cowplot::plot_grid(NULL,gg2+scale_x_continuous(breaks = c(0,5,10)),NULL,ncol = 3, rel_widths = c(0.75,2,0.25))
  if(i==1){
    fig3_supp[[i]] <- cowplot::plot_grid(top_row,gg,ncol = 1,rel_heights = c(1.2,4))
  } else {
    fig3_supp[[i]] <- cowplot::plot_grid(top_row,gg,ncol = 1,rel_heights = c(1,2))
  }
  
}

gg_sub1 <- cowplot::plot_grid(fig3_supp[[1]],NULL,fig3_supp[[2]],nrow = 3,rel_heights = c(5/3,.05,3.6/3),labels = c("a","","b"))
svg(filename = "supp_fig1.svg",width = 12,height = 18, pointsize = 12)
gg_sub1
dev.off()
cowplot::save_plot(filename = "supp_fig1.png",gg_sub1,base_width = 13,base_height = 20)

# Correlation Table -----------------------------------------------------------------------

df1$melted$prev <- df1$melted$value[df1$melted$variable=="mean_prev"]
df2$melted$prev <- df2$melted$value[df2$melted$variable=="mean_prev"]
melted <- rbind(df1$melted,df2$melted)
melted$sampled <- factor(melted$sampled,levels(melted$sampled))
melted$metric <- factor(melted$metric, levels = c("% Polygenomic","COI","% Unique","COU","mean_prev","iIBD","pIBD"))

text <- group_by(melted[melted$time>-1 & melted$time <9 & melted$prev>0 & melted$value>0,],metric,sampled) %>% 
  summarise(cor = round(cor(prev,value, method="kendall"),2))
write.table(spread(text[!text$metric %in% c("mean_prev",NA),],key=metric,value = cor),
            row.names = F, file = "table1.txt")

# Fig 4 -----------------------------------------------------------------------
fig4_data <- readRDS("fig4.rds")
fig4_dfa <- fig4_data$a
fig4_dfb <- fig4_data$b
fig4_dfc <- fig4_data$c


g1 <- ggplot(fig4_dfa[fig4_dfa$unphased==FALSE,],aes(x=ss,y=power,color=Year_PCR,shape = Year_PCR)) + geom_hline(yintercept=0.8,linetype="dashed") + geom_smooth(se=FALSE,span=0.5) + 
  xlab("Sample Size") + ylab("Power") + geom_point()  + scale_color_grey() + mytheme + facet_wrap(~metric,ncol = 4) + 
  ggtitle("Starting PCR Prevalence = 45%\n") + theme(legend.title = element_blank(), panel.spacing.x=unit(1, "lines"))

g2 <- ggplot(fig4_dfb[fig4_dfb$unphased==FALSE,],aes(x=ss,y=power,color=Year_PCR,shape = Year_PCR)) + geom_hline(yintercept=0.8,linetype="dashed") + geom_smooth(se=FALSE,span=0.5) + 
  xlab("Sample Size") + ylab("Power") + geom_point()  + scale_color_grey() + mytheme + facet_wrap(~metric,ncol = 4) + 
  ggtitle("Starting PCR Prevalence = 22.5%\n") + theme(legend.title = element_blank(), panel.spacing.x=unit(1, "lines"))

g3 <- ggplot(fig4_dfc[fig4_dfc$unphased==FALSE,],aes(x=ss,y=power,color=Year_PCR,shape = Year_PCR)) + geom_hline(yintercept=0.8,linetype="dashed") + geom_smooth(se=FALSE,span=0.5) + 
  xlab("Sample Size") + ylab("Power") + geom_point()  + scale_color_grey() + mytheme + facet_wrap(~metric,ncol = 4) + 
  ggtitle("Starting PCR Prevalence = 22.5%\n") + theme(legend.title = element_blank())

lb <- cowplot::ggdraw(grid::linesGrob(x=c(0.1,0.8),y = 0.5,gp = grid::gpar(col="#404040")))
gg4 <- cowplot::plot_grid(g1,lb,g2,lb,g3,nrow=5,rel_heights = c(1,.1,1,.1,1),labels = c("a","","b","","c"))
svg(filename = "fig4.svg", height = 12*1.2,width = 10*1.2)
gg4
dev.off()
cowplot::save_plot(filename = "fig4.png",gg4, base_height = 12*1.2,base_width = 10*1.2)


g1 <- ggplot(fig4_dfa[fig4_dfa$unphased==TRUE,],aes(x=ss,y=power,color=Year_PCR,shape = Year_PCR)) + geom_hline(yintercept=0.8,linetype="dashed") + geom_smooth(se=FALSE,span=0.5) + 
  xlab("Sample Size") + ylab("Power") + geom_point()  + scale_color_grey() + mytheme + facet_wrap(~metric,ncol = 4) + 
  ggtitle("Starting PCR Prevalence = 45%\n") + theme(legend.title = element_blank(), panel.spacing.x=unit(1, "lines"))

g2 <- ggplot(fig4_dfb[fig4_dfb$unphased==TRUE,],aes(x=ss,y=power,color=Year_PCR,shape = Year_PCR)) + geom_hline(yintercept=0.8,linetype="dashed") + geom_smooth(se=FALSE,span=0.5) + 
  xlab("Sample Size") + ylab("Power") + geom_point()  + scale_color_grey() + mytheme + facet_wrap(~metric,ncol = 4) + 
  ggtitle("Starting PCR Prevalence = 22.5%\n") + theme(legend.title = element_blank(), panel.spacing.x=unit(1, "lines"))

g3 <- ggplot(fig4_dfc[fig4_dfc$unphased==TRUE,],aes(x=ss,y=power,color=Year_PCR,shape = Year_PCR)) + geom_hline(yintercept=0.8,linetype="dashed") + geom_smooth(se=FALSE,span=0.5) + 
  xlab("Sample Size") + ylab("Power") + geom_point()  + scale_color_grey() + mytheme + facet_wrap(~metric,ncol = 4) + 
  ggtitle("Starting PCR Prevalence = 22.5%\n") + theme(legend.title = element_blank())

lb <- cowplot::ggdraw(grid::linesGrob(x=c(0.1,0.9),y = 0.5,gp = grid::gpar(col="grey")))
gg_sub2 <- cowplot::plot_grid(g1,lb,g2,lb,g3,nrow=5,rel_heights = c(1,.1,1,.1,1),labels = c("a","","b","","c"))
svg(filename = "supp_fig2.svg", height = 12*1.2,width = 10*1.2)
gg_sub2
dev.off()
cowplot::save_plot(filename = "supp_fig2.png",gg_sub2, base_height = 12*1.2,base_width = 10*1.2)

# Fig 5 -----------------------------------------------------------------------

fig5a_data <- readRDS("fig5a.rds")
fig5b_data <- readRDS("fig5b.rds")

gg5a <- ggplot(fig5a_data,aes(x=Observed.Microscopy2.10.PfPR, y = Predicted.Microscopy2.10.PfPR ,color = Origin,size=Origin)) + 
  geom_abline(slope = 1,intercept = 0) + 
  geom_point() + 
  geom_text(aes(label=Location),nudge_x = c(0.055,-0.035,0.03,-0.03,0.05), nudge_y=c(-0.002,0.01,-0.02,0.03,0.006)) + 
  mytheme  + scale_size_manual(name = "Origin of Data:", values = c(1.5,4)) + 
  theme(legend.position = "top",legend.text = element_text(size=14), legend.title = element_text(size=16)) +
  scale_color_manual(name = "Origin of Data:",values = c("#e62400","#002366")) + 
  labs(x=expression("Observed Microscopy Prevalence "["2-10"]),
       y=expression("Predicted Microscopy Prevalence  "["2-10"]))

gg5b <- ggplot(fig5b_data, 
                aes(y=value, x=variable, fill=meta_data )) + 
  geom_bar(stat="identity",position = position_dodge(width=0.8), color="Black", width=0.75) + 
  theme_bw() + scale_fill_brewer(palette = "Paired", name = "Meta Data:", labels = c("None", "Age", "Age and Clinical Status")) +
  mytheme + scale_x_discrete(labels = parse(text = levels(fig5b_data$variable))) + 
  ylab("Model Error") + xlab("Model Metric") + 
  theme( legend.text = element_text(size=14), legend.title = element_text(size=16),
        legend.key.width = unit(10,"mm"),legend.position = "top")

gg5 <- cowplot::plot_grid(gg5a, NULL, gg5b, rel_widths = c(1,0.1,1), nrow = 1,labels = c("a","","b"))
svg(filename = "fig5.svg",width = 12*1.2,height = 6*1.2)
gg5
dev.off()
cowplot::save_plot(filename = "fig5.png",gg5,base_width = 12*1.2,base_height = 6*1.2)

# Supp Fig 3 -----------------------------------------------------------------------

supp_fig3_data <- readRDS("supp_fig3.rds")

gg_sub3 <- ggplot(supp_fig3_data,aes(x=reorder(Variable,Variable.Importance),y=Variable.Importance, fill = Variable.Importance)) + 
  geom_bar(stat="identity") + coord_flip() + mytheme + xlab("") + 
  ylab("") +
  scale_fill_gradientn(name = "Variable\nImportance", 
                       colours = scales::seq_gradient_pal("#112842","#54b0f5", "Lab")(seq(0,1,length.out=9)))

svg(filename = "supp_fig3.svg",width = 8,height = 6)
gg_sub3
dev.off()
cowplot::save_plot(filename = "supp_fig3.png",gg_sub3,base_width = 8,base_height = 6)

# Supp Fig 4 -----------------------------------------------------------------------

supp_fig4_data <- readRDS("supp_fig4.rds")

gg_sub4 <- ggplot(supp_fig4_data,aes(x=Sample.Size, y=value,color=variable)) + 
  geom_vline(xintercept = 200, linetype="dashed") + 
  geom_line(lwd=1) + geom_point() + 
  scale_color_manual(labels = parse(text=c("1 - R^2","RMSE","MAE")),
                                                       values = gg_color_hue(3),
                                                       name = "") + 
  mytheme + xlab("Sample Size") + ylab("Model Error") 

svg(filename = "supp_fig4.svg",width = 6,height = 6)
gg_sub4
dev.off()
cowplot::save_plot(filename = "supp_fig4.png",gg_sub4,base_width = 6,base_height = 6)


# Supp Fig 5 -----------------------------------------------------------------------

supp_fig5_data <- readRDS("supp_fig5.rds")

gg_sub5 <- ggplot(supp_fig5,aes(x=age,color=status,y=y,ymin=ymin,ymax=ymax)) + 
  geom_errorbar(position = position_dodge(width=0.7),width=0.4,lwd=0.5) + 
  geom_point(position=position_dodge(width=0.7)) +
  ylim(c(0,10)) + facet_grid(.~admin) + ylab("COI") + xlab("Age Group") + 
  mytheme + 
  scale_y_continuous(breaks=seq(0,10,1)) +
  scale_color_discrete(name="Infection Status",labels=c("Asymptomatic","Clinical"))

svg(filename = "supp_fig5.svg", height = 6,width = 9)
gg_sub5
dev.off()
cowplot::save_plot(filename = "supp_fig5.png",gg_sub5, base_height = 6, base_width = 9)

# Supp Fig 7 -----------------------------------------------------------------------

supp_fig7_data <- readRDS("supp_fig7.rds")
supp_fig7_data$Source <- factor(supp_fig7_data$Source,levels=c("THE REAL McCOIL","MODEL OUTPUT"))

gg_sub7 <- ggplot(supp_fig7_data,aes(x=Age,y=COI)) + geom_point(color="blue", alpha=0.4) + 
  geom_smooth(se=FALSE,span=0.9,method="loess",color="red") +
  scale_x_continuous(breaks=c(0,20,40,60,80)) +
  facet_grid(~Source) + mytheme

svg("supp_fig7.svg",width = 7,height = 4)
gg_sub7
dev.off()
cowplot::save_plot("supp_fig7.png",gg_sub7,base_width = 7,base_height = 4)
