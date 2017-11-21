mytheme <- ggplot2::theme_classic() + ggplot2::theme_light() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 12, family = "Arial"),
                 axis.title.y = ggplot2::element_text(margin=ggplot2::margin(c(0,10)),size = 18, family = "Arial"),
                 axis.title.x = ggplot2::element_text(margin=ggplot2::margin(c(10)),size = 18, family = "Arial"),
                 legend.text = ggplot2::element_text(size = 16, family = "Arial"),
                 legend.title = ggplot2::element_text(size = 18, family = "Arial"),
                 plot.title = ggplot2::element_text(size = 18, family = "Arial",hjust = 0.5),
                 axis.line = ggplot2::element_line(),
                 panel.border = ggplot2::element_rect(),
                 strip.background = ggplot2::element_rect(fill="white"),
                 strip.text.x = ggplot2::element_text(colour = 'black', size = 14, family = "Arial")
  )

##afri
df_kin <- data.frame("Time"=2000:2015,"Prev"=as.numeric(MAP[which(MAP$Name=="Kinshasa"),5:20]))

kin_coi <- all_df %>>% ggplot() + geom_ribbon(mapping = aes(x=Time,ymin=Prev - PrevCI,ymax=Prev+PrevCI),color="blue",fill="blue",alpha=0.2) +
  geom_line(mapping = aes(x=Time,y=Prev),color="blue",linetype="dashed") + 
  geom_line(mapping = aes(x=Time,y=meanCOI*1/7),color="red",linetype="dashed") + 
  geom_ribbon(color="red",fill="red",alpha=0.1,mapping=aes(x=Time,ymin=meanCOI*1/7 - ciCOI*1/7,ymax=meanCOI*1/7+ciCOI*1/7)) +
  scale_x_continuous(name = "Year", breaks = 2000:2015) +
  scale_y_continuous(
    name  = "Micrscopy PfPR 2-10", 
    sec.axis = sec_axis(~ . *  7/1 , name = "Complexity of Infection"), 
    limits = c(0.1, 0.6)) + 
  mytheme + theme(panel.grid.minor = element_blank()) +  geom_line(data = df_kin,mapping=aes(x=Time,y=Prev),linetype="dashed",color="blue",size=1.5) 


kin_clon <- all_df %>>% ggplot() + geom_ribbon(mapping = aes(x=Time,ymin=Prev - PrevCI,ymax=Prev+PrevCI),color="blue",fill="blue",alpha=0.2) +
  geom_line(mapping = aes(x=Time,y=Prev),color="blue",linetype="dashed") + 
  geom_line(mapping = aes(x=Time,y=1-Clonality),color="red",linetype="dashed") + 
  geom_ribbon(color="red",fill="red",alpha=0.1,mapping=aes(x=Time,ymin=1-(Clonality - ClonalityCI),ymax=1-(Clonality + ClonalityCI))) +
  scale_x_continuous(name = "Year", breaks = 2000:2015) +
  scale_y_continuous(
    name  = "Micrscopy PfPR 2-10", 
    sec.axis = sec_axis(~ . *  1/1 , name = "Clonality"), 
    limits = c(0,1)) + 
  mytheme + theme(panel.grid.minor = element_blank()) +  geom_line(data = df_kin,mapping=aes(x=Time,y=Prev),linetype="dashed",color="blue",size=1.5) 

df_thies <- data.frame("Time"=seq(2000.001,2015.001,1.000),"Prev"=as.numeric(MAP[which(MAP$Name=="Thies"),5:20]))
thies_coi <- all_df_thies %>>% ggplot() + geom_ribbon(mapping = aes(x=Time,ymin=Prev - PrevCI,ymax=Prev+PrevCI),color="blue",fill="blue",alpha=0.2) +
  geom_line(mapping = aes(x=Time,y=Prev),color="blue",linetype="dashed") + 
  geom_line(mapping = aes(x=Time,y=meanCOI*1/5),color="red",linetype="dashed") + 
  geom_ribbon(color="red",fill="red",alpha=0.1,mapping=aes(x=Time,ymin=meanCOI*1/5 - ciCOI*1/5,ymax=meanCOI*1/5+ciCOI*1/5)) +
 
  scale_x_continuous(name = "Year", breaks = 2000:2015) +
  scale_y_continuous(
    name  = "Micrscopy PfPR 2-10", 
    sec.axis = sec_axis(~ . *  5/1 , name = "Complexity of Infection"), 
    limits = c(0, 0.4)) + 
  mytheme + theme(panel.grid.minor = element_blank()) +  geom_line(data = df_thies,mapping=aes(x=Time,y=Prev),linetype="dashed",color="blue",size=1.5) 


thies_clon <- all_df_thies %>>% ggplot() + geom_ribbon(mapping = aes(x=Time,ymin=Prev - PrevCI,ymax=Prev+PrevCI),color="blue",fill="blue",alpha=0.2) +
  geom_line(mapping = aes(x=Time,y=Prev),color="blue",linetype="dashed") + 
  geom_line(mapping = aes(x=Time,y=1-Clonality),color="red",linetype="dashed") + 
  geom_ribbon(color="red",fill="red",alpha=0.1,mapping=aes(x=Time,ymin=1-(Clonality - ClonalityCI),ymax=1-(Clonality + ClonalityCI))) +
  scale_x_continuous(name = "Year", breaks = 2000:2015) +
  scale_y_continuous(
    name  = "Micrscopy PfPR 2-10", 
    sec.axis = sec_axis(~ . *  1/1 , name = "Clonality"), 
    limits = c(0,1)) + 
  mytheme + theme(panel.grid.minor = element_blank()) +  geom_line(data = df_thies,mapping=aes(x=Time,y=Prev),linetype="dashed",color="blue",size=1.5) 



pdf(file = paste0(res_dir,"/newfigure1.pdf",sep=""),width=20,height=15,family = "Times")
gridExtra::grid.arrange(kin_coi,kin_clon,thies_coi,thies_clon,ncol=2)
dev.off()
embed_fonts(paste0(res_dir,"/newfigure1.pdf",sep=""), outfile=paste0(res_dir,"/newfigure1_embed.pdf",sep=""))


df_kouffo <- data.frame("Time"=seq(2000.001,2015.001,1.000),"Prev"=as.numeric(MAP[which(MAP$Name=="Couffo"),5:20]))
kouffo_coi <- all_df_kouffo %>>% ggplot() + geom_ribbon(mapping = aes(x=Time,ymin=Prev - PrevCI,ymax=Prev+PrevCI),color="blue",fill="blue",alpha=0.2) +
  geom_line(mapping = aes(x=Time,y=Prev),color="blue",linetype="dashed") + 
  geom_point(mapping = aes(x=Time,y=noframeCOI*1/10),color="red") + 
  #geom_ribbon(color="red",fill="red",alpha=0.1,mapping=aes(x=Time,ymin=noframeCOI*1/10 - ciCOI*1/10,ymax=noframeCOI*1/10+ciCOI*1/10)) +
    
   # geom_smooth(color="red",fill="red",alpha=0.1,mapping=aes(x=Time,noframeCOI*1/10),level=.99,span=1,method="lm") +
  
  scale_x_continuous(name = "Year", breaks = 2000:2015) +
  scale_y_continuous(
    name  = "Micrscopy PfPR 2-10", 
      sec.axis = sec_axis(~ . *  10/1 , name = "Complexity of Infection"), 
    limits = c(0.15, 0.75)) + 
  mytheme + theme(panel.grid.minor = element_blank()) +  geom_line(data = df_kouffo,mapping=aes(x=Time,y=Prev),linetype="dashed",color="blue",size=1.5) 
  

      kouffo_coi_grid <- all_df_kouffo %>>% ggplot() + mytheme + #geom_ribbon(mapping = aes(x=Time,ymin=Prev - PrevCI,ymax=Prev+PrevCI),color="blue",fill="blue",alpha=0.2) + mytheme
      #geom_line(mapping = aes(x=Time,y=Prev),color="blue",linetype="dashed") + 
      #geom_point(mapping = aes(x=Time,y=noframeCOI*1/10,col=norframegroup)) + 
      #geom_ribbon(color="red",fill="red",alpha=0.1,mapping=aes(x=Time,ymin=noframeCOI*1/10 - ciCOI*1/10,ymax=noframeCOI*1/10+ciCOI*1/10)) +
      
      #geom_smooth(alpha=0.1,mapping=aes(x=Time,y=noframeCOI*1/10,col=norframegroup,fill=norframegroup),span=0.75,level=.95) +
      facet_grid(~norframegroup) + 
      scale_x_continuous(name = "Year", breaks = c(2000,2005,2010,2015)) +
      scale_y_continuous(name  = "Micrscopy PfPR 2-10",sec.axis = sec_axis(~ . *  10/1 , name = "Complexity of Infection"), limits = c(0.1, 0.8)) + 
      theme(panel.grid.minor = element_blank(),
                      panel.spacing = unit(2, "lines"),legend.position='top',
            strip.background = element_blank(),
            strip.text.x = element_blank()) + 
      geom_line(data = df_kouffo,aes(x=Time,y=Prev),linetype="dashed",color="black",size=1.5) +
      scale_color_manual(breaks = c("0 - 2","2 - 5","5 - 10","10 - 20","20 - 100"),
                                  labels = c("0 - 2","2 - 5","5 - 10","10 - 20","20 - 100"),
                                  values = gg_color_hue(5),
                                  name = "COI by \nAge (Years)") 
  
  
  

  kouffo_clon <- all_df_kouffo %>>% ggplot() + geom_ribbon(mapping = aes(x=Time,ymin=Prev - PrevCI,ymax=Prev+PrevCI),color="blue",fill="blue",alpha=0.2) +
  geom_line(mapping = aes(x=Time,y=Prev),color="blue",linetype="dashed") + 
  #  geom_line(mapping = aes(x=Time,y=1-Clonality),color="red",linetype="dashed") + 
  # geom_ribbon(color="red",fill="red",alpha=0.1,mapping=aes(x=Time,ymin=1-(Clonality - ClonalityCI),ymax=1-(Clonality + ClonalityCI))) +
  scale_x_continuous(name = "Year", breaks = 2000:2015) +
  scale_y_continuous(
    name  = "Micrscopy PfPR 2-10", 
    #  sec.axis = sec_axis(~ . *  1/1 , name = "Clonality"), 
    limits = c(0,1)) + 
  mytheme + theme(panel.grid.minor = element_blank()) +  geom_line(data = df_kouffo,mapping=aes(x=Time,y=Prev),linetype="dashed",color="blue",size=1.5) 
