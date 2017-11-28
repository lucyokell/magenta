
wwhat2 <- what
what <- ghana_fit
COIs <- Convert_Barcode_Vectors(ghana_fit$populations_event_and_strains_List)
COIs$COI[is.na(COIs$COI)] <- 0

inf_states <- c("S","D","A","U","T","P")

Triplet_Model <- c(0.6116584564860428,0.4635139573070608,0.4641050903119869,0.40752052545156003,
  0.33986863711001647,0.3411165845648605,0.23376026272578,0.14311986863711001)
x <- 1:8
lm <- lm(Triplet_Model~x)
predict.lm(lm)
lin_Triplet <- predict.lm(lm)

df <- data.frame("COI"=COIs$COI,"Age"=what$population_List$Ages/365,"Status"=inf_states[what$population_List$Infection_States+1])
df$Age[df$Age==0] <- 1/365
df <- dplyr::mutate(df,Age_Bin = cut(Age,breaks = c(0,1,3,5,10,20,40,60,100)))
df <- dplyr::mutate(df,COI_Detected_Triplet_Model = COI*lin_Triplet[Age_Bin])
COIs <- Convert_Barcode_Vectors(what$populations_event_and_strains_List,sub_patents_included = FALSE)
df <- dplyr::mutate(df,COI_Detected_From_Model = COIs$COI)
df <- dplyr::mutate(df,COI_Detected_From_Model_Plus_Triple = COI_Detected_From_Model*lin_Triplet[Age_Bin])

sdf <- summarySE(df[df$COI>0,], measurevar="COI_Detected_Triplet_Model", groupvars=c("Age_Bin"))

windows()
names(sdf)[3] <- "COI"
sdf$COI <- sdf$COI*1.8
ggplot(totalsdf, aes(x=Age_Bin, y=COI, fill=total)) + 
  geom_bar(stat="identity", width=.75, position = "dodge") +
  geom_errorbar(aes(ymin=COI-ci, ymax=COI+ci),
                width=.5,                    # Width of the error bars
                position=position_dodge(.75)) + 
  scale_fill_discrete(name="COI Measure",labels=c("Detectable","Total")) + mytheme + xlab("Age")


wh <- Pipeline(EIR = PfPR_to_EIR_heuristic(PfPR_micro = as.numeric(MAP[which(MAP$Name=="Kinshasa"),5]),ft = admin_units_seasonal$ft[admin_units_seasonal$admin1=="Kinshasa"]),
               ft = admin_units_seasonal$ft[admin_units_seasonal$admin1=="Kinshasa"],
               years = 40,N = 10000,update_length = 30,country = "DRC",admin = "Kinshasa",use_historic_interventions = TRUE,num_het_brackets = 5,num_age_brackets = 20,
               geometric_age_brackets = TRUE,max_age = 100,use_odin = FALSE,human_only_full_save = TRUE,yearly_save = TRUE,human_yearly_save = TRUE,summary_saves_only = TRUE)

wh <- Pipeline(EIR = PfPR_to_EIR_heuristic(PfPR_micro = 0.4275,ft = admin_units_seasonal$ft[admin_units_seasonal$admin1=="Kouffo"]),
               ft = admin_units_seasonal$ft[admin_units_seasonal$admin1=="Kouffo"],years = 20,N = 1000,update_length = 10,country = "Benin",
               admin = "Kouffo",use_historic_interventions = TRUE,num_het_brackets = 5,num_age_brackets = 20,geometric_age_brackets = TRUE,
               max_age = 100,use_odin = FALSE,human_only_full_save = TRUE,yearly_save = TRUE,human_yearly_save = TRUE)
