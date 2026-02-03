#example running of BRTs using Lancetfish distribution from DSLL data
library(matrixStats)
library(fmsb)
library(parallel)
library(lubridate)
library(dplyr)
library(forcats)
library(ggplot2)
library(viridis)
Species_Exact_Names<-readRDS("~/Documents/Eco_Cast_Plus_2025/Species_Names_SSLL.rds")
setwd("~/Documents/Eco_Cast_Plus_2025/Partial_Plots/Abund/V2/For_Pub")
Var_Contribution_Files<-list.files("~/Documents/Eco_Cast_Plus_2025/Partial_Plots/Abund/V2/For_Pub/", pattern = "\\.csv$")
myfiles = lapply(Var_Contribution_Files, read.csv)

Species_DFs<-list()
for ( i in 1:length(myfiles)){
Species<-as.data.frame(t(myfiles[[i]]))
Species$Variables<-rownames(Species)
colnames(Species)<-c("Perc_Importance","Variables")
Species$Species_Name<-Species_Exact_Names[i]

Species_DFs[[i]]<-Species
}

Full_Var_Imp_Abund<-do.call("rbind",Species_DFs)
Full_Var_Imp_Abund<-Full_Var_Imp_Abund[Full_Var_Imp_Abund$Variables!="X",]

Mean_Var_Imp<-Full_Var_Imp_Abund %>%
  group_by(Variables) %>%
  summarize(Mean_Imp=sum(Perc_Importance, na.rm=T)/26)

Full_Var_Imp_Abund_Mean<-merge(Full_Var_Imp_Abund, Mean_Var_Imp, by=c("Variables"))


Ex<-Full_Var_Imp_Abund_Mean %>%
  mutate(Variable = fct_reorder(Variables, Mean_Imp))
setwd("~/Documents/Eco_Cast_Plus_2025/Partial_Plots/Var_Peaks/Abund/")
Var_Peak_Files<-list.files("~/Documents/Eco_Cast_Plus_2025/Partial_Plots/Var_Peaks/Abund/", pattern = "\\.csv$")
myfilespeak = lapply(Var_Peak_Files, read.csv)

Species_Peak_DFs<-list()
for ( i in 1:length(myfilespeak)){
  Species<-as.data.frame(myfilespeak[[i]])
  
  Species_Peak_DFs[[i]]<-Species
}

Full_Var_Imp_Abund_Peak<-do.call("rbind",Species_Peak_DFs)
Full_Var_Imp_Abund_Peak<-Full_Var_Imp_Abund_Peak[Full_Var_Imp_Abund_Peak$Variables!="X",]
Full_Var_Imp_Abund_Peak$Species_Name<-Full_Var_Imp_Abund_Peak$Species.Name
Var_Abund_Peaks_Full<- left_join(Full_Var_Imp_Abund, Full_Var_Imp_Abund_Peak, by=c("Species_Name","Variables"))
lunar_effects_Abund<-Var_Abund_Peaks_Full[Var_Abund_Peaks_Full$Variables=="lunar_rad",]
View(lunar_effects_Abund)
write.csv(lunar_effects_Abund, "~/Documents/Eco_Cast_Plus_2025/Manuscripts/Lunar_Var_Imp_Abund.csv")
#Ocean_vars<-Abund_Var_Imp[Abund_Var_Imp$Class=="Dynamic Oceanography",]
Full_Var_Imp_Abund_Peak<-do.call("rbind",Species_Peak_DFs)
Full_Var_Imp_Abund_Peak<-Full_Var_Imp_Abund_Peak[Full_Var_Imp_Abund_Peak$Variables!="X",]
Full_Var_Imp_Abund_Peak$Species_Name<-Full_Var_Imp_Abund_Peak$Species.Name
Var_Abund_Peaks_Full<- left_join(Full_Var_Imp_Abund, Full_Var_Imp_Abund_Peak, by=c("Species_Name","Variables"))
lunar_effects_Abund<-Var_Abund_Peaks_Full[Var_Abund_Peaks_Full$Variables=="julian_end",]
View(lunar_effects_Abund)
write.csv(lunar_effects_Abund, "~/Documents/Eco_Cast_Plus_2025/Manuscripts/DOY_Var_Imp_Abund.csv")


Ocean_Var_Range_Mins<-Full_Var_Imp_Abund_Peak %>% 
  group_by(Variables) %>%
  summarize(Range=max(Peak, na.rm=T)-min(Peak, na.rm=T), Min=min(Peak, na.rm=T), Mean=mean(Peak, na.rm=T))

Ocean_Var_Range_Mins$Range[Ocean_Var_Range_Mins$Range=="-Inf"]<-NA
Ocean_Var_Range_Mins$Min[Ocean_Var_Range_Mins$Min=="Inf"]<-NA
Ocean_Var_Range_Mins$Mean[Ocean_Var_Range_Mins$Mean=="Inf"]<-NA


Ocean_Vars_merged<-left_join(Var_Abund_Peaks_Full, Ocean_Var_Range_Mins, by=c("Variables"))

Ocean_Vars_merged$Normalized_Val<-((Ocean_Vars_merged$Peak-Ocean_Vars_merged$Min)/Ocean_Vars_merged$Range)-0.5
Ocean_Vars_merged$Norm_Mean<-(Ocean_Vars_merged$Peak-Ocean_Vars_merged$Mean)/Ocean_Vars_merged$Mean
Ocean_Vars_merged2<-merge(Ocean_Vars_merged, Mean_Var_Imp, by=c("Variables"))


Ex1<-Ocean_Vars_merged2 %>%
  mutate(Variable = fct_reorder(Variables, Mean_Imp)) 

Var_Names<-read.csv("~/Documents/Eco_Cast_Plus_2025/SSLL_Variable_Names.csv", header=TRUE)

Ex1<-merge(Ex1, Mean_Var_Imp, by=c("Variables"), all.x=TRUE)
Ex1<-merge(Ex1, Var_Names, by=c("Variables"), all.x=TRUE)


Ex1<-Ex1 %>%
  mutate(Variable_Names = fct_reorder(Variable_Name, Mean_Imp.x))
#png("Mean_Percent_Import_Dyn_Oce_Vals_v2_CPUE.png", height=8, width=8, units="in", res=300)
#dev.off()
Ex2<-Ex1[Ex1$Species_Name!="Lagocephalus lagocephalus",]

Ex2$Normalized_Val[Ex2$Variable_Names=="Lunar Cycle"]<-NA

png("~/Documents/Eco_Cast_Plus_2025/Figures/SSLL_CPUE_New_Method_Abund.png", height=8, width=9, units="in", res=300)
ggplot() +geom_point(data=Ex2 ,aes(y = Species_Name, x = Variable_Names,size=Perc_Importance,fill=Normalized_Val),pch=21)+scale_fill_gradient2(low="blue", high="red",mid="white", limits=c(-0.5,0.5), guide = guide_colourbar(title="Relative Maximum"))+xlab("Predictors")+ylab("Species")+theme_bw() +scale_x_discrete(limits=rev)+theme(axis.text.y = element_text(color="black",face="bold", size=14), axis.text.x = element_text(color="black",face="bold", angle = 45, size=12 , hjust = 1),legend.title = element_text(size = 14), legend.text = element_text(size = 14), axis.text = element_text(size = 14),      # Axis tick labels,legend.position = "none"
                                                                                                                                                                                                                                                                                                                                            axis.title = element_text(size = 16),     # Axis titles
                                                                                                                                                                                                                                                                                                                                            strip.text = element_text(size = 16))
dev.off()

T_S_Only_SSLL<-Ex2[Ex2$Variables=="SST" | Ex2$Variables=="SSS",]


T_S_Only_Wide<-reshape(T_S_Only_SSLL, idvar = "Species_Name", timevar = "Variables", direction = "wide")
library(ggrepel)
T_S_Only_Wide<-T_S_Only_Wide[!is.na(T_S_Only_Wide$Peak.SSS),]

T_S_Only_Wide$Spiciness<-swSpice(salinity=T_S_Only_Wide$Peak.SSS, temperature = T_S_Only_Wide$Peak.SST)
T_S_Only_Wide$Spiciness_Raw<-swSpice(salinity=(T_S_Only_Wide$Peak.SSS+T_S_Only_Wide$Error.SSS), temperature = (T_S_Only_Wide$Raw_Peak.SST+T_S_Only_Wide$Error.SST))

T_S_Only_Wide$Species_Name<-as.factor(T_S_Only_Wide$Species_Name)
T_S_Only_Wide1<-T_S_Only_Wide %>%
  mutate(Species_Name= fct_reorder(Species_Name, Spiciness))
T_S_Only_Wide1$Spice_Unc<-abs(T_S_Only_Wide1$Spiciness-T_S_Only_Wide1$Spiciness_Raw)

png("~/Documents/Eco_Cast_Plus_2025/Figures/Spiciness_Space_CPUE_New_Method_SSLL_Abund.png", height=7, width=8, units="in", res=300)
ggplot(data=T_S_Only_Wide1 ,aes(y =Species_Name , x = Spiciness)) +scale_y_discrete(limits=rev)+geom_point(fill="black",pch=21)+geom_errorbar(aes(xmin = Spiciness-Spice_Unc, xmax = Spiciness+Spice_Unc), width = 0.2)+xlab(expression(paste("Spiciness (kg ", m^{-3},")")))+ylab("Species")+theme_bw()+theme(axis.text.y = element_text(color="black",face="bold", size=14), axis.text.x = element_text(color="black",face="bold", size=12),legend.title = element_text(size = 14), legend.text = element_text(size = 14), axis.text = element_text(size = 14),      # Axis tick labels,legend.position = "none"
                                                                                                                                                                                                                                                                                                               axis.title = element_text(size = 16),     # Axis titles
                                                                                                                                                                                                                                                                                                               strip.text = element_text(size = 16))

dev.off() 













Ex1<-Ocean_Vars_merged2 %>%
  mutate(Variable = fct_reorder(Variables, Mean_Imp))

png("Mean_Percent_Import_Dyn_Oce_Vals_v2.png", height=8, width=8, units="in", res=300)
ggplot(data = Ex1) +
  aes(y = Variable, x = Perc_Importance, color=Normalized_Val) +geom_point(aes(fill=Normalized_Val),colour="black",pch=21, size=5)+scale_fill_gradient2(low="blue", high="red",mid="white", limits=c(-0.5,0.5), guide = guide_colourbar(title="Relative Maximum"))+xlab("Mean Perc. Importance")
dev.off()


ggplot(data = Ocean_Vars_merged) +
  aes(y = Variables, x = Normalized_Val, color=Perc_Importance) +geom_point()+scale_color_viridis()


Lunar_Fish<-Ex[Ex$Variable=="lunar_rad" & Ex$Perc_Importance>=3,]
  
  