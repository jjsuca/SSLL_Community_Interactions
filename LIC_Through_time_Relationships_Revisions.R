#Predicting certain species SDMs over the time series
library(matrixStats)
library(fmsb)
library(parallel)
library(lubridate)
library(DescTools)
library(ggplot2)  
library(viridis)
library(lubridate)
library(tidyr)
library(trend)
Env_Vars_SSLL_Quarter_Deg<-readRDS("~/Documents/LL_Data/Weekly_Quart_Deg_SSLL_GLORYS_Vars_Chl.rds")

setwd("~/Documents/Eco_Cast_Plus_2025")
source("BRT_Eval_Function_JJS.R")
Env_Vars_SSLL_Quarter_Deg<-readRDS("~/Documents/LL_Data/Weekly_Quart_Deg_SSLL_GLORYS_Vars_Chl.rds")

Used_preds<-c("bait_code","log_chla","julian_end", "year","lunar_rad","SST","SSH","SST_SD","SSH_SD","SSS","SSS_SD","U","V","Current_Speed","MLD", "lite_per_hook")

Preds<-which(colnames(Env_Vars_SSLL_Quarter_Deg) %in% Used_preds)
colnames(Env_Vars_SSLL_Quarter_Deg[Preds])
#read in models

Model_Names<-list.files(path="~/Documents/Eco_Cast_Plus_2025/SSLL_PA_Model")
Model_Fit_Metrics<-read.csv("~/Documents/Eco_Cast_Plus_2025/Partial_Plots/Abund/V2/Abund_Metrics_SSLL_Response_Update_2025_tc5.csv", header=TRUE)

Species_with_Good_Fit<-Model_Fit_Metrics$X[Model_Fit_Metrics$Mean_R2>=0.2]

Relevant_Species<-c("Caretta caretta","Carcharhinus longimanus","Coryphaena hippurus", "Isurus oxyrinchus","Kajikia audax", "Lampris guttatus", "Lepidocybium flavobrunneum" ,"Prionace glauca", "Pteroplatytrygon violacea","Ruvettus pretiosus", "Thunnus alalunga","Thunnus albacares","Thunnus obesus","Xiphias gladius") 

LIC_Species<-readRDS( "~/Documents/Eco_Cast_Plus_2025/Spatial_Correlations/LIC_Swrd_Sp_SSLL.rds")      
LIC_Species_Matix<-t(LIC_Species[,2,])
LIC_Species_Matix<-LIC_Species_Matix[1:13,]
rownames(LIC_Species_Matix)<-Relevant_Species[1:13]
LIC_Species<-as.data.frame(t(LIC_Species_Matix))
LIC_Species_Stacked<-stack(LIC_Species)
Dates<-unique(Env_Vars_SSLL_Quarter_Deg$Date_for_Pairing)
LIC_Species_Stacked$Date<-Dates
LIC_Species_Stacked$Month<-month(LIC_Species_Stacked$Date)
LIC_Species_Stacked$Year<-year(LIC_Species_Stacked$Date)

LIC_Species_Stacked$Quarter<-NA
LIC_Species_Stacked$Quarter[LIC_Species_Stacked$Month<4]<-1
LIC_Species_Stacked$Quarter[LIC_Species_Stacked$Month>=4 & LIC_Species_Stacked$Month<=6]<-2

LIC_Species_Stacked$Quarter[LIC_Species_Stacked$Month>=10]<-4
LIC_Species_Stacked$Quarter<-as.factor(as.character(LIC_Species_Stacked$Quarter))
LIC_Species_Stacked<-LIC_Species_Stacked[!is.na(LIC_Species_Stacked$Quarter),]
Species_Names<-readxl::read_excel("~/Documents/Eco_Cast_Plus_2025/DSLL/Data_Outputs/Model_Fits/Model_Skill_Table_DSLL.xlsx")
Ex3<-merge(LIC_Species_Stacked, Species_Names, by.x=c("ind"), by.y=c("Binomical Nomen."), all.x=TRUE)


png("~/Documents/Eco_Cast_Plus_2025/Figures/Revisions/LIC_By_Species_Good_Spp_Updated_Points_for_Pub_Revision.png", height=6, width=9, units="in", res=300)
ggplot(
  
  Ex3,
  
  aes(
    
    x = values,
    
    y = English_Name,
    
    fill = Year
    
  )
  
) +
  
  geom_jitter(
    
    shape = 21,
    
    color = "black",
    
    width = 0.015,
    
    height = 0.30,
    
    size = 1.25,
    
    alpha = 0.7
    
  )+
  
  facet_wrap(
    
    ~ Quarter,
    
    nrow = 1,
    
    labeller = labeller(
      
      Quarter = c(
        
        "1" = "Q1 (Jan–Mar)",
        
        "2" = "Q2 (Apr–Jun)",
        
        "4" = "Q4 (Oct–Dec)"
        
      )
      
    )
    
  ) +
  
  scale_y_discrete(limits = rev) +
  
  scale_fill_viridis_c(
    
    option = "plasma",
    
    name = "Year"
    
  ) +
  
  theme_bw() +
  
  theme(
    
    axis.text.y = element_text(
      
      face = "bold",
      
      color = "black"
      
    ),
    
    axis.title.x = element_text(
      
      face = "bold",
      
      color = "black"
      
    ),
    
    axis.text.x = element_text(
      
      face = "bold",
      
      color = "black"
      
    ),
    
    strip.text = element_text(
      
      face = "bold",
      
      size = 12
      
    ),
    
    legend.background = element_rect(fill = "white")
    
  ) +
  
  xlab("Local Index of Collocation") +
  
  ylab("")
dev.off()

ggplot(LIC_Species_Stacked, aes(y = values, x = Date,colour=ind)) +
  geom_point()+ theme_bw()+scale_color_viridis(discrete=TRUE,option="plasma")+ theme(axis.text.y= element_text(face="bold.italic", color="black"), axis.title.x =element_text(face="bold", color="black"), axis.text.x= element_text(face="bold", color="black"), legend.background = element_rect(fill="grey") )+
  ylab("")+xlab("Local Index of Collocation")


LIC_Species_Stacked$Year<-year(LIC_Species_Stacked$Date)

LIC_Annual_SPP_SSLL<-LIC_Species_Stacked %>%
  dplyr::group_by(Year, ind) %>%
  dplyr::summarize(LIC=mean(values), LIC_SD=sd(values))

ggplot(LIC_Annual_SPP_SSLL, aes(y = LIC, x = Year,colour=ind)) +
  geom_point()+ geom_path()+theme_bw()+scale_color_viridis(discrete=TRUE)+ theme(axis.text.y= element_text(face="bold.italic", color="black"), axis.title.x =element_text(face="bold", color="black"), axis.text.x= element_text(face="bold", color="black"), legend.background = element_rect(fill="grey") )+
  ylab("")+xlab("Local Index of Collocation")


PDO<-read.csv("~/Documents/Figure_Code/Data/PDO.csv", header=TRUE)
#Uku<-read.csv("Uku_Recruitment_2023.csv", header=TRUE)


PDO_Long<-PDO %>%
  pivot_longer(!Year, names_to = "Month", values_to = "PDO")
PDO_Long<-PDO_Long %>%
  mutate(Month = recode(Month,
                        January = 1,
                        February = 2,
                        March = 3,
                        April = 4,
                        May = 5,
                        June = 6,
                        July = 7,
                        August = 8,
                        September = 9,
                        October = 10,
                        November = 11,
                        December = 12
  ))



PDO_Fishing<-PDO_Long[PDO_Long$Month<6 | PDO_Long$Month>=11,]

PDO_Fishing_Annual<- PDO_Fishing %>%
  dplyr::group_by(Year) %>%
  dplyr:: summarise(PDO=mean(PDO))

NPGO<-read.csv("~/Documents/Figure_Code/Data/NPGO.csv", header=TRUE)
#Uku<-read.csv("Uku_Recruitment_2023.csv", header=TRUE)

NPGO_Long<-NPGO %>%
  pivot_longer(!Year, names_to = "Month", values_to = "NPGO")
NPGO_Long<-NPGO_Long %>%
  mutate(Month = recode(Month,
                        January = 1,
                        February = 2,
                        March = 3,
                        April = 4,
                        May = 5,
                        June = 6,
                        July = 7,
                        August = 8,
                        September = 9,
                        October = 10,
                        November = 11,
                        December = 12
  ))



NPGO_Fishing<-NPGO_Long[NPGO_Long$Month<6 | NPGO_Long$Month>=11,]

NPGO_Fishing_Annual<- NPGO_Fishing %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(NPGO=mean(NPGO))

ENSO<-read.csv("~/Documents/Alalaua/meiv2.csv", header=TRUE)
#Uku<-read.csv("Uku_Recruitment_2023.csv", header=TRUE)


ENSO_Long<-ENSO %>%
  pivot_longer(!Year, names_to = "Month", values_to = "ENSO")
ENSO_Long<-ENSO_Long %>%
  mutate(Month = recode(Month,
                        January = 1,
                        February = 2,
                        March = 3,
                        April = 4,
                        May = 5,
                        June = 6,
                        July = 7,
                        August = 8,
                        September = 9,
                        October = 10,
                        November = 11,
                        December = 12
  ))

ENSO_Fishing<-ENSO_Long[ENSO_Long$Month<6 | ENSO_Long$Month>=11,]

ENSO_Fishing_Annual<- ENSO_Fishing %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(ENSO=mean(ENSO))

Climate_Indices<-
  left_join(LIC_Annual_SPP_SSLL, NPGO_Fishing_Annual, by='Year') %>%
  left_join(., PDO_Fishing_Annual, by='Year') %>%
  left_join(., ENSO_Fishing_Annual, by='Year') 
  
  
Cors_NPGO <- Climate_Indices |> dplyr::group_by(ind) |> dplyr::summarize(Cor=round(cor(NPGO,LIC),3), p=round(cor.test(NPGO,LIC)$p.value,4))
Cors_PDO <- Climate_Indices |> dplyr::group_by(ind) |> dplyr::summarize(Cor=round(cor(PDO,LIC),3),p=round(cor.test(PDO,LIC)$p.value,4))
Cors_ENSO <- Climate_Indices |> dplyr::group_by(ind) |> dplyr::summarize(Cor=round(cor(ENSO,LIC),3),p=round(cor.test(ENSO,LIC)$p.value,4))


MK_Time <- Climate_Indices |> dplyr::group_by(ind) |> dplyr::summarize(Tau=mk.test(LIC)$estimates[3], p=mk.test(LIC)$p.value)

ggplot(Climate_Indices, aes(NPGO, LIC)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ ind)

Sig_Spp_NPGO<-Cors_NPGO$ind[Cors_NPGO$p<0.01]
Sig_Spp_PDO<-Cors_PDO$ind[Cors_PDO$p<0.01]
Sig_Spp_ENSO<-Cors_ENSO$ind[Cors_ENSO$p<0.01]
Sig_Spp_Time<-Cors_Time$ind[Cors_Time$p<0.01]

LIC_Annual_SPP_SSLL_Time <- Climate_Indices %>%
  
  filter(ind %in% c("Thunnus albacares",
                    
                    "Thunnus obesus",
                    
                    "Isurus oxyrinchus",
                    
                    "Ruvettus pretiosus")) %>%
  
  mutate(ind = factor(ind,
                      
                      levels = c("Isurus oxyrinchus",
                                 
                                 "Ruvettus pretiosus",
                                 
                                 "Thunnus obesus",
                                 
                                 "Thunnus albacares")))


png("~/Documents/Eco_Cast_Plus_2025/Figures/LIC_By_Species_w_Trends_Update_Non_viridis.png", height=5.5, width=11, units="in", res=300)
ggplot(LIC_Annual_SPP_SSLL_Time,
       
       aes(x = Year, y = LIC, colour = ind)) +
  
  geom_point() +
  
  geom_path() +
  
  geom_errorbar(aes(ymin = LIC - LIC_SD,
                    
                    ymax = LIC + LIC_SD),
                
                width = 0.2) +
  
  theme_bw() +
  
  scale_color_manual(
      values = c("#0074D9", "#B10DC9", "#85144b", "#FF4136"),
    
    labels = c(
      
      "Shortfin mako",
      
      "Oilfish",
      
      "Bigeye tuna",
      
      "Yellowfin tuna"
      
    )
    
  ) +
  
  xlab("Year") +
  
  ylab("Local Index of Collocation") +
  
  labs(colour = "Species") +
  
  theme(
    
    legend.title = element_text(size = 14),
    
    legend.text = element_text(size = 14),
    
    axis.text = element_text(size = 14),
    
    axis.title = element_text(size = 16),
    
    strip.text = element_text(size = 16)
    
  )


dev.off()


ggplot(LIC_Annual_SPP_SSLL_NPGO, aes(Year, LIC)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ ind, scales = "free",)
  
  