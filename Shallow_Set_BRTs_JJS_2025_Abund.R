#example running of BRTs using Lancetfish distribution from DSLL data
library(matrixStats)
library(fmsb)
library(parallel)
library(lubridate)

setwd("~/Documents/Eco_Cast_Plus_2025")
source("BRT_Eval_Function_JJS.R")

df<-readRDS("~/Documents/LL_Data/LL_Catch_Surface_GLORYS_ESA_Chl.rds")
df<-df[df$hks_per_flt<=8 & df$year>=2005 & df$year<=2024,]
df<-df[!is.na(df$hks_per_flt),]

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
df[is.nan(df)] <- NA
df$Random<-rnorm(nrow(df))
df$bait_code<-as.factor(df$bait_code)
df$year_end<-as.factor(as.character(df$year))
df$permit_num<-as.factor(df$permit_num)
df$year<-as.factor(as.numeric(df$year))
df$log_chla<-log(df$Chl)
df$julian_end<-yday(as.Date(substr(df$set_begin_datetime, 1, 10), format="%Y-%m-%d"))
df$Full_Soak_Time_Hrs<-df$soak_time
df$lite_per_hook<-df$num_lite_devices/df$num_hks_set

df<-df[df$soak_time>=10 & df$soak_time<=40,]
df<-df[df$num_hks_set>=500 & df$num_hks_set<=2000,]



df<-df[, !names(df) %in% c("Month","Date", "Year")]

df<-df[!is.na(df$`Coryphaena hippurus`),]
All_Taxa_Names<-colnames(df)[25:192]
Abundant_Taxa<-colSums(df[, colnames(df) %in% All_Taxa_Names])
df_PA<-df[, colnames(df) %in% All_Taxa_Names]
df_PA[df_PA>0]<-1

DF_by_SPP<-colSums(df_PA, na.rm=T)/(nrow(df_PA))

SPP_Half_Perc<-which(DF_by_SPP>=0.01)
Species<-DF_by_SPP[SPP_Half_Perc]
print(Species)


Species_to_keep<-SPP_Half_Perc+24
Species_Names<-colnames(df)[Species_to_keep]


Unused_Groups<-c("Actinopterygii", "Alopiidae","Elasmobranchii", "Istiophoridae","Isurus","Scombridae", "Phoebastria nigripes", "Phoebastria immutabilis","Remora remora")
Species_to_keep<-Species_to_keep[!(Species_Names) %in% Unused_Groups]
Species_Names<-colnames(df)[Species_to_keep]
Species_Exact_Names<-Species_Names
Species_to_keep_PA<-Species_to_keep-24

Species_to_Use_SSLL_Abund<-c(Species_to_keep, 108)
Species_to_Use_SSLL_PA<-c(Species_to_keep_PA, 84)


df_Half_Perc_Abund<-df[ ,c(1:24, Species_to_Use_SSLL_Abund, 193:ncol(df))]
df_Half_Perc_PA<-cbind(df[,1:24], df_PA[,Species_to_Use_SSLL_PA], df[,193:ncol(df)])

DSLL_BRT_Models_PA<-list()


#removing permit number 

Used_preds<-c("bait_code","log_chla","julian_end", "year","lunar_rad","Random","SST","SSH","SST_SD","SSH_SD","SSS","SSS_SD","U","V","Current_Speed","MLD", "lite_per_hook")


Preds<-which(colnames(df_Half_Perc_PA) %in% Used_preds)
colnames(df_Half_Perc_PA[Preds])
DSLL_BRT_Models<-list()

Species<-c(25:50)


Species_to_Model<-list()

Species_to_Model<-list()
for ( i in 1: 26){
  Pres<-df_Half_Perc_Abund[,Species[i]]>0
  df_pres<-df_Half_Perc_Abund[Pres,]
  df_pres$Log_Abund<-log((df_pres[,Species[i]]/df_pres$num_hks_set)*1000)
  K<-ncol(df_pres)
  Preds<-which(colnames(df_pres) %in% Used_preds)
  Df_ex<- df_pres[, c(Preds, K)]
  Species_to_Model[[i]]<-Df_ex
  
}


BRT_Eval<-function(df){
  BRT_DSLL_Est_1<-gbm.step(df, gbm.x = c(1:(ncol(df)-1)), gbm.y = ncol(df),tree.complexity=5,learning.rate = 0.005 , bag.fraction = 0.75,family="gaussian", n.folds=5)
  Ntrees<-ifelse(length(BRT_DSLL_Est_1$trees)<1000, 1000,length(BRT_DSLL_Est_1$trees))
  
  BRT_DSLL_Est<-fit.brt.n_eval_Balanced_Fixed(df, gbm.x = c(1:(ncol(df)-1)), gbm.y = ncol(df),lr=0.005, tc=5,family="gaussian",nt=Ntrees, bag.fraction = 0.75,25)
  var_tested<-names(df[,c(1:(ncol(df)-1))])
  
  Abund_Model<-BRT_DSLL_Est[[1]]
  iters=length(Abund_Model)
  percent_contrib<-NULL#list()
  for(q in 1:iters){                               
    sum1<-summary(Abund_Model[q][[1]]  , plot=F )
    sum2<-sum1[order(sum1[,1], levels = var_tested),]
    percent_contrib<-cbind(percent_contrib, sum2[,2])
    rownames(percent_contrib)<-sum1[order(sum1[,1], levels = var_tested),1]
  }
  
  
  Mean_Abund_Contributions<-as.data.frame(t(rowMeans(percent_contrib)))
  
  Predictors_to_Keep_Index<-which(Mean_Abund_Contributions>Mean_Abund_Contributions$Random)
  
  Predictors_to_Keep<-Mean_Abund_Contributions[,Predictors_to_Keep_Index]
  Reduced_Predictors<-which(colnames(df) %in% colnames(Predictors_to_Keep))
  
  BRT_DSLL_Est_Red<-fit.brt.n_eval_Balanced_Fixed(df, gbm.x = Reduced_Predictors, gbm.y =  ncol(df),lr=0.005, tc=5,family="gaussian",nt=Ntrees, bag.fraction = 0.75,50)
  return(BRT_DSLL_Est_Red)
  
}

#s<-Species_to_Model[[1]]
#BRT_Eval(Species_to_Model[[1]])
detectCores()

system.time(for (i in 1:25){
DSLL_1_13_BRTs_Abund<-list()
DSLL_1_13_BRTs_Abund[[i]]<-BRT_Eval(Species_to_Model[[i]])
saveRDS(DSLL_1_13_BRTs_Abund[[i]], paste0("~/Documents/Eco_Cast_Plus_2025/SSLL_Abund_Model_v2/",Species_Exact_Names[i],"_SSLL_Abund_BRT.rds"))
rm(DSLL_1_13_BRTs_Abund)
gc()
print(paste("Completed", i, "of 25"))
})

#system.time(DSLL_1_13_BRTs_PA<-parallel::mclapply(Species_to_Model[1:13], BRT_Eval, mc.cores = 8))


#saveRDS(DSLL_1_13_BRTs_PA, "PA_DSLL_BRT_Models_1_13_2025.rds")
