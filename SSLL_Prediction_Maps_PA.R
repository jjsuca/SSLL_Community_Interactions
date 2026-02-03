#Predicting certain species SDMs over the time series
library(matrixStats)
library(fmsb)
library(parallel)
library(lubridate)
library(DescTools)

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

matches <- sapply(Relevant_Species, function(pattern) grepl(pattern, Model_Names))

# Identify which main_vec elements contain any search_vec value
matching_indices <- rowSums(matches) > 0

# Get the matching elements
Relevant_Models <- Model_Names[matching_indices]


for (i in 1:length(Relevant_Models)){
  Model_Estimates<-matrix(, nrow=nrow(Env_Vars_SSLL_Quarter_Deg), ncol=50)
  SSLL_PA_Models<-readRDS(paste0("~/Documents/Eco_Cast_Plus_2025/SSLL_PA_Model/",Relevant_Models[i]))
  
  for (k in 1:50){
    Model_Estimates[,k]<-predict.gbm(SSLL_PA_Models[[1]][[k]], Env_Vars_SSLL_Quarter_Deg,
                                     n.trees=SSLL_PA_Models[[1]][[k]]$gbm.call$best.trees, type="response")
    print(paste("completed", k, "of 50"))}
  saveRDS(Model_Estimates, paste0("~/Documents/Eco_Cast_Plus_2025/Model_Predictions/PA/",Relevant_Species[i],"_SSLL_PA_Model_Est.rds"))
  print(paste("completed", Relevant_Species[i],i, "of", length(Relevant_Models)))
rm(Model_Estimates)
  
}


