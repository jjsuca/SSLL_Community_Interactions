#example running of BRTs using Lancetfish distribution from DSLL data
library(matrixStats)
library(fmsb)
library(parallel)
library(lubridate)

setwd("~/Documents/Eco_Cast_Plus_2025")
source("BRT_Eval_Function_JJS.R")

df<-readRDS("~/Documents/LL_Data/LL_Catch_Surface_GLORYS_ESA_Chl_Rugosity.rds")
df<-df[df$hks_per_flt<=8 & df$year>=2005 & df$year<=2024,]
df<-df[!is.na(df$hks_per_flt),]

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
df[is.nan(df)] <- NA
df$Random<-rnorm(nrow(df))
df$bait_code<-as.factor(df$bait_code)
df$year_end<-as.factor(as.character(df$year))
df$permit_num<-as.factor(df$permit_num)
df$year<-as.numeric(df$year)
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

Species_Names<-colnames(df_PA)[Species_to_Use_SSLL_PA]

df_Half_Perc_Abund<-df[ ,c(1:24, Species_to_Use_SSLL_Abund, 193:ncol(df))]
df_Half_Perc_PA<-cbind(df[,1:24], df_PA[,Species_to_Use_SSLL_PA], df[,193:ncol(df)])


# Used_preds<-c("HKS","permit_num","bait_code", "sst","log_chla","julian_end", "year_end","lunar_rad","sst_anom","SLA","Eddy_Classification","current.divergence","Random","PDO","MEI",  "bathymetry","EKE", "model_sss","chla.anom","mld","NPGO", "FLTLN_LEN","BRNCHLN_LEN", "HKS_FLT", "LITE_HK","Full_Soak_Time_Hrs")
# 
# 
# Preds<-which(colnames(df_Half_Perc_PA) %in% Used_preds)
# SSLL_BRT_Models<-list()
# 
Species_Exact_Names<-Species_Names

Model_Names<-list.files("~/Documents/Eco_Cast_Plus_2025/SSLL_Abund_Model_v2")
#
#Leatherback_Model<-Model_Names[14]
#Red_Model_Names<-Model_Names[-14]
#Model_Names<-c(Red_Model_Names, Leatherback_Model)

Model_Skill_Abund<-matrix(,length(Model_Names),4)#
setwd("~/Documents/Eco_Cast_Plus_2025/Figures/PDP_Revisions/Abund/")
for (p in 1:length(Model_Names)){#
  SSLL_Abund_Models<-readRDS(paste0("~/Documents/Eco_Cast_Plus_2025/SSLL_Abund_Model_v2/",Model_Names[p]))
 Model_Abund_Eval<-matrix(,50,2)
 if (is.null(SSLL_Abund_Models)==FALSE){
 Model_Evals_Abund<-unlist(unlist(SSLL_Abund_Models[[2]]))
   Model_Abund_Eval <- matrix(Model_Evals_Abund, nrow = 2, byrow = TRUE)
 Model_Skill_Abund[p,1]<-mean( Model_Abund_Eval[,1])
 Model_Skill_Abund[p,2]<-min( Model_Abund_Eval[,1])
 Model_Skill_Abund[p,3]<-mean( Model_Abund_Eval[,2])
 Model_Skill_Abund[p,4]<-min( Model_Abund_Eval[,2])
 }
 else if (is.null(SSLL_Abund_Models)==TRUE){
   Model_Skill_Abund[p,1]<-NA
   Model_Skill_Abund[p,2]<-NA
   Model_Skill_Abund[p,3]<-NA
   Model_Skill_Abund[p,4]<-NA
}
rownames(Model_Skill_Abund)<-Species_Exact_Names[1:25] #
colnames(Model_Skill_Abund)<-c("Mean_R2","Min_R2","Mean_RMSE","Max_RMSE")
#write.csv(Model_Skill_Abund,"Abund_Metrics_SSLL_Response_Update_240628_1_13.csv")



# #####Full model plot#########
`%nin%` = Negate(`%in%`)
   Model<-SSLL_Abund_Models
   if (is.null(Model)==FALSE){
 Fish_Models<-SSLL_Abund_Models[[1]]
 Species_Name<-Species_Exact_Names[p]
 var_tested<-Fish_Models[[1]]$var.names
 #var_tested<-colnames(Abund_Data_Reef_Fishes_MHI[Predictors])
 Fish_Models_Good<-Fish_Models#[Q]
 percent_contrib<-NULL#list()
 iters=length(Fish_Models_Good)
 part_plot<-list()
 part_plot<-list()
 percent_contrib<-NULL#list()
 Continuous_Preds<-which(var_tested %nin% c("year","permit_num","bait_code", "Eddy_Classification","HKS_FLT"))
 for(q in 1:iters){                                #this was 50
 mod<-Fish_Models_Good[q][[1]]
   ###
   part_plot1<-data.frame(row.names=1:100)
   for(x in c(Continuous_Preds)){ ###

       pp<-plot(mod ,var_tested[x],return.grid=T) ###
       part_plot1<-cbind(part_plot1, pp) ###
     }

#   ###
   part_plot[[q]]<-part_plot1 ###

   sum1<-summary(Fish_Models_Good[q][[1]]  , plot=F )
   sum2<-sum1[order(sum1[,1], levels = var_tested),]
   percent_contrib<-cbind(percent_contrib, sum2[,2])
   rownames(percent_contrib)<-sum1[order(sum1[,1], levels = var_tested),1]
 }
 All_percent_contribution<-cbind(rownames(percent_contrib), paste(round(rowMeans(percent_contrib),2), round(rowSds(percent_contrib),2), sep=" ± "))
 Combined_All_percent_contribution<-All_percent_contribution
#
#
 Mean_Abund_Contributions<-as.data.frame(t(rowMeans(percent_contrib)))
 write.csv( Mean_Abund_Contributions,paste0("Var_Contributions_",Species_Name,"_Abund_SSLL_2025_tc5v2.csv"))

 Abund_Predictors_Plot<- rbind(rep(max(Mean_Abund_Contributions),length(var_tested)) , rep(0,length(var_tested)) , Mean_Abund_Contributions)
 Abund_Predictors_Plot[]<-sapply(Abund_Predictors_Plot, as.numeric)
 par(mfrow=c(1,1))

 # png(paste0("Radar_Chart_",Species_Name,"_Abund_SSLL_2025_tc5v2.png"), height=6, width=6, units="in",res=300)
 #radarchart(Abund_Predictors_Plot,  pfcol=rgb(0.0,0.3,0.5,0.5), pcol=rgb(0.0,0.3,0.5,0.5), title=paste0(Species_Name,"_Abund"))
  #dev.off()
#
 All_percent_contribution<-cbind(rownames(percent_contrib), paste(round(rowMeans(percent_contrib),2), round(rowSds(percent_contrib),2), sep=" ± "))
#
 png(paste0("Partial_plots_",Species_Name,"_Abund_SSLL_2025_tc5_Revision.png"), height=18,width=14, res=300, units="in")
 par(mfrow=c(5,5))
 mn_part_plot<-list()
 for(y in c(Continuous_Preds)){
   id<-which(colnames(part_plot[[1]])==var_tested[y])
   all1<-NULL
   all2<-NULL
   for(z in 1:iters){											 #this was 50
     all1<-rbind(all1, cbind(c(part_plot[[z]][,id])))
     all2<-rbind(all2, cbind(c(part_plot[[z]][,id+1])))
   }
   all3 <- cbind(all1, all2)
   all1 <- all3[order(all3[,1]),]
   
   # Common x grid
   xgrid <- seq(min(all1[,1], na.rm=TRUE),
                max(all1[,1], na.rm=TRUE),
                length.out=200)
   
   # Interpolate each model's partial dependence onto common grid
   pred_mat <- matrix(NA, nrow=length(xgrid), ncol=iters)
   
   for(z in 1:iters){
     
     tmp_x <- part_plot[[z]][,id]
     tmp_y <- part_plot[[z]][,id+1]
     
     o <- order(tmp_x)
     
     pred_mat[,z] <- approx(
       x = tmp_x[o],
       y = tmp_y[o],
       xout = xgrid,
       rule = 2
     )$y
   }
   
   # Mean response across models
   mn <- rowMeans(pred_mat, na.rm=TRUE)
   
   # Envelope across models
   lower <- apply(pred_mat, 1, min, na.rm=TRUE)
   upper <- apply(pred_mat, 1, max, na.rm=TRUE)
   
   # Plot limits
   ylim_use <- range(c(lower, upper), na.rm=TRUE)
   
   plot(xgrid, mn,
        type="n",
        ylim=ylim_use,
        xlab=var_tested[y],
        ylab=paste0("f(", var_tested[y], ")"),
        cex.axis=1.2,
        cex.lab=1.2)
   
   # Shaded model range
   polygon(
     c(xgrid, rev(xgrid)),
     c(lower, rev(upper)),
     col=adjustcolor("grey70", alpha.f=0.5),
     border=NA
   )
   
   # Mean smooth line
   lines(xgrid, mn, lwd=2)
   
   rug(na.omit(unlist(df_Half_Perc_PA[var_tested[y]])))
   
   legend(
     "bottomright",
     paste(
       All_percent_contribution[
         which(All_percent_contribution[,1]==var_tested[y]),2
       ],
       "%",
       sep=" "
     ),
     bty="n",
     cex=1.4
   )}
dev.off()
 rm(Fish_Models)
 rm(Model)
 gc()
 }
rm(SSLL_Abund_Models)}

#write.csv(Model_Skill_Abund,"Abund_Metrics_SSLL_Response_Update_2025_tc5v2.csv")
