#install.packages('gamlss')
rm(list=ls())
library(glmnet)
library(dplyr);
library(ggplot2)
#library(ggseg)
# install.packages("ggplot2")
library(ggplot2)
library(stringr)
library(robustbase)
# setwd('C:/ZZZ/work/manuscript/Lifespan-main/test_zzz/aseg.vol.table')
#install.packages('gamlss')
library(gamlss)
library(ggplot2)
library(reshape2)


datapath='/media/ps/Seagate Basic/Lifespan_freesurfer_results/Github/Source-codes'# change
setwd(datapath)#set the filepath



source("100.common-variables.r")
source("101.common-functions.r")

source("ZZZ_function.R")
source("300.variables.r")
source("301.functions.r")

library(raincloudplots)
library(commonmark)
library(ggplot2)
#library(ggstatsplot)
library(lmerTest)
library(effectsize)
library(ggrain)
library(ggsci)
library(ggdist)
library(hrbrthemes)
library(lsr)
library(proxy)


datapath<-'/media/ps/Seagate Basic/Lifespan_freesurfer_results/Results/Step1_normative_models/V7_DK_updated_20260512'
savepath=paste0('/media/ps/Seagate Basic/Lifespan_freesurfer_results/HC_subtypes/results')

path<-list()

str_lab<-'global_and_local_feature'# change
path[[1]]<-paste0(datapath,'/aseg.vol.table');
path[[2]]<-paste0(datapath,'/lh.aparc.thickness.table');
path[[3]]<-paste0(datapath,'/rh.aparc.thickness.table');
path[[4]]<-paste0(datapath,'/lh.aparc.area.table');
path[[5]]<-paste0(datapath,'/rh.aparc.area.table');
path[[6]]<-paste0(datapath,'/lh.aparc.volume.table');
path[[7]]<-paste0(datapath,'/rh.aparc.volume.table');
path[[8]]<-paste0(datapath,'/Global_feature');



Exd_var<-c("Global_feature_meanCT2_lhMeanThickness_loop_our_model.rds",
           "Global_feature_meanCT2_lhVertex_loop_our_model.rds" ,
           "Global_feature_meanCT2_rhMeanThickness_loop_our_model.rds",
           "Global_feature_meanCT2_rhVertex_loop_our_model.rds",
           "Global_feature_totalSA2_lh_loop_our_model.rds",
           "Global_feature_totalSA2_rh_loop_our_model.rds",
           'aseg.vol.table_EstimatedTotalIntraCranialVol_loop_our_model.rds',
           "aseg.vol.table_cerebellum_WM_loop_our_model.rds" ,
           #"aseg.vol.table_cerebellum_total_loop_our_model.rds",
           "aseg.vol.table_cerebellum_GM_loop_our_model.rds",
           # "aseg.vol.table_Brain-Stem_loop_our_model.rds",
           "lh.aparc.thickness.table_eTIV_loop_our_model.rds",
           "rh.aparc.thickness.table_eTIV_loop_our_model.rds",
           "lh.aparc.area.table_eTIV_loop_our_model.rds",
           "rh.aparc.area.table_eTIV_loop_our_model.rds",
           "lh.aparc.volume.table_eTIV_loop_our_model.rds",
           "rh.aparc.volume.table_eTIV_loop_our_model.rds")


RDSfile_all<-NULL 
for(i in 1:length(path))
{
  feature_path=path[[i]];
  setwd(feature_path)
  myfile <- list.files()
  RDSfile_all0 <- myfile[grep(myfile,pattern ="our_model.rds$")]
  for(j in RDSfile_all0)
  {
    if(!(j %in% Exd_var))
    {
      #RDSfile_all0<-paste0(feature_path,'/',RDSfile_all0);
      RDSfile_all<-c(RDSfile_all,paste0(feature_path,'/',j));
    }
  }
}


feature_col<-NULL
results<-readRDS(RDSfile_all[1]);
#data<-results$all_data_original
data<-results$all_data_original[,c(1:16,56)]
for(RDSfile_loop in RDSfile_all)
{
  
  if(!(RDSfile_loop %in% Exd_var)){
    print(RDSfile_loop)
    RDSfile<-RDSfile_loop
    results<-readRDS(RDSfile)
    tem_data1<-results$Quant_data[results$i]
    tem_data1<-data.frame(tem_data1);
    
    all_data<-results$all_data
    all_data<-all_data[!is.na(tem_data1$Quant_score)&!is.infinite(tem_data1$Quant_score),]
    
    
    tem_data1<-tem_data1[!is.na(tem_data1$Quant_score)&!is.infinite(tem_data1$Quant_score),]
    tem_data1<-data.frame(tem_data1);
    
    
    rownames(tem_data1)<-rownames(results$all_data);
    index<-intersect(rownames(tem_data1),rownames(data));
    data<-cbind(data[index,],tem_data1[index,])
    colnames(data)[dim(data)[2]]<-results$i
    
    feature_col<-c(feature_col,results$i)
    
    print(dim(data))
  }
}  



HC_data<-data[data$Diagnosis=='HC'&data$Country=='China',] #Chinese populations

start_col=18;
end_col=dim(HC_data)[2]

df_scaled <- HC_data[,c(start_col:end_col)]



data_clean <- df_scaled[complete.cases(df_scaled), ]
data_clean <- data_clean[, sapply(data_clean, function(x) length(unique(x)) > 1)]

dist_mat <- proxy::dist(data_clean, method = "euclidean")


hc <- hclust(dist_mat, method = "ward.D2")


plot(hc, main = "Dendrogram of hierarchical clustering for the data", xlab = "samples", sub = "", cex = 0.6)


#Elbow Method）
wss <- function(hc, max_k = 10) {
  sapply(2:max_k, function(k) {
    clusters <- cutree(hc, k)
    sum(sapply(unique(clusters), function(c) {
      idx <- which(clusters == c)
      if(length(idx) <= 1) return(0)
      sum(dist(df_scaled[idx, ])^2) / length(idx)
    }))
  })
}

wss_values <- wss(hc, max_k = 10)
plot(2:10, wss_values, type = "b", pch = 19,
     xlab = "cluster k", ylab = "within-cluster sum of squares")



cluster_labels <- cutree(hc, k = 4)
table(cluster_labels)  


cluster_labels<-data.frame(cluster_labels)


ind<-intersect(rownames(data),rownames(cluster_labels));

data<-cbind(data[ind,],cluster_labels[ind,]);

colnames(data)[dim(data)[2]]<-'Cluster'

setwd(savepath);
saveRDS(data,'HC_data_cluster_euclidean.rds')

