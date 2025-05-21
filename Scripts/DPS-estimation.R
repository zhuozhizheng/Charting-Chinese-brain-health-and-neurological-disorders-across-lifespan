#install.packages('gamlss')
rm(list=ls())


#set your own directory
datapath='***/Source-codes/'#Change the directory where you save the Source-codes 
feature_path0='**/Diseases'; #Chnage the path you have saved your normative models in "Disease-application-normative-model.R"
savepath='**/Results/DPS_estimation';# Create and determine the directory where you would save the results  

library(glmnet)
library(dplyr);
library(ggplot2)
library(stringr)
library(robustbase)
library(gamlss)
library(reshape2)
library(raincloudplots)
library(commonmark)
library(lmerTest)
library(effectsize)
library(ggrain)
library(ggsci)
library(ggdist)
library(hrbrthemes)

setwd(datapath)#set the filepath
source("100.common-variables.r")
source("101.common-functions.r")

source("ZZZ_function.R")
source("300.variables.r")
source("301.functions.r")



path<-list()

str_lab<-'global_and_local_feature_all_for_ref_group_sample_no_site_match'# change
str_feature<-c('centile');


path[[1]]<-paste0(feature_path0,'aseg.vol.table');
path[[2]]<-paste0(feature_path0,'lh.aparc.thickness.table');
path[[3]]<-paste0(feature_path0,'rh.aparc.thickness.table');
path[[4]]<-paste0(feature_path0,'lh.aparc.area.table');
path[[5]]<-paste0(feature_path0,'rh.aparc.area.table');
path[[6]]<-paste0(feature_path0,'lh.aparc.volume.table');
path[[7]]<-paste0(feature_path0,'rh.aparc.volume.table');
path[[8]]<-paste0(feature_path0,'Global_feature');




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
           #"aseg.vol.table_Brain-Stem_loop_our_model.rds",
           "lh.aparc.thickness.table_eTIV_loop_our_model.rds",
           "rh.aparc.thickness.table_eTIV_loop_our_model.rds",
           "lh.aparc.area.table_eTIV_loop_our_model.rds",
           "rh.aparc.area.table_eTIV_loop_our_model.rds",
           "lh.aparc.volume.table_eTIV_loop_our_model.rds",
           "rh.aparc.volume.table_eTIV_loop_our_model.rds",
           
           "Global_feature_meanCT2_lhMeanThickness_loop_our_model.rds",
           "Global_feature_meanCT2_lhVertex_loop_our_model.rds" ,
           "Global_feature_meanCT2_rhMeanThickness_loop_our_model.rds",
           "Global_feature_meanCT2_rhVertex_loop_our_model.rds",
           "Global_feature_totalSA2_lh_loop_our_model.rds",
           "Global_feature_totalSA2_rh_loop_our_model.rds",
           'aseg.vol.table_EstimatedTotalIntraCranialVol_loop_our_model.rds',
           "aseg.vol.table_cerebellum_WM_loop_our_model.rds" ,
           #"aseg.vol.table_cerebellum_total_loop_our_model.rds",
           "aseg.vol.table_cerebellum_GM_loop_our_model.rds",
           #"aseg.vol.table_Brain-Stem_loop_our_model.rds",
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



results<-readRDS(RDSfile_all[1]);
#data<-results$all_data_original
data<-results$all_data_original[,c('Age','Sex','Site_ZZZ','Diagnosis')]

start_col=dim(data)[2]+1
str<-str_lab;

for(RDSfile_loop in RDSfile_all)
{
  
  if(!(RDSfile_loop %in% Exd_var)){
  print(RDSfile_loop)
  RDSfile<-RDSfile_loop
  results<-readRDS(RDSfile)
  
  
  if(!(str_feature %in% c('zscore','raw')))
  {
    print('centile')
    tem_data1<-results$Quant_data[results$i]# centile
    tem_data1<-data.frame(tem_data1);
    tem_data1[is.infinite(tem_data1[,1])|is.na(tem_data1[,1]),1]<-mean(tem_data1[!(is.infinite(tem_data1[,1])|is.na(tem_data1[,1])),1])
  }
  
  if(str_feature %in% c('zscore'))
  {
    print('zscore')
    tem_data1<-results$Zscore[results$i]# zscore
    tem_data1<-data.frame(tem_data1);
    tem_data1[is.infinite(tem_data1[,1])|is.na(tem_data1[,1]),1]<-mean(tem_data1[!(is.infinite(tem_data1[,1])|is.na(tem_data1[,1])),1])
  }
  
  
  if(str_feature %in% c('raw'))
  {
    print('raw')
    tem_data1<-results$all_data["tem_feature"]# raw
    tem_data1<-data.frame(tem_data1);
    tem_data1[is.infinite(tem_data1[,1])|is.na(tem_data1[,1]),1]<-mean(tem_data1[!(is.infinite(tem_data1[,1])|is.na(tem_data1[,1])),1])
  }
  
  
  
  tem_data1<-data.frame(tem_data1);
  index<-intersect(rownames(tem_data1),rownames(data));
  data<-cbind(data[index,],tem_data1[index,])
  colnames(data)[dim(data)[2]]<-results$i
  }
}  

tem_colnames<-colnames(data)[start_col:dim(data)[2]];
tem_colnames[tem_colnames=='total.surface.arrea']<-'Surface.area';
tem_colnames[tem_colnames=='Ventricles']<-'CSF';
tem_colnames[tem_colnames=='cerebellum.total']<-'Cerebellum';
tem_colnames[tem_colnames=='Brain.Stem']<-'Brainstem';
tem_colnames[tem_colnames=='mean.thickness']<-'Cortical.thickness';


for(i_sub in 1:length(tem_colnames))
{
  # if(str_detect(re_feature,'volume')|str_detect(re_feature,'area')|str_detect(re_feature,'thickness'))
  # {next}
  
  if(str_detect(tem_colnames[i_sub],'Right'))
  {
    
    tem_colnames[i_sub]<-paste0(tem_colnames[i_sub],'.volume');
    tem_colnames[i_sub]<-gsub('Right','rh',tem_colnames[i_sub]);
    tem_colnames[i_sub]<-tolower(tem_colnames[i_sub])
  }
  
  if(str_detect(tem_colnames[i_sub],'Left'))
  {
    
    tem_colnames[i_sub]<-paste0(tem_colnames[i_sub],'.volume');
    tem_colnames[i_sub]<-gsub('Left','lh',tem_colnames[i_sub]);
    tem_colnames[i_sub]<-tolower(tem_colnames[i_sub])
  }
  
}

colnames(data)[start_col:dim(data)[2]]<-tem_colnames



datazzz<-data;

data<-datazzz
data<-na.omit(data)
data$Sex<-as.character(data$Sex)
data[data$Sex=='Male','Sex']<-c(1);
data[data$Sex=='Female','Sex']<-c(0);
data$Sex<-as.numeric(data$Sex)


# data[data$Diagnosis=='SVD','Diagnosis']<-'CSVD';
# data[data$Diagnosis=='AQP4Pos_NMOSD','Diagnosis']<-'NMOSD';
# 
# disease<-c('MCI','AD',"PD","CSVD",
#            "MS","NMOSD")


 
 disease<-c('MCI','AD') # for example data


data1<-data[data$Diagnosis=='HC',];

for(i_dis in disease)
{
  data1<-rbind(data1,data[data$Diagnosis==i_dis,])
  
}

data<-data1;  




if (!(dir.exists(paste0(savepath,'/',str))))
{dir.create(paste0(savepath,'/',str))}

setwd(paste0(savepath,'/',str))


# data[data$Sex=='Female','Sex']<-1;
# data[data$Sex=='Male','Sex']<-0;
# data$Sex<-as.numeric(data$Sex);

tem_data1<-data;
list_eta<-NULL
num=0;
list_roc<-list()
list_roc_results<-list()
list_stats<-data.frame(matrix(NA,5,length(disease)));
colnames(list_stats)<-disease;
rownames(list_stats)<-c('median','Q25','Q75','Min','Max');
model<-list()
model_sum<-list()
feature_sum<-list()
Disease_risk<-list()
HC_risk<-list()
sum1<-list()
for(ii in disease[1:length(disease)])
{
  ii
  
  num=num+1;

# tem_data<-data[data$Diagnosis=='HC'|data$Diagnosis==disease[num]&
#                  data$Age>=min(data[data$Diagnosis==disease[num],'Age'],na.rm=TRUE)&
#                  data$Age<=max(data[data$Diagnosis==disease[num],'Age'],na.rm=TRUE),];
  
  tem_data<-data[data$Age>=min(data[data$Diagnosis==disease[num],'Age'],na.rm=TRUE)&
                 data$Age<=max(data[data$Diagnosis==disease[num],'Age'],na.rm=TRUE),];
  tem_data_HC1<-tem_data[tem_data$Diagnosis=='HC',]
  tem_data_cur_disease<-tem_data[tem_data$Diagnosis==disease[num],];
  
  tem_data_HC<-NULL;tem_data_HC<-data.frame(tem_data_HC);
  
  for(i_sub_HC in 1:dim(tem_data_HC1)[1])
  {
  if (tem_data_HC1$Site_ZZZ[i_sub_HC] %in% unique(tem_data_cur_disease$Site_ZZZ))
  {
    tem_data_HC<-rbind(tem_data_HC,tem_data_HC1[i_sub_HC,]);
  }
  }
  
  
  if(dim(tem_data_HC)[1]>dim(tem_data_cur_disease)[1])
  {
  sampled_indices <- sample(1:nrow(tem_data_HC), size = dim(tem_data_cur_disease)[1], replace = FALSE)
  tem_data_HC0<-tem_data_HC[sampled_indices, ]
  }
  
  if(dim(tem_data_HC)[1]<=dim(tem_data_cur_disease)[1])
  {
    tem_data_HC0<-tem_data_HC
  }
  
  tem_data_disease0<-tem_data[tem_data$Diagnosis!='HC',]
  
  tem_data<-rbind(tem_data_HC0,tem_data_disease0);
  tem_data[tem_data$Diagnosis!=disease[num],'Diagnosis']<-c('HC')
  
  
  tem_data<-na.omit(tem_data)
  tem_data$Diagnosis<-factor(tem_data$Diagnosis,levels=c('HC',disease[num]))
  library(MatchIt)
  #tem_data$Sex<-as.factor(tem_data$Sex)
  #match.it <- matchit(Diagnosis ~ Age + Sex+ Site_ZZZ, data = tem_data, method="nearest", ratio=1)
  match.it <- matchit(Diagnosis ~ Age + Sex, data = tem_data, method="nearest", ratio=1)
  #a <- summary(match.it)
  tem_data<-tem_data[c(rownames(match.it$match.matrix),match.it$match.matrix),]
  tem_data<-na.omit(tem_data)
  

tem_data$Diagnosis<-factor(tem_data$Diagnosis,
                           levels = c('HC',disease[num]));

tem_data0<-tem_data[,c('Diagnosis','Age','Sex',colnames(tem_data)[start_col:dim(tem_data)[2]])];

#tem_data0<-tem_data0[!is.na(tem_data0$Age)&!is.na(tem_data0$Sex)&!is.na(tem_data0[,4]),]

tem_data0<-na.omit(tem_data0);

tem_data0[tem_data0$Diagnosis=='HC','label']<-0;

tem_data0[tem_data0$Diagnosis==disease[num],'label']<-1;


tem_data0<-tem_data0[is.finite(rowSums(tem_data0[,4:dim(tem_data0)[2]])),]

tem_data0<-tem_data0[,2:dim(tem_data0)[2]]

tem_data0<-tem_data0[,sapply(tem_data0,is.numeric)]

#SVM classificaiton
library(caret)
library(e1071)
library(kernlab)

predict_value<-NULL
predict_value<-data.frame(predict_value)
predict_value_new<-NULL
predict_value_new<-data.frame(predict_value_new)

set.seed(123);

list_save<-list()
model_sel<-'svmRadial'
#spit data into training and testing ten-folds cross validation
k=10;kin<-5
len=nrow(tem_data0);
pool<-c(1:len);

for(i in c(1:k))
{
  seg=as.integer(len/k);
  print(i)
  if(i<k)
  {
    indexk<-sample(x=pool,size = seg,replace = F)
    pool<-setdiff(pool,indexk)
    patient_ID<-rownames(tem_data0)[indexk]
    
    testing<-tem_data0[indexk,]
    training<-tem_data0[-indexk,]
  }
  
  if(i==k)
  {
    indexk<-pool
    patient_ID<-rownames(tem_data0)[indexk]
    
    testing<-tem_data0[indexk,]
    training<-tem_data0[-indexk,]
  }
  
  y=training[,'label'];
  x=training[,c(1:(dim(testing)[2]-1))]
  
  
  #lasso regression for feature reduction for each fold
  cvfit<-cv.glmnet(as.matrix(x),as.numeric(y),family = "binomial", 
                   nfolds = 10,alpha = 1,nlambda=50,
                   gamma = c(0, 0.25, 0.5, 0.75, 1),relax=TRUE,
                   type.measure='auc',parallel=TRUE);
  #plot(cvfit)
  cvfit$lambda.min
  cvfit$lambda.1se
  
  l.coef2<-coef(cvfit$glmnet.fit,s=cvfit$lambda.min,exact=F)
  tem_feature<-l.coef2@Dimnames[[1]]
  sel_feature<-tem_feature[l.coef2@i+1];
  sel_feature<-sel_feature[c(2:length(sel_feature))]
  
  training<-training[,c('label',sel_feature)]
  testing<-testing[,c('label',sel_feature)]
 
  training$label<-factor(training$label)
  testing$label<-factor(testing$label)
  
  

  tuned_model <- e1071::tune.svm(label ~ .,
                          data = training,
                          kernel='radial', 
                          type='C-classification',#'eps-regression'
                          #degree=c(1,2,3),
                          gamma = c(0.01, 0.1, 1),
                          cost=c(0.01, 0.1, 1, 10, 100),
                          cross = kin)
  
  tuned_model$best.parameters
  
  
  svm_model <- svm(
    label ~ .,
    data = training,
    kernel='radial',
    type='C-classification',
    cost = tuned_model$best.parameters$cost,
    gamma = tuned_model$best.parameters$gamma,
    probability=TRUE
  )
  
  

  
  #use the model for testing data
  prediction <- predict(svm_model, newdata = testing,probability = TRUE)
  
  probabilities <- data.frame(attr(prediction, "probabilities"))
  
  prediction<-probabilities$X1
  

  prediction<-data.frame(prediction)
  rownames(prediction)<-rownames(testing)
  predict_value<-rbind(predict_value,prediction)
  print(predict_value)
  #use the model for new data
  data_new0<-data
  data_new0[is.na(data_new0)]<-0.5
  
  prediction<-predict(svm_model,newdata=data_new0,probability = TRUE);
  
  probabilities <- data.frame(attr(prediction, "probabilities"))
  
  prediction<-probabilities$X1
  
  prediction<-data.frame(prediction)
  rownames(prediction)<-rownames(data_new0)
  if(i==1){predict_value_new<-prediction}
  if(i>1){predict_value_new<-cbind(predict_value_new,prediction)}
  
  
  list_save[[paste0(model_sel,'_',ii,'_Select_feature_fold_',i)]]<-sel_feature
  list_save[[paste0(model_sel,'_',ii,'_Model_fold_',i)]]<-svm_model
  
}
  

saveRDS(list_save,paste0(ii,'_model_cross_validation.rds'))

ind<-intersect(rownames(tem_data0),rownames(predict_value))


list_save[[paste0(model_sel,'_feature_in_all_folds')]]<-
  intersect(list_save[[paste0(model_sel,'_',ii,'_Select_feature_fold_',1)]],
            list_save[[paste0(model_sel,'_',ii,'_Select_feature_fold_',2)]])

for(i in c(3:k))
{
  list_save[[paste0(model_sel,'_feature_in_all_folds')]]<-
    intersect(list_save[[paste0(model_sel,'_',ii,'_Select_feature_fold_',i)]],
              list_save[[paste0(model_sel,'_feature_in_all_folds')]])
}


tem_data0[ind,'predict_label']<-predict_value[ind,1]

  

ind<-intersect(rownames(data_new0),rownames(predict_value_new))

for(col in 1:dim(predict_value_new)[2])
{
colnames(predict_value_new)[col]<-c(paste0('predict_label',col))
}

predict_value_new1<-data.frame((predict_value_new))
predict_value_new1[is.na(predict_value_new1)]<-0
#predict_value_new1<-as.numeric(predict_value_new1)
for(i in 1:dim(predict_value_new1[1]))
{
  tem_predict_value_new1<-predict_value_new1[i,c(1:k)]
  predict_value_new1[i,'predict_label']<-mean(as.numeric(tem_predict_value_new1))
  
}

data_new0[ind,disease[num]]<-predict_value_new1[ind,'predict_label']


real_label<-tem_data0$label;
predict_label<-as.numeric(tem_data0$predict_label)

library(pROC)
list_roc[[disease[num]]]<-pROC::roc(real_label,predict_label,ci=TRUE,boot.n = 1000)

tem_roc<-list_roc[[disease[num]]]
Sen<-tem_roc$sensitivities[which.max(tem_roc$sensitivities+tem_roc$specificities-1)]
Spe<-tem_roc$specificities[which.max(tem_roc$sensitivities+tem_roc$specificities-1)]
AUC<-tem_roc$auc
AUC_low_CI<-tem_roc$ci[1]
AUC_high_CI<-tem_roc$ci[3]

thr<-tem_roc$thresholds[which.max(tem_roc$sensitivities+tem_roc$specificities-1)]

predict_value=ifelse(tem_roc$original.predictor>thr,1,0);

ACC=sum(predict_value==tem_roc$original.response)/length(predict_value)

conf_matrix <- confusionMatrix(as.factor(predict_value), as.factor(tem_roc$original.response))

list_roc_results[[disease[num]]]['Sen']<-Sen
list_roc_results[[disease[num]]]['Spe']<-Spe
list_roc_results[[disease[num]]]['Acc']<-ACC
list_roc_results[[disease[num]]]['AUC']<-AUC
list_roc_results[[disease[num]]]['AUC_low_CI']<-AUC_low_CI
list_roc_results[[disease[num]]]['AUC_high_CI']<-AUC_high_CI
list_roc_results[[disease[num]]]['features']<-list_save
list_roc_results[[disease[num]]]['conf_matrix']<-conf_matrix 


tem_data1<-data_new0

tem_data2<-tem_data1[tem_data1$Diagnosis=='HC',];

for(dis_z in c(disease))
{
  tem_data2<-rbind(tem_data2,tem_data1[tem_data1$Diagnosis==dis_z,])
}

HC_risk[[disease[num]]]<-tem_data2

tem_data2$Diagnosis<-factor(tem_data2$Diagnosis,levels=c('HC',disease))

rad_data<-aggregate(tem_data2[,disease[num]],by=list(type=tem_data2$Diagnosis),mean,na.rm=TRUE)

list_stats[1,disease[num]]<-quantile(tem_data2[tem_data2$Diagnosis==disease[num],disease[num]],0.5,na.rm=TRUE);
list_stats[2,disease[num]]<-quantile(tem_data2[tem_data2$Diagnosis==disease[num],disease[num]],0.25,na.rm=TRUE);
list_stats[3,disease[num]]<-quantile(tem_data2[tem_data2$Diagnosis==disease[num],disease[num]],0.75,na.rm=TRUE);
list_stats[4,disease[num]]<-quantile(tem_data2[tem_data2$Diagnosis==disease[num],disease[num]],0,na.rm=TRUE);
list_stats[5,disease[num]]<-quantile(tem_data2[tem_data2$Diagnosis==disease[num],disease[num]],1,na.rm=TRUE);


rad_data<-data.frame(rad_data)
rad_data1<-rad_data


boot_num=1000;
mean_value_boot<-data.frame(matrix(NA,1000,1))
for(dis in disease)
{

  tem_dis<-tem_data2[tem_data2$Diagnosis==dis,disease[num]];
  #mean_value[dis]<-mean(tem_dis,na.rm=TRUE);
  #sd_value[dis]<-sd(tem_dis,na.rm=TRUE);

  for(j in c(1:boot_num))
  {
    tem_dis1<-tem_dis[sample(1:length(tem_dis), length(tem_dis), replace = TRUE)]
    mean_value_boot[j,1]<-median(tem_dis1,na.rm=TRUE);
    #sd_value_boot[j]<-sd(tem_dis1,na.rm=TRUE);
  }

  sum1[ii]<-print(paste0(ii,'-',dis,':mean=',median(tem_dis,na.rm=TRUE),
                         '  quantile 25%=',quantile(tem_dis,0.25,na.rm=TRUE),
                         '  quantile 75%=',quantile(tem_dis,0.75,na.rm=TRUE),
               '  lower CI=',quantile(mean_value_boot[1,],0.025),
               '  higher CI=',quantile(mean_value_boot[1,],0.975)))
  # print(mean(tem_dis,na.rm=TRUE))
  # print(quantile(mean_value_boot[1,],0.025))
  # print(quantile(mean_value_boot[1,],0.975))
}


library(fmsb)
rownames(rad_data1)<-rad_data1$type
rad_data1[,'Max']<-(1)
rad_data1[,'Min']<-(0)
rad_data1[,'Value']<-rad_data1$x

print(rad_data1)
rad_data0 <- data.frame(t(rad_data1[,c("Max", "Min", "Value")]));

Disease_risk[[disease[num]]]<-rad_data1

png(filename = paste0(str_lab,'_',disease[num],'_model_plot.png'),
    width = 1480,
    height = 1480,
    units = "px",
    bg = "white",
    res = 300)
colnames(rad_data0)[dim(rad_data0)[2]]<-c('NMO')
p<-radarchart(rad_data0,
              axistype = 1,
              # Customize the polygon
              pcol = "grey", 
              pfcol = scales::alpha("grey", 0.5), 
              plwd = 1, plty = 1,
              # Customize the grid
              cglcol = "grey", cglty = 2, 
              cglwd = 0.8,
              # Customize the axis
              axislabcol = "grey",
              # Variable labels
              vlcex = 1,palcex = 1, 
              vlabels = colnames(rad_data0),
              caxislabels = c(0,rad_data0['Max',1]/4, 
                              rad_data0['Max',1]/4*2, 
                              rad_data0['Max',1]/4*3, 
                              rad_data0['Max',1]/4*4))
print(p)
dev.off()




}



png(filename = paste0(str,'_ROC_all_diseases_lasso_selection',".png"),
    width = 1480,
    height = 1000,
    units = "px",
    bg = "white",
    res = 300)

p<-ggroc(list_roc,linetype='solid',linewidth=1,alpha=0.5)+
  #geom_segment(aes(x=1,xend=0,y=1,yend=0),color='darkgrey',linetype='dashed')+
  ggsci::scale_color_npg()+
  guides(color = guide_legend(title = "Group"))+
  theme_classic()
print(p)

dev.off()


res<-list();
res$list_roc_results_lasso_selection<-list_roc_results
res$list_roc_lasso_selection<-list_roc
res$model<-model
#res$feature<-feature_sum
res$disease_risk_prediction<-Disease_risk
res$HC_risk_prediction<-HC_risk
res$list_stats<-list_stats
#res$list_save<-list_save
#res$predicted_all_diseases<-tem_data2

print(list_stats)
saveRDS(res,paste0(str,'_ROC_loop_lasso_selection.rds'));


for(dis in disease)
{
  tem_HC<-HC_risk[[dis]]
  tem_HC<-tem_HC[tem_HC$Diagnosis=='HC',];
  
  print(paste0(dis,dis,dis,dis))
  print(sum(tem_HC[,dis]>=0.95))
  print(sum(tem_HC[,dis]>=0.95)/dim(tem_HC)[1]);
  print(median(tem_HC[tem_HC[,dis]>=0.95,'Age'],na.rm=T))
  print(quantile(tem_HC[tem_HC[,dis]>=0.95,'Age'],0.25,na.rm=T))
  print(quantile(tem_HC[tem_HC[,dis]>=0.95,'Age'],0.75,na.rm=T))
  #print(sum(tem_HC[tem_HC[,dis]>=0.95,'Sex']))
  print((length(tem_HC[tem_HC[,dis]>=0.95,'Sex'])-sum(tem_HC[tem_HC[,dis]>=0.95,'Sex'])))
  print((length(tem_HC[tem_HC[,dis]>=0.95,'Sex'])-sum(tem_HC[tem_HC[,dis]>=0.95,'Sex']))/
          length(tem_HC[tem_HC[,dis]>=0.95,'Sex']))
  print(paste0("  "))
}


print(list_stats)
