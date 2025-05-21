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


datapath='C:/ZZZ/work/manuscript/Lifespan-main/'# change to your path
feature_path0='D:/Lifespan_freesurfer_results/Results/Step1_normative_models/V7_DK_no_CI/' #change to your path
savepath='D:/Lifespan_freesurfer_results/Results/Step5_prognosis/V7_DK/multivariate_analysis'#change to your path

#datapath='/mnt/sda/zzz/Lifespan-main/'# change
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

path<-list()

str_lab<-'global_and_local_feature'# change


str_feature<-c('centile');
# 
# feature_path0='D:/Lifespan_freesurfer_results/Results/Step1_normative_models/V7_DK_no_CI/'
# savepath='D:/Lifespan_freesurfer_results/Results/Step5_prognosis/V7_DK_zscore/multivariate_analysis'
# str_feature<-c('zscore');

# feature_path0='D:/Lifespan_freesurfer_results/Results/Step1_normative_models/V7_DK_no_CI/'
# savepath='D:/Lifespan_freesurfer_results/Results/Step5_prognosis/V7_DK_raw/multivariate_analysis'
# str_feature<-c('raw');


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
#data<-results$all_data[,1:17]
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

feature_col<-tem_colnames;

data_backup<-data;

data<-data_backup

add_data<-read.csv(paste0(datapath,'new_final_list1_update_20250325_update.csv'),header=TRUE);

rownames(add_data)<-paste0(add_data$Freesufer_Path2,add_data$Freesufer_Path3)

row_ind<-intersect(rownames(add_data),rownames(data))

data<-cbind(data[row_ind,],add_data[row_ind,])

data$Sex<-as.character(data$Sex)
data[data$Sex=='Male','Sex']<-c(1);
data[data$Sex=='Female','Sex']<-c(0);
data$Sex<-as.numeric(data$Sex)

data[data$Diagnosis=='SVD','Diagnosis']<-'CSVD';
data[data$Diagnosis=='AQP4Pos_NMOSD','Diagnosis']<-'NMOSD';

disease<-c("PD")

data1<-data[data$Diagnosis=='HC',];

for(i_dis in disease)
{
  data1<-rbind(data1,data[data$Diagnosis==i_dis,])
  
}

data<-data1;  


if (!(dir.exists(paste0(savepath,'/',str))))
{dir.create(paste0(savepath,'/',str))}

setwd(paste0(savepath,'/',str))


tem_data1<-data;
list_eta<-NULL
num=0;
list_roc<-list()
list_roc_results<-list()
model<-list()

setwd(savepath);



data$DBS_improvement_rate_off<-data$DBS_Improvement_difference_off/data$Preoperative_UPDRS_III_med_off
data[data$DBS_improvement_rate_off<0&!is.na(data$DBS_improvement_rate_off),'DBS_improvement_rate_off']<-NA
var_pre<-c("Preoperative_UPDRS_III_med_on" ,"Preoperative_med_improvement_rate" ,
           "Post_UPDRS_III_med_off","DBS_improvement_rate_off")

#feature_col<-c(feature_col,'Preoperative_UPDRS_III_med_off') #add baseline UPDRS as additional variable

data$Preoperative_UPDRS_III_med_on<-as.numeric(data$Preoperative_UPDRS_III_med_on)
data$Preoperative_med_improvement_rate<-as.numeric(data$Preoperative_med_improvement_rate)
data$Post_UPDRS_III_med_off<-as.numeric(data$Post_UPDRS_III_med_off)
data$DBS_improvement_rate_off<-as.numeric(data$DBS_improvement_rate_off)


results=list();
#model_sel_all<-c('svmRadial','glm','svmLinear',"svmPoly")
model_sel_all<-c('svmRadial')
#results=list();
list<-data.frame(matrix(0,length(var_pre)*length(c('HC',disease))*length(feature_col),25))

list_all<-data.frame(matrix(0,length(var_pre),15))

list_all_disease<-data.frame(matrix(0,length(var_pre)*length(c('HC',disease))*length(model_sel_all),19))

per_num=1000

list_clinical_measure_per<-data.frame(matrix(0,length(var_pre)*length(c('HC',disease))*per_num*length(model_sel_all),16))


colnames(list)<-c('clinical_feature','MRI_feature','disease','coefficient',
                  'pvalue_t','tvalue','degree1','degere2','Fvalue','R2','Adjust_R2',
                  'Eta2','Eta2_CI_low','Eta2_CI_high',
                  'Eta2_parital','Eta2_parital_CI_low','Eta2_parital_CI_high',
                  'rvalue','pvalue_r','MAE','MAPE','MSE','RMSE','MSLE','MedAE')

colnames(list_all)<-c('clinical_feature','pvalue_F',
                      'degree1','degree2','Fvalue','R2','Adjust_R2',
                      'rvalue','pvalue_r','MAE','MAPE','MSE','RMSE','MSLE','MedAE')

colnames(list_all_disease)<-c('clinical_feature','disease','pvalue_F',
                              'degree1','degree2','Fvalue','R2','Adjust_R2',
                              'rvalue','pvalue_r','rvalue_5%','rvalue_95%','MAE','MAPE','MSE','RMSE','MSLE','MedAE','Model')


colnames(list_clinical_measure_per)<-c('clinical_feature','disease','R2','Adjust_R2',
                                       'rvalue','pvalue_r','rvalue_5%','rvalue_95%','MAE','MAPE','MSE','RMSE','MSLE','MedAE','Per_index','Model')


num=0;num1=0;num2=0;num3=0
for(i_var in 1:length(var_pre))
{
  num1=num1+1;
  tem_data<-data[!is.na(data[,var_pre[i_var]]),];
  
  tem_data[tem_data$Diagnosis=='HC','G']<-'HC'
  tem_data[tem_data$Diagnosis!='HC','G']<-'Disease'
  
  tem_data$G<-factor(tem_data$G,levels=c('HC','Disease'))
  
  
  tem_data$Diagnosis<-factor(tem_data$Diagnosis,
                             levels = c('HC',disease));
  library(ppcor)
  #disease feature associations
  
  list_r<-data.frame(matrix(NA,length(c('HC',disease)),length(feature_col)))
  list_p<-data.frame(matrix(NA,length(c('HC',disease)),length(feature_col)))
  colnames(list_r)<-feature_col;
  colnames(list_p)<-feature_col;
  rownames(list_r)<-c('HC',disease)
  rownames(list_p)<-c('HC',disease)
  
  print('error1')
  for(i_dis in c('HC',disease))
  {
    
    for(i_feature in feature_col)
    {
      num=num+1;
      
      data_ass<-data[data$Diagnosis==i_dis,c(var_pre[i_var],i_feature)];
      data_ass<-na.omit(data_ass);if(dim(data_ass)[1]==0){next}
      tem_st<-cor.test(x=data_ass[,var_pre[i_var]],y=data_ass[,i_feature],method='pearson')
      
      data_ass_lm<-data_ass;colnames(data_ass_lm)[1:2]<-c('cli_feature','MR_feature');
      tem_lm<-lm(cli_feature~MR_feature,data=data_ass_lm);
      
      tem_eta2<-eta_squared(tem_lm,partial = FALSE);
      tem_eta2_par<-eta_squared(tem_lm);
      
      
      tem_sum<-summary(tem_lm)
      list[num,1]<-var_pre[i_var];
      list[num,2]<-i_feature;
      list[num,3]<-i_dis
      list[num,4]<-tem_sum$coefficients[2,1]
      list[num,5]<-tem_sum$coefficients[2,4]
      list[num,6]<-tem_sum$coefficients[2,3]
      list[num,7]<-tem_sum$df[1]-1
      list[num,8]<-tem_sum$df[2]
      list[num,9]<-tem_sum$fstatistic[1]
      list[num,10]<-tem_sum$r.squared
      list[num,11]<-tem_sum$adj.r.squared
      
      list[num,12]<-tem_eta2$Eta2[1];
      list[num,13]<-tem_eta2$CI_low[1];
      list[num,14]<-tem_eta2$CI_high[1];
      
      list[num,15]<-tem_eta2_par$Eta2[1] ;
      list[num,16]<-tem_eta2_par$CI_low[1];
      list[num,17]<-tem_eta2_par$CI_high[1];
      
      list[num,18]<-tem_st$estimate;
      list[num,19]<-tem_st$p.value;
      
      list[num,20]<-mean(abs(tem_lm$residuals))
      list[num,21]<-mean(abs(tem_lm$residuals)/abs(data_ass_lm$cli_feature))
      list[num,22]<-mean(abs(tem_lm$residuals**2))
      list[num,23]<-sqrt(mean(abs(tem_lm$residuals**2)))
      list[num,24]<-mean((log(1+data_ass_lm$cli_feature)-log(1+tem_lm$fitted.values))**2)
      list[num,25]<-median(abs(tem_lm$residuals))
      # list_rlist[num,3]<-tem_sum$adj.r.squared[num,i_feature]<-tem_st$estimate 
      # list_p[num,i_feature]<-tem_st$p.value
      
    }
  }
  
  
  results$signal_feature_lm<-list;
  
  print('error2')
  #for each disease
  for(i_dis in c('PD'))
  {
    
    if(i_dis!='SVD'&&var_pre[i_var]=="SVD_score"){next}
    
    #tem_data0<-tem_data[tem_data$Diagnosis==i_dis,c(var_pre[i_var],feature_col)];
    
    tem_data0<-tem_data[tem_data$Diagnosis==i_dis,c(feature_col,var_pre[i_var])];
    colnames(tem_data0)[dim(tem_data0)[2]]<-'label'
    #tem_data0[tem_data0$label<=0,'label']<-NA
    tem_data0<-na.omit(tem_data0);if(dim(tem_data0)[1]<10){next};
    
    #tem_data0$label<-tem_data0$label*100
    tem_data0_zz<-tem_data0;
    
    ###feature selection
    #SVR regression

    library(caret)
    library(e1071)
    library(kernlab)
    predict_value<-NULL
    
    set.seed(123);
    
    list_save<-list()
    
    
    for(model_sel in model_sel_all)
    {
    num2=num2+1;
    #model_sel<-'svmRadial'
    #model_sel<-'svmLinear'
    #model_sel<-"svmPoly"
    
    #model_sel<-'glm'
    #spit data into training and testing ten-folds cross validation
    k=10;kin<-5
    len=nrow(tem_data0);
    pool<-c(1:len);
    
    #lasso regression for feature reduction for each fold
    y=tem_data0[,'label'];
    x=tem_data0[,c(1:(dim(tem_data0)[2]-1))]
    cvfit<-cv.glmnet(as.matrix(x),as.numeric(y),family = "gaussian",
                     nfolds = 10,alpha = 1,nlambda=50,
                     gamma = c(0, 0.25, 0.5, 0.75, 1),relax=TRUE,
                     type.measure='mse',parallel=TRUE);

    plot(cvfit)
    cvfit$lambda.min
    cvfit$lambda.1se

    l.coef2<-coef(cvfit$glmnet.fit,s=1,exact=F)
    tem_feature<-l.coef2@Dimnames[[1]]
    sel_feature<-tem_feature[l.coef2@i+1];
    # 
    # 
    if(length(sel_feature)<=1)
    {sel_feature=colnames(x)}
    if(length(sel_feature)>1)
    {sel_feature<-sel_feature[c(2:length(sel_feature))]}

    # sel_feature=colnames(x)

    sel_feature<-list[list$disease==i_dis&
                      list$clinical_feature==var_pre[i_var]&
                        list$pvalue_t<=0.05,'MRI_feature']
    
    #SVR
    # library(doParallel)
    # cl <- makeCluster(12)
    # registerDoParallel(cl)
    
    ctrl <- trainControl(
      method = "cv",         
      number = kin,           
      #method = "LOOCV",
      savePredictions = TRUE,  
      verboseIter = TRUE      
      #allowParallel = TRUE
    )
    
    folds<-createFolds(y=tem_data0$label,k=k);
    
    predict_value<-NULL
    
    for(i in c(1:k))
    {
      
      print(i)
      
      testing<-tem_data0[folds[[i]],];
      training<-tem_data0[-folds[[i]],];
      # if(i<k)
      # {
      #   indexk<-sample(x=pool,size = seg,replace = F)
      #   pool<-setdiff(pool,indexk)
      #   patient_ID<-rownames(tem_data0)[indexk]
      # 
      #   testing<-tem_data0[indexk,]
      #   training<-tem_data0[-indexk,]
      #   #training<-tem_data0[sample(x=c(1:len),size = floor(len*0.9),
      #   #                          replace = F),]
      # }
      # 
      # if(i==k)
      # {
      #   indexk<-pool
      #   patient_ID<-rownames(tem_data0)[indexk]
      # 
      #   testing<-tem_data0[indexk,]
      #   training<-tem_data0[-indexk,]
      #   #training<-tem_data0[sample(x=c(1:len),size = floor(len*0.9),
      #   #                          replace = F),]
      # }
      # 
      #training<-tem_data0
      y=training[,'label'];
      x=training[,c(1:(dim(training)[2]-1))]
      
      
      #lasso regression for feature reduction for each fold
      cvfit<-cv.glmnet(as.matrix(x),as.numeric(y),family = "gaussian",
                       nfolds = 10,alpha = 1,nlambda=50,
                       gamma = c(0, 0.25, 0.5, 0.75, 1),relax=TRUE,
                       type.measure='auc',parallel=TRUE);

      #plot(cvfit)
      cvfit$lambda.min
      cvfit$lambda.1se

      l.coef2<-coef(cvfit$glmnet.fit,s=cvfit$lambda.min,exact=F)
      tem_feature<-l.coef2@Dimnames[[1]]
      sel_feature<-tem_feature[l.coef2@i+1];
      # 
      # 
      if(length(sel_feature)<=1)
      {sel_feature=colnames(x)}
      else
      {sel_feature<-sel_feature[c(2:length(sel_feature))]}

      #sel_feature=colnames(x)
      
      training<-training[,c('label',sel_feature)]
      testing<-testing[,c('label',sel_feature)]
      
if(model_sel=='svmLinear')
{
   
      tune_grid_linear <- expand.grid(C = 10^seq(-2, 2, length = 5))
      
      svm_model<- train(
        label ~ .,
        data = training,
        method = model_sel,    
        trControl = ctrl,
        tuneGrid = tune_grid_linear,
        family =  gaussian(link = "identity"),
        metric = "RMSE",
        preProcess = c("center", "scale")
      )
}
      
      
  if(model_sel=="svmPoly")
  {
    tune_grid_poly<- expand.grid(degree = c(2, 3),
                                 scale = c(0.1, 1),
                                 C = 2^seq(-5, 5, length = 5))
                                 
    svm_model= train(label ~ ., data = training, method = model_sel,
                 preProcess = c("center", "scale"),
                 tuneGrid = tune_grid_poly,
                 metric = "RMSE")
   }   
          
   if(model_sel=='svmRadial')
   {
      tune_grid_radial <- expand.grid(
        sigma = c(0.01, 0.1, 1),  
        C = 2^seq(-5, 15, by = 2)  
      )
      
      svm_model<- train(
        label ~ .,
        data = training,
        method = model_sel,    
        trControl = ctrl,
        tuneGrid = tune_grid_radial,
        family =  gaussian(link = "identity"),
        metric = "RMSE",
        preProcess = c("center", "scale")
      )
   }  
      
      
      if(model_sel=='glm')
      {
        
        svm_model<- glm(
          label ~ .,
          data = training,
          family =  gaussian(link = "identity")
        )
      }  
      # tuned_model <- e1071::tune.svm(label ~ .,
      #                                data = training,
      #                                kernel="radial",
      #                                type='eps-regression',
      #                                #degree=c(1,2,3),
      #                                gamma = 10^(-2:0),
      #                                cost=10^(-2:2),
      #                                epsilon = 10^(-2:0),
      #                                cross = 1)
      # 
      # tuned_model$best.parameters
      # 
      # 
      # svm_model <- svm(
      #   label ~ .,
      #   data = training,
      #   kernel="radial",
      #   type='eps-regression',
      #   cost = tuned_model$best.parameters$cost,
      #   gamma = tuned_model$best.parameters$gamma,
      #   epsilon=tuned_model$best.parameters$epsilon
      #   #epsilon=0.1
      # )
      
    
      
      #use the model for testing data
      prediction <- predict(svm_model, newdata = testing)
      
      prediction<-data.frame(prediction)
      rownames(prediction)<-rownames(testing)
      predict_value<-rbind(predict_value,prediction)
      
      
      #list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Select_feature_fold_',i)]]<-sel_feature
      list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Model_fold_',i)]]<-svm_model
      
    }
    
    #list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_feature')]]<-cvfit
    
    ind<-intersect(rownames(tem_data0),rownames(predict_value))
    
    
    # list_save[[paste0(model_sel,'_feature_in_all_folds')]]<-
    #   intersect(list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Select_feature_fold_',1)]],
    #             list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Select_feature_fold_',2)]])
    # 
    # for(i in c(3:k))
    # {
    #   list_save[[paste0(model_sel,'_feature_in_all_folds')]]<-
    #     intersect(list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Select_feature_fold_',i)]],
    #               list_save[[paste0(model_sel,'_feature_in_all_folds')]])
    # }
    # 
    # list_save[[paste0(model_sel,'_feature_in_all_folds')]]<-list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Select_feature_fold_',1)]]
    # 
    
    tem_data0[ind,'predict_label']<-predict_value[ind,1]
    
    
    real_label<-tem_data0$label;
    predict_label<-as.numeric(tem_data0$predict_label)
    
    # if(length(real_label[is.infinite(real_label)|is.na(real_label)])>0)
    # {
    #   real_label[is.infinite(real_label)|is.na(real_label)]<-
    #     mean(tem_data1[!(is.infinite(real_label)|is.na(real_label))])
    # }
    # 
    # if(length(predict_label[is.infinite(predict_label)|is.na(predict_label)])>0)
    # {
    #   predict_label[is.infinite(predict_label)|is.na(predict_label)]<-
    #     mean(tem_data1[!(is.infinite(predict_label)|is.na(predict_label))])
    # }
    # 
    
    # library(raincloudplots)
    # library(ggwordcloud)
    # setwd(savepath)
    # png(filename = paste0(i_dis,'_',var_pre[i_var],'_predictive_features_assocaition_plot.png'),
    #     width = 1480,
    #     height = 1480,
    #     units = "px",
    #     bg = "white",
    #     res = 300)
    # 
    # p<-ggplot(tem_list_fig, aes(label = MRI_feature,
    #                         size = -log(pvalue_r),color = rvalue)) +
    #   geom_text_wordcloud(eccentricity=1) +
    #   scale_size_area(max_size = 4)+
    #   theme_minimal()+
    #   scale_color_gradient2(limits=c(min(tem_list_fig$rvalue),max(tem_list_fig$rvalue)),low = "#00CCFF",mid="white",high = "#990000")
    # print(p)
    # dev.off()
    
    
    # setwd(savepath)
    # png(filename = paste0(i_dis,'_',var_pre[i_var],'_predictive_features_assocaition_plot.png'),
    #     width = 1480,
    #     height = 1480,
    #     units = "px",
    #     bg = "white",
    #     res = 300)
    # 
    # p<-ggplot(tem_data0, aes(x = label, y = predict_label)) +
    #   geom_point() +  
    #   geom_smooth(method = "lm", color = "darkred", fill = "gray", alpha = 0.3) +  # 添加线性拟合曲线和置信区间
    #   labs(x = "Actural measure", y = "Predicted measure") +
    #   theme_minimal()  
    # print(p)
    # dev.off()
    
    
    
    library(ggpubr)
    
    
    png(filename = paste0(i_dis,'_',var_pre[i_var],'_predictive_features_assocaition_plot_',model_sel,'.png'),
        width = 1000,
        height = 1000,
        units = "px",
        bg = "white",
        res = 300)
    
    p<-ggplot(tem_data0, aes(x = label, y = predict_label))+
      #geom_point(size=2)+
      geom_point(size=1,color='grey',shape=16)+
      geom_smooth(method = 'lm', formula = y ~ x, se = T,
                  linewidth=1,color='#990000')+
      stat_cor(data=tem_data0, method = "pearson")+
      labs(x=paste0('Real ',var_pre[i_var]),y=paste0('Predicted ',var_pre[i_var]))+
      theme_classic()
    print(p)
    dev.off()
    
    
    tem_cor<-cor.test(real_label,predict_label)
    
    ss_tot<-sum((real_label-mean(real_label))^2)
    ss_res<-sum((real_label-predict_label)^2)
    rsquared<-1-ss_res/ss_tot
    
  
    list_all_disease[num2,1]<-var_pre[i_var];
    list_all_disease[num2,2]<-i_dis;
    list_all_disease[num2,7]<-rsquared

    list_all_disease[num2,9]<-tem_cor$estimate;
    list_all_disease[num2,10]<-tem_cor$p.value;
    list_all_disease[num2,11]<-tem_cor$conf.int[1];
    list_all_disease[num2,12]<-tem_cor$conf.int[2];
    
    list_all_disease[num2,13]<-mean(abs(real_label-predict_label),na.rm = TRUE)
    list_all_disease[num2,14]<-mean(abs(real_label-predict_label)/abs(real_label),na.rm = TRUE)
    list_all_disease[num2,15]<-mean(abs((real_label-predict_label)**2),na.rm = TRUE)
    list_all_disease[num2,16]<-sqrt(mean(abs((real_label-predict_label)**2),na.rm = TRUE))
    list_all_disease[num2,17]<-mean((log(1+real_label)-log(1+predict_label))**2,na.rm = TRUE)
    list_all_disease[num2,18]<-median(abs((real_label-predict_label)),na.rm = TRUE)
    list_all_disease[num2,19]<-model_sel
    print(mean(abs(real_label-predict_label)/abs(real_label),na.rm = TRUE))
    
    }
  }
  
  
  
  
  print('error3')
  #for each disease for permutated data
  for(i_dis in c('PD'))
  {

    if(i_dis!='SVD'&&var_pre[i_var]=="SVD_score"){next}

    #tem_data0<-tem_data[tem_data$Diagnosis==i_dis,c(var_pre[i_var],feature_col)];

    tem_data00<-tem_data[tem_data$Diagnosis==i_dis,c(feature_col,var_pre[i_var])];
    colnames(tem_data00)[dim(tem_data00)[2]]<-'label'
    tem_data00<-na.omit(tem_data00);if(dim(tem_data00)[1]<10){next};
    tem_data00[tem_data00$label<=0,'label']<-mean(tem_data00[tem_data00$label>0,'label'])

    for(model_sel in model_sel_all)
    {
      

    for(per in 1:per_num)
    {
      num3=num3+1;
      tem_data0<-cbind(tem_data00[sample(nrow(tem_data00)),c(1:(dim(tem_data00)[2]-1))],tem_data00[,'label'])
      colnames(tem_data0)[dim(tem_data0)[2]]<-'label'


      ###feature selection

      #SVR 回归
      library(caret)
      library(e1071)
      library(kernlab)

      predict_value<-NULL
      predict_value_new<-NULL

      set.seed(123);

      list_save<-list()
      #model_sel<-'svmRadial'
      #spit data into training and testing ten-folds cross validation
      k=10;kin<-5
      len=nrow(tem_data0);
      pool<-c(1:len);
      
      folds<-createFolds(y=tem_data0$label,k=k);
      
      predict_value<-NULL
      

      for(i in c(1:k))
      {
        print(i)
        
        testing<-tem_data0[folds[[i]],];
        training<-tem_data0[-folds[[i]],];
     
        y=training[,'label'];
        x=training[,c(1:(dim(training)[2]-1))]
        
  
        
        #lasso regression for feature reduction for each fold
        cvfit<-cv.glmnet(as.matrix(x),as.numeric(y),family = "gaussian",
                         nfolds = 10,alpha = 1,nlambda=50,
                         gamma = c(0, 0.25, 0.5, 0.75, 1),relax=TRUE,
                         type.measure='mse',parallel=TRUE);

        #plot(cvfit)
        cvfit$lambda.min
        cvfit$lambda.1se

        l.coef2<-coef(cvfit$glmnet.fit,s=cvfit$lambda.min,exact=F)
        tem_feature<-l.coef2@Dimnames[[1]]
        sel_feature<-tem_feature[l.coef2@i+1];


        if(length(sel_feature)<=1)
        {sel_feature=colnames(x)}
        else
        {sel_feature<-sel_feature[c(2:length(sel_feature))]}
        
        #sel_feature=colnames(x)
        
        training<-training[,c('label',sel_feature)]
        testing<-testing[,c('label',sel_feature)]
        
        if(model_sel=='svmLinear')
        {
          
          tune_grid_linear <- expand.grid(C = 10^seq(-2, 2, length = 5))
          
          svm_model<- train(
            label ~ .,
            data = training,
            method = model_sel,    
            trControl = ctrl,
            tuneGrid = tune_grid_linear,
            family =  gaussian(link = "identity"),
            metric = "RMSE",
            preProcess = c("center", "scale")
          )
        }
        
        
        if(model_sel=="svmPoly")
        {
          tune_grid_poly<- expand.grid(degree = c(2, 3),
                                       scale = c(0.1, 1),
                                       C = 2^seq(-5, 5, length = 5))
          
          svm_model= train(label ~ ., data = training, method = model_sel,
                           preProcess = c("center", "scale"),
                           tuneGrid = tune_grid_poly,
                           metric = "RMSE")
        }   
        
        if(model_sel=='svmRadial')
        {
          tune_grid_radial <- expand.grid(
            sigma = c(0.01, 0.1, 1),  
            C = 2^seq(-5, 15, by = 2)  
          )
          
          svm_model<- train(
            label ~ .,
            data = training,
            method = model_sel,    
            trControl = ctrl,
            tuneGrid = tune_grid_radial,
            family =  gaussian(link = "identity"),
            metric = "RMSE",
            preProcess = c("center", "scale")
          )
        }  
        
        
        if(model_sel=='glm')
        {
          
          svm_model<- glm(
            label ~ .,
            data = training,
            family =  gaussian(link = "identity")
          )
        }  
  
        
        #use the model for testing data
        prediction <- predict(svm_model, newdata = testing)
        
        prediction<-data.frame(prediction)
        rownames(prediction)<-rownames(testing)
        predict_value<-rbind(predict_value,prediction)
        
        # list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Select_feature_fold_',i)]]<-sel_feature
        # list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Model_fold_',i)]]<-svm_model

      }

      ind<-intersect(rownames(tem_data0),rownames(predict_value))


      list_save[[paste0(model_sel,'_feature_in_all_folds')]]<-
        intersect(list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Select_feature_fold_',1)]],
                  list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Select_feature_fold_',2)]])

      for(i in c(3:k))
      {
        list_save[[paste0(model_sel,'_feature_in_all_folds')]]<-
          intersect(list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Select_feature_fold_',i)]],
                    list_save[[paste0(model_sel,'_feature_in_all_folds')]])
      }

      tem_data0[ind,'predict_label']<-predict_value[ind,1]


      real_label<-tem_data0$label;
      predict_label<-as.numeric(tem_data0$predict_label)



      if(length(real_label[is.infinite(real_label)|is.na(real_label)])>0)
      {
        real_label[is.infinite(real_label)|is.na(real_label)]<-
          mean(tem_data1[!(is.infinite(real_label)|is.na(real_label))])
      }

      if(length(predict_label[is.infinite(predict_label)|is.na(predict_label)])>0)
      {
        predict_label[is.infinite(predict_label)|is.na(predict_label)]<-
          mean(tem_data1[!(is.infinite(predict_label)|is.na(predict_label))])
      }


      tem_cor<-cor.test(real_label,predict_label)

      ss_tot<-sum((real_label-mean(real_label))^2)
      ss_res<-sum((real_label-predict_label)^2)
      rsquared<-1-ss_res/ss_tot

      #tem_sum<-tem_pre
      list_clinical_measure_per[num3,1]<-var_pre[i_var];
      list_clinical_measure_per[num3,2]<-i_dis;
      list_clinical_measure_per[num3,3]<-rsquared
      #list_clinical_measure_per[num3,4]<-tem_sum$adj.r.squared

      list_clinical_measure_per[num3,5]<-tem_cor$estimate;
      list_clinical_measure_per[num3,6]<-tem_cor$p.value;
      list_clinical_measure_per[num2,7]<-tem_cor$conf.int[1];
      list_clinical_measure_per[num2,8]<-tem_cor$conf.int[2];

      list_clinical_measure_per[num3,9]<-mean(abs(real_label-predict_label),na.rm = TRUE)
      list_clinical_measure_per[num3,10]<-mean(abs(real_label-predict_label)/abs(real_label),na.rm = TRUE)
      list_clinical_measure_per[num3,11]<-mean(abs((real_label-predict_label)**2),na.rm = TRUE)
      list_clinical_measure_per[num3,12]<-sqrt(mean(abs((real_label-predict_label)**2),na.rm = TRUE))
      list_clinical_measure_per[num3,13]<-mean((log(1+real_label)-log(1+predict_label))**2,na.rm = TRUE)
      list_clinical_measure_per[num3,14]<-median(abs((real_label-predict_label)),na.rm = TRUE)
      list_clinical_measure_per[num3,15]<-per
      list_clinical_measure_per[num3,16]<-model_sel
    }

  }
  }
  
}

setwd(savepath)
list[,'pFDR']<-p.adjust(list$pvalue_r,'fdr');
list_all_disease[,'pFDR']<-p.adjust(list_all_disease$pvalue_r,'fdr');
results$all_feature_svm_feature_disease_seperate<-list;
results$all_feature_svm_disease_seperate<-list_all_disease;
results$all_feature_svm_disease_seperate_per<-list_clinical_measure_per;

write.csv(list_all_disease,'clinical_association_all_diseases.csv')
saveRDS(results,'clinical_association.rds')


