#install.packages('gamlss')
rm(list=ls())


#set your own directory
datapath='***/Source-codes/'#Change the directory where you save the Source-codes 
feature_path0='**/Diseases'; #Chnage the path you have saved your normative models in "Disease-application-normative-model.R"
savepath='**/Results/Clinical_associaton';# Create and determine the directory where you would save the results  


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
library(caret)
library(e1071)
library(kernlab)
library(ppcor)


setwd(datapath)#set the filepath
source("100.common-variables.r")
source("101.common-functions.r")

source("ZZZ_function.R")
source("300.variables.r")
source("301.functions.r")


path<-list()

str_lab<-'global_and_local_feature_permutation'# change
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
    
    #feature_col<-c(feature_col,results$i)
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


# data[data$Diagnosis=='SVD','Diagnosis']<-'CSVD';
# data[data$Diagnosis=='AQP4Pos_NMOSD','Diagnosis']<-'NMOSD';
# 
# disease<-c('MCI','AD',"PD","CSVD",
#            "MS","NMOSD")

disease<-c('MCI','AD')

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




# data[,'Post_DBS_UPDRS_off']<-data$Post_UPDRS_III_med_off;
# data[,'Pre_UPDRS_III_med_off']<-data$Preoperative_UPDRS_III_med_off;
# 
# 
# var_pre<-c('MMSE','MOCA','BVMT','CVLT','SDMT','PASAT',
#            'Pre_UPDRS_III_med_off','EDSS')

var_pre<-c('MMSE','MOCA') # for example data


data$MOCA<-as.numeric(data$MOCA)
data$MMSE<-as.numeric(data$MMSE)
# data$BVMT<-as.numeric(data$BVMT)
# data$CVLT<-as.numeric(data$CVLT)
# data$SDMT<-as.numeric(data$SDMT)
# data$PASAT<-as.numeric(data$PASAT)
# data$Pre_UPDRS_III_med_off<-as.numeric(data$Pre_UPDRS_III_med_off)
# data$EDSS<-as.numeric(data$EDSS)

results=list();
model_sel_all<-c('svmRadial')
#results=list();
list<-data.frame(matrix(0,length(var_pre)*length(c('HC',disease))*length(feature_col),25))

list_all<-data.frame(matrix(0,length(var_pre),15))

list_all_disease<-data.frame(matrix(0,length(var_pre)*length(c('HC',disease))*length(model_sel_all),19))


per_num=1000

list_clinical_measure_per<-data.frame(matrix(0,length(var_pre)*length(c('HC',disease))*per_num,15))


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
                                       'rvalue','pvalue_r','rvalue_5%','rvalue_95%','MAE','MAPE','MSE','RMSE','MSLE','MedAE','Per_index')


num=0;num1=0;num2=0;num3=0
for(i_var in 1:length(var_pre))
{
  #num1=num1+1;
tem_data<-data[!is.na(data[,var_pre[i_var]]),];

tem_data[tem_data$Diagnosis=='HC','G']<-'HC'
tem_data[tem_data$Diagnosis!='HC','G']<-'Disease'

tem_data$G<-factor(tem_data$G,levels=c('HC','Disease'))


tem_data$Diagnosis<-factor(tem_data$Diagnosis,
                           levels = c('HC',disease));

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
    
    #data_ass <- data.frame(data_ass)
    #data_ass<-as.numeric(as.character(data_ass))
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


model_sel<-'svmRadial'
#model_sel<-'svmLinear'
#model_sel<-'glm'
#model_sel<-'svmPoly'
#for each disease
for(i_dis in c('HC',disease))
{
num2=num2+1;
if(i_dis!='SVD'&&var_pre[i_var]=="SVD_score"){next}
if(!(i_dis %in% c('HC','NMOSD','MS'))&&var_pre[i_var]=="BVMT"){next}
if(!(i_dis %in% c('HC','NMOSD','MS'))&&var_pre[i_var]=="CVLT"){next}
if(!(i_dis %in% c('HC','NMOSD','MS'))&&var_pre[i_var]=="SDMT"){next}
if(!(i_dis %in% c('NMOSD','MS'))&&var_pre[i_var]=="EDSS"){next}
if(!(i_dis %in% c('HC','NMOSD','MS'))&&var_pre[i_var]=="PASAT"){next}
if(!(i_dis %in% c('PD'))&&var_pre[i_var]=="Pre_UPDRS_III_med_off"){next}
#tem_data0<-tem_data[tem_data$Diagnosis==i_dis,c(var_pre[i_var],feature_col)];

tem_data0<-tem_data[tem_data$Diagnosis==i_dis,c(feature_col,var_pre[i_var])];
colnames(tem_data0)[dim(tem_data0)[2]]<-'label'
tem_data0<-na.omit(tem_data0);if(dim(tem_data0)[1]<10){next};
tem_data0[tem_data0$label<=0,'label']<-mean(tem_data0[tem_data0$label>0,'label'])
tem_data0_zz<-tem_data0;


#SVR
predict_value<-NULL
predict_value_new<-NULL

set.seed(123);

list_save<-list()

#spit data into training and testing ten-folds cross validation
k=10;kin<-5
#k=dim(tem_data0)[1];kin<-5
len=nrow(tem_data0);
pool<-c(1:len);
real_k=0;

# library(doParallel)
# cl <- makeCluster(12)
# registerDoParallel(cl)

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
)

folds<-createFolds(y=tem_data0$label,k=k);


for(i in c(1:k))
{
  
  print(i)
  
  testing<-tem_data0[folds[[i]],];
  training<-tem_data0[-folds[[i]],];

  y=training[,'label'];
  x=training[,c(1:(dim(testing)[2]-1))]
  
  
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
  
  
  #real_k=real_k+1
  
  if(length(sel_feature)<=1)
  {sel_feature=colnames(x)}
  else
  {sel_feature<-sel_feature[c(2:length(sel_feature))]}
  

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
      metric = "RMSE"
      #preProcess = c("center", "scale"),
    )
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
      metric = "RMSE"
      #preProcess = c("center", "scale"),
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
  
  
  if(model_sel=='glm')
  {
    
    svm_model<- glm(
      label ~ .,
      data = training,
      family =  gaussian(link = "identity")
    )
  }  
  
  #whichModel <- oneSE(svm_model$results, metric = "Accuracy", num=1000,maximize = TRUE)
  
  
  #use the model for testing data
  prediction <- predict(svm_model, newdata = testing)
  
  prediction<-data.frame(prediction)
  rownames(prediction)<-rownames(testing)
  predict_value<-rbind(predict_value,prediction)
  

  list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Select_feature_fold_',i)]]<-sel_feature
  list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Model_fold_',i)]]<-svm_model
  
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


setwd(savepath)
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

#tem_sum<-tem_pre
list_all_disease[num2,1]<-var_pre[i_var];
list_all_disease[num2,2]<-i_dis;
# list_all_disease[num2,3]<-pf(tem_sum$fstatistic[1], tem_sum$df[1]-1, tem_sum$df[2], lower.tail = FALSE)
# list_all_disease[num2,4]<-tem_sum$df[1]-1
# list_all_disease[num2,5]<-tem_sum$df[2]
# list_all_disease[num2,6]<-tem_sum$fstatistic[1]
list_all_disease[num2,7]<-rsquared
#list_all_disease[num2,8]<-tem_sum$adj.r.squared

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



# 
print('error3')
#for each disease for permutated data
for(i_dis in c('HC',disease))
{

  if(i_dis!='SVD'&&var_pre[i_var]=="SVD_score"){next}
  if(!(i_dis %in% c('HC','NMOSD','MS'))&&var_pre[i_var]=="BVMT"){next}
  if(!(i_dis %in% c('HC','NMOSD','MS'))&&var_pre[i_var]=="CVLT"){next}
  if(!(i_dis %in% c('HC','NMOSD','MS'))&&var_pre[i_var]=="SDMT"){next}
  if(!(i_dis %in% c('NMOSD','MS'))&&var_pre[i_var]=="EDSS"){next}
  if(!(i_dis %in% c('HC','NMOSD','MS'))&&var_pre[i_var]=="PASAT"){next}
  if(!(i_dis %in% c('PD'))&&var_pre[i_var]=="Pre_UPDRS_III_med_off"){next}

  #tem_data0<-tem_data[tem_data$Diagnosis==i_dis,c(var_pre[i_var],feature_col)];

  tem_data00<-tem_data[tem_data$Diagnosis==i_dis,c(feature_col,var_pre[i_var])];
  colnames(tem_data00)[dim(tem_data00)[2]]<-'label'
  tem_data00<-na.omit(tem_data00);if(dim(tem_data00)[1]<10){next};
  tem_data00[tem_data00$label<=0,'label']<-mean(tem_data00[tem_data00$label>0,'label'])

  for(per in 1:per_num)
  {
   num3=num3+1;
   tem_data0<-cbind(tem_data00[sample(nrow(tem_data00)),c(1:(dim(tem_data00)[2]-1))],tem_data00[,'label'])
   colnames(tem_data0)[dim(tem_data0)[2]]<-'label'

  #SVR regression

  predict_value<-NULL
  predict_value_new<-NULL

  set.seed(123);

  list_save<-list()
  model_sel<-'svmRadial'
  #spit data into training and testing ten-folds cross validation
  k=10;kin<-5
  len=nrow(tem_data0);
  pool<-c(1:len);

  
  # sel_feature<-list[list$disease==i_dis&
  #                     list$clinical_feature==var_pre[i_var]&
  #                     list$pvalue_t<=0.05,'MRI_feature']
  
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
  )
  
  folds<-createFolds(y=tem_data0$label,k=k);
  
  
  for(i in c(1:k))
  {

    print(i)
    
    testing<-tem_data0[folds[[i]],];
    training<-tem_data0[-folds[[i]],];
    
    y=training[,'label'];
    x=training[,c(1:(dim(testing)[2]-1))]
    
    
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
    
    
    #real_k=real_k+1
    
    if(length(sel_feature)<=1)
    {sel_feature=colnames(x)}
    else
    {sel_feature<-sel_feature[c(2:length(sel_feature))]}

    
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
        metric = "RMSE"
        #preProcess = c("center", "scale"),
      )
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
        metric = "RMSE"
        #preProcess = c("center", "scale"),
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
    
    
    if(model_sel=='glm')
    {
      
      svm_model<- glm(
        label ~ .,
        data = training,
        family =  gaussian(link = "identity")
      )
    }  
    
    #whichModel <- oneSE(svm_model$results, metric = "Accuracy", num=1000,maximize = TRUE)
    
    
    #use the model for testing data
    prediction <- predict(svm_model, newdata = testing)

    prediction<-data.frame(prediction)
    rownames(prediction)<-rownames(testing)
    predict_value<-rbind(predict_value,prediction)


    # list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Select_feature_fold_',i)]]<-sel_feature
    # list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Model_fold_',i)]]<-svm_model

  }

  ind<-intersect(rownames(tem_data0),rownames(predict_value))


  # list_save[[paste0(model_sel,'_feature_in_all_folds')]]<-
  #   intersect(list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Select_feature_fold_',1)]],
  #             list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Select_feature_fold_',2)]])

  # for(i in c(3:k))
  # {
  #   list_save[[paste0(model_sel,'_feature_in_all_folds')]]<-
  #     intersect(list_save[[paste0(model_sel,'_',var_pre[i_var],'_',i_dis,'_Select_feature_fold_',i)]],
  #               list_save[[paste0(model_sel,'_feature_in_all_folds')]])
  # }

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
  list_clinical_measure_per[num3,7]<-tem_cor$conf.int[1];
  list_clinical_measure_per[num3,8]<-tem_cor$conf.int[2];

  list_clinical_measure_per[num3,9]<-mean(abs(real_label-predict_label),na.rm = TRUE)
  list_clinical_measure_per[num3,10]<-mean(abs(real_label-predict_label)/abs(real_label),na.rm = TRUE)
  list_clinical_measure_per[num3,11]<-mean(abs((real_label-predict_label)**2),na.rm = TRUE)
  list_clinical_measure_per[num3,12]<-sqrt(mean(abs((real_label-predict_label)**2),na.rm = TRUE))
  list_clinical_measure_per[num3,13]<-mean((log(1+real_label)-log(1+predict_label))**2,na.rm = TRUE)
  list_clinical_measure_per[num3,14]<-median(abs((real_label-predict_label)),na.rm = TRUE)
  list_clinical_measure_per[num3,15]<-per


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

# #plot list_all_disease
# png(filename = paste0('heatmap_clinical_association_step_',model_sel,'.png'),
#     width = 2000,
#     height = 1480,
#     units = "px",
#     bg = "white",
#     res = 300)
# list_all_disease<-list_all_disease[!is.na(list_all_disease$disease),]
# list_all_disease$clinical_feature<-factor(list_all_disease$clinical_feature,levels=rev(c(var_pre)))
# list_all_disease$disease<-factor(list_all_disease$disease,c('HC',disease))
# heatmap=ggplot(data=list_all_disease,aes(x=disease,y=clinical_feature,fill=round(R2,2)))+
#   geom_tile()+
#   scale_fill_gradient2(limits=c(0,1),low='#4FBBD8',midpoint =0.5,high='#E84935',space='Lab')+
#   geom_text(aes(disease,clinical_feature,label=round(R2,2)),color='black',size=4)+
#   theme_minimal()+
#   theme(panel.grid.major = element_blank(),
#         axis.text.x = element_text(angle=45,vjust=1,size=10,hjust=1,color='black', family = "serif"),
#         axis.text.y = element_text(angle=0,vjust=0,size=10,hjust=1,color='black', family = "serif"))+
#   geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5),color='white',size=1)+
#   #geom_hline(yintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),color='white',size=1)
#   geom_hline(yintercept=c(1:12)-0.5,color='white',size=1)
# print(heatmap)
# dev.off()
# 
# png(filename = paste0('heatmap_clinical_association_adjust_R2_step_',model_sel,'.png'),
#     width = 3200,
#     height = 1880,
#     units = "px",
#     bg = "white",
#     res = 300)
# list_all_disease<-list_all_disease[!is.na(list_all_disease$disease),]
# list_all_disease$clinical_feature<-factor(list_all_disease$clinical_feature,levels=rev(c(var_pre)))
# list_all_disease$disease<-factor(list_all_disease$disease,c('HC',disease))
# heatmap=ggplot(data=list_all_disease,aes(x=disease,y=clinical_feature,fill=round(abs(Adjust_R2),2)))+
#   geom_tile()+
#   scale_fill_gradient2(limits=c(0,1),low='#4FBBD8',midpoint =0.35,high='#E84935',space='Lab')+
#   geom_text(aes(disease,clinical_feature,label=round(Adjust_R2,2)),color='black',size=4)+
#   theme_minimal()+
#   theme(panel.grid.major = element_blank(),
#         axis.text.x = element_text(angle=45,vjust=1,size=10,hjust=1,color='black', family = "serif"),
#         axis.text.y = element_text(angle=0,vjust=0,size=10,hjust=1,color='black', family = "serif"))+
#   geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5),color='white',size=1)+
#   #geom_hline(yintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),color='white',size=1)
#   geom_hline(yintercept=c(1:length(unique(list_all_disease$clinical_feature)))-0.5,color='white',size=1)
# print(heatmap)
# dev.off()
# 
# 
# #plot list_all_disease
# png(filename = paste0('heatmap_clinical_association_correlation_R_map_step_',model_sel,'.png'),
#     width = 2500,
#     height = 1500,
#     units = "px",
#     bg = "white",
#     res = 300)
# list_all_disease<-list_all_disease[!is.na(list_all_disease$disease),]
# list_all_disease$clinical_feature<-factor(list_all_disease$clinical_feature,levels=rev(c(var_pre)))
# list_all_disease$disease<-factor(list_all_disease$disease,c('HC',disease))
# heatmap=ggplot(data=list_all_disease,aes(x=disease,y=clinical_feature,fill=round(rvalue,2)))+
#   geom_tile()+
#   scale_fill_gradient2(limits=c(0,1),low='#4FBBD8',midpoint =0.5,high='#E84935',space='Lab')+
#   geom_text(aes(disease,clinical_feature,label=round(rvalue,2)),color='black',size=4)+
#   theme_minimal()+
#   theme(panel.grid.major = element_blank(),
#         axis.text=element_text(size=16),
#         axis.text.x = element_text(angle=0,vjust=0,size=12,hjust=0.5,color='black', family = "serif"),
#         axis.text.y = element_text(angle=0,vjust=0,size=12,hjust=1,color='black', family = "serif"))+
#   geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5),color='white',size=1)+
#   geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5),color='white',size=1)+
#   #geom_hline(yintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),color='white',size=1)
#   geom_hline(yintercept=c(1:12)-0.5,color='white',size=1)+
#   scale_x_discrete(labels = c("MCI" = "MCI","AD" = "AD", "PD" = "PD",
#                               "CSVD" = "CSVD","MS" = "MS",'NMOSD'="NMOSD"))+
#   scale_y_discrete(labels = c("relapses" = "Relapse Number","MOCA" = "MoCA", "BVLT"="BVMT","Post_DBS_UPRDS_off" = "UPRDS Receiving DBS",
#                             "Pre_UPDRS_III_med_on" = "UPRDS Receiving Medicine","Pre_UPDRS_III_med_off" = "UPRDS"))
# print(heatmap)
# dev.off()
# 
# png(filename = paste0('heatmap_clinical_association_correlation_mape_step_',model_sel,'.png'),
#     width = 2500,
#     height = 1500,
#     units = "px",
#     bg = "white",
#     res = 300)
# list_all_disease<-list_all_disease[!is.na(list_all_disease$disease),]
# list_all_disease$clinical_feature<-factor(list_all_disease$clinical_feature,levels=rev(c(var_pre)))
# list_all_disease$disease<-factor(list_all_disease$disease,c('HC',disease))
# heatmap=ggplot(data=list_all_disease,aes(x=disease,y=clinical_feature,fill=round(abs(MAPE),2)))+
#   geom_tile()+
#   scale_fill_gradient2(limits=c(0,2),low='#4FBBD8',midpoint =1,high='#E84935',space='Lab')+
#   geom_text(aes(disease,clinical_feature,label=round(MAPE,2)),color='black',size=4)+
#   theme_minimal()+
#   theme(panel.grid.major = element_blank(),
#         axis.text=element_text(size=16),
#         axis.text.x = element_text(angle=0,vjust=0,size=12,hjust=0.5,color='black', family = "serif"),
#         axis.text.y = element_text(angle=0,vjust=0,size=12,hjust=1,color='black', family = "serif"))+
#   geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5),color='white',size=1)+
#   #geom_hline(yintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),color='white',size=1)
#   geom_hline(yintercept=c(1:length(unique(list_all_disease$clinical_feature)))-0.5,color='white',size=1)+
#   scale_x_discrete(labels = c("MCI" = "MCI","AD" = "AD", "PD" = "PD",
#                               "CSVD" = "CSVD","MS" = "MS",'NMOSD'="NMOSD"))+
#   scale_y_discrete(labels = c("relapses" = "Relapse Number","MOCA" = "MoCA", "BVLT"="BVMT","Post_DBS_UPRDS_off" = "UPRDS Receiving DBS",
#                               "Pre_UPDRS_III_med_on" = "UPRDS Receiving Medicine","Pre_UPDRS_III_med_off" = "UPRDS"))
# print(heatmap)
# dev.off()



#plot list_all_disease for example data
png(filename = paste0('heatmap_clinical_association_step_',model_sel,'.png'),
    width = 2000,
    height = 1480,
    units = "px",
    bg = "white",
    res = 300)
list_all_disease<-list_all_disease[!is.na(list_all_disease$disease),]
list_all_disease$clinical_feature<-factor(list_all_disease$clinical_feature,levels=rev(c(var_pre)))
list_all_disease$disease<-factor(list_all_disease$disease,c('HC',disease))
heatmap=ggplot(data=list_all_disease,aes(x=disease,y=clinical_feature,fill=round(R2,2)))+
  geom_tile()+
  scale_fill_gradient2(limits=c(0,1),low='#4FBBD8',midpoint =0.5,high='#E84935',space='Lab')+
  geom_text(aes(disease,clinical_feature,label=round(R2,2)),color='black',size=4)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle=45,vjust=1,size=10,hjust=1,color='black', family = "serif"),
        axis.text.y = element_text(angle=0,vjust=0,size=10,hjust=1,color='black', family = "serif"))+
  geom_vline(xintercept=c(0.5,1.5,2.5),color='white',size=1)+
  #geom_hline(yintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),color='white',size=1)
  geom_hline(yintercept=c(1:12)-0.5,color='white',size=1)
print(heatmap)
dev.off()

png(filename = paste0('heatmap_clinical_association_adjust_R2_step_',model_sel,'.png'),
    width = 3200,
    height = 1880,
    units = "px",
    bg = "white",
    res = 300)
list_all_disease<-list_all_disease[!is.na(list_all_disease$disease),]
list_all_disease$clinical_feature<-factor(list_all_disease$clinical_feature,levels=rev(c(var_pre)))
list_all_disease$disease<-factor(list_all_disease$disease,c('HC',disease))
heatmap=ggplot(data=list_all_disease,aes(x=disease,y=clinical_feature,fill=round(abs(Adjust_R2),2)))+
  geom_tile()+
  scale_fill_gradient2(limits=c(0,1),low='#4FBBD8',midpoint =0.35,high='#E84935',space='Lab')+
  geom_text(aes(disease,clinical_feature,label=round(Adjust_R2,2)),color='black',size=4)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle=45,vjust=1,size=10,hjust=1,color='black', family = "serif"),
        axis.text.y = element_text(angle=0,vjust=0,size=10,hjust=1,color='black', family = "serif"))+
  geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5),color='white',size=1)+
  #geom_hline(yintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),color='white',size=1)
  geom_hline(yintercept=c(1:length(unique(list_all_disease$clinical_feature)))-0.5,color='white',size=1)
print(heatmap)
dev.off()


#plot list_all_disease
png(filename = paste0('heatmap_clinical_association_correlation_R_map_step_',model_sel,'.png'),
    width = 2500,
    height = 1500,
    units = "px",
    bg = "white",
    res = 300)
list_all_disease<-list_all_disease[!is.na(list_all_disease$disease),]
list_all_disease$clinical_feature<-factor(list_all_disease$clinical_feature,levels=rev(c(var_pre)))
list_all_disease$disease<-factor(list_all_disease$disease,c('HC',disease))
heatmap=ggplot(data=list_all_disease,aes(x=disease,y=clinical_feature,fill=round(rvalue,2)))+
  geom_tile()+
  scale_fill_gradient2(limits=c(0,1),low='#4FBBD8',midpoint =0.5,high='#E84935',space='Lab')+
  geom_text(aes(disease,clinical_feature,label=round(rvalue,2)),color='black',size=4)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        axis.text=element_text(size=16),
        axis.text.x = element_text(angle=0,vjust=0,size=12,hjust=0.5,color='black', family = "serif"),
        axis.text.y = element_text(angle=0,vjust=0,size=12,hjust=1,color='black', family = "serif"))+
  geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5),color='white',size=1)+
  geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5),color='white',size=1)+
  #geom_hline(yintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),color='white',size=1)
  geom_hline(yintercept=c(1:12)-0.5,color='white',size=1)+
  scale_x_discrete(labels = c("MCI" = "MCI","AD" = "AD", "PD" = "PD",
                              "CSVD" = "CSVD","MS" = "MS",'NMOSD'="NMOSD"))+
  scale_y_discrete(labels = c("relapses" = "Relapse Number","MOCA" = "MoCA", "BVLT"="BVMT","Post_DBS_UPRDS_off" = "UPRDS Receiving DBS",
                              "Pre_UPDRS_III_med_on" = "UPRDS Receiving Medicine","Pre_UPDRS_III_med_off" = "UPRDS"))
print(heatmap)
dev.off()

png(filename = paste0('heatmap_clinical_association_correlation_mape_step_',model_sel,'.png'),
    width = 2500,
    height = 1500,
    units = "px",
    bg = "white",
    res = 300)
list_all_disease<-list_all_disease[!is.na(list_all_disease$disease),]
list_all_disease$clinical_feature<-factor(list_all_disease$clinical_feature,levels=rev(c(var_pre)))
list_all_disease$disease<-factor(list_all_disease$disease,c('HC',disease))
heatmap=ggplot(data=list_all_disease,aes(x=disease,y=clinical_feature,fill=round(abs(MAPE),2)))+
  geom_tile()+
  scale_fill_gradient2(limits=c(0,2),low='#4FBBD8',midpoint =1,high='#E84935',space='Lab')+
  geom_text(aes(disease,clinical_feature,label=round(MAPE,2)),color='black',size=4)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        axis.text=element_text(size=16),
        axis.text.x = element_text(angle=0,vjust=0,size=12,hjust=0.5,color='black', family = "serif"),
        axis.text.y = element_text(angle=0,vjust=0,size=12,hjust=1,color='black', family = "serif"))+
  geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5),color='white',size=1)+
  #geom_hline(yintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),color='white',size=1)
  geom_hline(yintercept=c(1:length(unique(list_all_disease$clinical_feature)))-0.5,color='white',size=1)+
  scale_x_discrete(labels = c("MCI" = "MCI","AD" = "AD", "PD" = "PD",
                              "CSVD" = "CSVD","MS" = "MS",'NMOSD'="NMOSD"))+
  scale_y_discrete(labels = c("relapses" = "Relapse Number","MOCA" = "MoCA", "BVLT"="BVMT","Post_DBS_UPRDS_off" = "UPRDS Receiving DBS",
                              "Pre_UPDRS_III_med_on" = "UPRDS Receiving Medicine","Pre_UPDRS_III_med_off" = "UPRDS"))
print(heatmap)
dev.off()
