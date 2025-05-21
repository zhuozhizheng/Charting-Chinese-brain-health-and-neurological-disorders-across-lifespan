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



datapath='/mnt/sda/zzz/Lifespan-main/'#  change to your path
feature_path0='/mnt/sda/zzz/Lifespan_freesurfer_results/Results/Step1_normative_models/V7_DK_no_CI/'# change to your path
savepath='/mnt/sda/zzz/Lifespan_freesurfer_results/Results/Step6_disease_progression/V7_DK/multivariate_analysis_bootstrap_CI'# change to your path

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
str_lab<-'global_and_local_feature_permutation'# change
str_feature<-c('centile');
#feature_path0='/mnt/sda/zzz/Lifespan_freesurfer_results/Results/Step1_normative_models/V7_DK_no_CI/'# change to your path
#savepath='/mnt/sda/zzz/Lifespan_freesurfer_results/Results/Step6_disease_progression/V7_DK/multivariate_analysis_bootstrap_CI'# change to your path
#str_feature<-c('centile');


# feature_path0='/mnt/sda/zzz/Lifespan_freesurfer_results/Results/Step1_normative_models/V7_DK_no_CI_ZZZ/'
# savepath='/mnt/sda/zzz/Lifespan_freesurfer_results/Results/Step6_disease_progression/V7_DK_no_CI_ZZZ/multivariate_analysis_no_CI'
# str_feature<-c('centile');
# 
# feature_path0='/mnt/sda/zzz/Lifespan_freesurfer_results/Results/Step1_normative_models/V7_DK_no_CI/'
# savepath='/mnt/sda/zzz/Lifespan_freesurfer_results/Results/Step6_disease_progression/V7_DK_zscore/multivariate_analysis_no_CI'
# str_feature<-c('zscore');

# feature_path0='/mnt/sda/zzz/Lifespan_freesurfer_results/Results/Step1_normative_models/V7_DK_no_CI/'
# savepath='/mnt/sda/zzz/Lifespan_freesurfer_results/Results/Step6_disease_progression/V7_DK_raw/multivariate_analysis_no_CI'
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

disease<-c("MS","NMOSD")

data1<-data[data$Diagnosis=='HC',];

for(i_dis in disease)
{
  data1<-rbind(data1,data[data$Diagnosis==i_dis,])
  
}

data<-data1;  


# if (!(dir.exists(paste0(savepath,'/',str))))
# {dir.create(paste0(savepath,'/',str))}
# 
# setwd(paste0(savepath,'/',str))

setwd(paste0(savepath))

tem_data1<-data;
list_eta<-NULL
num=0;
list_roc<-list()
list_roc_results<-list()
model<-list()

setwd(savepath);


#var_pre<-c('EDSS_Progress','SPMS_conversion','Followup_relapse')

var_pre<-c('EDSS_Progress','SPMS_conversion')

data[data$EDSS_Progress==1&!is.na(data$EDSS_Progress),'EDSS_Progress']<-1;
data[data$EDSS_Progress==0&!is.na(data$EDSS_Progress),'EDSS_Progress']<-0;
data[data$SPMS_conversion==1&!is.na(data$SPMS_conversion),'SPMS_conversion']<-1;
data[data$SPMS_conversion==0&!is.na(data$SPMS_conversion),'SPMS_conversion']<-0;
# data[data$Followup_relapse>=1&!is.na(data$Followup_relapse),'Followup_relapse']<-1;
# data[data$Followup_relapse==0&!is.na(data$Followup_relapse),'Followup_relapse']<-0;



results=list();
#results=list();
list<-data.frame(matrix(0,length(var_pre)*length(c(disease))*length(feature_col),16))

list_all<-data.frame(matrix(0,length(var_pre),16))

colnames(list)<-c('clinical_feature','MRI_feature','disease','coef',
                  'exp(coef)1','se(coef)','z','Pr(>|z|)','exp(coef)',
                  'exp(-coef)','lower.95','upper.95',
                  'degree','concordance','Score(logrank)test','pvalue_Score(logrank)test')

colnames(list_all)<-c('clinical_feature','MRI_feature','disease','c-index')


boot_num=5000


list_bootstarp<-data.frame(matrix(0,length(var_pre)*boot_num,8))

colnames(list_bootstarp)<-c('clinical_feature','disease','c-index','per_index','Sur1','Sur2','Trans1','Trans2')
  
num=0;num1=0;boot_num0=0;
for(i_var in 1:length(var_pre))
{
  
  tem_data<-data[!is.na(data[,var_pre[i_var]]),];
  
  tem_data[tem_data$Diagnosis=='HC','G']<-'HC'
  tem_data[tem_data$Diagnosis!='HC','G']<-'Disease'
  
  tem_data$G<-factor(tem_data$G,levels=c('HC','Disease'))
  
  
  tem_data$Diagnosis<-factor(tem_data$Diagnosis,
                             levels = c('HC',disease));
  library(ppcor)
  #disease feature associations
  
 #univariate
  for(i_dis in c(disease))
  {
    
    for(i_feature in feature_col)
    {
      num=num+1;
      #print(i_feature)
      #data_ass<-data[data$Diagnosis==i_dis,c('Age','Sex','Followup_time_month','Duration_Month',var_pre[i_var],i_feature)];
      data_ass<-data[data$Diagnosis==i_dis,c('Followup_time_month',var_pre[i_var],i_feature)];
      
      data_ass<-na.omit(data_ass);if(dim(data_ass)[1]==0){next}
      
      library("survival")
      library("survminer")
      
      data_ass_lm<-data_ass;colnames(data_ass_lm)[2:3]<-c('cli_feature','MR_feature');
      
      
      #res.cox <- coxph(Surv(Followup_time_month, cli_feature) ~ Age + Sex + MR_feature+Duration_Month, data =  data_ass_lm)
      res.cox <- coxph(Surv(Followup_time_month, cli_feature) ~ MR_feature, data =  data_ass_lm)
      
      tem_sum<-summary(res.cox)
      
      
      list[num,1]<-var_pre[i_var];
      list[num,2]<-i_feature;
      list[num,3]<-i_dis
      list[num,4]<-tem_sum$coefficients[1,1]
      list[num,5]<-tem_sum$coefficients[1,2]
      list[num,6]<-tem_sum$coefficients[1,3]
      list[num,7]<-tem_sum$coefficients[1,4]
      list[num,8]<-tem_sum$coefficients[1,5]
      list[num,9]<-tem_sum$conf.int[1,1]
      list[num,10]<-tem_sum$conf.int[1,2]
      list[num,11]<-tem_sum$conf.int[1,3]
      list[num,12]<-tem_sum$conf.int[1,4]
      
      list[num,13]<-tem_sum$n
      list[num,14]<-tem_sum$concordance[1]
      list[num,15]<-tem_sum$sctest[1]
      list[num,16]<-tem_sum$sctest[3]
    }
  }
  

  
  #multivariate
  for(i_dis in c(disease))
  {
  num1=num1+1;
  #tem_data0<-tem_data[tem_data$Diagnosis==i_dis,c(var_pre[i_var],'Age','Sex','Followup_time_month','Duration_Month',feature_col,'Diagnosis')];
  tem_data0<-tem_data[tem_data$Diagnosis==i_dis,c(var_pre[i_var],'Followup_time_month',feature_col)];
  tem_data0<-na.omit(tem_data0);if(dim(tem_data0)[1]==0){next}
  tem_data0_zz<-tem_data0;
  #tem_data0<-tem_data0[,c(var_pre[i_var],'Age','Sex','Duration_Month',feature_col)];
  
  #list[,'pFDR']<-p.adjust(list$`Pr(>|z|)`,'fdr')
  
  # tem_list<-list[list$clinical_feature==var_pre[i_var]&
  #                  list$disease==i_dis&list$`Pr(>|z|)`<0.05&
  #                  list$`exp(coef)`<1,]
  
  # tem_list<-list[list$clinical_feature==var_pre[i_var]&
  #                  list$disease==i_dis&list$`Pr(>|z|)`<0.05,]
  
  if(i_dis=='NMOSD')
  {
  tem_list<-list[list$clinical_feature==var_pre[i_var]&
                   list$disease==i_dis&list$`Pr(>|z|)`<0.05,]
  }
  
  if(i_dis=='MS')
  {
    tem_list<-list[list$clinical_feature==var_pre[i_var]&
                     list$disease==i_dis&list$`Pr(>|z|)`<0.05,]
  }
  #tem_list<-tem_list[order(tem_list$`exp(coef)`),]
  
  tem_list<-tem_list[order(tem_list$`Pr(>|z|)`),]
  
  sel_feature<-tem_list$MRI_feature
  width1=1200;height1=720
  
  # if(var_pre[i_var]=='EDSS_Progress'&i_dis=='AQP4Pos_NMOSD')
  # {width1==1200;height1==720}
  # if(var_pre[i_var]=='EDSS_Progress'&i_dis=='MS')
  # {width1==1200;height1==100}
  # 
  # if(var_pre[i_var]=='SPMS_conversion'&i_dis=='MS')
  # {width1==1200;height1==720}
  # if(var_pre[i_var]=='Followup_relapse'&i_dis=='MS')
  # {width1==1200;height1==100}
  
  
  # png(filename = paste0(var_pre[i_var],'_',i_dis,'_clinical_prognosis_cox_error_plot_single_cox.png'), 
  #     width = width1,           
  #     height = height1,          
  #     units = "px",          
  #     bg = "white",          
  #     res = 300) 
  # tem_list$lower.95
  # 
  # tem_list<-tem_list[order(tem_list$`exp(coef)`),]
  # tem_list$MRI_feature<-factor(tem_list$MRI_feature,
  #                              levels = tem_list[order(tem_list$`exp(coef)`,decreasing = T),'MRI_feature'])
  # p<-ggplot(tem_list, aes(y = MRI_feature, x = exp(coef)) )+
  #   geom_pointrange(aes(xmin = lower.95, xmax = upper.95),
  #                   linewidth=1, size=0.6, color="black", shape=16)+
  #   scale_x_continuous(limits=c(0,1),
  #                      breaks=c(0,1))+
  #   geom_vline(xintercept = 1, linewidth=1,linetype = 2,alpha=1,color='grey')+
  #   #geom_vline(xintercept = 0.5, linewidth=0.6,linetype = 2,alpha=1,color='grey')+
  # theme_classic()
  # print(p)
  # dev.off()
  
  
  
  # sel_feature_col<-tem_list[,'MRI_feature']
  # 
  # tem_data0<-tem_data0[,c(var_pre[i_var],'Followup_time_month',sel_feature_col)];
  
  # tem_data0<-tem_data0[,c(var_pre[i_var],'Age','Sex',feature_col)];
  # 
  # X=as.matrix(tem_data0[,2:length(colnames(tem_data0))]);
  # Y=as.matrix(Surv(tem_data0_zz$Followup_time_month,tem_data0[,var_pre[i_var]]));
  # 
  # mod_cv <- cv.glmnet(x=X,
  #                     y=Y,
  #                     family="cox", nfolds = 10,
  #                     intercept = F, alpha=1)
  # 
  # best_lambda <- mod_cv$lambda.min
  # 
  # best_model <- glmnet(x=X,
  #                      y=Y,
  #                      family="cox",
  #                      alpha=1,lambda = best_lambda)
  # tem_coef<-data.frame(matrix(coef(best_model)))
  # 
  # rownames(tem_coef)<-c(colnames(tem_data0)[2:(dim(tem_data0)[2])]);
  # tem_coef
  # index_feature=which(tem_coef!=0); if(length(index_feature)==0){next}
  # 
  # index_feature<-index_feature[1:length(index_feature)]
  # sel_feature<-rownames(tem_coef)[index_feature]
  # 
  # sel_feature<-c(var_pre[i_var],'Followup_time_month',sel_feature);

  
  
  tem_data0<-tem_data0_zz[,c(var_pre[i_var],'Followup_time_month',sel_feature)]
  
  
  
  colnames(tem_data0)[1]<-c('cli_feature')
  
  
  res.cox <- coxph(Surv(Followup_time_month, cli_feature) ~ ., data =  tem_data0)
  tem_sum<-summary(res.cox)
  
  
  model1<-coxph(Surv(Followup_time_month, cli_feature)~.,data=tem_data0);
  
  library(autoReg)
  
  
  #install.packages("export")
  library(rrtable)
  library(export)
  
  model<-step(model1,direction="forward",scope=formula(Surv(Followup_time_month, cli_feature)~.),trace=0)
  
 
  
  saveRDS(model,paste0(i_dis,'_',var_pre[i_var],'.rds'));
  
  table2doc(model,file=paste0(i_dis,'_',var_pre[i_var]), digits = 1, digitspvals = 3, add.rownames =TRUE)
  
  model2<-autoReg(model1,uni=TRUE,threshold=0.05, final= TRUE) 
  saveRDS(model2,paste0(i_dis,'_',var_pre[i_var],'_autoReg.rds'));
  
  table2doc(model2,file=paste0(i_dis,'_',var_pre[i_var],'_autoReg'), digits = 1, digitspvals = 3, add.rownames =TRUE)
  
  #table2pptx(model)  #导出到ppt，可编辑数据
  
  #sel_feature1<-model$name[model$coefficients !=""]
  
  sel_feature1<-names(model$coefficients)
  
  sel_feature1<-str_remove(sel_feature1,'`')
  sel_feature1<-str_remove(sel_feature1,'`')
  
  tem_data0<-tem_data0_zz[,c(var_pre[i_var],'Followup_time_month',sel_feature1)]
  
  #tem_data0<-tem_data0_zz[,c(var_pre[i_var],'Followup_time_month',sel_feature)]
  
  
  
  
  
  
  #addd#addd#addd#addd#addd#addd#addd#addd#addd#addd#addd#addd#addd#addd#addd
  library(glmnet)
  colnames(tem_data0)[1]<-c('cli_feature')
  x <- model.matrix(~ . - Followup_time_month - cli_feature, data = tem_data0)
  y <- Surv(tem_data0$Followup_time_month, tem_data0$cli_feature)
  cvfit <- cv.glmnet(x, y, family = "cox")

  if(i_dis=='NMOSD')
  {
  l.coef2<-coef(cvfit$glmnet.fit,s=0.03,exact=F)
  }
  
  if(i_dis=='MS')
  
  {
    l.coef2<-coef(cvfit$glmnet.fit,s=cvfit$lambda.min,exact=F)
  }
  
  tem_feature<-l.coef2@Dimnames[[1]]
  sel_feature<-tem_feature[l.coef2@i+1];
  tem_data0<-tem_data0_zz[,c(var_pre[i_var],'Followup_time_month',sel_feature)]
  #addd#addd#addd#addd#addd#addd#addd#addd#addd#addd#addd#addd#addd#addd#addd
  
  
  tem_list_fig<-tem_list
  
  rownames(tem_list_fig)<-tem_list_fig$MRI_feature
  
  tem_list_fig<-tem_list_fig[sel_feature1,]
  # sel_feature_col_fig<-str_remove(tem_list_fig$MRI_feature,'`')
  # sel_feature_col_fig<-str_remove(sel_feature_col_fig,'`')
  
  
  library(raincloudplots)
  library(ggwordcloud)
  setwd(savepath)
  png(filename = paste0(i_dis,'_',var_pre[i_var],'_predictive_features_importantance_plot.png'),
      width = 1480,
      height = 1480,
      units = "px",
      bg = "white",
      res = 300)
  
  p<-ggplot(tem_list_fig, aes(label = MRI_feature,
                              size = -log(`Pr(>|z|)`),color = concordance)) +
    geom_text_wordcloud(eccentricity=1) +
    scale_size_area(max_size = 4)+
    theme_minimal()+
    scale_color_gradient(limits=c(0.5,1),low = "grey",high = "#990000")
  print(p)
  dev.off()
  
  
  
  
  

  
  colnames(tem_data0)[1]<-c('cli_feature')
  ###cross-validation
  library(caret)
  set.seed(7);
  folds<-createFolds(y=tem_data0$cli_feature,k=10);

  data2<-tem_data0
  cindex_results_test=NULL;
  for(i in 1:10)
  {
    fold_test<-tem_data0[folds[[i]],]
    fold_train<-tem_data0[-folds[[i]],]

    print(paste0('***fold',i,'***'))

    model<-coxph(Surv(Followup_time_month, cli_feature)~.,data=fold_train);

    Cin<-concordance(model,newdata=fold_test);

    cindex_results_test<-append(cindex_results_test,as.numeric(Cin$concordance));


    # fold_pre<-glm(label~.,family = 'gaussian',data=fold_train)
    fold_predict<-predict(model,type='survival',newdata=fold_test);
    #fold_predict<-predict(model,type='risk',newdata=fold_test);

    data2[folds[[i]],'Predict']<-fold_predict;

  }
  
  # ft3<-autoReg(coxmod,uni=TRUE,threshold=0.05, final= TRUE) 
  # myft(ft3)
  
  cindex_mean<-mean(cindex_results_test);
  # cindex_mean<-mean(cindex_results_test);
  # 
  # real<-tem_data0$label;
  # tem_data0[,'predict']<-data2$Predict
  # predict.<-data2$Predict
  # 
  # real[real==0]=1;
  
  
  # real<-real[!is.infinite(predict.)]
  # predict.<-predict.[!is.infinite(predict.)]
  # predict.<-predict.[!is.infinite(real)]
  # real<-real[!is.infinite(real)]
  
  list_all[num1,1]<-var_pre[i_var];
  list_all[num1,2]<-'feature';
  list_all[num1,3]<-i_dis
  list_all[num1,4]<-cindex_mean
  
  
  library("survminer")
 

  
  model<-coxph(Surv(Followup_time_month, cli_feature)~.,data=tem_data0);
  fold_predict<-predict(model,type='survival',newdata=tem_data0);
  data2[,'Predict']<-fold_predict;
  # 
  
  data2[data2$Predict<=mean(data2$Predict),'Subgroup']=0;
  data2[data2$Predict>mean(data2$Predict),'Subgroup']=1;
  
  # if(i_dis=='AQP4Pos_NMOSD'){
  # data2[data2$Predict<=median(data2$Predict),'Subgroup']=0;
  # data2[data2$Predict>median(data2$Predict),'Subgroup']=1;}


  
  fit <- survfit(Surv(Followup_time_month, cli_feature) ~ Subgroup, data = data2)
  list_sur<-list();
  summary(fit)

  summary(fit, times=5*12)$surv  # 生存概率
  1 - summary(fit, times=5*12)$surv  # 转化率(事件发生概率)
  
  list_sur[['sur_pro']]<-summary(fit, times=5*12)$surv
  list_sur[['trans_pro']]<-1-summary(fit, times=5*12)$surv
  
  saveRDS(summary(fit),paste0(var_pre[i_var],'_',i_dis,'_clinical_prognosis.rds'))
  saveRDS(list_sur,paste0(var_pre[i_var],'_',i_dis,'_transitrion_rate.rds'))
  
  png(filename = paste0(var_pre[i_var],'_',i_dis,'_clinical_prognosis.png'), 
      width = 1200,           
      height = 1200,          
      units = "px",          
      bg = "white",          
      res = 300) 
  
  p<-ggsurvplot(fit,
             pval = TRUE, #conf.int = TRUE,
             #risk.table = TRUE, # 添加风险表
             #risk.table.col = "strata", # 根据分层更改风险表颜色
             #linetype = "strata", # 根据分层更改线型
             surv.median.line = "hv", # 同时显示垂直和水平参考线
             #ggtheme = theme_bw(), # 更改ggplot2的主题
             legend.labs = c("High risk", "Low risk"),
             palette = c("#E84935","#4FBBD8"))#定义颜色
  print(p)
  dev.off()
  
  
  
  

  tem_data1<-tem_data0
  
  
  # boot_num=5000
  # 
  # 
  # list_bootstarp<-data.frame(matrix(0,length(var_pre)*boot_num,16))
  
  
  library(rsample)
  
  # 创建bootstrap重采样
  boot_samples <- bootstraps(tem_data1, times = boot_num)
 
  for(i_boot in 1:boot_num)
  {
  print(i_boot)
  boot_num0<-boot_num0+1
  
  tem_data0<-data.frame(boot_samples$splits[[i_boot]])
  tem_data0<-na.omit(tem_data0)
  
  # set.seed(sample(nrow(tem_data1),1,replace = TRUE));
  # tem_data0<-  tem_data1[sample(nrow(tem_data1), size =nrow(tem_data1) , replace = TRUE), ]
  ###cross-validation
  library(caret)
  
  folds<-createFolds(y=tem_data0$cli_feature,k=10);
  
  data2<-tem_data0
  cindex_results_test=NULL;
  
  div<-try(
  for(i in 1:10)
  {
    fold_test<-tem_data0[folds[[i]],]
    fold_train<-tem_data0[-folds[[i]],]
    
    #print(paste0('***fold',i,'***'))
    
    model<-coxph(Surv(Followup_time_month, cli_feature)~.,data=fold_train);
    
    Cin<-concordance(model,newdata=fold_test);
    
    cindex_results_test<-append(cindex_results_test,as.numeric(Cin$concordance));
    
    
    # fold_pre<-glm(label~.,family = 'gaussian',data=fold_train)
    fold_predict<-predict(model,type='survival',newdata=fold_test);
    #fold_predict<-predict(model,type='risk',newdata=fold_test);
    
    data2[folds[[i]],'Predict']<-fold_predict;
    
  })
  
  if(inherits(div, "try-error")) {
    next
  } 
  
  

  
  # ft3<-autoReg(coxmod,uni=TRUE,threshold=0.05, final= TRUE) 
  # myft(ft3)
  
  cindex_mean<-mean(cindex_results_test);
  # cindex_mean<-mean(cindex_results_test);
  # 
  # real<-tem_data0$label;
  # tem_data0[,'predict']<-data2$Predict
  # predict.<-data2$Predict
  # 
  # real[real==0]=1;
  
  
  # real<-real[!is.infinite(predict.)]
  # predict.<-predict.[!is.infinite(predict.)]
  # predict.<-predict.[!is.infinite(real)]
  # real<-real[!is.infinite(real)]
  

  
  list_bootstarp[boot_num0,1]<-var_pre[i_var]
  list_bootstarp[boot_num0,2]<-i_dis
  list_bootstarp[boot_num0,3]<-cindex_mean
  list_bootstarp[boot_num0,4]<-i_boot
  
  library("survminer")
  
  model<-coxph(Surv(Followup_time_month, cli_feature)~.,data=tem_data0);
  fold_predict<-predict(model,type='survival',newdata=tem_data0);
  data2[,'Predict']<-fold_predict;
  # 
  
  data2[data2$Predict<=mean(data2$Predict,na.rm=T)&!is.na(data2$Predict),'Subgroup']=0;
  data2[data2$Predict>mean(data2$Predict,na.rm=T)&!is.na(data2$Predict),'Subgroup']=1;
  
  # if(i_dis=='AQP4Pos_NMOSD'){
  # data2[data2$Predict<=median(data2$Predict),'Subgroup']=0;
  # data2[data2$Predict>median(data2$Predict),'Subgroup']=1;}
  
  fit <- survfit(Surv(Followup_time_month, cli_feature) ~ Subgroup, data = data2)
  
  summary(fit)
  
  tem_sum<-summary(fit, times=5*12)$surv  # 生存概率
  list_bootstarp[boot_num0,5]<-tem_sum[1]
  list_bootstarp[boot_num0,6]<-tem_sum[2]
  
  tem_sum<-1 - summary(fit, times=5*12)$surv  # 转化率(事件发生概率)
  
  list_bootstarp[boot_num0,7]<-tem_sum[1]
  list_bootstarp[boot_num0,8]<-tem_sum[2]

  
  }
  
  
  
  

  
  
  }
}

results$signal_feature_cox<-list;
results$all_feature_lm<-list_all;

saveRDS(results,'clinical_progosis_ms_nmo.rds')

saveRDS(list_bootstarp,'clinical_progosis_ms_nmo_bootstrap_CI.rds')