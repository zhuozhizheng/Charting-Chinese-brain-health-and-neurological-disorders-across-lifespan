rm(list=ls())
#set your own directory
datapath='E:/Lifespan_freesurfer_results/Github/Source-codes/'  # Change the directory where you save the Source-codes  
clinical_datapath='E:/Lifespan_freesurfer_results/Github/Datasets/Dataset-individual/Clinical_vars.csv'  # Change the directory where you save the clinical variables  
MR_datapath='E:/Lifespan_freesurfer_results/Github/Datasets/Dataset-individual/MR_measures.xlsx'  # Change the directory where you save the MR measures  
savepath1='E:/Lifespan_freesurfer_results/Github/Models/GAMLSS/DK' #Chnage the path you have saved your normative models in "Normative-model-fit.R"; here we used our models for example
savepath='E:/Lifespan_freesurfer_results/Github/Test_results/Individual';# Create and determine the directory where you would save the results  

setwd(datapath)
source("100.common-variables.r")
source("101.common-functions.r")
source("ZZZ_function.R")
source("300.variables.r")
source("301.functions.r")

library(dplyr)
library(stringr)

var<-c('Global.table',
       'aseg.vol.table',
       'lh.aparc.volume.table','rh.aparc.volume.table',
       'lh.aparc.thickness.table','rh.aparc.thickness.table',
       'lh.aparc.area.table','rh.aparc.area.table')



for(sheet in var) 
{ 
  
  setwd(datapath)
  MRI <- read_excel(MR_datapath,sheet=sheet)

  MRI<-as.data.frame(MRI)
  
  rownames(MRI)<-paste0(MRI$Freesurfer_Path2,MRI$Freesurfer_Path3)
  if(str_detect(sheet,'aparc'))
  {
    tem_feature<-colnames(MRI)[c(2:35)];
  }
  
  str=sheet;
  
  
  if(sheet=="aseg.vol.table")
  {
    MRI[,'cerebellum_WM']<-MRI$Left.Cerebellum.White.Matter+MRI$Right.Cerebellum.White.Matter
    MRI[,'cerebellum_GM']<-MRI$Left.Cerebellum.Cortex+MRI$Right.Cerebellum.Cortex
    MRI[,'cerebellum_total']<-MRI[,'cerebellum_WM']+MRI[,'cerebellum_GM'];
    MRI[,'CC']<-MRI$CC_Anterior+MRI$CC_Central+MRI$CC_Mid_Anterior+
      MRI$CC_Mid_Posterior+MRI$CC_Posterior
    
    tem_feature<-colnames(MRI)[c(6:9,12:14,16:17,24:31,69:71)];
    
    
    str=sheet;
    
  }
  
  if(sheet=='Global.table')
  {
    
    
    for(col in colnames(MRI)[2:19])
    {
      MRI[,col]<-as.numeric(MRI[,col])
    }
    
    for(i in 1:dim(MRI)[1])
    {
      MRI[i,'mean_thickness']<-
        (MRI[i,'lhMeanThickness']*MRI[i,'lhVertex']+
           MRI[i,'rhMeanThickness']*MRI[i,'rhVertex'])/(MRI[i,'lhVertex']+MRI[i,'rhVertex'])
      
      MRI[i,'total_surface_arrea']<-MRI[i,'lh_totaISA2']+MRI[i,'rh_totaISA2'];
    }
    
    tem_feature<-colnames(MRI)[c(2,3,4,5,13,23,24)];
    
    str='Global_feature';
    
  }
  
  
  if (!(dir.exists(paste0(savepath,'/',str))))
  {dir.create(paste0(savepath,'/',str))}
  
  setwd(paste0(savepath,'/',str))
  
  Z_data<-list();
  Quant_data<-list()
  
  
  for(i in tem_feature[1:length(tem_feature)])
  {
    
    print(i)
    setwd(paste0(savepath1,'/',str))
    rdsfile<-paste0(i,'_normative_model.rds')
    
    if(i=='Brain.Stem'){rdsfile='BrainStem_normative_model.rds'}
    if(i=='cerebellum_total'){rdsfile='Cerebellum_normative_model.rds'}
    if(i=='mean_thickness'){rdsfile='Cortical_thickness_normative_model.rds'}
    if(i=='total_surface_arrea'){rdsfile='Total_surface_area_normative_model.rds'}
    
    if(file.exists(rdsfile)){
      
      results<-readRDS(rdsfile)
    
    #for each feature, we should load the oirginal clinical information
    setwd(datapath)  
    data1<-read.csv(clinical_datapath,header=TRUE);
    data1$Site_ZZZ<-paste0(data1$Province,data1$Center,data1$Manufacturer)
    rownames(data1)<-paste0(data1$Freesufer_Path2,data1$Freesufer_Path3)
    
    setwd(paste0(savepath,'/',str))
   
    inter_row<-intersect(rownames(data1),rownames(MRI))
    data1=cbind(data1[inter_row,],MRI[inter_row,i])
  
    colnames(data1)[dim(data1)[2]]=c('tem_feature')
    
    
    data1$Site_ZZZ<-as.factor(data1$Site_ZZZ)
    
    data1$Sex<-as.factor(data1$Sex)
    
    data1$Sex<-factor(data1$Sex,levels=c('Female','Male'))
    
    
    data1<-data1[order(data1$Age),]
    data1[,'feature']<-data1$tem_feature
    
    

    data1<-data1[,c('Age','Sex','Site_ZZZ','feature','tem_feature','Diagnosis')]
    
    all_data1<-data1
    all_data1$Sex<-as.factor(all_data1$Sex)
    all_data1$Site_ZZZ<-as.factor(all_data1$Site_ZZZ)
    
  
 
    
    m2<-results$m2
    m0<-results$m0
   
    Z_score_sum<-NULL;
    Quant_score_sum<-NULL
    
    model1<-m2;
    
    
    
    tem_rnd<-matrix(0,dim(all_data1)[1],1);
    Model.Frame<-model.frame(formula = model1$mu.formula,data=all_data1);
    Model.Matrix<-model.matrix(model1$mu.formula,Model.Frame)
    Fit.fix<-matrix(model1$mu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
    if(!is.null(model1$mu.coefSmo[[1]]))
    {Fit.fix[length(Fit.fix)]=0;
    for(iz in 1:dim(all_data1)[1]){
      if(all_data1$Site_ZZZ[iz] %in% names(model1$mu.coefSmo[[1]]$coef)){
        tem_rnd[iz]<-model1$mu.coefSmo[[1]]$coef[as.character(all_data1$Site_ZZZ[iz])]
      }else{tem_rnd[iz]<-mean(model1$mu.coefSmo[[1]]$coef)}
    }
    }else{print('no random effects for this term')};
    
    
    mu<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))
    
    

    tem_rnd<-matrix(0,dim(all_data1)[1],1);
    Model.Frame<-model.frame(formula = model1$sigma.formula,data=all_data1);
    Model.Matrix<-model.matrix(model1$sigma.formula,Model.Frame)
    Fit.fix<-matrix(model1$sigma.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
    if(!is.null(model1$sigma.coefSmo[[1]]))
    {Fit.fix[length(Fit.fix)]=0;
    for(iz in 1:dim(all_data1)[1]){
      if(all_data1$Site_ZZZ[iz] %in% names(model1$sigma.coefSmo[[1]]$coef)){
        tem_rnd[iz]<-model1$sigma.coefSmo[[1]]$coef[as.character(all_data1$Site_ZZZ[iz])]
      }else{tem_rnd[iz]<-mean(model1$sigma.coefSmo[[1]]$coef)}
    }
    }else{print('no random effects for this term')};
    
    sigma<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))

    

    Model.Frame<-model.frame(formula = model1$nu.formula,data=all_data1);
    Model.Matrix<-model.matrix(model1$nu.formula,Model.Frame)
    Fit.fix<-matrix(model1$nu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
    if(!is.null(model1$nu.coefSmo[[1]]))
    {Fit.fix[length(Fit.fix)]=0;
    }else{print('no random effects for this term')};
    
    nu<-as.vector(Model.Matrix[,colnames(Model.Matrix)[1]]*Fit.fix[1])
    

    
    if(length(mu)!=dim(all_data1)[1])
    {
      print("Error, Please Check Data!!!")
    }
    
    
    Z_score_sum<-zzz_cent(obj=model1,type=c("z-scores"),mu=mu,sigma=sigma,nu=nu,
                          xname = 'Age',xvalues=all_data1$Age,yval=all_data1$tem_feature,
                          calibration=FALSE,lpar=3)
    
    Quant_score_sum<-zzz_cent(obj=model1,type=c("z-scores"),mu=mu,sigma=sigma,nu=nu,
                              xname = 'Age',xvalues=all_data1$Age,yval=all_data1$tem_feature,
                              calibration=FALSE,lpar=3,cdf=TRUE)
    
    
    Z_score_sum<-data.frame(Z_score_sum);
    colnames(Z_score_sum)<-c('Z_score');
    rownames(Z_score_sum)<-rownames(all_data1)
    
    Quant_score_sum<-data.frame(Quant_score_sum);
    colnames(Quant_score_sum)<-c('Quant_score');
    rownames(Quant_score_sum)<-rownames(all_data1)
    
    Z_data[[i]]<-Z_score_sum
    Quant_data[[i]]<-Quant_score_sum
    
    
    results_new<-list();
    results_new$m2<-m2
    results_new$Zscore<-Z_data
    results_new$Quant_data<-Quant_data
    results_new$data1<-data1
    results_new$all_data<-all_data1
    results_new$str<-str
    results_new$i<-i

    
    setwd(paste0(savepath,'/',str))
    
    saveRDS(results_new,paste0(str,'_',i,'_loop_our_model_individual.rds'))
    
  }
  
}

}

