rm(list=ls())

#set your own directory
datapath='E:/Lifespan_freesurfer_results/Github/Source-codes/'  # Change the directory where you save the Source-codes  
clinical_datapath='E:/Lifespan_freesurfer_results/Github/Datasets/Dataset-new/Clinical_vars.csv'  # Change the directory where you save the clinical variables  
MR_datapath='E:/Lifespan_freesurfer_results/Github/Datasets/Dataset-new/MR_measures.xlsx'  # Change the directory where you save the MR measures  
modelpath<-'E:/Lifespan_freesurfer_results/Github/Models/GAMLSS/DK' #Chnage the path you have saved your normative models in "Normative-model-fit.R"
savepath='E:/Lifespan_freesurfer_results/Github/Test_results/Newsite';# Create and determine the directory where you would save the results  



setwd(datapath)
source("100.common-variables.r")
source("101.common-functions.r")
source("ZZZ_function.R")
source("300.variables.r")
source("301.functions.r")

var<-c('Global.table',
       'aseg.vol.table',
       'lh.aparc.volume.table','rh.aparc.volume.table',
       'lh.aparc.thickness.table','rh.aparc.thickness.table',
       'lh.aparc.area.table','rh.aparc.area.table')


for(sheet in var) #for local feature
{ 
  

    setwd(datapath)
    MRI <- read_excel(MR_datapath,sheet=sheet)
    
    MRI<-MRI[MRI$Freesurfer_Path2!='Epilepsy_dicom_info_nii'&
               MRI$Freesurfer_Path2!='ET'&
               MRI$Freesurfer_Path2!='Guojibu_HC_nii',]
    MRI<-MRI[!is.na(MRI$Freesurfer_Path3),]
    MRI<-as.data.frame(MRI)
    
    rownames(MRI)<-paste0(MRI$Freesurfer_Path2,MRI$Freesurfer_Path3)
    
  
    tem_feature<-colnames(MRI)[c(2:35)];
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
  #tem_feature<-colnames(MRI)[c(5,13,23,24)];

  str='Global_feature';#str_lab='Global_feature';
  
}

  
if (!(dir.exists(paste0(savepath,'/',str))))
{dir.create(paste0(savepath,'/',str))}

setwd(paste0(savepath,'/',str))



Z_data<-list();
Quant_data<-list()


for(i in tem_feature[1:length(tem_feature)])
{

print(i)
setwd(paste0(modelpath,'/',str))

rdsfile<-paste0(i,'_normative_model.rds')

if(i=='Brain.Stem'){rdsfile='BrainStem_normative_model.rds'}
if(i=='cerebellum_total'){rdsfile='Cerebellum_normative_model.rds'}
if(i=='mean_thickness'){rdsfile='Cortical_thickness_normative_model.rds'}
if(i=='total_surface_arrea'){rdsfile='Total_surface_area_normative_model.rds'}

if(file.exists(rdsfile)){print('file exist');
  
  results<-readRDS(rdsfile)}
  
  
m2<-results$m2;
m0<-results$m0
#for each feature, we should load the oirginal clinical information
setwd(datapath)  
data1<-read.csv(clinical_datapath,header=TRUE);
data1$Site_ZZZ<-paste0(data1$Province,data1$Center,data1$Manufacturer)

library(dplyr)

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

#remove the extreme values
data1<-data1[!is.na(data1$tem_feature),]
data1<-data1[data1$feature>(mean(data1$feature)-3*sd(data1$feature))&
               data1$feature<(mean(data1$feature)+3*sd(data1$feature)),]


#select the columns
data1<-data1[,c('Age','Sex','Site_ZZZ','tem_feature','feature')]

data2<-data1;

#for mu
tem_rnd<-matrix(0,dim(data2)[1],1);
Model.Frame<-model.frame(formula = m2$mu.formula,data=data2);
Model.Matrix<-model.matrix(m2$mu.formula,Model.Frame)
Fit.fix<-matrix(m2$mu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))

if(!is.null(m2$mu.coefSmo[[1]]))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data2)[1]){
  tem_rnd[iz]<-0
}
}else{print('no random effects for this term')};

mufix<-as.vector(Model.Matrix %*% Fit.fix)


#for sigma
tem_rnd<-matrix(0,dim(data2)[1],1);
Model.Frame<-model.frame(formula = m2$sigma.formula,data=data2);
Model.Matrix<-model.matrix(m2$sigma.formula,Model.Frame)
Fit.fix<-matrix(m2$sigma.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
if(!is.null(m2$sigma.coefSmo[[1]]))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data2)[1]){
  tem_rnd[iz]<-0
}
}else{print('no random effects for this term')};

sigmafix<-as.vector(Model.Matrix %*% Fit.fix)

data2[,'mufix']<-mufix
data2[,'sigmafix']<-sigmafix

data2 <- na.omit(data2)

con=gamlss.control()
m2_cal<-gamlss(formula=feature~offset(mufix)+random(Site_ZZZ),
               sigma.formula = feature~offset(sigmafix)+random(Site_ZZZ),
           control=con,
           family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
           data=data2)


model_cal<-m2;
for(new_site in names(m2_cal$mu.coefSmo[[1]][1]$coef))
{
  if(!is.null(m2$mu.coefSmo[[1]][1]$coef))
  {
    model_cal$mu.coefSmo[[1]][1]$coef[new_site]<-m2_cal$mu.coefSmo[[1]][1]$coef[new_site]
  }
  if(!is.null(m2$sigma.coefSmo[[1]][1]$coef))
  {
    model_cal$sigma.coefSmo[[1]][1]$coef[new_site]<-m2_cal$sigma.coefSmo[[1]][1]$coef[new_site]
  }
  
}

setwd(paste0(savepath,'/',str))
saveRDS(model_cal, paste0(str,'_',i,'_loop_our_model_new_site_calibrated.rds'))
#plot calibrate models
#both male and female
#original model
model1<-m2;
num_length=5000
if(!is.null(model1$mu.coefSmo[[1]]$coef))
{
data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
} else
{
data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}
data4 <- do.call( what=expand.grid, args=data3 )

data4[,'feature']<-0
tem_rnd<-matrix(0,dim(data4)[1],1);
Model.Frame<-model.frame(formula = model1$mu.formula,data=data4);
Model.Matrix<-model.matrix(model1$mu.formula,Model.Frame)
Fit.fix<-matrix(model1$mu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
if(!is.null(model1$mu.coefSmo[[1]]))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data4)[1]){
  if(data4$Site_ZZZ[iz] %in% names(model1$mu.coefSmo[[1]]$coef)){
    tem_rnd[iz]<-model1$mu.coefSmo[[1]]$coef[as.character(data4$Site_ZZZ[iz])]
  }else{tem_rnd[iz]<-mean(model1$mu.coefSmo[[1]]$coef)}
}
}else{print('no random effects for this term')};


mu0<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))



if(!is.null(model1$sigma.coefSmo[[1]]$coef))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}
data4 <- do.call( what=expand.grid, args=data3 )

data4[,'feature']<-0
tem_rnd<-matrix(0,dim(data4)[1],1);
Model.Frame<-model.frame(formula = model1$sigma.formula,data=data4);
Model.Matrix<-model.matrix(model1$sigma.formula,Model.Frame)
Fit.fix<-matrix(model1$sigma.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
if(!is.null(model1$sigma.coefSmo[[1]]))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data4)[1]){
  if(data4$Site_ZZZ[iz] %in% names(model1$sigma.coefSmo[[1]]$coef)){
    tem_rnd[iz]<-model1$sigma.coefSmo[[1]]$coef[as.character(data4$Site_ZZZ[iz])]
  }else{tem_rnd[iz]<-mean(model1$sigma.coefSmo[[1]]$coef)}
}
}else{print('no random effects for this term')};

sigma0<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))



if(!is.null(model1$nu.coefSmo[[1]]$coef))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}

data4 <- do.call( what=expand.grid, args=data3 )

Model.Frame<-model.frame(formula = model1$nu.formula,data=data4);
Model.Matrix<-model.matrix(model1$nu.formula,Model.Frame)
Fit.fix<-matrix(model1$nu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
if(!is.null(model1$nu.coefSmo[[1]]))
{Fit.fix[length(Fit.fix)]=0;
}else{print('no random effects for this term')};

nu0<-as.vector(Model.Matrix[,colnames(Model.Matrix)[1]]*Fit.fix[1])



tem_par<-mu0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;

if(Seg>1)
{
for(Seg1 in c(2:Seg))
{
  par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
}
par=par/Seg
}
mu=par


tem_par<-sigma0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;
if(Seg>1)
{
for(Seg1 in c(2:Seg))
{
  par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
}
par=par/Seg
}
sigma=par


tem_par<-nu0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;
if(Seg>1)
{
for(Seg1 in c(2:Seg))
{
  par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
}
par=par/Seg
}
nu=par


p2<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
             cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
             calibration=FALSE,lpar=3)
p2[,'sigma']<-sigma


#calibrate model
model1<-model_cal;
num_length=5000
if(!is.null(model1$mu.coefSmo[[1]]$coef))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}
data4 <- do.call( what=expand.grid, args=data3 )

data4[,'feature']<-0
tem_rnd<-matrix(0,dim(data4)[1],1);
Model.Frame<-model.frame(formula = model1$mu.formula,data=data4);
Model.Matrix<-model.matrix(model1$mu.formula,Model.Frame)
Fit.fix<-matrix(model1$mu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
if(!is.null(model1$mu.coefSmo[[1]]))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data4)[1]){
  if(data4$Site_ZZZ[iz] %in% names(model1$mu.coefSmo[[1]]$coef)){
    tem_rnd[iz]<-model1$mu.coefSmo[[1]]$coef[as.character(data4$Site_ZZZ[iz])]
  }else{tem_rnd[iz]<-mean(model1$mu.coefSmo[[1]]$coef)}
}
}else{print('no random effects for this term')};


mu0<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))

#mu0 <- predict(model1, newdata = data4, type = "response", what = "mu")


if(!is.null(model1$sigma.coefSmo[[1]]$coef))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}
data4 <- do.call( what=expand.grid, args=data3 )

data4[,'feature']<-0
tem_rnd<-matrix(0,dim(data4)[1],1);
Model.Frame<-model.frame(formula = model1$sigma.formula,data=data4);
Model.Matrix<-model.matrix(model1$sigma.formula,Model.Frame)
Fit.fix<-matrix(model1$sigma.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
if(!is.null(model1$sigma.coefSmo[[1]]))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data4)[1]){
  if(data4$Site_ZZZ[iz] %in% names(model1$sigma.coefSmo[[1]]$coef)){
    tem_rnd[iz]<-model1$sigma.coefSmo[[1]]$coef[as.character(data4$Site_ZZZ[iz])]
  }else{tem_rnd[iz]<-mean(model1$sigma.coefSmo[[1]]$coef)}
}
}else{print('no random effects for this term')};

sigma0<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))
#sigma0 <- predict(model1, newdata = data4, type = "response", what = "sigma")


if(!is.null(model1$nu.coefSmo[[1]]$coef))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}

data4 <- do.call( what=expand.grid, args=data3 )

Model.Frame<-model.frame(formula = model1$nu.formula,data=data4);
Model.Matrix<-model.matrix(model1$nu.formula,Model.Frame)
Fit.fix<-matrix(model1$nu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
if(!is.null(model1$nu.coefSmo[[1]]))
{Fit.fix[length(Fit.fix)]=0;
}else{print('no random effects for this term')};

nu0<-as.vector(Model.Matrix[,colnames(Model.Matrix)[1]]*Fit.fix[1])

#nu0 <- predict(model1, newdata = data4, type = "response", what = "nu")

tem_par<-mu0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;

if(Seg>1)
{
  for(Seg1 in c(2:Seg))
  {
    par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
  }
  par=par/Seg
}
mu=par


tem_par<-sigma0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;
if(Seg>1)
{
  for(Seg1 in c(2:Seg))
  {
    par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
  }
  par=par/Seg
}
sigma=par


tem_par<-nu0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;
if(Seg>1)
{
  for(Seg1 in c(2:Seg))
  {
    par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
  }
  par=par/Seg
}
nu=par


p2_cal<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
             cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
             calibration=FALSE,lpar=3)
p2_cal[,'sigma']<-sigma





#site_specific model
model1<-model_cal;
num_length=5000
if(!is.null(model1$mu.coefSmo[[1]]$coef))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=unique(data1$Site_ZZZ))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}
data4 <- do.call( what=expand.grid, args=data3 )

data4[,'feature']<-0
tem_rnd<-matrix(0,dim(data4)[1],1);
Model.Frame<-model.frame(formula = model1$mu.formula,data=data4);
Model.Matrix<-model.matrix(model1$mu.formula,Model.Frame)
Fit.fix<-matrix(model1$mu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
if(!is.null(model1$mu.coefSmo[[1]]))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data4)[1]){
  if(data4$Site_ZZZ[iz] %in% names(model1$mu.coefSmo[[1]]$coef)){
    tem_rnd[iz]<-model1$mu.coefSmo[[1]]$coef[as.character(data4$Site_ZZZ[iz])]
  }else{tem_rnd[iz]<-mean(model1$mu.coefSmo[[1]]$coef)}
}
}else{print('no random effects for this term')};


mu0<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))




if(!is.null(model1$sigma.coefSmo[[1]]$coef))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=unique(data1$Site_ZZZ))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}
data4 <- do.call( what=expand.grid, args=data3 )

data4[,'feature']<-0
tem_rnd<-matrix(0,dim(data4)[1],1);
Model.Frame<-model.frame(formula = model1$sigma.formula,data=data4);
Model.Matrix<-model.matrix(model1$sigma.formula,Model.Frame)
Fit.fix<-matrix(model1$sigma.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
if(!is.null(model1$sigma.coefSmo[[1]]))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data4)[1]){
  if(data4$Site_ZZZ[iz] %in% names(model1$sigma.coefSmo[[1]]$coef)){
    tem_rnd[iz]<-model1$sigma.coefSmo[[1]]$coef[as.character(data4$Site_ZZZ[iz])]
  }else{tem_rnd[iz]<-mean(model1$sigma.coefSmo[[1]]$coef)}
}
}else{print('no random effects for this term')};

sigma0<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))



if(!is.null(model1$nu.coefSmo[[1]]$coef))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=unique(data1$Site_ZZZ))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}

data4 <- do.call( what=expand.grid, args=data3 )

Model.Frame<-model.frame(formula = model1$nu.formula,data=data4);
Model.Matrix<-model.matrix(model1$nu.formula,Model.Frame)
Fit.fix<-matrix(model1$nu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
if(!is.null(model1$nu.coefSmo[[1]]))
{Fit.fix[length(Fit.fix)]=0;
}else{print('no random effects for this term')};

nu0<-as.vector(Model.Matrix[,colnames(Model.Matrix)[1]]*Fit.fix[1])



tem_par<-mu0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;

if(Seg>1)
{
  for(Seg1 in c(2:Seg))
  {
    par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
  }
  par=par/Seg
}
mu=par


tem_par<-sigma0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;
if(Seg>1)
{
  for(Seg1 in c(2:Seg))
  {
    par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
  }
  par=par/Seg
}
sigma=par


tem_par<-nu0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;
if(Seg>1)
{
  for(Seg1 in c(2:Seg))
  {
    par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
  }
  par=par/Seg
}
nu=par


p2_site<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
                 cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
                 calibration=FALSE,lpar=3)
p2_site[,'sigma']<-sigma


library(reshape2);
colnames(p2)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma');

colnames(p2_cal)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma');

colnames(p2_site)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma');


if(!(sheet %in% c('lh.aparc.thickness.table','rh.aparc.thickness.table')))
   {
     scale1=10000;
     ylab1='×10^4 mm3';
   }

if(sheet %in% c('lh.aparc.thickness.table','rh.aparc.thickness.table'))
   {
     scale1=1;
     ylab1='mm';
   }



# png(filename = paste0(str,'_',i,'_all_without_sex_stratified.png'), 
#     width = 1480,           
#     height = 740,          
#     units = "px",          
#     bg = "white",          
#     res = 300)     
# 
# p3<-ggplot()+
#   #geom_line(data=mydata,aes(x=Age,y=Value,group=Percentile,color=Percentile))+
#   geom_point(data=data1[data1$Sex=='Female',],aes(x=Age,y=tem_feature/scale1),
#              colour=c('grey'),shape=16,size=3,alpha = 0.1)+
#   geom_point(data=data1[data1$Sex=='Male',],aes(x=Age,y=tem_feature/scale1),
#              colour=c('grey'),shape=17,size=3,alpha = 0.1)+
#   geom_line(data=p2,aes(x=Age,y=median/scale1),color=c('#262626'),linewidth=2,linetype=c('solid'))+
#   geom_line(data=p2_cal,aes(x=Age,y=median/scale1),color=c('#4FBBD8'),linewidth=2,linetype=c('dashed'))+
#   geom_line(data=p2_site,aes(x=Age,y=median/scale1),color=c('#E84935'),linewidth=2,linetype=c('dotted'))+
#   # geom_line(data=p2,aes(x=Age,y=lower99CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dashed'))+
#   # geom_line(data=p2,aes(x=Age,y=lower95CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dotted'))+
#   # geom_line(data=p2,aes(x=Age,y=upper95CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dotted'))+
#   # geom_line(data=p2,aes(x=Age,y=upper99CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dashed'))+
#   labs(title=paste0(i,' ',ylab1),x='',y='')+
#   theme_bw()+
#   theme(
#     axis.title = element_text(family = "serif",size=12,color = "black"),
#     axis.text.x = element_text(
#       size = 12,              
#       color = "black",          
#       family = "serif"     
#     ),
#     axis.text.y = element_text(
#       size = 10,              
#       color = "black",         
#       #face = "bold" ,          
#       family = "serif"
#     )
#   )+
#   # scale_x_log10(breaks = c(6,18,35,80),  
#   #               labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
#   scale_x_continuous(breaks = c(6,18,35,80),  
#                      labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
# 
# print(p3)  
# dev.off()




}

}



