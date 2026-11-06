rm(list=ls())

datapath='/media/ps/Seagate Basic/Lifespan_freesurfer_results/Github/Source-codes'# change
setwd(datapath)
source("100.common-variables.r")
source("101.common-functions.r")

source("ZZZ_function.R")
source("300.variables.r")
source("301.functions.r")
library(readxl)
library(ggplot2)
library(gamlss)

clinical_datapath='/media/ps/Seagate Basic/Lifespan_freesurfer_results/HC_subtypes/data/new_final_list1_update_20260120_update_with_WMH.csv'

savepath='/media/ps/Seagate Basic/Lifespan_freesurfer_results/HC_subtypes/models_subtypes';

MR_datapath='/media/ps/Seagate Basic/Lifespan_freesurfer_results/HC_subtypes/data/MR_Results.xlsx'

var<-c('Global.table',
       'aseg.vol.table',
       'lh.aparc.volume.table','rh.aparc.volume.table',
       'lh.aparc.thickness.table','rh.aparc.thickness.table',
       'lh.aparc.area.table','rh.aparc.area.table')



library(stringr)

data_ref<-readRDS('/media/ps/Seagate Basic/Lifespan_freesurfer_results/HC_subtypes/results/HC_data_cluster_euclidean.rds')


for(cluster in c(1:4))
{

for(sheet in var) #for local feature
{ 
  
  
  setwd(datapath)
  MRI <- read_excel(MR_datapath,sheet=sheet)
  
  MRI<-MRI[MRI$Freesurfer_Path2!='Epilepsy_dicom_info_nii'&
             MRI$Freesurfer_Path2!='ET'&
             MRI$Freesurfer_Path2!='Guojibu_HC_nii',]
  MRI<-MRI[!is.na(MRI$Freesurfer_Path3),]
  MRI<-as.data.frame(MRI)
  
  MRI<- MRI[!duplicated(paste0(MRI$Freesurfer_Path2, MRI$Freesurfer_Path3)), ]
  rownames(MRI) <- paste0(MRI$Freesurfer_Path2, MRI$Freesurfer_Path3)
  
  rownames(MRI)<-paste0(MRI$Freesurfer_Path2,MRI$Freesurfer_Path3)
  if(str_detect(sheet,'aparc'))
  {
    tem_feature<-colnames(MRI)[c(2:35)];
  }
  
  if(str_detect(sheet,'DU15NET'))
  {
    tem_feature<-colnames(MRI)[c(3:18)];
  }
  
  if(str_detect(sheet,'glasser'))
  {
    tem_feature<-colnames(MRI)[c(3:181)];
  }
  
  
  if(str_detect(sheet,'BN_Atlas'))
  {
    tem_feature<-colnames(MRI)[c(2:106)];
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
    #tem_feature<-colnames(MRI)[c(5,13,23,24)];
    
    str='Global_feature';#str_lab='Global_feature';
    
  }
  
  str_cluster<-paste0(str,'_',cluster)
  if (!(dir.exists(paste0(savepath,'/',str_cluster))))
  {dir.create(paste0(savepath,'/',str_cluster))}
  
  setwd(paste0(savepath,'/',str_cluster))
  
  
  setwd(datapath)
  data1<-read.csv(clinical_datapath,header=TRUE);
  #data1$Age<-data1$Age_y
  data1$Site_ZZZ<-paste0(data1$Province,data1$Center,data1$Manufacturer)
  data1[,'Euler']<-data1$euler_number_l+data1$euler_number_r;
  data1[is.na(data1$Euler),'Euler']<--2;
  
  Euler_bh<-data1[,'Euler'];
  
  Euler_bh<-Euler_bh[!is.na(Euler_bh)]
  median_Euler<-median(Euler_bh);
  low_Euler<-median_Euler-2*sd(Euler_bh[Euler_bh!=-2])
  
  
  
  Z_data<-list();
  Quant_data<-list()
  
  
  for(i in tem_feature[1:length(tem_feature)])
  {
    
    print(i)
    setwd(paste0(savepath,'/',str_cluster))
    if(file.exists(paste0(str,'_',i,'_loop_our_model.rds'))){print('file exist');next;print('file exist')};
    
    
    # for each feature, we should load the oirginal clinical information
    setwd(datapath)  
    data1<-read.csv(clinical_datapath,header=TRUE);
    data1$Site_ZZZ<-paste0(data1$Province,data1$Center,data1$Manufacturer)
    
    library(dplyr)
    site_count <- data1 %>%
      group_by(Site_ZZZ) %>%  
      summarise(count = n())  
    site_count<-site_count[order(site_count$count),]
    print(site_count)
    
    for(site in unique(site_count$Site_ZZZ))
    {
      if(site_count[site_count$Site_ZZZ==site,'count']<10)
      {
        data1[data1$Site_ZZZ==site,'Database_included']<-0
      }
      
    }
    #data1$Site_ZZZ<-paste0(data1$Province)
    
    data1<-data1[!is.na(data1$Freesufer_Path2),]
    data1<- data1[!duplicated(paste0(data1$Freesufer_Path2, data1$Freesufer_Path3)), ]
    
    rownames(data1) <- paste0(data1$Freesufer_Path2, data1$Freesufer_Path3)
    
    setwd(paste0(savepath,'/',str_cluster))
    
  
    
    
    inter_row<-intersect(rownames(data1),rownames(MRI))
    data1=cbind(data1[inter_row,],MRI[inter_row,i])
    
    
    colnames(data1)[dim(data1)[2]]=c('tem_feature')
    
    data1[,'Euler']<-data1$euler_number_l+data1$euler_number_r;
    data1[is.na(data1$Euler),'Euler']<--2;
    
 
    all_data<-data1[data1$Database_included==1&
                      !is.na(data1$baseline),]
    
    
    all_data <- all_data[!duplicated(paste0(all_data$Freesufer_Path2, all_data$Freesufer_Path3)), ]
    rownames(all_data) <- paste0(all_data$Freesufer_Path2, all_data$Freesufer_Path3)
    
    all_data_original<-all_data
    
    data1<-all_data
    
    data1=data1[!is.na(data1$tem_feature)&!is.na(data1$Data_baseline)&data1$Diagnosis=='HC'&data1$Age>=4&data1$Age<=85,]
    
    data1$Site_ZZZ<-as.factor(data1$Site_ZZZ)
    
    data1$Sex<-as.factor(data1$Sex)
    
    data1$Sex<-factor(data1$Sex,levels=c('Female','Male'))
    
    
    data1<-data1[order(data1$Age),]
    data1[,'feature']<-data1$tem_feature
    all_data[,'feature']<-all_data$tem_feature
    
    #remove the extreme values
    data1<-data1[!is.na(data1$tem_feature),]
    data1<-data1[data1$feature>(mean(data1$feature)-3*sd(data1$feature))&
                   data1$feature<(mean(data1$feature)+3*sd(data1$feature)),]
    
    data1<-data1[data1$feature>0,]
    #select the columns
    data1<-data1[,c('Age','Sex','Site_ZZZ','tem_feature','feature')]
    
  
    
    data_ref_cluster<-data_ref[data_ref$Cluster==cluster,]
    ind_cluster<-intersect(rownames(data1),rownames(data_ref_cluster));
    data1<-data1[ind_cluster,]
    
    #step1 choose the best fit
    
    list_par<-data.frame(matrix(0,3*3*2*2,1));
    
    con=gamlss.control()
    
    num=0;
    
    results_try<-try({
      for(i_poly in 1:3)
      {
        for(j_poly in 1:3)
        {
          for(i_rnd in 0:1)
          {
            for(j_rnd in 0:1)
            {
              num=num+1;
              list_par[num,1]<-i_poly
              list_par[num,2]<-j_poly
              list_par[num,3]<-i_rnd
              list_par[num,4]<-j_rnd
              
              #list_fit=fit_model(i_poly,j_poly,i_rnd,j_rnd,num)
              
            }}}}
      library(doParallel)
      library(foreach)
      
      
      cl<-makeCluster(10)
      registerDoParallel(cl)
      my_data<-foreach(num=1:dim(list_par)[1],
                       .combine=rbind,
                       .packages = c('gamlss')) %dopar% fit_model(num)
      stopCluster(cl)
      
      
      # my_data<-data.frame(matrix(0,3*3*2*2,9));
      # for(num in 1:dim(list_par)[1])
      # {
      # 
      #         list_fit1=fit_model(num)
      #         my_data[num,]<-list_fit1
      # }
      
      list_fit<-my_data
      print(list_fit)
      
      
      #fit using the bestfit npoly and random with lowest BIC
      model_ind<-which.min(list_fit$BIC);
      sel_mu_poly=list_fit$mu_poly[model_ind]
      sel_sigma_poly=list_fit$sigma_poly[model_ind]
      i_rnd=list_fit$mu_random[model_ind]
      j_rnd=list_fit$sigma_random[model_ind]
    })
    
    
    
    if(inherits(results_try,'try-error')) 
    {sel_mu_poly=2
    sel_sigma_poly=2
    i_rnd=1
    j_rnd=1
    con=gamlss.control(c.crit = 0.01, n.cyc = 2,autostep = FALSE)
    }
    
    
    con=gamlss.control(c.crit = 0.01, n.cyc = 3,autostep = FALSE)
    #data1<-data1_backup;
   # con=gamlss.control(c.crit = 0.01, n.cyc = 1,autostep = FALSE)
    
    m0<-best_fit(sel_mu_poly,sel_sigma_poly,i_rnd,j_rnd)
    
    
    if(i_rnd==1&j_rnd==1){
      m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                 sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                 control=con,
                 family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                 data=data1)}else if(i_rnd==1&j_rnd==0){
                   m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                              sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex,
                              control=con,
                              family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                              data=data1)}else if(i_rnd==0&j_rnd==1){
                                m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex,
                                           sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                                           control=con,
                                           family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                           data=data1)}else if(i_rnd==0&j_rnd==0){
                                             m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex,
                                                        sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex,
                                                        control=con,
                                                        family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                                        data=data1)}
    
    
    #for all population plot
    if(i_rnd==1&j_rnd==1){
      m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+random(Site_ZZZ),
                 sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+random(Site_ZZZ),
                 control=con,
                 family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                 data=data1)}else if(i_rnd==1&j_rnd==0){
                   m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+random(Site_ZZZ),
                              sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power)),
                              control=con,
                              family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                              data=data1)}else if(i_rnd==0&j_rnd==1){
                                m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power)),
                                           sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+random(Site_ZZZ),
                                           control=con,
                                           family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                           data=data1)}else if(i_rnd==0&j_rnd==0){
                                             m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power)),
                                                        sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power)),
                                                        control=con,
                                                        family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                                        data=data1)}
    
    
    
    
  
    #all population#all population#all population#all population#all population#all population
    #both male and female
    model1<-m3;
    
    num_length=5000
    if(!is.null(model1$mu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    mu0 <- predict(model1, newdata = data4, type = "response", what = "mu")
    
    if(!is.null(model1$sigma.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    sigma0 <- predict(model1, newdata = data4, type = "response", what = "sigma")
    
    
    if(!is.null(model1$nu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
    }
    
    data4 <- do.call( what=expand.grid, args=data3 )
    
    nu0 <- predict(model1, newdata = data4, type = "response", what = "nu")
    
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
    
    
    
    
    
    model1<-m2;
    
    
    
    if(!is.null(model1$mu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    mu0 <- predict(model1, newdata = data4, type = "response", what = "mu")
    
    if(!is.null(model1$sigma.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    sigma0 <- predict(model1, newdata = data4, type = "response", what = "sigma")
    
    
    if(!is.null(model1$nu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
    }
    
    data4 <- do.call( what=expand.grid, args=data3 )
    
    nu0 <- predict(model1, newdata = data4, type = "response", what = "nu")
    
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
    
    
    p2_all<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
                     cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
                     calibration=FALSE,lpar=3)
    p2_all[,'sigma']<-sigma
    
    
    
    #male
    model1<-m2;
    if(!is.null(model1$mu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Male"),Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Male")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    mu0 <- predict(model1, newdata = data4, type = "response", what = "mu")
    
    if(!is.null(model1$sigma.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Male"),Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Male")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    sigma0 <- predict(model1, newdata = data4, type = "response", what = "sigma")
    
    
    if(!is.null(model1$nu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Male"),Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Male")) 
    }
    
    data4 <- do.call( what=expand.grid, args=data3 )
    
    nu0 <- predict(model1, newdata = data4, type = "response", what = "nu")
    
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
    
    
    male_p2<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
                      cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
                      calibration=FALSE,lpar=3)
    male_p2[,'sigma']<-sigma
    
    
    #female
    
    model1<-m2;
    if(!is.null(model1$mu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female"),Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    mu0 <- predict(model1, newdata = data4, type = "response", what = "mu")
    
    if(!is.null(model1$sigma.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female"),Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female")) 
    }
    data4 <- do.call( what=expand.grid, args=data3 )
    
    sigma0 <- predict(model1, newdata = data4, type = "response", what = "sigma")
    
    
    if(!is.null(model1$nu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female"),Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
    } else
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female")) 
    }
    
    data4 <- do.call( what=expand.grid, args=data3 )
    
    nu0 <- predict(model1, newdata = data4, type = "response", what = "nu")
    
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
    
    female_p2<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
                        cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
                        calibration=FALSE,lpar=3)
    female_p2[,'sigma']<-sigma
    
    
    
    library(reshape2);
    colnames(p2)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma');
    mydata<-melt(p2,id='Age');colnames(mydata)<-c('Age','Percentile','Value')
    
    step_age<-(max(data1$Age)-min(data1$Age))/num_length
    dim(p2)[1]-1
    Grad_p2<-(p2$median[2:dim(p2)[1]]-p2$median[1:(dim(p2)[1]-1)])/step_age
    Grad_p2<-data.frame(c(Grad_p2,Grad_p2[dim(p2)[1]-1]));
    p2<-cbind(p2,Grad_p2)
    colnames(p2)[dim(p2)[2]]<-c('Gradient1')
    
    if(!(str_detect(sheet,'thickness')))
    {
      scale1=10000;
      ylab1='×10^4 mm3';
    }
    
    if(str_detect(sheet,'thickness'))
    {
      scale1=1;
      ylab1='mm';
    }
    
    png(filename = paste0(str,'_',i,'_all_without_sex_stratified_Gradient.png'), 
        width = 1480,           
        height = 740,          
        units = "px",          
        bg = "white",          
        res = 300)     
    
    p3<-ggplot()+
      #geom_line(data=mydata,aes(x=Age,y=Value,group=Percentile,color=Percentile))+
      # geom_point(data=data1[data1$Sex=='Female',],aes(x=Age,y=tem_feature),
      #            colour=c('#990000'),shape=16,size=3,alpha = 0.1)+
      # geom_point(data=data1[data1$Sex=='Male',],aes(x=Age,y=tem_feature),
      #            colour=c('#00CCFF'),shape=17,size=3,alpha = 0.1)+
      geom_line(data=p2,aes(x=Age,y=Gradient1/scale1),color=c('#262626'),linewidth=1,linetype=c('solid'))+
      labs(title=paste0(i,' ',ylab1),x='',y='')+
      theme_bw()+
      theme(
        axis.title = element_text(family = "serif",size=12,color = "black"),
        axis.text.x = element_text(
          size = 12,              
          color = "black",          
          family = "serif"      
        ),
        axis.text.y = element_text(
          size = 10,              
          color = "black",        
          #face = "bold" ,         
          family = "serif"
        )
      )+
      # scale_x_log10(breaks = c(6,18,35,80),  
      #               labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
      scale_x_continuous(breaks = c(6,18,35,80),  
                         labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
    
    
    print(p3)  
    dev.off()
    
    
    
    png(filename = paste0(str,'_',i,'_all_without_sex_stratified.png'), 
        width = 1480,           
        height = 740,          
        units = "px",          
        bg = "white",          
        res = 300)     
    
    p3<-ggplot()+
      #geom_line(data=mydata,aes(x=Age,y=Value,group=Percentile,color=Percentile))+
      geom_point(data=data1[data1$Sex=='Female',],aes(x=Age,y=tem_feature/scale1),
                 colour=c('#E84935'),shape=16,size=3,alpha = 0.1)+
      geom_point(data=data1[data1$Sex=='Male',],aes(x=Age,y=tem_feature/scale1),
                 colour=c('#4FBBD8'),shape=17,size=3,alpha = 0.1)+
      geom_line(data=p2,aes(x=Age,y=median/scale1),color=c('#262626'),linewidth=1,linetype=c('solid'))+
      geom_line(data=p2,aes(x=Age,y=lower99CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dashed'))+
      geom_line(data=p2,aes(x=Age,y=lower95CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dotted'))+
      geom_line(data=p2,aes(x=Age,y=upper95CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dotted'))+
      geom_line(data=p2,aes(x=Age,y=upper99CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dashed'))+
      labs(title=paste0(i,' ',ylab1),x='',y='')+
      theme_bw()+
      theme(
        axis.title = element_text(family = "serif",size=12,color = "black"),
        axis.text.x = element_text(
          size = 12,              
          color = "black",          
          family = "serif"     
        ),
        axis.text.y = element_text(
          size = 10,              
          color = "black",         
          #face = "bold" ,          
          family = "serif"
        )
      )+
      # scale_x_log10(breaks = c(6,18,35,80),  
      #               labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
      scale_x_continuous(breaks = c(6,18,35,80),  
                         labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
    
    print(p3)  
    dev.off()
    
    
    #for all population with sex stratified
    #female data
    colnames(female_p2)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma');
    mydata<-melt(female_p2,id='Age');colnames(mydata)<-c('Age','Percentile','Value')
    Female_p2<-female_p2;
    
    
    #male data
    colnames(male_p2)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma');
    mydata<-melt(male_p2,id='Age');colnames(mydata)<-c('Age','Percentile','Value')
    Male_p2<-male_p2;
    
    
    dim(Female_p2)[1]-1
    Grad_Female_p2<-(Female_p2$median[2:dim(Female_p2)[1]]-Female_p2$median[1:(dim(Female_p2)[1]-1)])/step_age
    Grad_Female_p2<-data.frame(c(Grad_Female_p2,Grad_Female_p2[dim(Female_p2)[1]-1]));
    Female_p2<-cbind(Female_p2,Grad_Female_p2)
    colnames(Female_p2)[dim(Female_p2)[2]]<-c('Gradient1')
    
    dim(Male_p2)[1]-1
    Grad_Male_p2<-(Male_p2$median[2:dim(Male_p2)[1]]-Male_p2$median[1:(dim(Male_p2)[1]-1)])/step_age
    Grad_Male_p2<-data.frame(c(Grad_Male_p2,Grad_Male_p2[dim(Male_p2)[1]-1]));
    Male_p2<-cbind(Male_p2,Grad_Male_p2)
    colnames(Male_p2)[dim(Male_p2)[2]]<-c('Gradient1')
    
    
    png(filename = paste0(str,'_',i,'_all_with_sex_stratified_Gradient.png'), 
        width = 1480,           
        height = 740,          
        units = "px",          
        bg = "white",          
        res = 300)     
    
    p3<-ggplot()+
      #geom_line(data=mydata,aes(x=Age,y=Value,group=Percentile,color=Percentile))+
      # geom_point(data=data1[data1$Sex=='Female',],aes(x=Age,y=tem_feature),
      #            colour=c('#990000'),shape=16,size=3,alpha = 0.1)+
      # geom_point(data=data1[data1$Sex=='Male',],aes(x=Age,y=tem_feature),
      #            colour=c('#00CCFF'),shape=17,size=3,alpha = 0.1)+
      geom_line(data=Female_p2,aes(x=Age,y=Gradient1/scale1),color=c('#E84935'),linewidth=1,linetype=c('solid'))+
      geom_line(data=Male_p2,aes(x=Age,y=Gradient1/scale1),color=c('#4FBBD8'),linewidth=1,linetype=c('solid'))+
      labs(title=paste0(i,' ',ylab1),x='',y='')+
      theme_bw()+
      theme(
        axis.title = element_text(family = "serif",size=12,color = "black"),
        axis.text.x = element_text(
          size = 12,              
          color = "black",          
          family = "serif"      
        ),
        axis.text.y = element_text(
          size = 10,              
          color = "black",         
          #face = "bold" ,         
          family = "serif"
        )
      )+
      # scale_x_log10(breaks = c(6,18,35,80),  
      #               labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
      scale_x_continuous(breaks = c(6,18,35,80),  
                         labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
    
    print(p3)  
    dev.off()
    
    
    png(filename = paste0(str,'_',i,'_all_with_sex_stratified_sigma.png'), 
        width = 1480,           
        height = 740,          
        units = "px",          
        bg = "white",          
        res = 300)     
    
    p3<-ggplot()+
      #geom_line(data=mydata,aes(x=Age,y=Value,group=Percentile,color=Percentile))+
      # geom_point(data=data1[data1$Sex=='Female',],aes(x=Age,y=tem_feature),
      #            colour=c('#990000'),shape=16,size=3,alpha = 0.1)+
      # geom_point(data=data1[data1$Sex=='Male',],aes(x=Age,y=tem_feature),
      #            colour=c('#00CCFF'),shape=17,size=3,alpha = 0.1)+
      geom_line(data=Female_p2,aes(x=Age,y=sigma),color=c('#E84935'),linewidth=1,linetype=c('solid'))+
      geom_line(data=Male_p2,aes(x=Age,y=sigma),color=c('#4FBBD8'),linewidth=1,linetype=c('solid'))+
      labs(title=paste0(i,' ',ylab1),x='',y='')+
      theme_bw()+
      theme(
        axis.title = element_text(family = "serif",size=12,color = "black"),
        axis.text.x = element_text(
          size = 12,              
          color = "black",         
          family = "serif"      
        ),
        axis.text.y = element_text(
          size = 10,              
          color = "black",         
          #face = "bold" ,          
          family = "serif"
        )
      )+
      # scale_x_log10(breaks = c(6,18,35,80),  
      #               labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
      scale_x_continuous(breaks = c(6,18,35,80),  
                         labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
    
    print(p3)  
    dev.off()
    
    
    
    
    png(filename = paste0(str,'_',i,'_all_with_sex_stratified.png'), 
        width = 1480,           
        height = 740,          
        units = "px",          
        bg = "white",          
        res = 300)     
    
    p3<-ggplot()+
      #geom_line(data=mydata,aes(x=Age,y=Value,group=Percentile,color=Percentile))+
      geom_point(data=data1[data1$Sex=='Female',],aes(x=Age,y=tem_feature/scale1),
                 colour=c('#E84935'),shape=16,size=3,alpha = 0.1)+
      geom_point(data=data1[data1$Sex=='Male',],aes(x=Age,y=tem_feature/scale1),
                 colour=c('#4FBBD8'),shape=17,size=3,alpha = 0.1)+
      geom_line(data=Female_p2,aes(x=Age,y=median/scale1),color=c('#E84935'),linewidth=1,linetype=c('solid'))+
      #geom_line(data=Female_p2,aes(x=Age,y=lower99CI),color=c('#990000'),linewidth=2,linetype=c('dashed'))+
      geom_line(data=Female_p2,aes(x=Age,y=lower95CI/scale1),color=c('#E84935'),linewidth=1,linetype=c('dotted'))+
      geom_line(data=Female_p2,aes(x=Age,y=upper95CI/scale1),color=c('#E84935'),linewidth=1,linetype=c('dotted'))+
      #geom_line(data=Female_p2,aes(x=Age,y=upper99CI),color=c('#990000'),linewidth=2,linetype=c('dashed'))+
      
      geom_line(data=Male_p2,aes(x=Age,y=median/scale1),color=c('#4FBBD8'),linewidth=1,linetype=c('solid'))+
      #geom_line(data=Male_p2,aes(x=Age,y=lower99CI),color=c('#00CCFF'),linewidth=2,linetype=c('dashed'))+
      geom_line(data=Male_p2,aes(x=Age,y=lower95CI/scale1),color=c('#4FBBD8'),linewidth=1,linetype=c('dotted'))+
      geom_line(data=Male_p2,aes(x=Age,y=upper95CI/scale1),color=c('#4FBBD8'),linewidth=1,linetype=c('dotted'))+
      #geom_line(data=Male_p2,aes(x=Age,y=upper99CI),color=c('#00CCFF'),linewidth=2,linetype=c('dashed'))+
      
      # geom_line(data=p2,aes(x=Age,y=median),color=c('#666666'),linewidth=1,linetype=c('solid'))+
      # #geom_line(data=Male_p2,aes(x=Age,y=lower99CI),color=c('#00CCFF'),linewidth=2,linetype=c('dashed'))+
      # geom_line(data=p2,aes(x=Age,y=lower95CI),color=c('#666666'),linewidth=1,linetype=c('dotted'))+
      # geom_line(data=p2,aes(x=Age,y=upper95CI),color=c('#666666'),linewidth=1,linetype=c('dotted'))+
      # #geom_line(data=Male_p2,aes(x=Age,y=upper99CI),color=c('#00CCFF'),linewidth=2,linetype=c('dashed'))+
      
      labs(title=paste0(i,' ',ylab1),x='',y='')+
      theme_bw()+
      theme(
        axis.title = element_text(family = "serif",size=12,color = "black"),
        axis.text.x = element_text(
          size = 12,              
          color = "black",          
          family = "serif"      
        ),
        axis.text.y = element_text(
          size = 10,             
          color = "black",         
          #face = "bold" ,          
          family = "serif"
        )
      )+
      # scale_x_log10(breaks = c(6,18,35,80),  
      #               labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
      scale_x_continuous(breaks = c(6,18,35,80),  
                         labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
    
    print(p3)  
    
    dev.off()
    
    #calculate the quantile for all cases including HC and Diseases
    
    Z_score_sum<-NULL;
    Quant_score_sum<-NULL
    
    #all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%
    #all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%
    #all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%
    #all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%#all_data%%
    
    all_data<-all_data[all_data$tem_feature!=''&!is.null(all_data$tem_feature)&!is.na(all_data$tem_feature)&!is.infinite(all_data$tem_feature),]
    all_data<-all_data[all_data$Age!=''&!is.null(all_data$Age)&!is.na(all_data$Age)&!is.infinite(all_data$Age),]
    all_data<-all_data[all_data$Site_ZZZ!=''&!is.null(all_data$Site_ZZZ)&!is.na(all_data$Site_ZZZ)&!is.infinite(all_data$Site_ZZZ),]
    all_data<-all_data[all_data$Sex!=''&!is.null(all_data$Sex)&!is.na(all_data$Sex)&!is.infinite(all_data$Sex),]
    
    all_data1<-all_data[,c('Age','Sex','Site_ZZZ','tem_feature')]
    
    all_data1$Sex<-as.factor(all_data1$Sex)
    all_data1$Site_ZZZ<-as.factor(all_data1$Site_ZZZ)
    
    
    model1<-m2;
    
    
    for(sub in 1:dim(all_data)[1])
    {
      
      if(!is.null(m2$mu.coefSmo[[1]]))
      {
        if(!(all_data1$Site_ZZZ[sub] %in% names(m2$mu.coefSmo[[1]]$coef)))
        {
          #all_data1$Site_ZZZ[sub]<-c('BeijingTiantanGE')
          all_data1$Site_ZZZ[sub]<-names(which.max(abs(m2$mu.coefSmo[[1]]$coef-mean(m2$mu.coefSmo[[1]]$coef))))
          
          #model1$mu.coefSmo[[1]]$coef[sub]=median(m2$mu.coefSmo[[1]]$coef)
        }
      }
      
      if(!is.null(m2$sigma.coefSmo[[1]]))
      {
        if(!(all_data1$Site_ZZZ[sub] %in% names(m2$sigma.coefSmo[[1]]$coef)))
        {
          #all_data1$Site_ZZZ[sub]<-c('BeijingTiantanGE')
          all_data1$Site_ZZZ[sub]<-names(which.max(abs(m2$sigma.coefSmo[[1]]$coef-mean(m2$sigma.coefSmo[[1]]$coef))))
          #model1$sigma.coefSmo[[1]]$coef[sub]=median(m2$sigma.coefSmo[[1]]$coef)
        }
      }
      
    }
    
    
    mu <- predict(model1, newdata = all_data1, type = "response", what = "mu")
    sigma <- predict(model1, newdata = all_data1, type = "response", what = "sigma")
    nu <- predict(model1, newdata = all_data1, type = "response", what = "nu")
    
    
    
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
    
    
    results<-list();
    results$Female_p2<-Female_p2
    results$Male_p2<-Male_p2
    results$p2<-p2
    results$peakage<-p2$Age[which.max(p2$median)]
    results$p2_all<-p2_all
    results$m2<-m2
    results$m0<-m0
    results$m3<-m3
    
    #results$list_fit<-list_fit
    
    results$Zscore<-Z_data
    results$Quant_data<-Quant_data
    results$data1<-data1
    results$all_data<-all_data1
    #results$all_data1<-all_data1
    results$str<-str
    results$i<-i
    results$all_data_original<-all_data_original
    
    
    saveRDS(results,paste0(str,'_',i,'_loop_our_model.rds'))
    
    
    # #five-fold cross-validation#five-fold cross-validation#five-fold cross-validation#five-fold cross-validation
    # #five-fold cross-validation#five-fold cross-validation#five-fold cross-validation#five-fold cross-validation
    # #five-fold cross-validation#five-fold cross-validation#five-fold cross-validation#five-fold cross-validation
    # #five-fold cross-validation#five-fold cross-validation#five-fold cross-validation#five-fold cross-validation
    # 
    # set.seed(123)
    # library('caret')
    # #data1<-data1[,c("Age","Sex","Site_ZZZ","feature"),]
    # data1$Site_ZZZ<-as.factor(data1$Site_ZZZ)
    # #folds<-createFolds(y=data1$Site_ZZZ,k=5,list=TRUE)
    # 
    # #sample by site
    # create_center_folds <- function(data, centers, k) {
    #   folds <- vector("list", k)
    #   for (center in unique(data[[centers]])) {
    #     center_data <- data[data[[centers]] == center, ]
    #     center_folds <- createFolds(center_data$Site_ZZZ, k = k, list = TRUE, returnTrain = FALSE)
    #     for (i in 1:k) {
    #       folds[[i]] <- c(folds[[i]], center_folds[[i]])
    #     }
    #   }
    #   return(folds)
    # }
    # 
    # k <- 5  
    # folds <- create_center_folds(data1, "Site_ZZZ", k)
    # 
    # 
    # Z_score_folds_HC1<-NULL
    # Quant_score_folds_HC1<-NULL
    # #con=gamlss.control(c.crit = 0.01, n.cyc = 5,autostep = FALSE)
    # 
    # for(i_fold in 1:5)
    # {
    #   train_data<-data1[-folds[[i_fold]],]
    #   test_data<-data1[folds[[i_fold]],]
    # 
    # if(i_rnd==1&j_rnd==1){
    #   m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
    #              sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
    #              control=con,
    #              family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
    #              data=train_data)}else if(i_rnd==1&j_rnd==0){
    #                m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
    #                           sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex,
    #                           control=con,
    #                           family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
    #                           data=train_data)}else if(i_rnd==0&j_rnd==1){
    #                             m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex,
    #                                        sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
    #                                        control=con,
    #                                        family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
    #                                        data=train_data)}else if(i_rnd==0&j_rnd==0){
    #                                          m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex,
    #                                                     sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex,
    #                                                     control=con,
    #                                                     family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
    #                                                     data=train_data)}
    # 
    # 
    # model1<-m2;
    # 
    # mu=NULL
    # sigma=NULL
    # nu=NULL
    # 
    # 
    # for(site in unique(test_data$Site_ZZZ))
    # {
    # 
    #   if(!is.null(m2$mu.coefSmo[[1]]))
    #   {
    #  
    #     for(sub in dim(test_data)[1])
    #     {
    #     if(!(test_data$Site_ZZZ[sub] %in% names(m2$mu.coefSmo[[1]]$coef)))
    #     {
    #       
    #       test_data$Site_ZZZ[sub]<-names(which.max(abs(m2$mu.coefSmo[[1]]$coef-mean(m2$mu.coefSmo[[1]]$coef))))
    #     } 
    #     }
    #     
    #   }
    # 
    #   if(!is.null(m2$sigma.coefSmo[[1]]))
    #   {
    #  
    #     for(sub in dim(test_data)[1])
    #     {
    #     if(!(test_data$Site_ZZZ[sub] %in% names(m2$sigma.coefSmo[[1]]$coef)))
    #     {
    #       test_data$Site_ZZZ[sub]<-names(which.max(abs(m2$sigma.coefSmo[[1]]$coef-mean(m2$sigma.coefSmo[[1]]$coef))))
    #     }  
    #     }
    #     
    #   }
    # }
    # 
    #   mu <- predict(model1, newdata = test_data, type = "response", what = "mu")
    #   sigma <- predict(model1, newdata = test_data, type = "response", what = "sigma")
    #   nu <- predict(model1, newdata = test_data, type = "response", what = "nu")
    # 
    # 
    # if(length(mu)!=dim(test_data)[1])
    # {
    #   print("Error, Please Check Data!!!")
    # }
    # 
    # 
    #   Z_score_folds_HC<-zzz_cent(obj=model1,type=c("z-scores"),mu=mu,sigma=sigma,nu=nu,
    #                         xname = 'Age',xvalues=test_data$Age,yval=test_data$feature,
    #                         calibration=FALSE,lpar=3)
    # 
    #   Z_score_folds_HC<-data.frame(Z_score_folds_HC);
    # 
    #   rownames(Z_score_folds_HC)<-data1$individual_ID[folds[[i_fold]]]
    # 
    #   Z_score_folds_HC1<-rbind(Z_score_folds_HC1,Z_score_folds_HC)
    # 
    # 
    #   Quant_score_folds_HC<-zzz_cent(obj=model1,type=c("z-scores"),mu=mu,sigma=sigma,nu=nu,
    #                             xname = 'Age',xvalues=test_data$Age,yval=test_data$feature,
    #                             calibration=FALSE,lpar=3,cdf=TRUE)
    #   Quant_score_folds_HC<-data.frame(Quant_score_folds_HC);
    # 
    #   rownames(Quant_score_folds_HC)<-data1$individual_ID[folds[[i_fold]]]
    #   Quant_score_folds_HC1<-rbind(Quant_score_folds_HC1,Quant_score_folds_HC)
    # 
    # }
    # 
    # #five-fold cross-validation#five-fold cross-validation#five-fold cross-validation#five-fold cross-validation
    # 
    # results$Z_score_folds_HC<-Z_score_folds_HC1;
    # results$Quant_score_folds_HC<-Quant_score_folds_HC1;
    # 
    # # results$Zscore[[i]][rownames(Z_score_folds_HC1),"Z_score"]<-Z_score_folds_HC1;
    # # results$Quant_data[[i]][rownames(Quant_score_folds_HC1),"Quant_score"]<-Quant_score_folds_HC1;
    # 
    # saveRDS(results,paste0(str,'_',i,'_loop_our_model.rds'))
    # 
    # 
    # 
    # 
    # #resample HCs for the robustness of fitted curves and peak ages
    # 
    # sample_num=1000
    # 
    # con=gamlss.control(c.crit = 0.01, n.cyc = 5,autostep = FALSE)
    # p2_sample=results$p2[,c('Age','median','sigma')];
    # p2_sample[,'Index']<-'Original'
    # peak_age_sample=data.frame(results$peakage);
    # colnames(peak_age_sample)<-'peakage';
    # peak_age_sample[,'Index']<-'Original'
    # 
    # 
    # for (i_sample in 1:sample_num)
    # {
    # 
    # results_try<-try({
    # print(i_sample)
    # data_resample<-data1[sample(nrow(data1), size=nrow(data1), replace = TRUE),]
    # 
    # #for all population plot
    # if(i_rnd==1&j_rnd==1){
    #   m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+random(Site_ZZZ),
    #              sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+random(Site_ZZZ),
    #              control=con,
    #              family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
    #              data=data_resample)}else if(i_rnd==1&j_rnd==0){
    #                m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+random(Site_ZZZ),
    #                           sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power)),
    #                           control=con,
    #                           family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
    #                           data=data_resample)}else if(i_rnd==0&j_rnd==1){
    #                             m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power)),
    #                                        sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+random(Site_ZZZ),
    #                                        control=con,
    #                                        family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
    #                                        data=data_resample)}else if(i_rnd==0&j_rnd==0){
    #                                          m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power)),
    #                                                     sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power)),
    #                                                     control=con,
    #                                                     family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
    #                                                     data=data_resample)}
    # 
    # 
    # 
    # 
    # model1<-m3;
    # 
    # num_length=5000
    # if(!is.null(model1$mu.coefSmo[[1]]$coef))
    # {
    #   data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
    # } else
    # {
    #   data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"))
    # }
    # data4 <- do.call( what=expand.grid, args=data3 )
    # 
    # mu0 <- predict(model1, newdata = data4, type = "response", what = "mu")
    # 
    # if(!is.null(model1$sigma.coefSmo[[1]]$coef))
    # {
    #   data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
    # } else
    # {
    #   data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"))
    # }
    # data4 <- do.call( what=expand.grid, args=data3 )
    # 
    # sigma0 <- predict(model1, newdata = data4, type = "response", what = "sigma")
    # 
    # 
    # if(!is.null(model1$nu.coefSmo[[1]]$coef))
    # {
    #   data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
    # } else
    # {
    #   data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"))
    # }
    # 
    # data4 <- do.call( what=expand.grid, args=data3 )
    # 
    # nu0 <- predict(model1, newdata = data4, type = "response", what = "nu")
    # 
    # tem_par<-mu0
    # par<-tem_par[1:num_length]
    # Seg=length(tem_par)/num_length;
    # 
    # if(Seg>1)
    # {
    #   for(Seg1 in c(2:Seg))
    #   {
    #     par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
    #   }
    #   par=par/Seg
    # }
    # mu=par
    # 
    # 
    # tem_par<-sigma0
    # par<-tem_par[1:num_length]
    # Seg=length(tem_par)/num_length;
    # if(Seg>1)
    # {
    #   for(Seg1 in c(2:Seg))
    #   {
    #     par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
    #   }
    #   par=par/Seg
    # }
    # sigma=par
    # 
    # 
    # tem_par<-nu0
    # par<-tem_par[1:num_length]
    # Seg=length(tem_par)/num_length;
    # if(Seg>1)
    # {
    #   for(Seg1 in c(2:Seg))
    #   {
    #     par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
    #   }
    #   par=par/Seg
    # }
    # nu=par
    # 
    # 
    # p2<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
    #              cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
    #              calibration=FALSE,lpar=3)
    # p2[,'sigma']<-sigma
    # 
    # 
    # colnames(p2)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma');
    # p2<-p2[,c('Age','median','sigma')]
    # p2[,'Index']<-paste0('resample',i_sample)
    # p2_sample<-rbind(p2_sample,p2);
    # 
    # 
    # peak_age=data.frame(p2$Age[which.max(p2$median)]);
    # colnames(peak_age)<-'peakage';
    # peak_age[,'Index']<-paste0('resample',i_sample)
    # peak_age_sample<-rbind(peak_age_sample,peak_age)
    # })
    # 
    # 
    # if(inherits(results_try,'try-error'))
    # {next}
    # 
    # }
    # 
    # 
    # lower_peak_age =  quantile(peak_age_sample$peakage, 0.025)
    # upper_peak_age =  quantile(peak_age_sample$peakage, 0.975)
    # median_peak_age = median(peak_age_sample$peakage)
    # 
    # results$p2_sample<-p2_sample
    # results$peak_age_sample<-peak_age_sample
    # 
    # results$CI<-c(median_peak_age,lower_peak_age,upper_peak_age)
    # 
    # saveRDS(results,paste0(str,'_',i,'_loop_our_model.rds'))
    # 
    # 
    # 
    # #calculate the CI of the fitted curves
    # p2_sample_CI<-NULL
    # 
    # for(i_sample in unique(p2_sample$Index))
    # {
    # 
    #   p2_sample_CI<-cbind(p2_sample_CI,p2_sample[p2_sample$Index==i_sample,'median'])
    # }
    # 
    # 
    # ci_data <- data.frame(
    #   Age = p2_sample$Age[p2_sample$Index=="Original"],
    #   lower = apply(p2_sample_CI, 1, function(x) quantile(x, 0.025)),
    #   upper = apply(p2_sample_CI, 1, function(x) quantile(x, 0.975)),
    #   median = p2_sample[p2_sample$Index=="Original",'median']
    # )
    # 
    # 
    # png(filename = paste0(str,'_',i,'_all_without_sex_stratified_resample_CI.png'),
    #     width = 1480,
    #     height = 740,
    #     units = "px",
    #     bg = "white",
    #     res = 300)
    # 
    # p3<-ggplot()+
    #   #geom_line(data=mydata,aes(x=Age,y=Value,group=Percentile,color=Percentile))+
    #   # geom_point(data=data1[data1$Sex=='Female',],aes(x=Age,y=tem_feature/scale1),
    #   #            colour=c('#E84935'),shape=16,size=3,alpha = 0.1)+
    #   # geom_point(data=data1[data1$Sex=='Male',],aes(x=Age,y=tem_feature/scale1),
    #   #            colour=c('#4FBBD8'),shape=17,size=3,alpha = 0.1)+
    # 
    #   geom_line(data=ci_data,aes(x=Age,y=median/scale1),color='black',linewidth=1,linetype=c('solid'))+
    #   geom_line(data=ci_data,aes(x=Age,y=lower/scale1),color='black',linewidth=0.5,linetype=c('dashed'))+
    #   geom_line(data=ci_data,aes(x=Age,y=upper/scale1),color='black',linewidth=0.5,linetype=c('dashed'))+
    #   labs(title=paste0(i,' ',ylab1),x='',y='')+
    #   theme_bw()+
    #   theme(
    #     legend.position = 'none',
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
    #   scale_x_log10(breaks = c(6,18,35,80),
    #                 labels = c("6 yr", "18 yr", "35 yr", "80 yr"))
    # 
    # print(p3)
    # dev.off()
    
    
  }
  
}

}

