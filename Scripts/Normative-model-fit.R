rm(list=ls())
#set your own directory
datapath='E:/Lifespan_freesurfer_results/Github/Source-codes/'  # Change the directory where you save the Source-codes  
clinical_datapath='E:/Lifespan_freesurfer_results/Github/Datasets/Dataset-norms/Clinical_vars.csv'  # Change the directory where you save the clinical variables  
MR_datapath='E:/Lifespan_freesurfer_results/Github/Datasets/Dataset-norms/MR_measures.xlsx'  # Change the directory where you save the MR measures  
savepath='E:/Lifespan_freesurfer_results/Github/Test_results/V7_DK' # Create and determine the directory where you would save the results  

#load source functions
setwd(datapath)
source("100.common-variables.r")
source("101.common-functions.r")
source("300.variables.r")
source("301.functions.r")
source("ZZZ_function.R")

library(gamlss)
library(readxl)
library(stringr)
library(ggplot2)
library(ggsci)

#select the sheetnames of the MR meaures you want to fit
var<-c('Global.table',
        'aseg.vol.table',
        'lh.aparc.volume.table','rh.aparc.volume.table',
        'lh.aparc.thickness.table','rh.aparc.thickness.table',
        'lh.aparc.area.table','rh.aparc.area.table')


#fit the normative references of global and local measures one-by-one 
for(sheet in var) 
{ 
  
 
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
    for(col in colnames(MRI)[2:13])
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
    
    tem_feature<-colnames(MRI)[c(2,3,4,5,13,17,18)];
     str='Global_feature';    
  }
  
  
  if (!(dir.exists(paste0(savepath,'/',str))))
  {dir.create(paste0(savepath,'/',str))}
  
  setwd(paste0(savepath,'/',str))
  
  
  data1<-read.csv(clinical_datapath,header=TRUE);
  data1$Site_ZZZ<-paste0(data1$Province,data1$Center,data1$Manufacturer)
  data1[,'Euler']<-data1$euler_number_l+data1$euler_number_r;
  Euler_bh<-data1[,'Euler'];
  median_Euler<-median(Euler_bh);
  low_Euler<-median_Euler-2*sd(Euler_bh)
  
  
  Z_data<-list();
  Quant_data<-list()
  
  for(i in tem_feature[1:length(tem_feature)]) ##fit the normative references of MR measures one-by-one 
  {
    
    print(i)
    setwd(paste0(savepath,'/',str))
    if(file.exists(paste0(str,'_',i,'_loop_our_model.rds'))){print('file exist');next;print('file exist')};
    
    
    # for each feature, we should load the oirginal clinical information

    data1<-read.csv(clinical_datapath,header=TRUE);
    data1$Site_ZZZ<-paste0(data1$Province,data1$Center,data1$Manufacturer)
           
    rownames(data1)<-paste0(data1$Freesufer_Path2,data1$Freesufer_Path3)
    
    
    setwd(paste0(savepath,'/',str))
    
    inter_row<-intersect(rownames(data1),rownames(MRI))
    data1=cbind(data1[inter_row,],MRI[inter_row,i])
    
    colnames(data1)[dim(data1)[2]]=c('tem_feature')
    all_data<-data1
    data1[,'Euler']<-data1$euler_number_l+data1$euler_number_r;
            
    rownames(all_data)<-paste0(all_data$Freesufer_Path2,all_data$Freesufer_Path3)
    all_data_original<-all_data
    
    data1<-all_data
    
    # data1=data1[!is.na(data1$tem_feature)&&data1$Euler>low_Euler&
    #                        data1$Age>=4&data1$Age<=85,]
    
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
    
    data1<-data1[,c('Age','Sex','Site_ZZZ','tem_feature','feature')]
    
    data1_backup<-data1;
    

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
                            
            }}}}
      library(doParallel)
      library(foreach)
      
      
      cl<-makeCluster(10)
      registerDoParallel(cl)
      my_data<-foreach(num=1:dim(list_par)[1],
                       .combine=rbind,
                       .packages = c('gamlss')) %dopar% fit_model(num)
      stopCluster(cl)
      
         
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
    }
    
    
    
    data1<-data1_backup;
    m0<-best_fit(sel_mu_poly,sel_sigma_poly,i_rnd,j_rnd)
    
   #for all population 
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
    
     
   
    #all population for both male and female (figure plot)
    num_length=5000 
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
    
    
    p2<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
                     cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
                     calibration=FALSE,lpar=3)
    p2[,'sigma']<-sigma
    
    
    
    #males
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
    
    
    #females
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
    



   #save main results
    results<-list();
    results$Female_p2<-Female_p2
    results$Male_p2<-Male_p2
    results$p2<-p2
    results$peakage<-p2$Age[which.max(p2$median)]
    results$m2<-m2
    results$m0<-m0
    results$list_fit<-list_fit
    
    results$data1<-data1
    results$str<-str
    results$i<-i
    results$all_data_original<-all_data_original
    
    
    saveRDS(results,paste0(str,'_',i,'_loop_our_model.rds'))
    
    
    #10-fold cross-validation  
    set.seed(123)
    library('caret')

    # Siet-stratified sample
    create_center_folds <- function(data, centers, k) {
    folds <- vector("list", k)
    for (center in unique(data[[centers]])) {
         center_data <- data[data[[centers]] == center, ]
         center_folds <- createFolds(center_data$Site_ZZZ, k = k, list = TRUE, returnTrain = FALSE)
         for (i in 1:k) {
           folds[[i]] <- c(folds[[i]], center_folds[[i]])
         }
       }
       return(folds)
     }
     
    k <- 10  
    folds <- create_center_folds(data1, "Site_ZZZ", k)
      
   Z_score_folds_HC1<-NULL
   Quant_score_folds_HC1<-NULL
    for(i_fold in 1:k)
     {
       train_data<-data1[-folds[[i_fold]],]
       test_data<-data1[folds[[i_fold]],]
     
     if(i_rnd==1&j_rnd==1){
       m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                  sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                  control=con,
                  family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                  data=train_data)}else if(i_rnd==1&j_rnd==0){
                    m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                               sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex,
                               control=con,
                               family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                               data=train_data)}else if(i_rnd==0&j_rnd==1){
                                 m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex,
                                            sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                                            control=con,
                                            family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                            data=train_data)}else if(i_rnd==0&j_rnd==0){
                                              m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex,
                                                         sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex,
                                                         control=con,
                                                         family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                                         data=train_data)}
     
     
     model1<-m2;
     mu=NULL
     sigma=NULL
     nu=NULL 
    
    for(site in unique(test_data$Site_ZZZ))
     {
     
       if(!is.null(m2$mu.coefSmo[[1]]))
       {
      
         for(sub in dim(test_data)[1])
         {
         if(!(test_data$Site_ZZZ[sub] %in% names(m2$mu.coefSmo[[1]]$coef)))
         {
           
           test_data$Site_ZZZ[sub]<-names(which.max(abs(m2$mu.coefSmo[[1]]$coef-mean(m2$mu.coefSmo[[1]]$coef))))
         } 
         }
         
       }
     
       if(!is.null(m2$sigma.coefSmo[[1]]))
       {
      
         for(sub in dim(test_data)[1])
         {
         if(!(test_data$Site_ZZZ[sub] %in% names(m2$sigma.coefSmo[[1]]$coef)))
         {
           test_data$Site_ZZZ[sub]<-names(which.max(abs(m2$sigma.coefSmo[[1]]$coef-mean(m2$sigma.coefSmo[[1]]$coef))))
         }  
         }
         
       }
     }
     
       mu <- predict(model1, newdata = test_data, type = "response", what = "mu")
       sigma <- predict(model1, newdata = test_data, type = "response", what = "sigma")
       nu <- predict(model1, newdata = test_data, type = "response", what = "nu")
      
     if(length(mu)!=dim(test_data)[1])
     {
      print("Error, Please Check Data!!!")
     }
     
     
       Z_score_folds_HC<-zzz_cent(obj=model1,type=c("z-scores"),mu=mu,sigma=sigma,nu=nu,
                             xname = 'Age',xvalues=test_data$Age,yval=test_data$feature,
                             calibration=FALSE,lpar=3)
     
       Z_score_folds_HC<-data.frame(Z_score_folds_HC);
     
       rownames(Z_score_folds_HC)<-data1$individual_ID[folds[[i_fold]]]
     
       Z_score_folds_HC1<-rbind(Z_score_folds_HC1,Z_score_folds_HC)
     
     
       Quant_score_folds_HC<-zzz_cent(obj=model1,type=c("z-scores"),mu=mu,sigma=sigma,nu=nu,
                                 xname = 'Age',xvalues=test_data$Age,yval=test_data$feature,
                                 calibration=FALSE,lpar=3,cdf=TRUE)
       Quant_score_folds_HC<-data.frame(Quant_score_folds_HC);
     
       rownames(Quant_score_folds_HC)<-data1$individual_ID[folds[[i_fold]]]
       Quant_score_folds_HC1<-rbind(Quant_score_folds_HC1,Quant_score_folds_HC)
     
     }
     
     results$Z_score_folds_HC<-Z_score_folds_HC1;
     results$Quant_score_folds_HC<-Quant_score_folds_HC1;
     saveRDS(results,paste0(str,'_',i,'_loop_our_model.rds'))
     
  }
  
}



