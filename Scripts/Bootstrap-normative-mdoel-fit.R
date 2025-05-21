
rm(list=ls())

#set your own directory
datapath='E:/Lifespan_freesurfer_results/Github/Source-codes/'  # Change the directory where you save the Source-codes  
feature_path0='E:/Lifespan_freesurfer_results/Github/Test_results/V7_DK/' #Chnage the path you have saved your normative models in "Normative-model-fit.R"
savepath=paste0(feature_path0,'bootstrap/') # Create and determine the directory where you would save the results  

setwd(datapath)
source("100.common-variables.r")
source("101.common-functions.r")
source("ZZZ_function.R")
source("300.variables.r")
source("301.functions.r")



library(stringr)
library(gamlss)
library(reshape2)
library(ggplot2)
library(ggsci)


path<-list()

path[[1]]<-paste0(feature_path0,'Global_feature');sheet='Global_feature'
#path[[1]]<-paste0(feature_path0,'aseg.vol.table');sheet='aseg.vol.table'
#path[[1]]<-paste0(feature_path0,'lh.aparc.thickness.table');sheet='lh.aparc.thickness.table'
#path[[1]]<-paste0(feature_path0,'rh.aparc.thickness.table');sheet='rh.aparc.thickness.table'
#path[[1]]<-paste0(feature_path0,'lh.aparc.area.table');sheet='lh.aparc.area.table'
#path[[1]]<-paste0(feature_path0,'rh.aparc.area.table');sheet='rh.aparc.area.table'
#path[[1]]<-paste0(feature_path0,'lh.aparc.volume.table');sheet='lh.aparc.volume.table'
#path[[1]]<-paste0(feature_path0,'rh.aparc.volume.table');sheet='rh.aparc.volume.table'



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
           #"aseg.vol.table_Brain.Stem_loop_our_model.rds",
           "lh.aparc.thickness.table_eTIV_loop_our_model.rds",
           "rh.aparc.thickness.table_eTIV_loop_our_model.rds",
           "lh.aparc.area.table_eTIV_loop_our_model.rds",
           "rh.aparc.area.table_eTIV_loop_our_model.rds",
           "lh.aparc.volume.table_eTIV_loop_our_model.rds",
           "rh.aparc.volume.table_eTIV_loop_our_model.rds",
           "lh.aparc.area.table_rh_WhiteSurfArea_area_loop_our_model.rds",
           "rh.aparc.area.table_rh_WhiteSurfArea_area_loop_our_model.rds",
           "lh.aparc.thickness.table_rh_MeanThickness_thickness_loop_our_model.rds",
           "rh.aparc.thickness.table_rh_MeanThickness_thickness_loop_our_model.rds")


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
RDSfile_all<-c(RDSfile_all,paste0(feature_path,'/',j));
}
}
}


for(RDSfile_loop in RDSfile_all)
{
  
  if(!(RDSfile_loop %in% Exd_var))
    {
    print(RDSfile_loop)
    RDSfile<-RDSfile_loop
    results<-readRDS(RDSfile)
    
    

    #resample HCs for the robustness of fitted curves and peak ages

    sample_num=1000

    con=gamlss.control(c.crit = 0.01, n.cyc = 5,autostep = FALSE)
    p2_sample=results$p2[,c('Age','median','sigma')];
    p2_sample[,'Index']<-'Original'
    peak_age_sample=data.frame(results$peakage);
    colnames(peak_age_sample)<-'peakage';
    peak_age_sample[,'Index']<-'Original'

    data1<-results$data1;
    m0=results$m0
    m2=results$m2
    str<-results$str
    i<-results$i
    if(file.exists(paste0(str,'_',i,'_all_without_sex_stratified_resample_CI.png')))
    {
      print('file exist')
      next}
    
    i_rnd=1;
    j_rnd=1
    if(is.null(m2$mu.coefSmo[[1]]$coef))
    {i_rnd=0}
    
    if(is.null(m2$sigma.coefSmo[[1]]$coef))
    {j_rnd=0}
      
    results1<-list()
    
    for (i_sample in 1:sample_num)
    {

    results_try<-try({
    print(i_sample)
    data_resample<-data1[sample(nrow(data1), size=nrow(data1), replace = TRUE),]

    #for all population plot
    if(i_rnd==1&j_rnd==1){
      m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                 sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                 control=con,
                 family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                 data=data_resample)}else if(i_rnd==1&j_rnd==0){
                   m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                              sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex,
                              control=con,
                              family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                              data=data_resample)}else if(i_rnd==0&j_rnd==1){
                                m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex,
                                           sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                                           control=con,
                                           family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                           data=data_resample)}else if(i_rnd==0&j_rnd==0){
                                             m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex,
                                                        sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex,
                                                        control=con,
                                                        family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                                        data=data_resample)}




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


    colnames(p2)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma');
    p2<-p2[,c('Age','median','sigma')]
    p2[,'Index']<-paste0('resample',i_sample)
    p2_sample<-rbind(p2_sample,p2);


    peak_age=data.frame(p2$Age[which.max(p2$median)]);
    colnames(peak_age)<-'peakage';
    peak_age[,'Index']<-paste0('resample',i_sample)
    peak_age_sample<-rbind(peak_age_sample,peak_age)
    })


    if(inherits(results_try,'try-error'))
    {next}

    }


    lower_peak_age =  quantile(peak_age_sample$peakage, 0.025)
    upper_peak_age =  quantile(peak_age_sample$peakage, 0.975)
    median_peak_age = median(peak_age_sample$peakage)

    results1$p2_sample<-p2_sample
    results$peak_age_sample<-peak_age_sample

    results$CI<-c(median_peak_age,lower_peak_age,upper_peak_age)

    
    
    #calculate the CI of the fitted curves
    p2_sample_CI<-NULL

    for(i_sample in unique(p2_sample$Index))
    {

      p2_sample_CI<-cbind(p2_sample_CI,p2_sample[p2_sample$Index==i_sample,'median'])
    }


    ci_data <- data.frame(
      Age = p2_sample$Age[p2_sample$Index=="Original"],
      lower = apply(p2_sample_CI, 1, function(x) quantile(x, 0.025)),
      upper = apply(p2_sample_CI, 1, function(x) quantile(x, 0.975)),
      median = p2_sample[p2_sample$Index=="Original",'median']
    )


    png(filename = paste0(str,'_',i,'_all_without_sex_stratified_resample_CI.png'),
        width = 1480,
        height = 740,
        units = "px",
        bg = "white",
        res = 300)

    p3<-ggplot()+
      #geom_line(data=mydata,aes(x=Age,y=Value,group=Percentile,color=Percentile))+
      # geom_point(data=data1[data1$Sex=='Female',],aes(x=Age,y=tem_feature/scale1),
      #            colour=c('#E84935'),shape=16,size=3,alpha = 0.1)+
      # geom_point(data=data1[data1$Sex=='Male',],aes(x=Age,y=tem_feature/scale1),
      #            colour=c('#4FBBD8'),shape=17,size=3,alpha = 0.1)+

      geom_line(data=ci_data,aes(x=Age,y=median/scale1),color='black',linewidth=1,linetype=c('solid'))+
      geom_line(data=ci_data,aes(x=Age,y=lower/scale1),color='black',linewidth=0.5,linetype=c('dashed'))+
      geom_line(data=ci_data,aes(x=Age,y=upper/scale1),color='black',linewidth=0.5,linetype=c('dashed'))+
      labs(title=paste0(i,' ',ylab1),x='',y='')+
      theme_bw()+
      theme(
        legend.position = 'none',
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
    
    
    results1$p2_sample_CI<-p2_sample_CI
    results$ci_data<-ci_data
    
   
    saveRDS(results,RDSfile)
    
    saveRDS(results1,paste0(savepath,str,'_',i,'_loop_our_model_bootstrap.rds'))
    

  }
}  


