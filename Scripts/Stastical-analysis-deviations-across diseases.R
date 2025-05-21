#install.packages('gamlss')
rm(list=ls())

library(gamlss)
library(ggplot2)
library(reshape2)


#set your own directory
datapath='***/Source-codes/'#Change the directory where you save the Source-codes 
feature_path0='**/Diseases'; #Chnage the path you have saved your normative models in "Disease-application-normative-model.R"
savepath='**/Results/Statistical_analyses_of_deviation_score_across_diseases';# Create and determine the directory where you would save the results  


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
library(coin)
library(bayestestR)
library(MatchIt)
library(dplyr)

tem_zzz_feature<-c('Global_feature','aseg.vol.table',
                   'lh.aparc.area.table','rh.aparc.area.table',
                   'lh.aparc.thickness.table','rh.aparc.thickness.table',
                   'lh.aparc.volume.table','rh.aparc.volume.table')



for(i_zzz_feature in tem_zzz_feature)
{

  feature_path=paste0(feature_path0,i_zzz_feature)# change
  
setwd(feature_path)
myfile <- list.files()
RDSfile_all <- myfile[grep(myfile,pattern ="our_model.rds$")]   
RDSfile_all


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

for(RDSfile_loop in RDSfile_all)
{
  
  if(RDSfile_loop %in% Exd_var){next}
 
  
  if (i_zzz_feature %in% c('aseg.vol.table','lh.aparc.volume.table','rh.aparc.volume.table','Global_feature'))
     {
      results0<-readRDS(paste0(feature_path0,'Global_feature/Global_feature_TCV_loop_our_model.rds'))
     }
 
     if (i_zzz_feature %in% c('lh.aparc.area.table','rh.aparc.area.table'))
        {
          results0<-readRDS(paste0(feature_path0,'Global_feature/Global_feature_total_surface_arrea_loop_our_model.rds'))
     }
     
     
     if (i_zzz_feature %in% c('lh.aparc.thickness.table','rh.aparc.thickness.table'))
     {
          results0<-readRDS(paste0(feature_path0,'Global_feature/Global_feature_mean_thickness_loop_our_model.rds'))
     }
     
    
     
  feature0<-results0$i
  
  data0<-results0$Quant_data[feature0];
  
  data0<-cbind(results0$all_data,data0)
  
  colnames(data0)[dim(data0)[2]]<-feature0
  
  
  
  
  RDSfile<-RDSfile_loop
  
  setwd(feature_path);
  results<-readRDS(RDSfile)

  
  disease<-c('MCI','AD',"PD","SVD",
             "MS","AQP4Pos_NMOSD")
  
  feature<-results$i
  
  data<-results$Quant_data[feature];
  data<-data.frame(data)
 
 
  ind_sel<-intersect(rownames(results$all_data_original),rownames(data))
  
  data<-cbind(results$all_data_original[ind_sel,],data)
  

  colnames(data)[dim(data)[2]]<-feature
  
  ind_sel<-intersect(rownames(data),rownames(data0))
  
  data<-cbind(data[ind_sel,],data0[ind_sel,feature0])
  
  colnames(data)[dim(data)[2]]<-'global_feature'

  data<-data[(data$Diagnosis=='HC'|data$Diagnosis=='MCI'|
               data$Diagnosis=='AD'|data$Diagnosis=='PD'|
               data$Diagnosis=='SVD'|data$Diagnosis=='MS'|
               data$Diagnosis=='AQP4Pos_NMOSD')&!is.na(data$Diagnosis),]
  
  if (!(i_zzz_feature %in% c('Global_feature'))&!(feature %in% c('Brain.Stem','cerebellum_total')))
  {
  res<-lm(data[,feature]~data[,"global_feature"])

  data[,feature]<-0.5+res$residuals
  }
  
  print(RDSfile_loop)
  data<-select(data,-c("global_feature"))
  
  
  #data<-data[data$Sex=="Male"&!is.na(data$Sex),]# 
  
  data1<-data[data$Diagnosis=='HC',];
  
  list<-data.frame(matrix(0,length(disease),36))
  list_pair<-data.frame(matrix(0,length(disease)*length(disease)-length(disease),37))
  
  colnames(list)<-c('feature','disease','pvalue','tvalue','degree','Fvalue','Adjust_R2',
                    'cohens_d','cohens_d_CI_low','cohens_d_CI_high',
                    'hedges_g','hedges_g_CI_low','hedges_g_CI_high',
                    'glass_delta','glass_delta_CI_low','glass_delta_CI_high',
                    'Eta2','Eta2_CI_low','Eta2_CI_high',
                    'OR','OR_CI_low','OR_CI_high',
                    'Overlap','Overlap_CI_low','Overlap_CI_high',
                    'rank_biserial','rank_biserial_CI_low','rank_biserial_CI_high',
                    'cliffs_delta','cliffs_delta_CI_low','cliffs_delta_CI_high',
                    'cohens_u1','cohens_u1_CI_low','cohens_u1_CI_high',
                    'wilcox_test_W','wilcox_test_p')
  
  colnames(list_pair)<-c('feature','disease1','disease2','pvalue','tvalue','degree','Fvalue','Adjust_R2',
                    'cohens_d','cohens_d_CI_low','cohens_d_CI_high',
                    'hedges_g','hedges_g_CI_low','hedges_g_CI_high',
                    'glass_delta','glass_delta_CI_low','glass_delta_CI_high',
                    'Eta2','Eta2_CI_low','Eta2_CI_high',
                    'OR','OR_CI_low','OR_CI_high',
                    'Overlap','Overlap_CI_low','Overlap_CI_high',
                    'rank_biserial','rank_biserial_CI_low','rank_biserial_CI_high',
                    'cliffs_delta','cliffs_delta_CI_low','cliffs_delta_CI_high',
                    'cohens_u1','cohens_u1_CI_low','cohens_u1_CI_high',
                    'wilcox_test_W','wilcox_test_p')
  

  data1<-data1[,c(1:18,dim(data1)[2])]
  data<-data[,c(1:18,dim(data)[2])]
  
  data0<-data1 
  
  num<-0
  
  for(i in disease)
  {
    num=num+1;
    
    data1<-data0;
    data1<-rbind(data1,data[data$Diagnosis==i,])
    
    
   
    tem_data1<-data1[data1$Age>=min(data1[data1$Diagnosis==i,'Age'],na.rm = TRUE)&
                       data1$Age<=max(data1[data1$Diagnosis==i,'Age'],na.rm = TRUE),]
   
    tem_data1<-na.omit(tem_data1)
    
    tem_data1$Diagnosis<-factor(tem_data1$Diagnosis,levels=c('HC',i))
    match.it <- matchit(Diagnosis ~ Age + Sex +Site_ZZZ, data = tem_data1, method="nearest", ratio=1)

    tem_data1<-tem_data1[c(rownames(match.it$match.matrix),match.it$match.matrix),]
    
    tem_data1<-na.omit(tem_data1)
  
    tem_data1$Diagnosis<-factor(tem_data1$Diagnosis,levels=c('HC',i))
    tem_stat<-lm(tem_data1[,feature]~tem_data1[,'Diagnosis'],data=tem_data1,na.action="na.omit",method = "qr")
    #tem_stat<-oneway_test(tem_data1[,feature]~tem_data1[,'Diagnosis'],data=tem_data1,
    #              distribution=approximate(nresample=1000)) #mento carlo simulation test 
    
    tem_stats<-summary(tem_stat)
    
    cohens_d(tem_data1[,feature]~tem_data1[,'Diagnosis'],data=tem_data1,
             distribution=approximate(nresample=1000))
    tem_d<-cohens_d(tem_data1[,feature]~tem_data1[,'Diagnosis'],data=tem_data1,na.action="na.omit",method = "qr")
    tem_g<-hedges_g(tem_data1[,feature]~tem_data1[,'Diagnosis'],data=tem_data1,na.action="na.omit",method = "qr")
    tem_delta<-glass_delta(tem_data1[,feature]~tem_data1[,'Diagnosis'],data=tem_data1,na.action="na.omit",method = "qr")
    
    model <- aov(tem_data1[,feature]~tem_data1[,'Diagnosis'],data=tem_data1,na.action="na.omit")
    
    tem_eta2<-eta_squared(model)
    
    G1<-tem_data1[tem_data1$Diagnosis=='HC',feature];
    G2<-tem_data1[tem_data1$Diagnosis==i,feature];
    
    tem_stats_nonpara<-wilcox.test(G1,G2,na.action="na.omit",paired=FALSE)
    
    tem_op<-p_overlap(G1,G2,parametric = FALSE)
    
    tem_rb<-rank_biserial(G1,G2,paired = FALSE);
    tem_cld<-cliffs_delta(G1,G2,paired = FALSE)
    tem_u1<-cohens_u1(G1,G2,parametric = TRUE)
    
    
    
    tem_data1[tem_data1$Diagnosis=='HC','label']<-0
    tem_data1[tem_data1$Diagnosis==i,'label']<-1
    
    tem_logs<-glm(tem_data1[,'label']~tem_data1[,feature],family=binomial(link='logit'),
                  data=tem_data1,na.action="na.omit")
    tem_logs_CI<-confint(tem_logs)
    
    
    list[num,1]<-RDSfile_loop
    list[num,2]<-i
    list[num,3]<-tem_stats$coefficients[2,3]
    list[num,4]<-tem_stats$coefficients[2,4]
    list[num,5]<-tem_stats$df[2]
    list[num,6]<-tem_stats$fstatistic[1]
    list[num,7]<-tem_stats$adj.r.squared
    list[num,8]<-tem_d$Cohens_d
    list[num,9]<-tem_d$CI_low
    list[num,10]<-tem_d$CI_high
    list[num,11]<-tem_g$Hedges_g
    list[num,12]<-tem_g$CI_low
    list[num,13]<-tem_g$CI_high
    
    list[num,14]<-tem_delta$Glass_delta
    list[num,15]<-tem_delta$CI_low
    list[num,16]<-tem_delta$CI_high
    
    list[num,17]<-tem_eta2$Eta2
    list[num,18]<-tem_eta2$CI_low
    list[num,19]<-tem_eta2$CI_high
    
    list[num,20]<-tem_logs$coefficients[2]
    list[num,21]<-tem_logs_CI[2,1]
    list[num,22]<-tem_logs_CI[2,2]
    
    list[num,23]<-tem_op$Overlap;
    list[num,24]<-tem_op$CI_low
    list[num,25]<-tem_op$CI_high
    
    list[num,26]<-tem_rb$r_rank_biserial;
    list[num,27]<-tem_rb$CI_low
    list[num,28]<-tem_rb$CI_high
    
    list[num,29]<-tem_cld$r_rank_biserial;
    list[num,30]<-tem_cld$CI_low
    list[num,31]<-tem_cld$CI_high
    
    list[num,32]<-tem_u1$Cohens_U1;
    list[num,33]<-tem_u1$CI_low
    list[num,34]<-tem_u1$CI_high
    
    list[num,35]<-tem_stats_nonpara$statistic
    list[num,36]<-tem_stats_nonpara$p.value
     #tem_delta<-glass_delta(tem_stats)
    
  }
  
  results$disease_comparison_global_regression<-list;
  
  
  saveRDS(results,RDSfile_loop)
  
  num1<-0
  
  for(i1 in disease)
  {
    for(i2 in disease)
    {
      if(i1!=i2){
    num1=num1+1;
    tem_data1<-data[data$Diagnosis==i1|data$Diagnosis==i2,];
   
    tem_data1$Diagnosis<-factor(tem_data1$Diagnosis,levels=c(i1,i2))
    
    tem_data1<-na.omit(tem_data1)
    #match.it <- matchit(Diagnosis ~ Age + Sex + Site_ZZZ, data = tem_data1, method="nearest", ratio=1)
    #match.it <- matchit(Diagnosis ~ Age + Sex, data = tem_data1, method="nearest", ratio=1)
    match.it <- matchit(Diagnosis ~ Sex+Site_ZZZ, data = tem_data1, method="nearest", ratio=1)
    #a <- summary(match.it)
    tem_data1<-tem_data1[c(rownames(match.it$match.matrix),match.it$match.matrix),]
    
    tem_data1<-na.omit(tem_data1)
   
    
    tem_stat<-lm(tem_data1[,feature]~tem_data1[,'Diagnosis'],data=tem_data1,na.action="na.omit",method = "qr")
    #tem_stat<-oneway_test(tem_data1[,feature]~tem_data1[,'Diagnosis'],data=tem_data1,
    #              distribution=approximate(nresample=1000)) #mento carlo simulation test 
    
    tem_stats<-summary(tem_stat)
    
    cohens_d(tem_data1[,feature]~tem_data1[,'Diagnosis'],data=tem_data1,
             distribution=approximate(nresample=1000))
    tem_d<-cohens_d(tem_data1[,feature]~tem_data1[,'Diagnosis'],data=tem_data1,na.action="na.omit",method = "qr")
    tem_g<-hedges_g(tem_data1[,feature]~tem_data1[,'Diagnosis'],data=tem_data1,na.action="na.omit",method = "qr")
    tem_delta<-glass_delta(tem_data1[,feature]~tem_data1[,'Diagnosis'],data=tem_data1,na.action="na.omit",method = "qr")
    
    model <- aov(tem_data1[,feature]~tem_data1[,'Diagnosis'],data=tem_data1,na.action="na.omit")
    
    tem_eta2<-eta_squared(model)
    
    tem_data1[tem_data1$Diagnosis==i1,'label']<-0
    tem_data1[tem_data1$Diagnosis==i2,'label']<-1
    
    tem_logs<-glm(tem_data1[,'label']~tem_data1[,feature],family=binomial(link='logit'),
                  data=tem_data1,na.action="na.omit")
    tem_logs_CI<-confint(tem_logs)
    
    
    
    G1<-tem_data1[tem_data1$Diagnosis==i1,feature];
    G2<-tem_data1[tem_data1$Diagnosis==i2,feature];
    
    tem_stats_nonpara<-wilcox.test(G1,G2,na.action="na.omit",paired=FALSE)
    
    tem_op<-p_overlap(G1,G2,parametric = FALSE)
    
    tem_rb<-rank_biserial(G1,G2,paired = FALSE);
    tem_cld<-cliffs_delta(G1,G2,paired = FALSE)
    tem_u1<-cohens_u1(G1,G2,parametric = TRUE)
    
  

    list_pair[num1,1]<-RDSfile_loop
    list_pair[num1,2]<-i1
    list_pair[num1,3]<-i2
    list_pair[num1,4]<-tem_stats$coefficients[2,3]
    list_pair[num1,5]<-tem_stats$coefficients[2,4]
    list_pair[num1,6]<-tem_stats$df[2]
    list_pair[num1,7]<-tem_stats$fstatistic[1]
    list_pair[num1,8]<-tem_stats$adj.r.squared
    list_pair[num1,9]<-tem_d$Cohens_d
    list_pair[num1,10]<-tem_d$CI_low
    list_pair[num1,11]<-tem_d$CI_high
    list_pair[num1,12]<-tem_g$Hedges_g
    list_pair[num1,13]<-tem_g$CI_low
    list_pair[num1,14]<-tem_g$CI_high
    
    list_pair[num1,15]<-tem_delta$Glass_delta
    list_pair[num1,16]<-tem_delta$CI_low
    list_pair[num1,17]<-tem_delta$CI_high
    
    list_pair[num1,18]<-tem_eta2$Eta2
    list_pair[num1,19]<-tem_eta2$CI_low
    list_pair[num1,20]<-tem_eta2$CI_high
    
    list_pair[num1,21]<-tem_logs$coefficients[2]
    list_pair[num1,22]<-tem_logs_CI[2,1]
    list_pair[num1,23]<-tem_logs_CI[2,2]
 
   
    list_pair[num1,24]<-tem_op$Overlap;
    list_pair[num1,25]<-tem_op$CI_low
    list_pair[num1,26]<-tem_op$CI_high
    
    list_pair[num1,27]<-tem_rb$r_rank_biserial;
    list_pair[num1,28]<-tem_rb$CI_low
    list_pair[num1,29]<-tem_rb$CI_high
    
    list_pair[num1,30]<-tem_cld$r_rank_biserial;
    list_pair[num1,31]<-tem_cld$CI_low
    list_pair[num1,32]<-tem_cld$CI_high
    
    list_pair[num1,33]<-tem_u1$Cohens_U1;
    list_pair[num1,34]<-tem_u1$CI_low
    list_pair[num1,35]<-tem_u1$CI_high
    
    list_pair[num1,36]<-tem_stats_nonpara$statistic
    list_pair[num1,37]<-tem_stats_nonpara$p.value
    
    #tem_delta<-glass_delta(tem_stats)
    }}}
  
  
  
  #results$disease_comparison<-list;
  
  results$disease_comparison_paired_global_regression<-list_pair;
  
  saveRDS(results,RDSfile_loop)
  
  
  
  data1<-data0
  for(i in disease)
  {
    data1<-rbind(data1,data[data$Diagnosis==i,])
  }
  
  data1$Diagnosis<-factor(data1$Diagnosis,levels=c('HC',disease))
  colnames(data1)[dim(data1)[2]]<-c('feature')
  # ggplot(data=data1,aes(Diagnosis, feature, fill = Diagnosis)) + 
  #   geom_rain(alpha = .5)+
  #   theme_classic()
 
  setwd(savepath)
  data1<-data1[data1$Diagnosis!='HC',]
  
  data1<-data1[!is.na(data1$Diagnosis),]

  
  png(filename = paste0(results$str,'_',results$i,'_between_group_comparison.png'), 
      width = 2500,           
      height = 1400,          
      units = "px",          
      bg = "white",          
      res = 300)  
 
  plot01<-ggplot(data=data1,aes(Diagnosis, feature)) + 
    # ggdist::stat_halfeye(aes(color=Diagnosis,fill=Diagnosis),alpha=0.8,adjust =1,position = 'identity',
    #                      width = 1, .width = c(0,1), justification = 0, point_colour = NA,density = "unbounded") + 
    
    ggdist::stat_halfeye(aes(color=Diagnosis,fill=Diagnosis),alpha=0.5,adjust =1,#position = 'dodge',
                         width = 0.5, .width = 0, scale=1,justification = 0, point_colour = NA,density = "unbounded") + 
    
    #geom_boxplot(aes(color=Diagnosis),width = .2) + 
    stat_summary(fun = "median", geom = "point",
                 shape=16,size=10,alpha=1,
                 aes(color = Diagnosis))+
    stat_summary(fun = "median", geom = "point",
                 shape=16,size=5,alpha=1,
                 color = 'white')+
    #geom_jitter(aes(color=Diagnosis),width = .05, alpha = .3) +
    ggsci::scale_color_npg()+
    ggsci::scale_fill_npg() +
    #coord_flip()+
    #ggdist::stat_dots(aes(color=Diagnosis,fill=Diagnosis),side = "left", dotsize = .5, justification = 1.1, binwidth = .1)+
    labs(title=paste0(results$str,'_',results$i),x='',y='')+
    geom_hline(aes(yintercept=0.5),linetype=c('solid'),linewidth=1,alpha=0.3)+
    scale_y_continuous(limits = c(0, 1),  
                       breaks = seq(0, 1, by = 0.5)) +
    
    #geom_text(data=error2,aes(x=group,y=values+9,label=label),color="red",fontface="bold")+
    # theme(panel.background = element_blank(),
    #       axis.text=element_text(size=16),
    #       axis.line = element_line(colour = "#000000",size = 0.2),
    #       axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))+
    theme(legend.position = 'none',
          panel.background = element_blank(),
          axis.text=element_text(size=16),
          axis.line = element_line(colour = "#000000",size = 0.2))+
    scale_x_discrete(labels = c("MCI" = "MCI","AD" = "AD", "PD" = "PD","SVD" = "CSVD","MS" = "MS",'AQP4Pos_NMOSD'="NMOSD"))+
    theme(
      axis.title = element_text(family = "serif",size=16,color = "black"),
      axis.text.x = element_text(
        size = 16,              
        color = "black",          
        family = "serif"      
      ),
      axis.text.y = element_text(
        size = 16,             
        color = "black",         
        #face = "bold" ,          
        family = "serif"
      )
    )
  
    #hrbrthemes::theme_ipsum(base_family = "Roboto Condensed") +
    # theme(
    #   plot.title = element_markdown(hjust = 0.5,vjust = .5,color = "black",
    #                                 size = 20, margin = margin(t = 1, b = 12)),
    #   plot.subtitle = element_markdown(hjust = 0,vjust = .5,size=15),
    #   plot.caption = element_markdown(face = 'bold',size = 12),
    #   legend.position = "none")
  print(plot01)
  dev.off()
}


}
  
  


