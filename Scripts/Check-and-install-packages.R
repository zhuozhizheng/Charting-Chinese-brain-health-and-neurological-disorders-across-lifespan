#check and install the main required packages
packages_list<-c("readxl","stringr","reshape2","MatchIt","dplyr","rsample",
                          "gamlss","lmerTest","effectsize","caret","coin","bayestestR","hrbrthemes",
                           "ppcor","e1071","kernlab","autoReg","survminer","commonmark","lsr",
                          "ggplot2","ggsci",  "raincloudplots","ggrain","ggdist","ggseg","ggpubrr",
                          "knitr")


for(i_package in packages_list)
{
if (!require(i_package, quietly = TRUE)) {
  install.packages(i_package)
}
}
