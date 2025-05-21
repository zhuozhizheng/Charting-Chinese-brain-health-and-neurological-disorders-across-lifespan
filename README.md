Charting Chinese brain health and neurological disorders across the lifespan
Table of contents
1. Introduction (#Introduction)
2. Datasets (#Datasets)
2.1. Example dataset for the normative reference estimation (#Dataset-norms)
2.2. Example disease dataset for the normative reference application (#Dataset-diseases)
2.3. Example individual data for the personalized application of normative reference (#Dataset-individual)
2.4. Example new dataset for the site-calibration of normative reference (#Dataset-new)
2.5. Supplementary/Extended datasets for a better interpretation of the normative reference application (#Dataset-extended)
3. Normative models and DPS models (#Models)
  3.1 Normative models by GAMLSS based on the Desikan–Killiany, HCP-MMP1.0, DU15NET, and Brainnetome parcellation atlases, as well as age censoring models based on Desikan–Killiany atlas (#GAMLSS)
  3.2 Normative models by HBR framework (#HBR) 
  3.3 Disease-specific classification models for DPS estimation (#DPS)
4. Scripts (#Scripts)
  4.1 Required R packages and installation (#Check-and-install-packages.R)
  4.2 Normative curve estimate and milestone determination (#Normative-model-fit.R)
  4.3 Bootstrap analysis of the estimated normative curves (#Bootstrap-normative-mdoel-fit.R)
  4.4 Apply normative model to disease datasets (#Disease-application-normative-model.R)
4.5 Apply normative model to an individual case (#Individual-application-normative-model.R)
4.6 Calibrate normative model by new dataset (#Calibration-normative-model-using-new-dataset.R)
  4.7 Statistical analyses of deviation score across diseases (#Stastical-analysis-deviations-across diseases)
4.8 Othe clinical task related scripts (#DPS-estimation.R; #Clincial-score-predcition.R; #PD-outcome-prediction.R; #Prognosis-risk-stratificaiton.R)
4.9. Individualized brain health report (#Individualized-brain-health-report.Rmd)
5. License(#License)

1. Introduction
The codes on Chinese normative reference construction and their downstream clinical applications. The example datasets could be found in file “Datasets”; the normative models (references) could be found in file “Models”; the main scripts could be found in file “Scripts”; the test outputs by the main scripts using example datasets are provided in file “Test_results”. For other detailed scripts on figure plots and comparison analyses, readers could contact the authors: zhuozhizheng@bjtth.org or liuyaou@bjtth.org. Please note that some source functions (#Source-codes) are from Bethlehem, R.A.I., Seidlitz, J., White, S.R. et al. Brain charts for the human lifespan. Nature 604, 525–533 (2022). https://doi.org/10.1038/s41586-022-04554-y.

This GitHub repository contains the main codes necessary to replicate the Chinese population-specific normative references and their clinical applications in neurological diseases. This repository does not contain the original datasets. We do not have permission to distribute some of the datasets included in the published study. However, many datasets are available upon request from the readers, when the corresponding You Liu liuyaou@bjtth.org received the data request.

2. Datasets
Example datasets for normative curve fitting, application, and calibration have been provided. Please note that these are demonstration datasets only, intended for script testing and validation purposes. They do not represent the original datasets used in our study and should not be utilized for research or clinical applications. Regarding the multicenter datasets for which we currently do not have redistribution rights, researchers interested in obtaining these datasets should contact the study authors directly (Zhizheng Zhuo at zhuozhizheng@bjtth.org or Yaou Liu at liuyaou@bjtth.org) to apply for data sharing permissions.

3. Models
Although we cannot share the all individual-level data (original datasets), we are able to share the outcome of the analysis, namely the normative models/fitted curves.

#GAMLSS contains the all the global and regional normative models fitted by GAMLSS based on the Desikan–Killiany, HCP-MMP1.0, DU15NET, and Brainnetome parcellation atlases. For HCP-MMP1.0, DU15NET, and Brainnetome Atlases, only regional cortical measures were currently released. In addition, age-censoring normative models (samples with ages of 18-70 years) based on Desikan–Killiany atlas are also provided. These age censoring normative references have narrow fitting confidence intervals, and would be more robust for a majority of studies and clinical applications for adults

#HBR contains normative models fitted by HBR framework. Please note currently only the global measures using different distribution functions are released. 

#DPS contains disease-specific classification models for disease propensity score (DPS) estimation for various neurological diseases. 

Note: Due to the space limits of the GitHub, we have uploaded these large files of all normative models into Zenodo. Zhuo, Z. (2025). Chinese brain structure normative models by GAMLSS (V1.0). Zenodo. https://doi.org/10.5281/zenodo.15447127; and Zhuo, Z. (2025). Chinese brain structure normative models by GAMLSS (V1.0). Zenodo. https://doi.org/10.5281/zenodo.15430370. You can assess these models through the links: https://zenodo.org/records/15447127; and https://zenodo.org/uploads/15430370. So if you want to test or use the normative models of DK, HCP-MMP, DU15NET, or Brainnetome, please download these models from Zenodo, and arrange these files as subfiles in file of “Models”, similar as the structure of GAMLSS-DK (we only provide the global normative models in GitHub repository due to the space limits). For HBR models and bootstrap files, you can direct download from Zenodo and put them into file of “Models” in your disk. 

4. Scripts
We have provided the main scripts for the normative reference development and application. For more scripts on the additional statistical analyses and comparative analyses, please contact the authors: zhuozhizheng@bjtth.org or liuyaou@bjtth.org. 
4.1 Required R packages and installation (#Check-and-install-packages)
All the R scripts are based on R version 4.3.2. You can download the R from https://cran.r-project.org/bin/windows/base/old/4.3.2/. For more information on R, please refer to https://www.r-project.org/.

#Check-and-install-packages.R: This script would automatically check and install the required R packages for the normative model fitting, application and plots. You just need run this script in your R. 

4.2 Normative curve estimate, milestone determination, and plot (#Normative-model-fit)
#Normative-model-fit.R: This script is for normative model fitting, peak age estimation and normative curve plot. In addition, 10-folds cross-validation for the deviation score calculation of HCs is also included in this script. Run this script, please update the following file paths
###
#set your own directory
datapath='***/Source-codes/' #Change the directory where you save the Source-codes  
clinical_datapath='**/Dataset-norms/Clinical_vars.csv' #Change the directory where you save the clinical variables
MR_datapath='**/Dataset-norms//MR_measures.xlsx' #Change the directory where you save the MR measures  
savepath='**/V7_DK' #Create and determine the directory where you would save the results  
###

4.3 Bootstrap analysis of the estimated normative curves
#Bootstrap-normative-mdoel-fit.R: As you have constructed the normative models in 4.2, and now need bootstrap analysis of the fitted normative curves and peak ages (here default 1000 bootstraps). You can use this script. 
Run this script, please update the following file paths.
###
#set your own directory
datapath='***/Source-codes/' #Change the directory where you save the Source-codes  
feature_path0='**/V7_DK/' #Chnage the path you have saved your normative models in "Normative-model-fit.R"
savepath=paste0(feature_path0,'bootstrap/') #Create and determine the directory where you would save the results 
###

4.4 Apply normative model to disease datasets
#Disease-application-normative-model.R: Using this script, you can calculate the deviation scores derived from the normative references constructed in 4.2. Please arrange your disease dataset similar as the Dataset-norms. Note that the new HCs which was not trained in the normative models could be arranged here accompanied with disease data for further statistical analysis. 
Run this script, please update the following file paths 
###
#set your own directory
datapath='***/Source-codes/' #Change the directory where you save the Source-codes  
clinical_datapath='**/Dataset-diseases/Clinical_vars.csv' #Change the directory where you save the clinical variables  
MR_datapath='**/Dataset-diseases/MR_measures.xlsx' #Change the directory where you save the MR measures  
savepath1='**/V7_DK' #Chnage the path you have saved your normative models in "Normative-model-fit.R"
savepath='**/Diseases'; #Create and determine the directory where you would save the results  
###

4.5 Apply normative model to an individual case
#Individual-application-normative-model.R: Using this script, you can calculate the individual deviation scores derived from the normative references constructed in 4.2. Please arrange your disease dataset similar as the Dataset-norms. In this setting, additional case1 and case2 should be added before the individual data, as the case 1 and case 2 have different sex and site information, which would avoid the potential error for the factor of sex and site included in the models. For a new individual data, you can just replace the individual data with your new data in the files of “Clinical_vars.csv” and “MR_measures.xlsx”. 

Run this script, please update the following file paths 
###
#set your own directory
datapath='***/Source-codes/' #Change the directory where you save the Source-codes  
clinical_datapath='**/Dataset-individual/Clinical_vars.csv' #Change the directory where you save the clinical variables  
MR_datapath='**/Dataset-individual/MR_measures.xlsx' #Change the directory where you save the MR measures  
savepath1='**/V7_DK' #Chnage the path you have saved your normative models in "Normative-model-fit.R"
savepath='**/Individual'; #Create and determine the directory where you would save the results  
###

4.6 Calibrate normative model by new dataset
#Calibration-normative-model-using-new-dataset.R: This script could help you calibrate the normative models using your own HC datasets, which were not unseen by the original models. This calibrate model could be used for your own disease data or other applications. 
Run this script, please update the following file paths 
###
#set your own directory
datapath='***/Source-codes/' #Change the directory where you save the Source-codes  
clinical_datapath='**/Dataset-new/Clinical_vars.csv' #Change the directory where you save the clinical variables  
MR_datapath='**/Dataset-new/MR_measures.xlsx' #Change the directory where you save the MR measures  
modelpath<-'**/V7_DK' #Chnage the path you have saved your normative models in "Normative-model-fit.R"
savepath='**/newsite'; #Create and determine the directory where you would save the results  
###

4.7 Statistical analyses of deviation score across diseases
#Stastical-analysis-deviations-across diseases: This script could be used for the between groups univariate analyses on the deviation scores. You can adjust the disease variable in the script to match your study. 
Run this script, please update the following file paths 
###
#set your own directory
datapath='***/Source-codes/'#Change the directory where you save the Source-codes 
feature_path0='**/Diseases'; #Chnage the path you have saved your normative models in "Disease-application-normative-model.R"
savepath='**/Results/Statistical_analyses_of_deviation_score_across_diseases'; #Create and determine the directory where you would save the results
###

4.8 Othe clinical task related scripts (#DPS-estimation.R; #Clincial-score-predcition.R; # PD-outcome-prediction.R; #Clinial-prognosis.R)
These scripts are used for specific clinical tasks. If you want to test and use these scripts on your own data, please change the directory similar to above scripts and revise related parameter setting in the script according to your purpose.

4.9. Individualized brain health report (#Individualized-brain-health-report.Rmd)
#Individualized-brain-health-report.Rmd: This script for individualized brain health report is written by R Markdown and current output report is a HTML file, which could be further printed as PDF file.
The data arrange is same as the Dataset-individual. Then you should change the directory of the source codes and data you store in your disk. 
Main changes have been labeled and commented in this script. You can change the directory and adjust the input data according to your datafile.

Run this script, please update the following file paths
###
datapath='E:/Lifespan_freesurfer_results/Github/Source-codes/' #Change the directory where you save the Source-codes 
datafile ='E:/Lifespan_freesurfer_results/Github/Datasets/Dataset-individual/MR_measures.xlsx' #Change the directory where you save the clinical variables  
clinicalfile='E:/Lifespan_freesurfer_results/Github/Datasets/Dataset-individual/Clinical_vars.csv' #Change the directory where you save the MR measures  
feature_path0='E:/Lifespan_freesurfer_results/Github/Models/GAMLSS/DK/' #Chnage the path you have saved your normative models in "Normative-model-fit.R"
model_path='E:/Lifespan_freesurfer_results/Github/Models/DPS' #Change the directory where you save the DPS models
savepath='E:/Lifespan_freesurfer_results/Github/Test_results/case_report'; #Create and determine the directory where you would save the results   
###

5.License(#License)
MIT License
Copyright (c) [2025] [Version V1.0]
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

