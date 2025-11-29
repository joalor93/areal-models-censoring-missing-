This repository contains codes used to develop the article "A Unified Spatiotemporal Framework for Modeling Censored and Missing Areal Responses".


To replicate results such as the CO concentratios application or some toy examples, you should install the auxiliarcpp and TreeSampleR packages. You can do this using the following lines: 

### auxiliarcpp
remotes::install_github("joalor93/areal-models-censoring-missing-", subdir = "auxiliarcpp")

### TreeSampleR
devtools::install_github('https://github.com/edrictam/TreeSampleR')


You would also find the following folders and files attached to them that will allow the replication of the results mentioned before. 

########################################################
##### 1.FOLDER: model_fitting  #########################
########################################################

It contains scripts used for running the models. The samples obtained for each model were saved as .rds archives, which we will describe in the folder "codes summarized". described below. 

The files in this directory are. 

** "beijing_inference_dagarar1.R" code for fitting the DAGAR-AR(1) model
** "beijing_inference_dagarar2.R" code for fitting the DAGAR-AR(2) model
** "beijing_inference_car.R" code for fitting the CAR model



########################################################
##### 2.FOLDER: codes summarized  ######################
########################################################

It contains the results obtained from the models fitted above. Its content is described below. 

** "summarized results.R". This code reproduces the results of Table 4 of the main manuscript. This table shows the estimated models for the CO concentrations in Beijing. It also shows different model criteria for model selection. The following archives are necessary to run this code, they are within the main folder too: 

*** "Estaciones_unidas_rushhour.csv". Contains the database with the response and covariates used to fit the model. 
*** "adj_matrix.rds". Adjacency matrix representing the neighborhood structure of the twelve stations.
*** "cad1ar1dagar_Beijing.rds". The first chain obtained for the DAGAR-AR(1) model. 
*** "cad2ar1dagar_Beijing.rds". The second chain obtained for the DAGAR-AR(1) model.
*** "cad3ar1dagar_Beijing.rds". The third chain obtained for the DAGAR-AR(1) model.
*** "cad1ar2dagar_Beijing.rds". The first chain obtained for the DAGAR-AR(2) model.
*** "car_beijing_2_40000.rds". A chain for the spatiotemporal CAR model (other chains were run, nonetheless, they were not used within this code )
*** "log_lik_ar1dagar.rds". Log likelihood associated to one chain of the DAGAR AR(1) model 
*** "log_lik_ar2dagar.rds". Log likelihood associated to one chain of the DAGAR AR(2) model
*** "log_lik_car.rds". Log likelihood associated to one chain of the CAR model
*** "results_pred_dens_2_ar1dagar.rds". Log predictive density of the test observations. It is associated to one chain of the DAGAR AR(1) model
*** "results_pred_dens_2_ar2dagar.rds". Log predictive density of the test observations. It is associated to one chain of the DAGAR AR(2) model
*** "results_pred_dens_car.rds" Log predictive density of the test observations. It is associated to one chain of the CAR model


** "prediction_graph_ar1dagar.R" . It contains the code to reproduce Figure S.5.16 (DAGAR prediction) of the supplementary material. The following archives are necessary to run this code. 

*** "Estaciones_unidas_rushhour.csv". Contains the database with the response and covariates used to fit the model. 
*** "adj_matrix.rds". Adjacency matrix representing the neighborhood structure of the twelve stations.
*** "ar1dagar_beijing_212_40000.rds". Obtained chains for the DAGAR-AR(1) model.
*** "results_pred_ar1dagar.rds". Predictive distribution for the DAGAR-AR(1) model. 


** "prediction_graph_car.R" . It contains the code to reproduce Figure 5 (CAR prediction) of the supplementary material. The following archives are necessary to run this code. 

*** "Estaciones_unidas_rushhour.csv". Contains the database with the response and covariates used to fit the model. 
*** "adj_matrix.rds". Adjacency matrix representing the neighborhood structure of the twelve stations.
*** "car_beijing_2_40000.rds". A chain for the spatiotemporal CAR model (other chains were run, nonetheless, they were not used within this code )
*** "results_pred_car.rds". Predictive distribution for the CAR model. 


** "assumptions_dagarar1.R". It contains the Gelman's convergence criteria, chains plotted, ACF as well as the posterior distributions of each parameter of the DAGAR - AR(1) model. The following archives are necessary to run this code. 

*** "cad1ar1dagar_Beijing.rds". The first chain obtained for the DAGAR-AR(1) model. 
*** "cad2ar1dagar_Beijing.rds". The second chain obtained for the DAGAR-AR(1) model.
*** "cad3ar1dagar_Beijing.rds". The third chain obtained for the DAGAR-AR(1) model.

** toy_examples.R Toy examples to illustrate the functioning of our routines for fitting areal models with censored/missing responses. 


