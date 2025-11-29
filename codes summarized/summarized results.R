######################################################
########## Lines to install the auxiliarcpp package##
####################################################


#setwd("C:/Users/ALEJANDRO/OneDrive/Documentos/Alejo/postdoc20242/censored_network_dagar/datasets/Beijing/codes sumarized")
#install.packages("auxiliarcpp_0.0.0.9000.tar.gz", type="source")



library(tmvtnorm)
library(doParallel)
library(doSNOW)
library(data.table)
library(dplyr)
library(ggplot2)
library(auxiliarcpp)


##### fixing the work directory #######################

setwd("C:/Users/ALEJANDRO/OneDrive/Documentos/Alejo/postdoc20242/censored_network_dagar/datasets/Beijing/codes sumarized/")

######################################################

###############################################
########## preparing the dataset #############
#############################################

data = fread(paste0("Estaciones_unidas_rushhour.csv"))
adj_matrix = readRDS("adj_matrix.rds")

adj_matcom = adj_matrix + t(adj_matrix)

datapre = data#%>%select("data", "station", "aggpm25","aggpm10" ,
#         "aggtemp","aggrain", "aggsumrain", "aggpres", "aggwsmp")


regions = c("Aotizhongxin", "Changping", "Dingling", "Dongsi",
            "Guanyuan", "Gucheng", "Huairou", "Nongzhanguan",
            "Shunyi", "Tiantan","Wanliu", "Wanshouxigong")




print(apply(datapre,2,function(x) sum(is.na(x))))

datapre2 = datapre%>%group_by("station")%>%
  mutate (aggtemp = if_else(is.na(TEMP),mean(TEMP,na.rm = TRUE),TEMP),
          aggrain = if_else(is.na(aggrain),mean(aggrain,na.rm = TRUE),aggrain),
          aggsumrain = if_else(is.na(aggsumrain),mean(aggsumrain,na.rm = TRUE),aggsumrain),
          aggpres = if_else(is.na(PRES),mean(PRES,na.rm = TRUE),PRES),
          aggwsmp = if_else(is.na(WSPM),mean(aggwsmp,na.rm = TRUE),WSPM))%>%ungroup()

#### prediction time without missing data "data>= "2017-01-24" & data<="2017-01-31"

datapre2 = datapre2%>%filter(data>= "2016-11-25" & data<="2017-02-25")

print(apply(datapre2,2,function(x) sum(is.na(x))))



ytotobs = datapre2$CO
indicens3 = as.numeric(is.na(ytotobs))
lower = rep(-Inf,sum(indicens3 ==1))
upper = rep(Inf,sum(indicens3 ==1))
ytotobs[indicens3 == 1] = mean(ytotobs,na.rm = TRUE)
xtotobs = cbind (1,datapre2$aggtemp - mean(datapre2$aggtemp),datapre2$aggwsmp - mean(datapre2$aggwsmp), (datapre2$aggpres - mean(datapre2$aggpres))/100)
y = log(ytotobs)
N_time_obs = length(ytotobs[datapre2$station == "Shunyi"])
lags = N_time_obs

###################################################################
##### obtained chains: DAGAR AR 1, DAGAR AR 2, CAR ###############
#################################################################

resar1dagar = readRDS("cad1ar1dagar_Beijing.rds")

resultsar2dagar = readRDS("ar2dagar_beijing_21_40000.rds")



resar2dagar = resultsar2dagar[[1]]

rescar =  readRDS("car_beijing_2_40000.rds")



thetatotar1dagar = cbind(resar1dagar$betaF,resar1dagar$sigma2F,resar1dagar$rhoF, resar1dagar$psiF, resar1dagar$phiF)

thetatotar2dagar = cbind(resar2dagar$betaF,resar2dagar$sigma2F,resar2dagar$rhoF, resar2dagar$psiF, resar2dagar$phi1F,resar2dagar$phi2F)

thetatotcar = cbind(rescar$betaF,rescar$sigma2F,rescar$rho_sF, rescar$rho_tF, rescar$rho_stF)


#######################################
######## summarized models ###########
#####################################

names1 = c(paste0("beta_",0:3), "sigma2", "rho", "phi", "rho*phi","tau2")
names2 = c(paste0("beta_",0:3), "sigma2", "rho", "phi_1", "phi_2", "rho*phi", "rho*phi2", "tau2")


#### DAGAR AR 1 ############
rhophi = resar1dagar$rhoF*resar1dagar$phiF
rhophiest = quantile(rhophi,c(0.025,0.5,0.975))

thetaestar1dagar = rbind(t(resar1dagar$betaest),resar1dagar$sigmaest,resar1dagar$rhoest,resar1dagar$phiest,rhophiest,resar1dagar$tau2est)
rownames(thetaestar1dagar) = names1
thetaestar1dagar[, c(1,2)] <- thetaestar1dagar[, c(2,1)]
colnames(thetaestar1dagar) = c( "50%","2.5%", "97.5%")
thetaestar1dagar

##### estimated values ###

round(thetaestar1dagar,3)

#### DAGAR AR 2 ############

rhophi1 = resar2dagar$rhoF*resar2dagar$phi1F
rhophi1est = quantile(rhophi1,c(0.025,0.5,0.975))

rhophi2 = resar2dagar$rhoF*resar2dagar$phi2F
rhophi2est = quantile(rhophi2,c(0.025,0.5,0.975))


thetaestar2dagar = rbind(t(resar2dagar$betaest),resar2dagar$sigmaest,resar2dagar$rhoest,resar2dagar$phi1est,resar2dagar$phi2est,
                         rhophi1est, rhophi2est, resar2dagar$tau2est)
rownames(thetaestar2dagar) = names2
thetaestar2dagar[, c(1,2)] <- thetaestar2dagar[, c(2,1)]
colnames(thetaestar2dagar) = c( "50%","2.5%", "97.5%")

##### estimated values ###

round(thetaestar2dagar,3)

#### CAR ###########

thetaestcar = rbind(t(rescar$betaest),rescar$sigmaest,rescar$rho_s_est,rescar$rho_t_est,rescar$rho_st_est)
names3 = c(paste0("beta_",0:3), "sigma2", "rho_s", "rho_t", "rho_st")
rownames(thetaestcar) = names3
thetaestcar[, c(1,2)] <- thetaestcar[, c(2,1)]
colnames(thetaestcar) = c( "50%","2.5%", "97.5%")

##### estimated values ###

round(thetaestcar,3)



##############################################
##### MODEL SELECTION CRITERIA (LOG-LIKELIHOODS)##############
#############################################


##### log_likelihoods from the chains #########

log_lik_ar1dagar = readRDS("log_lik_ar1dagar.rds")
log_lik_ar2dagar = readRDS("log_lik_ar2dagar.rds")
log_lik_car = readRDS("log_lik_car.rds")

###### expected AIC and expected BIC #######

dbarraar1dagar = -2*mean(log_lik_ar1dagar)
dbarraar2dagar = -2*mean(log_lik_ar2dagar)
dbarracar = -2*mean(log_lik_car)

meanthetaar1dagar = apply(thetatotar1dagar,2,mean)
meanbetaar1dagar = meanthetaar1dagar[1:ncol(xtotobs)]
meansigmaar1dagar = meanthetaar1dagar[ncol(xtotobs)+1]
meanrhoar1dagar =  meanthetaar1dagar[ncol(xtotobs)+2]
meanpsiar1dagar =  meanthetaar1dagar[ncol(xtotobs)+3]
meanphiar1dagar =  meanthetaar1dagar[ncol(xtotobs)+4]


meanthetaar2dagar = apply(thetatotar2dagar,2,mean)
meanbetaar2dagar = meanthetaar2dagar[1:ncol(xtotobs)]
meansigmaar2dagar = meanthetaar2dagar[ncol(xtotobs)+1]
meanrhoar2dagar =  meanthetaar2dagar[ncol(xtotobs)+2]
meanpsiar2dagar =  meanthetaar2dagar[ncol(xtotobs)+3]
meanphi1ar2dagar =  meanthetaar2dagar[ncol(xtotobs)+4]
meanphi2ar2dagar =  meanthetaar2dagar[ncol(xtotobs)+5]


meanthetacar = apply(thetatotcar,2,mean)
meanbetacar = meanthetacar[1:ncol(xtotobs)]
meansigmacar = meanthetacar[ncol(xtotobs)+1]
meanrhocar =  meanthetacar[ncol(xtotobs)+2]
meanpsicar =  meanthetacar[ncol(xtotobs)+3]
meanphicar =  meanthetacar[ncol(xtotobs)+4]





destar1dagar = -2*log_likelihood_ar1_fast_pkg_cpp(x = xtotobs,y = y ,cc = indicens3,lag = lags, adjmatinf = adj_matrix,
                                          theta =  c(meanbetaar1dagar,meansigmaar1dagar,meanrhoar1dagar,meanpsiar1dagar,
                                                     meanphiar1dagar), lower = as.vector(lower), upper = as.vector(upper),
                                          pkg = "auxiliarcpp")


destar2dagar = -2*log_likelihood_ar2_fast_pkg_cpp(x = xtotobs,y = y ,cc = indicens3,lag = lags, adjmatinf = adj_matrix,
                                                  theta =  c(meanbetaar2dagar,meansigmaar2dagar,meanrhoar2dagar,meanpsiar2dagar,
                                                             meanphi1ar2dagar,meanphi2ar2dagar), lower = as.vector(lower), upper = as.vector(upper),
                                                  pkg = "auxiliarcpp",fun_name = "spatimecovar_ar2_cpp")


destcar= -2*log_likelihood_car_fast_cpp(x = xtotobs,y = y ,cc = indicens3,lag = lags, adj_matcom = adj_matcom,
                                            theta =  meanthetacar, lower = as.matrix(lower), upper = as.matrix(upper))



##############################################
##### MODEL SELECTION CRITERIA ##############
#############################################


#### DAGAR AR1 ###########

DICar1dagar = 2*dbarraar1dagar - destar1dagar


EAICar1dagar = dbarraar1dagar + 2*ncol(thetatotar1dagar)

EBICar1dagar = dbarraar1dagar + 2*ncol(thetatotar1dagar)*log(nrow(thetatotar1dagar))

#### DAGAR AR2 ###########

DICar2dagar = 2*dbarraar2dagar - destar2dagar


EAICar2dagar = dbarraar2dagar + 2*ncol(thetatotar2dagar)

EBICar2dagar = dbarraar2dagar + 2*ncol(thetatotar2dagar)*log(nrow(thetatotar2dagar))


#### CAR ###########

DICcar = 2*dbarracar - destcar


EAICcar = dbarracar + 2*ncol(thetatotcar)

EBICcar = dbarracar + 2*ncol(thetatotcar)*log(nrow(thetatotcar))


pred_res_ar1_dagar = readRDS("results_pred_dens_2_ar1dagar.rds")
pred_res_ar2_dagar = readRDS("results_pred_dens_2_ar2dagar.rds")
pred_res_car = readRDS("results_pred_dens_car.rds")

pred1 = t(pred_res_ar1_dagar)
pred2 =t(pred_res_ar2_dagar)
pred3 = t(pred_res_car)


library(loo)

elpddagar = elpd(pred1)
elpddagar2 = elpd(pred2)
elpdcar = elpd(pred3)

ELPD = c(elpddagar$estimates[1,1], elpddagar2$estimates[1,1], elpdcar$estimates[1,1])


##### ALL models criteria #####

DIC = c(DICar1dagar,DICar2dagar,DICcar)

EAIC  = c(EAICar1dagar,EAICar2dagar,EAICcar)

EBIC  = c(EBICar1dagar,EBICar2dagar,EBICcar)


#################################
###### SUMMARIZED CRITERIA #####
################################

criteria = rbind(DIC,EAIC,EBIC,ELPD)

colnames(criteria) = c("DAGAR - AR(1)", "DAGAR - AR(2)", "CAR")

round(criteria,3)
