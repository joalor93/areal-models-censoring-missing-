
######################################################
########## Lines to install the auxiliarcpp package##
####################################################


#setwd("C:/Users/ALEJANDRO/OneDrive/Documentos/Alejo/postdoc20242/censored_network_dagar/datasets/Beijing/codes sumarized")
#install.packages("auxiliarcpp_0.0.0.9000.tar.gz", type="source")



library(mvtnorm)
library(invgamma)
library(auxiliarcpp)
library(tmvtnorm)
library(dplyr)

censlevel =0.15
misslevel = 0.05

m=4
N_time = 16
conn_graph = inclattice(m)
N_spa = nrow(conn_graph)

adj_matrix = conn_graph

adj_matrix[upper.tri(adj_matrix)] = 0


sigma2 = 2
beta = c(1,2,2.5)
rho = 0.8 ## spatial parameter rho
phi = 0.7## time parameter gamma
tau2 = 0.6
psi = tau2/sigma2


set.seed(123)
#########################################################
###### simuulating a  DAGAR Gaussian process #########
#######################################################

covtot = sigma2* spatimecovar_2_cpp(N_time,adj_matrix,rho,psi,phi) ### DAGAR covariance
xtot = cbind(1,rnorm(N_spa*N_time,4,2),rnorm(N_spa*N_time,2,2))
mu = xtot%*%beta
##mu = as.matrix(mu,t,N)
p =length(beta)
ydagar1 = as.vector(rmvnorm(1,mean = mu, sigma = covtot)) ### simulated response

#### generating censoring and missingness to the response #########
ydagarcen2 = quantile(ydagar1,probs = censlevel)
ydagarcen2 = ydagarcen2 + runif(1)
ydagar1[ydagar1<ydagarcen2] = ydagarcen2
indicens2 = 0 + (ydagar1 == ydagarcen2)
missobs = sample(which(indicens2 ==0), size = misslevel*length(ydagar1[indicens2==0]))
indimiss = rep(0,length(ydagar1))
indimiss[missobs] = 1
missimpvalue =  mean(ydagar1) #sample(c(-9999,9999),1)
ydagar1[indimiss==1] = missimpvalue

##### lower and upper limits #######
upper = ydagar1
upper[indicens2 ==1] = ydagarcen2
upper[indimiss ==1] = Inf
indicens3 = indicens2 + indimiss
upper = upper[indicens3==1]
lower = rep(-Inf,length(upper))

####### arguments to test the functions #########
thetaini = c(0.4,0.4,0.4)
lags = N_time
aprior = c(3.5,2,3.5,1)
bprior = c(2.5,4,2.5,28)
divproprho = 6
divproppsi = 6
divpropphi = 6
iter = 6000
burn = 600
thin = 12

divproprhosar =20
divproppsisar =20
divpropphisar = 20



#####Bayesian MCMC DAGAR ######

start <- Sys.time()

resar1dagar = bayesspatempcensauto_fix(ydagar1,xtot,thetaini,indicens3, iter,burn,thin,lags = N_time,
                                       adjmatinf  = adj_matrix, aprior = aprior, bprior = bprior,
                                       divproprho = divproprho, divproppsi = divproppsi, divpropphi = divpropphi,
                                       lower,upper)

print( Sys.time() - start)



#########################################################
###### simuulating a  SAR Gaussian process #############
#######################################################

covtotsar = sigma2*spatimecovarsar_2_cpp(lag = N_time, adjmat= conn_graph,
                                   rho = rho, psi = psi, phi = phi) ### SAR covariance matrix


ysar1 = as.vector(rmvnorm(1,mean = mu, sigma = covtotsar))### simulated response

#### generating censoring and missingness to the response #########
ysarcen2 = quantile(ysar1,probs = censlevel)
ysarcen2 = ysarcen2 + runif(1)
ysar1[ysar1<ysarcen2] = ysarcen2
indicens2 = 0 + (ysar1 == ysarcen2)
missobs = sample(which(indicens2 ==0), size = 0.05*length(ysar1[indicens2==0]))
indimiss = rep(0,length(ysar1))
indimiss[missobs] = 1
missimpvalue =  mean(ysar1)
ysar1[indimiss==1] = missimpvalue

##### lower and upper limits #######
upper = ysar1
upper[indicens2 ==1] = ysarcen2
upper[indimiss ==1] = Inf
indicens3 = indicens2 + indimiss
upper = upper[indicens3==1]
lower = rep(-Inf,length(upper))


#########################################################


#####Bayesian MCMC SAR ######

start <- Sys.time()

resar1sar = bayesspatempcensautosar_fix(ysar1,xtot,thetaini,indicens3, iter,burn,thin,lags = N_time,
                               adjmat  = conn_graph, aprior = aprior, bprior = bprior,
                               divproprho = divproprhosar, divproppsi = divproppsisar, divpropphi = divpropphisar,
                               lower,upper)

print( Sys.time() - start)




#######################################
######## summarized models ###########
#####################################

names1 = c(paste0("beta_",0:2), "sigma2", "rho", "phi", "rho*phi","tau2")

#### DAGAR AR 1 ############
rhophi = resar1dagar$rhoF*resar1dagar$phiF
rhophiest = quantile(rhophi,c(0.025,0.5,0.975))

thetaestar1dagar = rbind(t(resar1dagar$betaest),resar1dagar$sigmaest,resar1dagar$rhoest,resar1dagar$phiest,rhophiest,resar1dagar$tau2est)
rownames(thetaestar1dagar) = names1
thetaestar1dagar[, c(1,2)] <- thetaestar1dagar[, c(2,1)]
colnames(thetaestar1dagar) = c( "50%","2.5%", "97.5%")


#### SAR ############


rhophi = resar1sar$rhoF*resar1sar$phiF
rhophiest = quantile(rhophi,c(0.025,0.5,0.975))


thetaestar1sar = rbind(t(resar1sar$betaest),resar1sar$sigmaest,resar1sar$rhoest,resar1sar$phiest,rhophiest,
                         resar1sar$tau2est)
rownames(thetaestar1sar) = names1
thetaestar1sar[, c(1,2)] <- thetaestar1sar[, c(2,1)]
colnames(thetaestar1sar) = c( "50%","2.5%", "97.5%")

#########################################
##### ESTIMATED MODELS #################
#######################################

### DAGAR model ######
round(thetaestar1dagar,3)

### SAR model ######
round(thetaestar1sar,3)

#####acceptance probability
###DAGAR ####
resar1dagar$probacc
####SAR#######
resar1sar$probacc

####TRUE VALUES

#sigma2 = 2
#beta = c(1,2,2.5)
#rho = 0.8 ## spatial parameter rho
#phi = 0.7## time parameter gamma
#tau2 = 0.6



