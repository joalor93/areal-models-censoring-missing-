library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)
library(invgamma)
library(rstan)
library(auxiliarcpp)
library(tmvtnorm)
library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)
library(data.table)
library(dplyr)
library(ggplot2)
library(auxiliarcpp)


#install.packages("/home/alejandro/postdoc2024/auxiliarcpp_0.0.0.9000.tar.gz", repos = NULL, type="source")


###################
library(parallel)
library(foreach)
library(doParallel)

##########################################
#### Avoiding issues in ubuntu###########
########################################

Sys.setenv(
  OMP_NUM_THREADS = "1",          # OpenMP-using code
  OPENBLAS_NUM_THREADS = "1",     # OpenBLAS (Ubuntu default)
  MKL_NUM_THREADS = "1",          # if you ever link MKL
  VECLIB_MAXIMUM_THREADS = "1"    # if on macOS Accelerate (not Ubuntu, harmless)
)

if ("RhpcBLASctl" %in% rownames(installed.packages())) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}
if ("data.table" %in% rownames(installed.packages())) data.table::setDTthreads(1)


n.cores <- 3 ##### number of cores (running three parallel chains)

cl <- parallel::makeCluster(n.cores, type = "PSOCK")


clusterCall(cl, Sys.setenv,
            OMP_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1",
            MKL_NUM_THREADS = "1",
            VECLIB_MAXIMUM_THREADS = "1"
)
clusterEvalQ(cl, {
  if ("RhpcBLASctl" %in% rownames(installed.packages())) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }
  if ("data.table" %in% rownames(installed.packages())) data.table::setDTthreads(1)
  NULL
})


doParallel::registerDoParallel(cl)




#### Reading the dataset ######

setwd("/home/alejandro/postdoc2024/censoring_dagar/Beijing")
data = fread(paste0("Estaciones_unidas_rushhour.csv"))
adj_matrix = readRDS("adj_matrix.rds")
source("auxiliar_functions.R")


#### complete adjacency matrix
adj_matcom = adj_matrix

datapre = data

regions = c("Aotizhongxin", "Changping", "Dingling", "Dongsi",
            "Guanyuan", "Gucheng", "Huairou", "Nongzhanguan",
            "Shunyi", "Tiantan","Wanliu", "Wanshouxigong")


#### Imputing the missing data with the total mean (we use them as initial values)
datapre2 = datapre%>%group_by("station")%>%
  mutate (aggtemp = if_else(is.na(TEMP),mean(TEMP,na.rm = TRUE),TEMP),
          aggrain = if_else(is.na(aggrain),mean(aggrain,na.rm = TRUE),aggrain),
          aggsumrain = if_else(is.na(aggsumrain),mean(aggsumrain,na.rm = TRUE),aggsumrain),
          aggpres = if_else(is.na(PRES),mean(PRES,na.rm = TRUE),PRES),
          aggwsmp = if_else(is.na(WSPM),mean(aggwsmp,na.rm = TRUE),WSPM))%>%ungroup()


######### Selecting the dataset to analyze
datapre2 = datapre2%>%filter(data>= "2016-11-25" & data<="2017-02-25")


###########################################################
##### Defining the input to estimate the DAGAR AR(1) model
##########################################################

ytotobs = datapre2$CO
indicens3 = as.numeric(is.na(ytotobs))
lower = rep(-Inf,sum(indicens3 ==1))
upper = rep(Inf,sum(indicens3 ==1))
ytotobs[indicens3 == 1] = mean(ytotobs,na.rm = TRUE)
xtotobs = cbind (1,datapre2$aggtemp - mean(datapre2$aggtemp),datapre2$aggwsmp - mean(datapre2$aggwsmp), (datapre2$aggpres - mean(datapre2$aggpres))/100)

thetaini = c(0.4,0.2,0.4,0.4)
iter = 40000
burn = 8000
thin = 80
N_time_obs = length(ytotobs[datapre2$station == "Shunyi"])
lags = N_time_obs
adjmatinf = adj_matrix
aprior = c(3.5,2,3.5,3.5)
bprior = c(2.5,4,2.5,2.5)

divproprho = 230
divproppsi = 230
divpropphi = 230


##### Creating the cluster
my.cluster <- makeCluster(n.cores)


n_sim =3 #### runing three chains in parallel

result <- foreach(j = 1:n_sim,.packages = c("Rcpp",
                                            "RcppArmadillo",
                                            "mvtnorm",
                                            "tmvtnorm",
                                            "parallel",
                                            "doParallel",
                                            "invgamma",
                                            "auxiliarcpp",
                                            "Matrix",
                                            "parallel",
                                            "doParallel",
                                            "foreach",
                                            "doSNOW","dplyr")
                                              ) %dopar% {


                                              set.seed(12345*j)
##### routine to run the DAGAR AR(1) model

           res2 = bayesspatempcensauto_fix(y = log(ytotobs),xtot = xtotobs,thetaini = thetaini[1:3],indicens = indicens3, iter = iter,burn = burn,thin = thin,lags = N_time_obs,
                              adjmatinf = adj_matrix, aprior = aprior[1:3], bprior = bprior[1:3],
                              divproprho = divproprho, divproppsi = divproppsi, divpropphi = divpropphi,lower = lower, upper = upper)

           saveRDS(res2, paste0("ar1dagar_beijing_212_40000_",j,iter,burn,thin,".rds"))
return(res2)

}

stopCluster(cl)
