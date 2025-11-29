loglike = function(listres,X, y){
 ## listres = thetaFval[[2]];X=xobs;
  p = ncol(X)
  N = length(y)
  theta = listres$thetaFval
  adjmatinf = listres$adjmatinf
  beta = theta[1:p]
  sigma2 = theta[p+1]
  rho = theta[p+2]
  psi = theta[p+3]
  mufix = X%*%beta
  varcova = sigma2*vdagar_cpp(adjmatinf,rho,psi)
  invvarcova = solve(varcova)
  sigmamu = invvarcova%*%(y-mufix)
  loglik = 0
  for(i in 1:N){
    sigmai = 1/sqrt(invvarcova[i,i])
    mui = y[i] - (sigmamu[i]*sigmai*sigmai)
    loglik[i] = dnorm(y[i],mui,sigmai,log = TRUE)
  }
  return(loglik)
}




loglikespatemp = function(listres,Xmat, y,lags){
  p = ncol(Xmat)
  N = length(y)
  theta = listres$thetaFval
  adjmatinf = listres$adjmatinf
  beta = theta[1:p]
  sigma2 = theta[p+1]
  rho = theta[p+2]
  psi = theta[p+3]
  phi = theta[p+4]
  mufix = Xmat%*%beta
  ##varcova = sigma2*vdagar_cpp(adjmatinf,rho,psi)
  varcova = sigma2*spatimecovar_2_cpp(lag = lags,adjmatinf = adjmatinf,rho = rho,
                               psi = psi,phi = phi)
  invvarcova = solve(varcova)
  sigmamu = invvarcova%*%(y-mufix)
  loglik = 0
  for(i in 1:N){
    sigmai = 1/sqrt(invvarcova[i,i])
    mui = y[i] - (sigmamu[i]*sigmai*sigmai)
    loglik[i] = dnorm(y[i],mui,sigmai,log = TRUE)
  }
  return(loglik)
}




loglikespatemp_fix = function(thetaFval,Xmat, y,lags,adjmatinf){
  #thetaFval1 =thetaFval[1,]
  #Xmat = xobs
  #y = yobs
  #lags = N_time_obs
  #adjmatinf = adjm1
  p = ncol(Xmat)
  N = length(y)
  theta = thetaFval
  beta = theta[1:p]
  sigma2 = theta[p+1]
  rho = theta[p+2]
  psi = theta[p+3]
  phi = theta[p+4]
  mufix = Xmat%*%beta
  ##varcova = sigma2*vdagar_cpp(adjmatinf,rho,psi)
  varcova = sigma2*spatimecovar_2_cpp(lag = lags,adjmatinf = adjmatinf,rho = rho,
                                      psi = psi,phi = phi)
  invvarcova = solve(varcova)
  sigmamu = invvarcova%*%(y-mufix)
  loglik = 0
  for(i in 1:N){
    #i = 1
    sigmai = 1/sqrt(invvarcova[i,i])
    mui = y[i] - (sigmamu[i]*sigmai*sigmai)
    loglik[i] = dnorm(y[i],mui,sigmai,log = TRUE)
  }
  return(loglik)
}


loglikespatempcens_fix = function(thetaFval,Xmat, ytotobs,lags,adjmatinf,indicens){
  #thetaFval1 =thetaFval[1,]
  #Xmat = xobs
  #y = yobs
  #lags = N_time_obs
  #adjmatinf = adjm1
  p = ncol(Xmat)
  N = length(ytotobs)
  theta = thetaFval
  beta = theta[1:p]
  sigma2 = theta[p+1]
  rho = theta[p+2]
  psi = theta[p+3]
  phi = theta[p+4]
  ycensF = theta[(p+5):length(theta)]
  ytotobs[indicens ==1] = ycensF

  mufix = Xmat%*%beta
  ##varcova = sigma2*vdagar_cpp(adjmatinf,rho,psi)
  varcova = sigma2*spatimecovar_2_cpp(lag = lags,adjmatinf = adjmatinf,rho = rho,
                                      psi = psi,phi = phi)
  invvarcova = solve(varcova)
  sigmamu = invvarcova%*%(ytotobs-mufix)
  loglik = 0
  for(i in 1:N){
    #i = 1
    sigmai = 1/sqrt(invvarcova[i,i])
    mui = ytotobs[i] - (sigmamu[i]*sigmai*sigmai)
    loglik[i] = dnorm(ytotobs[i],mui,sigmai,log = TRUE)
  }
  return(loglik)
}

#' @export
spatimecovarcar_2 = function(lag,adjmat,rho,psi,phi){
matrix1 = spatimecovarcar_2_cpp(lag = lag,adjmatinf = adjmat,rho = rho,psi = psi,
                     phi = phi)
return(matrix1)
}

#' @export
spatimecovarsar_2 = function(lag,adjmatinf,rho,psi,phi){
  matrix1 = spatimecovarsar_2_cpp(lag = lag,adjmat = adjmatinf,rho = rho,psi = psi,
                                  phi = phi)
  return(matrix1)
}

#' @export
spatimecovarsar_21 = function(lag,adjmatinf,rho,psi,phi){
  matrix1 = spatimecovarsar_21_cpp(lag = lag,adjmatinf = adjmatinf,rho = rho,psi = psi,
                                  phi = phi)
  return(matrix1)
}


#' @export
spatimecovarsar_22 = function(lag,adjmatinf,rho,psi,phi){
  matrix1 = spatimecovarsar_21_cpp(lag = lag,adjmatinf = adjmatinf,rho = rho,psi = psi,
                                   phi = phi)
  return(matrix1)
}


#' @export
varcovspacar_2 = function(adjmat,rho){
  matrix1 = varcovspacar_cpp(adjmatinf = adjmat, rho = rho)
  return(matrix1)
}





#' @export
bayestrandgrap = function(y,xobs,thetaini,iter,burn,thin,
                          weight_mat, aprior = c(4,2), bprior = c(2,4),
                          divproprho = 18, divproppsi = 20){

  #adjsamini = primSpanningTree(n)
  #adjsaminilist = generate_weighted_spanning_tree_aldous_broder_optimized_2(weight_mat)
  #adjsamini = as.matrix(adjsaminilist$adjacency_matrix)

  #y= y; xobs = xobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;
  #weight_mat = weight_mat;aprior = c(4,2); bprior = c(2,4);
  #divproprho = 9; divproppsi = 9

  prob_mat = weight_mat/rowSums(weight_mat)
  start = which(prob_mat == max(prob_mat),arr.ind = TRUE)
  start = start[1,1]
  adjsamini = fast_cover(weight_mat,start = start,threshold = 1000)
  rhoF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  psiF = rep(0,iter)
  WF = list()
  thetaFval = list()
  rhoF[1] = thetaini[1]
  psiF[1] = thetaini[2]
  WF[[1]] = adjsamini
  logpw = rep(0,iter)
  #logpw[[1]] =  adjsaminilist$log_probability
  proba = prob_mat[adjsamini ==1]
  logpw[1] = sum(log(proba[proba!=0]))
  ##WF[[1]][upper.tri(WF[[1]])]=0

  for(i in 2:iter){
    #i = 2
    #adjsamcandauxlist = generate_weighted_spanning_tree_aldous_broder_optimized_2(weight_mat)
    adjsamcandaux = fast_cover(weight_mat,start = start,threshold = 1000)
    #adjsamcandaux = as.matrix(adjsamcandauxlist$adjacency_matrix)
    adjsamlastaux = WF[[i-1]]

    adjsamcand = adjsamcandaux
    adjsamlast = adjsamlastaux

    probacand = prob_mat[adjsamcand ==1]
    logpwcand = sum(log(probacand[probacand!=0]))

    #logpwcand =   adjsamcandauxlist$log_probability
    logpwlast = logpw[[i-1]]

    ##logpriorwcan = sum(log(weight_mat[adjsamcand == 1]))
    ##logpriorlast = sum(log(weight_mat[adjsamcand == 1]))
    #adjsamcand = primSpanningTree(n)
    #adjsamlast = WF[[i-1]]

    ##adjsamcand = adjsamcandaux
    ##adjsamlast = adjsamlastaux

    adjsamcand[upper.tri(adjsamcand)]=0
    adjsamlast[upper.tri(adjsamlast)]=0

    thetalast = c(rhoF[i-1],psiF[i-1])
    cons = 1
    muprop1 = min(0.99,thetalast[1]/cons)
    muprop2 = min(0.99,thetalast[2]/cons)
    muprop = c(muprop1,muprop2)
    sigmaprop1 = muprop[1]*(1-muprop[1])/divproprho
    sigmaprop2 = muprop[2]*(1-muprop[2])/divproppsi

    #sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop = c(sigmaprop1,sigmaprop2)
    omega = runif(1,0.99,1)
    theta1cand = min(0.99,rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    theta2cand = min(0.99,rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    thetacand = c(theta1cand,theta2cand)
    #thetacand

    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))

    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))

    lognum = log(posteriorvar3_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y, a = aprior, b = bprior)) + logpwcand
    logden = log(posteriorvar3_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y, a = aprior, b= bprior)) + logpwlast


    lognum = lognum + log(dbeta_rep_cpp(thetalast[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetalast[2],omega*muprop[2],sigmaprop[2])) + logpwlast
    logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) + logpwcand



    #acand =
    #num = num*dbeta(thetalast[1],5,2)*dbeta(thetalast[2],2,5)
    #den = den*dbeta(thetacand[1],5,2)*dbeta(thetacand[2],2,5)


    logunif= log(runif(1))
    #numden[i] = num/den

    if(logunif<lognum - logden)
    {
      count = count +1
      rhoF[i] <- thetacand[1]
      psiF[i] = thetacand[2]
      WF[[i]] = adjsamcandaux
      logpw[[i]] = logpwcand

    }else
    {
      rhoF[i] <- thetalast[1]
      psiF[i] = thetalast[2]
      WF[[i]] = adjsamlastaux
      logpw[[i]] = logpwlast
    }
   print(c(i,rhoF[i],psiF[i],count/i))
    #gc()
  }

  gc()

  rhoburn=rhoF[(burn+1):iter]
  rhoval= rhoburn[seq((burn+1),iter-burn,thin)]

  psiburn=psiF[(burn+1):iter]
  psival= psiburn[seq((burn+1),iter-burn,thin)]


  WFburn=WF[(burn+1):iter]
  WFval= WFburn[seq((burn+1),iter-burn,thin)]

  logpwburn=WF[(burn+1):iter]
  logpwval= WFburn[seq((burn+1),iter-burn,thin)]


  probacc = sum(count)/iter



  betaFval = matrix(0,length(rhoval),ncol(xobs))
  sigmaFval = 0
  thetaFval = list()
  n = length(y)
  p = ncol(xobs)
  for (i in 1:length(rhoval)){
    adjmatinf = WFval[[i]]
    adjmatinf[upper.tri(adjmatinf)] = 0
    rhovalnum = rhoval[i]
    psivalnum = psival[[i]]
    vdag = vdagar_cpp(adjmatinf = adjmatinf,rho = rhovalnum,psi = psivalnum)
    vdaginv = solve(vdag)
    vbetaest = solve(t(xobs)%*%vdaginv%*%xobs)
    betaest = vbetaest%*%t(xobs)%*%vdaginv%*%y
    S2 = t(y-xobs%*%betaest)%*%vdaginv%*%(y-xobs%*%betaest)
    betaFval[i,] = rmvnorm(1,betaest,vbetaest)
    sigmaFval[i] = rinvgamma(1,0.5*(n-p), S2/2)
    thetaval = c(betaFval[i,],sigmaFval[i],rhovalnum,psivalnum)
    thetaFval[[i]] = list(thetaFval =thetaval,adjmatinf = adjmatinf)
  }

  log_lik = rbind()
  for(i in 1:length(rhoval)){
  log_lik = rbind(log_lik,loglike(thetaFval[[i]],xobs,y))
  }

    #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  rhoest = quantile(rhoval,probs = c(0.025,0.5,0.975))
  psiest = quantile(psival,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigmaFval,probs = c(0.025,0.5,0.975))

  result = list(betaF = betaFval, sigmaF = sigmaFval, rhoF = rhoval, psiF = psival,
                thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
                psiest = psiest,probacc = probacc,log_lik = log_lik)

  return(result)
}




#' @export
bayesspatemp = function(y,xobs,thetaini,iter,burn,thin, lags,
                          weight_mat, aprior = c(4,2,4), bprior = c(2,4,2),
                          divproprho = 9, divproppsi = 9, divpropphi = 9,
                        apriorprop = c(4,2,4),bpriorprop = c(2,4,2)){

  #adjsamini = primSpanningTree(n)
  #adjsaminilist = generate_weighted_spanning_tree_aldous_broder_optimized_2(weight_mat)
  #adjsamini = as.matrix(adjsaminilist$adjacency_matrix)

 # y= y; xobs = xobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;
  #weight_mat = weight_mat;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time

  prob_mat = weight_mat/rowSums(weight_mat)
  start = which(prob_mat == max(prob_mat),arr.ind = TRUE)
  start = start[1,1]
  adjsamini = fast_cover(weight_mat,start = start,threshold = 1000)
  rhoF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  psiF = rep(0,iter)
  phiF = rep(0,iter)
  WF = list()
  thetaFval = list()
  rhoF[1] = thetaini[1]
  psiF[1] = thetaini[2]
  phiF[1] = thetaini[3]
  WF[[1]] = adjsamini
  logpw = rep(0,iter)
  #logpw[[1]] =  adjsaminilist$log_probability
  proba = prob_mat[adjsamini ==1]
  logpw[1] = sum(log(proba[proba!=0]))
  ##WF[[1]][upper.tri(WF[[1]])]=0

  for(i in 2:iter){
   # i = 2
    #adjsamcandauxlist = generate_weighted_spanning_tree_aldous_broder_optimized_2(weight_mat)
    adjsamcandaux = fast_cover(weight_mat,start = start,threshold = 1000)
    #adjsamcandaux = as.matrix(adjsamcandauxlist$adjacency_matrix)
    adjsamlastaux = WF[[i-1]]

    adjsamcand = adjsamcandaux
    adjsamlast = adjsamlastaux

    probacand = prob_mat[adjsamcand ==1]
    logpwcand = sum(log(probacand[probacand!=0]))

    #logpwcand =   adjsamcandauxlist$log_probability
    logpwlast = logpw[[i-1]]

    ##logpriorwcan = sum(log(weight_mat[adjsamcand == 1]))
    ##logpriorlast = sum(log(weight_mat[adjsamcand == 1]))
    #adjsamcand = primSpanningTree(n)
    #adjsamlast = WF[[i-1]]

    ##adjsamcand = adjsamcandaux
    ##adjsamlast = adjsamlastaux

    adjsamcand[upper.tri(adjsamcand)]=0
    adjsamlast[upper.tri(adjsamlast)]=0

    thetalast = c(rhoF[i-1],psiF[i-1],phiF[i-1])
   # cons = 1
    #muprop1 = min(0.99,thetalast[1]/cons)
    #muprop2 = min(0.99,thetalast[2]/cons)
    #muprop3 = min(0.99,thetalast[3]/cons)

    #muprop = c(muprop1,muprop2,muprop3)
    #sigmaprop1 = muprop[1]*(1-muprop[1])/divproprho
    #sigmaprop2 = muprop[2]*(1-muprop[2])/divproppsi
    #sigmaprop3 = muprop[3]*(1-muprop[3])/divpropphi

    #sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    #sigmaprop = c(sigmaprop1,sigmaprop2,sigmaprop3)
    #omega = runif(1,0.99,1)
  #  theta1cand = min(0.99,rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
   # theta2cand = min(0.99,rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
  #  theta3cand = min(0.99,rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))


   #apriorprop = aprior + 0.5
   #bpriorprop = bprior + 0.5
   # thetacand = c(theta1cand,theta2cand,theta3cand)
    thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))
    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))
   # a = posteriorvar3_cpp(thetacand,adjmatinf = adj_matrix,lags,X,ytot,a=c(1,1,1),b =c(1,1,1))
    # posterior_spatemp_cpp(thetacand,adjmatinf = adj_matrix,lag = N_time,X=xobs,y=y,a=c(1,1,1),b =c(1,1,1))
    #posterior_spatempver_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))

    lognum = posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior) + logpwcand
    logden = posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjsamlast,lag = lags, X = xobs,y = y, a = aprior, b= bprior) + logpwlast


   ## lognum = lognum + log(dbeta_rep_cpp(thetalast[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetalast[2],omega*muprop[2],sigmaprop[2])) +
     ##       log(dbeta_rep_cpp(thetalast[3],omega*muprop[3],sigmaprop[3])) ##+ logpwlast

    ##logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
      ##       log(dbeta_rep_cpp(thetacand[3],omega*muprop[3],sigmaprop[3])) ##+ logpwcand

     lognum = lognum + log(dbeta(thetalast[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetalast[2],apriorprop[2],bpriorprop[2])) +
         log(dbeta(thetalast[3],apriorprop[3],bpriorprop[3])) + logpwlast

    logden = logden + log(dbeta(thetacand[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetacand[2],apriorprop[2],bpriorprop[2])) +
      log(dbeta(thetacand[3],apriorprop[3],bpriorprop[3])) + logpwcand


    #acand =
    #num = num*dbeta(thetalast[1],5,2)*dbeta(thetalast[2],2,5)
    #den = den*dbeta(thetacand[1],5,2)*dbeta(thetacand[2],2,5)


    logunif= log(runif(1))
    #numden[i] = num/den

    if(logunif<lognum - logden)
    {
      count = count +1
      rhoF[i] <- thetacand[1]
      psiF[i] = thetacand[2]
      phiF[i] = thetacand[3]
      WF[[i]] = adjsamcandaux
     # logpw[[i]] = logpwcand

    }else
    {
      rhoF[i] <- thetalast[1]
      psiF[i] = thetalast[2]
      phiF[i] = thetalast[3]
      WF[[i]] = adjsamlastaux
      #logpw[[i]] = logpwlast
    }
    print(c(i,rhoF[i],psiF[i],phiF[i],count/i))
    #gc()
  }

  gc()

  rhoburn=rhoF[(burn+1):iter]
  rhoval= rhoburn[seq((burn+1),iter-burn,thin)]

  psiburn=psiF[(burn+1):iter]
  psival= psiburn[seq((burn+1),iter-burn,thin)]

  phiburn=phiF[(burn+1):iter]
  phival= phiburn[seq((burn+1),iter-burn,thin)]


  WFburn=WF[(burn+1):iter]
  WFval= WFburn[seq((burn+1),iter-burn,thin)]

  #logpwburn=WF[(burn+1):iter]
  #logpwval= WFburn[seq((burn+1),iter-burn,thin)]


  probacc = sum(count)/iter



  betaFval = matrix(0,length(rhoval),ncol(xobs))
  sigmaFval = 0
  thetaFval = list()
  n = length(y)
  p = ncol(xobs)
  for (i in 1:length(rhoval)){
    #i = 1
    adjmatinf = WFval[[i]]
    adjmatinf[upper.tri(adjmatinf)] = 0
    rhovalnum = rhoval[i]
    psivalnum = psival[[i]]
    phivalnum = phival[[i]]
    vdag =spatimecovar_2_cpp(lag = lags,adjmatinf = adjmatinf,rho = rhovalnum,
                       psi = psivalnum,phi = phivalnum)
    vdaginv = solve(vdag)
    vbetaest = solve(t(xobs)%*%vdaginv%*%xobs)
    betaest = vbetaest%*%t(xobs)%*%vdaginv%*%y
    S2 = t(y-xobs%*%betaest)%*%vdaginv%*%(y-xobs%*%betaest)
    betaFval[i,] = rmvnorm(1,betaest,vbetaest)
    sigmaFval[i] = rinvgamma(1,0.5*(n-p), S2/2)
    thetaval = c(betaFval[i,],sigmaFval[i],rhovalnum,psivalnum,phivalnum)
    thetaFval[[i]] = list(thetaFval =thetaval,adjmatinf = adjmatinf)
  }

  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  rhoest = quantile(rhoval,probs = c(0.025,0.5,0.975))
  psiest = quantile(psival,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigmaFval,probs = c(0.025,0.5,0.975))

  result = list(betaF = betaFval, sigmaF = sigmaFval, rhoF = rhoval, psiF = psival,
                thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
                psiest = psiest, probacc = probacc,log_lik = log_lik)

  return(result)
}



#' @export
bayesspatempauto = function(y,xobs,thetaini,iter,burn,thin, lags,
                        weight_mat, aprior = c(4,2,4), bprior = c(2,4,2),
                        divproprho = 9, divproppsi = 9, divpropphi = 9){

  #adjsamini = primSpanningTree(n)
  #adjsaminilist = generate_weighted_spanning_tree_aldous_broder_optimized_2(weight_mat)
  #adjsamini = as.matrix(adjsaminilist$adjacency_matrix)

  # y= y; xobs = xobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;
  #weight_mat = weight_mat;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time

  prob_mat = weight_mat/rowSums(weight_mat)
  start = which(prob_mat == max(prob_mat),arr.ind = TRUE)
  start = start[1,1]
  adjsamini = fast_cover(weight_mat,start = start,threshold = 1000)
  rhoF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  psiF = rep(0,iter)
  phiF = rep(0,iter)
  WF = list()
  thetaFval = list()
  rhoF[1] = thetaini[1]
  psiF[1] = thetaini[2]
  phiF[1] = thetaini[3]
  WF[[1]] = adjsamini
  logpw = rep(0,iter)
  #logpw[[1]] =  adjsaminilist$log_probability
  proba = prob_mat[adjsamini ==1]
  logpw[1] = sum(log(proba[proba!=0]))
  ##WF[[1]][upper.tri(WF[[1]])]=0

  for(i in 2:iter){
    # i = 2
    #adjsamcandauxlist = generate_weighted_spanning_tree_aldous_broder_optimized_2(weight_mat)
    adjsamcandaux = fast_cover(weight_mat,start = start,threshold = 1000)
    #adjsamcandaux = as.matrix(adjsamcandauxlist$adjacency_matrix)
    adjsamlastaux = WF[[i-1]]

    adjsamcand = adjsamcandaux
    adjsamlast = adjsamlastaux

    probacand = prob_mat[adjsamcand ==1]
    logpwcand = sum(log(probacand[probacand!=0]))

    #logpwcand =   adjsamcandauxlist$log_probability
    logpwlast = logpw[[i-1]]

    ##logpriorwcan = sum(log(weight_mat[adjsamcand == 1]))
    ##logpriorlast = sum(log(weight_mat[adjsamcand == 1]))
    #adjsamcand = primSpanningTree(n)
    #adjsamlast = WF[[i-1]]

    ##adjsamcand = adjsamcandaux
    ##adjsamlast = adjsamlastaux

    adjsamcand[upper.tri(adjsamcand)]=0
    adjsamlast[upper.tri(adjsamlast)]=0

    thetalast = c(rhoF[i-1],psiF[i-1],phiF[i-1])
    cons = 1
    muprop1 = min(0.99,thetalast[1]/cons)
    muprop2 = min(0.99,thetalast[2]/cons)
    muprop3 = min(0.99,thetalast[3]/cons)

    muprop = c(muprop1,muprop2,muprop3)
    sigmaprop1 = muprop[1]*(1-muprop[1])/divproprho
    sigmaprop2 = muprop[2]*(1-muprop[2])/divproppsi
    sigmaprop3 = muprop[3]*(1-muprop[3])/divpropphi

#    sigmaprop1 = muprop[1]*(1-muprop[1])/50
 #   sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop = c(sigmaprop1,sigmaprop2,sigmaprop3)
    omega = runif(1,0.99,1)
    theta1cand = min(0.99,rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    theta2cand = min(0.99,rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    theta3cand = min(0.99,rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))

    #if(theta1cand> 0.95){
      #theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
     # theta1cand = runif(1,0.1,0.6)
    #}

    if(theta2cand > 0.6){
      theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
      theta2cand = runif(1,0.1,0.6)
    }


    #if(theta3cand > 0.95){
      #theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
     # theta3cand = runif(1,0.1,0.6)
    #}



    if(theta1cand<0.02){
      #theta1cand =max(runif(1,0.1,0.2),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
      theta1cand = runif(1,0.1,0.2)
    }

    if(theta2cand<0.02){
      #theta2cand =max(runif(1,0.1,0.2),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
      theta2cand = runif(1,0.1,0.2)
    }


    if(theta3cand<0.02){
      #theta3cand =max(runif(1,0.1,0.2),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
      theta3cand = runif(1,0.1,0.2)
    }


    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
     thetacand = c(theta1cand,theta2cand,theta3cand)
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))
    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))
    # a = posteriorvar3_cpp(thetacand,adjmatinf = adj_matrix,lags,X,ytot,a=c(1,1,1),b =c(1,1,1))
    # posterior_spatemp_cpp(thetacand,adjmatinf = adj_matrix,lag = N_time,X=xobs,y=y,a=c(1,1,1),b =c(1,1,1))
    #posterior_spatempver_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    lognum = posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    logden = posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjsamlast,lag = lags, X = xobs,y = y, a = aprior, b= bprior)

    #if(is.na(num)){
     # num = runif(1)
    #}

    #if(is.na(den)){
     # den = runif(1)
    #}

     lognum = lognum #+ logpwcand
     logden = logden #+ logpwlast
    #lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)) + logpwcand
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjsamlast,lag = lags, X = xobs,y = y, a = aprior, b= bprior)) + logpwlast


     lognum = lognum + log(dbeta_rep_cpp(thetalast[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetalast[2],omega*muprop[2],sigmaprop[2])) +
           log(dbeta_rep_cpp(thetalast[3],omega*muprop[3],sigmaprop[3])) #+ logpwlast

    logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
           log(dbeta_rep_cpp(thetacand[3],omega*muprop[3],sigmaprop[3])) #+ logpwcand

   # lognum = lognum + log(dbeta(thetalast[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetalast[2],apriorprop[2],bpriorprop[2])) +
    #  log(dbeta(thetalast[3],apriorprop[3],bpriorprop[3])) + logpwlast

    #logden = logden + log(dbeta(thetacand[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetacand[2],apriorprop[2],bpriorprop[2])) +
     # log(dbeta(thetacand[3],apriorprop[3],bpriorprop[3])) + logpwcand


    #acand =
    #num = num*dbeta(thetalast[1],5,2)*dbeta(thetalast[2],2,5)
    #den = den*dbeta(thetacand[1],5,2)*dbeta(thetacand[2],2,5)


    logunif= log(runif(1))
    #numden[i] = num/den

    if(logunif<lognum - logden)
    {
      count = count +1
      rhoF[i] <- thetacand[1]
      psiF[i] = thetacand[2]
      phiF[i] = thetacand[3]
      WF[[i]] = adjsamcandaux
      # logpw[[i]] = logpwcand

    }else
    {
      rhoF[i] <- thetalast[1]
      psiF[i] = thetalast[2]
      phiF[i] = thetalast[3]
      WF[[i]] = adjsamlastaux
      #logpw[[i]] = logpwlast
    }
    print(c(i,rhoF[i],psiF[i],phiF[i],count/i))
    #gc()
  }

  gc()

  rhoburn=rhoF[(burn+1):iter]
  rhoval= rhoburn[seq((burn+1),iter-burn,thin)]

  psiburn=psiF[(burn+1):iter]
  psival= psiburn[seq((burn+1),iter-burn,thin)]

  phiburn=phiF[(burn+1):iter]
  phival= phiburn[seq((burn+1),iter-burn,thin)]


  WFburn=WF[(burn+1):iter]
  WFval= WFburn[seq((burn+1),iter-burn,thin)]

  #logpwburn=WF[(burn+1):iter]
  #logpwval= WFburn[seq((burn+1),iter-burn,thin)]


  probacc = sum(count)/iter



  betaFval = matrix(0,length(rhoval),ncol(xobs))
  sigmaFval = 0
  thetaFval = list()
  n = length(y)
  p = ncol(xobs)
  for (i in 1:length(rhoval)){
    #i = 1
    adjmatinf = WFval[[i]]
    adjmatinf[upper.tri(adjmatinf)] = 0
    rhovalnum = rhoval[i]
    psivalnum = psival[[i]]
    phivalnum = phival[[i]]
    vdag =spatimecovar_2_cpp(lag = lags,adjmatinf = adjmatinf,rho = rhovalnum,
                             psi = psivalnum,phi = phivalnum)
    vdaginv = solve(vdag)
    vbetaest = solve(t(xobs)%*%vdaginv%*%xobs)
    betaest = vbetaest%*%t(xobs)%*%vdaginv%*%y
    S2 = t(y-xobs%*%betaest)%*%vdaginv%*%(y-xobs%*%betaest)
    sigmaFval[i] = rinvgamma(1,0.5*(n-p), S2/2)
    betaFval[i,] = rmvnorm(1,betaest, sigmaFval[i]*vbetaest)
    thetaval = c(betaFval[i,],sigmaFval[i],rhovalnum,psivalnum,phivalnum)
    thetaFval[[i]] = list(thetaFval =thetaval,adjmatinf = adjmatinf)
  }

  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  rhoest = quantile(rhoval,probs = c(0.025,0.5,0.975))
  psiest = quantile(psival,probs = c(0.025,0.5,0.975))
  phiest = quantile(phival,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigmaFval,probs = c(0.025,0.5,0.975))

  result = list(betaF = betaFval, sigmaF = sigmaFval, rhoF = rhoval, psiF = psival,
                thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
                psiest = psiest,phiest = phiest, probacc = probacc,log_lik = log_lik)

  return(result)
}


#' @export
generating_tree_subtrees = function(subtrees){

  ############# getting a tree from the subtreees ########
 # num_vertices <- sapply(subtrees, vcount)
  #w <- order(num_vertices, decreasing = TRUE)
  w = sample(1:length(subtrees),length(subtrees))
  # OPTION 2: Or use a function to pick nodes, e.g., highest degree node
  vertex = rbind()
  for(i in 2:length(subtrees)){
    # set.seed(rnorm(1))
    node_a <- names(sample(degree(subtrees[[w[i-1]]]),1))
    node_b <- names(sample(degree(subtrees[[w[i]]]),1))
    vertex = rbind(vertex,c(node_a,node_b))

  }




  # Merge all subtrees into one disconnected graph
  forest <- subtrees[[1]]
  for (i in 2:length(subtrees)) {
    forest <- forest %u% subtrees[[i]]
  }


  forest <- add_edges(forest, vertex)


  return(forest)

}


#' @export
bayesspatempauto_2 = function(y,xobs,thetaini,iter,burn,thin, lags,
                            subtrees, aprior = c(4,2,4), bprior = c(2,4,2),
                            divproprho = 9, divproppsi = 9, divpropphi = 9){

  #adjsamini = primSpanningTree(n)
  #adjsaminilist = generate_weighted_spanning_tree_aldous_broder_optimized_2(weight_mat)
  #adjsamini = as.matrix(adjsaminilist$adjacency_matrix)

  # y= y; xobs = xobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;
  #weight_mat = weight_mat;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time

  #prob_mat = weight_mat/rowSums(weight_mat)
  #start = which(prob_mat == max(prob_mat),arr.ind = TRUE)
  #start = start[1,1]
  #adjsamini = fast_cover(weight_mat,start = start,threshold = 1000)
  treerand = generating_tree_subtrees(subtrees)
  adjsamini = as.matrix(as_adjacency_matrix(treerand))
  rhoF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  psiF = rep(0,iter)
  phiF = rep(0,iter)
  WF = list()
  thetaFval = list()
  rhoF[1] = thetaini[1]
  psiF[1] = thetaini[2]
  phiF[1] = thetaini[3]
  WF[[1]] = adjsamini
 # logpw = rep(0,iter)
  #logpw[[1]] =  adjsaminilist$log_probability
  #proba = prob_mat[adjsamini ==1]
  #logpw[1] = sum(log(proba[proba!=0]))
  ##WF[[1]][upper.tri(WF[[1]])]=0

  for(i in 2:iter){
    # i = 2
    #adjsamcandauxlist = generate_weighted_spanning_tree_aldous_broder_optimized_2(weight_mat)
    #adjsamcandaux = fast_cover(weight_mat,start = start,threshold = 1000)
    treerand1 = generating_tree_subtrees(subtrees)
    adjsamcandaux = as.matrix(as_adjacency_matrix(treerand1))

    #adjsamcandaux = as.matrix(adjsamcandauxlist$adjacency_matrix)
    adjsamlastaux = WF[[i-1]]

    adjsamcand = adjsamcandaux
    adjsamlast = adjsamlastaux

    #probacand = prob_mat[adjsamcand ==1]
    #logpwcand = sum(log(probacand[probacand!=0]))

    #logpwcand =   adjsamcandauxlist$log_probability
    #logpwlast = logpw[[i-1]]

    ##logpriorwcan = sum(log(weight_mat[adjsamcand == 1]))
    ##logpriorlast = sum(log(weight_mat[adjsamcand == 1]))
    #adjsamcand = primSpanningTree(n)
    #adjsamlast = WF[[i-1]]

    ##adjsamcand = adjsamcandaux
    ##adjsamlast = adjsamlastaux

    adjsamcand[upper.tri(adjsamcand)]=0
    adjsamlast[upper.tri(adjsamlast)]=0

    thetalast = c(rhoF[i-1],psiF[i-1],phiF[i-1])
    cons = 1
    muprop1 = min(0.99,thetalast[1]/cons)
    muprop2 = min(0.99,thetalast[2]/cons)
    muprop3 = min(0.99,thetalast[3]/cons)

    muprop = c(muprop1,muprop2,muprop3)
    sigmaprop1 = muprop[1]*(1-muprop[1])/divproprho
    sigmaprop2 = muprop[2]*(1-muprop[2])/divproppsi
    sigmaprop3 = muprop[3]*(1-muprop[3])/divpropphi

    #    sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #   sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop = c(sigmaprop1,sigmaprop2,sigmaprop3)
    omega = runif(1,0.99,1)
    theta1cand = min(0.99,rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    theta2cand = min(0.99,rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    theta3cand = min(0.99,rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))

    #if(theta1cand> 0.95){
    #theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    # theta1cand = runif(1,0.1,0.6)
    #}

    if(theta2cand > 0.6){
      theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
      theta2cand = runif(1,0.1,0.6)
    }


    #if(theta3cand > 0.95){
    #theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
    # theta3cand = runif(1,0.1,0.6)
    #}



    if(theta1cand<0.08){
      #theta1cand =max(runif(1,0.1,0.2),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
      theta1cand = runif(1,0.1,0.2)
    }

    if(theta2cand<0.08){
      #theta2cand =max(runif(1,0.1,0.2),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
      theta2cand = runif(1,0.1,0.2)
    }


    if(theta3cand<0.08){
      #theta3cand =max(runif(1,0.1,0.2),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
      theta3cand = runif(1,0.1,0.2)
    }


    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    thetacand = c(theta1cand,theta2cand,theta3cand)
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))
    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))
    # a = posteriorvar3_cpp(thetacand,adjmatinf = adj_matrix,lags,X,ytot,a=c(1,1,1),b =c(1,1,1))
    # posterior_spatemp_cpp(thetacand,adjmatinf = adj_matrix,lag = N_time,X=xobs,y=y,a=c(1,1,1),b =c(1,1,1))
    #posterior_spatempver_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    lognum = posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    logden = posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjsamlast,lag = lags, X = xobs,y = y, a = aprior, b= bprior)

    #if(is.na(num)){
    # num = runif(1)
    #}

    #if(is.na(den)){
    # den = runif(1)
    #}

    lognum = lognum #+ logpwcand
    logden = logden #+ logpwlast
    #lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)) + logpwcand
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjsamlast,lag = lags, X = xobs,y = y, a = aprior, b= bprior)) + logpwlast


    lognum = lognum + log(dbeta_rep_cpp(thetalast[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetalast[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_cpp(thetalast[3],omega*muprop[3],sigmaprop[3])) #+ logpwlast

    logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_cpp(thetacand[3],omega*muprop[3],sigmaprop[3])) #+ logpwcand

    # lognum = lognum + log(dbeta(thetalast[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetalast[2],apriorprop[2],bpriorprop[2])) +
    #  log(dbeta(thetalast[3],apriorprop[3],bpriorprop[3])) + logpwlast

    #logden = logden + log(dbeta(thetacand[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetacand[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetacand[3],apriorprop[3],bpriorprop[3])) + logpwcand


    #acand =
    #num = num*dbeta(thetalast[1],5,2)*dbeta(thetalast[2],2,5)
    #den = den*dbeta(thetacand[1],5,2)*dbeta(thetacand[2],2,5)


    logunif= log(runif(1))
    #numden[i] = num/den

    if(logunif<lognum - logden)
    {
      count = count +1
      rhoF[i] <- thetacand[1]
      psiF[i] = thetacand[2]
      phiF[i] = thetacand[3]
      WF[[i]] = adjsamcandaux
      # logpw[[i]] = logpwcand

    }else
    {
      rhoF[i] <- thetalast[1]
      psiF[i] = thetalast[2]
      phiF[i] = thetalast[3]
      WF[[i]] = adjsamlastaux
      #logpw[[i]] = logpwlast
    }
    print(c(i,rhoF[i],psiF[i],phiF[i],count/i))
    #gc()
  }

  gc()

  rhoburn=rhoF[(burn+1):iter]
  rhoval= rhoburn[seq((burn+1),iter-burn,thin)]

  psiburn=psiF[(burn+1):iter]
  psival= psiburn[seq((burn+1),iter-burn,thin)]

  phiburn=phiF[(burn+1):iter]
  phival= phiburn[seq((burn+1),iter-burn,thin)]


  WFburn=WF[(burn+1):iter]
  WFval= WFburn[seq((burn+1),iter-burn,thin)]

  #logpwburn=WF[(burn+1):iter]
  #logpwval= WFburn[seq((burn+1),iter-burn,thin)]


  probacc = sum(count)/iter



  betaFval = matrix(0,length(rhoval),ncol(xobs))
  sigmaFval = 0
  thetaFval = list()
  n = length(y)
  p = ncol(xobs)
  for (i in 1:length(rhoval)){
    #i = 1
    adjmatinf = WFval[[i]]
    adjmatinf[upper.tri(adjmatinf)] = 0
    rhovalnum = rhoval[i]
    psivalnum = psival[[i]]
    phivalnum = phival[[i]]
    vdag =spatimecovar_2_cpp(lag = lags,adjmatinf = adjmatinf,rho = rhovalnum,
                             psi = psivalnum,phi = phivalnum)
    vdaginv = solve(vdag)
    vbetaest = solve(t(xobs)%*%vdaginv%*%xobs)
    betaest = vbetaest%*%t(xobs)%*%vdaginv%*%y
    S2 = t(y-xobs%*%betaest)%*%vdaginv%*%(y-xobs%*%betaest)
    sigmaFval[i] = rinvgamma(1,0.5*(n-p), S2/2)
    betaFval[i,] = rmvnorm(1,betaest, sigmaFval[i]*vbetaest)
    thetaval = c(betaFval[i,],sigmaFval[i],rhovalnum,psivalnum,phivalnum)
    thetaFval[[i]] = list(thetaFval =thetaval,adjmatinf = adjmatinf)
  }

  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  rhoest = quantile(rhoval,probs = c(0.025,0.5,0.975))
  psiest = quantile(psival,probs = c(0.025,0.5,0.975))
  phiest = quantile(phival,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigmaFval,probs = c(0.025,0.5,0.975))

  result = list(betaF = betaFval, sigmaF = sigmaFval, rhoF = rhoval, psiF = psival,
                thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
                psiest = psiest,phiest = phiest, probacc = probacc,log_lik = log_lik)

  return(result)
}












#' @export
bayesspatemp_fix = function(y,xobs,thetaini,iter,burn,thin, lags,
                        adjmatinf, aprior = c(4,2,4), bprior = c(2,4,2),
                        divproprho = 9, divproppsi = 9, divpropphi = 9){

  # y= yobs; xobs = xobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;
  #adjmatinf = adjm1;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time_obs





  rhoF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  psiF = rep(0,iter)
  phiF = rep(0,iter)
  rhoF[1] = thetaini[1]
  psiF[1] = thetaini[2]
  phiF[1] = thetaini[3]

  for(i in 2:iter){
    # i = 2
    adjmatinf[upper.tri(adjmatinf)]=0

    thetalast = c(rhoF[i-1],psiF[i-1],phiF[i-1])

    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5

    cons = 1
    muprop1 = min(0.99,thetalast[1]/cons)
    muprop2 = min(0.99,thetalast[2]/cons)
    muprop3 = min(0.99,thetalast[3]/cons)

    muprop = c(muprop1,muprop2,muprop3)
    sigmaprop1 = muprop[1]*(1-muprop[1])/divproprho
    sigmaprop2 = muprop[2]*(1-muprop[2])/divproppsi
    sigmaprop3 = muprop[3]*(1-muprop[3])/divpropphi
    sigmaprop3 = abs(sigmaprop3)

    #    sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #   sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop = c(sigmaprop1,sigmaprop2,sigmaprop3)
    omega = runif(1,0.99,1)
    theta1cand = min(0.99,rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    theta2cand = min(0.99,rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    theta3cand = min(0.99,rbeta_rep_ab(1,omega*muprop[3],sigmaprop[3],-1,1))

    thetacand = c(theta1cand, theta2cand, theta3cand)


    lognum = posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    logden = posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior)


    lognum = lognum + log(dbeta_rep_cpp(thetalast[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetalast[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_ab(thetalast[3],omega*muprop[3],sigmaprop[3],-1,1))

    logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_ab(thetacand[3],omega*muprop[3],sigmaprop[3],-1,1))


    logunif= log(runif(1))

    if(logunif<lognum - logden)
    {
      count = count +1
      rhoF[i] <- thetacand[1]
      psiF[i] = thetacand[2]
      phiF[i] = thetacand[3]

    }else
    {
      rhoF[i] <- thetalast[1]
      psiF[i] = thetalast[2]
      phiF[i] = thetalast[3]
    }
    print(c(i,rhoF[i],psiF[i],phiF[i],count/i))
    #gc()
  }

  gc()

  rhoburn=rhoF[(burn+1):iter]
  rhoval= rhoburn[seq((burn+1),iter-burn,thin)]

  psiburn=psiF[(burn+1):iter]
  psival= psiburn[seq((burn+1),iter-burn,thin)]

  phiburn=phiF[(burn+1):iter]
  phival= phiburn[seq((burn+1),iter-burn,thin)]

  probacc = sum(count)/iter



  betaFval = matrix(0,length(rhoval),ncol(xobs))
  sigmaFval = 0
  nmc = length(rhoval)
  n = length(y)
  p = ncol(xobs)
  thetaFval = matrix(0,nmc,p+4)

  for (i in 1:length(rhoval)){
    #i = 1
    rhovalnum = rhoval[i]
    psivalnum = psival[[i]]
    phivalnum = phival[[i]]
    vdag =spatimecovar_2_cpp(lag = lags,adjmatinf = adjmatinf,rho = rhovalnum,
                             psi = psivalnum,phi = phivalnum)
    vdaginv = solve(vdag)
    vbetaest = solve(t(xobs)%*%vdaginv%*%xobs)
    betaest = vbetaest%*%t(xobs)%*%vdaginv%*%y
    S2 = t(y-xobs%*%betaest)%*%vdaginv%*%(y-xobs%*%betaest)
    betaFval[i,] = rmvnorm(1,betaest,vbetaest)
    sigmaFval[i] = rinvgamma(1,0.5*(n-p), S2/2)
    thetaval = c(betaFval[i,],sigmaFval[i],rhovalnum,psivalnum,phivalnum)
    thetaFval[i,] = thetaval
  }

  #tau2F = sigma2F*psiF
  thetaF = cbind(rhoF,psiF,phiF)


  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)

  log_lik = apply(thetaFval,1,loglikespatemp_fix,Xmat = xobs,y,lags = lags,adjmatinf = adjmatinf)
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  rhoest = quantile(rhoval,probs = c(0.025,0.5,0.975))
  psiest = quantile(psival,probs = c(0.025,0.5,0.975))
  phiest = quantile(phival,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigmaFval,probs = c(0.025,0.5,0.975))

  result = list(betaF = betaFval, sigmaF = sigmaFval, rhoF = rhoval, psiF = psival, thetaF = thetaF,
                thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
                psiest = psiest, phiest = phiest,probacc = probacc,log_lik = log_lik)

  return(result)
}







#' @export
bayesspatempauto_fix = function(y,xobs,thetaini,iter,burn,thin, lags,
                            adjmatinf, aprior = c(4,2,4), bprior = c(2,4,2),
                            divproprho = 9, divproppsi = 9, divpropphi = 9){

  # y= yobs; xobs = xobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;
  #adjmatinf = adjm1;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time_obs





  rhoF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  psiF = rep(0,iter)
  phiF = rep(0,iter)
  rhoF[1] = thetaini[1]
  psiF[1] = thetaini[2]
  phiF[1] = thetaini[3]

  for(i in 2:iter){
    # i = 2
    adjmatinf[upper.tri(adjmatinf)]=0

    thetalast = c(rhoF[i-1],psiF[i-1],phiF[i-1])

    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))

    #lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior))
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior))


    #lognum = lognum + log(dbeta(thetalast[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetalast[2],apriorprop[2],bpriorprop[2])) +
     # log(dbeta(thetalast[3],apriorprop[3],bpriorprop[3]))

    #logden = logden + log(dbeta(thetacand[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetacand[2],apriorprop[2],bpriorprop[2])) +
     # log(dbeta(thetacand[3],apriorprop[3],bpriorprop[3]))

    cons = 1
    muprop1 = min(0.99,thetalast[1]/cons)
    muprop2 = min(0.99,thetalast[2]/cons)
    muprop3 = min(0.99,thetalast[3]/cons)

    muprop = c(muprop1,muprop2,muprop3)
    sigmaprop1 = muprop[1]*(1-muprop[1])/divproprho
    sigmaprop2 = muprop[2]*(1-muprop[2])/divproppsi
    sigmaprop3 = muprop[3]*(1-muprop[3])/divpropphi

    #    sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #   sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop = c(sigmaprop1,sigmaprop2,sigmaprop3)
    omega = runif(1,0.99,1)
    theta1cand = min(0.99,rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    theta2cand = min(0.99,rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    theta3cand = min(0.99,rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))


    #if(theta1cand> 0.93){
      #theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
     # theta1cand = runif(1,0.1,0.6)
    #}

    #if(theta2cand > 0.93){
      #theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
     # theta2cand = runif(1,0.1,0.6)
    #}


    #if(theta3cand > 0.93){
      #theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
     # theta3cand = runif(1,0.1,0.6)
    #}



    if(theta1cand<0.05){
      #theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
      theta1cand = runif(1,0.1,0.2)
    }

    if(theta2cand<0.05){
      #theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
      theta2cand = runif(1,0.1,0.2)
    }


    if(theta3cand<0.05){
      #theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
      theta3cand = runif(1,0.1,0.2)
    }


    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    thetacand = c(theta1cand,theta2cand,theta3cand)
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))
    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))
    # a = posteriorvar3_cpp(thetacand,adjmatinf = adj_matrix,lags,X,ytot,a=c(1,1,1),b =c(1,1,1))
    # posterior_spatemp_cpp(thetacand,adjmatinf = adj_matrix,lag = N_time,X=xobs,y=y,a=c(1,1,1),b =c(1,1,1))
    #posterior_spatempver_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))

    lognum = posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    logden = posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior)

   # lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior)) #+ logpwcand
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior)) #+ logpwlast



    #if(is.na(num)){
     # num = runif(1)
    #}

    #if(is.na(den)){
     # den = runif(1)
    #}

    #lognum = log(num) #+ logpwcand
    #logden = log(den) #+ logpwlast

    lognum = lognum + log(dbeta_rep_cpp(thetalast[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetalast[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_cpp(thetalast[3],omega*muprop[3],sigmaprop[3])) #+ logpwlast

    logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_cpp(thetacand[3],omega*muprop[3],sigmaprop[3])) #+ logpwcand



    logunif= log(runif(1))

    if(logunif<lognum - logden)
    {
      count = count +1
      rhoF[i] <- thetacand[1]
      psiF[i] = thetacand[2]
      phiF[i] = thetacand[3]

    }else
    {
      rhoF[i] <- thetalast[1]
      psiF[i] = thetalast[2]
      phiF[i] = thetalast[3]
    }
    print(c(i,rhoF[i],psiF[i],phiF[i],count/i))
    #gc()
  }

  gc()

  rhoburn=rhoF[(burn+1):iter]
  rhoval= rhoburn[seq((burn+1),iter-burn,thin)]

  psiburn=psiF[(burn+1):iter]
  psival= psiburn[seq((burn+1),iter-burn,thin)]

  phiburn=phiF[(burn+1):iter]
  phival= phiburn[seq((burn+1),iter-burn,thin)]

  probacc = sum(count)/iter



  betaFval = matrix(0,length(rhoval),ncol(xobs))
  sigmaFval = 0
  nmc = length(rhoval)
  n = length(y)
  p = ncol(xobs)
  thetaFval = matrix(0,nmc,p+4)

  for (i in 1:length(rhoval)){
    #i = 1
    rhovalnum = rhoval[i]
    psivalnum = psival[[i]]
    phivalnum = phival[[i]]
    vdag =spatimecovar_2_cpp(lag = lags,adjmatinf = adjmatinf,rho = rhovalnum,
                             psi = psivalnum,phi = phivalnum)
    vdaginv = solve(vdag)
    vbetaest = solve(t(xobs)%*%vdaginv%*%xobs)
    betaest = vbetaest%*%t(xobs)%*%vdaginv%*%y
    S2 = t(y-xobs%*%betaest)%*%vdaginv%*%(y-xobs%*%betaest)
    sigmaFval[i] = rinvgamma(1,0.5*(n-p), S2/2)
    betaFval[i,] = rmvnorm(1,betaest,sigmaFval[i]*vbetaest)
    thetaval = c(betaFval[i,],sigmaFval[i],rhovalnum,psivalnum,phivalnum)
    thetaFval[i,] = thetaval
  }

  tau2val = sigmaFval*psival
  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)


  ##### descomentar despues (otro paper)
  #log_lik = apply(thetaFval,1,loglikespatemp_fix,Xmat = xobs,y,lags = lags,adjmatinf = adjmatinf)
###################################


  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  rhoest = quantile(rhoval,probs = c(0.025,0.5,0.975))
  psiest = quantile(psival,probs = c(0.025,0.5,0.975))
  tau2est = quantile(tau2val,probs = c(0.025,0.5,0.975))
  phiest = quantile(phival,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigmaFval,probs = c(0.025,0.5,0.975))

  #### descomentar despues (otro paper)

  #result = list(betaF = betaFval, sigmaF = sigmaFval, rhoF = rhoval, psiF = psival,tau2F = tau2val,
   #             thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
    #            psiest = psiest,phiest = phiest, tau2est = tau2est, probacc = probacc,log_lik = log_lik)

  result = list(betaF = betaFval, sigmaF = sigmaFval, rhoF = rhoval, psiF = psival,tau2F = tau2val,
                thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
                psiest = psiest,phiest = phiest, tau2est = tau2est, probacc = probacc)

  return(result)
}




#' @export
bayesspatempautosar_fix = function(y,xobs,thetaini,iter,burn,thin, lags,
                                adjmat, aprior = c(4,2,4), bprior = c(2,4,2),
                                divproprho = 9, divproppsi = 9, divpropphi = 9){

  # y= yobs; xobs = xobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;
  #adjmatinf = adjm1;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time_obs





  rhoF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  psiF = rep(0,iter)
  phiF = rep(0,iter)
  rhoF[1] = thetaini[1]
  psiF[1] = thetaini[2]
  phiF[1] = thetaini[3]

  for(i in 2:iter){
    # i = 2
   # adjmatinf[upper.tri(adjmatinf)]=0

    thetalast = c(rhoF[i-1],psiF[i-1],phiF[i-1])

    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))

    #lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior))
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior))


    #lognum = lognum + log(dbeta(thetalast[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetalast[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetalast[3],apriorprop[3],bpriorprop[3]))

    #logden = logden + log(dbeta(thetacand[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetacand[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetacand[3],apriorprop[3],bpriorprop[3]))

    cons = 1
    muprop1 = min(0.99,thetalast[1]/cons)
    muprop2 = min(0.99,thetalast[2]/cons)
    muprop3 = min(0.99,thetalast[3]/cons)

    muprop = c(muprop1,muprop2,muprop3)
    sigmaprop1 = muprop[1]*(1-muprop[1])/divproprho
    sigmaprop2 = muprop[2]*(1-muprop[2])/divproppsi
    sigmaprop3 = muprop[3]*(1-muprop[3])/divpropphi

    #    sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #   sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop = c(sigmaprop1,sigmaprop2,sigmaprop3)
    omega = runif(1,0.99,1)
    theta1cand = min(0.99,rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    theta2cand = min(0.99,rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    theta3cand = min(0.99,rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))


    #if(theta1cand> 0.93){
    #theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    # theta1cand = runif(1,0.1,0.6)
    #}

    #if(theta2cand > 0.93){
    #theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    # theta2cand = runif(1,0.1,0.6)
    #}


    #if(theta3cand > 0.93){
    #theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
    # theta3cand = runif(1,0.1,0.6)
    #}




    if(theta1cand<0.05){
      #theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
      theta1cand = runif(1,0.1,0.2)
    }

    if(theta2cand<0.05){
      #theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
      theta2cand = runif(1,0.1,0.2)
    }


    if(theta3cand<0.05){
      #theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
      theta3cand = runif(1,0.1,0.2)
    }


    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    thetacand = c(theta1cand,theta2cand,theta3cand)
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))
    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))
    # a = posteriorvar3_cpp(thetacand,adjmatinf = adj_matrix,lags,X,ytot,a=c(1,1,1),b =c(1,1,1))
    # posterior_spatemp_cpp(thetacand,adjmatinf = adj_matrix,lag = N_time,X=xobs,y=y,a=c(1,1,1),b =c(1,1,1))
    #posterior_spatempver_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))

    lognum = posterior_spatempsar_cpp(theta = thetacand,adjmat = adjmat, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    logden = posterior_spatempsar_cpp(theta = thetalast,adjmat = adjmat,lag = lags, X = xobs,y = y, a = aprior, b= bprior)

    # lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior)) #+ logpwcand
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior)) #+ logpwlast



    #if(is.na(num)){
    # num = runif(1)
    #}

    #if(is.na(den)){
    # den = runif(1)
    #}

    #lognum = log(num) #+ logpwcand
    #logden = log(den) #+ logpwlast

    lognum = lognum + log(dbeta_rep_cpp(thetalast[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetalast[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_cpp(thetalast[3],omega*muprop[3],sigmaprop[3])) #+ logpwlast

    logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_cpp(thetacand[3],omega*muprop[3],sigmaprop[3])) #+ logpwcand



    logunif= log(runif(1))

    if(logunif<lognum - logden)
    {
      count = count +1
      rhoF[i] <- thetacand[1]
      psiF[i] = thetacand[2]
      phiF[i] = thetacand[3]

    }else
    {
      rhoF[i] <- thetalast[1]
      psiF[i] = thetalast[2]
      phiF[i] = thetalast[3]
    }
    print(c(i,rhoF[i],psiF[i],phiF[i],count/i))
    #gc()
  }

  gc()

  rhoburn=rhoF[(burn+1):iter]
  rhoval= rhoburn[seq((burn+1),iter-burn,thin)]

  psiburn=psiF[(burn+1):iter]
  psival= psiburn[seq((burn+1),iter-burn,thin)]

  phiburn=phiF[(burn+1):iter]
  phival= phiburn[seq((burn+1),iter-burn,thin)]

  probacc = sum(count)/iter



  betaFval = matrix(0,length(rhoval),ncol(xobs))
  sigmaFval = 0
  nmc = length(rhoval)
  n = length(y)
  p = ncol(xobs)
  thetaFval = matrix(0,nmc,p+4)

  for (i in 1:length(rhoval)){
    #i = 1
    rhovalnum = rhoval[i]
    psivalnum = psival[[i]]
    phivalnum = phival[[i]]
    vdag =spatimecovarsar_2_cpp(lag = lags,adjmat = adjmat,rho = rhovalnum,
                             psi = psivalnum,phi = phivalnum)
    vdaginv = solve(vdag)
    vbetaest = solve(t(xobs)%*%vdaginv%*%xobs)
    betaest = vbetaest%*%t(xobs)%*%vdaginv%*%y
    S2 = t(y-xobs%*%betaest)%*%vdaginv%*%(y-xobs%*%betaest)
    sigmaFval[i] = rinvgamma(1,0.5*(n-p), S2/2)
    betaFval[i,] = rmvnorm(1,betaest,sigmaFval[i]*vbetaest)
    thetaval = c(betaFval[i,],sigmaFval[i],rhovalnum,psivalnum,phivalnum)
    thetaFval[i,] = thetaval
  }

  tau2val = sigmaFval*psival
  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)


  ##### descomentar despues (otro paper)
  #log_lik = apply(thetaFval,1,loglikespatemp_fix,Xmat = xobs,y,lags = lags,adjmatinf = adjmatinf)
  ###################################


  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  rhoest = quantile(rhoval,probs = c(0.025,0.5,0.975))
  psiest = quantile(psival,probs = c(0.025,0.5,0.975))
  tau2est = quantile(tau2val,probs = c(0.025,0.5,0.975))
  phiest = quantile(phival,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigmaFval,probs = c(0.025,0.5,0.975))

  #### descomentar despues (otro paper)

  #result = list(betaF = betaFval, sigmaF = sigmaFval, rhoF = rhoval, psiF = psival,tau2F = tau2val,
  #             thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
  #            psiest = psiest,phiest = phiest, tau2est = tau2est, probacc = probacc,log_lik = log_lik)

  result = list(betaF = betaFval, sigmaF = sigmaFval, rhoF = rhoval, psiF = psival,tau2F = tau2val,
                thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
                psiest = psiest,phiest = phiest, tau2est = tau2est, probacc = probacc)

  return(result)
}

















#' @export
dbeta_rep_ab <- function(x, mu, sigma2, a, b) {
  if (a >= b) stop("Require a < b.")
  if (mu <= a || mu >= b) stop("mu must be in (a, b).")
  if (sigma2 <= 0) stop("sigma2 must be > 0.")

  width <- b - a

  # rescale to (0,1)
  z    <- (x - a) / width
  mu_s <- (mu - a) / width
  s2_s <- sigma2 / (width^2)

  if (mu_s <= 0 || mu_s >= 1) stop("Scaled mu not in (0,1).")
  if (s2_s >= mu_s * (1 - mu_s)) {
    stop("Invalid sigma2: must be less than (mu-a)*(b-mu).")
  }

  # recover alpha, beta
  cons  <- ((1 - mu_s) / s2_s) - (1 / mu_s)
  alpha <- cons * (mu_s^2)
  beta  <- cons * mu_s * (1 - mu_s)

  if (alpha <= 0 || beta <= 0) stop("alpha, beta must be > 0.")
  # density
  ifelse(x <= a | x >= b, 0, dbeta(z, alpha, beta) / width)
}


#' @export
rbeta_rep_ab <- function(n, mu, sigma2, a, b) {
  if (a >= b) stop("Require a < b.")
  if (mu <= a || mu >= b) stop("mu must be in (a, b).")
  if (sigma2 <= 0) stop("sigma2 must be > 0.")

  width <- b - a
  mu_s  <- (mu - a) / width
  s2_s  <- sigma2 / (width^2)

  if (mu_s <= 0 || mu_s >= 1) stop("Scaled mu not in (0,1).")
  if (s2_s >= mu_s * (1 - mu_s)) {
    stop("sigma2 must be < (mu-a)*(b-mu).")
  }

  # Recover alpha, beta from (mu*, s2*)
  cons  <- ((1 - mu_s) / s2_s) - (1 / mu_s)
  alpha <- cons * (mu_s^2)
  beta  <- cons * mu_s * (1 - mu_s)
  if (alpha <= 0 || beta <= 0) stop("alpha,beta must be > 0 (check mu,sigma2).")

  a + width * rbeta(n, alpha, beta)
}




#' @export
bayesspatempcensauto_fix = function(y,xtot,thetaini,indicens, iter,burn,thin, lags,
                                adjmatinf, aprior = c(4,2,4), bprior = c(2,4,2),
                                divproprho = 9, divproppsi = 9, divpropphi = 9,lower,upper){

  # y= ytotobs; xtot = xtotobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;indicens = indicens3;
  #adjmatinf = adj_matrix;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time_obs




  yobs = y[indicens == 0]
  xobs = xtot[indicens ==0, ]
  rhoF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  psiF = rep(0,iter)
  phiF = rep(0,iter)
  p = ncol(xtot)
  betaF = matrix(0,iter,p)
  sigma2F = 0
  n = length(y)
  ycensF = matrix(0,iter,sum(indicens==1))
  rhoF[1] = thetaini[1]
  psiF[1] = thetaini[2]
  phiF[1] = thetaini[3]
  #corrmatini = spatimecovar_2_cpp(lag = lags,adjmatinf =adjmatinf, rho = rhoF[1], psi = psiF[1],
   #                              phi = phiF[1])
  #invcorrmatini = solve(corrmatini)
  #betaini = solve(t(xtot)%*%invcorrmatini%*%xtot)%*%t(xtot)%*%invcorrmatini%*%y

  #sigma2ini = (t(y-xtot%*%betaini)%*%corrmatini%*%(y-xtot%*%betaini))/(n-p)
  #betaF[1,] = betaini
  #sigma2F[1] = as.numeric(sigma2ini)

  #covmatini = sigma2F[1]*corrmatini

  ############################################################
  #covobsobsini = covmatini[indicens ==0, indicens ==0]
  #covcenscensini = covmatini[indicens ==1, indicens ==1]
  #covcensobsini = covmatini[indicens ==1, indicens ==0]
  #covobscensini = covmatini[indicens ==0, indicens ==1]
  ############################################################

  #yobs = y[indicens ==0]
  #xobserv = xtot[indicens ==0,]
  #xcens = xtot[indicens ==1,]
  #covobsobsiniinv = solve(covobsobsini)
  #mucensini = xcens%*%betaini + covcensobsini%*%covobsobsiniinv%*%(yobs - xobserv%*%betaini)
  #Sigmacensini = covcenscensini - covcensobsini%*%covobsobsiniinv%*%covobscensini
  ycensF[1,] = y[indicens ==1]
  #rtmvnorm(1,mean =as.vector(mucensini), sigma = Sigmacensini,lower = lower,upper = upper,algorithm = "gibbs",thinning = 2)
  ycom = y
  ycom[indicens ==1] = ycensF[1,]
  ycom2 = matrix(0,iter,)
  adjmatinf[upper.tri(adjmatinf)]=0


  for(i in 2:iter){
   #i = 2
    thetalast = c(rhoF[i-1],psiF[i-1],phiF[i-1])

    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))

    #lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior))
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior))


    #lognum = lognum + log(dbeta(thetalast[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetalast[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetalast[3],apriorprop[3],bpriorprop[3]))

    #logden = logden + log(dbeta(thetacand[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetacand[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetacand[3],apriorprop[3],bpriorprop[3]))
    cons = 1
    muprop1 = min(0.99,thetalast[1]/cons)
    muprop2 = min(0.99,thetalast[2]/cons)
    muprop3 = min(0.99,thetalast[3]/cons)

    muprop = c(muprop1,muprop2,muprop3)
    sigmaprop1 = muprop[1]*(1-muprop[1])/divproprho
    sigmaprop2 = muprop[2]*(1-muprop[2])/divproppsi
    sigmaprop3 = muprop[3]*(1-muprop[3])/divpropphi
    sigmaprop3 = abs(sigmaprop3)

    #    sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #   sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop = c(sigmaprop1,sigmaprop2,sigmaprop3)
    omega = runif(1,0.99,1)
    theta1cand = min(0.99,rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    theta2cand = min(0.99,rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    #theta3cand = min(0.99,rbeta_rep_ab(1,omega*muprop[3],sigmaprop[3],-1,1))
    theta3cand = min(0.99,rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
    #theta3cand = rnorm(1,muprop[3],sigmaprop[3])

    #if(theta1cand > 0.95){
     # theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
     #theta1cand = runif(1,0.1,0.6)
    #}

    #if(theta2cand > 0.95){
     # theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    #  theta2cand = runif(1,0.1,0.45)
    #}


    #if(theta3cand> 0.95){
     # theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
      #theta3cand = runif(1,0.1,0.6)
    #}

   # if(theta3cand<0.05){
      #theta3cand =max(runif(1,0.05,1),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
    #  theta3cand = runif(1,0.05,0.3)
    #}



#    if(theta1cand<0.05){
 #     theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
  #    theta1cand = runif(1,0.1,0.2)
  #  }

  # if(theta2cand<0.03){
    #theta2cand =max(runif(1,0.03,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
   #   theta2cand = runif(1,0.03,0.15)
  #  }


    #if(theta3cand<0.05){
     #theta3cand =max(runif(1,0.05,0.15),rbeta_rep_ab(1,omega*muprop[3],sigmaprop[3],-1,1))
     #theta3cand = runif(1,0.1,0.2)
    #}


    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    thetacand = c(theta1cand,theta2cand,theta3cand)

    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))
    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))
    # a = posteriorvar3_cpp(thetacand,adjmatinf = adj_matrix,lags,X,ytot,a=c(1,1,1),b =c(1,1,1))
    # posterior_spatemp_cpp(thetacand,adjmatinf = adj_matrix,lag = N_time,X=xobs,y=y,a=c(1,1,1),b =c(1,1,1))
    #posterior_spatempver_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))

    auxnum = posterior_spatempcens_2_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xtot,y = ycom, a = aprior, b = bprior)
    auxden = posterior_spatempcens_2_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xtot,y = ycom, a = aprior, b= bprior)

    #num = auxnum$posterior
    #den = auxden$posterior

    # lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior)) #+ logpwcand
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior)) #+ logpwlast



  #  if(is.na(num)){
   #   num = runif(1)
    #}

    #if(is.na(den)){
     # den = runif(1)
    #}

    lognum = auxnum$posterior_log #+ logpwcand
    logden = auxden$posterior_log #+ logpwlast

    lognum = lognum + log(dbeta_rep_cpp(thetalast[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetalast[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_cpp(thetalast[3],omega*muprop[3],sigmaprop[3])) #+ logpwlast
   # (x, mu, sigma2, a, b)

    logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_cpp(thetacand[3],omega*muprop[3],sigmaprop[3])) #+ logpwcand
   # logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
    # log(dnorm(thetacand[3],omega*muprop[3],sigmaprop[3])) #+ logpwcand



    logunif= log(runif(1))
    n = length(y)
    p = ncol(xtot)

    if(logunif<lognum - logden)
    {
      count = count +1
      rhoF[i] <- thetacand[1]
      psiF[i] = thetacand[2]
      phiF[i] = thetacand[3]
      corrmat = auxnum$varcovspatemp
      betaest = auxnum$betaest
      varbeta = auxnum$varbeta
      s2est = auxnum$S2

    }else
    {
      rhoF[i] <- thetalast[1]
      psiF[i] = thetalast[2]
      phiF[i] = thetalast[3]
      corrmat = auxden$varcovspatemp
      betaest = auxden$betaest
      varbeta = auxden$varbeta
      s2est = auxden$S2
    }

    sigma2F[i] = rinvgamma(1,0.5*(n-p), (n-p)*s2est/2)
    betaF[i,] = rmvnorm(1,betaest,sigma2F[i]*varbeta)
    covmat = sigma2F[i] *corrmat
    ############################################################
    covobsobs = covmat[indicens ==0, indicens ==0]
    covcenscens = covmat[indicens ==1, indicens ==1]
    covcensobs = covmat[indicens ==1, indicens ==0]
    covobscens = covmat[indicens ==0, indicens ==1]
    ############################################################
    yobs = y[indicens ==0]
    xobserv = xtot[indicens ==0,]
    xcens = xtot[indicens ==1,]
    covobsobsinv = solve(covobsobs)
    mucens = xcens%*%betaF[i,] + covcensobs%*%covobsobsinv%*%(yobs - xobserv%*%betaF[i,])
    Sigmacens = covcenscens - covcensobs%*%covobsobsinv%*%covobscens
    ycensF[i,] = rtmvnorm(1,mean =as.vector(mucens), sigma = Sigmacens,lower = lower,upper = upper,algorithm = "gibbs",thinning = 4)
    ycom[indicens ==1] = ycensF[i,]
    #print(c(i,sigma2F[i],rhoF[i],psiF[i],phiF[i],count/i))
    cat(sprintf("\rIteration: %d", i))
    #gc()
  }
  cat("\n")
  gc()
  tau2F = sigma2F*psiF
  thetaF = cbind(betaF,sigma2F,rhoF,psiF,phiF,tau2F,ycensF)

  rhoburn=rhoF[(burn+1):iter]
 # rhoval= rhoburn[seq((burn+1),iter-burn,thin)]
  rhoval= rhoburn[seq(1,iter-burn,thin)]

  psiburn=psiF[(burn+1):iter]
  psival= psiburn[seq(1,iter-burn,thin)]

  phiburn=phiF[(burn+1):iter]
  phival= phiburn[seq(1,iter-burn,thin)]

  sigma2Fburn=sigma2F[(burn+1):iter]
  sigma2Fval= sigma2Fburn[seq(1,iter-burn,thin)]

  tau2val = sigma2Fval*psival

  betaFburn=betaF[(burn+1):iter,]
  betaFval= betaFburn[seq(1,iter-burn,thin),]

  ycensFburn = ycensF[(burn+1):iter,]
  ycensFval = ycensFburn[seq(1,iter-burn,thin),]

  probacc = sum(count)/iter

  nmc = length(rhoval)
  thetaFval = cbind(betaFval,sigma2Fval,rhoval,psival,phival,ycensFval)


  rhoest = quantile(rhoval,probs = c(0.025,0.5,0.975))
  psiest = quantile(psival,probs = c(0.025,0.5,0.975))
  tau2est = quantile(tau2val,probs = c(0.025,0.5,0.975))
  phiest = quantile(phival,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigma2Fval,probs = c(0.025,0.5,0.975))


  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)


  #### maybe later
  #log_lik = apply(thetaFval,1,loglikespatempcens_fix,Xmat = xtot,ytotobs = y,lags = lags,adjmatinf = adjmatinf,indicens = indicens)
  ######


  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  #result = list(betaF = betaFval, sigma2F = sigma2Fval, rhoF = rhoval, psiF = psival,phiF = phival,tau2F = tau2val,
   #            ycensF = ycensFval, thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
    #          psiest = psiest,phiest = phiest, tau2est = tau2est, probacc = probacc,log_lik = log_lik)

   result = list(betaF = betaFval, sigma2F = sigma2Fval, rhoF = rhoval, psiF = psival,phiF = phival,tau2F = tau2val,
                ycensF = ycensFval, thetaF =thetaF, thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
               psiest = psiest,phiest = phiest, tau2est = tau2est, probacc = probacc)


    return(result)
}





#' @export
bayesspatempcensautosar_fix = function(y,xtot,thetaini,indicens, iter,burn,thin, lags,
                                    adjmat, aprior = c(4,2,4), bprior = c(2,4,2),
                                    divproprho = 9, divproppsi = 9, divpropphi = 9,lower,upper){

  # y= ytotobs; xtot = xtotobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;indicens = indicens3;
  #adjmatinf = adj_matrix;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time_obs




  yobs = y[indicens == 0]
  xobs = xtot[indicens ==0, ]
  rhoF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  psiF = rep(0,iter)
  phiF = rep(0,iter)
  p = ncol(xtot)
  betaF = matrix(0,iter,p)
  sigma2F = 0
  n = length(y)
  ycensF = matrix(0,iter,sum(indicens==1))
  rhoF[1] = thetaini[1]
  psiF[1] = thetaini[2]
  phiF[1] = thetaini[3]
  #corrmatini = spatimecovar_2_cpp(lag = lags,adjmatinf =adjmatinf, rho = rhoF[1], psi = psiF[1],
  #                              phi = phiF[1])
  #invcorrmatini = solve(corrmatini)
  #betaini = solve(t(xtot)%*%invcorrmatini%*%xtot)%*%t(xtot)%*%invcorrmatini%*%y

  #sigma2ini = (t(y-xtot%*%betaini)%*%corrmatini%*%(y-xtot%*%betaini))/(n-p)
  #betaF[1,] = betaini
  #sigma2F[1] = as.numeric(sigma2ini)

  #covmatini = sigma2F[1]*corrmatini

  ############################################################
  #covobsobsini = covmatini[indicens ==0, indicens ==0]
  #covcenscensini = covmatini[indicens ==1, indicens ==1]
  #covcensobsini = covmatini[indicens ==1, indicens ==0]
  #covobscensini = covmatini[indicens ==0, indicens ==1]
  ############################################################

  #yobs = y[indicens ==0]
  #xobserv = xtot[indicens ==0,]
  #xcens = xtot[indicens ==1,]
  #covobsobsiniinv = solve(covobsobsini)
  #mucensini = xcens%*%betaini + covcensobsini%*%covobsobsiniinv%*%(yobs - xobserv%*%betaini)
  #Sigmacensini = covcenscensini - covcensobsini%*%covobsobsiniinv%*%covobscensini
  ycensF[1,] = y[indicens ==1]
  #rtmvnorm(1,mean =as.vector(mucensini), sigma = Sigmacensini,lower = lower,upper = upper,algorithm = "gibbs",thinning = 2)
  ycom = y
  ycom[indicens ==1] = ycensF[1,]
  ycom2 = matrix(0,iter,)
  #adjmatinf[upper.tri(adjmatinf)]=0


  for(i in 2:iter){
    #i = 2
    thetalast = c(rhoF[i-1],psiF[i-1],phiF[i-1])

    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))

    #lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior))
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior))


    #lognum = lognum + log(dbeta(thetalast[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetalast[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetalast[3],apriorprop[3],bpriorprop[3]))

    #logden = logden + log(dbeta(thetacand[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetacand[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetacand[3],apriorprop[3],bpriorprop[3]))
    cons = 1
    muprop1 = min(0.965,thetalast[1]/cons)
    muprop2 = min(0.965,thetalast[2]/cons)
    muprop3 = min(0.965,thetalast[3]/cons)

    muprop = c(muprop1,muprop2,muprop3)
    sigmaprop1 = abs(muprop[1]*(1-muprop[1])/divproprho)
    sigmaprop1 = min(abs((muprop[1]-(-0.97))*(0.97-muprop[1])-0.005),sigmaprop1)
    sigmaprop1 = abs(sigmaprop1)
    sigmaprop2 = abs(muprop[2]*(1-muprop[2])/divproppsi)
    sigmaprop3 = abs(muprop[3]*(1-
                                  muprop[3])/divpropphi)
    sigmaprop3 = min(abs((muprop[3]-(-0.97))*(0.97-muprop[3])-0.005),sigmaprop3)
    sigmaprop3 = abs(sigmaprop3)

    #    sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #   sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop = c(sigmaprop1,sigmaprop2,sigmaprop3)
    omega = runif(1,0.99,1)
    theta1cand = min(0.99,rbeta_rep_ab(1,omega*muprop[1],sigmaprop[1],-0.97,0.97))
    theta2cand = min(0.99,rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    theta3cand = min(0.99,rbeta_rep_ab(1,omega*muprop[3],sigmaprop[3],-0.97,0.97))
    #theta3cand = rnorm(1,muprop[3],sigmaprop[3])

    #if(theta1cand > 0.95){
    # theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    #theta1cand = runif(1,0.1,0.6)
    #}

    #if(theta2cand > 0.95){
    # theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    #  theta2cand = runif(1,0.1,0.45)
    #}


    #if(theta3cand> 0.95){
    # theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
    #theta3cand = runif(1,0.1,0.6)
    #}


    if(theta1cand<0.05){
      theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
      theta1cand = runif(1,0.1,0.2)
    }

    if(theta2cand<0.05){
      theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
      theta2cand = runif(1,0.1,0.2)
    }


    if(theta3cand<0.05){
      theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
      theta3cand = runif(1,0.1,0.2)
    }


    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    thetacand = c(theta1cand,theta2cand,theta3cand)

    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))
    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))
    # a = posteriorvar3_cpp(thetacand,adjmatinf = adj_matrix,lags,X,ytot,a=c(1,1,1),b =c(1,1,1))
    # posterior_spatemp_cpp(thetacand,adjmatinf = adj_matrix,lag = N_time,X=xobs,y=y,a=c(1,1,1),b =c(1,1,1))
    #posterior_spatempver_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))

    auxnum = posterior_spatempcenssar_2_cpp(theta = thetacand,adjmat = adjmat, lag = lags, X = xtot,y = ycom, a = aprior, b = bprior)
    auxden = posterior_spatempcenssar_2_cpp(theta = thetalast,adjmat = adjmat,lag = lags, X = xtot,y = ycom, a = aprior, b= bprior)

    #num = auxnum$posterior
    #den = auxden$posterior

    # lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior)) #+ logpwcand
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior)) #+ logpwlast



    #  if(is.na(num)){
    #   num = runif(1)
    #}

    #if(is.na(den)){
    # den = runif(1)
    #}

    lognum = auxnum$posterior_log #+ logpwcand
    logden = auxden$posterior_log #+ logpwlast

    lognum = lognum + log(dbeta_rep_ab(thetalast[1],omega*muprop[1],sigmaprop[1],-1,1)) + log(dbeta_rep_cpp(thetalast[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_ab(thetalast[3],omega*muprop[3],sigmaprop[3],-1,1)) #+ logpwlast
    # (x, mu, sigma2, a, b)

    logden = logden + log(dbeta_rep_ab(thetacand[1],omega*muprop[1],sigmaprop[1],-1,1)) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_ab(thetacand[3],omega*muprop[3],sigmaprop[3],-1,1)) #+ logpwcand
    # logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
    # log(dnorm(thetacand[3],omega*muprop[3],sigmaprop[3])) #+ logpwcand



    logunif= log(runif(1))
    n = length(y)
    p = ncol(xtot)

    if(logunif<lognum - logden)
    {
      count = count +1
      rhoF[i] <- thetacand[1]
      psiF[i] = thetacand[2]
      phiF[i] = thetacand[3]
      corrmat = auxnum$varcovspatemp
      betaest = auxnum$betaest
      varbeta = auxnum$varbeta
      s2est = auxnum$S2

    }else
    {
      rhoF[i] <- thetalast[1]
      psiF[i] = thetalast[2]
      phiF[i] = thetalast[3]
      corrmat = auxden$varcovspatemp
      betaest = auxden$betaest
      varbeta = auxden$varbeta
      s2est = auxden$S2
    }

    sigma2F[i] = rinvgamma(1,0.5*(n-p), (n-p)*s2est/2)
    betaF[i,] = rmvnorm(1,betaest,sigma2F[i]*varbeta)
    covmat = sigma2F[i] *corrmat
    ############################################################
    covobsobs = covmat[indicens ==0, indicens ==0]
    covcenscens = covmat[indicens ==1, indicens ==1]
    covcensobs = covmat[indicens ==1, indicens ==0]
    covobscens = covmat[indicens ==0, indicens ==1]
    ############################################################
    yobs = y[indicens ==0]
    xobserv = xtot[indicens ==0,]
    xcens = xtot[indicens ==1,]
    covobsobsinv = solve(covobsobs)
    mucens = xcens%*%betaF[i,] + covcensobs%*%covobsobsinv%*%(yobs - xobserv%*%betaF[i,])
    Sigmacens = covcenscens - covcensobs%*%covobsobsinv%*%covobscens
    ycensF[i,] = rtmvnorm(1,mean =as.vector(mucens), sigma = Sigmacens,lower = lower,upper = upper,algorithm = "gibbs",thinning = 4)
    ycom[indicens ==1] = ycensF[i,]
    #print(c(i,sigma2F[i],rhoF[i],psiF[i],phiF[i],count/i))
    cat(sprintf("\rIteration: %d", i))
    #gc()
  }
  cat("\n")
  gc()

  tau2F = sigma2F*psiF
  thetaF = cbind(betaF,sigma2F,rhoF,psiF,phiF,tau2F,ycensF)

  rhoburn=rhoF[(burn+1):iter]
  rhoval= rhoburn[seq(1,iter-burn,thin)]

  psiburn=psiF[(burn+1):iter]
  psival= psiburn[seq(1,iter-burn,thin)]

  phiburn=phiF[(burn+1):iter]
  phival= phiburn[seq(1,iter-burn,thin)]

  sigma2Fburn=sigma2F[(burn+1):iter]
  sigma2Fval= sigma2Fburn[seq(1,iter-burn,thin)]

  tau2val = sigma2Fval*psival

  betaFburn=betaF[(burn+1):iter,]
  betaFval= betaFburn[seq(1,iter-burn,thin),]

  ycensFburn = ycensF[(burn+1):iter,]
  ycensFval = ycensFburn[seq(1,iter-burn,thin),]

  probacc = sum(count)/iter

  nmc = length(rhoval)
  thetaFval = cbind(betaFval,sigma2Fval,rhoval,psival,phival,ycensFval)


  rhoest = quantile(rhoval,probs = c(0.025,0.5,0.975))
  psiest = quantile(psival,probs = c(0.025,0.5,0.975))
  tau2est = quantile(tau2val,probs = c(0.025,0.5,0.975))
  phiest = quantile(phival,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigma2Fval,probs = c(0.025,0.5,0.975))


  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)


  #### maybe later
  #log_lik = apply(thetaFval,1,loglikespatempcens_fix,Xmat = xtot,ytotobs = y,lags = lags,adjmatinf = adjmatinf,indicens = indicens)
  ######


  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  #result = list(betaF = betaFval, sigma2F = sigma2Fval, rhoF = rhoval, psiF = psival,phiF = phival,tau2F = tau2val,
  #            ycensF = ycensFval, thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
  #          psiest = psiest,phiest = phiest, tau2est = tau2est, probacc = probacc,log_lik = log_lik)

  result = list(betaF = betaFval, sigma2F = sigma2Fval, rhoF = rhoval, psiF = psival,phiF = phival,tau2F = tau2val,
                ycensF = ycensFval,thetaF = thetaF,  thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
                psiest = psiest,phiest = phiest, tau2est = tau2est, probacc = probacc)


  return(result)
}





#' @export
bayesspatempcensauto_ar2_fix = function(y,xtot,thetaini,indicens, iter,burn,thin, lags,
                                    adjmatinf, aprior = c(4,2,4,4), bprior = c(2,4,2,2),
                                    divproprho = 9, divproppsi = 9, divpropphi = 9,lower,upper){

  # y= ytotobs; xtot = xtotobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;indicens = indicens3;
  #adjmatinf = adj_matrix;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time_obs




  yobs = y[indicens == 0]
  xobs = xtot[indicens ==0, ]
  rhoF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  psiF = rep(0,iter)
  phi1F = rep(0,iter)
  phi2F = rep(0,iter)
  p = ncol(xtot)
  betaF = matrix(0,iter,p)
  sigma2F = 0
  n = length(y)
  ycensF = matrix(0,iter,sum(indicens==1))
  rhoF[1] = thetaini[1]
  psiF[1] = thetaini[2]
  phi1F[1] = thetaini[3]
  phi2F[1] = thetaini[4]
  #corrmatini = spatimecovar_2_cpp(lag = lags,adjmatinf =adjmatinf, rho = rhoF[1], psi = psiF[1],
  #                              phi = phiF[1])
  #invcorrmatini = solve(corrmatini)
  #betaini = solve(t(xtot)%*%invcorrmatini%*%xtot)%*%t(xtot)%*%invcorrmatini%*%y

  #sigma2ini = (t(y-xtot%*%betaini)%*%corrmatini%*%(y-xtot%*%betaini))/(n-p)
  #betaF[1,] = betaini
  #sigma2F[1] = as.numeric(sigma2ini)

  #covmatini = sigma2F[1]*corrmatini

  ############################################################
  #covobsobsini = covmatini[indicens ==0, indicens ==0]
  #covcenscensini = covmatini[indicens ==1, indicens ==1]
  #covcensobsini = covmatini[indicens ==1, indicens ==0]
  #covobscensini = covmatini[indicens ==0, indicens ==1]
  ############################################################

  #yobs = y[indicens ==0]
  #xobserv = xtot[indicens ==0,]
  #xcens = xtot[indicens ==1,]
  #covobsobsiniinv = solve(covobsobsini)
  #mucensini = xcens%*%betaini + covcensobsini%*%covobsobsiniinv%*%(yobs - xobserv%*%betaini)
  #Sigmacensini = covcenscensini - covcensobsini%*%covobsobsiniinv%*%covobscensini
  ycensF[1,] = y[indicens ==1]
  #rtmvnorm(1,mean =as.vector(mucensini), sigma = Sigmacensini,lower = lower,upper = upper,algorithm = "gibbs",thinning = 2)
  ycom = y
  ycom[indicens ==1] = ycensF[1,]
  ycom2 = matrix(0,iter,)
  adjmatinf[upper.tri(adjmatinf)]=0


  for(i in 2:iter){
    #i = 2
    thetalast = c(rhoF[i-1],psiF[i-1],phi1F[i-1], phi2F[i-1])

    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))

    #lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior))
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior))


    #lognum = lognum + log(dbeta(thetalast[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetalast[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetalast[3],apriorprop[3],bpriorprop[3]))

    #logden = logden + log(dbeta(thetacand[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetacand[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetacand[3],apriorprop[3],bpriorprop[3]))
    cons = 1
    muprop1 = min(0.99,thetalast[1]/cons)
    muprop2 = min(0.99,thetalast[2]/cons)
    muprop3 = min(0.99,thetalast[3]/cons)
    muprop4 = min(0.99,thetalast[4]/cons)

    muprop = c(muprop1,muprop2,muprop3,muprop4)

    sigmaprop1 = muprop[1]*(1-muprop[1])/divproprho
    sigmaprop2 = muprop[2]*(1-muprop[2])/divproppsi
    #sigmaprop3 = muprop[3]*(1-muprop[3])/divpropphi
    #sigmaprop3 = abs(sigmaprop3)
    sigmaprop4 = muprop[4]*(1-muprop[4])/divpropphi
    sigmaprop4 = abs(sigmaprop4)

    #    sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #   sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop3 = 0
    sigmaprop = c(sigmaprop1,sigmaprop2,sigmaprop3,sigmaprop4)
    omega = runif(1,0.99,1)
    theta1cand = min(0.99,rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    theta2cand = min(0.99,rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    theta4cand = min(0.99,rbeta_rep_ab(1,omega*muprop[4],sigmaprop[4],-1,1))
    a= theta4cand -1
    b = 1-theta4cand
    muprop3 = if_else(muprop[3]<=a | muprop[3]>=b, (a+b)/2, muprop[3])
    sigmaprop3 = min(muprop3*(1-muprop3)/divpropphi, abs((muprop3-a)*(b-muprop3)-0.01))
    sigmaprop3 = if_else(sigmaprop3<=0,0.00002,sigmaprop3)
    #sigmaprop3 = max(0.002,abs(sigmaprop3))

    theta3cand = rbeta_rep_ab(1,muprop3,sigmaprop3,a,b)

    #theta3cand = rbeta_rep_ab(1,omega*muprop[3],sigmaprop[3],theta4cand-1,1-theta4cand)

    #theta3cand = rnorm(1,muprop[3],sigmaprop[3])

    #if(theta1cand > 0.95){
    # theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    #theta1cand = runif(1,0.1,0.6)
    #}

    #if(theta2cand > 0.95){
    # theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    #  theta2cand = runif(1,0.1,0.45)
    #}


    #if(theta3cand> 0.95){
    # theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
    #theta3cand = runif(1,0.1,0.6)
    #}



    #    if(theta1cand<0.05){
    #     theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    #    theta1cand = runif(1,0.1,0.2)
    #  }

    #    if(theta2cand<0.1){
    #   theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    #    theta2cand = runif(1,0.1,0.2)
    # }


    #if(theta3cand<0.05){
    #theta3cand =max(runif(1,0.05,0.15),rbeta_rep_ab(1,omega*muprop[3],sigmaprop[3],-1,1))
    #theta3cand = runif(1,0.1,0.2)
    #}


    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    thetacand = c(theta1cand,theta2cand,theta3cand,theta4cand)

    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))
    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))
    # a = posteriorvar3_cpp(thetacand,adjmatinf = adj_matrix,lags,X,ytot,a=c(1,1,1),b =c(1,1,1))
    # posterior_spatemp_cpp(thetacand,adjmatinf = adj_matrix,lag = N_time,X=xobs,y=y,a=c(1,1,1),b =c(1,1,1))
    #posterior_spatempver_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))

    auxnum = posterior_spatempcens_ar2_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xtot,y = ycom, a = aprior, b = bprior)
    auxden = posterior_spatempcens_ar2_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xtot,y = ycom, a = aprior, b= bprior)

    #um = auxnum$posterior
    #den = auxden$posterior

    # lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior)) #+ logpwcand
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior)) #+ logpwlast



    #  if(is.na(num)){
    #   num = runif(1)
    #}

    #if(is.na(den)){
    # den = runif(1)
    #}

    lognum = auxnum$posterior_log #+ logpwcand
    logden = auxden$posterior_log #+ logpwlast

    lognum = lognum + log(dbeta_rep_cpp(thetalast[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetalast[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_ab(thetalast[4],omega*muprop[4],sigmaprop[4],-1,1)) + log(dbeta_rep_ab(thetalast[3],omega*muprop3,sigmaprop3,thetalast[4]-1,1 - thetalast[4]))
    # (x, mu, sigma2, a, b)

    logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_ab(thetacand[4],omega*muprop[4],sigmaprop[4],-1,1)) + log(dbeta_rep_ab(thetacand[3],omega*muprop3,sigmaprop3,thetacand[4]-1,1 - thetacand[4]))
    # logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
    # log(dnorm(thetacand[3],omega*muprop[3],sigmaprop[3])) #+ logpwcand



    logunif= log(runif(1))
    n = length(y)
    p = ncol(xtot)

    if(logunif<lognum - logden)
    {
      count = count +1
      rhoF[i] <- thetacand[1]
      psiF[i] = thetacand[2]
      phi1F[i] = thetacand[3]
      phi2F[i] = thetacand[4]
      corrmat = auxnum$varcovspatemp
      betaest = auxnum$betaest
      varbeta = auxnum$varbeta
      s2est = auxnum$S2

    }else
    {
      rhoF[i] <- thetalast[1]
      psiF[i] = thetalast[2]
      phi1F[i] = thetalast[3]
      phi2F[i] = thetalast[4]
      corrmat = auxden$varcovspatemp
      betaest = auxden$betaest
      varbeta = auxden$varbeta
      s2est = auxden$S2
    }

    sigma2F[i] = rinvgamma(1,0.5*(n-p), (n-p)*s2est/2)
    betaF[i,] = rmvnorm(1,betaest,sigma2F[i]*varbeta)
    covmat = sigma2F[i] *corrmat
    ############################################################
    covobsobs = covmat[indicens ==0, indicens ==0]
    covcenscens = covmat[indicens ==1, indicens ==1]
    covcensobs = covmat[indicens ==1, indicens ==0]
    covobscens = covmat[indicens ==0, indicens ==1]
    ############################################################
    yobs = y[indicens ==0]
    xobserv = xtot[indicens ==0,]
    xcens = xtot[indicens ==1,]
    covobsobsinv = solve(covobsobs)
    mucens = xcens%*%betaF[i,] + covcensobs%*%covobsobsinv%*%(yobs - xobserv%*%betaF[i,])
    Sigmacens = covcenscens - covcensobs%*%covobsobsinv%*%covobscens
    ycensF[i,] = rtmvnorm(1,mean =as.vector(mucens), sigma = Sigmacens,lower = lower,upper = upper,algorithm = "gibbs",thinning = 4)
    ycom[indicens ==1] = ycensF[i,]
    print(c(i,sigma2F[i],rhoF[i],psiF[i],phi1F[i],phi2F[i], count/i))
    #gc()
  }

  gc()

  tau2F = sigma2F*psiF
  thetaF = cbind(betaF,sigma2F,rhoF,psiF,phi1F, phi2F, tau2F,ycensF)


  rhoburn=rhoF[(burn+1):iter]
  rhoval= rhoburn[seq(1,iter-burn,thin)]

  psiburn=psiF[(burn+1):iter]
  psival= psiburn[seq(1,iter-burn,thin)]

  phi1burn=phi1F[(burn+1):iter]
  phi1val= phi1burn[seq(1,iter-burn,thin)]

  phi2burn=phi2F[(burn+1):iter]
  phi2val= phi2burn[seq(1,iter-burn,thin)]

  sigma2Fburn=sigma2F[(burn+1):iter]
  sigma2Fval= sigma2Fburn[seq(1,iter-burn,thin)]

  tau2val = sigma2Fval*psival

  betaFburn=betaF[(burn+1):iter,]
  betaFval= betaFburn[seq(1,iter-burn,thin),]

  ycensFburn = ycensF[(burn+1):iter,]
  ycensFval = ycensFburn[seq(1,iter-burn,thin),]

  probacc = sum(count)/iter

  nmc = length(rhoval)
  thetaFval = cbind(betaFval,sigma2Fval,rhoval,psival,phi1val,phi2val,ycensFval)


  rhoest = quantile(rhoval,probs = c(0.025,0.5,0.975))
  psiest = quantile(psival,probs = c(0.025,0.5,0.975))
  tau2est = quantile(tau2val,probs = c(0.025,0.5,0.975))
  phi1est = quantile(phi1val,probs = c(0.025,0.5,0.975))
  phi2est = quantile(phi2val,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigma2Fval,probs = c(0.025,0.5,0.975))


  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)


  #### maybe later
  #log_lik = apply(thetaFval,1,loglikespatempcens_fix,Xmat = xtot,ytotobs = y,lags = lags,adjmatinf = adjmatinf,indicens = indicens)
  ######


  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  #result = list(betaF = betaFval, sigma2F = sigma2Fval, rhoF = rhoval, psiF = psival,phiF = phival,tau2F = tau2val,
  #            ycensF = ycensFval, thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
  #          psiest = psiest,phiest = phiest, tau2est = tau2est, probacc = probacc,log_lik = log_lik)

  result = list(betaF = betaFval, sigma2F = sigma2Fval, rhoF = rhoval, psiF = psival,phi1F = phi1val, phi2F = phi2val, tau2F = tau2val,
                ycensF = ycensFval, thetaF = thetaF, thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
                psiest = psiest,phi1est = phi1est, phi2est = phi2est, tau2est = tau2est, probacc = probacc)


  return(result)
}





#' @export
bayesspatempcensautosar_ar2_fix = function(y,xtot,thetaini,indicens, iter,burn,thin, lags,
                                        adjmat, aprior = c(4,2,4,4), bprior = c(2,4,2,2),
                                        divproprho = 9, divproppsi = 9, divpropphi = 9,lower,upper){

  # y= ytotobs; xtot = xtotobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;indicens = indicens3;
  #adjmatinf = adj_matrix;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time_obs




  yobs = y[indicens == 0]
  xobs = xtot[indicens ==0, ]
  rhoF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  psiF = rep(0,iter)
  phi1F = rep(0,iter)
  phi2F = rep(0,iter)
  p = ncol(xtot)
  betaF = matrix(0,iter,p)
  sigma2F = 0
  n = length(y)
  ycensF = matrix(0,iter,sum(indicens==1))
  rhoF[1] = thetaini[1]
  psiF[1] = thetaini[2]
  phi1F[1] = thetaini[3]
  phi2F[1] = thetaini[4]
  #corrmatini = spatimecovar_2_cpp(lag = lags,adjmatinf =adjmatinf, rho = rhoF[1], psi = psiF[1],
  #                              phi = phiF[1])
  #invcorrmatini = solve(corrmatini)
  #betaini = solve(t(xtot)%*%invcorrmatini%*%xtot)%*%t(xtot)%*%invcorrmatini%*%y

  #sigma2ini = (t(y-xtot%*%betaini)%*%corrmatini%*%(y-xtot%*%betaini))/(n-p)
  #betaF[1,] = betaini
  #sigma2F[1] = as.numeric(sigma2ini)

  #covmatini = sigma2F[1]*corrmatini

  ############################################################
  #covobsobsini = covmatini[indicens ==0, indicens ==0]
  #covcenscensini = covmatini[indicens ==1, indicens ==1]
  #covcensobsini = covmatini[indicens ==1, indicens ==0]
  #covobscensini = covmatini[indicens ==0, indicens ==1]
  ############################################################

  #yobs = y[indicens ==0]
  #xobserv = xtot[indicens ==0,]
  #xcens = xtot[indicens ==1,]
  #covobsobsiniinv = solve(covobsobsini)
  #mucensini = xcens%*%betaini + covcensobsini%*%covobsobsiniinv%*%(yobs - xobserv%*%betaini)
  #Sigmacensini = covcenscensini - covcensobsini%*%covobsobsiniinv%*%covobscensini
  ycensF[1,] = y[indicens ==1]
  #rtmvnorm(1,mean =as.vector(mucensini), sigma = Sigmacensini,lower = lower,upper = upper,algorithm = "gibbs",thinning = 2)
  ycom = y
  ycom[indicens ==1] = ycensF[1,]
  ycom2 = matrix(0,iter,)
 # adjmatinf[upper.tri(adjmatinf)]=0


  for(i in 2:iter){
    #i = 2
    thetalast = c(rhoF[i-1],psiF[i-1],phi1F[i-1], phi2F[i-1])

    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))

    #lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior))
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior))


    #lognum = lognum + log(dbeta(thetalast[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetalast[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetalast[3],apriorprop[3],bpriorprop[3]))

    #logden = logden + log(dbeta(thetacand[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetacand[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetacand[3],apriorprop[3],bpriorprop[3]))
    cons = 1
    muprop1 = min(0.965,thetalast[1]/cons)
    muprop2 = min(0.965,thetalast[2]/cons)
    muprop3 = min(0.965,thetalast[3]/cons)
    muprop4 = min(0.965,thetalast[4]/cons)

    muprop = c(muprop1,muprop2,muprop3,muprop4)

    sigmaprop1 = abs(muprop[1]*(1-muprop[1])/divproprho)
    sigmaprop1 = min(abs((muprop[1]-(-0.97))*(0.97-muprop[1])-0.005),sigmaprop1)
    sigmaprop1 = abs(sigmaprop1)
    sigmaprop2 = muprop[2]*(1-muprop[2])/divproppsi
    sigmaprop2 = abs(sigmaprop2)
    #sigmaprop3 = abs(sigmaprop3)
    sigmaprop4 = muprop[4]*(1-muprop[4])/divpropphi
    sigmaprop4 = abs(sigmaprop4)

    #    sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #   sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop3 = 0.00002
    sigmaprop = c(sigmaprop1,sigmaprop2,sigmaprop3,sigmaprop4)
    omega = 1
    theta1cand = min(0.99,rbeta_rep_ab(1,omega*muprop[1],sigmaprop[1],-1,1))
    theta2cand = min(0.99,rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    theta4cand = min(0.99,rbeta_rep_ab(1,omega*muprop[4],sigmaprop[4],-1,1))
    a= (theta4cand -1) + 0.03
    b = ((1-theta4cand))-0.03
    muprop3 = min(0.965,thetalast[3]/cons)
    muprop3 = ifelse(muprop[3]<=a |muprop[3]>=b, (a+b)/2,muprop3)
    sigmaprop3 = abs(muprop3*(1-muprop3)/divpropphi)
    sigmaprop3 = min(0.95*abs((muprop3-a)*(b-muprop3)),sigmaprop3)
    sigmaprop3 = abs(sigmaprop3)
    sigmaprop3 = max(sigmaprop3,0.0000005)
    print(c(sigmaprop3, (muprop3-a)*(b-muprop3)))
    theta3cand = rbeta_rep_ab(1,omega*muprop3,sigmaprop3,a,b)

    #theta3cand = rnorm(1,muprop[3],sigmaprop[3])

    #if(theta1cand > 0.95){
    # theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    #theta1cand = runif(1,0.1,0.6)
    #}

    #if(theta2cand > 0.95){
    # theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    #  theta2cand = runif(1,0.1,0.45)
    #}


    #if(theta3cand> 0.95){
    # theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
    #theta3cand = runif(1,0.1,0.6)
    #}



    #    if(theta1cand<0.05){
    #     theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    #    theta1cand = runif(1,0.1,0.2)
    #  }

    #    if(theta2cand<0.1){
    #   theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    #    theta2cand = runif(1,0.1,0.2)
    # }


    #if(theta3cand<0.05){
    #theta3cand =max(runif(1,0.05,0.15),rbeta_rep_ab(1,omega*muprop[3],sigmaprop[3],-1,1))
    #theta3cand = runif(1,0.1,0.2)
    #}


    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    thetacand = c(theta1cand,theta2cand,theta3cand,theta4cand)

    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))
    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))
    # a = posteriorvar3_cpp(thetacand,adjmatinf = adj_matrix,lags,X,ytot,a=c(1,1,1),b =c(1,1,1))
    # posterior_spatemp_cpp(thetacand,adjmatinf = adj_matrix,lag = N_time,X=xobs,y=y,a=c(1,1,1),b =c(1,1,1))
    #posterior_spatempver_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))

    auxnum = posterior_spatempcenssar_ar2_cpp(theta = thetacand,adjmat = adjmat, lag = lags, X = xtot,y = ycom, a = aprior, b = bprior)
    auxden = posterior_spatempcenssar_ar2_cpp(theta = thetalast,adjmat = adjmat,lag = lags, X = xtot,y = ycom, a = aprior, b= bprior)

    #um = auxnum$posterior
    #den = auxden$posterior

    # lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior)) #+ logpwcand
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior)) #+ logpwlast



    #  if(is.na(num)){
    #   num = runif(1)
    #}

    #if(is.na(den)){
    # den = runif(1)
    #}
    aest = thetalast[4]-1
    best = 1- thetalast[4]
    sigmaprop3est =  if_else(sigmaprop3>= (muprop3-aest)*(best-muprop3),0.000002,sigmaprop3)
    lognum = auxnum$posterior_log #+ logpwcand
    logden = auxden$posterior_log #+ logpwlast

    lognum = lognum + log(dbeta_rep_ab(thetalast[1],omega*muprop[1],sigmaprop[1],-1,1)) + log(dbeta_rep_cpp(thetalast[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_ab(thetalast[4],omega*muprop[4],sigmaprop[4],-1,1)) + log(dbeta_rep_ab(thetalast[3],omega*muprop3,sigmaprop3est,aest,best))
    # (x, mu, sigma2, a, b)

    logden = logden + log(dbeta_rep_ab(thetacand[1],omega*muprop[1],sigmaprop[1],-1,1)) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_ab(thetacand[4],omega*muprop[4],sigmaprop[4],-1,1)) + log(dbeta_rep_ab(thetacand[3],omega*muprop3,sigmaprop3,thetacand[4]-1,1 - thetacand[4]))
    # logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
    # log(dnorm(thetacand[3],omega*muprop[3],sigmaprop[3])) #+ logpwcand



    logunif= log(runif(1))
    n = length(y)
    p = ncol(xtot)

    if(logunif<lognum - logden)
    {
      count = count +1
      rhoF[i] <- thetacand[1]
      psiF[i] = thetacand[2]
      phi1F[i] = thetacand[3]
      phi2F[i] = thetacand[4]
      corrmat = auxnum$varcovspatemp
      betaest = auxnum$betaest
      varbeta = auxnum$varbeta
      s2est = auxnum$S2

    }else
    {
      rhoF[i] <- thetalast[1]
      psiF[i] = thetalast[2]
      phi1F[i] = thetalast[3]
      phi2F[i] = thetalast[4]
      corrmat = auxden$varcovspatemp
      betaest = auxden$betaest
      varbeta = auxden$varbeta
      s2est = auxden$S2
    }

    sigma2F[i] = rinvgamma(1,0.5*(n-p), (n-p)*s2est/2)
    betaF[i,] = rmvnorm(1,betaest,sigma2F[i]*varbeta)
    covmat = sigma2F[i] *corrmat
    ############################################################
    covobsobs = covmat[indicens ==0, indicens ==0]
    covcenscens = covmat[indicens ==1, indicens ==1]
    covcensobs = covmat[indicens ==1, indicens ==0]
    covobscens = covmat[indicens ==0, indicens ==1]
    ############################################################
    yobs = y[indicens ==0]
    xobserv = xtot[indicens ==0,]
    xcens = xtot[indicens ==1,]
    covobsobsinv = solve(covobsobs)
    mucens = xcens%*%betaF[i,] + covcensobs%*%covobsobsinv%*%(yobs - xobserv%*%betaF[i,])
    Sigmacens = covcenscens - covcensobs%*%covobsobsinv%*%covobscens
    ycensF[i,] = rtmvnorm(1,mean =as.vector(mucens), sigma = Sigmacens,lower = lower,upper = upper,algorithm = "gibbs",thinning = 4)
    ycom[indicens ==1] = ycensF[i,]
    print(c(i,sigma2F[i],rhoF[i],psiF[i],phi1F[i],phi2F[i], count/i))
    #gc()
  }

  gc()

  tau2F = sigma2F*psiF
  thetaF = cbind(betaF,sigma2F,rhoF,psiF,phi1F, phi2F,tau2F,ycensF)


  rhoburn=rhoF[(burn+1):iter]
  rhoval= rhoburn[seq(1,iter-burn,thin)]

  psiburn=psiF[(burn+1):iter]
  psival= psiburn[seq(1,iter-burn,thin)]

  phi1burn=phi1F[(burn+1):iter]
  phi1val= phi1burn[seq(1,iter-burn,thin)]

  phi2burn=phi2F[(burn+1):iter]
  phi2val= phi2burn[seq(1,iter-burn,thin)]

  sigma2Fburn=sigma2F[(burn+1):iter]
  sigma2Fval= sigma2Fburn[seq(1,iter-burn,thin)]

  tau2val = sigma2Fval*psival

  betaFburn=betaF[(burn+1):iter,]
  betaFval= betaFburn[seq(1,iter-burn,thin),]

  ycensFburn = ycensF[(burn+1):iter,]
  ycensFval = ycensFburn[seq(1,iter-burn,thin),]

  probacc = sum(count)/iter

  nmc = length(rhoval)
  thetaFval = cbind(betaFval,sigma2Fval,rhoval,psival,phi1val,phi2val,ycensFval)


  rhoest = quantile(rhoval,probs = c(0.025,0.5,0.975))
  psiest = quantile(psival,probs = c(0.025,0.5,0.975))
  tau2est = quantile(tau2val,probs = c(0.025,0.5,0.975))
  phi1est = quantile(phi1val,probs = c(0.025,0.5,0.975))
  phi2est = quantile(phi2val,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigma2Fval,probs = c(0.025,0.5,0.975))


  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)


  #### maybe later
  #log_lik = apply(thetaFval,1,loglikespatempcens_fix,Xmat = xtot,ytotobs = y,lags = lags,adjmatinf = adjmatinf,indicens = indicens)
  ######


  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  #result = list(betaF = betaFval, sigma2F = sigma2Fval, rhoF = rhoval, psiF = psival,phiF = phival,tau2F = tau2val,
  #            ycensF = ycensFval, thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
  #          psiest = psiest,phiest = phiest, tau2est = tau2est, probacc = probacc,log_lik = log_lik)

  result = list(betaF = betaFval, sigma2F = sigma2Fval, rhoF = rhoval, psiF = psival,phi1F = phi1val, phi2F = phi2val, tau2F = tau2val,
                ycensF = ycensFval, thetaF = thetaF,thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
                psiest = psiest,phi1est = phi1est, phi2est = phi2est, tau2est = tau2est, probacc = probacc)


  return(result)
}















#' @export
bayesspatempcensautocar_fix = function(y,xtot,thetaini,indicens, iter,burn,thin, lags,
                                    adjmat, aprior = c(4,2,4), bprior = c(2,4,2),
                                    divproprho = 9, divproppsi = 9, divpropphi = 9,lower,upper){

  # y= ytotobs; xtot = xtotobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;indicens = indicens3;
  #adjmatinf = adj_matrix;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time_obs




  yobs = y[indicens == 0]
  xobs = xtot[indicens ==0, ]
  rhoF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  psiF = rep(0,iter)
  phiF = rep(0,iter)
  p = ncol(xtot)
  betaF = matrix(0,iter,p)
  sigma2F = 0
  n = length(y)
  ycensF = matrix(0,iter,sum(indicens==1))
  rhoF[1] = thetaini[1]
  psiF[1] = thetaini[2]
  phiF[1] = thetaini[3]
  #corrmatini = spatimecovar_2_cpp(lag = lags,adjmatinf =adjmatinf, rho = rhoF[1], psi = psiF[1],
  #                              phi = phiF[1])
  #invcorrmatini = solve(corrmatini)
  #betaini = solve(t(xtot)%*%invcorrmatini%*%xtot)%*%t(xtot)%*%invcorrmatini%*%y

  #sigma2ini = (t(y-xtot%*%betaini)%*%corrmatini%*%(y-xtot%*%betaini))/(n-p)
  #betaF[1,] = betaini
  #sigma2F[1] = as.numeric(sigma2ini)

  #covmatini = sigma2F[1]*corrmatini

  ############################################################
  #covobsobsini = covmatini[indicens ==0, indicens ==0]
  #covcenscensini = covmatini[indicens ==1, indicens ==1]
  #covcensobsini = covmatini[indicens ==1, indicens ==0]
  #covobscensini = covmatini[indicens ==0, indicens ==1]
  ############################################################

  #yobs = y[indicens ==0]
  #xobserv = xtot[indicens ==0,]
  #xcens = xtot[indicens ==1,]
  #covobsobsiniinv = solve(covobsobsini)
  #mucensini = xcens%*%betaini + covcensobsini%*%covobsobsiniinv%*%(yobs - xobserv%*%betaini)
  #Sigmacensini = covcenscensini - covcensobsini%*%covobsobsiniinv%*%covobscensini
  ycensF[1,] = y[indicens ==1]
  #rtmvnorm(1,mean =as.vector(mucensini), sigma = Sigmacensini,lower = lower,upper = upper,algorithm = "gibbs",thinning = 2)
  ycom = y
  ycom[indicens ==1] = ycensF[1,]
  ycom2 = matrix(0,iter,)
 ## adjmatinf[upper.tri(adjmatinf)]=0
  adjmatinf = adjmat


  for(i in 2:iter){
    #i = 2
    thetalast = c(rhoF[i-1],psiF[i-1],phiF[i-1])

    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))

    #lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior))
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior))


    #lognum = lognum + log(dbeta(thetalast[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetalast[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetalast[3],apriorprop[3],bpriorprop[3]))

    #logden = logden + log(dbeta(thetacand[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetacand[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetacand[3],apriorprop[3],bpriorprop[3]))
    cons = 1
    muprop1 = min(0.99,thetalast[1]/cons)
    muprop2 = min(0.99,thetalast[2]/cons)
    muprop3 = min(0.99,thetalast[3]/cons)

    muprop = c(muprop1,muprop2,muprop3)
    sigmaprop1 = muprop[1]*(1-muprop[1])/divproprho
    sigmaprop2 = muprop[2]*(1-muprop[2])/divproppsi
    sigmaprop3 = muprop[3]*(1-muprop[3])/divpropphi

    #    sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #   sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop = c(sigmaprop1,sigmaprop2,sigmaprop3)
    omega = runif(1,0.99,1)
    theta1cand = min(0.99,rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    theta2cand = min(0.99,rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    theta3cand = min(0.99,rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))


    #if(theta1cand > 0.95){
    # theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    #theta1cand = runif(1,0.1,0.6)
    #}

    #if(theta2cand > 0.95){
    # theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    #  theta2cand = runif(1,0.1,0.45)
    #}


    #if(theta3cand> 0.95){
    # theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
    #theta3cand = runif(1,0.1,0.6)
    #}



    if(theta1cand<0.05){
      theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
      theta1cand = runif(1,0.1,0.2)
    }

    if(theta2cand<0.05){
      theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
      theta2cand = runif(1,0.1,0.2)
    }


    if(theta3cand<0.05){
      theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
      theta3cand = runif(1,0.1,0.2)
    }


    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    thetacand = c(theta1cand,theta2cand,theta3cand)

    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))
    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))
    # a = posteriorvar3_cpp(thetacand,adjmatinf = adj_matrix,lags,X,ytot,a=c(1,1,1),b =c(1,1,1))
    # posterior_spatemp_cpp(thetacand,adjmatinf = adj_matrix,lag = N_time,X=xobs,y=y,a=c(1,1,1),b =c(1,1,1))
    #posterior_spatempver_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))

    auxnum = posterior_spatempcenscar_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xtot,y = ycom, a = aprior, b = bprior)
    auxden = posterior_spatempcenscar_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xtot,y = ycom, a = aprior, b= bprior)

    #num = auxnum$posterior
    #den = auxden$posterior

    # lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior)) #+ logpwcand
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior)) #+ logpwlast



    #  if(is.na(num)){
    #   num = runif(1)
    #}

    #if(is.na(den)){
    # den = runif(1)
    #}

    lognum = auxnum$posterior_log #+ logpwcand
    logden = auxden$posterior_log #+ logpwlast

    lognum = lognum + log(dbeta_rep_cpp(thetalast[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetalast[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_cpp(thetalast[3],omega*muprop[3],sigmaprop[3])) #+ logpwlast

    logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2])) +
      log(dbeta_rep_cpp(thetacand[3],omega*muprop[3],sigmaprop[3])) #+ logpwcand



    logunif= log(runif(1))
    n = length(y)
    p = ncol(xtot)

    if(logunif<lognum - logden)
    {
      count = count +1
      rhoF[i] <- thetacand[1]
      psiF[i] = thetacand[2]
      phiF[i] = thetacand[3]
      corrmat = auxnum$varcovspatemp
      betaest = auxnum$betaest
      varbeta = auxnum$varbeta
      s2est = auxnum$S2

    }else
    {
      rhoF[i] <- thetalast[1]
      psiF[i] = thetalast[2]
      phiF[i] = thetalast[3]
      corrmat = auxden$varcovspatemp
      betaest = auxden$betaest
      varbeta = auxden$varbeta
      s2est = auxden$S2
    }

    sigma2F[i] = rinvgamma(1,0.5*(n-p), (n-p)*s2est/2)
    betaF[i,] = rmvnorm(1,betaest,sigma2F[i]*varbeta)
    covmat = sigma2F[i] *corrmat
    ############################################################
    covobsobs = covmat[indicens ==0, indicens ==0]
    covcenscens = covmat[indicens ==1, indicens ==1]
    covcensobs = covmat[indicens ==1, indicens ==0]
    covobscens = covmat[indicens ==0, indicens ==1]
    ############################################################
    yobs = y[indicens ==0]
    xobserv = xtot[indicens ==0,]
    xcens = xtot[indicens ==1,]
    covobsobsinv = solve(covobsobs)
    mucens = xcens%*%betaF[i,] + covcensobs%*%covobsobsinv%*%(yobs - xobserv%*%betaF[i,])
    Sigmacens = covcenscens - covcensobs%*%covobsobsinv%*%covobscens
    ycensF[i,] = rtmvnorm(1,mean =as.vector(mucens), sigma = Sigmacens,lower = lower,upper = upper,algorithm = "gibbs")
    ycom[indicens ==1] = ycensF[i,]
    print(c(i,sigma2F[i],rhoF[i],psiF[i],phiF[i],count/i))
    #gc()
  }

  gc()

  rhoburn=rhoF[(burn+1):iter]
  rhoval= rhoburn[seq((burn+1),iter-burn,thin)]

  psiburn=psiF[(burn+1):iter]
  psival= psiburn[seq((burn+1),iter-burn,thin)]

  phiburn=phiF[(burn+1):iter]
  phival= phiburn[seq((burn+1),iter-burn,thin)]

  sigma2Fburn=sigma2F[(burn+1):iter]
  sigma2Fval= sigma2Fburn[seq((burn+1),iter-burn,thin)]

  tau2val = sigma2Fval*psival

  betaFburn=betaF[(burn+1):iter,]
  betaFval= betaFburn[seq((burn+1),iter-burn,thin),]

  ycensFburn = ycensF[(burn+1):iter,]
  ycensFval = ycensFburn[seq((burn+1),iter-burn,thin),]

  probacc = sum(count)/iter

  nmc = length(rhoval)
  thetaFval = cbind(betaFval,sigma2Fval,rhoval,psival,phival,ycensFval)


  rhoest = quantile(rhoval,probs = c(0.025,0.5,0.975))
  psiest = quantile(psival,probs = c(0.025,0.5,0.975))
  tau2est = quantile(tau2val,probs = c(0.025,0.5,0.975))
  phiest = quantile(phival,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigma2Fval,probs = c(0.025,0.5,0.975))


  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)


  #### maybe later
  log_lik = apply(thetaFval,1,loglikespatempcens_fix,Xmat = xtot,ytotobs = y,lags = lags,adjmatinf = adjmatinf,indicens = indicens)
  ######


  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  result = list(betaF = betaFval, sigma2F = sigma2Fval, rhoF = rhoval, psiF = psival,phiF = phival,tau2F = tau2val,
                ycensF = ycensFval, thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
                psiest = psiest,phiest = phiest, tau2est = tau2est, probacc = probacc,log_lik = log_lik)

  return(result)
}






#' @export
bayesspatempcensautocar_2_fix = function(y,xtot,thetaini,indicens, iter,burn,thin, lags,
                                       adjmat, aprior = c(4,2,4,4), bprior = c(2,4,2,4),
                                       alim = c(-10,-10,-10,0), blim = c(10,10,10,1),
                                       divproprho_s = 9, divproprho_t = 9, divproprho_st = 9, divproppsi = 9, lower,upper){

  # y= ytotobs; xtot = xtotobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;indicens = indicens3;
  #adjmatinf = adj_matrix;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time_obs




  yobs = y[indicens == 0]
  xobs = xtot[indicens ==0, ]
  rho_sF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  rho_tF = rep(0,iter)
  rho_stF = rep(0,iter)
  psiF = rep(0,iter)
  p = ncol(xtot)
  betaF = matrix(0,iter,p)
  sigma2F = 0
  n = length(y)
  ycensF = matrix(0,iter,sum(indicens==1))
  rho_sF[1] = thetaini[1]
  rho_tF[1] = thetaini[2]
  rho_stF[1] = thetaini[3]
  psiF[1] = thetaini[4]
  #corrmatini = spatimecovar_2_cpp(lag = lags,adjmatinf =adjmatinf, rho = rhoF[1], psi = psiF[1],
  #                              phi = phiF[1])
  #invcorrmatini = solve(corrmatini)
  #betaini = solve(t(xtot)%*%invcorrmatini%*%xtot)%*%t(xtot)%*%invcorrmatini%*%y

  #sigma2ini = (t(y-xtot%*%betaini)%*%corrmatini%*%(y-xtot%*%betaini))/(n-p)
  #betaF[1,] = betaini
  #sigma2F[1] = as.numeric(sigma2ini)

  #covmatini = sigma2F[1]*corrmatini

  ############################################################
  #covobsobsini = covmatini[indicens ==0, indicens ==0]
  #covcenscensini = covmatini[indicens ==1, indicens ==1]
  #covcensobsini = covmatini[indicens ==1, indicens ==0]
  #covobscensini = covmatini[indicens ==0, indicens ==1]
  ############################################################

  #yobs = y[indicens ==0]
  #xobserv = xtot[indicens ==0,]
  #xcens = xtot[indicens ==1,]
  #covobsobsiniinv = solve(covobsobsini)
  #mucensini = xcens%*%betaini + covcensobsini%*%covobsobsiniinv%*%(yobs - xobserv%*%betaini)
  #Sigmacensini = covcenscensini - covcensobsini%*%covobsobsiniinv%*%covobscensini
  ycensF[1,] = y[indicens ==1]
  #rtmvnorm(1,mean =as.vector(mucensini), sigma = Sigmacensini,lower = lower,upper = upper,algorithm = "gibbs",thinning = 2)
  ycom = y
  ycom[indicens ==1] = ycensF[1,]
  ycom2 = matrix(0,iter,)
  ## adjmatinf[upper.tri(adjmatinf)]=0
  adjmatinf = adjmat
  Wt <- abs(outer(1:lags, 1:lags, "-")) == 1
  Wt = 1*Wt


  for(i in 2:iter){
    #i = 2
    thetalast = c(rho_sF[i-1],rho_tF[i-1],rho_stF[i-1], psiF[i-1])

    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))

    #lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior))
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior))


    #lognum = lognum + log(dbeta(thetalast[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetalast[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetalast[3],apriorprop[3],bpriorprop[3]))

    #logden = logden + log(dbeta(thetacand[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetacand[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetacand[3],apriorprop[3],bpriorprop[3]))
    cons = 1
    muprop1 = min(0.99,thetalast[1]/cons)
    muprop2 = min(0.99,thetalast[2]/cons)
    muprop3 = min(0.99,thetalast[3]/cons)
    muprop4 = min(0.99,thetalast[4]/cons)

    muprop = c(muprop1,muprop2,muprop3,muprop4)

    #sigmaprop1 = abs(muprop[1]*(1-muprop[1]))/divproprho_s
    #sigmaprop1 = min(abs((muprop[1]-alim[1])*(blim[1] - muprop[1])+0.000015),sigmaprop1)
    #sigmaprop1 = if_else(sigmaprop1< (muprop[1]-alim[1])*(blim[1] - muprop[1]), sigmaprop1,0.000005 )

    #sigmaprop2 = abs(muprop[2]*(1-abs(muprop[2])))/divproprho_t
    sigmaprop1 = abs((muprop[1]-alim[1])*(blim[1] - muprop[1]))/divproprho_s
    sigmaprop2 = abs((muprop[2]-alim[2])*(blim[2] - muprop[2]))/divproprho_t
    sigmaprop3 = abs((muprop[3]-alim[3])*(blim[3] - muprop[3]))/divproprho_st

    #sigmaprop2 = min(abs((muprop[2]-alim[2])*(blim[2] - muprop[2])+0.000015),sigmaprop2)
    #sigmaprop2 = if_else(sigmaprop2< (muprop[2]-alim[2])*(blim[2] - muprop[2]), sigmaprop2,0.000005 )
    #sigmaprop2 = min((muprop[2]-alim[2])*(blim[2] - muprop[2])+0.000015,sigmaprop2)
    #sigmaprop2 = abs(sigmaprop2)
    #sigmaprop3 = abs(muprop[3]*(1-muprop[3]))/divproprho_st
    #sigmaprop3 = min(abs((muprop[3]-alim[3])*(blim[3] - muprop[3])+0.000015),sigmaprop3)

    #sigmaprop3 = if_else(sigmaprop3< (muprop[3]-alim[3])*(blim[3] - muprop[3]), sigmaprop3,0.000005 )

    sigmaprop4 = abs(muprop[4]*(1-muprop[4]))/divproppsi

    print(c(sigmaprop2, (muprop[2]-alim[2])*(blim[2] - muprop[2])))
    #    sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #   sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop = c(sigmaprop1,sigmaprop2,sigmaprop3, sigmaprop4)
    omega = runif(1,0.99,1)
    theta1cand = rbeta_rep_ab(1,muprop[1],sigmaprop[1],alim[1],blim[1])
    theta2cand = rbeta_rep_ab(1,muprop[2],sigmaprop[2],alim[2],blim[2])
    theta3cand = rbeta_rep_ab(1,muprop[3],sigmaprop[3],alim[3],blim[3])
    theta4cand = rbeta_rep_ab(1,muprop[4],sigmaprop[4],alim[4],blim[4])


    if(theta1cand <= alim[1]){
    theta1cand =max(alim[1] +0.00005 ,rbeta_rep_ab(1,muprop[1],sigmaprop[1],alim[1],blim[1]))
   # theta1cand = runif(1,0.1,0.6)
    }

    if(theta2cand <= alim[2]){
      theta2cand = max(alim[2] +0.00005 ,rbeta_rep_ab(1,muprop[2],sigmaprop[2],alim[2],blim[2]))
      # theta1cand = runif(1,0.1,0.6)
    }


    if(theta3cand <= alim[3]){
      theta3cand =max(alim[3] +0.005 ,rbeta_rep_ab(1,muprop[3],sigmaprop[3],alim[3],blim[3]))
      # theta1cand = runif(1,0.1,0.6)
    }

    if(theta4cand <= alim[4]){
      theta4cand =max(alim[4] +0.005 ,rbeta_rep_ab(1,muprop[4],sigmaprop[4],alim[4],blim[4]))
      # theta1cand = runif(1,0.1,0.6)
    }


    if(theta1cand >= blim[1]){
      theta1cand =min(blim[1] - 0.005 ,rbeta_rep_ab(1,muprop[1],sigmaprop[1],alim[1],blim[1]))
      # theta1cand = runif(1,0.1,0.6)
    }

    if(theta2cand >= blim[2]){
      theta2cand =min(blim[2] -0.00005 ,rbeta_rep_ab(1,muprop[2],sigmaprop[2],alim[2],blim[2]))
      # theta1cand = runif(1,0.1,0.6)
    }


    if(theta3cand >= blim[3]){
      theta3cand =min(blim[3] - 0.00005 ,rbeta_rep_ab(1,muprop[3],sigmaprop[3],alim[3],blim[3]))
      # theta1cand = runif(1,0.1,0.6)
    }

    if(theta4cand >= blim[4]){
      theta4cand =min(blim[4] -0.00005 ,rbeta_rep_ab(1,muprop[4],sigmaprop[4],alim[4],blim[4]))
      # theta1cand = runif(1,0.1,0.6)
    }

    #if(theta3cand> 0.95){
    # theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
    #theta3cand = runif(1,0.1,0.6)
    #}



    #if(theta1cand<0.05){
     # theta1cand =max(runif(1,0.05,0.15),rbeta_rep_ab(1,omega*muprop[1],sigmaprop[1],alim[1],blim[1]))
    #  theta1cand = runif(1,0.1,0.2)
    #}

    #if(theta2cand<0.05){
     # theta2cand =max(runif(1,0.05,0.15),rbeta_rep_ab(1,omega*muprop[2],sigmaprop[2],alim[2],blim[2]))
    #  theta2cand = runif(1,0.1,0.2)
    #}


    #if(theta3cand<0.05){
     # theta3cand =max(runif(1,0.05,0.15),rbeta_rep_ab(1,omega*muprop[3],sigmaprop[3],alim[3],blim[3]))
    #  theta3cand = runif(1,0.1,0.2)
    #}

    #if(theta4cand<0.05){
     # theta4cand =max(runif(1,0.05,0.15),rbeta_rep_ab(1,omega*muprop[4],sigmaprop[4],alim[4],blim[4]))
    #  theta4cand = runif(1,0.1,0.2)
    #}

    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    thetacand = c(theta1cand,theta2cand,theta3cand, theta4cand)

    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))
    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))
    # a = posteriorvar3_cpp(thetacand,adjmatinf = adj_matrix,lags,X,ytot,a=c(1,1,1),b =c(1,1,1))
    # posterior_spatemp_cpp(thetacand,adjmatinf = adj_matrix,lag = N_time,X=xobs,y=y,a=c(1,1,1),b =c(1,1,1))
    #posterior_spatempver_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))

    auxnum = posterior_spatempcenscar2_cpp(theta = thetacand,adjmatinf = adjmatinf, Wt = Wt, X = xtot,y = ycom, a = aprior, b = bprior, alim = alim, blim = blim)
    auxden = posterior_spatempcenscar2_cpp(theta = thetalast,adjmatinf = adjmatinf, Wt = Wt, X = xtot,y = ycom, a = aprior, b= bprior, alim = alim, blim = blim)

    #num = auxnum$posterior
    #den = auxden$posterior

    # lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior)) #+ logpwcand
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior)) #+ logpwlast



    #  if(is.na(num)){
    #   num = runif(1)
    #}

    #if(is.na(den)){
    # den = runif(1)
    #}

    lognum = auxnum$posterior_log #+ logpwcand
    logden = auxden$posterior_log #+ logpwlast

    lognum = lognum + log(dbeta_rep_ab(thetalast[1],muprop[1],sigmaprop[1],alim[1],blim[1]))+
                      log(dbeta_rep_ab(thetalast[2],muprop[2],sigmaprop[2],alim[2],blim[2])) +
                      log(dbeta_rep_ab(thetalast[3],muprop[3],sigmaprop[3],alim[3],blim[3])) +
                      log(dbeta_rep_ab(thetalast[4],muprop[4],sigmaprop[4],alim[4],blim[4]))

    logden = logden + log(dbeta_rep_ab(thetacand[1],muprop[1],sigmaprop[1], alim[1],blim[1])) +
                    + log(dbeta_rep_ab(thetacand[2],muprop[2],sigmaprop[2], alim[2],blim[2])) +
                    + log(dbeta_rep_ab(thetacand[3],muprop[3],sigmaprop[3], alim[3],blim[3])) +
                    + log(dbeta_rep_ab(thetacand[4],muprop[4],sigmaprop[4], alim[4],blim[4]))

    logunif= log(runif(1))
    n = length(y)
    p = ncol(xtot)

    if(logunif<lognum - logden)
    {
      count = count +1
      rho_sF[i] <- thetacand[1]
      rho_tF[i] = thetacand[2]
      rho_stF[i] = thetacand[3]
      psiF[i] = thetacand[4]
      corrmat = auxnum$varcovspatemp
      betaest = auxnum$betaest
      varbeta = auxnum$varbeta
      s2est = auxnum$S2

    }else
    {
      rho_sF[i] =  thetalast[1]
      rho_tF[i] =  thetalast[2]
      rho_stF[i] = thetalast[3]
      psiF[i] = thetalast[4]

      corrmat = auxden$varcovspatemp
      betaest = auxden$betaest
      varbeta = auxden$varbeta
      s2est = auxden$S2
    }

    sigma2F[i] = rinvgamma(1,0.5*(n-p), (n-p)*s2est/2)
    betaF[i,] = rmvnorm(1,betaest,sigma2F[i]*varbeta)
    covmat = sigma2F[i] *corrmat
    ############################################################
    covobsobs = covmat[indicens ==0, indicens ==0]
    covcenscens = covmat[indicens ==1, indicens ==1]
    covcensobs = covmat[indicens ==1, indicens ==0]
    covobscens = covmat[indicens ==0, indicens ==1]
    ############################################################
    yobs = y[indicens ==0]
    xobserv = xtot[indicens ==0,]
    xcens = xtot[indicens ==1,]
    covobsobsinv = solve(covobsobs)
    mucens = xcens%*%betaF[i,] + covcensobs%*%covobsobsinv%*%(yobs - xobserv%*%betaF[i,])
    Sigmacens = covcenscens - covcensobs%*%covobsobsinv%*%covobscens
    ycensF[i,] = rtmvnorm(1,mean =as.vector(mucens), sigma = Sigmacens,lower = lower,upper = upper,algorithm = "gibbs")
    ycom[indicens ==1] = ycensF[i,]
    print(c(i,sigma2F[i],rho_sF[i],rho_tF[i],rho_stF[i], psiF[i], count/i))
    #gc()
  }

  gc()

  rho_s_burn=rho_sF[(burn+1):iter]
  rho_s_val= rho_s_burn[seq((burn+1),iter-burn,thin)]

  rho_t_burn=rho_tF[(burn+1):iter]
  rho_t_val= rho_t_burn[seq((burn+1),iter-burn,thin)]

  rho_st_burn=rho_stF[(burn+1):iter]
  rho_st_val= rho_st_burn[seq((burn+1),iter-burn,thin)]

  psiburn= psiF[(burn+1):iter]
  psival=  psiburn[seq((burn+1),iter-burn,thin)]


  sigma2Fburn=sigma2F[(burn+1):iter]
  sigma2Fval= sigma2Fburn[seq((burn+1),iter-burn,thin)]

  tau2val = sigma2Fval*psival

  betaFburn=betaF[(burn+1):iter,]
  betaFval= betaFburn[seq((burn+1),iter-burn,thin),]

  ycensFburn = ycensF[(burn+1):iter,]
  ycensFval = ycensFburn[seq((burn+1),iter-burn,thin),]

  probacc = sum(count)/iter

  nmc = length(rho_s_val)
  thetaFval = cbind(betaFval,sigma2Fval,rho_s_val,rho_t_val,rho_st_val,psival,ycensFval)


  rho_s_est = quantile(rho_s_val,probs = c(0.025,0.5,0.975))
  rho_t_est = quantile(rho_t_val,probs = c(0.025,0.5,0.975))
  rho_st_est = quantile(rho_st_val,probs = c(0.025,0.5,0.975))
  psiest = quantile(psival,probs = c(0.025,0.5,0.975))
  tau2est = quantile(tau2val,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigma2Fval,probs = c(0.025,0.5,0.975))


  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)


  #### maybe later
 # log_lik = apply(thetaFval,1,loglikespatempcens_fix,Xmat = xtot,ytotobs = y,lags = lags,adjmatinf = adjmatinf,indicens = indicens)
  ######


  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  result = list(betaF = betaFval, sigma2F = sigma2Fval, rho_sF = rho_s_val, rho_tF = rho_t_val,rho_stF = rho_st_val, psiF = psival, tau2F = tau2val,
                ycensF = ycensFval, thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rho_s_est = rho_s_est,
                rho_t_est = rho_t_est, rho_st_est = rho_st_est, psiest = psiest, tau2est = tau2est, probacc = probacc)

  return(result)
}






#' @export
bayesspatempcensautocar_zero_fix = function(y,xtot,thetaini,indicens, iter,burn,thin, lags,
                                         adjmat, aprior = c(4,2,4,4), bprior = c(2,4,2,4),
                                         alim = c(-10,-10,-10,0), blim = c(10,10,10,1),
                                         divproprho_s = 9, divproprho_t = 9, divproprho_st = 9, lower,upper){

  # y= ytotobs; xtot = xtotobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;indicens = indicens3;
  #adjmatinf = adj_matrix;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time_obs




  yobs = y[indicens == 0]
  xobs = xtot[indicens ==0, ]
  rho_sF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  rho_tF = rep(0,iter)
  rho_stF = rep(0,iter)
  p = ncol(xtot)
  betaF = matrix(0,iter,p)
  sigma2F = 0
  n = length(y)
  ycensF = matrix(0,iter,sum(indicens==1))
  rho_sF[1] = thetaini[1]
  rho_tF[1] = thetaini[2]
  rho_stF[1] = thetaini[3]
  #corrmatini = spatimecovar_2_cpp(lag = lags,adjmatinf =adjmatinf, rho = rhoF[1], psi = psiF[1],
  #                              phi = phiF[1])
  #invcorrmatini = solve(corrmatini)
  #betaini = solve(t(xtot)%*%invcorrmatini%*%xtot)%*%t(xtot)%*%invcorrmatini%*%y

  #sigma2ini = (t(y-xtot%*%betaini)%*%corrmatini%*%(y-xtot%*%betaini))/(n-p)
  #betaF[1,] = betaini
  #sigma2F[1] = as.numeric(sigma2ini)

  #covmatini = sigma2F[1]*corrmatini

  ############################################################
  #covobsobsini = covmatini[indicens ==0, indicens ==0]
  #covcenscensini = covmatini[indicens ==1, indicens ==1]
  #covcensobsini = covmatini[indicens ==1, indicens ==0]
  #covobscensini = covmatini[indicens ==0, indicens ==1]
  ############################################################

  #yobs = y[indicens ==0]
  #xobserv = xtot[indicens ==0,]
  #xcens = xtot[indicens ==1,]
  #covobsobsiniinv = solve(covobsobsini)
  #mucensini = xcens%*%betaini + covcensobsini%*%covobsobsiniinv%*%(yobs - xobserv%*%betaini)
  #Sigmacensini = covcenscensini - covcensobsini%*%covobsobsiniinv%*%covobscensini
  ycensF[1,] = y[indicens ==1]
  #rtmvnorm(1,mean =as.vector(mucensini), sigma = Sigmacensini,lower = lower,upper = upper,algorithm = "gibbs",thinning = 2)
  ycom = y
  ycom[indicens ==1] = ycensF[1,]
  ycom2 = matrix(0,iter,)
  ## adjmatinf[upper.tri(adjmatinf)]=0
  adjmatinf = adjmat
  Wt <- abs(outer(1:lags, 1:lags, "-")) == 1
  Wt = 1*Wt


  for(i in 2:iter){
    #i = 2
    thetalast = c(rho_sF[i-1],rho_tF[i-1],rho_stF[i-1])

    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))

    #lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior))
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior))


    #lognum = lognum + log(dbeta(thetalast[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetalast[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetalast[3],apriorprop[3],bpriorprop[3]))

    #logden = logden + log(dbeta(thetacand[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetacand[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetacand[3],apriorprop[3],bpriorprop[3]))
    cons = 1
    muprop1 = min(0.99,thetalast[1]/cons)
    muprop2 = min(0.99,thetalast[2]/cons)
    muprop3 = min(0.99,thetalast[3]/cons)

    muprop = c(muprop1,muprop2,muprop3)

    #sigmaprop1 = abs(muprop[1]*(1-muprop[1]))/divproprho_s
    #sigmaprop1 = min(abs((muprop[1]-alim[1])*(blim[1] - muprop[1])+0.000015),sigmaprop1)
    #sigmaprop1 = if_else(sigmaprop1< (muprop[1]-alim[1])*(blim[1] - muprop[1]), sigmaprop1,0.000005 )

    #sigmaprop2 = abs(muprop[2]*(1-abs(muprop[2])))/divproprho_t
    sigmaprop1 = abs((muprop[1]-alim[1])*(blim[1] - muprop[1]))/divproprho_s
    sigmaprop2 = abs((muprop[2]-alim[2])*(blim[2] - muprop[2]))/divproprho_t
    sigmaprop3 = abs((muprop[3]-alim[3])*(blim[3] - muprop[3]))/divproprho_st

    #sigmaprop2 = min(abs((muprop[2]-alim[2])*(blim[2] - muprop[2])+0.000015),sigmaprop2)
    #sigmaprop2 = if_else(sigmaprop2< (muprop[2]-alim[2])*(blim[2] - muprop[2]), sigmaprop2,0.000005 )
    #sigmaprop2 = min((muprop[2]-alim[2])*(blim[2] - muprop[2])+0.000015,sigmaprop2)
    #sigmaprop2 = abs(sigmaprop2)
    #sigmaprop3 = abs(muprop[3]*(1-muprop[3]))/divproprho_st
    #sigmaprop3 = min(abs((muprop[3]-alim[3])*(blim[3] - muprop[3])+0.000015),sigmaprop3)

    #sigmaprop3 = if_else(sigmaprop3< (muprop[3]-alim[3])*(blim[3] - muprop[3]), sigmaprop3,0.000005 )


    print(c(sigmaprop2, (muprop[2]-alim[2])*(blim[2] - muprop[2])))
    #    sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #   sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop = c(sigmaprop1,sigmaprop2,sigmaprop3)
    omega = runif(1,0.99,1)
    theta1cand = rbeta_rep_ab(1,muprop[1],sigmaprop[1],alim[1],blim[1])
    theta2cand = rbeta_rep_ab(1,muprop[2],sigmaprop[2],alim[2],blim[2])
    theta3cand = rbeta_rep_ab(1,muprop[3],sigmaprop[3],alim[3],blim[3])


    if(theta1cand <= alim[1]){
      theta1cand =max(alim[1] +0.00005 ,rbeta_rep_ab(1,muprop[1],sigmaprop[1],alim[1],blim[1]))
      # theta1cand = runif(1,0.1,0.6)
    }

    if(theta2cand <= alim[2]){
      theta2cand = max(alim[2] +0.00005 ,rbeta_rep_ab(1,muprop[2],sigmaprop[2],alim[2],blim[2]))
      # theta1cand = runif(1,0.1,0.6)
    }


    if(theta3cand <= alim[3]){
      theta3cand =max(alim[3] +0.00005 ,rbeta_rep_ab(1,muprop[3],sigmaprop[3],alim[3],blim[3]))
      # theta1cand = runif(1,0.1,0.6)
    }


    if(theta1cand >= blim[1]){
      theta1cand =min(blim[1] - 0.00005 ,rbeta_rep_ab(1,muprop[1],sigmaprop[1],alim[1],blim[1]))
      # theta1cand = runif(1,0.1,0.6)
    }

    if(theta2cand >= blim[2]){
      theta2cand =min(blim[2] -0.00005 ,rbeta_rep_ab(1,muprop[2],sigmaprop[2],alim[2],blim[2]))
      # theta1cand = runif(1,0.1,0.6)
    }


    if(theta3cand >= blim[3]){
      theta3cand =min(blim[3] - 0.00005 ,rbeta_rep_ab(1,muprop[3],sigmaprop[3],alim[3],blim[3]))
      # theta1cand = runif(1,0.1,0.6)
    }

    #if(theta3cand> 0.95){
    # theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
    #theta3cand = runif(1,0.1,0.6)
    #}



    #if(theta1cand<0.05){
    # theta1cand =max(runif(1,0.05,0.15),rbeta_rep_ab(1,omega*muprop[1],sigmaprop[1],alim[1],blim[1]))
    #  theta1cand = runif(1,0.1,0.2)
    #}

    #if(theta2cand<0.05){
    # theta2cand =max(runif(1,0.05,0.15),rbeta_rep_ab(1,omega*muprop[2],sigmaprop[2],alim[2],blim[2]))
    #  theta2cand = runif(1,0.1,0.2)
    #}


    #if(theta3cand<0.05){
    # theta3cand =max(runif(1,0.05,0.15),rbeta_rep_ab(1,omega*muprop[3],sigmaprop[3],alim[3],blim[3]))
    #  theta3cand = runif(1,0.1,0.2)
    #}

    #if(theta4cand<0.05){
    # theta4cand =max(runif(1,0.05,0.15),rbeta_rep_ab(1,omega*muprop[4],sigmaprop[4],alim[4],blim[4]))
    #  theta4cand = runif(1,0.1,0.2)
    #}

    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    thetacand = c(theta1cand,theta2cand,theta3cand)

    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))
    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))
    # a = posteriorvar3_cpp(thetacand,adjmatinf = adj_matrix,lags,X,ytot,a=c(1,1,1),b =c(1,1,1))
    # posterior_spatemp_cpp(thetacand,adjmatinf = adj_matrix,lag = N_time,X=xobs,y=y,a=c(1,1,1),b =c(1,1,1))
    #posterior_spatempver_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))

    auxnum = posterior_spatempcenscar_zero_cpp(theta = thetacand,adjmatinf = adjmatinf, Wt = Wt, X = xtot,y = ycom, a = aprior, b = bprior, alim = alim, blim = blim)
    auxden = posterior_spatempcenscar_zero_cpp(theta = thetalast,adjmatinf = adjmatinf, Wt = Wt, X = xtot,y = ycom, a = aprior, b= bprior, alim = alim, blim = blim)

    #num = auxnum$posterior
    #den = auxden$posterior

    # lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior)) #+ logpwcand
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior)) #+ logpwlast



    #  if(is.na(num)){
    #   num = runif(1)
    #}

    #if(is.na(den)){
    # den = runif(1)
    #}

    lognum = auxnum$posterior_log #+ logpwcand
    logden = auxden$posterior_log #+ logpwlast

    lognum = lognum + log(dbeta_rep_ab(thetalast[1],muprop[1],sigmaprop[1],alim[1],blim[1]))+
      log(dbeta_rep_ab(thetalast[2],muprop[2],sigmaprop[2],alim[2],blim[2])) +
      log(dbeta_rep_ab(thetalast[3],muprop[3],sigmaprop[3],alim[3],blim[3]))


    logden = logden + log(dbeta_rep_ab(thetacand[1],muprop[1],sigmaprop[1], alim[1],blim[1])) +
             log(dbeta_rep_ab(thetacand[2],muprop[2],sigmaprop[2], alim[2],blim[2])) +
             log(dbeta_rep_ab(thetacand[3],muprop[3],sigmaprop[3], alim[3],blim[3]))


    logunif= log(runif(1))
    n = length(y)
    p = ncol(xtot)

    if(logunif<lognum - logden)
    {
      count = count +1
      rho_sF[i] <- thetacand[1]
      rho_tF[i] = thetacand[2]
      rho_stF[i] = thetacand[3]
      corrmat = auxnum$varcovspatemp
      betaest = auxnum$betaest
      varbeta = auxnum$varbeta
      s2est = auxnum$S2

    }else
    {
      rho_sF[i] =  thetalast[1]
      rho_tF[i] =  thetalast[2]
      rho_stF[i] = thetalast[3]

      corrmat = auxden$varcovspatemp
      betaest = auxden$betaest
      varbeta = auxden$varbeta
      s2est = auxden$S2
    }

    sigma2F[i] = rinvgamma(1,0.5*(n-p), (n-p)*s2est/2)
    betaF[i,] = rmvnorm(1,betaest,sigma2F[i]*varbeta)
    covmat = sigma2F[i] *corrmat
    ############################################################
    covobsobs = covmat[indicens ==0, indicens ==0]
    covcenscens = covmat[indicens ==1, indicens ==1]
    covcensobs = covmat[indicens ==1, indicens ==0]
    covobscens = covmat[indicens ==0, indicens ==1]
    ############################################################
    yobs = y[indicens ==0]
    xobserv = xtot[indicens ==0,]
    xcens = xtot[indicens ==1,]
    covobsobsinv = solve(covobsobs)
    mucens = xcens%*%betaF[i,] + covcensobs%*%covobsobsinv%*%(yobs - xobserv%*%betaF[i,])
    Sigmacens = covcenscens - covcensobs%*%covobsobsinv%*%covobscens
    ycensF[i,] = rtmvnorm(1,mean =as.vector(mucens), sigma = Sigmacens,lower = lower,upper = upper,algorithm = "gibbs")
    ycom[indicens ==1] = ycensF[i,]
    print(c(i,sigma2F[i],rho_sF[i],rho_tF[i],rho_stF[i], count/i))
    #gc()
  }

  gc()

  rho_s_burn=rho_sF[(burn+1):iter]
  rho_s_val= rho_s_burn[seq((burn+1),iter-burn,thin)]

  rho_t_burn=rho_tF[(burn+1):iter]
  rho_t_val= rho_t_burn[seq((burn+1),iter-burn,thin)]

  rho_st_burn=rho_stF[(burn+1):iter]
  rho_st_val= rho_st_burn[seq((burn+1),iter-burn,thin)]


  sigma2Fburn=sigma2F[(burn+1):iter]
  sigma2Fval= sigma2Fburn[seq((burn+1),iter-burn,thin)]

  betaFburn=betaF[(burn+1):iter,]
  betaFval= betaFburn[seq((burn+1),iter-burn,thin),]

  ycensFburn = ycensF[(burn+1):iter,]
  ycensFval = ycensFburn[seq((burn+1),iter-burn,thin),]

  probacc = sum(count)/iter

  nmc = length(rho_s_val)
  thetaFval = cbind(betaFval,sigma2Fval,rho_s_val,rho_t_val,rho_st_val,ycensFval)


  rho_s_est = quantile(rho_s_val,probs = c(0.025,0.5,0.975))
  rho_t_est = quantile(rho_t_val,probs = c(0.025,0.5,0.975))
  rho_st_est = quantile(rho_st_val,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigma2Fval,probs = c(0.025,0.5,0.975))


  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)


  #### maybe later
  # log_lik = apply(thetaFval,1,loglikespatempcens_fix,Xmat = xtot,ytotobs = y,lags = lags,adjmatinf = adjmatinf,indicens = indicens)
  ######


  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  result = list(betaF = betaFval, sigma2F = sigma2Fval, rho_sF = rho_s_val, rho_tF = rho_t_val,rho_stF = rho_st_val,
                ycensF = ycensFval, thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rho_s_est = rho_s_est,
                rho_t_est = rho_t_est, rho_st_est = rho_st_est, probacc = probacc)

  return(result)
}

















#' @export
bayesspatempcensauto_fix_zero = function(y,xtot,thetaini,indicens, iter,burn,thin, lags,
                                    adjmatinf, aprior = c(4,4), bprior = c(2,2),
                                    divproprho = 9, divpropphi = 9,lower,upper){

  # y= yobs; xobs = xobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;
  #adjmatinf = adjm1;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time_obs




  yobs = y[indicens == 0]
  xobs = xtot[indicens ==0, ]
  rhoF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  psiF = rep(0,iter)
  phiF = rep(0,iter)
  p = ncol(xtot)
  betaF = matrix(0,iter,p)
  sigma2F = 0
  n = length(y)
  ycensF = matrix(0,iter,sum(indicens==1))
  rhoF[1] = thetaini[1]
  phiF[1] = thetaini[2]
 # corrmatini = spatimecovar_2_cpp(lag = lags,adjmatinf =adjmatinf, rho = rhoF[1], psi = 0,
  #                                phi = phiF[1])
  #invcorrmatini = solve(corrmatini)
  #betaini = solve(t(xtot)%*%invcorrmatini%*%xtot)%*%t(xtot)%*%invcorrmatini%*%y

  #sigma2ini = (t(y-xtot%*%betaini)%*%corrmatini%*%(y-xtot%*%betaini))/(n-p)
 # betaF[1,] = betaini
  #sigma2F[1] = as.numeric(sigma2ini)

  #covmatini = sigma2F[1]*corrmatini

  ############################################################
  #covobsobsini = covmatini[indicens ==0, indicens ==0]
  #covcenscensini = covmatini[indicens ==1, indicens ==1]
  #covcensobsini = covmatini[indicens ==1, indicens ==0]
  #covobscensini = covmatini[indicens ==0, indicens ==1]
  ############################################################

  #yobs = y[indicens ==0]
  #xobserv = xtot[indicens ==0,]
  #xcens = xtot[indicens ==1,]
  #covobsobsiniinv = solve(covobsobsini)
#  mucensini = xcens%*%betaini + covcensobsini%*%covobsobsiniinv%*%(yobs - xobserv%*%betaini)
 # Sigmacensini = covcenscensini - covcensobsini%*%covobsobsiniinv%*%covobscensini
  ycensF[1,] = y[indicens ==1]
  #rtmvnorm(1,mean =as.vector(mucensini), sigma = Sigmacensini,lower = lower,upper = upper,algorithm = "gibbs",thinning = 2)
  ycom = y
  ycom[indicens ==1] = ycensF[1,]
  ycom2 = matrix(0,iter,)
  adjmatinf[upper.tri(adjmatinf)]=0


  for(i in 2:iter){
    #i = 2
    thetalast = c(rhoF[i-1],phiF[i-1])

    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))

    #lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior))
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior))


    #lognum = lognum + log(dbeta(thetalast[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetalast[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetalast[3],apriorprop[3],bpriorprop[3]))

    #logden = logden + log(dbeta(thetacand[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetacand[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetacand[3],apriorprop[3],bpriorprop[3]))

    cons = 1
    muprop1 = min(0.99,thetalast[1]/cons)
    muprop2 = min(0.99,thetalast[2]/cons)

    muprop = c(muprop1,muprop2)
    sigmaprop1 = muprop[1]*(1-muprop[1])/divproprho
    sigmaprop2 = muprop[2]*(1-muprop[2])/divproppsi

    #    sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #   sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop = c(sigmaprop1,sigmaprop2)
    omega = runif(1,0.99,1)
    theta1cand = min(0.99,rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    theta2cand = min(0.99,rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))


    #if(theta1cand > 0.95){
    # theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    #theta1cand = runif(1,0.1,0.6)
    #}

    #if(theta2cand > 0.95){
    # theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    #  theta2cand = runif(1,0.1,0.45)
    #}


    #if(theta3cand> 0.95){
    # theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
    #theta3cand = runif(1,0.1,0.6)
    #}



    if(theta1cand<0.05){
      theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
      theta1cand = runif(1,0.1,0.2)
    }

    if(theta2cand<0.02){
      theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
      theta2cand = runif(1,0.1,0.2)
    }


    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    thetacand = c(theta1cand,theta2cand)

    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))
    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))
    # a = posteriorvar3_cpp(thetacand,adjmatinf = adj_matrix,lags,X,ytot,a=c(1,1,1),b =c(1,1,1))
    # posterior_spatemp_cpp(thetacand,adjmatinf = adj_matrix,lag = N_time,X=xobs,y=y,a=c(1,1,1),b =c(1,1,1))
    #posterior_spatempver_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))

    auxnum = posterior_spatempcens_zero_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xtot,y = ycom, a = aprior, b = bprior)
    auxden = posterior_spatempcens_zero_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xtot,y = ycom, a = aprior, b= bprior)

    #num = auxnum$posterior
    #den = auxden$posterior

    # lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior)) #+ logpwcand
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior)) #+ logpwlast



    #  if(is.na(num)){
    #   num = runif(1)
    #}

    #if(is.na(den)){
    # den = runif(1)
    #}

    lognum = auxnum$posterior_log #+ logpwcand
    logden = auxden$posterior_log #+ logpwlast

    lognum = lognum + log(dbeta_rep_cpp(thetalast[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetalast[2],omega*muprop[2],sigmaprop[2]))

    logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_cpp(thetacand[2],omega*muprop[2],sigmaprop[2]))



    logunif= log(runif(1))
    n = length(y)
    p = ncol(xtot)

    if(logunif<lognum - logden)
    {
      count = count +1
      rhoF[i] <- thetacand[1]
      phiF[i] = thetacand[2]
      corrmat = auxnum$varcovspatemp
      betaest = auxnum$betaest
      varbeta = auxnum$varbeta
      s2est = auxnum$S2

    }else
    {
      rhoF[i] <- thetalast[1]
      phiF[i] = thetalast[2]
      corrmat = auxden$varcovspatemp
      betaest = auxden$betaest
      varbeta = auxden$varbeta
      s2est = auxden$S2
    }

    sigma2F[i] = rinvgamma(1,0.5*(n-p), (n-p)*s2est/2)
    betaF[i,] = rmvnorm(1,betaest,sigma2F[i]*varbeta)
    covmat = sigma2F[i] *corrmat
    ############################################################
    covobsobs = covmat[indicens ==0, indicens ==0]
    covcenscens = covmat[indicens ==1, indicens ==1]
    covcensobs = covmat[indicens ==1, indicens ==0]
    covobscens = covmat[indicens ==0, indicens ==1]
    ############################################################
    yobs = y[indicens ==0]
    xobserv = xtot[indicens ==0,]
    xcens = xtot[indicens ==1,]
    covobsobsinv = solve(covobsobs)
    mucens = xcens%*%betaF[i,] + covcensobs%*%covobsobsinv%*%(yobs - xobserv%*%betaF[i,])
    Sigmacens = covcenscens - covcensobs%*%covobsobsinv%*%covobscens
    ycensF[i,] = rtmvnorm(1,mean =as.vector(mucens), sigma = Sigmacens,lower = lower,upper = upper,algorithm = "gibbs",thinning = 2)
    ycom[indicens ==1] = ycensF[i,]
    print(c(i,sigma2F[i],rhoF[i],phiF[i],count/i))
    #gc()
  }

  gc()

  rhoburn=rhoF[(burn+1):iter]
  rhoval= rhoburn[seq((burn+1),iter-burn,thin)]

  phiburn=phiF[(burn+1):iter]
  phival= phiburn[seq((burn+1),iter-burn,thin)]

  sigma2Fburn=sigma2F[(burn+1):iter]
  sigma2Fval= sigma2Fburn[seq((burn+1),iter-burn,thin)]

  betaFburn=betaF[(burn+1):iter,]
  betaFval= betaFburn[seq((burn+1),iter-burn,thin),]

  ycensFburn = ycensF[(burn+1):iter,]
  ycensFval = ycensFburn[seq((burn+1),iter-burn,thin),]

  probacc = sum(count)/iter

  nmc = length(rhoval)
  thetaFval = cbind(betaFval,sigma2Fval,rhoval,phival,ycensFval)


  rhoest = quantile(rhoval,probs = c(0.025,0.5,0.975))
  phiest = quantile(phival,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigma2Fval,probs = c(0.025,0.5,0.975))


  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)


  #### maybe later
  ##log_lik = apply(thetaFval,1,loglikespatemp_fix,Xmat = xobs,y,lags = lags,adjmatinf = adjmatinf)
  ######


  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  result = list(betaF = betaFval, sigma2F = sigma2Fval, rhoF = rhoval, phiF = phival,
                ycensF = ycensFval, thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
                phiest = phiest, probacc = probacc)#,log_lik = log_lik)

  return(result)
}




#' @export
bayesspatempcensauto_ar2_fix_zero = function(y,xtot,thetaini,indicens, iter,burn,thin, lags,
                                         adjmatinf, aprior = c(4,4,4), bprior = c(2,2,2),
                                         divproprho = 9, divpropphi = 9,lower,upper){

  # y= yobs; xobs = xobs;thetaini = thetaini; iter = iter;
  #burn = burn; thin = thin;
  #adjmatinf = adjm1;aprior = c(4,2,4); bprior = c(2,4,2);
  #divproprho = 9; divproppsi = 9; divpropphi = 9; lags = N_time_obs

  yobs = y[indicens == 0]
  xobs = xtot[indicens ==0, ]
  rhoF = rep(0,iter)
  numden = 0
  count = 0
  count1 = 0
  psiF = rep(0,iter)
  phi1F = rep(0,iter)
  phi2F = rep(0,iter)
  p = ncol(xtot)
  betaF = matrix(0,iter,p)
  sigma2F = 0
  n = length(y)
  ycensF = matrix(0,iter,sum(indicens==1))
  rhoF[1] = thetaini[1]
  phi1F[1] = thetaini[2]
  phi2F[1] = thetaini[3]
  # corrmatini = spatimecovar_2_cpp(lag = lags,adjmatinf =adjmatinf, rho = rhoF[1], psi = 0,
  #                                phi = phiF[1])
  #invcorrmatini = solve(corrmatini)
  #betaini = solve(t(xtot)%*%invcorrmatini%*%xtot)%*%t(xtot)%*%invcorrmatini%*%y

  #sigma2ini = (t(y-xtot%*%betaini)%*%corrmatini%*%(y-xtot%*%betaini))/(n-p)
  # betaF[1,] = betaini
  #sigma2F[1] = as.numeric(sigma2ini)

  #covmatini = sigma2F[1]*corrmatini

  ############################################################
  #covobsobsini = covmatini[indicens ==0, indicens ==0]
  #covcenscensini = covmatini[indicens ==1, indicens ==1]
  #covcensobsini = covmatini[indicens ==1, indicens ==0]
  #covobscensini = covmatini[indicens ==0, indicens ==1]
  ############################################################

  #yobs = y[indicens ==0]
  #xobserv = xtot[indicens ==0,]
  #xcens = xtot[indicens ==1,]
  #covobsobsiniinv = solve(covobsobsini)
  #  mucensini = xcens%*%betaini + covcensobsini%*%covobsobsiniinv%*%(yobs - xobserv%*%betaini)
  # Sigmacensini = covcenscensini - covcensobsini%*%covobsobsiniinv%*%covobscensini
  ycensF[1,] = y[indicens ==1]
  #rtmvnorm(1,mean =as.vector(mucensini), sigma = Sigmacensini,lower = lower,upper = upper,algorithm = "gibbs",thinning = 2)
  ycom = y
  ycom[indicens ==1] = ycensF[1,]
  ycom2 = matrix(0,iter,)
  adjmatinf[upper.tri(adjmatinf)]=0


  for(i in 2:iter){
  #  i = 2
    thetalast = c(rhoF[i-1],phi1F[i-1],phi2F[i-1])

    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))

    #lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior))
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior))


    #lognum = lognum + log(dbeta(thetalast[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetalast[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetalast[3],apriorprop[3],bpriorprop[3]))

    #logden = logden + log(dbeta(thetacand[1],apriorprop[1],bpriorprop[1])) + log(dbeta(thetacand[2],apriorprop[2],bpriorprop[2])) +
    # log(dbeta(thetacand[3],apriorprop[3],bpriorprop[3]))

    cons = 1
    muprop1 = min(0.99,thetalast[1]/cons)
    muprop2 = min(0.99,thetalast[2]/cons)
    muprop3 = min(0.99,thetalast[3]/cons)

    muprop = c(muprop1,muprop2,muprop3)
    sigmaprop1 = muprop[1]*(1-muprop[1])/divproprho
    sigmaprop2 = muprop[2]*(1-muprop[2])/divproppsi
    sigmaprop2 = abs(sigmaprop2)
    sigmaprop3 = muprop[3]*(1-muprop[3])/divproppsi
    sigmaprop3 = abs(sigmaprop3)

    #    sigmaprop1 = muprop[1]*(1-muprop[1])/50
    #   sigmaprop2 = muprop[2]*(1-muprop[2])/50

    #sigmaprop1 = 1/120
    #sigmaprop2 = 1/120
    sigmaprop = c(sigmaprop1,sigmaprop2,sigmaprop3)
    omega = runif(1,0.99,1)
    theta1cand = min(0.99,rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    theta3cand = min(0.99,rbeta_rep_ab(1,omega*muprop[3],sigmaprop[3],-1,1))
    theta2cand = min(0.99,rbeta_rep_ab(1,omega*muprop[2],sigmaprop[2],theta3cand-1,1-theta3cand))

    #if(theta1cand > 0.95){
    # theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    #theta1cand = runif(1,0.1,0.6)
    #}

    #if(theta2cand > 0.95){
    # theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
    #  theta2cand = runif(1,0.1,0.45)
    #}


    #if(theta3cand> 0.95){
    # theta3cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[3],sigmaprop[3]))
    #theta3cand = runif(1,0.1,0.6)
    #}



  #  if(theta1cand<0.05){
   #   theta1cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[1],sigmaprop[1]))
    #  theta1cand = runif(1,0.1,0.2)
    #}

    #if(theta2cand<0.02){
     # theta2cand =max(runif(1,0.05,0.15),rbeta_rep_cpp(omega*muprop[2],sigmaprop[2]))
      #theta2cand = runif(1,0.1,0.2)
    #}


    #apriorprop = aprior + 0.5
    #bpriorprop = bprior + 0.5
    thetacand = c(theta1cand,theta2cand, theta3cand)

    #thetacand = c(rbeta(1,apriorprop[1],bpriorprop[1]),rbeta(1,apriorprop[2],bpriorprop[2]),rbeta(1,apriorprop[3],bpriorprop[3]))
    #thetacand = c(rbeta(1,5,2),rbeta(1,2,5))

    #thetacand = c(runif(1,min = 0.2,max =  max(1,1)),
    #             runif(1,min = 0.05,max = 0.4))
    # a = posteriorvar3_cpp(thetacand,adjmatinf = adj_matrix,lags,X,ytot,a=c(1,1,1),b =c(1,1,1))
    # posterior_spatemp_cpp(thetacand,adjmatinf = adj_matrix,lag = N_time,X=xobs,y=y,a=c(1,1,1),b =c(1,1,1))
    #posterior_spatempver_cpp(theta = thetacand,adjmatinf = adjsamcand, lag = lags, X = xobs,y = y, a = aprior, b = bprior)
    #num = posteriorvar2_cpp(theta = thetacand,adjmatinf = adjsamcand,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))
    #den = posteriorvar2_cpp(theta = thetalast,adjmatinf = adjsamlast,X = xobs,y = y,mutheta = c(0.5,0.5),sigmatheta = c(1/12,1/12))

    auxnum = posterior_spatempcens_ar2_zero_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xtot,y = ycom, a = aprior, b = bprior)
    auxden = posterior_spatempcens_ar2_zero_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xtot,y = ycom, a = aprior, b= bprior)

    #num = auxnum$posterior
    #den = auxden$posterior

    # lognum = log(posterior_spatemp_cpp(theta = thetacand,adjmatinf = adjmatinf, lag = lags, X = xobs,y = y, a = aprior, b = bprior)) #+ logpwcand
    #logden = log(posterior_spatemp_cpp(theta = thetalast,adjmatinf = adjmatinf,lag = lags, X = xobs,y = y, a = aprior, b= bprior)) #+ logpwlast



    #  if(is.na(num)){
    #   num = runif(1)
    #}

    #if(is.na(den)){
    # den = runif(1)
    #}

    lognum = auxnum$posterior_log #+ logpwcand
    logden = auxden$posterior_log #+ logpwlast

    lognum = lognum + log(dbeta_rep_cpp(thetalast[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_ab(thetalast[2],omega*muprop[2],sigmaprop[2],thetalast[3]-1,1-thetalast[3])) +
      log(dbeta_rep_ab(thetalast[3],omega*muprop[3],sigmaprop[3],-1,1))

    logden = logden + log(dbeta_rep_cpp(thetacand[1],omega*muprop[1],sigmaprop[1])) + log(dbeta_rep_ab(thetacand[2],omega*muprop[2],sigmaprop[2], thetacand[3]-1, 1-thetacand[3])) +
      log(dbeta_rep_ab(thetacand[3],omega*muprop[3],sigmaprop[3], -1, 1))



    logunif= log(runif(1))
    n = length(y)
    p = ncol(xtot)

    if(logunif<lognum - logden)
    {
      count = count +1
      rhoF[i] <- thetacand[1]
      phi1F[i] = thetacand[2]
      phi2F[i] = thetacand[3]
      corrmat = auxnum$varcovspatemp
      betaest = auxnum$betaest
      varbeta = auxnum$varbeta
      s2est = auxnum$S2

    }else
    {
      rhoF[i] <- thetalast[1]
      phi1F[i] = thetalast[2]
      phi2F[i] = thetalast[3]
      corrmat = auxden$varcovspatemp
      betaest = auxden$betaest
      varbeta = auxden$varbeta
      s2est = auxden$S2
    }

    sigma2F[i] = rinvgamma(1,0.5*(n-p), (n-p)*s2est/2)
    betaF[i,] = rmvnorm(1,betaest,sigma2F[i]*varbeta)
    covmat = sigma2F[i] *corrmat
    ############################################################
    covobsobs = covmat[indicens ==0, indicens ==0]
    covcenscens = covmat[indicens ==1, indicens ==1]
    covcensobs = covmat[indicens ==1, indicens ==0]
    covobscens = covmat[indicens ==0, indicens ==1]
    ############################################################
    yobs = y[indicens ==0]
    xobserv = xtot[indicens ==0,]
    xcens = xtot[indicens ==1,]
    covobsobsinv = solve(covobsobs)
    mucens = xcens%*%betaF[i,] + covcensobs%*%covobsobsinv%*%(yobs - xobserv%*%betaF[i,])
    Sigmacens = covcenscens - covcensobs%*%covobsobsinv%*%covobscens
    ycensF[i,] = rtmvnorm(1,mean =as.vector(mucens), sigma = Sigmacens,lower = lower,upper = upper,algorithm = "gibbs",thinning = 2)
    ycom[indicens ==1] = ycensF[i,]
    print(c(i,sigma2F[i],rhoF[i],phi1F[i], phi2F[i], count/i))
    #gc()
  }

  gc()

  rhoburn=rhoF[(burn+1):iter]
  rhoval= rhoburn[seq((burn+1),iter-burn,thin)]

  phi1burn=phi1F[(burn+1):iter]
  phi1val= phi1burn[seq((burn+1),iter-burn,thin)]

  phi2burn=phi2F[(burn+1):iter]
  phi2val= phi2burn[seq((burn+1),iter-burn,thin)]


  sigma2Fburn=sigma2F[(burn+1):iter]
  sigma2Fval= sigma2Fburn[seq((burn+1),iter-burn,thin)]

  betaFburn=betaF[(burn+1):iter,]
  betaFval= betaFburn[seq((burn+1),iter-burn,thin),]

  ycensFburn = ycensF[(burn+1):iter,]
  ycensFval = ycensFburn[seq((burn+1),iter-burn,thin),]

  probacc = sum(count)/iter

  nmc = length(rhoval)
  thetaFval = cbind(betaFval,sigma2Fval,rhoval,phi1val, phi2val, ycensFval)


  rhoest = quantile(rhoval,probs = c(0.025,0.5,0.975))
  phi1est = quantile(phi1val,probs = c(0.025,0.5,0.975))
  phi2est = quantile(phi2val,probs = c(0.025,0.5,0.975))
  betaest = apply(betaFval,2,quantile,probs = c(0.025,0.5,0.975))
  sigmaest = quantile(sigma2Fval,probs = c(0.025,0.5,0.975))


  #start_time <- Sys.time()
  #log_lik = rbind()
  #for(i in 1:length(rhoval)){
  #  log_lik = rbind(log_lik,loglikespatemp(thetaFval[[i]],xobs,y,lags = lags))
  #}
  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = sapply(thetaFval,loglikespatemp,Xmat = xobs,y,lags = lags)


  #### maybe later
  ##log_lik = apply(thetaFval,1,loglikespatemp_fix,Xmat = xobs,y,lags = lags,adjmatinf = adjmatinf)
  ######


  #end_time <- Sys.time()
  #elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #start_time <- Sys.time()
  #log_lik = t(sapply(thetaFval,loglikespatemp_cpp,Xmat = xobs,y,lags = lags))
  #end_time <- Sys.time()
  ##elapsed_time1 <- end_time - start_time
  #print(elapsed_time1)

  #thetaFval = cbind(betaFval,sigmaFval,rhoval,psival)

  result = list(betaF = betaFval, sigma2F = sigma2Fval, rhoF = rhoval, phi1F = phi1val, phi2F = phi2val,
                ycensF = ycensFval, thetaFval = thetaFval, betaest = betaest, sigmaest = sigmaest, rhoest = rhoest,
                phi1est = phi1est, phi2est = phi2est, probacc = probacc)#,log_lik = log_lik)

  return(result)
}


matrixb2 = function(adjmatinf,rho){
  ##adjmatinf = adjm1
  nneigh = rowSums(adjmatinf)
  b =  rho / (1 + nneigh * rho^2)
  for (i in 1:length(nneigh)){
    adjmatinf[i,][adjmatinf[i,]==1] = b[i]
  }
  return(adjmatinf)
}





matrixF <- function(nneigh, rho) {
  N <- length(nneigh)
  tau <- (1 + (nneigh - 1) * rho^2) / (1 - rho^2)
  Fmat <- diag(tau)
  return(Fmat)
}

matrixFinv <- function(nneigh, rho) {
  N <- length(nneigh)
  tau <- (1 + (nneigh - 1) * rho^2) / (1 - rho^2)
  Fmatinv <- diag(1/tau)
  return(Fmatinv)
}


varcovdagar_2 = function (adjmatinf,rho,sigma2,tau2){
  nneigh = rowSums(adjmatinf)
  matB = matrixb2(adjmatinf,rho = rho)
  matFinv = sigma2*matrixFinv(nneigh,rho)
  I = diag(1,nrow(matFinv),nrow(matFinv))
  inv_I_minus_B = solve(I-matB)
  return (inv_I_minus_B%*%matFinv%*%t(inv_I_minus_B) + tau2*I)
}

#' @export
pred_dist_prior = function(listres,interval=FALSE,predind=NULL,X,yobs){
##listres = res
##  rhoval = listres$rhoF
##  psival = listres$psiF
##  betaval = listres$betaF
##  sigmaval = listres$sigmaF
##  WFval = listres$thetaFval[[2]]
  h1 = 0
  #X = xobs
  #yobs = y
  p = ncol(X)
  #listres = res
  thetaaux = listres$thetaFval
  for ( j in 1:length(thetaaux)){
    # i = 1
    #j = 1
   # listres = res
    thetaaux2 = thetaaux[[j]]$thetaFval
    rhoaux = thetaaux2[p+2]
    psiaux = thetaaux2[p+3]
    adjmatinfaux = listres$thetaFval[[j]]$adjmatinf
    h1[j] = auxiliarcpp::posteriorvar3_cpp(theta = c(rhoaux,psiaux),adjmatinf = adjmatinfaux,X = xobs,y = y, a = c(4,2), b = c(2,4))
    #h2[i] = posteriorvar3_cpp(theta = c(median(rhoval),median(psival)),adjmatinf = WFval[[i]],X = xobs,y = y, a = c(4,2), b = c(2,4))
      }

  posmax= which(h1 == max(h1))
  posmax=posmax[1]

  thetaFmax = thetaaux[[posmax]]

  psi = thetaFmax$thetaFval[p+3]
  rho = thetaFmax$thetaFval[p+2]
  sigma2 = thetaFmax$thetaFval[p+1]
  tau2 = psi*sigma2
  beta = thetaFmax$thetaFval[1:p]
  adjmatinf = thetaFmax$adjmatinf
  adjmatinf[upper.tri(adjmatinf)] = 0

  varcov = varcovdagar_2(adjmatinf,rho,sigma2,tau2)
  if(is.null(predind)){
    muobs = X%*%beta
    if(interval == FALSE){
      return(muobs)
    }else{
      q = qnorm(0.975)
      li = muobs - q*diag(varcov)
      ls = muobs + q*diag(varcov)
      return(cbind(muobs,li,ls))
    }
  }else{
    obs = (predind ==0)
    pred = (predind==1)
    xobs = X[obs,]
    muobs = xobs%*%beta
    xpred = X[pred,]
    varcovobsinv = solve(varcov[obs,obs])
    mupred = xpred%*%beta + varcov[pred,obs]%*%varcovobsinv%*%(yobs-muobs)
    if(interval == FALSE){
      return(mupred)
    }else{
      sigmapred = varcov[pred,pred] - varcov[pred,obs]%*%varcovobsinv%*%varcov[obs,pred]
      q = qnorm(0.975)
      li  = mupred - q*diag(sigmapred)
      ls  = mupred + q*diag(sigmapred)
      return(cbind(mupred,li,ls))
    }
  }
}


#' @export
pred_dist_prior_ind = function(listparam,interval=FALSE,predind=NULL,Xmat,yobs){
  ##listres = res
  ##  rhoval = listres$rhoF
  ##  psival = listres$psiF
  ##  betaval = listres$betaF
  ##  sigmaval = listres$sigmaF
  ##  WFval = listres$thetaFval[[2]]
  #listparam = h[[1]]$res$thetaFval[[1]]
  #X = xobs
  #yobs = y
  p = ncol(Xmat)
  #listres = res
  thetachain = listparam$thetaFval

  psi = thetachain[p+3]
  rho = thetachain[p+2]
  sigma2 = thetachain[p+1]
  tau2 = psi*sigma2
  beta = thetachain[1:p]
  adjmatinf = listparam$adjmatinf
  adjmatinf[upper.tri(adjmatinf)] = 0

  varcov = varcovdagar_2(adjmatinf,rho,sigma2,tau2)
  if(is.null(predind)){
    muobs = Xmat%*%beta
    if(interval == FALSE){
      return(muobs)
    }else{
      q = qnorm(0.975)
      li = muobs - q*diag(varcov)
      ls = muobs + q*diag(varcov)
      return(cbind(muobs,li,ls))
    }
  }else{
    obs = (predind ==0)
    pred = (predind==1)
    xobs = Xmat[obs,]
    muobs = xobs%*%beta
    xpred = Xmat[pred,]
    varcovobsinv = solve(varcov[obs,obs])
    mupred = xpred%*%beta + varcov[pred,obs]%*%varcovobsinv%*%(yobs-muobs)
    if(interval == FALSE){
      return(mupred)
    }else{
      sigmapred = varcov[pred,pred] - varcov[pred,obs]%*%varcovobsinv%*%varcov[obs,pred]
      q = qnorm(0.975)
      li  = mupred - q*diag(sigmapred)
      ls  = mupred + q*diag(sigmapred)
      return(cbind(mupred,li,ls))
    }
  }
}

#### spatiotemporal predictions ###########################

#### computing the posterior of rho,psi, phi, for a particular sample

posterior_sam = function(listreswithin,Xmat,yobs,aprior,bprior,lags_obs){
  #listres = thetaFval[[1]]
  #Xmat = xobs
  p = ncol(Xmat)
  thetaaux = listreswithin$thetaFval
  rhoaux = thetaaux[p+2]
  psiaux = thetaaux[p+3]
  phiaux = thetaaux[p+4]
  adjmatinfaux = listreswithin$adjmatinf
  theta = c(rhoaux,psiaux,phiaux)
  h = posterior_spatemp_cpp(theta,adjmatinf =adjmatinfaux,lag = lags_obs,X = Xmat,y = yobs,a = aprior,b =bprior)
  return(h)
  }


### example
#listreswithin = thetaFval[[1]]
#Xmat = xobs
#yobs = y
#aprior = aprior
#bprior = bprior
#lags = lags

#posterior_sam(listreswithin,Xmat,yobs,aprior,bprior,lags)

#' @export
pred_dist_spatemp = function(listres,interval=FALSE,predind=NULL,Xmat,yobs,lags,lags_obs=NULL,aprior,bprior){

 #listres =res
 #interval=FALSE
 #predind=indipred
 #Xmat = xtot
 #yobs = yobs
 #lags = N_time
 #aprior = c(4,2,4)
 #bprior = c(2,4,2)
 #lags_obs = N_time_obs


  if(is.null(predind)){
    thetaaux = listres$thetaFval
    h1 = sapply(X=thetaaux, posterior_sam, Xmat = Xmat,yobs= yobs, a = aprior,b = bprior,lags_obs = lags)
    p = ncol(Xmat)
    posmax= which(h1 == max(h1))
    posmax=posmax[1]

    thetaFmax = thetaaux[[posmax]]

    phi = thetaFmax$thetaFval[p+4]
    psi = thetaFmax$thetaFval[p+3]
    rho = thetaFmax$thetaFval[p+2]
    sigma2 = thetaFmax$thetaFval[p+1]
    ##tau2 = psi*sigma2
    beta = thetaFmax$thetaFval[1:p]
    adjmatinf = thetaFmax$adjmatinf
    adjmatinf[upper.tri(adjmatinf)] = 0

    # varcov = varcovdagar_2(adjmatinf,rho,sigma2,tau2)
    varcov = sigma2*spatimecovar_2_cpp(lag = lags, adjmatinf = adjmatinf,rho,psi,phi)
    muobs = Xmat%*%beta
    if(interval == FALSE){
      return(muobs)
    }else{
      q = qnorm(0.975)
      li = muobs - q*diag(varcov)
      ls = muobs + q*diag(varcov)
      return(cbind(muobs,li,ls))
    }
  }else{
    obs = (predind ==0)
    pred = (predind==1)
    xobs = Xmat[obs,]
    xpred = Xmat[pred,]
    thetaaux = listres$thetaFval
    h1 = sapply(X=thetaaux, posterior_sam, Xmat = xobs,yobs= yobs, a = aprior,b = bprior,lags_obs = lags_obs)
    p = ncol(Xmat)
    posmax= which(h1 == max(h1))
    posmax=posmax[1]

    thetaFmax = thetaaux[[posmax]]

    phi = thetaFmax$thetaFval[p+4]
    psi = thetaFmax$thetaFval[p+3]
    rho = thetaFmax$thetaFval[p+2]
    sigma2 = thetaFmax$thetaFval[p+1]
    ##tau2 = psi*sigma2
    beta = thetaFmax$thetaFval[1:p]
    adjmatinf = thetaFmax$adjmatinf
    adjmatinf[upper.tri(adjmatinf)] = 0

    # varcov = varcovdagar_2(adjmatinf,rho,sigma2,tau2)
    muobs = xobs%*%beta
    varcov = sigma2*spatimecovar_2_cpp(lag = lags, adjmatinf = adjmatinf,rho,psi,phi)



    varcovobsinv = solve(varcov[obs,obs])
    mupred = xpred%*%beta + varcov[pred,obs]%*%varcovobsinv%*%(yobs-muobs)
    if(interval == FALSE){
      return(mupred)
    }else{
      sigmapred = varcov[pred,pred] - varcov[pred,obs]%*%varcovobsinv%*%varcov[obs,pred]
      q = qnorm(0.975)
      li  = mupred - q*diag(sigmapred)
      ls  = mupred + q*diag(sigmapred)
      return(cbind(mupred,li,ls))
    }
  }
}


### example#######
#listres = list()
#listres$thetaFval = thetaFval
#pred_dist_spatemp(listres,interval=FALSE,predind=NULL,Xmat = xobs,yobs=y,lags=lags)


#' @export
pred_dist_spatemp_ind = function(listparam,predind=NULL,Xmat,yobs,lags){
  ##listres = res
  ##  rhoval = listres$rhoF
  ##  psival = listres$psiF
  ##  betaval = listres$betaF
  ##  sigmaval = listres$sigmaF
  ##  WFval = listres$thetaFval[[2]]
  #listparam = h[[1]]$res$thetaFval[[1]]
  #X = xobs
  #yobs = y
  p = ncol(Xmat)
  #listres = res
  thetachain = listparam$thetaFval
  phi = thetachain[p+4]
  psi = thetachain[p+3]
  rho = thetachain[p+2]
  sigma2 = thetachain[p+1]
  ##tau2 = psi*sigma2
  beta = thetachain[1:p]
  adjmatinf = listparam$adjmatinf
  adjmatinf[upper.tri(adjmatinf)] = 0

  ##varcov = varcovdagar_2(adjmatinf,rho,sigma2,tau2)
  varcov = sigma2*spatimecovar_2_cpp(lag = lags, adjmatinf = adjmatinf,rho,psi,phi)
  if(is.null(predind)){
    sample = t(rmvnorm(1,Xmat%*%beta,diag(1,nrow(Xmat))))
  }else{
    obs = (predind ==0)
    pred = (predind==1)
    xobs = Xmat[obs,]
    muobs = xobs%*%beta
    xpred = Xmat[pred,]
    varcovobsinv = solve(varcov[obs,obs])
    mupred = xpred%*%beta + varcov[pred,obs]%*%varcovobsinv%*%(yobs-muobs)
    sigmapred = varcov[pred,pred] - varcov[pred,obs]%*%varcovobsinv%*%varcov[obs,pred]
    sample = t(rmvnorm(1,mupred,sigmapred))
  }
return(sample)
  }


#' @export
pred_dist_spatempcens_ind = function(theta,predind=NULL,Xmat,ytotobs,lags,indicens,adjmatinf){
  ##listres = res
  ##  rhoval = listres$rhoF
  ##  psival = listres$psiF
  ##  betaval = listres$betaF
  ##  sigmaval = listres$sigmaF
  ##  WFval = listres$thetaFval[[2]]
  #listparam = h[[1]]$res$thetaFval[[1]]
  #X = xobs
  #yobs = y
  #Xmat =xtotobs
  #theta = res$thetaFval[1,]

  p = ncol(Xmat)
  #listres = res
  #thetachain = listparam$thetaFval
  phi = theta[p+4]
  psi = theta[p+3]
  rho = theta[p+2]
  sigma2 = theta[p+1]
  ##tau2 = psi*sigma2
  beta = theta[1:p]
  ycensF = theta[(p+5):length(theta)]
  ytotobs[indicens ==1] = ycensF
 # adjmatinf = listparam$adjmatinf
  adjmatinf[upper.tri(adjmatinf)] = 0

  ##varcov = varcovdagar_2(adjmatinf,rho,sigma2,tau2)
  varcov = sigma2*spatimecovar_2_cpp(lag = lags, adjmatinf = adjmatinf,rho,psi,phi)
  if(is.null(predind)){
    sample = t(rmvnorm(1,Xmat%*%beta,diag(1,nrow(Xmat))))
  }else{
    obs = (predind ==0)
    pred = (predind==1)
    xobs = Xmat[obs,]
    muobs = xobs%*%beta
    xpred = Xmat[pred,]
    varcovobsinv = solve(varcov[obs,obs])
    mupred = xpred%*%beta + varcov[pred,obs]%*%varcovobsinv%*%(ytotobs-muobs)
    sigmapred = varcov[pred,pred] - varcov[pred,obs]%*%varcovobsinv%*%varcov[obs,pred]
    sample = t(rmvnorm(1,mupred,sigmapred))
  }
  return(sample)
}



#' @export
pred_dist_spatempcens_ar2_ind = function(theta,predind=NULL,Xmat,ytotobs,lags,indicens,adjmatinf){
  ##listres = res
  ##  rhoval = listres$rhoF
  ##  psival = listres$psiF
  ##  betaval = listres$betaF
  ##  sigmaval = listres$sigmaF
  ##  WFval = listres$thetaFval[[2]]
  #listparam = h[[1]]$res$thetaFval[[1]]
  #X = xobs
  #yobs = y
  #Xmat =xtotobs
  #theta = res$thetaFval[1,]

  p = ncol(Xmat)
  #listres = res
  #thetachain = listparam$thetaFval
  phi2 = theta[p+5]
  phi1 = theta[p+4]
  psi = theta[p+3]
  rho = theta[p+2]
  sigma2 = theta[p+1]
  ##tau2 = psi*sigma2
  beta = theta[1:p]
  ycensF = theta[(p+6):length(theta)]
  ytotobs[indicens ==1] = ycensF
  # adjmatinf = listparam$adjmatinf
  adjmatinf[upper.tri(adjmatinf)] = 0

  ##varcov = varcovdagar_2(adjmatinf,rho,sigma2,tau2)
  varcov = sigma2*spatimecovar_ar2_cpp(lag = lags, adjmatinf = adjmatinf,rho,psi,phi1,phi2)
  if(is.null(predind)){
    sample = t(rmvnorm(1,Xmat%*%beta,diag(1,nrow(Xmat))))
  }else{
    obs = (predind ==0)
    pred = (predind==1)
    xobs = Xmat[obs,]
    muobs = xobs%*%beta
    xpred = Xmat[pred,]
    varcovobsinv = solve(varcov[obs,obs])
    mupred = xpred%*%beta + varcov[pred,obs]%*%varcovobsinv%*%(ytotobs-muobs)
    sigmapred = varcov[pred,pred] - varcov[pred,obs]%*%varcovobsinv%*%varcov[obs,pred]
    sample = t(rmvnorm(1,mupred,sigmapred))
  }
  return(sample)
}



#' @export
pred_dist_spatempcens_car_ind = function(theta,predind=NULL,Xmat,ytotobs,lags,indicens,adj_matcom){
  ##listres = res
  ##  rhoval = listres$rhoF
  ##  psival = listres$psiF
  ##  betaval = listres$betaF
  ##  sigmaval = listres$sigmaF
  ##  WFval = listres$thetaFval[[2]]
  #listparam = h[[1]]$res$thetaFval[[1]]
  #X = xobs
  #yobs = y
  #Xmat =xtotobs
  #theta = res$thetaFval[1,]

  p = ncol(Xmat)
  #listres = res
  #thetachain = listparam$thetaFval
  rho_st = theta[p+4]
  rho_t = theta[p+3]
  rho_s = theta[p+2]
  sigma2 = theta[p+1]
  ##tau2 = psi*sigma2
  beta = theta[1:p]
  ycensF = theta[(p+5):length(theta)]
  ytotobs[indicens ==1] = ycensF
  # adjmatinf = listparam$adjmatinf
  Wt <- abs(outer(1:lags, 1:lags, "-")) == 1
  Wt = 1*Wt

  varcov =sigma2 * spatimecovarcar_zero_cpp(adjmatinf = adj_matcom,Wt =Wt, rho_s = rho_s,rho_t = rho_t, rho_st = rho_st)


  ##varcov = varcovdagar_2(adjmatinf,rho,sigma2,tau2)
  if(is.null(predind)){
    sample = t(rmvnorm(1,Xmat%*%beta,diag(1,nrow(Xmat))))
  }else{
    obs = (predind ==0)
    pred = (predind==1)
    xobs = Xmat[obs,]
    muobs = xobs%*%beta
    xpred = Xmat[pred,]
    varcovobsinv = solve(varcov[obs,obs])
    mupred = xpred%*%beta + varcov[pred,obs]%*%varcovobsinv%*%(ytotobs-muobs)
    sigmapred = varcov[pred,pred] - varcov[pred,obs]%*%varcovobsinv%*%varcov[obs,pred]
    sample = t(rmvnorm(1,mupred,sigmapred))
  }
  return(sample)
}












#' @export
pred_dist_spatempcens_sar_ar2_ind = function(theta,predind=NULL,Xmat,ytotobs,lags,indicens,adjmatinf){
  ##listres = res
  ##  rhoval = listres$rhoF
  ##  psival = listres$psiF
  ##  betaval = listres$betaF
  ##  sigmaval = listres$sigmaF
  ##  WFval = listres$thetaFval[[2]]
  #listparam = h[[1]]$res$thetaFval[[1]]
  #X = xobs
  #yobs = y
  #Xmat =xtotobs
  #theta = res$thetaFval[1,]

  p = ncol(Xmat)
  #listres = res
  #thetachain = listparam$thetaFval
  phi2 = theta[p+5]
  phi1 = theta[p+4]
  psi = theta[p+3]
  rho = theta[p+2]
  sigma2 = theta[p+1]
  ##tau2 = psi*sigma2
  beta = theta[1:p]
  ycensF = theta[(p+6):length(theta)]
  ytotobs[indicens ==1] = ycensF
  # adjmatinf = listparam$adjmatinf
  #adjmatinf[upper.tri(adjmatinf)] = 0

  ##varcov = varcovdagar_2(adjmatinf,rho,sigma2,tau2)
  varcov = sigma2*spatimecovarsar_2_ar2_cpp(lag = lags, adjmatinf = adjmatinf,rho,psi,phi1,phi2)
  if(is.null(predind)){
    sample = t(rmvnorm(1,Xmat%*%beta,diag(1,nrow(Xmat))))
  }else{
    obs = (predind ==0)
    pred = (predind==1)
    xobs = Xmat[obs,]
    muobs = xobs%*%beta
    xpred = Xmat[pred,]
    varcovobsinv = solve(varcov[obs,obs])
    mupred = xpred%*%beta + varcov[pred,obs]%*%varcovobsinv%*%(ytotobs-muobs)
    sigmapred = varcov[pred,pred] - varcov[pred,obs]%*%varcovobsinv%*%varcov[obs,pred]
    sample = t(rmvnorm(1,mupred,sigmapred))
  }
  return(sample)
}



#' @export
pred_dist_spatempcens_sar_ind = function(theta,predind=NULL,Xmat,ytotobs,lags,indicens,adjmatinf){
  ##listres = res
  ##  rhoval = listres$rhoF
  ##  psival = listres$psiF
  ##  betaval = listres$betaF
  ##  sigmaval = listres$sigmaF
  ##  WFval = listres$thetaFval[[2]]
  #listparam = h[[1]]$res$thetaFval[[1]]
  #X = xobs
  #yobs = y
  #Xmat =xtotobs
  #theta = res$thetaFval[1,]

  p = ncol(Xmat)
  #listres = res
  #thetachain = listparam$thetaFval
  phi = theta[p+4]
  psi = theta[p+3]
  rho = theta[p+2]
  sigma2 = theta[p+1]
  ##tau2 = psi*sigma2
  beta = theta[1:p]
  ycensF = theta[(p+5):length(theta)]
  ytotobs[indicens ==1] = ycensF
  # adjmatinf = listparam$adjmatinf
 # adjmatinf[upper.tri(adjmatinf)] = 0

  ##varcov = varcovdagar_2(adjmatinf,rho,sigma2,tau2)
  varcov = sigma2*spatimecovarsar_21_cpp(lag = lags, adjmatinf = adjmatinf,rho,psi,phi)
  if(is.null(predind)){
    sample = t(rmvnorm(1,Xmat%*%beta,diag(1,nrow(Xmat))))
  }else{
    obs = (predind ==0)
    pred = (predind==1)
    xobs = Xmat[obs,]
    muobs = xobs%*%beta
    xpred = Xmat[pred,]
    varcovobsinv = solve(varcov[obs,obs])
    mupred = xpred%*%beta + varcov[pred,obs]%*%varcovobsinv%*%(ytotobs-muobs)
    sigmapred = varcov[pred,pred] - varcov[pred,obs]%*%varcovobsinv%*%varcov[obs,pred]
    sample = t(rmvnorm(1,mupred,sigmapred))
  }
  return(sample)
}






#' @export
pred_dist_spatemp_ind_fix = function(theta,adjmatinf,predind=NULL,Xmat,yobs,lags){
  p = ncol(Xmat)
  #listres = res
  thetachain = theta
  phi = thetachain[p+4]
  psi = thetachain[p+3]
  rho = thetachain[p+2]
  sigma2 = thetachain[p+1]
  ##tau2 = psi*sigma2
  beta = thetachain[1:p]
  #adjmatinf = listparam$adjmatinf
  adjmatinf[upper.tri(adjmatinf)] = 0

  ##varcov = varcovdagar_2(adjmatinf,rho,sigma2,tau2)
  varcov = sigma2*spatimecovar_2_cpp(lag = lags, adjmatinf = adjmatinf,rho,psi,phi)
  if(is.null(predind)){
    sample = t(rmvnorm(1,Xmat%*%beta,diag(1,nrow(Xmat))))
  }else{
    obs = (predind ==0)
    pred = (predind==1)
    xobs = Xmat[obs,]
    muobs = xobs%*%beta
    xpred = Xmat[pred,]
    varcovobsinv = solve(varcov[obs,obs])
    mupred = xpred%*%beta + varcov[pred,obs]%*%varcovobsinv%*%(yobs-muobs)
    sigmapred = varcov[pred,pred] - varcov[pred,obs]%*%varcovobsinv%*%varcov[obs,pred]
    sample = t(rmvnorm(1,mupred,sigmapred))
  }
  return(sample)
}



#' @export
pred_dist_spatemp_sar_ind_fix = function(theta,adjmatinf,predind=NULL,Xmat,yobs,lags){
  p = ncol(Xmat)
  #listres = res
  thetachain = theta
  phi = thetachain[p+4]
  psi = thetachain[p+3]
  rho = thetachain[p+2]
  sigma2 = thetachain[p+1]
  ##tau2 = psi*sigma2
  beta = thetachain[1:p]
  #adjmatinf = listparam$adjmatinf
  adjmatinf[upper.tri(adjmatinf)] = 0

  ##varcov = varcovdagar_2(adjmatinf,rho,sigma2,tau2)
  varcov = sigma2*spatimecovarsar_21_cpp(lag = lags, adjmatinf = adjmatinf,rho,psi,phi)
  if(is.null(predind)){
    sample = t(rmvnorm(1,Xmat%*%beta,diag(1,nrow(Xmat))))
  }else{
    obs = (predind ==0)
    pred = (predind==1)
    xobs = Xmat[obs,]
    muobs = xobs%*%beta
    xpred = Xmat[pred,]
    varcovobsinv = solve(varcov[obs,obs])
    mupred = xpred%*%beta + varcov[pred,obs]%*%varcovobsinv%*%(yobs-muobs)
    sigmapred = varcov[pred,pred] - varcov[pred,obs]%*%varcovobsinv%*%varcov[obs,pred]
    sample = t(rmvnorm(1,mupred,sigmapred))
  }
  return(sample)
}














#' @export
log_likelihood_ar1 = function(x, y, cc, lag, adjmatinf, theta,lower, upper){
  p = ncol(x)
  beta = theta[1:p]
  sigma2 = theta[p+1]
  rho = theta[p+2]
  psi = theta[p+3]
  phi = theta[p+4]
  #ycensf = theta[p+5:(length(theta))]
  eta = x%*%beta
  yobs = y[cc == 0]
  #ycens = y[cc = 1]
  covmat =sigma2 * spatimecovar_2_cpp(lag = lag, adjmatinf = adjmatinf, rho, psi, phi)
  covobsobs = covmat[cc ==0, cc ==0]
  covcenscens = covmat[cc ==1, cc ==1]
  covcensobs = covmat[cc ==1, cc ==0]
  covobscens = covmat[cc ==0, cc ==1]
  ############################################################
  covobsobsinv = solve(covobsobs)
  mucens = eta[cc == 1] + covcensobs%*%covobsobsinv%*%(yobs - eta[cc ==0])
  Sigmacens = covcenscens - covcensobs%*%covobsobsinv%*%covobscens
  like = dmvnorm(rep(0,length(yobs)), mean = eta[cc == 0], sigma = covobsobs, log = TRUE)
  + log(pmvnorm(lower = lower, upper = upper, mean = as.numeric(mucens), sigma = Sigmacens))
  return(like)
}


#' @export
inclattice=function(m){
  n=m^2
  Minc=matrix(0,n,n)
  for(i in 1:(m-1))	for(j in 1:(m-1)) Minc[(i-1)*m+j,(i-1)*m+j+1]=Minc[(i-1)*m+j,i*m+j]=1
  for(i in 1:(m-1)) Minc[(i-1)*m+m,i*m+m]=1
  for(j in 1:(m-1)) Minc[(m-1)*m+j,(m-1)*m+j+1]=1
  Minc+t(Minc)
}
