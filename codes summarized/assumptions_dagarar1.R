library(coda)
library(tidyverse)
library(patchwork)
library(ggpubr)

setwd("C:/Users/ALEJANDRO/OneDrive/Documentos/Alejo/postdoc20242/censored_network_dagar/datasets/Beijing/codes sumarized/")

#results = readRDS("ar1dagar_beijing_212_40000.rds")

cad1 = readRDS("cad1ar1dagar_Beijing.rds" )

cad2 = readRDS("cad2ar1dagar_Beijing.rds" )

cad3 = readRDS("cad3ar1dagar_Beijing.rds" )

results = list(cad1 = cad1, cad2 = cad2, cad3 = cad3)

iter = 40000
burn = 10000
thin = 50
chains = c(1,2,3)

p = 4

######### chain 1 ########################

cad1beta = results[[chains[1]]]$thetaF[,1:p]
cad1betaburn = cad1beta[burn:iter,]
cad1betaval = cad1betaburn[seq(1,iter-burn,thin), ]

cad1sigma2 = results[[chains[1]]]$thetaF[,(p+1)]
cad1sigma2burn = cad1sigma2[burn:iter]
cad1sigma2val = cad1sigma2burn[seq(1,iter-burn,thin)]

cad1rho = results[[chains[1]]]$thetaF[,(p+2)]
cad1rhoburn = cad1rho[burn:iter]
cad1rhoval = cad1rhoburn[seq(1,iter-burn,thin)]


cad1phi = results[[chains[1]]]$thetaF[,(p+4)]
cad1phiburn = cad1phi[burn:iter]
cad1phival = cad1phiburn[seq(1,iter-burn,thin)]


cad1tau2 = results[[chains[1]]]$thetaF[,(p+5)]
cad1tau2burn = cad1tau2[burn:iter]
cad1tau2val = cad1tau2burn[seq(1,iter-burn,thin)]

cad1 = list(betaF = cad1betaval,sigma2F = cad1sigma2val, rhoF = cad1rhoval, phiF = cad1phival,
            tau2F = cad1tau2val)


######### chain 2 ########################

cad2beta = results[[chains[2]]]$thetaF[,1:p]
cad2betaburn = cad2beta[burn:iter,]
cad2betaval = cad2betaburn[seq(1,iter-burn,thin), ]

cad2sigma2 = results[[chains[2]]]$thetaF[,(p+1)]
cad2sigma2burn = cad2sigma2[burn:iter]
cad2sigma2val = cad2sigma2burn[seq(1,iter-burn,thin)]

cad2rho = results[[chains[2]]]$thetaF[,(p+2)]
cad2rhoburn = cad2rho[burn:iter]
cad2rhoval = cad2rhoburn[seq(1,iter-burn,thin)]


cad2phi = results[[chains[2]]]$thetaF[,(p+4)]
cad2phiburn = cad2phi[burn:iter]
cad2phival = cad2phiburn[seq(1,iter-burn,thin)]


cad2tau2 = results[[chains[2]]]$thetaF[,(p+5)]
cad2tau2burn = cad2tau2[burn:iter]
cad2tau2val = cad2tau2burn[seq(1,iter-burn,thin)]

cad2 = list(betaF = cad2betaval,sigma2F = cad2sigma2val, rhoF = cad2rhoval, phiF = cad2phival,
            tau2F = cad2tau2val)


######### chain 3 ########################

cad3beta = results[[chains[3]]]$thetaF[,1:p]
cad3betaburn = cad3beta[burn:iter,]
cad3betaval = cad3betaburn[seq(1,iter-burn,thin), ]

cad3sigma2 = results[[chains[3]]]$thetaF[,(p+1)]
cad3sigma2burn = cad3sigma2[burn:iter]
cad3sigma2val = cad3sigma2burn[seq(1,iter-burn,thin)]

cad3rho = results[[chains[3]]]$thetaF[,(p+2)]
cad3rhoburn = cad3rho[burn:iter]
cad3rhoval = cad3rhoburn[seq(1,iter-burn,thin)]


cad3phi = results[[chains[3]]]$thetaF[,(p+4)]
cad3phiburn = cad3phi[burn:iter]
cad3phival = cad3phiburn[seq(1,iter-burn,thin)]


cad3tau2 = results[[chains[3]]]$thetaF[,(p+5)]
cad3tau2burn = cad3tau2[burn:iter]
cad3tau2val = cad3tau2burn[seq(1,iter-burn,thin)]

cad3 = list(betaF = cad3betaval,sigma2F = cad3sigma2val, rhoF = cad3rhoval, phiF = cad3phival,
            tau2F = cad3tau2val)






##### transforming each parameter chain into an mcmc object
chain1beta0 = mcmc(cad1$betaF[,1])
chain2beta0 = mcmc(cad2$betaF[,1])
chain3beta0 = mcmc(cad3$betaF[,1])

chain1beta1 = mcmc(cad1$betaF[,2])
chain2beta1 = mcmc(cad2$betaF[,2])
chain3beta1 = mcmc(cad3$betaF[,2])

chain1beta2 = mcmc(cad1$betaF[,3])
chain2beta2 = mcmc(cad2$betaF[,3])
chain3beta2 = mcmc(cad3$betaF[,3])



chain1beta3 = mcmc(cad1$betaF[,4])
chain2beta3 = mcmc(cad2$betaF[,4])
chain3beta3 = mcmc(cad3$betaF[,4])


chain1sigma = mcmc(cad1$sigma2F)
chain2sigma = mcmc(cad2$sigma2F)
chain3sigma = mcmc(cad3$sigma2F)

chain1rho = mcmc(cad1$rhoF)
chain2rho = mcmc(cad2$rhoF)
chain3rho = mcmc(cad3$rhoF)

chain1phi = mcmc(cad1$phiF)
chain2phi = mcmc(cad2$phiF)
chain3phi = mcmc(cad2$phiF)

chain1tau2 = mcmc(cad1$tau2F)
chain2tau2 = mcmc(cad2$tau2F)
chain3tau2 = mcmc(cad2$tau2F)


####### joining the two chains ######
listmcmcbeta0 = list(chain1beta0,chain2beta0,chain3beta0)
listmcmcbeta0 = mcmc.list(listmcmcbeta0)

listmcmcbeta1 = list(chain1beta1,chain2beta1,chain3beta1)
listmcmcbeta1 = mcmc.list(listmcmcbeta1)

listmcmcbeta2 = list(chain1beta2,chain2beta2, chain3beta2)
listmcmcbeta2 = mcmc.list(listmcmcbeta2)

listmcmcbeta3 = list(chain1beta3,chain2beta3, chain3beta3)
listmcmcbeta3 = mcmc.list(listmcmcbeta3)

listmcmcsigma = list(chain1sigma,chain2sigma,chain3sigma)
listmcmcsigma = mcmc.list(listmcmcsigma)

listmcmcrho = list(chain1rho,chain2rho,chain3rho)
listmcmcrho = mcmc.list(listmcmcrho)

listmcmcphi = list(chain1phi,chain2phi,chain3phi)
listmcmcphi = mcmc.list(listmcmcphi)

listmcmctau2 = list(chain1tau2,chain2tau2,chain3tau2)
listmcmctau2 = mcmc.list(listmcmctau2)


###########################################

###### gelman 's convergence criteria ####################

####beta0###
gelman.diag(listmcmcbeta0)
####beta1###
gelman.diag(listmcmcbeta1)
####beta2###
gelman.diag(listmcmcbeta2)
####beta3###
gelman.diag(listmcmcbeta3)
####sigma2###
gelman.diag(listmcmcsigma)
####tau2###
gelman.diag(listmcmctau2)
####rho###
gelman.diag(listmcmcrho)
####gamma (phi)###
gelman.diag(listmcmcphi)

###########################################

####### preparing the data ##########################

####chain 1 ###################
data1 = cbind(cad1$betaF,cad1$sigma2F,cad1$rhoF,cad1$phiF,cad1$tau2F)
data1 = data.frame(data1)
nchain = nrow(data1)
data1 = data.frame(data1,'Chain 1',1:nchain)
colnames(data1) = c(paste0('beta_',0:3),'sigma2','rho','phi','tau2','chain','Index')

####chain 2 ###################
data2 = data.frame(cad2$betaF,cad2$sigma2F,cad2$rhoF,cad2$phiF,cad2$tau2F)
nchain = nrow(data2)
data2 = data.frame(data2,'Chain 2',1:nchain)
colnames(data2) = c(paste0('beta_',0:3),'sigma2','rho','phi','tau2','chain','Index')

####chain 3 ###################
data3 = data.frame(cad3$betaF,cad3$sigma2F,cad3$rhoF,cad3$phiF,cad3$tau2F)
nchain = nrow(data3)
data3 = data.frame(data3,'Chain 3',1:nchain)
colnames(data3) = c(paste0('beta_',0:3),'sigma2','rho','phi','tau2','chain','Index')


###### joining all the chains #######
datatot = rbind(data1,data2,data3)


labelsb = c("Chain 1", "Chain 2", "Chain 3")

###### chain beta0 #######

g1 = datatot%>%ggplot(aes(x=Index,y=beta_0,color=factor(chain)))+ geom_line() +
  theme_bw() +
  scale_colour_grey(name = "Chain",start = 0.3,end = 0.7, labels = labelsb) +
  ylab(expression(beta[0])) + xlab("Index")

ggsave("beta0chain.eps",plot = g1, width = 6, height = 6, dpi = 600)



###### chain beta1 #######
g2 = datatot%>%ggplot(aes(x=Index,y=beta_1,color=factor(chain)))+ geom_line() +
  theme_bw() +
  scale_colour_grey(name = "Chain",start = 0.3,end = 0.7, labels = labelsb) +
  ylab(expression(beta[1])) + xlab("Index")

ggsave("beta1chain.eps",plot = g2, width = 6, height = 6, dpi = 600)


###### chain beta2 #######
g3 = datatot%>%ggplot(aes(x=Index,y=beta_2,color=factor(chain)))+ geom_line() +
  theme_bw() +
  scale_colour_grey(name = "Chain",start = 0.3,end = 0.7, labels = labelsb) +
  ylab(expression(beta[2])) + xlab("Index")

ggsave("beta2chain.eps",plot = g3, width = 6, height = 6, dpi = 600)


###### chain beta3 #######
g4 = datatot%>%ggplot(aes(x=Index,y=beta_3,color=factor(chain)))+ geom_line() +
  theme_bw() +
  scale_colour_grey(name = "Chain",start = 0.3,end = 0.7, labels = labelsb) +
  ylab(expression(beta[3])) + xlab("Index")

ggsave("beta3chain.eps",plot = g4, width = 6, height = 6, dpi = 600)



##### chains betas #########
betaschains = ggarrange(g1,g2,g3,g4,common.legend = TRUE,legend = "bottom")



###### chain sigma2 #######
g5 = datatot%>%ggplot(aes(x=Index,y=sigma2,color=factor(chain)))+ geom_line() +
  theme_bw() +
  scale_colour_grey(name = "Chain",start = 0.3,end = 0.7, labels = labelsb) +
  ylab(expression(sigma^2)) + xlab("Index")

ggsave("sigma2chain.eps",plot = g5, width = 6, height = 6, dpi = 600)


###### chain rho #######
g6 = datatot%>%ggplot(aes(x=Index,y=rho,color=factor(chain)))+ geom_line() +
  theme_bw() +
  scale_colour_grey(name = "Chain",start = 0.3,end = 0.7, labels = labelsb) +
  ylab(expression(rho)) + xlab("Index")

ggsave("rhochain.eps",plot = g6, width = 6, height = 6, dpi = 600)



###### chain gamma #######
g7 = datatot%>%ggplot(aes(x=Index,y=phi,color=factor(chain)))+ geom_line() +
  theme_bw() +
  scale_colour_grey(name = "Chain",start = 0.3,end = 0.7, labels = labelsb) +
  ylab(expression(gamma)) + xlab("Index")

ggsave("gammachain.eps",plot = g7, width = 6, height = 6, dpi = 600)


###### chain tau2 #######
g8 = datatot%>%ggplot(aes(x=Index,y=tau2,color=factor(chain)))+ geom_line() +
  theme_bw() +
  scale_colour_grey(name = "Chain",start = 0.3,end = 0.7, labels = labelsb) +
  ylab(expression(tau^2)) + xlab("Index")

ggsave("tau2chain.eps",plot = g8, width = 6, height = 6, dpi = 600)



##### chains covariance parameters #######
covchains = ggarrange(g5,g6,g7,g8,common.legend = TRUE,legend = "bottom")

covchains

##### acf plots for all the parameters
library(ggfortify)

##### ACF beta0 #####
p1 <- autoplot(acf(cad1$betaF[,1], plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.col=1)

p1 = p1+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(beta[0]))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )

p1

ggsave("acfbeta0.eps",plot = p1, width = 6, height = 6, dpi = 600)



##### ACF beta1 ########
p2 <- autoplot(acf(cad1$betaF[,2], plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.col=1)

p2 = p2+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(beta[1]))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()+ theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )

p2

ggsave("acfbeta1.eps",plot = p2, width = 6, height = 6, dpi = 600)



##### ACF beta2 ########
p3 <- autoplot(acf(cad1$betaF[,3], plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.col=1)

p3 = p3+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(beta[2]))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )


p3

ggsave("acfbeta2.eps",plot = p3, width = 6, height = 6, dpi = 600)



##### ACF beta3 ########
p4 <- autoplot(acf(cad1$betaF[,4], plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.col=1)

p4 = p4+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(beta[3]))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )

p4


ggsave("acfbeta3.eps",plot = p4, width = 6, height = 6, dpi = 600)



betasacf = ggarrange(p1,p2,p3,p4, common.legend = TRUE, legend = "bottom")


##### ACF sigma2 ########
p5 <- autoplot(acf(cad1$sigma2F, plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.col=1)

p5 = p5+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(sigma^2))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )


p5

ggsave("acfsigma2.eps",plot = p5, width = 6, height = 6, dpi = 600)


##### ACF rho ########
p6 <- autoplot(acf(cad1$rhoF, plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.col=1)

p6 = p6+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(rho))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )

p6

ggsave("acfrho.eps",plot = p6, width = 6, height = 6, dpi = 600)



##### ACF gamma ########
p7 <- autoplot(acf(cad1$phiF, plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.col=1)

p7 = p7+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(gamma))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )


p7

ggsave("acfgamma.eps",plot = p7, width = 6, height = 6, dpi = 600)


##### ACF tau2 ########
p8 <- autoplot(acf(cad1$tau2F, plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.col=1)

p8 = p8+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(tau^2))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )


p8

ggsave("acftau2.eps",plot = p8, width = 6, height = 6, dpi = 600)




covaacf = ggarrange(p5,p6,p7,p8, common.legend = TRUE, legend = "bottom")


###### posterior densities beta######

#### density beta0 ########

beta0post = ggplot(datatot, aes(x = beta_0)) +
  # Smooth density curve (darker gray, thinner line)
  geom_histogram(aes(y = after_stat(count / sum(count))),
    fill = "gray50",     # darker gray fill
    color = "black",     # black border
    alpha = 0.8,         # more transparent shading
    linewidth = 0.3
  ) +
  labs(
    x = expression(beta[0]),
    y = "Frequency"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )


beta0post

ggsave("beta0pos.eps",plot = beta0post, width = 6, height = 6, dpi = 600)



#### density beta1 ########

beta1post = ggplot(datatot, aes(x = beta_1)) +
  # Smooth density curve (darker gray, thinner line)
  geom_histogram(aes(y = after_stat(count / sum(count))),
    fill = "gray50",     # darker gray fill
    color = "black",     # black border
    alpha = 0.8,         # more transparent shading
    linewidth = 0.3
  ) +
  labs(
    x = expression(beta[1]),
    y = "Frequency"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )

beta1post

ggsave("beta1pos.eps",plot = beta1post, width = 6, height = 6, dpi = 600)


#### density beta2 ########
beta2post = ggplot(datatot, aes(x = beta_2)) +
  # Smooth density curve (darker gray, thinner line)
  geom_histogram(aes(y = after_stat(count / sum(count))),
    fill = "gray50",     # darker gray fill
    color = "black",     # black border
    alpha = 0.8,         # more transparent shading
    linewidth = 0.3
  ) +
  labs(
    x = expression(beta[2]),
    y = "Frequency"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text( hjust = 0.5, size = 23),
    axis.title = element_text( size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )

beta2post

ggsave("beta2pos.eps",plot = beta2post, width = 6, height = 6, dpi = 600)


#### density beta3 ########

beta3post = ggplot(datatot, aes(x = beta_3)) +
  # Smooth density curve (darker gray, thinner line)
  geom_histogram(aes(y = after_stat(count / sum(count))),
    fill = "gray50",     # darker gray fill
    color = "black",     # black border
    alpha = 0.8,         # more transparent shading
    linewidth = 0.3
  ) +
  labs(
    x = expression(beta[3]),
    y = "Frequency"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text( hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )

ggsave("beta3pos.eps",plot = beta3post, width = 6, height = 6, dpi = 600)


##### density beta
betaposterior = ggarrange(beta0post,beta1post,beta2post,beta3post)



#### density sigma2 ########
sigma2post = ggplot(datatot, aes(x = sigma2)) +
  # Smooth density curve (darker gray, thinner line)
  geom_histogram(aes(y = after_stat(count / sum(count))),
    fill = "gray50",     # darker gray fill
    color = "black",     # black border
    alpha = 0.8,         # more transparent shading
    linewidth = 0.3
  ) +
  labs(
    x = expression(sigma^2),
    y = "Frequency"
  )  +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text( hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )

sigma2post

ggsave("sigma2pos.eps",plot = sigma2post, width = 6, height = 6, dpi = 600)


#### density rho ########
rhopost = ggplot(datatot, aes(x = rho)) +
  # Smooth density curve (darker gray, thinner line)
   geom_histogram(aes(y = after_stat(count / sum(count))),
    fill = "gray50",     # darker gray fill
    color = "black",     # black border
    alpha = 0.8,         # more transparent shading
    linewidth = 0.3
  ) +
  labs(
    x = expression(rho),
    y = "Frequency"
  )  +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text( hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )

rhopost

ggsave("rhopos.eps",plot = rhopost, width = 6, height = 6, dpi = 600)


#### density gamma ########

phipost = ggplot(datatot, aes(x = phi)) +
  # Smooth density curve (darker gray, thinner line)
  geom_histogram(aes(y = after_stat(count / sum(count))),
    fill = "gray50",     # darker gray fill
    color = "black",     # black border
    alpha = 0.8,         # more transparent shading
    linewidth = 0.3
  ) +
  labs(
    x = expression(gamma),
    y = "Frequency"
  )  +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text( hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )

phipost

ggsave("gammapos.eps",plot = phipost, width = 6, height = 6, dpi = 600)


#### density tau2 ########

tau2post = ggplot(datatot, aes(x = tau2)) +
  # Smooth density curve (darker gray, thinner line)
  geom_histogram(aes(y = after_stat(count / sum(count))),
    fill = "gray50",     # darker gray fill
    color = "black",     # black border
    alpha = 0.8,         # more transparent shading
    linewidth = 0.3
  ) +
  labs(
    x = expression(tau^2),
    y = "Frequency"
  )  +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text( hjust = 0.5, size = 23),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )

tau2post

ggsave("tau2pos.eps",plot = tau2post, width = 6, height = 6, dpi = 600)



covaposteriori = ggarrange(sigma2post,rhopost,phipost,tau2post)
ggsave("covaposteriori.eps",plot = covaposteriori, width = 10, height = 10, dpi = 600)

