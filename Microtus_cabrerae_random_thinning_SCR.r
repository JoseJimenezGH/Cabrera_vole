#==============================================================================#
#                                                                              #
#                    MULTISESSION RANDOM THINNING-SCR MODEL                    #
#                           Microtus cabrerae                                  #
#                             Jose Jimenez                                     #
#                         25/04/2023 10:06:36                                  #
#                                                                              #
#==============================================================================#
#rm(list=ls())

library(sf)
library(raster)
library(nimble)
library(basicMCMCplots)
library(coda)
library(lattice)
library(secr)
library(jagsUI)
library(makeJAGSmask)
library(rgdal)
library(scrbook)
library(mcmcOutput)
library(MCMCvis)
library(tidyverse)
library(ggplot2)
library(terra)
library(fields)

setwd('C:/Users/Usuario/OneDrive/40 Proyecto Topillo Cabrera/05 R code/DEF')
source("SCR_Functions.R")


# SITE 2
#========
# Elevation
elev1 <- raster("./gisData/Sitio2.tif")
plot(elev1, asp=TRUE)
ras1<-stack(elev1)
names(ras1)<-c('Elev')
ras1<-scale(ras1)

# Identified samples
mc.ch1 <- read.capthist("./dataBase/capt2.txt", "./dataBase/trap2.txt", detector='count', cov="sex", noccasions=1)
summary(mc.ch1)
cabrera1<-aperm(mc.ch1,c(1,3,2))

M1<-200
a1<-attributes(mc.ch1)
sex1<-a1$covariates[]
sex1.indicator <- as.numeric(as.factor(sex1$sex))
sex1<- sex1.indicator

traplocs1<-as.matrix(secr::traps(mc.ch1))
X1<-data.matrix(traplocs1)
rownames(X1)<-1:469
colnames(X1)<-c("X","Y")

eff<-read.table("./dataBase/Eff2.txt", header=TRUE)[,2]
effs1<-(eff-mean(eff))/sd(eff)

nind1<-dim(cabrera1)[1]
J1<-nrow(X1)
y1<-apply(cabrera1,c(1,2),sum)
y1[y1>0]<-1
Yaug1<-array(0,c(M1,J1))
Yaug1[1:nind1,]<-y1
SEX.a <- c(sex1 - 1, rep(NA, (M1-nind1)))
KT1<-rep(1,J1)

# Bundle data
mymask1 <- convertRaster(ras1, as.data.frame(X1))
str(mymask1)

# Non-ID samples
mc.ch1ni <- read.capthist("./dataBase/nnid2.txt", "./dataBase/trap2.txt", detector='count', noccasions=1)
summary(mc.ch1ni)
cabrera1ni<-aperm(mc.ch1ni,c(1,3,2))
nnid1<-apply(cabrera1ni,c(2,3),sum)

################################################################################

# SITE 4
#========
# Elevation
elev2 <- raster("./gisData/Sitio4.tif")
plot(elev2, asp=TRUE)
ras2<-stack(elev2)
names(ras2)<-c('Elev')
ras2<-scale(ras2)

# Identified samples
mc.ch2 <- read.capthist("./dataBase/capt4.txt", "./dataBase/trap4.txt", detector='count', cov="sex", noccasions=1)
summary(mc.ch2)
cabrera2<-aperm(mc.ch2,c(1,3,2))

M2<-200
a2<-attributes(mc.ch2)
sex2<-a2$covariates[]
sex2.indicator <- as.numeric(as.factor(sex2$sex))
sex2<- sex2.indicator

traplocs2<-as.matrix(secr::traps(mc.ch2))
X2<-data.matrix(traplocs2)
rownames(X2)<-1:932
colnames(X2)<-c("X","Y")

# Sampling effort
eff<-read.table("./dataBase/Eff4.txt", header=TRUE)[,2]
effs2<-(eff-mean(eff))/sd(eff)

nind2<-dim(cabrera2)[1]
J2<-nrow(X2)
y2<-apply(cabrera2,c(1,2),sum)
Yaug2<-array(0,c(M2,J2))
Yaug2[1:nind2,]<-y2
SEX.b <- c(sex2 - 1, rep(NA, (M2-nind2)))

# Bundle data
mymask2 <- convertRaster(ras2, as.data.frame(X2))
str(mymask2)

# Non-ID samples
mc.ch2ni <- read.capthist("./dataBase/nnid4.txt", "./dataBase/trap4.txt", detector='count', noccasions=1)
summary(mc.ch2ni)
cabrera2ni<-aperm(mc.ch2ni,c(1,3,2))
nnid2<-apply(cabrera2ni,c(2,3),sum)

##############################  Capture plots   ################################

all<-vect("./gisData/TODOS.shp")
ni2<-vect("./gisData/nnid2.shp")
ni4<-vect("./gisData/nnid4.shp")
traps2<-vect("./gisData/traps2.shp")
traps4<-vect("./gisData/traps4.shp")

dev.new(width=9.979167,height=5.541667)
plot(elev1, asp=TRUE)
plot(traps2, add=TRUE, border='grey60', lwd=0.1)
points(all)
points(ni2,col="red")
sbar(100, xy="bottomright", divs=4, cex=.8, ticks=TRUE)
text(181731.7,4391544, "m")

dev.new(width=7.760417,height=9.104167)
plot(elev2, asp=TRUE)
plot(traps4, add=TRUE, border='grey60', lwd=0.1)
points(all)
points(ni4,col="red")
sbar(100, xy="bottomright", divs=4, cex=.8, ticks=TRUE)
text(182947.3,4390104, "m")

################################################################################

## define the model
code <- nimbleCode({

  for(i in 1:nSites){
    psi[i] ~ dunif(0,1)
    beta01[i] ~ dunif(-10, 10)          # Elevation
    beta02[i] ~ dunif(-10, 10)          # Elevation (^2)
    alpha1.pS[i] ~ dunif(-10, 10)       # Effort
    alpha2.pS[i] ~ dunif(-10, 10)       # Effort (^2)
    psi.sex[i] ~ dunif(0,1)
    lp0.sex[i] ~ dunif(-10, 10)
  }

  id.prob1 ~ dunif(0,1)
  id.prob2 ~ dunif(0,1)

  for(i in 1:nSex){
    sigma[i] ~ dunif(0,10)             # scaled sigma
    sigmaR[i]<-sigma[i]*pixelWidth     # real sigma
  }

  # Site 2
  #=========
  # Detection probability covaried with effort and sex (as a random effect)
  for(j in 1:J1){
    log(p0S1[1,j]) <-  alpha1.pS[1]*Eff1[j] + alpha2.pS[1]*(Eff1[j])^2 + lp0.sex[1]
    log(p0S1[2,j]) <-  alpha1.pS[1]*Eff1[j] + alpha2.pS[1]*(Eff1[j])^2 - lp0.sex[1]
  }

  # Spatial covariates for density
  for(i in 1:(upperLimit1[1]-1)) {
    for(j in 1:(upperLimit1[2]-1)) {
      log(lam1[i, j]) <- beta01[1]*Elev1[i,j]+beta02[1]*(Elev1[i,j])^2
    }
  }
  # convert 'lam01' to 0 for non-habitat
  lam01[(1:(upperLimit1[1]-1)), (1:(upperLimit1[2]-1))] <- lam1[(1:(upperLimit1[1]-1)), (1:(upperLimit1[2]-1))]  * habMat1[(1:(upperLimit1[1]-1)), (1:(upperLimit1[2]-1))]
  # 'probs1' must sum to 1
  probs1[(1:(upperLimit1[1]-1)), (1:(upperLimit1[2]-1))]  <- lam01[(1:(upperLimit1[1]-1)), (1:(upperLimit1[2]-1))]  / sum(lam01[(1:(upperLimit1[1]-1)), (1:(upperLimit1[2]-1))] )

  for (i in 1:M1){                       # augmentation
    z1[i] ~ dbern(psi[1])                # it's the individual real?
    SEX.a[i]~ dbern(psi.sex[1])          # Potential males
    SEX2.a[i]<-SEX.a[i] + 1
    males1[i] <- z1[i] * SEX.a[i]	       # Realized males
    sexfemale1[i] <- 1-SEX.a[i]	         # Potential females
    females1[i] <- z1[i] * sexfemale1[i] # Realized females
    S1[i, 1] ~ dunif(1, upperLimit1[1])  # uniform priors for the activity centres for each individual
    S1[i, 2] ~ dunif(1, upperLimit1[2])
    negLogDen1[i] <- -log(probs1[trunc(S1[i,1]), trunc(S1[i,2])]) # zeros trick
    zeros1[i] ~ dpois(negLogDen1[i])
    Dsq1[i,1:J1] <- (S1[i,1]-trapMat1[1:J1,1])^2 + (S1[i,2]-trapMat1[1:J1,2])^2
    lamS1[i,1:J1] <- p0S1[SEX2.a[i],1:J1] * exp(-Dsq1[i,1:J1]/(2*sigma[SEX2.a[i]]^2))*z1[i]
    y.fullS1[i,1:J1] ~ dPoissonVector(lamS1[i,1:J1])  # Vectorized Poisson

    for (j in 1:J1) {
      y.obs1[i,j] ~ dbin(id.prob1, y.fullS1[i,j])  # Id. samples over total (yfullS1) with prob=id.prob1
    }
  }

  # Site 4
  #=========
  # Detection probability covaried with effort and sex (as a random effect)
  for(j in 1:J2){
    log(p0S2[1,j]) <-  alpha1.pS[2]*Eff2[j] + alpha2.pS[2]*(Eff2[j])^2 + lp0.sex[2]
    log(p0S2[2,j]) <-  alpha1.pS[2]*Eff2[j] + alpha2.pS[2]*(Eff2[j])^2 - lp0.sex[2]
  }

  # Spatial covariates for density
  for(i in 1:(upperLimit2[1]-1)) {
    for(j in 1:(upperLimit2[2]-1)) {
      log(lam2[i, j]) <- beta01[2]*Elev2[i,j]+beta02[2]*(Elev2[i,j])^2
    }
  }
  # convert 'lam02' to 0 for non-habitat
  lam02[(1:(upperLimit2[1]-1)), (1:(upperLimit2[2]-1))] <- lam2[(1:(upperLimit2[1]-1)), (1:(upperLimit2[2]-1))]  * habMat2[(1:(upperLimit2[1]-1)), (1:(upperLimit2[2]-1))]
  # 'probs2' must sum to 1
  probs2[(1:(upperLimit2[1]-1)), (1:(upperLimit2[2]-1))]  <- lam02[(1:(upperLimit2[1]-1)), (1:(upperLimit2[2]-1))]  / sum(lam02[(1:(upperLimit2[1]-1)), (1:(upperLimit2[2]-1))] )

  for (i in 1:M2){                       # augmentation
    z2[i] ~ dbern(psi[2])                # it's the individual real?
    SEX.b[i]~ dbern(psi.sex[2])          # Potential males
    SEX2.b[i]<-SEX.b[i] + 1
    males2[i] <- z2[i] * SEX.b[i]	     # Realized males
    sexfemale2[i] <- 1-SEX.b[i]	         # Potential females
    females2[i] <- z2[i] * sexfemale2[i] # Realized females
    S2[i, 1] ~ dunif(1, upperLimit2[1])  # uniform priors for the activity centres for each individual
    S2[i, 2] ~ dunif(1, upperLimit2[2])
    negLogDen2[i] <- -log(probs2[trunc(S2[i,1]), trunc(S2[i,2])]) # zeros trick
    zeros2[i] ~ dpois(negLogDen2[i])
    Dsq2[i,1:J2] <- (S2[i,1]-trapMat2[1:J2,1])^2 + (S2[i,2]-trapMat2[1:J2,2])^2
    lamS2[i,1:J2] <- p0S2[SEX2.b[i],1:J2] * exp(-Dsq2[i,1:J2]/(2*sigma[SEX2.b[i]]^2))*z2[i]
    y.fullS2[i,1:J2] ~ dPoissonVector(lamS2[i,1:J2])  # Vectorized Poisson

    for (j in 1:J2) {
      y.obs2[i,j] ~ dbin(id.prob2, y.fullS2[i,j])  # Id. samples over total (yfullS2) with prob=id.prob2
    }
  }

  N1 <- sum(z1[1:M1])               # Realized number of individuals
  Nmales1 <- sum(males1[1:M1]) 	    # Realized number of males
  Nfemales1 <- sum(females1[1:M1])  # Realized number of females
  SR1 <- Nmales1 / N1			    # Male sex ratio
  D1<- N1/area1
  N2 <- sum(z2[1:M2])               # Realized number of individuals
  Nmales2 <- sum(males2[1:M2]) 	    # Realized number of males
  Nfemales2 <- sum(females2[1:M2])  # Realized number of females
  SR2 <- Nmales2 / N2			    # Male sex ratio
  D2<- N2/area2

})

# Vectorized Poisson distribution
dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lambda = double(1),
  log = integer(0, default = 0)) {
    J <- length(x)
    ans <- 0.0
    for(j in 1:J)
      ans <- ans + dpois(x[j], lambda[j], 1)
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rPoissonVector  <- nimbleFunction(
  run = function(n = integer(), lambda = double(1)) {
    J <- length(lambda)
    ans<- numeric(J)
    for(j in 1:J)
      ans[j] <- rpois(1, lambda[j])
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(
  dPoissonVector = list(
    BUGSdist = "dPoissonVector(lambda)",
    Rdist = "dPoissonVector(lambda)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = double(1)', 'lambda = double(1)'))
))




# Data
SEX.a <- c(sex1 - 1, rep(NA, (M1-nind1)))
SEX.b <- c(sex2 - 1, rep(NA, (M2-nind2)))
str(data  <-   list(# SITE 2
                    #===================================================
                    y.obs1 = Yaug1,                  # ID individuals
                    zeros1=rep(0,M1),                # zero trick
                    Eff1=effs1,                      # effort
                    Elev1=mymask1$covMat,            # elevation data
                    SEX.a=SEX.a,                     # sex
                    habMat1=mymask1$habMat,          # habitat available
                    trapMat1=mymask1$trapMat,        # trap
                    # SITE 4
                    #===================================================
                    y.obs2 = Yaug2,                  # ID individuals
                    zeros2=rep(0,M2),                # zero trick
                    Eff2=effs2,                      # effort
                    Elev2=mymask2$covMat,            # elevation data
                    SEX.b=SEX.b,                     # sex
                    habMat2=mymask2$habMat,          # habitat available
                    trapMat2=mymask2$trapMat))       # trap

# Constants
J1 <- nrow(X1)
J2 <- nrow(X2)
nnidd1<-as.numeric(apply(nnid1,1,sum))
nnidd2<-as.numeric(apply(nnid2,1,sum))
str(constants<-list(nSex=2,                          # number of sexes
                    nSites=2,                        # number of sites
                    pixelWidth=mymask1$pixelWidth,   # pixel size
                    # SITE 2
                    #===================================================
                    M1=M1,                           # data augmentation
                    J1=469,                          # number of traps
                    area1=0.8512,                    # area (ha)
                    nnidd1=nnidd1,                   # Non-ID data
                    upperLimit1=mymask1$upperLimit,  # upper limit
                    # SITE 4
                    #===================================================
                    M2=M2,                           # data augmentation
                    J2=932,                          # number of traps
                    area2=1.4876,                    # area (ha)
                    nnidd2=nnidd2,                   # Non-ID data
                    upperLimit2=mymask2$upperLimit)) # upper limit

## Inits
load('./dataBase/InitsMultElev.RData')

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
Rmodel$initializeInfo()
Rmodel$calculate()
Cmodel <- compileNimble(Rmodel)
params<-c('N1','Nmales1','Nfemales1','D1','SR1',
          'N2','Nmales2','Nfemales2','D2','SR2',
          'psi.sex','psi',
          'alpha1.pS','alpha2.pS','lp0.sex',
          'beta01','beta02',
          'SEX.a','SEX.b',
          'sigmaR','S1','z1','S2','z2')

conf<-configureMCMC(Rmodel, monitors=params)

#### --- load custom samplers for random thinning SCR  --- ####
# replace with new sampler for y (sample without replacement with sum to n[j,k] constraint)
conf$removeSampler("y.fullS1")     # to deal with non-ID data from Site 2
for(j in 1:J1){
  conf$addSampler(target = paste(paste("y.fullS1[1:",M1,", ",j,"]"), sep=""),
                  type = 'IDSampler1',control = list(nnidd1 = nnidd1[j], j=j, M1=M1),
                  silent = TRUE)
  }
conf$removeSampler("y.fullS2")     # to deal with non-ID data from Site 4
for(j in 1:J2){
  conf$addSampler(target = paste(paste("y.fullS2[1:",M2,", ",j,"]"), sep=""),
                  type = 'IDSampler2',control = list(nnidd2 = nnidd2[j], j=j, M2=M2),
                  silent = TRUE)
  }

MCMC <- buildMCMC(conf)
CompMCMC <- compileNimble(MCMC, project = Rmodel)

## Execute MCMC algorithm and extract samples
## Run model
nb=10000       # Burnin
ni=50000 +nb   # Iters
nc=3           # Chains

start.time2<-Sys.time()
outNim <- runMCMC(CompMCMC, niter = ni , nburnin = nb , nchains = nc, inits=inits,
                  setSeed = TRUE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
end.time<-Sys.time()
end.time-start.time2 # post-compilation run time

save(outNim, file="outNim.RData")

summary(mcmcOutput(outNim[,c('N1','Nmales1','Nfemales1','D1','SR1','psi[1]','psi.sex[1]',
                             'N2','Nmales2','Nfemales2','D2','SR2','psi[2]','psi.sex[2]',
                             'sigmaR[1]','sigmaR[2]',
                             'alpha1.pS[1]','alpha2.pS[1]','lp0.sex[1]',
                             'alpha1.pS[2]','alpha2.pS[2]','lp0.sex[2]',
                             'beta01[1]','beta02[1]','beta01[2]','beta02[2]')]))

# MCMC values from mcmc.list object ‘outNim[, c("N1", "Nmales1", "Nfemales1", "D1", "SR1", "psi[1]", ’ MCMC values from mcmc.list object ‘    "psi.sex[1]", "N2", "Nmales2", "Nfemales2", "D2", "SR2", ’ MCMC values from mcmc.list object ‘    "psi[2]", "psi.sex[2]", "sigmaR[1]", "sigmaR[2]", "alpha1.pS[1]", ’ MCMC values from mcmc.list object ‘    "alpha2.pS[1]", "lp0.sex[1]", "alpha1.pS[2]", "alpha2.pS[2]", ’ MCMC values from mcmc.list object ‘    "lp0.sex[2]", "beta01[1]", "beta02[1]", "beta01[2]", "beta02[2]")]’ 
# The object has 26 nodes with 50000 draws for each of 3 chains.
# l95 and u95 are the limits of a 95% Highest Density Credible Interval.
# Rhat is the estimated potential scale reduction factor:
        # largest is 1.00; NONE are greater than 1.10.
# MCEpc is the Monte Carlo standard error as a percentage of the posterior SD:
        # largest is 2.6%; NONE are greater than 5%.

                # mean     sd  median     l95     u95  Rhat MCEpc
# N1           119.020 11.742 118.000  96.000 141.000 0.978 1.965
# Nmales1       40.668  5.668  40.000  30.000  51.000 0.977 1.209
# Nfemales1     78.352 11.481  77.000  56.000 100.000 0.977 1.840
# D1           139.826 13.794 138.628 112.782 165.648 0.978 1.965
# SR1            0.344  0.049   0.342   0.250   0.439 1.000 1.319
# psi[1]         0.594  0.067   0.591   0.465   0.729 0.998 1.820
# psi.sex[1]     0.346  0.064   0.344   0.224   0.472 0.999 1.337
# N2           136.064 14.791 134.000 107.000 164.000 0.982 2.560
# Nmales2       41.396  6.552  41.000  30.000  54.000 1.000 1.571
# Nfemales2     94.668 14.326  93.000  69.000 123.000 1.000 2.345
# D2            91.466  9.943  90.078  72.600 110.917 0.982 2.560
# SR2            0.306  0.049   0.303   0.214   0.403 0.999 1.543
# psi[2]         0.679  0.080   0.673   0.525   0.838 0.999 2.441
# psi.sex[2]     0.309  0.062   0.305   0.190   0.431 1.000 1.493
# sigmaR[1]      2.177  0.139   2.171   1.908   2.451 1.001 2.409
# sigmaR[2]      2.752  0.195   2.744   2.379   3.138 1.001 2.268
# alpha1.pS[1]   0.387  0.098   0.386   0.201   0.584 1.000 1.149
# alpha2.pS[1]  -0.183  0.058  -0.181  -0.295  -0.071 0.999 1.174
# lp0.sex[1]     0.008  0.139   0.008  -0.263   0.277 1.001 1.758
# alpha1.pS[2]   0.587  0.105   0.586   0.380   0.790 1.000 1.210
# alpha2.pS[2]  -0.195  0.047  -0.194  -0.286  -0.103 1.000 1.188
# lp0.sex[2]    -0.018  0.140  -0.018  -0.291   0.260 1.002 1.817
# beta01[1]     -0.863  0.380  -0.830  -1.632  -0.173 1.000 2.265
# beta02[1]     -1.295  0.470  -1.262  -2.226  -0.417 1.000 2.313
# beta01[2]      0.239  0.186   0.235  -0.131   0.603 1.000 1.741
# beta02[2]     -0.184  0.187  -0.159  -0.551   0.149 1.000 2.441
