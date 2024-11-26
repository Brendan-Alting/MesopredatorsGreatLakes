#####Repeating discrete space nimble model with my data: 

library(readxl)
library(tidyverse)
library(sp)
library(sf)
library(nimble)
library(MCMCvis)
library(ggmcmc)
library(gridExtra)

source("R Scripts/3.0 functions for 3.1.r") # required for nimble model

#-----------------------------------
borderhuh <- st_read("Raw Data/studyareaclipped.shp")

bordersimply <- st_buffer(borderhuh, 10) #do by small take too long
bordersimply <- st_union(bordersimply, dTolerance = 10)
bordersimply <- st_buffer(bordersimply, 150) #now do more

cams<- read.csv("Raw Data/camlocs_MLNP.csv")
cam.type <- as.integer(cams$cam.type == "Baited") +0 # 0 for trail, 1 for baited

camlocs<- cams %>% dplyr::select(x=Easting,y=Northing)
camssf <- st_as_sf(camlocs, coords = c("x","y"), crs = st_crs(borderhuh))

camsbufferfox <- st_buffer(camssf, dist = 3000)
camsbufferdingo <- st_buffer(camssf, dist = 8000)
camsbufferquoll <- st_buffer(camssf, dist = 4000)
camsbuffergoanna <- st_buffer(camssf, dist = 3000)

camsbufferfox <- st_union(camsbufferfox, dTolerance = 10)
camsbufferdingo <- st_union(camsbufferdingo, dTolerance = 10)
camsbufferquoll <- st_union(camsbufferquoll,dTolerance = 10)
camsbuffergoanna <- st_union(camsbuffergoanna, dTolerance = 10)


borderfox <- st_intersection(camsbufferfox, bordersimply)
borderdingo <- st_intersection(camsbufferdingo, bordersimply)
borderquoll <- st_intersection(camsbufferquoll, bordersimply)
bordergoanna <- st_intersection(camsbuffergoanna, bordersimply)



###begin with ...

foxobs<- read.csv("Raw Data/foxdetmlnp.csv")
dingoobs <- read.csv("Raw Data/dingodetmlnpnopupnoun.csv")
quollobs <- read.csv("Raw Data/quolldetmlnp.csv")
goannaobs <- read.csv("Raw Data/goannadetmlnp.csv")



#------------------------------------
foxobs<- as.matrix(foxobs)
dingoobs <- as.matrix(dingoobs)
quollobs <- as.matrix(quollobs)
goannaobs <-as.matrix(goannaobs)

cam.type <- as.integer(cams$cam.type == "Baited") +0 # 0 for trail, 1 for baited


#---------------------------------------------

camlocs<- cams %>% dplyr::select(x=Easting,y=Northing)

eps<- 700 # pixel size of state space


#-----------------------------------
# Assemble data ---- FOR FOXES


cam.regionfox<- make_grid(borderfox, cell_diameter=eps, overlap="centre",xy=camlocs)

hex.cent<- st_centroid(cam.regionfox)
X<- as.matrix(camlocs)/1000
grid<- st_coordinates(hex.cent)/1000
npix<- nrow(grid)
pix_area<- as.numeric(st_area(cam.regionfox)[1]/1e6)
Area<- as.numeric(st_area(borderfox)/1e6)




M<-200  # data augmentation - unmarked
n <- rowSums(foxobs, na.rm = TRUE)
get.k <- function(x) length(x[!is.na(x)])
K <- apply(foxobs, 1, get.k)
J <- nrow(foxobs)



constants <- list(M=M,J=J,npix=npix,pix_area=pix_area, A=Area)

data <- list(n=n, grid=grid, X=X, K=K,cam.type = cam.type)


#================================================


code <- nimbleCode({
  
  # Priors
  sigma ~ dgamma(5, 7)
  mean.p ~ dunif(0, 1) # Detection intercept on prob. scale
  alpha0 <- logit(mean.p) # Detection intercept (baseline)
  alpha1 ~ dnorm(0,1) # Effect of bait
  beta0 ~ dnorm(0, sd=2) ##############
  
  for(j in 1:npix) {
    log(muu[j])<- beta0 + log(pix_area)
    ppu[j]<- muu[j]/ENu
  }
  ENu<- sum(muu[1:npix])
  psiu <- ENu/M
  
  # Unmarked part
  for(i in 1:M) {
    z[i] ~ dbern(psiu)
    s[i] ~ dcat(ppu[1:npix])
    g[i,1]<- grid[s[i], 1]
    g[i,2]<- grid[s[i], 2]
    for(j in 1:J) {
      d2[i, j] <- (g[i, 1] - X[j, 1])^2 + (g[i, 2] - X[j, 2])^2
      prob[i, j] <- exp(-d2[i, j] / (2 * sigma^2)) * z[i]
    }
  }
  for (j in 1:J) {
    logit_g0[j] <- alpha0 + alpha1 * cam.type[j]
    g0_eff[j] <- exp(logit_g0[j]) / (1 + exp(logit_g0[j]))
    Ptrap[j] <- 1 - prod(1 - prob[1:M, j] * g0_eff[j])
  }
  n[1:J] ~ dbin_by_row(Ptrap[1:J], K[1:J])
  
  Nu<- sum(z[1:M])
  D<- Nu/A	
})

#------------------------------------------------
# Initial values
sst <- sample(1:npix, M, replace = TRUE)
zust<- rbinom(M,1,0.2)

inits1 <- list(sigma=0.7, alpha1= -2, mean.p=0.25, s=sst, z=zust,beta0 =-1)
             

## create the model object
Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits1, check = FALSE)
Rmcmc<- compileNimble(Rmodel, showCompilerOutput = F)

conf <- configureMCMC(Rmodel)

## add a binary state sampler for each wm and w node
conf$removeSamplers(c("z"), print = FALSE)
Nodes <- Rmodel$expandNodeNames("z")
for(Node in Nodes) conf$addSampler(target = Node, type = "binary", print=FALSE)

conf$removeSamplers(c('alpha1','beta0'), print=FALSE)
conf$addSampler(target=c('alpha1','beta0'), type='AF_slice')

conf$resetMonitors()
conf$addMonitors(c("sigma","beta0","mean.p","psiu","alpha1","D","Nu","s","z","g"))

Cmodel <- buildMCMC(conf)

Cmcmc <- compileNimble(Cmodel, project = Rmodel, resetFunctions = T)

ni<- 12000 #iterations  
nb<- 2000 #burn in
nc<- 3 #chains
nt<- 1 #thinning
#--------------------


inits = function(){inits1}

againszfox<- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, thin = nt, inits = inits,  
               samplesAsCodaMCMC = TRUE)

saveRDS(againszfox, "Derived Data/DiscreteSpaceFox.rds")







#-----------------------------------
# Assemble data ---- FOR Dingoes



cam.regiondingo<- make_grid(borderdingo, cell_diameter=eps, overlap="centre",xy=camlocs)

hex.cent<- st_centroid(cam.regiondingo)
X<- as.matrix(camlocs)/1000
grid<- st_coordinates(hex.cent)/1000
npix<- nrow(grid)
pix_area<- as.numeric(st_area(cam.regiondingo)[1]/1e6)
Area<- as.numeric(st_area(borderdingo)/1e6)


M<-200  # data augmentation - unmarked
n <- rowSums(dingoobs, na.rm = TRUE)
get.k <- function(x) length(x[!is.na(x)])
K <- apply(dingoobs, 1, get.k)
J <- nrow(dingoobs)




constants <- list(M=M,J=J,npix=npix,pix_area=pix_area, A=Area)

data <- list(n=n, grid=grid, X=X, K=K,cam.type = cam.type)


#================================================


code <- nimbleCode({
  
  # Priors
  sigma ~ dgamma(3, 2)
  mean.p ~ dunif(0, 1) # Detection intercept on prob. scale
  alpha0 <- logit(mean.p) # Detection intercept (baseline)
  alpha1 ~ dnorm(0,1) # Effect of bait
  beta0 ~ dnorm(0, sd=2) ##############
  
  for(j in 1:npix) {
    log(muu[j])<- beta0 + log(pix_area)
    ppu[j]<- muu[j]/ENu
  }
  ENu<- sum(muu[1:npix])
  psiu <- ENu/M
  
  # Unmarked part
  for(i in 1:M) {
    z[i] ~ dbern(psiu)
    s[i] ~ dcat(ppu[1:npix])
    g[i,1]<- grid[s[i], 1]
    g[i,2]<- grid[s[i], 2]
    for(j in 1:J) {
      d2[i, j] <- (g[i, 1] - X[j, 1])^2 + (g[i, 2] - X[j, 2])^2
      prob[i, j] <- exp(-d2[i, j] / (2 * sigma^2)) * z[i]
    }
  }
  for (j in 1:J) {
    logit_g0[j] <- alpha0 + alpha1 * cam.type[j]
    g0_eff[j] <- exp(logit_g0[j]) / (1 + exp(logit_g0[j]))
    Ptrap[j] <- 1 - prod(1 - prob[1:M, j] * g0_eff[j])
  }
  n[1:J] ~ dbin_by_row(Ptrap[1:J], K[1:J])
  
  Nu<- sum(z[1:M])
  D<- Nu/A	
})

#------------------------------------------------
# Initial values
sst <- sample(1:npix, M, replace = TRUE)
zust<- rbinom(M,1,0.5)

inits1 <- list(sigma=1.5, alpha1=-2, mean.p=0.3, s=sst, z=zust,beta0 =-1)


## create the model object
Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits1, check = FALSE)
Rmcmc<- compileNimble(Rmodel, showCompilerOutput = F)

conf <- configureMCMC(Rmodel)

## add a binary state sampler for each wm and w node
conf$removeSamplers(c("z"), print = FALSE)
Nodes <- Rmodel$expandNodeNames("z")
for(Node in Nodes) conf$addSampler(target = Node, type = "binary", print=FALSE)

conf$removeSamplers(c('alpha1','beta0'), print=FALSE)
conf$addSampler(target=c('alpha1','beta0'), type='AF_slice')

conf$resetMonitors()
conf$addMonitors(c("sigma","beta0","mean.p","psiu","alpha1","D","Nu","s","z","g"))

Cmodel <- buildMCMC(conf)

Cmcmc <- compileNimble(Cmodel, project = Rmodel, resetFunctions = T)

ni<- 12000  
nb<- 2000
nc<- 3
nt<- 1
#--------------------

#--------------------

inits = function(){inits1}


againszdingo<- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, thin = nt, inits = inits,  
                       samplesAsCodaMCMC = TRUE)



saveRDS(againszdingo, "Derived Data/DiscreteSpaceDingonopupnountest.rds")





#-----------------------------------
# Assemble data ---- FOR Quolls

cam.regionquoll<- make_grid(borderquoll, cell_diameter=eps, overlap="centre",xy=camlocs)

hex.cent<- st_centroid(cam.regionquoll)
X<- as.matrix(camlocs)/1000
grid<- st_coordinates(hex.cent)/1000
npix<- nrow(grid)
pix_area<- as.numeric(st_area(cam.regionquoll)[1]/1e6)
Area<- as.numeric(st_area(borderquoll)/1e6)


M<-200  # data augmentation - unmarked
n <- rowSums(quollobs, na.rm = TRUE)
get.k <- function(x) length(x[!is.na(x)])
K <- apply(quollobs, 1, get.k)
J <- nrow(quollobs)



constants <- list(M=M,J=J,npix=npix,pix_area=pix_area, A=Area)

data <- list(n=n, grid=grid, X=X, K=K,cam.type = cam.type)


#================================================


code <- nimbleCode({
  
  # Priors
  sigma ~ dgamma(5, 7)
  mean.p ~ dunif(0, 1) # Detection intercept on prob. scale
  alpha0 <- logit(mean.p) # Detection intercept (baseline)
  alpha1 ~ dnorm(0,1) # Effect of bait
  beta0 ~ dnorm(0, sd=2) ##############
  
  for(j in 1:npix) {
    log(muu[j])<- beta0 + log(pix_area)
    ppu[j]<- muu[j]/ENu
  }
  ENu<- sum(muu[1:npix])
  psiu <- ENu/M
  
  # Unmarked part
  for(i in 1:M) {
    z[i] ~ dbern(psiu)
    s[i] ~ dcat(ppu[1:npix])
    g[i,1]<- grid[s[i], 1]
    g[i,2]<- grid[s[i], 2]
    for(j in 1:J) {
      d2[i, j] <- (g[i, 1] - X[j, 1])^2 + (g[i, 2] - X[j, 2])^2
      prob[i, j] <- exp(-d2[i, j] / (2 * sigma^2)) * z[i]
    }
  }
  for (j in 1:J) {
    logit_g0[j] <- alpha0 + alpha1 * cam.type[j]
    g0_eff[j] <- exp(logit_g0[j]) / (1 + exp(logit_g0[j]))
    Ptrap[j] <- 1 - prod(1 - prob[1:M, j] * g0_eff[j])
  }
  n[1:J] ~ dbin_by_row(Ptrap[1:J], K[1:J])
  
  Nu<- sum(z[1:M])
  D<- Nu/A	
})

#------------------------------------------------
# Initial values
sst <- sample(1:npix, M, replace = TRUE)
zust<- rbinom(M,1,0.5)

inits1 <- list(sigma=0.5, alpha1=2, mean.p=0.05, s=sst, z=zust,beta0 =-1)


## create the model object
Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits1, check = FALSE)
Rmcmc<- compileNimble(Rmodel, showCompilerOutput = F)

conf <- configureMCMC(Rmodel)

## add a binary state sampler for each wm and w node
conf$removeSamplers(c("z"), print = FALSE)
Nodes <- Rmodel$expandNodeNames("z")
for(Node in Nodes) conf$addSampler(target = Node, type = "binary", print=FALSE)

conf$removeSamplers(c('alpha1','beta0'), print=FALSE)
conf$addSampler(target=c('alpha1','beta0'), type='AF_slice')

conf$resetMonitors()
conf$addMonitors(c("sigma","beta0","mean.p","psiu","alpha1","D","Nu","s","z","g"))

Cmodel <- buildMCMC(conf)

Cmcmc <- compileNimble(Cmodel, project = Rmodel, resetFunctions = T)

ni<- 12000  # for testing - run more iterations to improve convergence
nb<- 2000 # same
nc<- 3
nt<- 1
#--------------------


inits = function(){inits1}

againszquoll<- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, thin = nt, inits = inits,  
                       samplesAsCodaMCMC = TRUE)

saveRDS(againszquoll, "Derived Data/DiscreteSpaceQuoll.rds")




#-----------------------------------
# Assemble data ---- FOR Goannas
eps<- 1195 # pixel size of state space  ###Works!


cam.regiongoanna<- make_grid(bordergoanna, cell_diameter=eps, overlap="centre",xy=camlocs)

hex.cent<- st_centroid(cam.regiongoanna)
X<- as.matrix(camlocs)/1000
grid<- st_coordinates(hex.cent)/1000
npix<- nrow(grid)
pix_area<- as.numeric(st_area(cam.regiongoanna)[1]/1e6)
Area<- as.numeric(st_area(bordergoanna)/1e6)





M<-800  # data augmentation - unmarked
n <- rowSums(goannaobs, na.rm = TRUE)
get.k <- function(x) length(x[!is.na(x)])
K <- apply(goannaobs, 1, get.k)
J <- nrow(goannaobs)



constants <- list(M=M,J=J,npix=npix,pix_area=pix_area, A=Area)

data <- list(n=n, grid=grid, X=X, K=K,cam.type = cam.type)


#================================================


code <- nimbleCode({
  
  # Priors
  sigma ~ dgamma(2,4)
  mean.p ~ dunif(0, 1) # Detection intercept on prob. scale
  alpha0 <- logit(mean.p) # Detection intercept (baseline)
  alpha1 ~ dnorm(0,1) # Effect of bait
  beta0 ~ dnorm(0, sd=2) ##############
  
  for(j in 1:npix) {
    log(muu[j])<- beta0 + log(pix_area)
    ppu[j]<- muu[j]/ENu
  }
  ENu<- sum(muu[1:npix])
  psiu <- ENu/M
  
  # Unmarked part
  for(i in 1:M) {
    z[i] ~ dbern(psiu)
    s[i] ~ dcat(ppu[1:npix])
    g[i,1]<- grid[s[i], 1]
    g[i,2]<- grid[s[i], 2]
    for(j in 1:J) {
      d2[i, j] <- (g[i, 1] - X[j, 1])^2 + (g[i, 2] - X[j, 2])^2
      prob[i, j] <- exp(-d2[i, j] / (2 * sigma^2)) * z[i]
    }
  }
  for (j in 1:J) {
    logit_g0[j] <- alpha0 + alpha1 * cam.type[j]
    g0_eff[j] <- exp(logit_g0[j]) / (1 + exp(logit_g0[j]))
    Ptrap[j] <- 1 - prod(1 - prob[1:M, j] * g0_eff[j])
  }
  n[1:J] ~ dbin_by_row(Ptrap[1:J], K[1:J])
  
  Nu<- sum(z[1:M])
  D<- Nu/A	
})

#------------------------------------------------
# Initial values
sst <- sample(1:npix, M, replace = TRUE)
zust<- rbinom(M,1,0.5)

inits1 <- list(sigma=0.4, alpha1=2, mean.p=0.05, s=sst, z=zust,beta0 =-1)


## create the model object
Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits1, check = FALSE)
Rmcmc<- compileNimble(Rmodel, showCompilerOutput = F)

conf <- configureMCMC(Rmodel)

## add a binary state sampler for each wm and w node
conf$removeSamplers(c("z"), print = FALSE)
Nodes <- Rmodel$expandNodeNames("z")
for(Node in Nodes) conf$addSampler(target = Node, type = "binary", print=FALSE)

conf$removeSamplers(c('alpha1','beta0'), print=FALSE)
conf$addSampler(target=c('alpha1','beta0'), type='AF_slice')

conf$resetMonitors()
conf$addMonitors(c("sigma","beta0","mean.p","psiu","alpha1","D","Nu","s","z","g"))

Cmodel <- buildMCMC(conf)

Cmcmc <- compileNimble(Cmodel, project = Rmodel, resetFunctions = T)

ni<- 12000  # for testing - run more iterations to improve convergence
nb<- 2000 # same
nc<- 3
nt<- 1
#--------------------


inits = function(){inits1}

againszgoanna<- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, thin = nt, inits = inits,  
                       samplesAsCodaMCMC = TRUE)

#start 1:24 Friyay- each bar equals 1/56th of the entire chain. So 1 bar time = 1/168th of total time. 
#Do minutes * 168 = est length. 

saveRDS(againszgoanna, "Derived Data/DiscreteSpaceGoanna.rds")



#end