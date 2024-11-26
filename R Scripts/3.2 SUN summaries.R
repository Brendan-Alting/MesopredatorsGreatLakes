#3.2. Getting summary statistics and evaluation plots from bayesian models. We need to calculate posterior mode for density first.
library(MCMCvis)
library(coda)
library(ggmcmc)
library(gridExtra)

foxesfit<- readRDS("Derived Data/DiscreteSpaceFox.rds") 
dingoesfit<- readRDS("Derived Data/DiscreteSpaceDingo.rds") 
dingoesfitnopupnoun<- readRDS("Derived Data/DiscreteSpaceDingonopupnountest.rds") 
quollsfit<- readRDS("Derived Data/DiscreteSpaceQuoll.rds") 
goannasfit <- readRDS("Derived Data/DiscreteSpaceGoanna.rds")

#lets get the posterior distribution first.

#fox 

foxsamples <- as.matrix(foxesfit)[,"D"]

foxd <- density(foxsamples)
mode_foxd <- foxd$x[which.max(foxd$y)]

mode_foxd

#dingo 

dingosamples <- as.matrix(dingoesfitnopupnoun)[,"D"]

dingod <- density(dingosamples)
mode_dingod <- dingod$x[which.max(dingod$y)]

mode_dingod

#quoll

quollsamples <- as.matrix(quollsfit)[,"D"]

quolld <- density(quollsamples)
mode_quolld <- quolld$x[which.max(quolld$y)]

mode_quolld

#goanna

goannasamples <- as.matrix(goannasfit)[,"D"]

goannad <- density(goannasamples)
mode_goannad <- goannad$x[which.max(goannad$y)]

mode_goannad


#diagnostic plots

#fox

foxsummary <- data.frame(MCMCsummary(foxesfit, c("alpha1","mean.p","psiu","sigma","Nu","D")))
write.csv(foxsummary, file = "Derived Data/Supplementary/Fox Summary.csv")

foxS<- ggs(foxesfit)

p1fox<- ggs_density(foxS, family="Nu") + xlim(0, 200)
p2fox<- ggs_density(foxS, family="D") + xlim(0, 1.1)
p4fox<- ggs_density(foxS, family="sigma") + xlim(0, 1)
p6fox<-ggs_density(foxS, family="mean.p") + xlim(0, 0.75)

grid.arrange(p1fox, p2fox,p4fox,p6fox, ncol=2, nrow=2)


dp1fox<- ggs_traceplot(foxS, family="Nu")
dp2fox<- ggs_traceplot(foxS, family="D")
dp4fox<- ggs_traceplot(foxS, family="sigma")
dp5fox<- ggs_traceplot(foxS, family="psiu")
dp6fox<- ggs_traceplot(foxS, family="mean.p")
dp3fox<- ggs_traceplot(foxS,family ="beta0")



grid.arrange(dp1fox, dp2fox,dp3fox,dp5fox, dp4fox,dp6fox, ncol=3, nrow=2)

#dingo

dingosummarynopup <- data.frame(MCMCsummary(dingoesfitnopupnoun, c("alpha1","mean.p","psiu","sigma","Nu","D")))
write.csv(dingosummarynopup, file = "Derived Data/Supplementary/dingo Summary3.csv")

dingoS<- ggs(dingoesfitnopupnoun)

p1dingo<- ggs_density(dingoS, family="Nu") + xlim(0, 220)
p2dingo<- ggs_density(dingoS, family="D") + xlim(0, 1)
p4dingo<- ggs_density(dingoS, family="sigma") + xlim(0, 3)
p6dingo<-ggs_density(dingoS, family="mean.p") + xlim(0, 0.75)

grid.arrange(p1dingo, p2dingo,p4dingo,p6dingo, ncol=2, nrow=2)


dp1dingo<- ggs_traceplot(dingoS, family="Nu")
dp2dingo<- ggs_traceplot(dingoS, family="D")
dp4dingo<- ggs_traceplot(dingoS, family="sigma")
dp5dingo<- ggs_traceplot(dingoS, family="psiu")
dp6dingo<- ggs_traceplot(dingoS, family="mean.p")
dp3dingo<- ggs_traceplot(dingoS,family ="beta0")



grid.arrange(dp1dingo, dp2dingo,dp3dingo,dp5dingo, dp4dingo,dp6dingo, ncol=3, nrow=2)

#Quoll

quollsummary <- data.frame(MCMCsummary(quollsfit, c("alpha1","mean.p","psiu","sigma","Nu","D")))
write.csv(quollsummary, file = "Derived Data/Supplementary/quoll Summary.csv")

quollS<- ggs(quollsfit)

p1quoll<- ggs_density(quollS, family="Nu") + xlim(0, 200)
p2quoll<- ggs_density(quollS, family="D") + xlim(0, 1.2)
p4quoll<- ggs_density(quollS, family="sigma") + xlim(0, 1)
p6quoll<-ggs_density(quollS, family="mean.p") + xlim(0, 0.15)

grid.arrange(p1quoll, p2quoll,p4quoll,p6quoll, ncol=2, nrow=2)


dp1quoll<- ggs_traceplot(quollS, family="Nu")
dp2quoll<- ggs_traceplot(quollS, family="D")
dp4quoll<- ggs_traceplot(quollS, family="sigma")
dp5quoll<- ggs_traceplot(quollS, family="psiu")
dp6quoll<- ggs_traceplot(quollS, family="mean.p")
dp3quoll<- ggs_traceplot(quollS,family ="beta0")



grid.arrange(dp1quoll, dp2quoll,dp3quoll,dp5quoll, dp4quoll,dp6quoll, ncol=3, nrow=2)

#goanna

goannasummary <- data.frame(MCMCsummary(goannasfit, c("alpha1","mean.p","psiu","sigma","Nu","D")))
write.csv(goannasummary, file = "Derived Data/Supplementary/goanna Summary.csv")

goannaS<- ggs(goannasfit)

p1goanna<- ggs_density(goannaS, family="Nu") + xlim(0, 600)
p2goanna<- ggs_density(goannaS, family="D") + xlim(0, 3.2)
p4goanna<- ggs_density(goannaS, family="sigma") + xlim(0.25, 0.85)
p6goanna<-ggs_density(goannaS, family="mean.p") + xlim(0, 0.06)

grid.arrange(p1goanna, p2goanna,p4goanna,p6goanna, ncol=2, nrow=2)


dp1goanna<- ggs_traceplot(goannaS, family="Nu")
dp2goanna<- ggs_traceplot(goannaS, family="D")
dp4goanna<- ggs_traceplot(goannaS, family="sigma")
dp5goanna<- ggs_traceplot(goannaS, family="psiu")
dp6goanna<- ggs_traceplot(goannaS, family="mean.p")
dp3goanna<- ggs_traceplot(goannaS,family ="beta0")



grid.arrange(dp1goanna, dp2goanna,dp3goanna,dp5goanna, dp4goanna,dp6goanna, ncol=3, nrow=2)
