#Run SECR models for quolls

library(secr)
options(scipen=999)

#First check detection function 

nullHRquoll <- secr.fit(quollsdethist, 
                   mask = correctmaskclipped,
                   link = "log",
                   detectfn = "HR",
                   ncores = 4,
                   trace = TRUE,
                   model = list(D ~1, sigma ~1, g0 ~ 1),
                   details = list(fastproximity = FALSE))

nullEXquoll <- secr.fit(quollsdethist, 
                   mask = correctmaskclipped,
                   link = "log",
                   detectfn = "EX",
                   ncores = 4,
                   trace = TRUE,
                   model = list(D ~1, sigma ~1, g0 ~ 1),
                   details = list(fastproximity = FALSE))

nullHNquoll <- secr.fit(quollsdethist, 
                   mask = correctmaskclipped,
                   link = "log",
                   detectfn = "HN",
                   ncores = 4,
                   trace = TRUE,
                   model = list(D ~1, sigma ~1, g0 ~ 1),
                   details = list(fastproximity = FALSE))

#and then full models. 

fullHRquoll <- secr.fit(quollsdethist, 
                        mask = correctmaskclipped,
                        link = "log",
                        detectfn = "HR",
                        ncores = 4,
                        trace = TRUE,
                        model = list(D ~scaledurb+scaledcamp, sigma ~1, g0 ~ TrailType+bk),
                        details = list(fastproximity = FALSE))

fullEXquoll <- secr.fit(quollsdethist, 
                        mask = correctmaskclipped,
                        link = "log",
                        detectfn = "EX",
                        ncores = 4,
                        trace = TRUE,
                        model = list(D ~scaledurb+scaledcamp, sigma ~1, g0 ~ TrailType+bk),
                        details = list(fastproximity = FALSE))

fullHNquoll <- secr.fit(quollsdethist, 
                        mask = correctmaskclipped,
                        link = "log",
                        detectfn = "HN",
                        ncores = 4,
                        trace = TRUE,
                        model = list(D ~scaledurb+scaledcamp, sigma ~1, g0 ~ TrailType+bk),
                        details = list(fastproximity = FALSE))


quollmodelsnull <- secrlist(nullEXquoll,nullHRquoll, nullHNquoll)
quollmodelsfull <- secrlist(fullEXquoll,fullHRquoll,fullHNquoll)
AIC(quollmodelsnull)
AIC(quollmodelsfull)


initialsigmaquoll <- RPSV(quollsdethist,CC=TRUE)
esaPlot(quollmodelsnull,max.buffer = 8*initialsigmaquoll)
esaPlot(quollmodelsfull,max.buffer = 16*initialsigmaquoll)

fullEXquoll
#HR doesn't stabilise. So lets use EX as best function, as it at least stabilises. 

#repeat for dingoes

nullHRdingo <- secr.fit(dingoesdethist, 
                        mask = correctmaskclipped,
                        link = "log",
                        detectfn = "HR",
                        ncores = 4,
                        trace = TRUE,
                        model = list(D ~1, sigma ~1, g0 ~ 1),
                        details = list(fastproximity = FALSE))

nullEXdingo <- secr.fit(dingoesdethist, 
                        mask = correctmaskclipped,
                        link = "log",
                        detectfn = "EX",
                        ncores = 4,
                        trace = TRUE,
                        model = list(D ~1, sigma ~1, g0 ~ 1),
                        details = list(fastproximity = FALSE))

nullHNdingo <- secr.fit(dingoesdethist, 
                        mask = correctmaskclipped,
                        link = "log",
                        detectfn = "HN",
                        ncores = 4,
                        trace = TRUE,
                        model = list(D ~1, sigma ~1, g0 ~ 1),
                        details = list(fastproximity = FALSE))

#followed by the full models 

fullHRdingo <- secr.fit(dingoesdethist, 
                          mask = correctmaskclipped,
                          link = "log",
                          detectfn = "HR",
                          ncores = 4,
                          trace = TRUE,
                          model = list(D ~scaledurb+scaledcamp, sigma ~1, g0 ~ TrailType+bk),
                          details = list(fastproximity = FALSE))

fullEXdingo <- secr.fit(dingoesdethist, 
                          mask = correctmaskclipped,
                          link = "log",
                          detectfn = "EX",
                          ncores = 4,
                          trace = TRUE,
                          model = list(D ~scaledurb+scaledcamp, sigma ~1, g0 ~ TrailType+bk),
                          details = list(fastproximity = FALSE))

fullHNdingo <- secr.fit(dingoesdethist, 
                          mask = correctmaskclipped,
                          link = "log",
                          detectfn = "HN",
                          ncores = 4,
                          trace = TRUE,
                          model = list(D ~scaledurb+scaledcamp, sigma ~1, g0 ~ TrailType+bk),
                          details = list(fastproximity = FALSE))



#check best
dingomodelsnull <- secrlist(nullEXdingo,nullHRdingo, nullHNdingo)
dingomodelsfull <- secrlist(fullHNdingo,fullEXdingo,fullHRdingo)
AIC(dingomodelsnull)
AIC(dingomodelsfull)


initialsigmadingo <- RPSV(dingoesdethist,CC=TRUE)
esaPlot(dingomodelsnull,max.buffer = 8*initialsigmadingo)
esaPlot(dingomodelsfull)
esaPlot(fullHNdingo)

#check best 
fullHNdingo

#all density estimates exactly the same-v good sign!






#Get density surface
DensitySurfacedingoes <- predictDsurface(fullHNdingo, cl.D = TRUE, alpha = 0.05 )
DensitySurfacequolls <- predictDsurface(fullEXquoll, cl.D = TRUE, alpha = 0.05 )

#Extract density estimates by summing predicted values from across mask 

#density estimates
DensityEstdingoes <- c(mean(covariates(DensitySurfacedingoes)$D.0*100))
LCIsdingoes <- c(mean(covariates(DensitySurfacedingoes)$lcl.0*100))
UCIsdingoes <- c(mean(covariates(DensitySurfacedingoes)$ucl.0*100))

DensityEstquolls <- c(mean(covariates(DensitySurfacequolls)$D.0*100))
LCIsquolls <- c(mean(covariates(DensitySurfacequolls)$lcl.0*100))
UCIsquolls <- c(mean(covariates(DensitySurfacequolls)$ucl.0*100))


Densities <- data.frame(
  Density = c(DensityEstdingoes,DensityEstquolls),
  LCI = c(LCIsdingoes,LCIsquolls),
  UCI = c(UCIsdingoes,UCIsquolls),
  Species = c("Dingo","Quoll"))

write.csv(Densities, "Derived Data/DensityEstimatess.csv", row.names = FALSE)




#End. 

