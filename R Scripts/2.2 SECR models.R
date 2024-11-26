#Run SECR models for quolls

library(secr)

#First check detection function 

nullHRquoll <- secr.fit(quollsdethist, 
                   mask = maskclippedforpaperquoll,
                   link = "log",
                   detectfn = "HR",
                   ncores = 4,
                   trace = TRUE,
                   model = list(D ~1, sigma ~1, g0 ~ 1),
                   details = list(fastproximity = FALSE))

nullEXquoll <- secr.fit(quollsdethist, 
                   mask = maskclippedforpaperquoll,
                   link = "log",
                   detectfn = "EX",
                   ncores = 4,
                   trace = TRUE,
                   model = list(D ~1, sigma ~1, g0 ~ 1),
                   details = list(fastproximity = FALSE))

nullHNquoll <- secr.fit(quollsdethist, 
                   mask = maskclippedforpaperquoll,
                   link = "log",
                   detectfn = "HN",
                   ncores = 4,
                   trace = TRUE,
                   model = list(D ~1, sigma ~1, g0 ~ 1),
                   details = list(fastproximity = FALSE))

#check best
AIC(nullHRquoll,nullEXquoll,nullHNquoll)

#HR wins, use for full model


#repeat for dingoes

nullHRdingo <- secr.fit(dingoesdethist, 
                        mask = maskclippedforpaperdingo,
                        link = "log",
                        detectfn = "HR",
                        ncores = 4,
                        trace = TRUE,
                        model = list(D ~1, sigma ~1, g0 ~ 1),
                        details = list(fastproximity = FALSE))

nullEXdingo <- secr.fit(dingoesdethist, 
                        mask = maskclippedforpaperdingo,
                        link = "log",
                        detectfn = "EX",
                        ncores = 4,
                        trace = TRUE,
                        model = list(D ~1, sigma ~1, g0 ~ 1),
                        details = list(fastproximity = FALSE))

nullHNdingo <- secr.fit(dingoesdethist, 
                        mask = maskclippedforpaperdingo,
                        link = "log",
                        detectfn = "HN",
                        ncores = 4,
                        trace = TRUE,
                        model = list(D ~1, sigma ~1, g0 ~ 1),
                        details = list(fastproximity = FALSE))

#check best
AIC(nullHRdingo,nullEXdingo,nullHNdingo)

#EX wins, use for full model


quollSECR <- secr.fit(quollsdethist, 
                      mask = maskclippedforpaperquoll,
                      link = "log",
                      detectfn = "HR",
                      ncores = 4,
                      trace = TRUE,
                      model = list(D ~1, sigma ~1, g0 ~ CamType),
                      details = list(fastproximity = FALSE))

dingoSECR <- secr.fit(dingoesdethist, 
                      mask = maskclippedforpaperdingo,
                      link = "log",
                      detectfn = "EX",
                      ncores = 4,
                      trace = TRUE,
                      model = list(D ~1, sigma ~1, g0 ~ CamType),
                      details = list(fastproximity = FALSE))

options(scipen=999)
quollSECR
dingoSECR
#End. 

#Use these estimates to compare to the SUN models. 