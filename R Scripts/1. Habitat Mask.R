#Creating habitat mask
library(secr)
library(sf)
library(sp)
library(dplyr)
library(alphahull)




study_area_clipped <- sf::read_sf("Raw Data/studyareaclipped.shp")



trapsall <- read.csv(file = "Raw Data/utmpredsall.csv")
trapsnocovs <- trapsall%>%
  select(-CamType)

traps <- read.traps(data = trapsnocovs, detector = "proximity")







#Make excessively large mask, for covariates. Use trapsfile we made, and then the study area clipped. 
#this is going to take a while. It's a huge area w small resolution 

maskclippedforpaperdingo <- make.mask(traps, 
                                 type = 'trapbuffer', 
                                 poly = study_area_clipped, 
                                 buffer = 8000, 
                                 keep.poly = F, 
                                 nx = 64, 
                                 ny = 64)

maskclippedforpaperquoll<- make.mask(traps, 
                                     type = 'trapbuffer', 
                                     poly = study_area_clipped, 
                                     buffer = 4000, 
                                     keep.poly = F, 
                                     nx = 64, 
                                     ny = 64)


#mask for each species complete. 
