#Creating habitat masks
library(secr)
library(sf)
library(sp)
library(dplyr)
library(alphahull)



study_area_clipped <- sf::read_sf("Raw Data/finalhabitatmapaug2023cut.shp", crs = 32756)

#simplify habitat categories

study_area_clipped<- study_area_clipped%>%
  mutate(SimpleVeg = case_when(
    vegForm %in% c("Large Campground", "Campground", "Golf Course", "Cleared", "Grasslands", "Sand","Caravan Park","Urban")~ "Cleared_Or_Open",
    vegForm %in% c("Forested Wetlands","Saline Wetlands", "Freshwater Wetlands","Heathlands")~ "Wetland_Or_Heathland",
    vegForm %in% c("Wet Sclerophyll Forests (Grassy sub-formation)", "Wet Sclerophyll Forests (Shrubby sub-formation)","Dry Sclerophyll Forests (Shrub/grass sub-formation)","Dry Sclerophyll Forests (Shrubby sub-formation)", "Dry Sclerophyll Sandmined (Shrubby sub-formation)","Rainforests")~"Forest", TRUE~NA))



trapsall <- read.csv(file = "Raw Data/utmpredsall.csv")


trapsnocovs <- trapsall%>%
  select(-CamType)

traps <- read.traps(data = trapsnocovs, detector = "proximity")







#Make excessively large mask, for covariates. Use trapsfile we made, and then the study area clipped. 
#this is going to take a while. It's a huge area w small resolution 

correctmask <-make.mask(traps, 
              type = 'trapbuffer', 
              poly = study_area_clipped, 
              buffer = 16000, 
              keep.poly = F, 
              nx = 121, 
              ny = 121)




###This is to cut the points in the lakes, and some on the other side of the lake that shouldn't be in there (not connected because across waterbody) 
correctmaskclipped <- correctmask%>%
  filter(!(x >= 451741.26 & x <= 455162.12 &
             y >= 6416000.44 & y <= 6417229.43)&
           !(x < 439765.29 & y > 6413338.22))


#Now we add distance to urban covariates, and habitat. 


###Use our new coords
HPD <- read.csv(file = "Raw Data/CoordsCovsARS.csv", header = TRUE) %>%  
  st_as_sf(coords = c("x", "y"), crs = 32756) 

HPSUrban <- HPD %>%
  filter(Habitat == "Urban")

HPSCampgroud <- HPD%>%
  filter(!Habitat =="Urban",
         !Name == "Jimmys Beach ",
         !Name == "Reflections Hawks Nest",
         !Name == "Pacific Palms Holiday Park",
         !Name == "The Ruins Campground")


HPSCampgroud <- HPSCampgroud%>%
  mutate(x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2])%>%
  as.data.frame()%>%
  dplyr::select(x,y)

HPSUrban <- HPSUrban%>%
  mutate(x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2])%>%
  as.data.frame()%>%
  dplyr::select(x,y)


distsurbcorrectdingoes<- nedist(HPSUrban, correctmaskclipped,correctmaskclipped)

distscampgroundcorrectdingoes <- nedist(HPSCampgroud, correctmaskclipped,correctmaskclipped)


distsurbcorrectdingoes <- (data.frame(apply(distsurbcorrectdingoes, 2, min, na.rm = TRUE)))

distscampgroundcorrectdingoes <- (data.frame(apply(distscampgroundcorrectdingoes, 2, min, na.rm = TRUE)))


distsurbcorrectdingoes <- distsurbcorrectdingoes%>%
  rename(urbdist = apply.distsurbcorrectdingoes..2..min..na.rm...TRUE.)


distscampgroundcorrectdingoes <- distscampgroundcorrectdingoes%>%
  rename(campgrounddist = apply.distscampgroundcorrectdingoes..2..min..na.rm...TRUE.)


distsurbcorrectdingoes <- as.data.frame(distsurbcorrectdingoes)
distscampgroundcorrectdingoes <- as.data.frame(distscampgroundcorrectdingoes)


which(is.infinite(distsurbcorrectdingoes$urbdist))
which(is.infinite(distscampgroundcorrectdingoes$campgrounddist))


#some infinite values (couldn't calculate path). Look visually, and take the point nearest plus a little bit. 
distscampgroundcorrectdingoes <- distscampgroundcorrectdingoes %>%
  mutate(
    campgrounddist= case_when(
      row_number() == 6 ~ 4710,
      TRUE ~ campgrounddist 
    )
  )
##

distsurbcorrectdingoes <- distsurbcorrectdingoes %>%
  mutate(
    urbdist= case_when(
      row_number() == 6 ~ 1192,
      TRUE ~ urbdist  
    )
  )
##


maskdummy <- st_as_sf(correctmaskclipped, coords = c('x','y'),crs = st_crs(study_area_clipped))



maskdummy$urbdists <- distsurbcorrectdingoes$urbdist
maskdummy$campgrounddists <- distscampgroundcorrectdingoes$campgrounddist

maskdummy$scaledurb <- scale(maskdummy$urbdists)
maskdummy$scaledcamp <- scale(maskdummy$campgrounddists)

maskdummy$point_id <- 1:nrow(maskdummy)

covariates(correctmaskclipped) <- NULL

correctmaskclipped <- addCovariates(correctmaskclipped,maskdummy)


#plot to check .
plot(correctmaskclipped, covariate = 'scaledurb',legend = FALSE)
plot(correctmaskclipped, covariate = 'scaledcamp',legend = FALSE)




####DONE 

