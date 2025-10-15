####2. Now getting spatial detection history and camera trap operation files. This will be used for secr models. 

library(camtrapR)
library(secr)
library(sf)
library(dplyr)
library(zoo)

#read in both files with camera information, and quoll independent events
camop<- read.csv(file = "Raw Data/CamOpMatrixCompleteRevise.csv", header = T)
trapsall <- read.csv(file = "Raw Data/utmpredsall.csv", header = T)

#study area clipped is from 1. habitat mask
traps_sf <- st_as_sf(trapsall, coords = c("x", "y"), crs = st_crs(study_area_clipped))
traps_with_habitat <- st_join(traps_sf, study_area_clipped)
trapsall$Habitat <- traps_with_habitat$SimpleVeg


trapsall$TrailType <- camop$TrailType



#####remove unidentifiable from the dataframes made in 0.00cleaning, as that was total counts. 
dingo30min <- dingorawretry%>%
  filter(!Individual == "Unidentifiable")

quoll30min <- quollrawretry%>%
  filter(!Individual == "Unidentifiable")


#make camera operation matrix

camopsecr <- cameraOperation(CTtable = camop,
                                stationCol = "Trap",
                                cameraCol = "Camera",
                                setupCol = "Setup_date",
                                retrievalCol = "Retrieval_date",
                                hasProblems = TRUE,
                                camerasIndependent = FALSE,
                                byCamera = FALSE,
                                allCamsOn = FALSE,
                                writecsv = FALSE,
                                occasionStartTime = 0,
                                dateFormat = "dmy")

camopsecr<- camopsecr[, -c(1)]
camopsecr <- camopsecr[, -c(91:98)]


#First need to sort files into  correct order: 
row_names <- rownames(camopsecr)
# Extract the prefix and numeric part from row names
prefix <- gsub("[0-9]", "", row_names)
numeric_part <- as.numeric(gsub("\\D", "", row_names))

# Create a sorting index based on prefix and numeric part
sort_index <- order(prefix, numeric_part)

# Reorder the dataframe based on the sorting index
camopsecr <- camopsecr[sort_index, , drop = FALSE]

#Now can finally make spatial detection history for dingoes and quolls

quollsdethist <- spatialDetectionHistory(recordTableIndividual = quoll30min,
                                         camOp = camopsecr,
                                         CTtable = trapsall,
                                         output = "binary",
                                         stationCol = "Trap",
                                         Xcol = "x", 
                                         Ycol = "y",
                                         species = "Quoll",
                                         stationCovariateCols = c("TrailType","Habitat"),
                                         individualCol = "Individual",
                                         timeZone = "Australia/Sydney",
                                         recordDateTimeCol = "DateTime", 
                                         recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                         occasionLength = 1,
                                         maxNumberDays = 90,
                                         day1 = "survey",
                                         includeEffort = TRUE)

dingoesdethist <- spatialDetectionHistory(recordTableIndividual = dingo30min,
                                         camOp = camopsecr,
                                         CTtable = trapsall,
                                         output = "binary",
                                         stationCol = "Trap",
                                         Xcol = "x", 
                                         Ycol = "y",
                                         species = "Dingo",
                                         stationCovariateCols = c("TrailType","Habitat"),
                                         individualCol = "Individual",
                                         timeZone = "Australia/Sydney",
                                         recordDateTimeCol = "DateTime", 
                                         recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                         occasionLength = 1,
                                         maxNumberDays = 90,
                                         day1 = "survey",
                                         includeEffort = TRUE)

#Done. Now should have all elements for a secr model. 

