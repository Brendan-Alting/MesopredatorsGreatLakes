####4. Now getting spatial detection history and camera trap operation files. 

library(camtrapR)
library(secr)
library(sf)
library(dplyr)
library(zoo)

#read in both files with camera information, and quoll independent events
camop<- read.csv(file = "Raw Data/CamOpMatrixComplete.csv", header = T)
trapsall <- read.csv(file = "Raw Data/utmpredsall.csv", header = T)
quoll15min <- read.csv(file = "Raw Data/quoll15min.csv", header = T)  
dingo15min <- read.csv(file = "Raw Data/dingo15minnopupnoun.csv",header=T)


#remove unidentifiable
dingo15min <- dingo15min %>%
  filter(!Individual=="Unidentifiable",
          !Trap=="PS4",
           !Trap=="PS7",
           !Trap=="PS8",
           !Trap=="PS10",
           !Trap=="PS22")
         


#make camera operation matrix

camopsecr <- cameraOperation(CTtable = camop,
                                stationCol = "Station",
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

quollsdethist <- spatialDetectionHistory(recordTableIndividual = quoll15min,
                                         camOp = camopsecr,
                                         CTtable = trapsall,
                                         output = "binary",
                                         stationCol = "Trap",
                                         Xcol = "x", 
                                         Ycol = "y",
                                         species = "Quoll",
                                         stationCovariateCols = c("CamType"),
                                         individualCol = "Individual",
                                         timeZone = "Australia/Sydney",
                                         recordDateTimeCol = "DateTime", 
                                         recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                         occasionLength = 1,
                                         maxNumberDays = 90,
                                         day1 = "survey",
                                         includeEffort = TRUE)

dingoesdethist <- spatialDetectionHistory(recordTableIndividual = dingo15min,
                                         camOp = camopsecr,
                                         CTtable = trapsall,
                                         output = "binary",
                                         stationCol = "Trap",
                                         Xcol = "x", 
                                         Ycol = "y",
                                         species = "Dingo",
                                         stationCovariateCols = c("CamType"),
                                         individualCol = "Individual",
                                         timeZone = "Australia/Sydney",
                                         recordDateTimeCol = "DateTime", 
                                         recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                         occasionLength = 1,
                                         maxNumberDays = 90,
                                         day1 = "survey",
                                         includeEffort = TRUE)

#Done. Now should have all elements for a secr model. 