####Running temporal analysis, with traps assigned to correct locations of risk. 

library(overlap)
library(circular)
library(dplyr)
library(lubridate)

#First just get a bit of cleaning done. 
quolltemporal <- read.csv(file = "Raw Data/quoll15min.csv",header =T)
dingotemporal <- read.csv(file = "Raw Data/dingo15minnopup.csv",header =T)
foxtemporal <- read.csv(file = "Raw Data/fox15min.csv",header=T)
goannatemporal <- read.csv(file = "Raw Data/goanna15min.csv",header=T)



quoll <- quolltemporal %>%
  group_by(Trap,DateTime)%>%
  summarise(detections = n())

quoll$DateTime <- as.POSIXct(quoll$DateTime)

dingo <- dingotemporal %>%
  group_by(Trap,DateTime)%>%
  summarise(detections = n())

dingo$DateTime <- as.POSIXct(dingo$DateTime)

fox <- foxtemporal %>%
  group_by(Trap,DateTime)%>%
  summarise(detections = n())

fox$DateTime <- as.POSIXct(fox$DateTime)

goanna <- goannatemporal %>%
  group_by(Trap,DateTime)%>%
  summarise(detections = n())

goanna$DateTime <-as.POSIXct(goanna$DateTime)

quoll$Species <- "Quoll"
dingo$Species <- "Dingo"
fox$Species <- "Fox"
goanna$Species <- "Goanna"

speciesdetections <- rbind(quoll,dingo,fox,goanna)



##Join with traps file so can separate into high and low risk  

#read in trapsall 
trapstemporal <- read.csv(file = "Derived Data/riskratingtraps.csv")


speciesdetections <- speciesdetections %>%
  left_join(select(trapstemporal, Trap, closest_risk), by = "Trap")%>%
  rename(DateTimeOriginal = DateTime)

#lets get summary stats 

summaryrisk <- speciesdetections%>%
  group_by(Species)%>%
  count(closest_risk)


#For ovelap and circular need to make into radians for easier analysis. 


hours <- hour(speciesdetections$DateTimeOriginal)
minutes <- minute(speciesdetections$DateTimeOriginal)
seconds <- second(speciesdetections$DateTimeOriginal)

radians <- (hours * 3600 + minutes * 60 + seconds) * (2 * pi / (24 * 3600))
speciesdetections$radians <- radians 

HRdetections <- speciesdetections %>% filter(risk_category == "High Risk")
LRdetections <- speciesdetections %>% filter(risk_category == "Low Risk")



####HRS SEPARATED
detHRquoll <- HRdetections[HRdetections$Species == "Quoll",]
detHRfox <- HRdetections[HRdetections$Species == "Fox",]
detHRdingo <- HRdetections[HRdetections$Species == "Dingo",]
detHRgoanna <- HRdetections[HRdetections$Species == "Goanna",]

#LR SEPARATED
detLRquoll <- LRdetections[LRdetections$Species == "Quoll",]
detLRfox <- LRdetections[LRdetections$Species == "Fox",]
detLRdingo <- LRdetections[LRdetections$Species == "Dingo",]
detLRgoanna <- LRdetections[LRdetections$Species == "Goanna",]

####Get radians

quollradsHR <- unlist(detHRquoll$radians)
foxradsHR <- unlist(detHRfox$radians)
dingoradsHR <- unlist(detHRdingo$radians)
goannaradsHR <- unlist(detHRgoanna$radians)

quollradsLR <- unlist(detLRquoll$radians)
foxradsLR <- unlist(detLRfox$radians)
dingoradsLR <- unlist(detLRdingo$radians)
goannaradsLR <- unlist(detLRgoanna$radians)



####Now we will calculate the temporal overlap stuff: 

##Quoll fox
quollfoxoverlapHR <- overlapEst(quollradsHR, foxradsHR, type = "Dhat1")
quollfoxoverlapHR

quollfoxoverlapHRboot <- bootstrap(quollradsHR, foxradsHR, 10000, type = "Dhat1")
BSQuollfoxHR <- mean(quollfoxoverlapHRboot)
BSQuollfoxHR

quollfoxoverlapLR <- overlapEst(quollradsLR, foxradsLR, type = "Dhat4")
quollfoxoverlapLR

quollfoxoverlapLRboot <- bootstrap(quollradsLR, foxradsLR, 10000, type = "Dhat4")
BSQuollfoxLR <- mean(quollfoxoverlapLRboot)
BSQuollfoxLR


confintsquollfoxHR <- as.data.frame(bootCIlogit(quollfoxoverlapHR, quollfoxoverlapHRboot))
confintsquollfoxHR$UseState <- "HR"
confintsquollfoxHR <- confintsquollfoxHR["basic0", , drop = FALSE]
confintsquollfoxHR$Overlap <- quollfoxoverlapHR

confintsquollfoxLR <- as.data.frame(bootCIlogit(quollfoxoverlapLR, quollfoxoverlapLRboot))
confintsquollfoxLR$UseState <- "LR"
confintsquollfoxLR <- confintsquollfoxLR["basic0", , drop = FALSE]
confintsquollfoxLR$Overlap <- quollfoxoverlapLR

confintsquollfox <- rbind(confintsquollfoxHR, confintsquollfoxLR)
confintsquollfox$Comparison <- "Quoll Fox"

#Quoll Dingo 
quolldingooverlapHR <- overlapEst(quollradsHR, dingoradsHR, type = "Dhat1")
quolldingooverlapHR

quolldingooverlapHRboot <- bootstrap(quollradsHR, dingoradsHR, 10000, type = "Dhat1")
BSQuolldingoHR <- mean(quolldingooverlapHRboot)
BSQuolldingoHR

quolldingooverlapLR <- overlapEst(quollradsLR, dingoradsLR, type = "Dhat4")
quolldingooverlapLR

quolldingooverlapLRboot <- bootstrap(quollradsLR, dingoradsLR, 10000, type = "Dhat4")
BSQuolldingoLR <- mean(quolldingooverlapLRboot)
BSQuolldingoLR



confintsquolldingoHR <- as.data.frame(bootCIlogit(quolldingooverlapHR, quolldingooverlapHRboot))
confintsquolldingoHR$UseState <- "HR"
confintsquolldingoHR <- confintsquolldingoHR["basic0", , drop = FALSE]
confintsquolldingoHR$Overlap <- quolldingooverlapHR

confintsquolldingoLR <- as.data.frame(bootCIlogit(quolldingooverlapLR, quolldingooverlapLRboot))
confintsquolldingoLR$UseState <- "LR"
confintsquolldingoLR <- confintsquolldingoLR["basic0", , drop = FALSE]
confintsquolldingoLR$Overlap <- quolldingooverlapLR


confintsquolldingo <- rbind(confintsquolldingoHR, confintsquolldingoLR)
confintsquolldingo$Comparison <- "Quoll Dingo"


###Fox dingo
foxoverlapHR <- overlapEst(foxradsHR, dingoradsHR, type = "Dhat1")
foxoverlapHR

foxoverlapHRboot <- bootstrap(foxradsHR, dingoradsHR, 10000, type = "Dhat1")
BSfoxHR <- mean(foxoverlapHRboot)
BSfoxHR

foxoverlapLR <- overlapEst(foxradsLR, dingoradsLR, type = "Dhat4")
foxoverlapLR

foxoverlapLRboot <- bootstrap(foxradsLR, dingoradsLR, 10000, type = "Dhat4")
BSfoxLR <- mean(foxoverlapLRboot)
BSfoxLR



confintsfoxdingoHR <- as.data.frame(bootCIlogit(foxoverlapHR, foxoverlapHRboot))
confintsfoxdingoHR$UseState <- "HR"
confintsfoxdingoHR <- confintsfoxdingoHR["basic0", , drop = FALSE]
confintsfoxdingoHR$Overlap <- foxoverlapHR

confintsfoxdingoLR <- as.data.frame(bootCIlogit(foxoverlapLR, foxoverlapLRboot))
confintsfoxdingoLR$UseState <- "LR"
confintsfoxdingoLR <- confintsfoxdingoLR["basic0", , drop = FALSE]
confintsfoxdingoLR$Overlap <- foxoverlapLR



confintsfoxdingo <- rbind(confintsfoxdingoHR, confintsfoxdingoLR)
confintsfoxdingo$Comparison <- "Fox Dingo"

#######Goanna fox

###Fox goanna
foxgoannaoverlapHR <- overlapEst(foxradsHR, goannaradsHR, type = "Dhat1")
foxgoannaoverlapHR

foxgoannaoverlapHRboot <- bootstrap(foxradsHR, goannaradsHR, 10000, type = "Dhat1")
BSfoxgoannaHR <- mean(foxgoannaoverlapHRboot)
BSfoxgoannaHR

foxgoannaoverlapLR <- overlapEst(foxradsLR, goannaradsLR, type = "Dhat4")
foxgoannaoverlapLR

foxgoannaoverlapLRboot <- bootstrap(foxradsLR, goannaradsLR, 10000, type = "Dhat4")
BSfoxgoannaLR <- mean(foxgoannaoverlapLRboot)
BSfoxgoannaLR



confintsfoxgoannaHR <- as.data.frame(bootCIlogit(foxgoannaoverlapHR, foxgoannaoverlapHRboot))
confintsfoxgoannaHR$UseState <- "HR"
confintsfoxgoannaHR <- confintsfoxgoannaHR["basic0", , drop = FALSE]
confintsfoxgoannaHR$Overlap <- foxgoannaoverlapHR

confintsfoxgoannaLR <- as.data.frame(bootCIlogit(foxgoannaoverlapLR, foxgoannaoverlapLRboot))
confintsfoxgoannaLR$UseState <- "LR"
confintsfoxgoannaLR <- confintsfoxgoannaLR["basic0", , drop = FALSE]
confintsfoxgoannaLR$Overlap <- foxgoannaoverlapLR



confintsfoxgoanna <- rbind(confintsfoxgoannaHR, confintsfoxgoannaLR)
confintsfoxgoanna$Comparison <- "Fox Monitor"

################Dingo Goanna



###dingo goanna
dingogoannaoverlapHR <- overlapEst(dingoradsHR, goannaradsHR, type = "Dhat4")
dingogoannaoverlapHR

dingogoannaoverlapHRboot <- bootstrap(dingoradsHR, goannaradsHR, 10000, type = "Dhat4")
BSdingogoannaHR <- mean(dingogoannaoverlapHRboot)
BSdingogoannaHR

dingogoannaoverlapLR <- overlapEst(dingoradsLR, goannaradsLR, type = "Dhat4")
dingogoannaoverlapLR

dingogoannaoverlapLRboot <- bootstrap(dingoradsLR, goannaradsLR, 10000, type = "Dhat4")
BSdingogoannaLR <- mean(dingogoannaoverlapLRboot)
BSdingogoannaLR



confintsdingogoannaHR <- as.data.frame(bootCIlogit(dingogoannaoverlapHR, dingogoannaoverlapHRboot))
confintsdingogoannaHR$UseState <- "HR"
confintsdingogoannaHR <- confintsdingogoannaHR["basic0", , drop = FALSE]
confintsdingogoannaHR$Overlap <- dingogoannaoverlapHR

confintsdingogoannaLR <- as.data.frame(bootCIlogit(dingogoannaoverlapLR, dingogoannaoverlapLRboot))
confintsdingogoannaLR$UseState <- "LR"
confintsdingogoannaLR <- confintsdingogoannaLR["basic0", , drop = FALSE]
confintsdingogoannaLR$Overlap <- dingogoannaoverlapLR



confintsdingogoanna <- rbind(confintsdingogoannaHR, confintsdingogoannaLR)
confintsdingogoanna$Comparison <- "Dingo Monitor"



###Quoll goanna
quollgoannaoverlapHR <- overlapEst(quollradsHR, goannaradsHR, type = "Dhat1")
quollgoannaoverlapHR

quollgoannaoverlapHRboot <- bootstrap(quollradsHR, goannaradsHR, 10000, type = "Dhat1")
BSquollgoannaHR <- mean(quollgoannaoverlapHRboot)
BSquollgoannaHR

quollgoannaoverlapLR <- overlapEst(quollradsLR, goannaradsLR, type = "Dhat4")
quollgoannaoverlapLR

quollgoannaoverlapLRboot <- bootstrap(quollradsLR, goannaradsLR, 10000, type = "Dhat4")
BSquollgoannaLR <- mean(quollgoannaoverlapLRboot)
BSquollgoannaLR



confintsquollgoannaHR <- as.data.frame(bootCIlogit(quollgoannaoverlapHR, quollgoannaoverlapHRboot))
confintsquollgoannaHR$UseState <- "HR"
confintsquollgoannaHR <- confintsquollgoannaHR["basic0", , drop = FALSE]
confintsquollgoannaHR$Overlap <- quollgoannaoverlapHR

confintsquollgoannaLR <- as.data.frame(bootCIlogit(quollgoannaoverlapLR, quollgoannaoverlapLRboot))
confintsquollgoannaLR$UseState <- "LR"
confintsquollgoannaLR <- confintsquollgoannaLR["basic0", , drop = FALSE]
confintsquollgoannaLR$Overlap <- quollgoannaoverlapLR


####Now combine all together into a dataframe, for plotting. 


confintsquollgoanna <- rbind(confintsquollgoannaHR, confintsquollgoannaLR)
confintsquollgoanna$Comparison <- "Quoll Monitor"


confintsalloverlap <- rbind(confintsfoxdingo, confintsquolldingo, confintsquollfox, confintsfoxgoanna, confintsdingogoanna,confintsquollgoanna)

confintsalloverlap$Comparison <- factor(confintsalloverlap$Comparison, levels = c("Fox Dingo", "Quoll Dingo", "Quoll Fox", "Fox Monitor", "Quoll Monitor", "Dingo Monitor"))


library(circular)
####Now watson wheeler tests of homogeneity of angles stuff. 

#quoll-fox 
QFHR <- watson.wheeler.test(list(quollradsHR, foxradsHR)) #quoll fox LR
QFHR$Combo <- "QFHR"
QFLR <- watson.wheeler.test(list(quollradsLR, foxradsLR))
QFLR$Combo <- "QFLR"
#quoll-dingo
QDHR <- watson.wheeler.test(list(quollradsHR, dingoradsHR))
QDHR$Combo <- "QDHR"
QDLR <- watson.wheeler.test(list(quollradsLR, dingoradsLR))
QDLR$Combo <- "QDLR"
##fox-dingo
FDHR <- watson.wheeler.test(list(foxradsHR, dingoradsHR))
FDHR$Combo <- "FDHR"
FDLR <- watson.wheeler.test(list(foxradsLR, dingoradsLR))
FDLR$Combo <- "FDLR"
#Fox-monitor
FMHR <- watson.wheeler.test(list(foxradsHR, goannaradsHR))
FMHR$Combo <- "FMHR"
FMLR <- watson.wheeler.test(list(foxradsLR, goannaradsLR))
FMLR$Combo <- "FMLR"
#Dingo-monitor
MDHR <- watson.wheeler.test(list(goannaradsHR, dingoradsHR))
MDHR$Combo <- "MDHR"
MDLR <- watson.wheeler.test(list(goannaradsLR, dingoradsLR))
MDLR$Combo <- "MDLR"
#Quoll-monitor
QMHR <- watson.wheeler.test(list(quollradsHR, goannaradsHR))
QMHR$Combo <- "QMHR"
QMLR <- watson.wheeler.test(list(quollradsLR, goannaradsLR))
QMLR$Combo <- "QMLR"

#Combine results
mergedresults <- data.frame(rbind(QFHR,QFLR,QDHR,QDLR,FDHR,FDLR,FMHR,FMLR,MDHR,MDLR,QMHR,QMLR))
mergedresults
#done. 

#Prop.test




###now plot trends


#7.2- Plotting temporal effects

####First we will plot the overlap: 

png("Figures/Overlap In LR and HR.png", width = 12, height = 13, res= 300, units = "in")

par(mfrow = c(2,1))
###Overlap plots for HR
densityPlot(foxradsHR, xscale = 24,main = "", add = FALSE, xcenter = c("midnight"),col = "orange", lwd = 4, extend = NULL)
densityPlot(goannaradsHR, xscale = 24, add = TRUE, xcenter = c("midnight"), col = "#0072B2", lwd =1, extend = NULL)
densityPlot(quollradsHR, xscale = 24, add = TRUE, xcenter = c("midnight"), col = "black", lwd = 4, extend = NULL)
densityPlot(dingoradsHR, xscale = 24, add = TRUE, xcenter = c("midnight"), col = "#CC79A7", lwd =4, extend = NULL)

legend("left", legend = c("Fox (n = 62)", "Quoll (n = 53)", "Dingo (n = 457)", "Lace \nMonitor (n = 310)", ""), col = c("orange", "black", "#CC79A7", "#0072B2", "white"), lwd = 3)

mtext("a) High Risk", side = 3, line = 1, adj = 0, cex = 1.5) 


###Overlap plots for LReral
densityPlot(foxradsLR, xscale = 24, add = FALSE, xcenter = c("midnight"), col = "orange",main = "", lwd = 4, extend = NULL)
densityPlot(quollradsLR, xscale = 24, add = TRUE, xcenter = c("midnight"), col = "black", lwd = 4, extend = NULL)
densityPlot(dingoradsLR, xscale = 24, add = TRUE, xcenter = c("midnight"), col = "#CC79A7", lwd =4, extend = NULL)
densityPlot(goannaradsLR, xscale = 24, add = TRUE, xcenter = c("midnight"), col = "#0072B2", lwd =1, extend = NULL)

mtext("b) Low Risk", side = 3, line = 1, adj = 0, cex = 1.5)  # Label for the LR plot

legend("left", legend = c("Fox (n = 182)", "Quoll (n = 125)", "Dingo (n = 143)", "Lace\nMonitor (n = 401)", ""), col = c("orange", "black", "#CC79A7", "#0072B2", "white"), lwd = 3)

dev.off()


#Finished plotting that. 

#Next step is plot the overlap ones 

############plotting all overlaps together############ 
overlapallplot <- ggplot(confintsalloverlap, aes(x = Comparison, y = Overlap, color = UseState)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), size = 1.2, 
                position = position_dodge(width = 0.5), lwd = 1, width = 0.2) +
  scale_color_manual(values = c("HR" = "orange", "LR" = "#0072B2"),
                     labels = c("HR" = "High Risk", "LR" = "Low Risk")) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.5, 
                  position = position_dodge(width = 0.5)) +
  labs(x = "Species Dyad", y = "Temporal Overlap Coefficient", color = "Dingo Risk \n     Zone") +
  ylim(0,1) +
  theme_minimal()+
  theme(text = element_text(size = 18),  # Adjust this value for the desired text size
        axis.title.x = element_text(size = 25, margin = margin(t=10)),  # Adjust this value for the desired x-axis label size
        axis.title.y = element_text(size = 25),  # Adjust this value for the desired y-axis label size
        legend.title = element_text(size = 25),  # Adjust this value for the desired legend title size
        legend.text = element_text(size = 20),
        panel.grid  = element_blank(),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.line.x = element_line(color = "black", size = 0.5), # Add x-axis line
        axis.line.y = element_line(color = "black", size = 0.5))

overlapallplot

png("Figures/Overlap Plots.png", width = 12, height = 10, res= 600, units = "in")

overlapallplot 

dev.off()






