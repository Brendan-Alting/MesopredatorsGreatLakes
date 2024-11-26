library(spatstat)
library(MCMCvis)
library(sf)
library(tidyverse)
library(SpatialKDE)

#read border sf
borderhuh <- st_read("Raw Data/studyareaclipped.shp")
bordersimply <- st_buffer(borderhuh, 10,crs =28356)
st_crs(bordersimply)<-(28356)
bordersimply$geometry <- bordersimply$geometry / 1000

owinborder <- as.owin(bordersimply)



# camera locations
cam.locs<- read.csv("Raw Data/camlocs_MLNP.csv",header=T)


cam.locs$Easting<- cam.locs$Easting/1000 # km
cam.locs$Northing<- cam.locs$Northing/1000 # km

#extract cam type
cam.type <- cam.locs$cam.type

cam.region<- ripras(cam.locs$Easting,cam.locs$Northing) #convex polygon around locations

cam.region<- dilation(cam.region, 0)  # add 0km buffer
cam.region <- intersect.owin(cam.region, owinborder)


#---------------------
# Fox data
#

nsim<- 100
r=seq(0,2,0.1)


DingoFox.res<- matrix(NA, nrow=length(r),ncol=nsim)
DingoQuoll.res<- matrix(NA, nrow=length(r),ncol=nsim)
DingoGoanna.res<- matrix(NA, nrow=length(r),ncol=nsim)
FoxQuoll.res<- matrix(NA, nrow=length(r),ncol=nsim)
FoxGoanna.res <- matrix(NA, nrow=length(r),ncol=nsim)
QuollGoanna.res <-matrix(NA, nrow=length(r),ncol=nsim)




foxesfit<- readRDS("Derived Data/DiscreteSpaceFox.rds") # Use SPA_nimble.r to generate posteriors
dingoesfitnopup<- readRDS("Derived Data/DiscreteSpaceDingonopupnountest.rds") # use SPA_nimble.r to generate posteriors
quollsfit<- readRDS("Derived Data/DiscreteSpaceQuoll.rds") # use SPA_nimble.r to generate posteriors
goannasfit <- readRDS("Derived Data/DiscreteSpaceGoanna.rds") # use SPA_nimble.r to generate posteriors


##First define a function for jittering - the data was in a grid cell, so we can imagine we are randomising into an  grid cell,with centre x and y.Add or substract a random value from the middle, with length of the 1/2 side length of squre. Then the points will be equally random in each grid cell, filling the study area.


jitter_points <- function(Sx, Sy, grid_side_length) {
  # Generate random offsets within the grid cell
  offset_x <- runif(length(Sx), min = -grid_side_length/2, max = grid_side_length/2)
  offset_y <- runif(length(Sy), min = -grid_side_length/2, max = grid_side_length/2)
  
  # Apply the offsets to the original coordinates
  jittered_x <- Sx + offset_x
  jittered_y <- Sy + offset_y
  
  return(list(x = jittered_x, y = jittered_y))
}


#1.1 (1.1km for goanna. )
#0.65 (650m side length for dingo fox quoll)


jittered_data <- list()

for (i in 1:nsim) {
  cat("Doing simulation ", i, "\n")
  
  # Dingo data
  out <- MCMCpstr(dingoesfitnopup, c("g", "z"), type = "chains")
  nsamp <- sample(0000:30000, size = 50, replace = FALSE)
  
  Sx <- out$g[, 1, nsamp]
  Sy <- out$g[, 2, nsamp]
  z <- out$z[, nsamp]
  
  dingo.x <- Sx[z == 1]
  dingo.y <- Sy[z == 1]
  
  ok <- inside.owin(dingo.x, dingo.y, cam.region)
  dingo.x <- dingo.x[ok]
  dingo.y <- dingo.y[ok]
  
  dingo_jittered <- jitter_points(dingo.x, dingo.y, grid_side_length = 0.65)
  dingo.x <- dingo_jittered$x
  dingo.y <- dingo_jittered$y
  cat("Applied jitter to Dingo data", i, "\n")
  
  # Fox data
  out <- MCMCpstr(foxesfit, c("g", "z"), type = "chains")
  
  Sx <- out$g[, 1, nsamp]
  Sy <- out$g[, 2, nsamp]
  z <- out$z[, nsamp]
  
  foxes.x <- Sx[z == 1]
  foxes.y <- Sy[z == 1]
  
  ok <- inside.owin(foxes.x, foxes.y, cam.region)
  foxes.x <- foxes.x[ok]
  foxes.y <- foxes.y[ok]
  
  foxes_jittered <- jitter_points(foxes.x, foxes.y, grid_side_length = 0.65)
  foxes.x <- foxes_jittered$x
  foxes.y <- foxes_jittered$y
  cat("Applied jitter to Fox data", i, "\n")
  
  # Quoll data
  out <- MCMCpstr(quollsfit, c("g", "z"), type = "chains")
  
  Sx <- out$g[, 1, nsamp]
  Sy <- out$g[, 2, nsamp]
  z <- out$z[, nsamp]
  
  quoll.x <- Sx[z == 1]
  quoll.y <- Sy[z == 1]
  
  ok <- inside.owin(quoll.x, quoll.y, cam.region)
  quoll.x <- quoll.x[ok]
  quoll.y <- quoll.y[ok]
  
  quoll_jittered <- jitter_points(quoll.x, quoll.y, grid_side_length = 0.65)
  quoll.x <- quoll_jittered$x
  quoll.y <- quoll_jittered$y
  cat("Applied jitter to Quoll data", i, "\n")
  
  # Goanna data
  out <- MCMCpstr(goannasfit, c("g", "z"), type = "chains")
  
  Sx <- out$g[, 1, nsamp]
  Sy <- out$g[, 2, nsamp]
  z <- out$z[, nsamp]
  
  goanna.x <- Sx[z == 1]
  goanna.y <- Sy[z == 1]
  
  ok <- inside.owin(goanna.x, goanna.y, cam.region)
  goanna.x <- goanna.x[ok]
  goanna.y <- goanna.y[ok]
  
  goanna_jittered <- jitter_points(goanna.x, goanna.y, grid_side_length = 1.1)
  goanna.x <- goanna_jittered$x
  goanna.y <- goanna_jittered$y
  cat("Applied jitter to Goanna data", i, "\n")
  
  # Combine all data
  X <- c(dingo.x, foxes.x, goanna.x, quoll.x)
  Y <- c(dingo.y, foxes.y, goanna.y, quoll.y)
  M <- factor(c(rep("D", length(dingo.x)), rep("F", length(foxes.x)), rep("G", length(goanna.x)), rep("Q", length(quoll.x))))
  
  # Save the results (optional)
  # List of jittered points for later use if needed
  jittered_data[[i]] <- list(X = X, Y = Y, M = M)
}

saveRDS(jittered_data, "Derived Data/jitteredpoints50.rds")#also have jitteredpoints.rds, which is 100. #also have jitteredpoints10.rds, which is 10. 

# Mark Connect Process
for (i in 1:nsim) {
  cat("Doing mark connect for simulation ", i, "\n")
  
  # Use previously jittered points (assuming they were saved in a list)
  jittered_points <- jittered_data[[i]]
  X <- jittered_points$X
  Y <- jittered_points$Y
  M <- jittered_points$M
  
  hr.ppp <- ppp(X, Y, marks = M, window = cam.region)
  hr.ppp <- as.ppp(hr.ppp)
  
  # Dingo and Fox mark connect
  DingoFox <- markconnect(hr.ppp, "D", "F", normalise = TRUE, r = r, correction = "translate")
  cat("DingoFoxMarkConnectDone", i, "\n")
  
  # Dingo and Quoll mark connect
  DingoQuoll <- markconnect(hr.ppp, "D", "Q", normalise = TRUE, r = r, correction = "translate")
  cat("DingoQuollMarkConnectDone", i, "\n")
  
  # Dingo and Goanna mark connect
  DingoGoanna <- markconnect(hr.ppp, "D", "G", normalise = TRUE, r = r, correction = "translate")
  cat("DingoGoannaMarkConnectDone", i, "\n")
  
  FoxQuoll <- markconnect(hr.ppp, "F", "Q", normalise=TRUE, r=r,correction = "translate")
  cat("FoxQuollMarkConnectDone",i,"\n")
  
  FoxGoanna <- markconnect(hr.ppp, "F", "G", normalise=TRUE, r=r,correction = "translate")
  cat("FoxGoannaMarkConnectDone",i,"\n")
  
  QuollGoanna <- markconnect(hr.ppp, "Q", "G", normalise=TRUE, r=r,correction = "translate")
  cat("AllMarkConnectDone",i,"\n")
  
  DingoFox.res[,i] <- DingoFox$trans
  DingoQuoll.res[,i] <- DingoQuoll$trans
  DingoGoanna.res[,i] <- DingoGoanna$trans
  FoxQuoll.res[,i] <- FoxQuoll$trans
  FoxGoanna.res[,i] <- FoxGoanna$trans
  QuollGoanna.res[,i] <- QuollGoanna$trans
}

###putting this ousside for a sec


DingoFox<- apply(DingoFox.res,1,mean)
DingoFox.cl<- apply(DingoFox.res,1,quantile, c(0.025,0.975))

DingoQuoll<- apply(DingoQuoll.res,1,mean)
DingoQuoll.cl<- apply(DingoQuoll.res,1,quantile, c(0.025,0.975))

DingoGoanna<- apply(DingoGoanna.res,1,mean)
DingoGoanna.cl<- apply(DingoGoanna.res,1,quantile, c(0.025,0.975))

FoxQuoll <- apply(FoxQuoll.res,1,mean)
FoxQuoll.cl <- apply(FoxQuoll.res,1,quantile,c(0.025,0.975))

FoxGoanna <- apply(FoxGoanna.res,1,mean)
FoxGoanna.cl <- apply(FoxGoanna.res,1,quantile,c(0.025,0.975))

QuollGoanna<- apply(QuollGoanna.res,1,mean)
QuollGoanna.cl <- apply(QuollGoanna.res,1,quantile,c(0.025,0.975))


DingoFox.results<- data.frame(r=r,pr=DingoFox,lcl=DingoFox.cl[1,],ucl=DingoFox.cl[2,],Dyad="Dingo-Fox")

DingoQuoll.results<- data.frame(r=r,pr=DingoQuoll,lcl=DingoQuoll.cl[1,],ucl=DingoQuoll.cl[2,],Dyad="Dingo-Quoll")

DingoGoanna.results<- data.frame(r=r,pr=DingoGoanna,lcl=DingoGoanna.cl[1,],ucl=DingoGoanna.cl[2,],Dyad="Dingo-Monitor")

FoxQuoll.results<- data.frame(r=r,pr=FoxQuoll,lcl=FoxQuoll.cl[1,],ucl=FoxQuoll.cl[2,],Dyad="Fox-Quoll")

FoxGoanna.results<- data.frame(r=r,pr=FoxGoanna,lcl=FoxGoanna.cl[1,],ucl=FoxGoanna.cl[2,],Dyad="Fox-Monitor")

QuollGoanna.results<- data.frame(r=r,pr=QuollGoanna,lcl=QuollGoanna.cl[1,],ucl=QuollGoanna.cl[2,],Dyad="Quoll-Monitor")

####now we've got the six results. Next step is to write these DIRECTLY into the MC.Results df. Then we delete all of this and pretend that we did it all in this. 

#then delete this. 




























################################################################
########################FINISH DELETE ABOVE#####################
################################################################















mc.results<- bind_rows(DingoFox.results,DingoQuoll.results,DingoGoanna.results,FoxQuoll.results,FoxGoanna.results,QuollGoanna.results)

mc.results$Dyad <- factor(mc.results$Dyad, levels = c("Fox-Quoll","Quoll-Monitor","Fox-Monitor", "Dingo-Fox","Dingo-Monitor","Dingo-Quoll"))

mc.results <- mc.results%>%
  mutate(Type = case_when(Dyad == "Fox-Quoll"~"Meso-Meso",
                          Dyad == "Quoll-Monitor"~"Meso-Meso",
                          Dyad == "Fox-Monitor"~"Meso-Meso",
                          Dyad == "Dingo-Fox"~"Apex-Meso",
                          Dyad == "Dingo-Monitor"~"Apex-Meso",
                          Dyad == "Dingo-Quoll"~"Apex-Meso"))

markcorrplot <- mc.results %>%
  ggplot(aes(r, pr, colour=Dyad)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=lcl, ymax=ucl, colour=Dyad), alpha=0.25, fill="grey40", linetype="blank") +
  ylab(expression(paste(italic(p[i][j](r))))) +
  xlab(expression(paste(italic(r),~ (km)))) +
  geom_hline(aes(yintercept = 1), linetype=2, size=1) +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")) +  # Distinct colors
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")) +
  labs(color= "Species Dyad")+
  theme_bw() +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size =12),
    axis.text = element_text(size =14),
    axis.title = element_text(size = 14)
    )

markcorrplot



png("Figures/MarkCorr100.png", width = 8, height = 5, res= 300, units = "in")

markcorrplot
dev.off()


write.csv(mc.results, file = "Derived Data/mcresults50.csv")


##plot density results of the hr.ppp- this is in the plotting script.




####Now here we plot the densities of activity centres of each sp.




nsamp<- sample(0:30000, size=500, replace=F) # 500 posterior samples


out<- MCMCpstr(foxesfit, c("g","z"), type="chains")


Sx<- out$g[,1,nsamp]
Sy<- out$g[,2,nsamp]
z<- out$z[,nsamp]

fox.x<- Sx[z==1]
fox.y<- Sy[z==1]




#also get a ppp for plotting
Fox <- factor(rep("Fox", length(fox.x)))


foxppp <- ppp(fox.x,fox.y,marks = Fox,window = cam.region)
foxppp <- as.ppp(foxppp)

###########Repeat dingoes


out<- MCMCpstr(dingoesfitnopup, c("g","z"), type="chains")


Sx<- out$g[,1,nsamp]
Sy<- out$g[,2,nsamp]
z<- out$z[,nsamp]

dingo.x<- Sx[z==1]
dingo.y<- Sy[z==1]



#also get a ppp for plotting
Dingo <- factor(rep("Dingo", length(dingo.x)))

dingoppp <- ppp(dingo.x,dingo.y,M = "Dingo",window = cam.region)
dingoppp <- as.ppp(dingoppp)

#repeat for quolls


out<- MCMCpstr(quollsfit, c("g","z"), type="chains")


Sx<- out$g[,1,nsamp]
Sy<- out$g[,2,nsamp]
z<- out$z[,nsamp]

quoll.x<- Sx[z==1]
quoll.y<- Sy[z==1]



#also get a ppp for plotting
Quoll <- factor(rep("Quoll", length(quoll.x)))

quollppp <- ppp(quoll.x,quoll.y,M = "Quoll",window = cam.region)
quollppp <- as.ppp(quollppp)

#repeat for goannas

out<- MCMCpstr(goannasfit, c("g","z"), type="chains")


Sx<- out$g[,1,nsamp]
Sy<- out$g[,2,nsamp]
z<- out$z[,nsamp]

goanna.x<- Sx[z==1]
goanna.y<- Sy[z==1]





Goanna <- factor(rep("Goanna", length(goanna.x)))

goannappp <- ppp(goanna.x,goanna.y,marks  = Goanna,window = cam.region)
goannappp <- as.ppp(goannappp)

####NOW have all 4 dataframes.





#####Now we'll plot the individual hr. ppp objects as smoothed density surfaces

png("Figures/DensityPlots.png", width = 12, height = 13, res= 300, units = "in")
par(mfrow=c(2,2), mar = c(1, 1, 1, 1))

plot(density.ppp(dingoppp,sigma =1),main="",cex.axis=2)
plot(density.ppp(foxppp,sigma =1),main="",cex.axis=2)
plot(density.ppp(quollppp,sigma =1),main="",cex.axis=2)
plot(density.ppp(goannappp,sigma =1),main="",cex.axis=2)

dev.off()

























###lets classify some areas as locations of dingo risk or not. 

nsamp<- sample(000:30000, size=30000, replace=F) # every posterior sample


out<- MCMCpstr(dingoesfitnopup, c("g","z"), type="chains")


Sx<- out$g[,1,nsamp]
Sy<- out$g[,2,nsamp]
z<- out$z[,nsamp]

dingo.x<- Sx[z==1]
dingo.y<- Sy[z==1]
dingo.x0<- Sx[z==0]
dingo.y0<- Sy[z==0]




dingotest <- data.frame(cbind(x = dingo.x,y = dingo.y,z=1))
dingotest0<- data.frame(cbind(x=dingo.x0,y=dingo.y0,z=0))
dingotestboth <- rbind(dingotest,dingotest0)

dingotestboth <- dingotestboth %>%
  mutate(x = x * 1000, y = y * 1000) %>%
  group_by(x, y) %>%
  summarise(countdingo = sum(z))

dingotestsf <- dingotestboth %>%
  st_as_sf(coords = c("x", "y"),crs=28356)



dingotestboth <- dingotestboth %>%
  mutate(risk_category = cut(countdingo, breaks = c(0,3212,20259), labels = c("Low Risk","High Risk")))

dingotestboth$risk_category <- factor(dingotestboth$risk_category,levels =c("Low Risk","High Risk"))



####Lets determine if the cameras are closer or further to an area of risk. 
trapstemp <- read.csv(file = "Raw Data/utmpredsall.csv")

trapstemp$risk_category <- NA

n_nearest <-5

for (i in 1:nrow(trapstemp)) {
  cam_x <- trapstemp$x[i]
  cam_y <- trapstemp$y[i]
  
  # Calculate Euclidean distances to all points in dingotestboth
  distances <- sqrt((dingotestboth$x - cam_x)^2 + (dingotestboth$y - cam_y)^2)
  
  # Find the indices of the five closest points
  nearest_idx <- order(distances)[1:n_nearest]
  
  # Extract the risk categories of these nearest points
  nearest_risks <- dingotestboth$risk_category[nearest_idx]
  
  # Remove NA values from nearest_risks for checking
  nearest_risks_non_na <- na.omit(nearest_risks)
  
  # Check if any of the nearest points have a high-risk category
  if ("High Risk" %in% nearest_risks_non_na) {
    trapstemp$closest_risk[i] <- "High Risk"
  } else {
    trapstemp$closest_risk[i] <- "Low Risk"
  }
}



trapstemp <- trapstemp %>%
  mutate(closest_risk = case_when(
    closest_risk == 1 ~ "Low Risk",
    closest_risk == 2 ~ "High Risk",
    is.na(closest_risk) ~ "Low Risk",
    TRUE ~ as.character(closest_risk)  # Keep original value if it doesn't match 1 or 2
  ))


write.csv(trapstemp, file = "Derived Data/riskratingtraps.csv",row.names=F)

trapstemp$Trap<- factor(trapstemp$Trap)

ggplot(dingotestboth, aes(x = x, y = y)) +
  geom_point(aes(color = risk_category, size = countdingo), alpha = 0.7) +
  scale_color_manual(values = c("Low Risk" = "green", "High Risk" = "red")) +
  scale_size_continuous(name = "Dingo Count") +
  labs(title = "Spatial Distribution of Dingo Risk Categories",
       x = "Easting",
       y = "Northing",
       color = "Risk Category") +
  theme_minimal() +
  theme(legend.position = "right")+
  geom_point(data =trapstemp, aes(x=x,y=y,fill = Trap),shape =21,fill = "black",size =3)+
  geom_text(data = trapstemp, aes(x = x, y = y, label = Trap), vjust = -1.5,hjust=1.5, color = "black",size =2)
  




#End









