#plotting secr covariate responses
#take results from 2.3 

library(secr)
library(ggplot2)
library(tidyverse)
library(patchwork)

trail_levels <- unique(covariates(traps(dingoesdethist))$TrailType)


###FIRST DINGOES
###get dataframes 
#urb
newdatadensityploturb <- expand.grid(
  scaledurb = seq(min(covariates(correctmaskclipped)$scaledurb), max(covariates(correctmaskclipped)$scaledurb), length.out = 100),
  scaledcamp = mean(covariates(correctmaskclipped)$scaledcamp),
  TrailType = factor("4wdtrack", levels = trail_levels),
  bk = factor("TRUE", levels = c("TRUE","FALSE")))

#camp
newdatadensityplotcamp <- expand.grid(
  scaledcamp = seq(min(covariates(correctmaskclipped)$scaledcamp), max(covariates(correctmaskclipped)$scaledcamp), length.out = 100),
  scaledurb = mean(covariates(correctmaskclipped)$scaledurb),
  TrailType = factor("4wdtrack", levels = trail_levels),
  bk = factor("TRUE", levels = c("TRUE","FALSE")))

#then predict
#urb
predurbandingo <- predict(fullHNdingo, newdata = newdatadensityploturb, parameter = "D")

#camp
predcampdingo <- predict(fullHNdingo, newdata = newdatadensityplotcamp, parameter = "D")


#get these values into nice frames for plotting
#urban
predicted_values_dingourb <- unlist(sapply(predurbandingo, "[", "D","estimate"))*100
lower_bound_dingourb <- unlist(sapply(predurbandingo, "[", "D","lcl"))*100
upper_bound_dingourb <- unlist(sapply(predurbandingo, "[", "D","ucl"))*100
urbdist_value <- seq(min(covariates(correctmaskclipped)$scaledurb), max(covariates(correctmaskclipped)$scaledurb), length.out = 100)
plot_data_dingourb <- cbind.data.frame(predicted_values_dingourb, lower_bound_dingourb, upper_bound_dingourb, urbdist_value)

#camp
predicted_values_dingocamp <- unlist(sapply(predcampdingo, "[", "D","estimate"))*100
lower_bound_dingocamp <- unlist(sapply(predcampdingo, "[", "D","lcl"))*100
upper_bound_dingocamp <- unlist(sapply(predcampdingo, "[", "D","ucl"))*100
campdist_value <- seq(min(covariates(correctmaskclipped)$scaledcamp), max(covariates(correctmaskclipped)$scaledcamp), length.out = 100)
plot_data_dingocamp <- cbind.data.frame(predicted_values_dingocamp, lower_bound_dingocamp, upper_bound_dingocamp, campdist_value)


####also need to get the unscaled values for plotting: 

urb_mean <- mean(covariates(correctmaskclipped)$urb, na.rm = TRUE)
urb_sd   <- sd(covariates(correctmaskclipped)$urb, na.rm = TRUE)
camp_mean <- mean(covariates(correctmaskclipped)$camp, na.rm = TRUE)
camp_sd   <- sd(covariates(correctmaskclipped)$camp, na.rm = TRUE)

plot_data_dingourb$unscaled_x <- (plot_data_dingourb$urbdist_value * urb_sd) + urb_mean
plot_data_dingourb$covariate <- "Distance from urban centre (m)"
plot_data_dingourb$species <- "Dingo"

plot_data_dingocamp$unscaled_x <- (plot_data_dingocamp$campdist_value * camp_sd) + camp_mean
plot_data_dingocamp$covariate <- "Distance from campground (m)"
plot_data_dingocamp$species <- "Dingo"



##Quolls
###get dataframes 
#urb
newdatadensityploturb <- expand.grid(
  scaledurb = seq(min(covariates(correctmaskclipped)$scaledurb), max(covariates(correctmaskclipped)$scaledurb), length.out = 100),
  scaledcamp = mean(covariates(correctmaskclipped)$scaledcamp),
  TrailType = factor("4wdtrack", levels = trail_levels),
  bk = factor("TRUE", levels = c("TRUE","FALSE")))

#camp
newdatadensityplotcamp <- expand.grid(
  scaledcamp = seq(min(covariates(correctmaskclipped)$scaledcamp), max(covariates(correctmaskclipped)$scaledcamp), length.out = 100),
  scaledurb = mean(covariates(correctmaskclipped)$scaledurb),
  TrailType = factor("4wdtrack", levels = trail_levels),
  bk = factor("TRUE", levels = c("TRUE","FALSE")))

#then predict
#urb
predurbanquoll <- predict(fullEXquoll, newdata = newdatadensityploturb, parameter = "D")

#camp
predcampquoll <- predict(fullEXquoll, newdata = newdatadensityplotcamp, parameter = "D")


#get these values into nice frames for plotting
#urban
predicted_values_quollurb <- unlist(sapply(predurbanquoll, "[", "D","estimate"))*100
lower_bound_quollurb <- unlist(sapply(predurbanquoll, "[", "D","lcl"))*100
upper_bound_quollurb <- unlist(sapply(predurbanquoll, "[", "D","ucl"))*100
urbdist_value <- seq(min(covariates(correctmaskclipped)$scaledurb), max(covariates(correctmaskclipped)$scaledurb), length.out = 100)
plot_data_quollurb <- cbind.data.frame(predicted_values_quollurb, lower_bound_quollurb, upper_bound_quollurb, urbdist_value)

#camp
predicted_values_quollcamp <- unlist(sapply(predcampquoll, "[", "D","estimate"))*100
lower_bound_quollcamp <- unlist(sapply(predcampquoll, "[", "D","lcl"))*100
upper_bound_quollcamp <- unlist(sapply(predcampquoll, "[", "D","ucl"))*100
campdist_value <- seq(min(covariates(correctmaskclipped)$scaledcamp), max(covariates(correctmaskclipped)$scaledcamp), length.out = 100)
plot_data_quollcamp <- cbind.data.frame(predicted_values_quollcamp, lower_bound_quollcamp, upper_bound_quollcamp, campdist_value)


#finally plot
#urban
ggplot(plot_data_quollurb, aes(x = urbdist_value, y = predicted_values_quollurb)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_bound_quollurb, ymax = upper_bound_quollurb), alpha = 0.2)

#camp
ggplot(plot_data_quollcamp, aes(x = campdist_value, y = predicted_values_quollcamp)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_bound_quollcamp, ymax = upper_bound_quollcamp), alpha = 0.2)


plot_data_quollurb$unscaled_x <- (plot_data_quollurb$urbdist_value * urb_sd) + urb_mean
plot_data_quollurb$covariate <- "Distance from urban centre (m)"
plot_data_quollurb$species <- "Quoll"


plot_data_quollcamp$unscaled_x <- (plot_data_quollcamp$campdist_value * camp_sd) + camp_mean
plot_data_quollcamp$covariate <- "Distance from campground (m)"
plot_data_quollcamp$species <- "Quoll"




###plot final: 

plot_secr_all <- bind_rows(
  plot_data_dingourb %>% rename(pred = predicted_values_dingourb,
                                lower = lower_bound_dingourb,
                                upper = upper_bound_dingourb,
                                distance_scale_value = urbdist_value,
                                distance_value = unscaled_x),
  plot_data_dingocamp %>% rename(pred = predicted_values_dingocamp,
                                 lower = lower_bound_dingocamp,
                                 upper = upper_bound_dingocamp,
                                 distance_scale_value = campdist_value,
                                 distance_value = unscaled_x
  ),
  plot_data_quollurb %>% rename(pred = predicted_values_quollurb,
                                lower = lower_bound_quollurb,
                                upper = upper_bound_quollurb,
                                distance_scale_value = urbdist_value,
                                distance_value = unscaled_x),
  plot_data_quollcamp %>% rename(pred = predicted_values_quollcamp,
                                 lower = lower_bound_quollcamp,
                                 upper = upper_bound_quollcamp,
                                 distance_scale_value = campdist_value,
                                 distance_value = unscaled_x)
)

plot_secr_all$covariate <- factor(plot_secr_all$covariate,
                                  levels = c("Distance from campground (m)",
                                             "Distance from urban centre (m)"))

#finally plot

#first dingo
dingosecrplots <- ggplot(plot_secr_all %>% filter(species == "Dingo"),
                         aes(x = distance_value, y = pred)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  facet_wrap(~ covariate, scales = "free", nrow = 1, strip.position = "bottom") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    legend.position = "none"
  ) +
  xlab("") +
  ylab(expression(paste("Density (individuals/", km^2, ")"))) +
  geom_text(
    data = plot_secr_all %>% filter(species == "Dingo") %>% distinct(covariate),
    aes(x = -Inf, y = Inf, label = c("(b) Dingo", "(a) Dingo")),
    hjust = -0.3, vjust = 1.6, size = 5, fontface = "bold",
    inherit.aes = FALSE
  )

dingosecrplots

#quoll
quollsecrplots <- ggplot(plot_secr_all %>% filter(species == "Quoll"),
                         aes(x = distance_value, y = pred)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  facet_wrap(~ covariate, scales = "free", nrow = 1, strip.position = "bottom") +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 17),
    axis.text.y = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    legend.position = "none",
    panel.spacing = unit(1, "lines")
  ) +
  xlab("") +
  ylab(expression(paste("Density (individuals/", km^2, ")"))) +
  geom_text(
    data = plot_secr_all %>% filter(species == "Quoll") %>% distinct(covariate),
    aes(x = -Inf, y = Inf, label = c("(d) Quoll", "(c) Quoll")),
    hjust = -0.3, vjust = 1.6, size = 5, fontface = "bold",
    inherit.aes = FALSE
  )

quollsecrplots

combinedplot <- (dingosecrplots/quollsecrplots)
combinedplot

#and save
ggsave(
  filename = "Figures/combined_density_effects.png",
  plot = combinedplot,
  width = 10,   
  height = 10,  
  dpi = 200     
)

