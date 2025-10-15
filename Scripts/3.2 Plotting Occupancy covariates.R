#3.2 plotting occupancy models results

library(ggplot2)
library(gridExtra)

trail_levels <- levels(factor(site_covs$TrailType))

urbanpredict <- expand.grid(
  scaled_urbdist = seq(min(site_covs$scaled_urbdist), max(site_covs$scaled_urbdist), length.out = 100),
  scaled_campdist = mean(site_covs$scaled_campdist),
  TrailType = factor("4wdtrack", levels = trail_levels))

#camp
camppredict <- expand.grid(
  scaled_campdist = seq(min(site_covs$scaled_campdist), max(site_covs$scaled_campdist), length.out = 100),
  scaled_urbdist = mean(site_covs$scaled_urbdist),
  TrailType = factor("4wdtrack", levels = trail_levels))

pred_fox_urb <- predict(m2_fox, type = "state", newdata = urbanpredict)
pred_fox_camp <- predict(m2_fox, type = "state", newdata = camppredict)
pred_goanna_urb <- predict(m2_goanna, type = "state", newdata = urbanpredict)
pred_goanna_camp <- predict(m2_goanna, type = "state", newdata = camppredict)


#and make a sequence of the unscaled covariates. 

urbanpredict<- urbanpredict%>%
  mutate(distance_value = seq(min(site_covs$urbdist),max(site_covs$urbdist),length.out=100),
         covariate = "Distance from urban centre (m)")%>%
  select(-c(scaled_campdist,TrailType))%>%
  rename(distance_scale_value = scaled_urbdist)


camppredict <- camppredict%>%
  mutate(distance_value = seq(min(site_covs$campdist),max(site_covs$campdist),length.out=100),
         covariate = "Distance from campground (m)")%>%
  select(-c(scaled_urbdist,TrailType))%>%
  rename(distance_scale_value = scaled_campdist)

#add names for each predicted so we know which is which
pred_fox_urb <- pred_fox_urb%>%
  mutate(species = "Fox")
pred_fox_camp <- pred_fox_camp%>%
  mutate(species = "Fox")
pred_goanna_urb <- pred_goanna_urb%>%
  mutate(species = "Goanna")
pred_goanna_camp <- pred_goanna_camp%>%
  mutate(species = "Goanna")


pred_fox_urb_plus_dists <- bind_cols(pred_fox_urb,urbanpredict)
pred_fox_camp_plus_dists <- bind_cols(pred_fox_camp,camppredict)
pred_goanna_urb_plus_dists <- bind_cols(pred_goanna_urb,urbanpredict)
pred_goanna_camp_plus_dists <- bind_cols(pred_goanna_camp,camppredict)

plot_occ_all <- bind_rows(pred_fox_urb_plus_dists,pred_fox_camp_plus_dists,pred_goanna_urb_plus_dists,pred_goanna_camp_plus_dists)


#first lets do lil foxes. 
foxoccplots <- ggplot(plot_occ_all %>% filter(species == "Fox"),
                      aes(x = distance_value, y = Predicted)) +
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
  ylab(expression(paste("Expected abundance at site"))) +
  geom_text(
    data = plot_occ_all %>% filter(species == "Fox") %>% distinct(covariate),
    aes(x = -Inf, y = Inf, label = c("(b) Fox", "(a) Fox")),
    hjust = -0.3, vjust = 1.6, size = 5, fontface = "bold",
    inherit.aes = FALSE
  )

foxoccplots


###now for oannagays

goannaoccplots <- ggplot(plot_occ_all %>% filter(species == "Goanna"),
                         aes(x = distance_value, y = Predicted)) +
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
  ylab(expression(paste("Expected abundance at site"))) +
  geom_text(
    data = plot_occ_all %>% filter(species == "Goanna") %>% distinct(covariate),
    aes(x = -Inf, y = Inf, label = c("(d) Lace monitor", "(c) Lace monitor")),
    hjust = -0.2, vjust = 1.6, size = 5, fontface = "bold",
    inherit.aes = FALSE
  )

goannaoccplots

combinedoccplot <- (foxoccplots/goannaoccplots)
combinedoccplot

ggsave(
  filename = "Figures/combined_occupancy_effects.png",  
  plot = combinedoccplot,
  width = 10,   
  height = 10,  
  dpi = 200     
)
