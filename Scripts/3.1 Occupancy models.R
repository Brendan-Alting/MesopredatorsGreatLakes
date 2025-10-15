#3.1 run occupancy models for foxes and goannas.
library(tidyverse)
library(unmarked)
library(sf)

###start by cleaning fox and goanna data: 



foxrawretry <- read.csv(file = "Raw Data/fox15min.csv",header = T)
goannarawretry <- read.csv(file = "Raw Data/goanna15min.csv",header = T)


#also just read in other data. 
camsretry<- read.csv("Raw Data/CamOpMatrixCompleteRevise.csv")

camdata <- st_read("Raw Data/CamsCovs.shp")
#####  covariates for later.
camsretry_updated <- camsretry %>%
  left_join(camdata %>% select(Trap, campdist,urbdist), by = "Trap")



ps_cams <- grep("^PS", camsretry_updated$Trap, value = TRUE)




foxrawretry$DateTime <- as.POSIXct(foxrawretry$DateTime, format = "%Y-%m-%d %H:%M:%S")

foxrawretry <- foxrawretry %>%
  filter(Trap %in% ps_cams) %>%
  arrange(Trap, DateTime) %>%
  group_by(Trap) %>%
  mutate(time_diff = as.numeric(difftime(DateTime, lag(DateTime), units = "mins")),
         keep_row = is.na(time_diff) | time_diff > 30 |No_Adults_Seq>1) %>%
  filter(keep_row) %>%
  ungroup()%>%
  mutate(Day= as.Date(DateTime))

goannarawretry$DateTime <- as.POSIXct(goannarawretry$DateTime, format = "%Y-%m-%d %H:%M:%S")

goannarawretry <- goannarawretry %>%
  filter(Trap %in% ps_cams) %>%
  arrange(Trap, DateTime) %>%
  group_by(Trap) %>%
  mutate(time_diff = as.numeric(difftime(DateTime, lag(DateTime), units = "mins")),
         keep_row = is.na(time_diff) | time_diff > 30 |No_Adults_Seq>1) %>%
  filter(keep_row) %>%
  ungroup()%>%
  mutate(Day= as.Date(DateTime))


##########OCCUPANCY UNMARKED.
####Get into day format. 

fox_daily <- foxrawretry %>%
  group_by(Trap, Day) %>%
  summarise(det = n(), .groups = "drop")


all_days <- seq(min(foxrawretry$Day), max(foxrawretry$Day))


# Expand to full Trap × day grid
fox_matrix <- expand.grid(Trap = ps_cams, Day = all_days) %>%
  left_join(fox_daily, by = c("Trap", "Day")) %>%
  mutate(det = ifelse(is.na(det), 0, det)) %>%
  pivot_wider(names_from = Day, values_from = det)

# Make matrix with rownames = Trap
rownames(fox_matrix) <- fox_matrix$Trap
fox_matrix <- as.matrix(fox_matrix[,-1])



goanna_daily <- goannarawretry %>%
  group_by(Trap, Day) %>%
  summarise(det = n(), .groups = "drop")



# Expand to full Trap × day grid
goanna_matrix <- expand.grid(Trap = ps_cams, Day = all_days) %>%
  left_join(goanna_daily, by = c("Trap", "Day")) %>%
  mutate(det = ifelse(is.na(det), 0, det)) %>%
  pivot_wider(names_from = Day, values_from = det)

# Make matrix with rownames = Trap
rownames(goanna_matrix) <- goanna_matrix$Trap
goanna_matrix <- as.matrix(goanna_matrix[,-1])






site_covs <- camsretry_updated %>%
  filter(Trap %in% ps_cams) %>%
  arrange(match(Trap, rownames(fox_matrix))) %>%
  select(Trap, campdist, urbdist,TrailType) %>%
  mutate(
    scaled_campdist = scale(campdist),
    scaled_urbdist  = scale(urbdist)
  ) %>%
  column_to_rownames("Trap")


####calculate the camera effort.

effort <- camsretry_updated %>%
  mutate(
    Setup_date     = as.Date(Setup_date, format = "%d/%m/%Y"),
    Retrieval_date = as.Date(Retrieval_date, format = "%d/%m/%Y"),
    Problem1_from  = as.Date(Problem1_from, format = "%d/%m/%Y"),
    Problem1_to    = as.Date(Problem1_to,   format = "%d/%m/%Y"),
    Problem2_from  = as.Date(Problem2_from, format = "%d/%m/%Y"),
    Problem2_to    = as.Date(Problem2_to,   format = "%d/%m/%Y")
  )

effort_daily <- expand.grid(Trap = effort$Trap, Day = all_days) %>%
  left_join(effort %>% 
              select(Trap, Setup_date, Retrieval_date, 
                     Problem1_from, Problem1_to, 
                     Problem2_from, Problem2_to),
            by = "Trap") %>%
  mutate(
    effort = ifelse(
      Day >= Setup_date & Day <= Retrieval_date &
        !( !is.na(Problem1_from) & !is.na(Problem1_to) & Day >= Problem1_from & Day <= Problem1_to ) &
        !( !is.na(Problem2_from) & !is.na(Problem2_to) & Day >= Problem2_from & Day <= Problem2_to ),
      1, 0
    )
  ) %>%
  select(Trap, Day, effort)

effortPS <- effort_daily%>%
  filter(grepl("^PS",Trap))

effortmatrix <- effortPS %>%
  pivot_wider(names_from = Day,values_from = effort)

obs_covs_effort <- as.matrix(effortmatrix[ , -1])
rownames(obs_covs_effort) <- effortmatrix$Trap

#make non functioning occasions 'nas'.
fox_matrixwithnas <- fox_matrix
fox_matrixwithnas[obs_covs_effort==0] <- NA

goanna_matrixwithnas <- goanna_matrix
goanna_matrixwithnas[obs_covs_effort==0] <- NA

#########So now should have stuff for occupancy. 























#####

umf_fox <- unmarkedFrameOccu(
  y = fox_matrixwithnas,
  siteCovs = site_covs
)

umf_goanna<- unmarkedFrameOccu(
  y = goanna_matrixwithnas,
  siteCovs = site_covs
)



###First mods null
m0_fox <- occu(~1 ~1, data = umf_fox)
m0_goanna <- occu(~1 ~1, data = umf_goanna)

summary(m0_fox)
summary(m0_goanna)



#models with all variables: 
m2_fox <- occuRN(~ TrailType ~ scaled_urbdist+scaled_campdist, data = umf_fox)
m2_goanna <- occuRN(~ TrailType ~ scaled_urbdist + scaled_campdist, data = umf_goanna)

sumfox <- summary(m2_fox)
sumgoanna <- summary(m2_goanna)

#first we'll collate fox data

fox_coef <- rbind(data.frame(Parameter = rownames(sumfox$state), sumfox$state, Type = "state"),
                  data.frame(Parameter = rownames(sumfox$det),   sumfox$det,   Type = "det")
)



fox_ci_state <- confint(m2_fox, type = "state")
fox_ci_det   <- confint(m2_fox, type = "det")

fox_ci <- rbind(
  data.frame(Parameter = rownames(fox_ci_state), LCI = fox_ci_state[,1], UCI = fox_ci_state[,2], Type = "state"),
  data.frame(Parameter = rownames(fox_ci_det),   LCI = fox_ci_det[,1],   UCI = fox_ci_det[,2],   Type = "det")
)

fox_ci$Parameter <- gsub("^(lam|p)\\(|\\)$", "", fox_ci$Parameter)

# Merge CIs
fox_results <- merge(fox_coef, fox_ci, by = c("Parameter","Type"), all.x = TRUE)
fox_results$Species <- "Fox"



#and repeat for goannas

goanna_coef <- rbind(data.frame(Parameter = rownames(sumgoanna$state), sumgoanna$state, Type = "state"),
                  data.frame(Parameter = rownames(sumgoanna$det),   sumgoanna$det,   Type = "det")
)

goanna_ci_state <- confint(m2_goanna, type = "state")
goanna_ci_det   <- confint(m2_goanna, type = "det")

goanna_ci <- rbind(
  data.frame(Parameter = rownames(goanna_ci_state), LCI = goanna_ci_state[,1], UCI = goanna_ci_state[,2], Type = "state"),
  data.frame(Parameter = rownames(goanna_ci_det),   LCI = goanna_ci_det[,1],   UCI = goanna_ci_det[,2],   Type = "det")
)

goanna_ci$Parameter <- gsub("^(lam|p)\\(|\\)$", "", goanna_ci$Parameter)


# Merge CIs
goanna_results <- merge(goanna_coef, goanna_ci, by = c("Parameter","Type"), all.x = TRUE)
goanna_results$Species <- "Goanna"


model_results <- rbind(fox_results, goanna_results)
model_results <- model_results%>%
  filter(!(Parameter == "(Intercept)"))

write.csv(model_results, "Derived Data/OccupancyResults.csv", row.names = FALSE)





