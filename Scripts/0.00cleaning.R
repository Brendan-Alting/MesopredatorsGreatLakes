library(tidyverse)
#cleaning for analysis file. 

#need to get species sums rather than presence absence. 

dingorawretry <- read.csv(file = "Raw Data/dingo15minnopup.csv",header = T)
quollrawretry <- read.csv(file = "Raw Data/quoll15min.csv",header = T)


#need to clean the dingo/quoll df a bit. Some unidentifiables are same as the detected individuals. 


quollrawretry$DateTime <- as.POSIXct(quollrawretry$DateTime, format = "%Y-%m-%d %H:%M:%S")

quollrawretry <- quollrawretry %>%
  arrange(Trap, DateTime) %>%
  group_by(Trap) %>%
  mutate(time_diff = as.numeric(difftime(DateTime, lag(DateTime), units = "mins")),
         keep_row = is.na(time_diff) | time_diff > 30) %>%
  filter(keep_row) %>%
  ungroup()


dingorawretry$DateTime <- as.POSIXct(dingorawretry$DateTime, format = "%Y-%m-%d %H:%M:%S")

dingorawretry <- dingorawretry %>%
  arrange(Trap, DateTime) %>%
  group_by(Trap) %>%
  mutate(time_diff = as.numeric(difftime(DateTime, lag(DateTime), units = "mins")),
         keep_row = is.na(time_diff) | time_diff > 30 |No_Adults_Seq>1) %>%
  filter(keep_row) %>%
  ungroup()

