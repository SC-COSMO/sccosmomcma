# 01_soc_distancing.R
#------------------------------------------------------------------------------#
# This script processes Google Mobility data to estimate where effective       #
# changes on contacts occurred for Mexico City and the State of Mexico.        #
#                                                                              #
# Author:                                                                      #
#     Valeria Gracia Olvera, <valeria.gracia@cide.edu>                         #
#------------------------------------------------------------------------------#

rm(list = ls())


# Load libraries  ---------------------------------------------------------

library(epitools)
library(dampack)
library(stringi)
library(countrycode)
library(usmap)
library(sccosmoData)
#library(data.table)
library(tidyverse)
library(stringi)
library(lubridate)
library(segmented)
library(gridExtra)
library(grid)
library(doParallel)
library(foreach)
library(dampack)



# Load functions, data and set variables ----------------------------------

get_cache(reinit = TRUE)

# Locations in Google Mobility Data
v_variables <- c("retail_recreation", 
                 "grocery_and_pharmacy", 
                 "parks", "transit", 
                 "workplaces", "residential"
)

# Vector of states
v_states <- c(CMX = "Mexico City",
                    MEX = "State of Mexico")

# Number of states
n_states <- length(v_states)

# # ### Google mobility data for social distancing from (sccosmo-Data v0.27)
# df_interventions_loc_raw <- get_interventions(country = "Mexico",
#                                               state   = "",
#                                               county  = "")

# save(df_interventions_loc_raw, file = "data/df_interventions_loc_raw_2020-12-01.rda")

# Load data
load("data/df_interventions_loc_raw_2020-12-01.rda")

df_interventions_loc <- df_interventions_loc_raw %>%
  group_by(country, state) %>%
  arrange(date) %>%
  mutate(date_0 = date - date[1]) %>%
  ungroup() %>%
  filter(state %in% v_states)

# Mobility for each state
df_mob_state <- df_interventions_loc %>%
  filter(state %in% v_states #& date <= "2020-12-07"
  ) %>%
  select(state, date, retail_recreation:residential) %>%
  pivot_longer(cols = retail_recreation:residential,
               names_to = "Location", values_to = "Mobility") %>%
  mutate(date_0    = as.numeric(date - date[1]),
         s_code    = ifelse(state == "Mexico City", "CMX", "MEX")) %>%
  select(state, s_code, date, date_0, Location, Mobility)

# Order data set
df_mob_state <- df_mob_state[order(df_mob_state$state,df_mob_state$date),]

# Mobility by state for all outcomes except "residential" and "workplaces"
df_mob_state_mean_all <- df_mob_state %>%
  filter(Location != "residential" & Location != "workplaces") %>%
  group_by(date, state) %>%
  summarise(mean_mob = mean(Mobility)) %>%
  ungroup() %>%
  mutate(date_0    = as.numeric(date - date[1]),
         Location  = "All",
         predicted = NA) %>%
  select(state, date, date_0, Location, mean_mob, predicted)


# All ---------------------------------------------------------------------

df_NPIs_long <- data.frame(NULL)

for(state_i in v_states){
  # Select state
  print(state_i)
  abbrev_state <- names(v_states)[which(v_states == state_i)]
  
  # for(n_change_points in 4:8){
  
  if(state_i == "State of Mexico"){
    n_change_points <- 6 # 7 AIC 2343.61
  }else if(state_i == "Mexico City"){
    n_change_points <- 7 # 7 AIC 2320.84
  }
  #----------------------------------------------------------------------------#
  #### 3. Changepoint estimation                                           ####
  #----------------------------------------------------------------------------#
  
  fit_lm <-   lm(mean_mob ~ 1 + date_0,
                 data = df_mob_state_mean_all[df_mob_state_mean_all$state == state_i,])   # intercept-only model
  fit_segmented <-  segmented(obj = fit_lm, seg.Z = ~ 1 + date_0, 
                              # psi = list(date_0 = c(8,14,32,36,62,64,130,230)) #, 
                              npsi = n_change_points, 
                              control = seg.control(n.boot=500, seed = 54321)
  )
  
  print(paste0("AIC model for ",state_i," with ",n_change_points," changepoints ", AIC(fit_segmented)))
  
  # summary(fit_segmented)
  # plot(fit_segmented)
  # points(df_temp)
  # lines.segmented(fit_segmented)
  # points.segmented(fit_segmented)
  # print(min(df_temp$date))
  # print(fit_segmented$psi[,2])
  
  df_mob_state_mean_all$predicted[df_mob_state_mean_all$state == state_i] <- predict(fit_segmented)
  # }
  
  cp_estimates <- as.Date(fit_segmented$psi[,2]+as.Date("2020-02-15"))
  
  df_NPIs_long <- rbind.data.frame(df_NPIs_long,
                                   data.frame(state = state_i,
                                              name = paste0("Date_NPI", 1:n_change_points),
                                              cp_estimates, row.names = NULL))
}


ggplot(df_mob_state_mean_all, aes(x = date)) +
  geom_point(aes(y = mean_mob, colour = "#666666"), alpha = 0.5) +
  geom_line(aes(y = mean_mob, colour = "#666666"), linetype = "dashed", alpha = 0.5) +
  geom_line(aes(y = predicted, colour = "firebrick2"), size = 0.9) +
  facet_wrap(~state) +
  ylab("Mean mobility variation (%)") +
  xlab("") +
  scale_x_date(date_labels = "%m/%d",
               breaks = number_ticks(12)) +
  theme_bw(base_size = 12) +
  scale_colour_manual(name="",
                      # labels map onto colors and pretty labels
                      values=c("#666666", "firebrick2"),
                      labels=c("Observed","Estimated")) +
  theme(axis.text.x = element_text(angle = 60, 
                                   hjust = 1),
        strip.background = element_rect(fill   = "transparent",
                                        colour = "transparent"),
        legend.position="bottom", 
        legend.title=element_blank())

# Save plot
ggsave(file = "figs/Mobility_traj_NPIs_MCMA_2020-12-01.png", width = 14, height = 7)


## Long wide data.frame with cp
df_NPIs <- df_NPIs_long %>%
  pivot_wider(values_from = cp_estimates, names_from = name) %>%
  rbind.data.frame(data.frame(
    state = "MCMA",
    Date_NPI1 = median(df_NPIs$Date_NPI1, na.rm = TRUE),
    Date_NPI2 = median(df_NPIs$Date_NPI2, na.rm = TRUE),
    Date_NPI3 = median(df_NPIs$Date_NPI3, na.rm = TRUE),
    Date_NPI4 = median(df_NPIs$Date_NPI4, na.rm = TRUE),
    Date_NPI5 = median(df_NPIs$Date_NPI5, na.rm = TRUE),
    Date_NPI6 = median(df_NPIs$Date_NPI6, na.rm = TRUE),
    Date_NPI7 = median(df_NPIs$Date_NPI7, na.rm = TRUE)
  ))

save(df_NPIs, file = "data/df_NPIs_2020-12-01.RData")
save(df_mob_state_mean_all, file = "data/df_mob_MCMA_obs_estim_2020-12-01.RData")
