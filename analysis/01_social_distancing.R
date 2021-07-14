#01_soc_distancing.R
#------------------------------------------------------------------------------#
# This script processes mobility data from Google's COVID-19 Community Mobili- #
# ty Reports https://www.google.com/covid19/mobility/ to estimate the time     #
# points at which there was a change in mobility (and on contacts) in Mexico   #
# City and the State of Mexico.                                                #
#                                                                              #
# Author:                                                                      #
#       Jeremy D Goldhaber-Fiebert, PhD                                        #
#       Valeria Gracia Olvera, <valeria.gracia@cide.edu>                       #
#------------------------------------------------------------------------------#

rm(list = ls())

# Load libraries  ---------------------------------------------------------

library(tidyverse)
library(segmented)
# library(dampack)

# Load functions, data and set variables  ---------------------------------
# Target functions
source("R/03_target_functions.R")

# Locations for mobility data
v_variables <- c("retail_recreation", 
                 "grocery_and_pharmacy", 
                 "parks", "transit", 
                 "workplaces", "residential"
)

# States of interest
v_states <- c(CMX = "Mexico City",
              MEX = "State of Mexico")

# Number of states
n_states <- length(v_states)

# Date chop
date_chop <- 2

# Targets
l_targets_MCMA <- gen_targets(v_states_calib = "MCMA",
                              n_time_stamp   = "2020-12-13",
                              n_date_last    = "2020-12-07")

# Load and wrangle data ---------------------------------------------------

n_time_stamp_data <- "2020-12-29"

# Google and Foursquare data for social distancing from sccosmo-Data v0.28
load(paste0("data-raw/df_interventions_loc_raw_",n_time_stamp_data,".rda"))

df_interventions_loc <- df_interventions_loc_raw %>%
  group_by(country, state) %>%
  arrange(date) %>%
  mutate(date_0 = date - date[1]) %>%
  ungroup() %>%
  filter(state %in% v_states & date <= "2020-12-13"
  )

# Run piecewise estimation ------------------------------------------------
# Number of changepoints
n_change_points <- 6

# Create array to store results
a_breakpoint_locs <- array(0, 
                           dim = c(n_states, 
                                   n_change_points, 
                                   length(v_variables)),
                           dimnames = list(names(v_states), 
                                           1:n_change_points, 
                                           v_variables))

# Empty object to store results
df_all_out <- data.frame(NULL)

# Initialize outcomes of interest
k <- 1

# Start iterations
for(c_var in v_variables) { # c_var = "grocery_and_pharmacy"
  df_out <- data.frame(NULL)
  for(i in 1:n_states){
    state_i <- v_states[i]
    
    df_temp <- df_interventions_loc %>% filter(state == state_i)
    
    df_out <- tryCatch (
      {
        fit_lm =  lm(formula = eval(parse(text = paste0(c_var, " ~ 1 + date_0"))), data = df_temp)   # intercept-only model
        fit_segmented = segmented(obj = fit_lm, seg.Z = ~date_0, npsi = n_change_points, 
                                  control = seg.control(n.boot=500, seed = 54321))  # n_change_points along date
        
        a_breakpoint <- array(0, dim = c(1, n_change_points))
        
        if (is.null(fit_segmented$psi[,2])) {
          a_breakpoint[1, 1] <- NA
          a_breakpoint[1, 2] <- NA
          a_breakpoint[1, 3] <- NA
          a_breakpoint[1, 4] <- NA
          a_breakpoint[1, 5] <- NA
          a_breakpoint[1, 6] <- NA
        } else {
          a_breakpoint[1, 1] <- fit_segmented$psi[1,2]
          a_breakpoint[1, 2] <- fit_segmented$psi[2,2]
          a_breakpoint[1, 3] <- fit_segmented$psi[3,2]
          a_breakpoint[1, 4] <- fit_segmented$psi[4,2]
          a_breakpoint[1, 5] <- fit_segmented$psi[5,2]
          a_breakpoint[1, 6] <- fit_segmented$psi[6,2]
        }
        
        df_out_temp <- data.frame(a_breakpoint)
        df_out_temp$i <- i
        df_out_temp$k <- k
        df_out_temp$state <- state_i
        df_out <- rbind.data.frame(df_out, df_out_temp)
      }, error=function(cond) {
        message("Here's the original error message:")
        message(cond)
        skip_to_next <- TRUE
        df_out_temp <- data.frame(array(NA, dim = c(1, n_change_points)))
        df_out_temp$i <- i
        df_out_temp$k <- k
        df_out <- rbind.data.frame(df_out, df_out_temp)
      }
      
    )
  }
  df_all_out <- rbind.data.frame(df_all_out, df_out)
  k <- k + 1
  df_all_out
} 

kk <- 1
for(c_var in v_variables) {
  ii <- 1
  for (state_i in v_states) { 
    v_temp <- as.numeric(df_all_out %>% filter(k == kk , i == ii) %>% select(starts_with("X")))
    a_breakpoint_locs[ii,,kk] <- v_temp
    ii <- ii + 1
  }
  kk <- kk + 1
}

# Create variable to indicate to which CP each date belongs
df_interventions_loc_mod <- df_interventions_loc %>%
  mutate(intervention_period = 0)

i <- 1
for (state_i in v_states){
  print(state_i)
  for(j in 1:n_change_points) {
    
    print(paste0("    ", round(median(a_breakpoint_locs[i,j,], na.rm = TRUE))))
    df_interventions_loc_mod <- df_interventions_loc_mod %>%
      mutate(intervention_period = replace(intervention_period,
                                           state == state_i & 
                                             date_0 >= round(median(a_breakpoint_locs[i,j,], na.rm = TRUE)),
                                           j))
    
  } 
  i <- i + 1
}

# Create data.frame with the first date for the intervention period
df_interventions_loc_NPIs <- df_interventions_loc_mod %>%
  group_by(country, state, intervention_period) %>%
  arrange(date) %>%
  filter(date == date[1]) %>%
  ungroup() %>%
  select(country, state, date, intervention_period) %>%
  pivot_wider(values_from = date, names_prefix = "Date_INT", names_from = intervention_period) %>%
  drop_na()

# Summarize intervention dates for MCMA
df_NPIs <- rbind.data.frame(df_interventions_loc_NPIs,
                            data.frame(country = "Mexico",
                                       state   = "MCMA",
                                       Date_INT0 = mean(df_interventions_loc_NPIs$Date_INT0),
                                       Date_INT1 = mean(df_interventions_loc_NPIs$Date_INT1),
                                       Date_INT2 = mean(df_interventions_loc_NPIs$Date_INT2),
                                       Date_INT3 = mean(df_interventions_loc_NPIs$Date_INT3),
                                       Date_INT4 = mean(df_interventions_loc_NPIs$Date_INT4),
                                       Date_INT5 = mean(df_interventions_loc_NPIs$Date_INT5),
                                       Date_INT6 = mean(df_interventions_loc_NPIs$Date_INT6)
                            ))

# Save file
# save(df_NPIs, file = paste0("data/df_NPIs_",n_time_stamp_data,".RData"))

# Targets vs NPIs ---------------------------------------------------------

# Load data
load(paste0("data/df_NPIs_",n_time_stamp_data,".RData"))

# State
state_i <- "MCMA"

# Filter data
df_temp_plot <- l_targets_MCMA$cases_inc

n_time_stamp <- unique(df_temp_plot$time_stamp)
abbrev_state <- unique(df_temp_plot$abbrev_state)

ggplot(df_temp_plot, aes(x = Date, y = value,
                         ymin = lb, ymax = ub,
                         color = type, shape = type)) +
  geom_point(size = 3) + 
  geom_errorbar() +
  geom_vline(xintercept = as.numeric(as.Date(df_NPIs$Date_INT1[df_NPIs$state == "MCMA"])),
             colour = "blue")+
  geom_vline(xintercept = as.numeric(as.Date(df_NPIs$Date_INT2[df_NPIs$state == "MCMA"])),
             colour = "blue")+
  geom_vline(xintercept = as.numeric(as.Date(df_NPIs$Date_INT3[df_NPIs$state == "MCMA"])),
             colour = "blue")+
  geom_vline(xintercept = as.numeric(as.Date(df_NPIs$Date_INT4[df_NPIs$state == "MCMA"])),
             colour = "blue")+
  geom_vline(xintercept = as.numeric(as.Date(df_NPIs$Date_INT5[df_NPIs$state == "MCMA"])),
             colour = "blue")+
  geom_vline(xintercept = as.numeric(as.Date(df_NPIs$Date_INT6[df_NPIs$state == "MCMA"])),
             colour = "blue")+
  scale_shape_manual(values = c(1)) +
  scale_color_manual(values = c("black")) +
  scale_x_date(breaks = number_ticks(16), 
               date_labels = "%b/%d") +
  scale_y_continuous(breaks = number_ticks(8)) +
  labs(title = "Targets vs NPIs",
       subtitle = state_i,
       x = "", y = "Incident cases") +
  theme_bw(base_size = 16) +
  theme(legend.position = "",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 10))

# Save plot
ggsave(paste0("figs/03_targets_NPIs_", abbrev_state,"_",n_time_stamp_data,".png"),
       width = 12, height = 8)
