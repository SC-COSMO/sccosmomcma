#01_social_distancing.R
################################################################################
# This script processes Google and Foursquare data for social distancing for   #
# Mexico City and State of Mexico.                                              #
################################################################################

rm(list = ls())

#----------------------------------------------------------------------------#
#### 1. Load libraries                                                   ####
#----------------------------------------------------------------------------#

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


#----------------------------------------------------------------------------#
#### 2. Load functions, data and set variables                           ####
#----------------------------------------------------------------------------#
source("R/03_target_functions.R")

get_cache(reinit = TRUE)

date_chop <- 2
v_variables <- c("retail_recreation", "grocery_and_pharmacy", 
                 "parks", 
                 "transit", "workplaces", "residential"
)
n_change_points <- 8#6
v_states_calib <- c("Mexico City", "State of Mexico")
n_states <- length(v_states_calib)

# load("data/l_targets_2020-12-13.RData")   # targets generated

l_targets_all <- gen_targets(v_states_calib = v_states_calib,
                             n_time_stamp   = "2020-12-13",
                             n_date_last    = "2020-12-07")

### Google and Foursquare data for social distancing
# load("data-raw/mobility-google-mexico_2020-09-01.RData")
df_interventions_loc_raw <-get_interventions(country = "Mexico", 
                                             state   = "", 
                                             county  = "")

df_interventions_loc <- df_interventions_loc_raw %>%
  group_by(country, state) %>%
  arrange(date) %>%
  mutate(date_0 = date - date[1]) %>%
  ungroup() %>%
  filter(state %in% v_states_calib #& date <= "2020-10-25"
         )

#----------------------------------------------------------------------------#
#### 3. Changepoint estimation                                           ####
#----------------------------------------------------------------------------#
a_breakpoint_locs <- array(0, dim = c(n_states, 
                                      n_change_points, 
                                      length(v_variables)))

num_reserved_for_me <- 1
registerDoParallel(max(1, min(detectCores()-num_reserved_for_me,n_states+1)))
df_all_out <- c()
k <- 1
for(c_var in v_variables) { # c_var = "grocery_and_pharmacy"
  df_out <- c()
  opts <- list(attachExportEnv = TRUE)
  df_out <- foreach(parallel_i = 1:n_states,
                    .combine=rbind,
                    .export = ls(globalenv()),
                    .packages=c("segmented",
                                "tidyverse",
                                "dplyr",
                                "lubridate",
                                "stringr",
                                "stringi"
                    )
                    ,
                    .options.snow = opts) %dopar%  { #parallel_i = 1
                      
                      i <- parallel_i
                      state_i <- v_states_calib[i]
                      #                            print(state_i)
                      #                            print(c_var)
                      df_temp <- df_interventions_loc %>% filter(state == state_i)
                      #                              skip_to_next <- FALSE
                      df_out <- tryCatch (
                        {
                          fit_lm =  lm(formula = eval(parse(text = paste0(c_var, " ~ 1 + date_0"))), data = df_temp)   # intercept-only model
                          fit_segmented = segmented(obj = fit_lm, seg.Z = ~date_0, npsi = n_change_points, 
                                                    control = seg.control(n.boot=500, seed = 54321))  # Five change points along date
                          # # summary(fit_segmented)
                          # # plot(fit_segmented)
                          # # points(df_temp)
                          # # lines.segmented(fit_segmented)
                          # # points.segmented(fit_segmented)  
                          # print(min(df_temp$date))
                          # print(fit_segmented$psi[,2])
                          a_breakpoint <- array(0, dim = c(1, 
                                                           n_change_points))
                          if (#skip_to_next == TRUE | 
                            is.null(fit_segmented$psi[,2])) {
                            a_breakpoint[1, 1] <- NA
                            a_breakpoint[1, 2] <- NA
                            a_breakpoint[1, 3] <- NA
                            a_breakpoint[1, 4] <- NA
                            a_breakpoint[1, 5] <- NA
                            a_breakpoint[1, 6] <- NA
                            a_breakpoint[1, 7] <- NA
                            a_breakpoint[1, 8] <- NA
                          } else {
                            a_breakpoint[1, 1] <- fit_segmented$psi[1,2]
                            a_breakpoint[1, 2] <- fit_segmented$psi[2,2]
                            a_breakpoint[1, 3] <- fit_segmented$psi[3,2]
                            a_breakpoint[1, 4] <- fit_segmented$psi[4,2]
                            a_breakpoint[1, 5] <- fit_segmented$psi[5,2]
                            a_breakpoint[1, 6] <- fit_segmented$psi[6,2]
                            a_breakpoint[1, 7] <- fit_segmented$psi[7,2]
                            a_breakpoint[1, 8] <- fit_segmented$psi[8,2]
                          }
                          df_out <- data.frame(a_breakpoint)
                          df_out$i <- i
                          df_out$k <- k
                          df_out$state <- state_i
                          df_out <- df_out
                          df_out
                        }, error=function(cond) {
                          message("Here's the original error message:")
                          message(cond)
                          skip_to_next <- TRUE
                          df_out <- data.frame(array(NA, dim = c(1, 
                                                                 n_change_points)))
                          df_out$i <- i
                          df_out$k <- k
                          df_out <- df_out
                          df_out
                        }
                        
                      )
                      
                      df_out <- df_out
                      df_out
                    }
  df_all_out <- bind_rows(df_all_out, df_out)
  k <- k + 1
  df_all_out
} 
stopImplicitCluster()

kk <- 1
for(c_var in v_variables) {
  ii <- 1
  for (state_i in v_states_calib) { 
    v_temp <- as.numeric(df_all_out %>% filter(k == kk , i == ii) %>% select(starts_with("X")))
    a_breakpoint_locs[ii,,kk] <- v_temp
    ii <- ii + 1
  }
  kk <- kk + 1
}

df_interventions_loc_mod <- df_interventions_loc %>%
  mutate(intervention_period = 0)

df_interventions_mean <- df_interventions_loc_mod

i <- 1
for (state_i in v_states_calib){
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


df_interventions_loc_NPIs <- df_interventions_loc_mod %>%
  group_by(country, state, intervention_period) %>%
  arrange(date) %>%
  filter(date == date[1]) %>%
  ungroup() %>%
  select(country, state, date, intervention_period) %>%
  pivot_wider(values_from = date, names_prefix = "Date_INT", names_from = intervention_period) %>%
  mutate(flag_late_soc_int = (Date_INT8 > (max(l_targets_all$cases_inc$Date) - date_chop - 2*7))) %>%
  drop_na()

### Save file
save(df_interventions_loc_NPIs, file = paste0("data/df_interventions_loc_",n_time_stamp,".RData"))


### Plot
df_interventions_traj  <- df_interventions_loc %>%
  select(state, date, retail_recreation:residential) %>%
  pivot_longer(cols = retail_recreation:residential,
               names_to = "Location", values_to = "Mobility") %>%
  group_by(week = cut(date, "week"), Location, state) %>% 
  summarise(mean_mob = mean(Mobility)) %>%
  ungroup() %>%
  left_join(df_interventions_loc_NPIs, by = "state")

ggplot(df_interventions_traj, 
       aes(x = week, y = mean_mob,
           color = as.factor(Location),
           group = as.factor(Location))) +
  geom_point() + geom_line() + 
  xlab("Week") + 
  ylab("Mobility variation (%)") +
  theme_bw(base_size = 12) + 
  geom_vline(aes(xintercept = week(Date_INT2)-6),
             data = df_interventions_traj,
             colour = "#393F3E",
             linetype = "twodash") +
  geom_vline(aes(xintercept = week(Date_INT4)-6),
             data = df_interventions_traj,
             colour = "#393F3E",
             linetype = "dashed") +
  geom_vline(aes(xintercept = week(Date_INT5)-6),
             data = df_interventions_traj,
             colour = "#393F3E",
             linetype = "twodash") +
  geom_vline(aes(xintercept = week(Date_INT6)-6),
             data = df_interventions_traj,
             colour = "#393F3E",
             linetype = "dashed") +
  geom_vline(aes(xintercept = week(Date_INT7)-6),
             data = df_interventions_traj,
             colour = "#393F3E",
             linetype = "twodash") +
  geom_vline(aes(xintercept = week(Date_INT8)-6),
             data = df_interventions_traj,
             colour = "#393F3E",
             linetype = "dashed") +
  facet_wrap(~state) +
  scale_colour_discrete(labels=c("Grocery and pharmacy","Parks",
                                 "Residential","Retail and recreation",
                                 "Transit","Workplaces"),) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position="bottom", legend.title=element_blank(),
        text = element_text(size = 16))

ggsave(file = "figs/Mobility_traj_NPIs_MCMA.png", width = 14, height = 8)
ggsave(file = "figs/Mobility_traj_NPIs_MCMA.pdf", width = 14, height = 8)


#----------------------------------------------------------------------------#
#### 4. Plot NPIs vs targets                                              ####
#----------------------------------------------------------------------------#

for(state_i in v_states_calib){ # state_i <- "Aguascalientes"
  print(state_i)
  
  df_alltargets_state <- l_targets_all$df_all_targets %>%
    filter(state == state_i)
  
  n_time_stamp <- unique(df_targets$cases$time_stamp)
  abbrev_state <- unique(df_alltargets_state$abbrev_state)
  
  gg_targets <- ggplot(df_alltargets_state, aes(x = Date, y = value,
                                                ymin = lb, ymax = ub,
                                                color = type, shape = type)) +
    facet_wrap( ~ series , scales = "free_y") +
    geom_point(size = 3) + # for target: shape = 8
    geom_errorbar() +
    geom_vline(xintercept = as.numeric(as.Date(df_interventions_loc_NPIs$Date_INT0[df_interventions_loc_NPIs$state == state_i])),
               colour = "blue")+
    geom_vline(xintercept = as.numeric(as.Date(df_interventions_loc_NPIs$Date_INT1[df_interventions_loc_NPIs$state == state_i])),
               colour = "blue")+
    geom_vline(xintercept = as.numeric(as.Date(df_interventions_loc_NPIs$Date_INT2[df_interventions_loc_NPIs$state == state_i])),
               colour = "blue")+
    geom_vline(xintercept = as.numeric(as.Date(df_interventions_loc_NPIs$Date_INT3[df_interventions_loc_NPIs$state == state_i])),
               colour = "blue")+
    geom_vline(xintercept = as.numeric(as.Date(df_interventions_loc_NPIs$Date_INT4[df_interventions_loc_NPIs$state == state_i])),
               colour = "blue")+
    geom_vline(xintercept = as.numeric(as.Date(df_interventions_loc_NPIs$Date_INT5[df_interventions_loc_NPIs$state == state_i])),
               colour = "blue")+
    geom_vline(xintercept = as.numeric(as.Date(df_interventions_loc_NPIs$Date_INT6[df_interventions_loc_NPIs$state == state_i])),
               colour = "blue")+
    geom_vline(xintercept = as.numeric(as.Date(df_interventions_loc_NPIs$Date_INT7[df_interventions_loc_NPIs$state == state_i])),
               colour = "blue")+
    geom_vline(xintercept = as.numeric(as.Date(df_interventions_loc_NPIs$Date_INT8[df_interventions_loc_NPIs$state == state_i])),
               colour = "blue")+
    scale_shape_manual(values = c(1)) +
    scale_color_manual(values = c("black")) +
    scale_x_date(breaks = number_ticks(16)) +
    scale_y_continuous(breaks = number_ticks(8)) +
    # labs(title = “Cumulative total infectious cases of COVID-19 Mexico”,
    #      subtitle = “Days since first case in each state”,
    #      x = “Days since first case”, y = “Cumulative cases”) +
    theme_bw(base_size = 16) +
    theme(legend.position = "",
          axis.text.x = element_text(angle = 60, hjust = 1, size = 10))
  ggsave(paste0("figs/03_targets_NPIs_", abbrev_state, "_",n_time_stamp,".pdf"), 
         plot = gg_targets,
         width = 12, height = 8)
}

