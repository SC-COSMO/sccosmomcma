#' Generate cumulative mortality risk
#'
#' \code{get_cum_mort_risk} generates a 30- and 60-day cumulative mortality risk
#' of Covid-19 for MCMA. 
#' 
#' @param save_data Logical. Saves data in a csv file if TRUE.
#' @return 
#' A data.frame with cumulative mortality risk expected by 30 and 60 days.
#' @export
get_cum_mort_risk <- function(save_data = TRUE){

# *****************************************************************************
#### 02_Load_data ####
# *****************************************************************************
#Load confirmed 

load("data-raw/ssa_covid_surv_pos_MCMA_2020-12-13.RData")

#Rename data
cov_dic <- ssa_covid_surv_pos

remove(ssa_covid_surv_pos)

# *****************************************************************************
#### 03_Initial_setup ####
# *****************************************************************************
#Create variable with the date of the database

n_time_stamp <- unique(cov_dic$time_stamp)


# *****************************************************************************
#### 03_Wrangle_data ####
# *****************************************************************************

# *****************************************************************************
##### 03.01_Main_data_set #### 
# *****************************************************************************

#Create variable with the day count from the first patient confirmed with the virus
#For MCMA, according to the data reported by the official data, the first case was
#reported on February 24, 2020.

cov_ZM_dic <- cov_dic %>% 
  mutate(day_dx = as.numeric(difftime(cov_dic$date_dx, 
                                      lubridate::ymd( "2020-02-24")),
                             unit = c("days")))

# *****************************************************************************
#### 04 Cox models to get Cumulative Mortality Risk (CMR) for the MCMA ####
# *****************************************************************************


##### 04.01 Wich model we need to choose ##### 

#We develop different models using natural splines and penalised splines to see
#wich one had the best fit. Here are the results using a Bayesian information 
#criterion(BIC) and the Akaike information criterion (BIC):


###### 04.01.01 Survival models######

#mod.cov_10k_ZM <- coxph(Surv(time_dx, dead_ind) ~ ns(day_dx, 10),data = cov_ZM_dic)
#
#mod.cov_7k_ZM <- coxph(Surv(time_dx, dead_ind) ~ ns(day_dx, 7),data = cov_ZM_dic)
#
#mod.cov_5k_ZM <- coxph(Surv(time_dx, dead_ind) ~ ns(day_dx, 5),data = cov_ZM_dic)
#
#
#mod.cov_10k_ZM_ps <- coxph(Surv(time_dx, dead_ind) ~ pspline(day_dx, 10), 
#                           data = cov_ZM_dic)
#
#mod.cov_7k_ZM_ps <- coxph(Surv(time_dx, dead_ind) ~ pspline(day_dx, 7), 
#                          data = cov_ZM_dic)
#
#mod.cov_5k_ZM_ps <- coxph(Surv(time_dx, dead_ind) ~ pspline(day_dx, 5), 
#                          data = cov_ZM_dic)
#

##### 04.02 Penalised splines with 10 degrees of splines ####

#This is our model to get the Cumulative Mortality Risk. Here, we can see that we 
#that we use the function coxph from the package "survival" to fit a survival model

#We can see our survival object "Surv" with the variable "time_dx" that indicates
#the day that the patient was confirmed with the virus starting the count from
#the first confirmed case in MCMA. The variable dead_ind is categorical where 
#dead_ind equals to 1 if the patient died and 0 otherwise. Then we applied our
#penalised spline to the variable day_dx with 10 degrees of splines.

mod.cov_10k_ZM_ps <- coxph(Surv(time_dx, dead_ind) ~ pspline(day_dx, 10), 
                           data = cov_ZM_dic)


##### 04.03 Bayesian information criterion (BIC) ##### 

###### 04.03.01 Polynomial splines ######

#BIC(mod.cov_10k_ZM_ps, mod.cov_7k_ZM_ps, mod.cov_5k_ZM_ps)

#                         df      BIC 
#mod.cov_10k_ZM_ps 10.139884 701424.6 <- This one had the best fit
#mod.cov_7k_ZM_ps   7.140368 701445.2 
#mod.cov_5k_ZM_ps   5.098348 701480.8 

#******************************************************************************
####### 04.03.01_A Best fit: Penalised spline with 10 degrees of spline ####### 
#******************************************************************************

###### 04.03.02 Natural splines ######

#BIC(mod.cov_10k_ZM, mod.cov_7k_ZM, mod.cov_5k_ZM)

#               df      BIC
#mod.cov_10k_ZM 10 701469.4 <- This one had the best fit
#mod.cov_7k_ZM   7 701514.4
#mod.cov_5k_ZM   5 701522.8

#******************************************************************************
####### 04.03.02_B Best fit: Natural spline with 10 knots #######
#******************************************************************************

##### 04.04 Akaike information criterion (AIC)#####

###### 04.04.01 Polynomial splines ######

#AIC(mod.cov_10k_ZM_ps, mod.cov_7k_ZM_ps, mod.cov_5k_ZM_ps)

#                         df      AIC
#mod.cov_10k_ZM_ps 10.139884 701341.0 <- This one had the best fit
#mod.cov_7k_ZM_ps   7.140368 701386.4
#mod.cov_5k_ZM_ps   5.098348 701438.8

#******************************************************************************
####### 04.04.01_A Best fit: Penalised spline with 10 degrees of spline ####### 
#******************************************************************************


###### 04.04.02 Natural splines ######

#AIC(mod.cov_10k_ZM, mod.cov_7k_ZM, mod.cov_5k_ZM)

#               df      AIC
#mod.cov_10k_ZM 10 701387.0 <- This one had the best fit
#mod.cov_7k_ZM   7 701456.7
#mod.cov_5k_ZM   5 701481.6

#******************************************************************************
####### 04.04.02_B Best fit: Natural spline with 10 knots #######
#******************************************************************************

#We can see that the model that had the best fit, according to the AIC and BIC is
#Penalised spline with 10 degrees of spline


# *****************************************************************************
#### 05 New dataframes for the models ####
# *****************************************************************************

#Expand grid to create a new data frame with all the possible combinations of our
#dataset

#We create a variable day_dx developing a sequence starting from the day 1 to the 
#last day of confirmed reported in the official data

df_ND_MR_ZM <- with(cov_ZM_dic, expand.grid(
  day_dx = seq(1, max(day_dx))))

#We create a data frame grouping by the date of dead summarising the number of 
#dates by date

group_death <- cov_ZM_dic %>% group_by(date_dead) %>% 
  summarise(death = sum(dead_ind)) %>% rename(`date`=date_dead) 

#We add the date to our dataframe starting from the first day that a patient was 
#confirmed with the virus till the last day of actualization of the data.

#We put a -7 because of the for the adjustment of any data within that period by 
#the government

df_ND_MR_ZM <- df_ND_MR_ZM %>%
  mutate(date = seq(as.Date("2020-02-24"),as.Date(n_time_stamp) - 7,by = 1))

# *****************************************************************************
#### 06 Objects of the models to get the Cumulative Mortality Risk ####
# *****************************************************************************

#We get the survfit of our model
survMR_10k_ZM_ps <- survfit(mod.cov_10k_ZM_ps, newdata = df_ND_MR_ZM)

#We get the Mortality Risk expected from our models
sum_surv_MR_10k_ZM_30_ps <- summary(survMR_10k_ZM_ps, times = 30)   # day 30
sum_surv_MR_10k_ZM_60_ps <- summary(survMR_10k_ZM_ps, times = 60)   # day 60

# *****************************************************************************
#### 07 Observed CMR ####
# *****************************************************************************

df_MR_obs <- cov_ZM_dic %>% 
  mutate(dead_30d = ifelse(dead_ind & ((date_dead - date_dx) <= 30), 
                           yes = 1,  # If individual die 30 days after detection 
                           no  = 0),
         dead_60d = ifelse(dead_ind & ((date_dead - date_dx) <= 60), 
                           yes = 1,  # If individual die 60 days after detection 
                           no  = 0)) %>%
  group_by(date_dx) %>%
  summarise(idx  = n(), 
            dead = sum(dead_ind), 
            dead_30 = sum(dead_30d),
            dead_60 = sum(dead_60d)) %>%
  complete(date_dx = seq.Date(as.Date("2020-02-24"), 
                              as.Date(n_time_stamp) - 7, 
                              by="day"), # Create a sequence of dates
           fill = list(idx  = 0, 
                       dead = 0, 
                       dead_30 = 0, 
                       dead_60 = 0)) %>% # Fill the dates without cases with zero
  ungroup() %>% 
  mutate(MRObs_30 = ifelse(idx != 0, dead_30/idx, 0),
         MRObs_60 = ifelse(idx != 0, dead_60/idx, 0)) %>% 
  rename(date = date_dx)

df_MR_obs_long <- df_MR_obs %>%
  pivot_longer(cols = starts_with("MRObs_"), names_to = "MR", values_to = "value",
               names_prefix = "MRObs_") %>%
  mutate(MR = case_when(MR == 30 ~ "Mortality Risk Day 30",
                        MR == 60 ~ "Mortality Risk Day 60"),
         type = "Observed",
         day_dx = as.numeric(difftime(date, 
                                      lubridate::ymd( "2020-02-24")),
                             unit = c("days")),
         lb = NA,
         ub = NA) %>%
  select(day_dx, date, type, MR, value, lb, ub)
  

# *****************************************************************************
#### 08 Final data.frame with IC ####
# *****************************************************************************

##We create a final LONG data.frame using the information that we got in our summary
#
#For 30 days
df_MR_10k_30d_ZM <- cbind.data.frame(df_ND_MR_ZM,
                                      data.frame(
                                        type  = "Estimated", 
                                        MR    = "Mortality Risk Day 30",
                                        value = 1 - t(sum_surv_MR_10k_ZM_30_ps$surv),
                                        lb    = 1 - t(sum_surv_MR_10k_ZM_30_ps$upper),
                                        ub    = 1 - t(sum_surv_MR_10k_ZM_30_ps$lower),
                                        check.names = FALSE)
                                     )

#For 60 days
df_MR_10k_60d_ZM <- cbind.data.frame(df_ND_MR_ZM,
                                     data.frame(
                                       type  = "Estimated", 
                                       MR    = "Mortality Risk Day 60",
                                       value = 1 - t(sum_surv_MR_10k_ZM_60_ps$surv),
                                       lb    = 1 - t(sum_surv_MR_10k_ZM_60_ps$upper),
                                       ub    = 1 - t(sum_surv_MR_10k_ZM_60_ps$lower),
                                       check.names = FALSE)
                                     )

# Add observed data
df_MR_30_60_days <- rbind(df_MR_10k_30d_ZM,df_MR_10k_60d_ZM,df_MR_obs_long)

if(save_data){
  # Save a csv file with 30-day and 60-day CMR 
  write.csv(df_MR_30_60_days, 
            file      = paste0("data/MortalityRisk_MCMA_",n_time_stamp,".csv"), 
            row.names = FALSE)
}

return(df_MR_30_60_days)
}

#' Plot expected cumulative mortality risk
#'
#' \code{plot_mort_risk} plot expected cumulative mortality risk of Covid-19 for MCMA. 
#' 
#' @param df_CMR Data.frame. Data of cumulative mortality risk generated by 
#' \code{get_cum_mort_risk}.
#' @param n_MR_day Numeric. Day of expected mortality risk. Options: 30, 60.
#' @param save_plot Logical. Saves the plot if TRUE.
#' @return 
#' A plot a 30- or 60-day cumulative mortality risk.
#' @export
plot_mort_risk <- function(df_CMR,
                           n_MR_day  = 30, 
                           save_plot = FALSE){

  
  if(n_MR_day == 30){
    MR_day <- "Mortality Risk Day 30"
  }else if(n_MR_day == 60){
    MR_day <- "Mortality Risk Day 60"
  }

  # Filter data
  df_plot <- df_CMR %>%
    filter(MR == MR_day)
  
  df_plot$type <- as.factor(df_plot$type)
  
  # Vector of colors
  v_names_colours <- levels(df_plot$type)
  
  if(MR_day == "Mortality Risk Day 30"){
    v_colors_colour <- c("#0F57CA", "#a13023")
    names(v_colors_colour) <- v_names_colours
    
    # Vector of colors to fill
    v_colors_fill <- c("#0F57CA", NA)
    names(v_colors_fill) <- v_names_colours
    
    # Shape
    v_shape <- c(NA, 16)
    names(v_shape) <- v_names_colours
    
    # Line
    v_linetype <- c("solid", "blank")
    names(v_linetype) <- v_names_colours
    
  }else if(MR_day == "Mortality Risk Day 60"){
    v_colors_colour <- c("#4A9B6E", "#a13023")
    names(v_colors_colour) <- v_names_colours
    
    # Vector of colors to fill
    v_colors_fill <- c("#4A9B6E", NA)
    names(v_colors_fill) <- v_names_colours
    
    # Shape
    v_shape <- c(NA, 16)
    names(v_shape) <- v_names_colours
    
    # Line
    v_linetype <- c("solid", "blank")
    names(v_linetype) <- v_names_colours
  }
  
  df_plot_obs <- subset(df_plot, type == "Observed")
  df_plot_hat <- subset(df_plot, type == "Estimated")
  
  plot_MR <- ggplot(df_plot_obs) +
         geom_point(mapping = aes(x = date, 
                                  y = value, 
                                  colour = type,
                                  fill = type),
                    shape = 16,
                    size = 1.5) +
           geom_line(data = df_plot_hat,
                     aes(x = date, 
                         y = value, 
                         colour = type),
                     size = 0.9) +
           geom_ribbon(data = df_plot_hat,
                       aes(x = date,
                           y = value,
                           ymin = lb,
                           ymax = ub,
                           fill = type),
                       alpha = 0.4,
                       colour = NA) +
           scale_color_manual(values = v_colors_colour) +
           scale_fill_manual(values =  v_colors_fill)  +
           guides(fill = guide_legend(title = ""),
                  color = guide_legend(title = "",
                                       ncol = 2,
                                       override.aes = list(
                                         shape    = v_shape,
                                         linetype = v_linetype,
                                         color    = v_colors_colour,
                                         fill     = v_colors_fill))) +
    scale_x_date("",
                 date_labels = "%m/%d", #Change the format of the date in the x axis
                 breaks = number_ticks(12)
    ) + 
    scale_y_continuous("Cumulative mortality risk", 
                       breaks = number_ticks(12)) +
    theme(text = element_text(size = 13),
          plot.caption = element_text(hjust = 0,
                                      size = 12),
          axis.text.x = element_text(angle = 60,
                                     hjust = 1,
                                     colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.background = element_rect(fill = "white",
                                          colour = "gray",
                                          size = 0.15,
                                          linetype = "solid"),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      size = 0.7),
          legend.position = "bottom",
          legend.text = element_text(size = 10),
          legend.justification = "center",
          legend.direction = "horizontal",
          legend.title = element_blank(),
          legend.key = element_rect(fill   = "transparent",
                                    colour = "transparent")) +
    labs(title = MR_day)
  
  if(save_plot){
    if(n_MR_day == 30){
      ggsave(paste0("figs/MortRisk_30day_MCMA.png"), plot = plot_MR,
             width = 10, height = 7)
    }else if(n_MR_day == 60){
      ggsave(paste0("figs/MortRisk_60day_MCMA.png"), plot = plot_MR,
             width = 10, height = 7)
    }
  }
  
  return(plot_MR)
}





