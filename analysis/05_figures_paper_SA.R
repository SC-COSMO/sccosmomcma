# 05_figures_SM_paper.R
#-------------------------------------------------------------------------# 
# This script creates the figures presented for the sensitivity analysis  # 
# presented in the MCMA paper.                                            #
#                                                                         # 
# Authors:                                                                #
#     - Valeria Gracia Olvera, MsC, <valeria.gracia@cide.edu>             #
#-------------------------------------------------------------------------# 

rm(list=ls()) # empty the environment


# Load packages and set constants -----------------------------------------

# Libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(dampack)
library(ggpubr)
library(grid)

# Dates
n_date_projs <- "2020-12-13"
n_date_obs <- "2020-12-13"
n_time_stamp_hosp <- "2021-03-07"

# # Combine data.frames for both scenarios
# df_out_total_prob <- data.frame(NULL)
# df_out_total_prob_all <- data.frame(NULL)
# 
# for(n_proj_type in c("SQ", "SA")){
# 
#   load(paste0("output/05_projections_probabilistic_all_",n_proj_type,"_MCMA_", n_date_projs, ".RData"))
#   load(paste0("output/05_projections_probabilistic_",n_proj_type,"_MCMA_", n_date_projs, ".RData"))
# 
#   df_out_total_prob <- rbind.data.frame(df_out_total_prob,
#                                         df_out_mex_total_prob)
# 
#   df_out_total_prob_all <- rbind.data.frame(df_out_total_prob_all,
#                                             df_out_mex_total_prob_all)
# 
# }
# 
# save(df_out_total_prob_all, file =paste0("output/05_projections_probabilistic_all_MCMA_",n_date_projs,".RData"))
# save(df_out_total_prob, file = paste0("output/05_projections_probabilistic_MCMA_",n_date_projs,".RData"))


# Load data ---------------------------------------------------------------

# Projections
load(paste0("output/05_projections_probabilistic_all_MCMA_",n_date_projs,".RData"))
load(paste0("output/05_projections_probabilistic_MCMA_",n_date_projs,".RData"))

# Observed data
load("data/l_targets_2020-12-13.RData")

# Hospitalization data
# load(paste0("data/df_hosp_MCMA.RData"))
load(paste0("data/df_hosp_MCMA_2020-12-21.RData"))


# Plot specifications  ----------------------------------------------------

# Names of outcomes
v_outcomes <- c(Hosp_tot      = "Total hospitalizations",
                CDR_prop      = "CDR proportion",
                Dcov_tot      = "COVID deaths",
                DXCumtot      = "Cumulative cases",
                InfCumtot     = "Cumulative infections",
                InfCumprop    = "Cumulative infections proportion",
                DXtot         = "Detected cases",
                DXDcov_tot    = "Cumulative deaths",
                NonICU_tot    = "Hospitalizations without ventilator",
                ICU_tot       = "Hospitalizations with ventilator",
                DXDIncCov_tot = "Incident deaths",
                DXInctot      = "Incident cases",
                InfInctot     = "Incident infections",
                Inftot        = "Prevalent infections",
                Rt            = "R effective",
                Rec_tot       = "Recovered prevalence",
                Rec_prop      = "Recovered prevalence proportion")

# State
state_i <- unique(df_out_total_prob$state)

# Abbreviation for state
abbrev_state <- unique(df_out_total_prob$s_code)

# Data time stamp
n_time_stamp <- as.Date(unique(df_out_total_prob$time_stamp))

# Last date
n_date_last <- as.Date("2020-12-07")

# Max date
max_date <- as.Date("2021-03-07")

# NPI date
date_NPI <- as.Date("2020-03-23")

# NPI update
date_NPI_update <- as.Date("2021-01-10")

# Holidays dates
n_holidays_start <- as.Date("2020-12-24")

# Text size
text_size <- 24
caption_size <- 20
annotation_size <- 17
bracket_size <- 1.1
label_bracket_size <- 6.8
arrow_size <- 3

# X axis formating
angle_x_axis <- 90
label_x_axis <- ""

# Date breaks
v_dates_breaks <- c(as.Date("2020-02-15"), as.Date("2020-03-01"), 
                    as.Date("2020-03-15"), as.Date("2020-04-01"), 
                    as.Date("2020-04-15"), as.Date("2020-05-01"),
                    as.Date("2020-05-15"), as.Date("2020-06-01"), 
                    as.Date("2020-06-15"), as.Date("2020-07-01"), 
                    as.Date("2020-07-15"), as.Date("2020-08-01"),
                    as.Date("2020-08-15"), as.Date("2020-09-01"), 
                    as.Date("2020-09-15"), as.Date("2020-10-01"),
                    as.Date("2020-10-15"), as.Date("2020-11-01"),
                    as.Date("2020-11-15"), as.Date("2020-12-01"), 
                    as.Date("2020-12-15"), as.Date("2021-01-01"), 
                    as.Date("2021-01-15"), as.Date("2021-02-01"),
                    as.Date("2021-02-15"), as.Date("2021-03-01"))

# Set colors
jet.colors <- colorRampPalette(c("black", "#00007F", "blue", "#007FFF",
                                 "cyan", "#7FFF7F", "yellow", "#FF7F00",
                                 "red", "#7F0000"))
color_map  <-  jet.colors(100)

# Function 
# Number of ticks for plots: dampack package
number_ticks <- function(n) {
  function(limits) {
    pretty(limits, n + 1)
  }
}


# Wrangle data ------------------------------------------------------------

# Set vector of interventions and labels
v_interv_names <- c(BaseCase           = "Policy A. Physical distancing: status quo; Schooling: not in-person",
                    IncreaseSD         = "Policy B. Physical distancing: +26% (SA: +24%) compared to status quo; Schooling: not in-person",
                    IncreaseSDSchoolSD = "Policy C. Physical distancing: +26% (SA: +24%) compared to status quo; Schooling: in-person",
                    SchoolSD           = "Policy D. Physical distancing: status quo; Schooling: in-person")

# Rename interventions
df_out_total_prob$Intervention <- as.character(df_out_total_prob$Intervention)
df_out_total_prob_all$Intervention <- as.character(df_out_total_prob_all$Intervention)

for(n_interv in names(v_interv_names)){ 
  interv_type <- v_interv_names[which(names(v_interv_names) == n_interv)]
  
  df_out_total_prob$Intervention[df_out_total_prob$intervention_type == n_interv] <- interv_type
  df_out_total_prob_all$Intervention[df_out_total_prob_all$intervention_type == n_interv] <- interv_type
}

df_out_total_prob$Intervention <- factor(df_out_total_prob$Intervention,
                                         levels = v_interv_names,
                                         ordered = TRUE)

df_out_total_prob_all$Intervention <- factor(df_out_total_prob_all$Intervention,
                                             levels = v_interv_names,
                                             ordered = TRUE)

# Projection type
v_proj_labels <- c(SQ = "Base-case",
                   SA = "Lower Covid-19 risk in children")

# Projection type names
df_out_total_prob$proj_type_label <- ""
df_out_total_prob_all$proj_type_label <- ""
for(n_proj_type in c("SQ", "SA")){ 
  df_out_total_prob$proj_type_label[df_out_total_prob$proj_type == n_proj_type] <- v_proj_labels[n_proj_type]
  df_out_total_prob_all$proj_type_label[df_out_total_prob_all$proj_type == n_proj_type] <- v_proj_labels[n_proj_type]
}

df_out_total_prob$proj_type_label <- factor(df_out_total_prob$proj_type_label,
                                            levels = v_proj_labels,
                                            ordered = TRUE)

df_out_total_prob_all$proj_type_label <- factor(df_out_total_prob_all$proj_type_label,
                                                levels = v_proj_labels,
                                                ordered = TRUE)

# Create variable with base case labels
v_basecase_names <- c(StatusQuo = "End-of-year Holidays: status quo",
                      Holidays  = "End-of-year Holidays: higher social contacts")

df_out_total_prob$BaseCase <- ""
df_out_total_prob_all$BaseCase <- ""

for(n_basecase in names(v_basecase_names)){ 
  basecase_type <- v_basecase_names[which(names(v_basecase_names) == n_basecase)]
  df_out_total_prob$BaseCase[df_out_total_prob$BaseCase_type == n_basecase] <- basecase_type
  df_out_total_prob_all$BaseCase[df_out_total_prob_all$BaseCase_type == n_basecase] <- basecase_type
}

df_out_total_prob$BaseCase <- factor(df_out_total_prob$BaseCase,
                                     levels = v_basecase_names,
                                     ordered = TRUE)

df_out_total_prob_all$BaseCase <- factor(df_out_total_prob_all$BaseCase,
                                         levels = v_basecase_names,
                                         ordered = TRUE)

# Rename observed outcomes
df_out_total_prob$Outcome <- as.character(df_out_total_prob$Outcome)
df_out_total_prob$Outcome[df_out_total_prob$Outcome == "Cumulative detected cases"] <- "Cumulative cases"
df_out_total_prob$Outcome[df_out_total_prob$Outcome == "Incident detected cases"] <- "Incident cases"
df_out_total_prob$Outcome[df_out_total_prob$Outcome == "Detected COVID deaths"] <- "Cumulative deaths" 
df_out_total_prob$Outcome[df_out_total_prob$Outcome == "Incident COVID19 deaths infections"] <- "Incident deaths"

df_out_total_prob_all$Outcome <- as.character(df_out_total_prob_all$Outcome)
df_out_total_prob_all$Outcome[df_out_total_prob_all$Outcome == "Cumulative detected cases"] <- "Cumulative cases"
df_out_total_prob_all$Outcome[df_out_total_prob_all$Outcome == "Incident detected cases"] <- "Incident cases"
df_out_total_prob_all$Outcome[df_out_total_prob_all$Outcome == "Detected COVID deaths"] <- "Cumulative deaths" 
df_out_total_prob_all$Outcome[df_out_total_prob_all$Outcome == "Incident COVID19 deaths infections"] <- "Incident deaths"


## Offset hospitalizations ------------------------------------------------

offset_flag <- TRUE

# Modify data applying an offset
n_date_offset <- "2020-12-07"

# Outcomes of interest
v_hosp_outcomes <- c("Total hospitalizations",
                     "Hospitalizations without ventilator",
                     "Hospitalizations with ventilator")

if(offset_flag){
  for(n_proj_type in c("SQ", "SA")){
    for(n_hosp_outcome in c(v_hosp_outcomes)){
      denominator <- df_hosp_MCMA$hosp_tot[df_hosp_MCMA$dates == n_date_offset] 
      
      offset_hosp_ratio <- mean(df_out_total_prob_all$value[df_out_total_prob_all$dates == n_date_offset & 
                                                              df_out_total_prob_all$BaseCase_type == "StatusQuo" &
                                                              df_out_total_prob_all$proj_type == n_proj_type &
                                                              df_out_total_prob_all$intervention_type == "BaseCase" &
                                                              df_out_total_prob_all$Outcome == n_hosp_outcome])/denominator
      
      
      # Modify data 
      df_out_total_prob_all$value[df_out_total_prob_all$proj_type == n_proj_type &
                                    df_out_total_prob_all$Outcome == n_hosp_outcome] <- df_out_total_prob_all$value[df_out_total_prob_all$proj_type == n_proj_type &
                                                                                                                      df_out_total_prob_all$Outcome == n_hosp_outcome]/offset_hosp_ratio
      df_out_total_prob_all$lb[df_out_total_prob_all$proj_type == n_proj_type &
                                 df_out_total_prob_all$Outcome == n_hosp_outcome] <- df_out_total_prob_all$lb[df_out_total_prob_all$proj_type == n_proj_type &
                                                                                                                df_out_total_prob_all$Outcome == n_hosp_outcome]/offset_hosp_ratio
      df_out_total_prob_all$ub[df_out_total_prob_all$proj_type == n_proj_type &
                                 df_out_total_prob_all$Outcome == n_hosp_outcome] <- df_out_total_prob_all$ub[df_out_total_prob_all$proj_type == n_proj_type &
                                                                                                                df_out_total_prob_all$Outcome == n_hosp_outcome]/offset_hosp_ratio
      
    }
  }
}

# Covid-19 cases and deaths observed data
df_covid_obs_data <- rbind(l_targets$cases,
                           l_targets$cases_inc,
                           l_targets$deaths,
                           l_targets$deaths_inc) %>%
  ungroup() %>%
  filter(state %in% unique(df_out_total_prob$state)) %>%
  rename(Outcome = Target,
         s_code = abbrev_state,
         dates = Date,
         sd = se) %>%
  select(-country, -county, -var_outcome, -series, -population, -Date0) %>%
  mutate(type = "Observed",
         BaseCase_type = "Observed",
         Intervention = "Observed",
         intervention_type = "Observed",
         median = NA)

df_covid_obs_data <- df_covid_obs_data[,c("state", "s_code", "type", "BaseCase_type", "Intervention",
                                          "intervention_type", "Outcome", "dates", "value", "median",
                                          "sd", "lb", "ub","time_stamp")]

df_covid_obs_data$Outcome <- as.character(df_covid_obs_data$Outcome)
df_covid_obs_data$Outcome[df_covid_obs_data$Outcome == "Incident confirmed infections"] <- "Incident cases"
df_covid_obs_data$Outcome[df_covid_obs_data$Outcome == "Cumulative confirmed infections"] <- "Cumulative cases"
df_covid_obs_data$Outcome[df_covid_obs_data$Outcome == "Cumulative COVID19 deaths infections"] <-"Cumulative deaths" 
df_covid_obs_data$Outcome[df_covid_obs_data$Outcome == "Incident COVID19 deaths infections"] <- "Incident deaths"


# Weekly data
df_covid_obs_data_week <- rbind(l_targets$cases,
                                l_targets$cases_inc,
                                l_targets$deaths,
                                l_targets$deaths_inc) %>%
  ungroup() %>%
  filter(state %in% unique(df_out_total_prob$state) &
           Date <= n_date_last) %>%
  rename(Outcome = Target,
         s_code = abbrev_state,
         dates = Date,
         sd = se) %>%
  select(-country, -county, -var_outcome, -series, -population, -Date0) %>%
  mutate(type = "Observed",
         BaseCase_type = "Observed",
         Intervention = "Observed",
         intervention_type = "Observed",
         median = NA) %>%
  group_by(state, type, Outcome, DateLast, time_stamp, s_code,
           BaseCase_type, Intervention, intervention_type, 
           week = cut.Date(dates, breaks = "week")) %>%
  summarise(value_week = mean(value),
            lb_week = mean(lb),
            ub_week = mean(ub),
            sd_week = mean(sd)) %>%
  rename(value = value_week,
         lb = lb_week,
         ub = ub_week,
         sd = sd_week,
         dates = week) %>%
  mutate(median = NA) %>%
  ungroup()


df_covid_obs_data_week <- df_covid_obs_data_week[,c("state", "s_code", "type", "BaseCase_type", "Intervention",
                                                    "intervention_type", "Outcome", "dates", "value", "median",
                                                    "sd", "lb", "ub","time_stamp")]

df_covid_obs_data_week$Outcome <- as.character(df_covid_obs_data_week$Outcome)
df_covid_obs_data_week$Outcome[df_covid_obs_data_week$Outcome == "Incident confirmed infections"] <- "Incident cases"
df_covid_obs_data_week$Outcome[df_covid_obs_data_week$Outcome == "Cumulative confirmed infections"] <- "Cumulative cases"
df_covid_obs_data_week$Outcome[df_covid_obs_data_week$Outcome == "Cumulative COVID19 deaths infections"] <-"Cumulative deaths" 
df_covid_obs_data_week$Outcome[df_covid_obs_data_week$Outcome == "Incident COVID19 deaths infections"] <- "Incident deaths"


# Figure: comparison between approaches ---------------------------------

# Interventions
select_intervention <- c("SchoolSD",
                         "IncreaseSD",
                         "IncreaseSDSchoolSD",
                         "BaseCase")

# Outcomes
select_outcome <- c("Incident cases",
                    "Incident deaths")

# Scenarios
select_basecase <- c("Holidays",
                     "StatusQuo")

# Create data.frame
df_figX <-  df_out_total_prob_all %>%
  filter((dates > "2020-12-07" & dates <= "2021-03-07") &
           Outcome %in% select_outcome &
           BaseCase_type %in% select_basecase &
           intervention_type %in% select_intervention) %>%
  group_by(county, type, BaseCase, BaseCase_type, Intervention, intervention_type, 
           Outcome, proj_type, proj_type_label, simulation) %>%
  summarise(cum_value = sum(value)) %>%
  mutate(time_stamp = n_time_stamp) %>%
  rename(state = county,
         value = cum_value) %>%
  ungroup()


# Summarise data
df_figX_summ <- df_figX %>%
  group_by(state, type, BaseCase, BaseCase_type, Intervention, intervention_type, 
           Outcome, proj_type, proj_type_label) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            median = quantile(value, probs = 0.50, names = FALSE),
            Q1 = quantile(value, probs = 0.25, names = FALSE),
            Q3 = quantile(value, probs = 0.75, names = FALSE),
            IQR = Q3 - Q1,
            min = Q1 - 1.5*IQR,
            max = Q3 + 1.5*IQR,
            lb = quantile(value, probs = 0.025, names = FALSE),
            ub = quantile(value, probs = 0.975, names = FALSE)) %>%
  rename(value = mean)

# Save data
save(df_figX, file = "figs/figs_SA/data_frames/df_figX_Diff_SQ_SA.RData")


## Plot specifications ----------------------------------------------------

# Set vector of interventions and labels
v_interv_labels <- c(BaseCase           = "Policy A",
                     IncreaseSD         = "Policy B",
                     IncreaseSDSchoolSD = "Policy C",
                     SchoolSD           = "Policy D" 
)

# Caption
plot_caption <- ""
for(intv in unique(df_figX$Intervention)){
  plot_caption <- paste0(plot_caption, paste0(intv,"\n"))
}

# Colors
v_names_colours <- c(as.character(levels(df_figX$proj_type_label)))
v_colors_colour <- c("#0F57CA", "#008450")
names(v_colors_colour) <- v_names_colours

# Vector of colors to fill
v_colors_fill <- v_colors_colour
names(v_colors_fill) <- v_names_colours


## Plot -------------------------------------------------------------------
l_plot_figX <- vector(mode = "list", length = 4)
i <- 1

for(n_outcome in select_outcome){
  for(n_policy in unique(df_figX$intervention_type)){
    
    abbrev_plot <- names(v_outcomes)[which(v_outcomes == n_outcome)]
    
    # Plot's title
    if(n_outcome == "Incident cases"){
      gg_title <- "Cumulative cases from 12/07 until 03/07"
    }else if(n_outcome == "Incident deaths"){
      gg_title <- "Cumulative deaths from 12/07 until 03/07"
    }
    
    # Plot subtitle
    gg_subtitle <- v_interv_names[which(names(v_interv_names) == n_policy)]
    
    # Plot
    l_plot_figX[[i]] <- ggplot(subset(df_figX, Outcome == n_outcome & 
                                        intervention_type == n_policy), 
                               aes(x     = proj_type_label, 
                                   y     = value,
                                   color = proj_type_label)) +
      geom_boxplot(size = 1) +
      facet_wrap( ~ BaseCase, ncol = 2, scales = "fixed") +
      scale_colour_manual("", values = v_colors_colour) +
      scale_x_discrete("", labels = NULL) +
      scale_y_continuous("", 
                         breaks = number_ticks(5),
                         labels = scales::comma) + # trans = "log"
      labs(subtitle = gg_subtitle) +
      theme(text = element_text(size = text_size),
            plot.caption = element_text(hjust = 0,
                                        size = 19),
            axis.text.x = element_text(angle = 0, 
                                       hjust = 0.5, 
                                       vjust = 0.5,
                                       colour = "black"),
            axis.text.y = element_text(colour = "black"),
            panel.background = element_rect(fill = "white", 
                                            colour = "gray", 
                                            size = 0.15, 
                                            linetype = "solid"),
            panel.border = element_rect(colour = "black", 
                                        fill = NA, 
                                        size = 0.7), 
            strip.background = element_rect(fill   = "transparent",
                                            colour = "transparent"),
            legend.position = "top", 
            legend.justification = "center", 
            legend.direction = "horizontal", 
            legend.title = element_blank(),
            legend.text = element_text(size = 22),
            legend.key = element_rect(fill   = "transparent", 
                                      colour = "transparent",
                                      size   = unit(4, "cm")))
    
    if(n_policy == "BaseCase"){
      l_plot_figX[[i]] <- l_plot_figX[[i]] +
        labs(title = gg_title)
    }
    
    i <- i + 1 
  }
  
  i <- 1
  
  plot_figX <- ggarrange(l_plot_figX[[i]],
                         NULL,
                         l_plot_figX[[i+1]],
                         NULL,
                         l_plot_figX[[i+2]],
                         NULL,
                         l_plot_figX[[i+3]],
                         ncol = 1,
                         common.legend = T,
                         legend = "bottom",
                         heights = c(1,0.05, 1,0.05,1,0.05,1),
                         # labels = c("Policy A","",
                         #            "Policy B","", 
                         #            "Policy C","", 
                         #            "Policy D"),
                         hjust = 0,vjust = 1.5,
                         font.label = list(size = 26, face = "bold"))
  
  # Save plot
  ggsave(file = paste0("figs/figs_SA/figX_Diff_",abbrev_plot,"_SA_SQ.jpg"),
         width = 18, height = 24, dpi = 300)
  
}


# FIGURE 1: cases and deaths ----------------------------------------------

# Vector of outcomes of interest
v_outcome <- c("Incident cases",
               "Incident deaths",
               "Cumulative cases",
               "Cumulative deaths")

# Wrangle projections
df_outcome_int <- df_out_total_prob_all %>%
  ungroup() %>%
  filter(county == state_i &
           Outcome %in% v_outcome &
           intervention_type == "BaseCase" &
           dates <= max_date &
           proj_type == "SA") %>%
  group_by(county, type, BaseCase, BaseCase_type, Intervention, 
           intervention_type, Outcome, dates, proj_type, proj_type_label) %>%
  summarise(mean = mean(value),
            median = quantile(value, probs = 0.5, names = FALSE),
            sd = sd(value),
            lb = quantile(value, probs = 0.025, names = FALSE),
            ub = quantile(value, probs = 0.975, names = FALSE)) %>%
  rename(value = mean,
         state = county) %>%
  mutate(time_stamp = n_time_stamp) %>%
  select(state, type, BaseCase, BaseCase_type, Intervention, 
         intervention_type, Outcome, dates, value, lb, ub, time_stamp, 
         proj_type, proj_type_label) %>%
  ungroup()

df_outcome_int$Intervention <- "Model-predicted"

# Observed data
df_outcome_obs <- rbind.data.frame(
  subset(df_covid_obs_data_week, Intervention == "Observed") %>%
    mutate(BaseCase_type   = "Holidays",
           BaseCase        = v_basecase_names["Holidays"],
           proj_type       = "SA",
           proj_type_label = "Lower Covid-19 risk in children"),
  subset(df_covid_obs_data_week, Intervention == "Observed") %>%
    mutate(BaseCase_type   = "StatusQuo",
           BaseCase        = v_basecase_names["StatusQuo"],
           proj_type       = "SA",
           proj_type_label = "Lower Covid-19 risk in children")) %>%
  select(state, type, BaseCase, BaseCase_type, Intervention, 
         intervention_type, Outcome, dates, value, lb, ub, time_stamp, 
         proj_type, proj_type_label)

df_outcome_obs$Intervention <- as.factor(df_outcome_obs$Intervention)

# Create and save data.frame
df_fig1_CasesDeaths <- rbind.data.frame(df_outcome_int,
                                        df_outcome_obs) 

# Save data
save(df_fig1_CasesDeaths,
     file = paste0("figs/figs_SA/data_frames/df_fig1_CasesDeaths.RData"))


## Plot specifications ----------------------------------------------------

# Colors
v_names_colours_obs <- as.character(levels(df_fig1_CasesDeaths$Intervention))
v_colors_colour <- rev(c("#a13023","#659583"))
names(v_colors_colour) <- v_names_colours_obs

# Vector of colors to fill
v_colors_fill <- rev(c(NA,"#659583"))
names(v_colors_fill) <- v_names_colours_obs

# Shape
v_shape <- rev(c(16,NA))
names(v_shape) <- v_names_colours_obs

# Line
v_linetype <- rev(c("blank","solid"))
names(v_linetype) <- v_names_colours_obs

# For annotation_custom
npi_date <- grobTree(textGrob("Date of NPI implementation", 
                              x = 0.12, 
                              y = 0.97, 
                              hjust = 0,
                              gp = gpar(col = "#393F3E", 
                                        fontsize = annotation_size)))

calib_day <- grobTree(textGrob("Last day used for calibration", 
                               x = 0.756, 
                               y = 0.97, 
                               hjust = 0,
                               gp = gpar(col = "#393F3E", 
                                         fontsize = annotation_size)))

## Plot -------------------------------------------------------------------

# Save plots in list
l_plots_fig1 <- vector(mode = "list", length = 4)
i <- 1

for(n_proj_type in c("SA")){ # n_proj_type = "SA"
  
  for(outcome in v_outcome){
    
    # Abbreviation for outcome to save the plot
    abbrev_outcome <- names(v_outcomes)[which(v_outcomes == outcome)]
    
    # Filter data to outcome and type of data
    df_obs_temp <- subset(df_fig1_CasesDeaths, 
                          Outcome == outcome & 
                            Intervention == "Observed" &
                            proj_type == n_proj_type)
    df_int_temp <- subset(df_fig1_CasesDeaths,  
                          Outcome == outcome & 
                            Intervention == "Model-predicted" &
                            proj_type == n_proj_type)
    
    # Bracket position
    y_end <- max(df_obs_temp$ub,df_int_temp$ub)
    
    # Plot
    l_plots_fig1[[i]] <- ggplot(data = df_obs_temp)  +
      geom_point(mapping = aes(x = dates, 
                               y = value, 
                               colour = Intervention,
                               fill = Intervention),
                 shape = 16,
                 size = 3,
                 alpha = 0.6) +
      geom_line(data = df_int_temp,
                aes(x = dates, 
                    y = value, 
                    colour = Intervention),
                size = 1.1) + #0.9) +
      geom_ribbon(data = df_int_temp,
                  aes(x = dates,
                      y = value,
                      ymin = lb,
                      ymax = ub,
                      fill = Intervention),
                  alpha = 0.4,
                  colour = NA) +
      facet_wrap(~BaseCase, ncol = 2) +
      scale_color_manual(values = rep(v_colors_colour,2)) +
      scale_fill_manual(values =  rep(v_colors_fill, 2))  +
      guides(fill = guide_legend(title = ""),
             color = guide_legend(title = "",
                                  ncol = 2,
                                  override.aes = list(shape    = v_shape,
                                                      linetype = v_linetype,
                                                      color    = v_colors_colour,
                                                      fill     = v_colors_fill))) +
      geom_vline(xintercept = c(as.numeric(as.Date("2020-03-23")), 
                                as.numeric(n_date_last)),
                 show.legend = TRUE,
                 size = 0.5, 
                 linetype = rep(c("dashed", "twodash"), 2),
                 color = "#393F3E") +
      annotation_custom(npi_date) +
      geom_bracket(data = NULL, xmin = n_date_last,
                   xmax = max_date,
                   y.position = y_end*1.04,
                   label = "Projections",
                   label.size = label_bracket_size,
                   size = bracket_size) +
      scale_x_date(label_x_axis,
                   date_labels = "%b/%d", #Change the format of the date in the x axis
                   breaks = v_dates_breaks
      ) + 
      scale_y_continuous(outcome, 
                         breaks = number_ticks(12),
                         labels = scales::comma) + # trans = "log"
      theme(text = element_text(size = text_size),
            axis.text.x = element_text(angle = angle_x_axis, 
                                       hjust = 0.5, 
                                       vjust = 0.5,
                                       colour = "black"),
            axis.text.y = element_text(colour = "black"),
            panel.background = element_rect(fill = "white", 
                                            colour = "gray", 
                                            size = 0.15, 
                                            linetype = "solid"),
            panel.border = element_rect(colour = "black", 
                                        fill = NA, 
                                        size = 0.7), 
            strip.background = element_rect(fill   = "transparent",
                                            colour = "transparent"),
            legend.position = "bottom", 
            legend.justification = "center", 
            legend.text = element_text(size = 22),
            legend.direction = "horizontal", 
            legend.title = element_blank(), 
            legend.key = element_rect(fill   = "transparent", 
                                      colour = "transparent",
                                      size   = unit(3, "cm")))  + 
      coord_cartesian(ylim = c(0, y_end*1.1))
    
    i <- i + 1
    
  }
  
  i <- 1
  
  # Combine graphs
  plot_inc <- ggarrange(l_plots_fig1[[i]], 
                        NULL,
                        l_plots_fig1[[i+1]],
                        ncol = 1,
                        common.legend = T, 
                        legend = "none",
                        heights = c(1,0.05, 1),
                        labels = c("A","","B"), 
                        hjust = 0,vjust = 1.5,
                        font.label = list(size = 26, face = "bold"))
  # Cumulative
  plot_cum <- ggarrange(l_plots_fig1[[i+2]], 
                        NULL,
                        l_plots_fig1[[i+3]],
                        ncol = 1,
                        common.legend = T, 
                        legend = "bottom",
                        heights = c(1,0.05, 1),
                        labels = c("C","","D"), 
                        hjust = 0,vjust = 1.5,
                        font.label = list(size = 26, face = "bold"))
  
  plot_fig1 <- ggarrange(plot_inc,
                         plot_cum,
                         ncol = 1)
  
  # Save plot
  ggsave(paste0("figs/figs_SA/fig1_CasesDeaths_",n_proj_type,".jpg"),
         width = 18, height = 32, dpi = 300)
  
}

# FIGURE 2: Hospitalizations ----------------------------------------------

# Vector of outcomes of interest
outcome <- "Total hospitalizations"

# Offset flag
offset_flag <- TRUE

# Hosp projected data
df_Hosp_proj <- df_out_total_prob_all %>%
  ungroup() %>%
  filter(county == state_i &
           intervention_type == "BaseCase" &
           Outcome == outcome &
           dates <= max_date &
           proj_type == "SA") %>%
  group_by(county, type, BaseCase, BaseCase_type, Intervention, 
           intervention_type, Outcome, dates, proj_type, proj_type_label) %>%
  summarise(mean = round(mean(value),0),
            median = round(quantile(value, probs = 0.5, names = FALSE),0),
            sd = round(sd(value),0),
            lb = round(quantile(value, probs = 0.025, names = FALSE),0),
            ub = round(quantile(value, probs = 0.975, names = FALSE),0)) %>%
  rename(value = mean,
         state = county) %>%
  mutate(time_stamp = n_time_stamp,
         cap_tot = NA) %>%
  select(state, type, BaseCase, BaseCase_type, Intervention, 
         intervention_type, Outcome, dates, value, lb, ub, cap_tot, time_stamp,
         proj_type, proj_type_label)

df_Hosp_proj$Intervention <- "Model-predicted"

# Observed data
df_hosp_tot <- df_hosp_MCMA %>%
  filter(dates <= n_date_last)
df_hosp_tot$dates <- as.Date(df_hosp_tot$dates)
df_hosp_tot <- df_hosp_tot[order(df_hosp_tot$dates,decreasing = F),]

df_hosp_tot <- df_hosp_tot %>%
  mutate(state = state_i,
         type = "Observed",
         Outcome = outcome,
         Intervention = "Observed",
         intervention_type = "Observed",
         lb = NA,
         ub = NA, time_stamp = as.Date(n_date_obs)
  ) %>%
  filter(dates <= n_time_stamp_hosp) %>%
  rename(value = hosp_tot) 

df_hosp_tot$Intervention <- as.factor(df_hosp_tot$Intervention)

# Create data.frame with two set of base-case
df_hosp_observed <- rbind.data.frame(df_hosp_tot %>% mutate(BaseCase_type   = "Holidays",
                                                            BaseCase        = v_basecase_names["Holidays"],
                                                            proj_type       = "SA",
                                                            proj_type_label = "Lower Covid-19 risk in children"),
                                     df_hosp_tot %>% mutate(BaseCase_type   = "StatusQuo",
                                                            BaseCase        = v_basecase_names["StatusQuo"],
                                                            proj_type       = "SA",
                                                            proj_type_label = "Lower Covid-19 risk in children")) %>%
  select(state, type, BaseCase, BaseCase_type, Intervention, 
         intervention_type, Outcome, dates, value, lb, ub, cap_tot, time_stamp,
         proj_type, proj_type_label
  )

# Bind data.frames
df_fig2_HospTot <- rbind.data.frame(df_Hosp_proj,
                                    df_hosp_observed)

# Save
if(offset_flag){
  save(df_fig2_HospTot, 
       file = paste0("figs/figs_SA/data_frames/df_fig2_HospTot_offset.RData"))
  
}else{
  save(df_fig2_HospTot, 
       file = paste0("figs/figs_SA/data_frames/df_fig2_HospTot.RData"))
}


## Plot specifications ----------------------------------------------------

# Order interventions
v_interv <- as.character(unique(df_fig2_HospTot$Intervention))
df_fig2_HospTot$Intervention <- ordered(df_fig2_HospTot$Intervention,
                                        rev(v_interv))

# Colors
v_names_colours_obs <- as.character(levels(df_fig2_HospTot$Intervention))
v_colors_colour <- c("#a13023", "#659583")
names(v_colors_colour) <- v_names_colours_obs

# Vector of colors to fill
v_colors_fill <- c("#a13023","#659583")
names(v_colors_fill) <- v_names_colours_obs

# Line
v_linetype <- c("solid","solid")
names(v_linetype) <- v_names_colours_obs

# For annotation_custom
npi_date <- grobTree(textGrob("Date of NPI implementation", 
                              x = 0.12, 
                              y = 0.97, 
                              hjust = 0,
                              gp = gpar(col = "#393F3E", 
                                        fontsize = annotation_size)))

calib_day <- grobTree(textGrob("Last day used for calibration", 
                               x = 0.675, 
                               y = 0.97, 
                               hjust = 0,
                               gp = gpar(col = "#393F3E", 
                                         fontsize = annotation_size)))

## Plot -------------------------------------------------------------------

# List of plots
l_plot_fig2 <- vector(mode = "list", length = 2)
i <- 1

for(n_proj_type in c("SA")){
  
  # Abbreviation for outcome to save the plot
  abbrev_outcome <- names(v_outcomes)[which(v_outcomes == outcome)]
  
  # Filter outcome of interest
  df_obs_temp <- subset(df_fig2_HospTot, 
                        Intervention == "Observed" &
                          proj_type == n_proj_type)
  
  df_int_temp <- subset(df_fig2_HospTot, 
                        Intervention == "Model-predicted" &
                          proj_type == n_proj_type)
  
  # Bracket position
  y_end <- max(df_obs_temp$value,df_int_temp$ub)
  
  # Plot
  l_plot_fig2[[i]] <- ggplot(data = df_obs_temp) +
    geom_area(position = "stack",
              aes(x = dates, 
                  y = value, 
                  colour = Intervention, 
                  fill = Intervention),
              alpha = 0.4) +
    geom_line(data = df_int_temp,#subset(df_outcome_int, dates >= max(df_hosp_tot$dates)),
              aes(x = dates, 
                  y = value, 
                  colour = Intervention),
              linetype = "solid") +
    geom_ribbon(data = df_int_temp, #subset(df_outcome_int, dates >= max(df_hosp_tot$dates)),
                aes(x = dates,
                    y = value,
                    ymin = lb,
                    ymax = ub,
                    fill = Intervention),
                alpha = 0.4,
                colour = NA) +
    geom_line(data = df_obs_temp,
              aes(x = dates,
                  y = cap_tot),
              size = 1.2,
              color = "#393F3E") +
    facet_wrap(~BaseCase, ncol = 2) +
    scale_color_manual(values = rep(v_colors_colour, 2)) +
    scale_fill_manual(values =  rep(v_colors_fill, 2))  +
    # scale_color_manual(values = v_colors_colour) + 
    # scale_fill_manual(values =  v_colors_fill)  +
    guides(fill = guide_legend(title = "",
                               reverse = T),
           color = guide_legend(title = "", 
                                ncol = 2,
                                override.aes = list(color    = v_colors_colour,
                                                    fill     = v_colors_fill),
                                reverse = T)) +
    geom_vline(xintercept = c(as.numeric(as.Date("2020-03-23")), as.numeric(n_date_last)),
               show.legend = TRUE,
               size = 0.5, 
               linetype = rep(c("dashed", "twodash"), 2),
               # linetype = c("dashed", "twodash"),
               color = "#393F3E") +
    annotation_custom(npi_date) +
    # annotation_custom(calib_day) +
    annotate("text",
             label  = paste0("Current capacity: ",
                             last(df_obs_temp$cap_tot)," beds"),
             x      = as.Date("2020-09-13"),
             y      = last(df_obs_temp$cap_tot) + 2000,
             size   = 6,
             colour = "#393F3E") +
    geom_bracket(data = NULL, xmin = n_date_last,
                 xmax = max_date,
                 y.position = y_end*1.04,
                 label = "Projections",
                 label.size = label_bracket_size,
                 size = bracket_size) +
    scale_x_date(label_x_axis,
                 date_labels = "%b/%d", #Change the format of the date in the x axis
                 breaks = v_dates_breaks,
    ) + 
    scale_y_continuous(outcome, 
                       breaks = number_ticks(8),
                       labels = scales::comma) + # trans = "log"
    # labs(title    = "COVID-19 hospitalizations with ventilator in MCMA, Mexico" # , caption  = gg_caption
    #      ) + 
    theme(text = element_text(size = text_size),
          plot.caption = element_text(hjust = 0,
                                      size = caption_size), 
          axis.text.x = element_text(angle = 90, 
                                     hjust = 0.5, 
                                     vjust = 0.5,
                                     colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.background = element_rect(fill = "white", 
                                          colour = "gray", 
                                          size = 0.15, 
                                          linetype = "solid"),
          panel.border = element_rect(colour = "black", 
                                      fill = NA, 
                                      size = 0.7), 
          strip.background = element_rect(fill   = "transparent",
                                          colour = "transparent"),
          legend.position = "bottom", 
          legend.justification = "center", 
          legend.direction = "vertical", 
          legend.text = element_text(size = 22),
          legend.title = element_blank(), 
          legend.key = element_rect(fill   = "transparent", 
                                    colour = "transparent",
                                    size = unit(3, "cm")))  + 
    coord_cartesian(ylim = c(0, y_end*1.1))
  
  # Save plot
  if(offset_flag){
    ggsave(paste0("figs/figs_SA/fig2_HospTot_offset_",n_proj_type,".jpg"),
           width = 16, height = 9, dpi = 300)
  }else{
    ggsave(paste0("figs/figs_SA/fig2_HospTot_",n_proj_type,".jpg"),
           width = 16, height = 9, dpi = 300)
  }
  
  i <- i + 1
  
}

# FIGURE 3: interventions cases and deaths --------------------------------

# Outcomes
v_outcome     <-  c("Incident deaths", 
                    "Incident cases")


# Interventions
select_intervention <- c("SchoolSD", 
                         "IncreaseSD",
                         "IncreaseSDSchoolSD",
                         "BaseCase")

# Data.frame
df_fig3_CasesDeaths_interv <- df_out_total_prob_all %>%
  ungroup() %>%
  filter(county == state_i &
           Outcome %in% v_outcome &
           intervention_type %in% select_intervention &
           dates > n_date_last &
           dates <= max_date &
           proj_type == "SA") %>%
  group_by(county, type, BaseCase, BaseCase_type, Intervention, 
           intervention_type, Outcome, dates, proj_type, proj_type_label) %>%
  summarise(mean = mean(value),
            median = quantile(value, probs = 0.5, names = FALSE),
            sd = sd(value),
            lb = quantile(value, probs = 0.025, names = FALSE),
            ub = quantile(value, probs = 0.975, names = FALSE)) %>%
  rename(value = mean)

# Set vector of interventions and labels
v_interv_names <- c(BaseCase           = "Policy A. Physical distancing: status quo; Schooling: not in-person",
                    IncreaseSD         = "Policy B. Physical distancing: +24% compared to status quo; Schooling: not in-person",
                    IncreaseSDSchoolSD = "Policy C. Physical distancing: +24% compared to status quo; Schooling: in-person",
                    SchoolSD           = "Policy D. Physical distancing: status quo; Schooling: in-person")

# Rename interventions
df_fig3_CasesDeaths_interv$Intervention <- as.character(df_fig3_CasesDeaths_interv$Intervention)

for(n_interv in names(v_interv_names)){ 
  interv_type <- v_interv_names[which(names(v_interv_names) == n_interv)]
  
  df_fig3_CasesDeaths_interv$Intervention[df_fig3_CasesDeaths_interv$intervention_type == n_interv] <- interv_type
}

df_fig3_CasesDeaths_interv$Intervention <- factor(df_fig3_CasesDeaths_interv$Intervention,
                                                  levels = v_interv_names,
                                                  ordered = TRUE)

# Order outcomes
v_outcom <- as.character(unique(df_fig3_CasesDeaths_interv$Outcome))
df_fig3_CasesDeaths_interv$Outcome <- ordered(df_fig3_CasesDeaths_interv$Outcome,
                                              v_outcom[c(2,1)])

# Save data.frame
save(df_fig3_CasesDeaths_interv,
     file = paste0("figs/figs_SA/data_frames/df_fig3_CasesDeaths_interv.RData"))


## Plot specifications ----------------------------------------------------

# Colors
v_names_colours <- c(as.character(levels(df_fig3_CasesDeaths_interv$Intervention)))
v_colors_colour <- c("#0F57CA", "#008450", "#FF9900", "#B81D13")
names(v_colors_colour) <- v_names_colours

# Vector of colors to fill
v_colors_fill <- v_colors_colour
names(v_colors_fill) <- v_names_colours

# Vector of colors to shape
v_colors_shape <- c(15, 16, 17, 18)
names(v_colors_shape) <- v_names_colours

# Bracket position
y_end <- max(df_fig3_CasesDeaths_interv$value)

NPI_update <- grobTree(textGrob("Policy start date",
                                x = 0.14,
                                y = 0.97,
                                hjust = 0,
                                gp = gpar(col = "#393F3E",
                                          fontsize = annotation_size)))

## Plot -------------------------------------------------------------------

# List of plots
i <- 1
l_plots_fig3 <- vector(mode = "list", length = 4)

for(n_proj_type in c("SA")){
  
  # Plot incident cases
  l_plots_fig3[[i]] <- ggplot(
    data = subset(df_fig3_CasesDeaths_interv, Outcome == "Incident cases" & 
                    proj_type == n_proj_type),
    aes(x = dates,
        y = value,
        colour = Intervention,
        shape = Intervention)) +
    # geom_bar(stat = "identity", position = position_dodge(width = 5),
    #          width = 3) +
    geom_line(size = 0.9) + 
    # geom_point(size = 3) +
    facet_wrap(~BaseCase, ncol = 2) +
    scale_color_manual(values = rep(v_colors_colour, 2)) +
    scale_shape_manual(values = rep(v_colors_shape, 2)) +
    geom_vline(xintercept = c(as.numeric(date_NPI_update)),
               show.legend = TRUE,
               size = 0.5, 
               linetype = c("dashed"),
               # linetype = c("dashed", "twodash"),
               color = "#393F3E") +
    annotation_custom(NPI_update) +
    scale_x_date(label_x_axis,
                 date_labels = "%b/%d", #Change the format of the date in the x axis
                 breaks = "week") +
    scale_y_continuous("Cases",#"Difference to case base scenario",
                       breaks = number_ticks(8),
                       # labels = scales::percent
                       labels = scales::comma) + # trans = "log"
    theme(text = element_text(size = text_size),
          plot.caption = element_text(hjust = 0,
                                      size = caption_size+2),
          axis.text.x = element_text(angle = 90,
                                     hjust = 0.5,
                                     vjust = 0.5,
                                     colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.background = element_rect(fill = "white",
                                          colour = "gray",
                                          size = 0.15,
                                          linetype = "solid"),
          strip.background = element_rect(fill   = "transparent",
                                          colour = "transparent"),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      size = 0.7),
          legend.position = "bottom", # c(0.44,0.5),
          legend.text = element_text(size = 20),
          legend.justification = "center",
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.key = element_rect(fill   = "transparent",
                                    colour = "transparent"))
  
  # Plot incident deaths
  l_plots_fig3[[i+1]] <- ggplot(
    data = subset(df_fig3_CasesDeaths_interv, Outcome == "Incident deaths" &
                    proj_type == n_proj_type),
    aes(x = dates,
        y = value,
        colour = Intervention,
        shape = Intervention)) +
    geom_line(size = 0.9) + 
    facet_wrap(~BaseCase, ncol = 2) +
    scale_color_manual(values = rep(v_colors_colour, 2)) +
    scale_shape_manual(values = rep(v_colors_shape, 2)) +
    geom_vline(xintercept = c(as.numeric(date_NPI_update)),
               show.legend = TRUE,
               size = 0.5, 
               linetype = c("dashed"),
               # linetype = c("dashed", "twodash"),
               color = "#393F3E") +
    annotation_custom(NPI_update) +
    scale_x_date(label_x_axis,
                 date_labels = "%b/%d", #Change the format of the date in the x axis
                 breaks = "week") +
    scale_y_continuous("Deaths",#"Difference to case base scenario",
                       breaks = number_ticks(8),
                       # labels = scales::percent
                       labels = scales::comma) + # trans = "log"
    theme(text = element_text(size = text_size),
          plot.caption = element_text(hjust = 0,
                                      size = caption_size+2),
          axis.text.x = element_text(angle = 90,
                                     hjust = 0.5,
                                     vjust = 0.5,
                                     colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.background = element_rect(fill = "white",
                                          colour = "gray",
                                          size = 0.15,
                                          linetype = "solid"),
          strip.background = element_rect(fill   = "transparent",
                                          colour = "transparent"),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      size = 0.7),
          legend.position = "bottom", # c(0.44,0.5),
          legend.text = element_text(size = 20),
          legend.justification = "center",
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.key = element_rect(fill   = "transparent",
                                    colour = "transparent"))
  
  
  ## Combine graphs ---------------------------------------------------------
  
  plot_fig3 <- ggarrange(l_plots_fig3[[i]],
                         NULL,
                         l_plots_fig3[[i+1]],
                         ncol = 1,
                         common.legend = T,
                         legend = "bottom", heights = c(1, 0.05, 1),
                         labels = c("A","","B"), hjust = 0,vjust = 1.5,
                         font.label = list(size = 26, face = "bold"))
  
  # Save plot
  ggsave(paste0("figs/figs_SA/fig3_CasesDeaths_interv_",n_proj_type,".jpg"),
         width = 16, height = 18, dpi = 300)
  
  i <- i + 2
}


# FIGURE 4: interventions total hospitalizations --------------------------

# Interventions
select_intervention <- c("SchoolSD", 
                         "IncreaseSD",
                         "IncreaseSDSchoolSD",
                         "BaseCase")

# Offset flag
offset_flag <- TRUE

# Data.frame
df_hosp_all_INT_fig4 <-  df_out_total_prob_all %>%
  filter(Outcome == "Total hospitalizations" &
           intervention_type %in% select_intervention)

df_hosp_all_INT_fig4 <- df_hosp_all_INT_fig4[order(df_hosp_all_INT_fig4$dates, decreasing = F),]

# Observed data
df_hosp_fig4 <- df_hosp_MCMA %>%
  filter(dates <= n_date_last)
df_hosp_fig4$dates <- as.Date(df_hosp_fig4$dates)
df_hosp_fig4 <- df_hosp_fig4[order(df_hosp_fig4$dates,decreasing = F),]

df_hosp_fig4 <- df_hosp_fig4 %>%
  mutate(Intervention = "Observed") %>%
  filter(dates <= n_date_last)

df_hosp_fig4$Intervention <- as.factor(df_hosp_fig4$Intervention)

# Data.frame
df_fig4_HospTot_interv <- df_hosp_all_INT_fig4 %>%
  filter(county == state_i &
           intervention_type %in% select_intervention &
           proj_type == "SA") %>%
  group_by(county, type, BaseCase, BaseCase_type, Intervention, 
           intervention_type, Outcome, dates, proj_type, proj_type_label) %>%
  summarise(mean = round(mean(value),0),
            median = round(quantile(value, probs = 0.5, names = FALSE),0),
            sd = sd(value),
            lb = round(quantile(value, probs = 0.025, names = FALSE),0),
            ub = round(quantile(value, probs = 0.975, names = FALSE),0)) %>%
  rename(value = mean, 
         state = county
  ) %>%
  filter(dates >= n_date_last) %>%
  mutate(cap_tot = unique(df_hosp_observed$cap_tot[df_hosp_observed$dates==n_date_last]),
         time_stamp = n_time_stamp)

df_fig4_HospTot_interv$dates <- as.Date(df_fig4_HospTot_interv$dates)

# Set vector of interventions and labels
v_interv_names <- c(BaseCase           = "Policy A. Physical distancing: status quo; Schooling: not in-person",
                    IncreaseSD         = "Policy B. Physical distancing: +24% compared to status quo; Schooling: not in-person",
                    IncreaseSDSchoolSD = "Policy C. Physical distancing: +24% compared to status quo; Schooling: in-person",
                    SchoolSD           = "Policy D. Physical distancing: status quo; Schooling: in-person")

# Rename interventions
df_fig4_HospTot_interv$Intervention <- as.character(df_fig4_HospTot_interv$Intervention)

for(n_interv in names(v_interv_names)){ 
  interv_type <- v_interv_names[which(names(v_interv_names) == n_interv)]
  
  df_fig4_HospTot_interv$Intervention[df_fig4_HospTot_interv$intervention_type == n_interv] <- interv_type
}

df_fig4_HospTot_interv$Intervention <- factor(df_fig4_HospTot_interv$Intervention,
                                                  levels = v_interv_names,
                                                  ordered = TRUE)

# Save
if(offset_flag){
  save(df_fig4_HospTot_interv,
       file = paste0("figs/figs_SA/data_frames/df_fig4_HospTot_interv_offset.RData"))
}else{
  save(df_fig4_HospTot_interv,
       file = paste0("figs/figs_SA/data_frames/df_fig4_HospTot_interv.RData"))
}


## Plot specifications ----------------------------------------------------

# Colors
v_names_colours <- c(as.character(levels(df_fig4_HospTot_interv$Intervention)))
v_colors <- c("#0F57CA", "#008450", "#FF9900", "#B81D13")

# Vector of colors to colour
v_colors_colour <- v_colors
names(v_colors_colour) <- v_names_colours

# Vector of colors to fill
v_colors_fill <- v_colors
names(v_colors_fill) <- v_names_colours

# Vector of colors to shape
v_colors_shape <- c(15, 16, 17, 18)
names(v_colors_shape) <- v_names_colours

# Annotations
NPI_update <- grobTree(textGrob("Policy start date",
                                x = 0.14,
                                y = 0.97,
                                hjust = 0,
                                gp = gpar(col = "#393F3E",
                                          fontsize = annotation_size)))

## Plot -------------------------------------------------------------------


# List of plots
l_plot_fig4 <- vector(mode = "list", length = 2)
i <- 1

for(n_proj_type in c("SA")){
  
  # Plot
  l_plot_fig4[[i]] <- ggplot(
    data = subset(df_fig4_HospTot_interv, proj_type == n_proj_type),
    aes(x = dates,
        y = value,
        colour = Intervention)) +
    # geom_bar(stat = "identity", position = position_dodge(width = 5),
    #          width = 3) +
    geom_line(size = 0.9) + 
    # geom_point(size = 3) +
    geom_line(aes(x = dates,
                  y = cap_tot),
              size = 1.2,
              color = "#393F3E" #, linetype = "twodash"
    ) +
    facet_wrap(~BaseCase, ncol = 2) +
    scale_color_manual(values = rep(v_colors_colour, 2)) +
    # scale_shape_manual(values = rep(v_colors_shape, 2)) +
    geom_vline(xintercept = c(as.numeric(date_NPI_update)),
               show.legend = TRUE,
               size = 0.5, 
               linetype = c("dashed"),
               # linetype = c("dashed", "twodash"),
               color = "#393F3E") +
    annotation_custom(NPI_update) +
    annotate("text",
             label  = paste0("Current capacity: ",
                             unique(df_hosp_observed$cap_tot[df_hosp_observed$dates==n_date_last])," beds"),
             x      = as.Date("2020-12-27"),
             y      = unique(df_hosp_observed$cap_tot[df_hosp_observed$dates==n_date_last]) + 1000,
             size   = 6,
             colour = "#393F3E") +
    # scale_fill_manual(values = rep(v_colors_fill, 2)) +
    # scale_color_manual(values = v_colors_colour) +
    # scale_fill_manual(values = v_colors_fill) +
    # guides(fill = guide_legend(title = ""),
    #        color = guide_legend(title = "",
    #                             ncol = 1,
    #                             override.aes = list(color = v_colors_colour,
    #                                                 fill  = v_colors_fill))) +
    scale_x_date(label_x_axis,
                 date_labels = "%b/%d", #Change the format of the date in the x axis
                 breaks = "week") +
    scale_y_continuous("Total hospitalizations",
                       breaks = number_ticks(8),
                       # labels = scales::percent
                       labels = scales::comma) + # trans = "log"
    #  labs(# title    = "COVID-19 incident cases in MCMA, Mexico",
    # caption  = gg_caption
    #       ) +
    theme(text = element_text(size = text_size),
          plot.caption = element_text(hjust = 0,
                                      size = caption_size),
          axis.text.x = element_text(angle = 90,
                                     hjust = 0.5,
                                     vjust = 0.5,
                                     colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.background = element_rect(fill = "white",
                                          colour = "gray",
                                          size = 0.15,
                                          linetype = "solid"),
          strip.background = element_rect(fill   = "transparent",
                                          colour = "transparent"),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      size = 0.7),
          legend.position = "bottom", # c(0.44,0.5),
          legend.text = element_text(size = 20),
          legend.justification = "center",
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.key = element_rect(fill   = "transparent",
                                    colour = "transparent"))
  
  # Save plot
  if(offset_flag){
    ggsave(paste0("figs/figs_SA/fig4_HospTot_interv_offset_",n_proj_type,".jpg"),
           width = 16, height = 9, dpi = 300)
  }else{
    ggsave(paste0("figs/figs_SA/fig4_HospTot_interv_",n_proj_type,".jpg"),
           width = 16, height = 9, dpi = 300)
  }
  
  i <- i + 1
  
}


# FIGURE 5: exceed hosp capacity ------------------------------------------

# Interventions
select_intervention <- c("SchoolSD", 
                         "IncreaseSD",
                         "IncreaseSDSchoolSD",
                         "BaseCase")

# Vector of outcomes of interest
outcome <- "Total hospitalizations"

# Hospital capacity
v_hosp_capacity <- c(Hosp_tot   = 9667,
                     NonICU_tot = 7008,
                     ICU_tot    = 2659)

# Interventions data
df_fig5_ExceedHospCap <- df_out_total_prob_all %>%
  filter(Outcome == outcome &
           intervention_type %in% select_intervention &
           dates >= n_date_last  &
           proj_type == "SA") %>%
  mutate(cap_tot = df_hosp_MCMA$cap_tot[df_hosp_MCMA$dates == n_date_last]) %>%
  mutate(exceed_cap = ifelse(value > cap_tot, 1, 0)) %>%
  group_by(county, type, BaseCase, BaseCase_type, Intervention, intervention_type, 
           Outcome, dates, proj_type, proj_type_label) %>%
  mutate(exceed_cap_count = sum(exceed_cap),
         exceed_cap_prop  = exceed_cap_count*100/1000) %>%
  summarise(value  = mean(exceed_cap_prop),
            median = quantile(exceed_cap_prop, probs = 0.5, names = FALSE),
            sd = sd(exceed_cap_prop),
            lb = quantile(exceed_cap_prop, probs = 0.025, names = FALSE),
            ub = quantile(exceed_cap_prop, probs = 0.975, names = FALSE)) %>%
  rename(state = county)


# Rename intervention_label
df_fig5_ExceedHospCap$intervention_label <- ""
for(n_interv in names(v_interv_labels)){ 
  df_fig5_ExceedHospCap$intervention_label[df_fig5_ExceedHospCap$intervention_type == n_interv] <- v_interv_labels[n_interv]
}

df_fig5_ExceedHospCap$intervention_label <- factor(df_fig5_ExceedHospCap$intervention_label,
                                                   levels = v_interv_labels,
                                                   ordered = TRUE)

# Set vector of interventions and labels
v_interv_names <- c(BaseCase           = "Policy A. Physical distancing: status quo; Schooling: not in-person",
                    IncreaseSD         = "Policy B. Physical distancing: +24% compared to status quo; Schooling: not in-person",
                    IncreaseSDSchoolSD = "Policy C. Physical distancing: +24% compared to status quo; Schooling: in-person",
                    SchoolSD           = "Policy D. Physical distancing: status quo; Schooling: in-person")

# Rename interventions
df_fig5_ExceedHospCap$Intervention <- as.character(df_fig5_ExceedHospCap$Intervention)

for(n_interv in names(v_interv_names)){ 
  interv_type <- v_interv_names[which(names(v_interv_names) == n_interv)]
  
  df_fig5_ExceedHospCap$Intervention[df_fig5_ExceedHospCap$intervention_type == n_interv] <- interv_type
}

df_fig5_ExceedHospCap$Intervention <- factor(df_fig5_ExceedHospCap$Intervention,
                                              levels = v_interv_names,
                                              ordered = TRUE)

# Save
if(offset_flag){
  save(df_fig5_ExceedHospCap,
       file = paste0("figs/figs_SA/data_frames/df_fig5_ExceedHospCap_",abbrev_outcome,"_offset.RData"))
}else{
  save(df_fig5_ExceedHospCap,
       file = paste0("figs/figs_SA/data_frames/df_fig5_ExceedHospCap_",abbrev_outcome,".RData"))
}

# Plot specifications -----------------------------------------------------

# Caption

# Set vector of interventions and labels
v_interv_labels <- c(BaseCase           = "Policy A",
                     IncreaseSD         = "Policy B",
                     IncreaseSDSchoolSD = "Policy C",
                     SchoolSD           = "Policy D" 
)

plot_caption <- ""
for(intv in unique(df_fig5_ExceedHospCap$Intervention)){
  plot_caption <- paste0(plot_caption, paste0(intv,"\n"))
}

## Plot -------------------------------------------------------------------

# List of plots
l_plots_fig5 <- vector(mode = "list", length = 2)
i <- 1

for(n_proj_type in c("SA")){
  
  # Abbrev outcome
  abbrev_outcome <- names(v_outcomes)[which(v_outcomes == outcome)]
  
  # Plot
  l_plots_fig5[[i]] <- ggplot(data = subset(df_fig5_ExceedHospCap, proj_type == n_proj_type), 
                              aes(x    = dates, 
                                  y    = intervention_label, 
                                  fill = value)) +
    geom_raster() +
    labs(y = "", fill = "Probability of\nexceeding hospital\ncapacity (%)") +
    scale_fill_gradientn(colours = color_map,
                         breaks = c(0.0, 25.0, 50.0, 75.0, 100.0),
                         labels = c("0", "25", "50", "75", "100"),
                         limits = c(0,100)) +  # Hawre's jet map
    facet_wrap( ~ BaseCase, ncol = 1) +
    scale_x_date(label_x_axis,
                 date_labels = "%b/%d", #Change the format of the date in the x axis
                 # breaks = seq.Date(from = first(df_hosp_all_INT_daily$week),
                 #                   to = last(df_hosp_all_INT_daily$week),
                 #                   by = "week")
                 breaks = "week") + 
    labs(caption = plot_caption,
         title = outcome) +
    theme(text = element_text(size = text_size),
          plot.caption = element_text(hjust = 0,
                                      size = 19), 
          axis.text.x = element_text(angle = 90, 
                                     hjust = 0.5, 
                                     vjust = 0.5,
                                     colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.background = element_rect(fill = "white", 
                                          colour = "gray", 
                                          size = 0.15, 
                                          linetype = "solid"),
          panel.border = element_rect(colour = "black", 
                                      fill = NA, 
                                      size = 0.7), 
          strip.background = element_rect(fill   = "transparent",
                                          colour = "transparent"),
          legend.position = "right", 
          legend.title = element_text(size = 18),
          legend.justification = "center", 
          legend.direction = "vertical", 
          legend.key.width = unit(0.9, "cm"),
          legend.key.height = unit(1.2, "cm"),
          legend.title.align = 0,
          legend.key = element_rect(fill   = "transparent", 
                                    colour = "transparent"))
  
  if(offset_flag){
    ggsave(paste0("figs/figs_SA/fig5_ExceedHospCap_",abbrev_outcome,"_",n_proj_type,"_offset.jpg"),
           width = 16, height = 11, dpi = 300)
  }else{
    ggsave(paste0("figs/figs_SA/fig5_ExceedHospCap_",abbrev_outcome,"_",n_proj_type,".jpg"),
           width = 16, height = 11, dpi = 300)
  }
  
  i <- i + 1
}


# Figure S4: Effective Reproduction Number --------------------------------

## Status quo -------------------------------------------------------------

# Outcome
outcome     <- "R effective"
outcome_eng <- "R effective"

# Data.frame
df_outcome_Rt <- df_out_total_prob %>%
  ungroup() %>%
  filter(state == state_i &
           Outcome == outcome_eng &
           dates <= max_date & 
           proj_type == "SA")

# Set vector of interventions and labels
v_interv_names <- c(BaseCase           = "Policy A. Physical distancing: status quo; Schooling: not in-person",
                    IncreaseSD         = "Policy B. Physical distancing: +24% compared to status quo; Schooling: not in-person",
                    IncreaseSDSchoolSD = "Policy C. Physical distancing: +24% compared to status quo; Schooling: in-person",
                    SchoolSD           = "Policy D. Physical distancing: status quo; Schooling: in-person")

# Rename interventions
df_outcome_Rt$Intervention <- as.character(df_outcome_Rt$Intervention)

for(n_interv in names(v_interv_names)){ 
  interv_type <- v_interv_names[which(names(v_interv_names) == n_interv)]
  
  df_outcome_Rt$Intervention[df_outcome_Rt$intervention_type == n_interv] <- interv_type
}

df_outcome_Rt$Intervention <- factor(df_outcome_Rt$Intervention,
                                     levels = v_interv_names,
                                     ordered = TRUE)

# Create data.frame for BaseCase
df_Rt_BaseCase <- subset(df_outcome_Rt, intervention_type == "BaseCase")
df_Rt_BaseCase$Intervention <- "Model-projected"

# Create data.frame for Interventions
df_RTt_INT <- subset(df_outcome_Rt, intervention_type %in% c("SchoolSD", "IncreaseSD","IncreaseSDSchoolSD"))

# Save data.frame
df_figS4_Rt <- rbind.data.frame(df_Rt_BaseCase,
                                df_RTt_INT)
save(df_figS4_Rt,
     file = paste0("figs/figs_SA/data_frames/df_figS4_Rt.RData"))

# List of plots
l_plots_S4 <- vector(mode = "list", length = 4)
i <- 1

for(n_proj_type in c("SA")){
  
  ### Plot: SQ --------------------------------------------------------------
  
  # Abbreviation for outcome to save the plot
  abbrev_outcome <- names(v_outcomes)[which(v_outcomes == outcome)]
  
  # Name of colors
  v_names_colours <- as.character(levels(df_Rt_BaseCase$BaseCase))
  
  # Colors
  v_colors_colour <- c("#659583", "#a13023" )
  names(v_colors_colour) <- v_names_colours
  
  # Vector of colors to fill
  v_colors_fill <- v_colors_colour
  names(v_colors_fill) <- v_names_colours
  
  # Line
  v_linetype <- c("solid", "solid")
  names(v_linetype) <- v_names_colours
  
  # Bracket position
  y_end <- max(df_Rt_BaseCase$ub) 
  
  # For annotation_custom
  npi_date <- grobTree(textGrob("Date of NPI implementation", 
                                x = 0.115, 
                                y = 0.97, 
                                hjust = 0,
                                gp = gpar(col = "#393F3E", 
                                          fontsize = annotation_size)))
  
  calib_day <- grobTree(textGrob("Last day used for calibration", 
                                 x = 0.675, 
                                 y = 0.97, 
                                 hjust = 0,
                                 gp = gpar(col = "#393F3E", 
                                           fontsize = annotation_size)))
  
  # Plot
  l_plots_S4[[i]] <- ggplot(data = subset(df_Rt_BaseCase, proj_type == n_proj_type)) +
    geom_line(aes(x = dates, 
                  y = value, 
                  colour = BaseCase),
              size = 0.9) +
    geom_ribbon(aes(x = dates,
                    y = value,
                    ymin = lb,
                    ymax = ub,
                    fill = BaseCase),
                alpha = 0.4,
                colour = NA) +
    scale_color_manual(values = v_colors_colour) + 
    scale_fill_manual(values =  v_colors_fill)  +
    guides(fill = guide_legend(title = ""),
           color = guide_legend(title = "", 
                                ncol = 1,
                                override.aes = list(#shape    = v_shape, 
                                  linetype = v_linetype,
                                  color    = v_colors_colour,
                                  fill     = v_colors_fill))) +
    geom_vline(xintercept = c(as.numeric(as.Date("2020-03-23")), as.numeric(n_date_last)),
               show.legend = TRUE,
               size = 0.5, 
               linetype = c("dashed", "twodash"),
               color = "#393F3E") +
    geom_hline(yintercept = 1, linetype = 3) +
    annotation_custom(npi_date) +
    scale_x_date(label_x_axis,
                 date_labels = "%b/%d", #Change the format of the date in the x axis
                 breaks = v_dates_breaks
    ) + 
    scale_y_continuous(outcome_eng, 
                       breaks = number_ticks(8)#,
                       #labels = scales::comma
    ) + # trans = "log"
    # labs(title    = "COVID-19 effective reproduction number in MCMA, Mexico" #, caption  = gg_caption
    # ) +
    theme(text = element_text(size = text_size),
          # plot.caption = element_text(hjust = 0,
          #                             size = caption_size), 
          axis.text.x = element_text(angle = angle_x_axis, 
                                     hjust = 0.5, 
                                     vjust = 0.5,
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
          legend.justification = "center", 
          legend.direction = "vertical", 
          legend.title = element_blank(), 
          legend.key = element_rect(fill   = "transparent", 
                                    colour = "transparent")) 
  
  
  ## Interventions ----------------------------------------------------------
  
  # Abbreviation for outcome to save the plot
  abbrev_outcome <- names(v_outcomes)[which(v_outcomes == outcome)]
  
  # Colors
  v_names_colours <- c(as.character(levels(df_RTt_INT$Intervention)[-1]))
  
  # Colours
  v_colors_colour <- c(#"#659583",  
    "#008450", "#FF9900", "#B81D13")
  names(v_colors_colour) <- v_names_colours
  
  # Vector of colors to fill
  v_colors_fill <- v_colors_colour
  names(v_colors_fill) <- v_names_colours
  
  # # Vector of colors to shape
  # v_colors_shape <- c(#8, 
  #   15, 16, 17, 18)
  # names(v_colors_shape) <- v_names_colours
  
  # End
  df_RTt_INT_temp <- subset(df_RTt_INT, dates >= n_date_last & proj_type == n_proj_type)
  y_end <- max(df_RTt_INT_temp$ub)
  y_start <- min(df_RTt_INT_temp$lb)
  
  # For annotation_custom
  scenarios_start <- grobTree(textGrob("Policy start date", 
                                       x = 0.145, 
                                       y = 0.97, 
                                       hjust = 0,
                                       gp = gpar(col = "#393F3E", 
                                                 fontsize = annotation_size)))
  
  # calib_day <- grobTree(textGrob("Last day used for calibration", 
  #                                x = 0.675, 
  #                                y = 0.97, 
  #                                hjust = 0,
  #                                gp = gpar(col = "#393F3E", 
  #                                          fontsize = annotation_size)))
  
  
  ### Plot: Interventions ---------------------------------------------------
  
  l_plots_S4[[i+1]] <- ggplot(data = df_RTt_INT_temp) + 
    geom_line(aes(x = dates,
                  y = value,
                  colour = Intervention),
              size = 0.9) +
    geom_ribbon(aes(x = dates,
                    y = value,
                    ymin = lb,
                    ymax = ub,
                    fill = Intervention),
                alpha = 0.4,
                colour = NA) +
    facet_wrap(~BaseCase, ncol = 2) +
    scale_color_manual(values = rep(v_colors_colour, 2)) +
    scale_fill_manual(values =  rep(v_colors_fill, 2))  +
    guides(fill = guide_legend(title = ""),
           color = guide_legend(title = "", 
                                ncol = 1,
                                override.aes = list(color = v_colors_colour,
                                                    fill  = v_colors_fill #, shape = v_colors_shape
                                ))) +
    geom_vline(xintercept = as.numeric(date_NPI_update),
               show.legend = TRUE,
               size = 0.5, 
               linetype = "dashed",
               color = "#393F3E") +
    geom_hline(yintercept = 1, linetype = 3) +
    annotation_custom(scenarios_start) +
    scale_x_date(label_x_axis,
                 date_labels = "%b/%d", #Change the format of the date in the x axis
                 breaks = v_dates_breaks[20:26]
    ) + 
    scale_y_continuous(outcome_eng, 
                       breaks = number_ticks(8)#,
                       #labels = scales::comma
    ) + # trans = "log"
    # labs(title    = "COVID-19 effective reproduction number in MCMA, Mexico" #, caption  = gg_caption
    # ) +
    theme(text = element_text(size = text_size),
          # plot.caption = element_text(hjust = 0,
          #                             size = caption_size), 
          axis.text.x = element_text(angle = 0, 
                                     hjust = 0.5, 
                                     vjust = 0.5,
                                     colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.background = element_rect(fill = "white", 
                                          colour = "gray", 
                                          size = 0.15, 
                                          linetype = "solid"),
          panel.border = element_rect(colour = "black", 
                                      fill = NA, 
                                      size = 0.7), 
          strip.background = element_rect(fill   = "transparent",
                                          colour = "transparent"),
          legend.position = "bottom",
          legend.text = element_text(size = 16),
          legend.justification = "center", 
          legend.direction = "vertical", 
          legend.title = element_blank(), 
          legend.key = element_rect(fill   = "transparent", 
                                    colour = "transparent")) +
    coord_cartesian(ylim = c(y_start*0.90, y_end*1.1))
  
  
  ## Combine two graphs -----------------------------------------------------
  
  plot_figS4 <- ggarrange(l_plots_S4[[i]],
                          NULL,
                          l_plots_S4[[i+1]],
                          ncol = 1,
                          heights = c(1, 0.1, 1),
                          labels = c("A","","B"), hjust = 0,vjust = 1.5,label.x = c(0,0,0), label.y = c(1,1,1),
                          font.label = list(size = 26, face = "bold"))
  
  ggsave(paste0("figs/figs_SA/figS4_Rt_BaseCase_INT_",n_proj_type,".jpg"), 
         width = 16, height = 18, dpi = 300)
  
  i <- i + 2
}

# FIGURE S6: Cumulative proportion of infections being detected -----------

# Compute proportions
df_outcome_all_detec_inf <- df_out_total_prob_all %>%
  ungroup() %>%
  filter(county == state_i &
           Outcome == "Cumulative cases" & 
           BaseCase_type == "Holidays" &
           intervention_type == "BaseCase" & 
           proj_type == "SA")

df_outcome_all_total_inf <- df_out_total_prob_all %>%
  ungroup() %>%
  filter(county == state_i &
           Outcome == "Cumulative infections" & 
           BaseCase_type == "Holidays" &
           intervention_type == "BaseCase" & 
           proj_type == "SA")

df_outcome_all_detec_inf$Intervention <- "Model-projected"
df_outcome_all_total_inf$Intervention <- "Model-projected"

# Compute proportion
df_prop_DetecInf <- data.frame(df_outcome_all_total_inf[,c("type","BaseCase_type", "Intervention",
                                                           "intervention_type", "Outcome", "dates",
                                                           "proj_type", "proj_type_label")],
                               value = df_outcome_all_detec_inf$value/df_outcome_all_total_inf$value
)

df_figS6_DetecInf_prop <- df_prop_DetecInf %>%
  group_by(type, BaseCase_type, Intervention, intervention_type, Outcome, dates,
           proj_type, proj_type_label) %>%
  summarise(mean = mean(value),
            median = quantile(value, probs = 0.5, names = FALSE),
            sd = sd(value),
            lb = quantile(value, probs = 0.025, names = FALSE),
            ub = quantile(value, probs = 0.975, names = FALSE)) %>%
  rename(value = mean)

# Intervention as factor
df_figS6_DetecInf_prop$Intervention <- as.factor(df_figS6_DetecInf_prop$Intervention)

# Save
save(df_figS6_DetecInf_prop,
     file = paste0("figs/figs_SA/data_frames/df_figS6_DetecInf_prop.RData"))

# Abbreviation for outcome to save the plot
abbrev_outcome <- "DetecInf_prop"

# Vector of colors
v_names_colours <- as.character(levels(df_figS6_DetecInf_prop$Intervention))
v_colors_colour <- c("#659583")
names(v_colors_colour) <- v_names_colours

# Vector of colors to fill
v_colors_fill <- c("#659583")
names(v_colors_fill) <- v_names_colours

# Line
v_linetype <- c("solid")
names(v_linetype) <- v_names_colours

# Bracket position
y_end <- max(df_figS6_DetecInf_prop$ub)


# For annotation_custom
npi_date <- grobTree(textGrob("Date of NPI implementation",
                              x = 0.119,
                              y = 0.97,
                              hjust = 0,
                              gp = gpar(col = "#393F3E",
                                        fontsize = annotation_size)))

calib_day <- grobTree(textGrob("Last day used for calibration",
                               x = 0.675,
                               y = 0.97,
                               hjust = 0,
                               gp = gpar(col = "#393F3E",
                                         fontsize = annotation_size)))

## Plot -------------------------------------------------------------------
for(n_proj_type in c("SA")){
  
  ggplot(data = subset(df_figS6_DetecInf_prop, proj_type == n_proj_type)) +
    geom_line(aes(x = dates,
                  y = value,
                  colour = Intervention),
              size = 0.9) +
    geom_ribbon(aes(x = dates,
                    y = value,
                    ymin = lb,
                    ymax = ub,
                    fill = Intervention),
                alpha = 0.4,
                colour = NA) +
    scale_color_manual("", values = v_colors_colour) +
    scale_fill_manual("", values =  v_colors_fill) +
    geom_vline(xintercept = c(as.numeric(as.Date("2020-03-23")), as.numeric(n_date_last)),
               show.legend = TRUE,
               size = 0.5,
               linetype = c("dashed", "twodash"),
               color = "#393F3E") +
    annotation_custom(npi_date) +
    # annotation_custom(calib_day) +
    scale_x_date(label_x_axis,
                 date_labels = "%b/%d", #Change the format of the date in the x axis
                 breaks = v_dates_breaks
    ) +
    scale_y_continuous("",
                       breaks = number_ticks(8),
                       labels = scales::percent) + # trans = "log"
    # labs(title    = "COVID-19 total infections in MCMA, Mexico" #, caption  = gg_caption
    # ) +
    theme(text = element_text(size = text_size),
          axis.text.x = element_text(angle = angle_x_axis,
                                     hjust = 0.5,
                                     vjust = 0.5,
                                     colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.background = element_rect(fill = "white", 
                                          colour = "gray", 
                                          size = 0.15, 
                                          linetype = "solid"),
          panel.border = element_rect(colour = "black", 
                                      fill = NA, 
                                      size = 0.7), 
          strip.background = element_rect(fill   = "transparent",
                                          colour = "transparent"),
          legend.position = "none",
          legend.justification = "center",
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.key = element_rect(fill   = "transparent",
                                    colour = "transparent")) 
  
  # Save plot
  ggsave(paste0("figs/figs_SA/figS6_DetecInf_prop_",n_proj_type,".jpg"),
         width = 12, height = 9, dpi = 300)
}


# ## Plot: facet ------------------------------------------------------------
# 
# ggplot(data = df_figS6_DetecInf_prop) +
#   geom_line(aes(x = dates,
#                 y = value,
#                 colour = Intervention),
#             size = 0.9) +
#   geom_ribbon(aes(x = dates,
#                   y = value,
#                   ymin = lb,
#                   ymax = ub,
#                   fill = Intervention),
#               alpha = 0.4,
#               colour = NA) +
#   facet_wrap(~proj_type_label) +
#   scale_color_manual("", values = rep(v_colors_colour, 2)) +
#   scale_fill_manual("", values =  rep(v_colors_fill,2)) +
#   geom_vline(xintercept = c(as.numeric(as.Date("2020-03-23")), as.numeric(n_date_last)),
#              show.legend = TRUE,
#              size = 0.5,
#              linetype = rep(c("dashed", "twodash"),2),
#              color = "#393F3E") +
#   annotation_custom(npi_date) +
#   # annotation_custom(calib_day) +
#   scale_x_date(label_x_axis,
#                date_labels = "%b/%d", #Change the format of the date in the x axis
#                breaks = v_dates_breaks
#   ) +
#   scale_y_continuous("",
#                      breaks = number_ticks(8),
#                      labels = scales::percent) + # trans = "log"
#   # labs(title    = "COVID-19 total infections in MCMA, Mexico" #, caption  = gg_caption
#   # ) +
#   theme(text = element_text(size = text_size),
#         axis.text.x = element_text(angle = angle_x_axis,
#                                    hjust = 0.5,
#                                    vjust = 0.5,
#                                    colour = "black"),
#         axis.text.y = element_text(colour = "black"),
#         panel.background = element_rect(fill = "white", 
#                                         colour = "gray", 
#                                         size = 0.15, 
#                                         linetype = "solid"),
#         panel.border = element_rect(colour = "black", 
#                                     fill = NA, 
#                                     size = 0.7), 
#         strip.background = element_rect(fill   = "transparent",
#                                         colour = "transparent"),
#         legend.position = "none",
#         legend.justification = "center",
#         legend.direction = "vertical",
#         legend.title = element_blank(),
#         legend.key = element_rect(fill   = "transparent",
#                                   colour = "transparent")) 
# 
# # Save plot
# ggsave(paste0("figs/figs_SA/figS6_DetecInf_prop_SA_SQ.jpg"),
#        width = 18, height = 9, dpi = 300)



# FIGURE S5: Cumulative proportion of population ever been infected -------

# Outcome
outcome <- outcome_eng <- "Recovered prevalence proportion"

# Data.frame
df_figS5_Rec_prop <- df_out_total_prob %>%
  ungroup() %>%
  filter(state == state_i &
           Outcome == outcome_eng & BaseCase_type == "Holidays" &
           intervention_type == "BaseCase" &
           dates <= max_date & 
           proj_type == "SA")

# Set type as factor
df_figS5_Rec_prop$type <- as.factor(df_figS5_Rec_prop$type)

# Save
save(df_figS5_Rec_prop,
     file = paste0("figs/figs_SA/data_frames/df_figS5_Rec_prop.RData"))


# Plot specifications -----------------------------------------------------

# Abbreviation for outcome to save the plot
abbrev_outcome <- names(v_outcomes)[which(v_outcomes == outcome)]

# Name of colors
v_names_colours <- as.character(levels(df_figS5_Rec_prop$type))

# Colours
v_colors_colour <- c("#a13023")
names(v_colors_colour) <- v_names_colours

# Vector of colors to fill
v_colors_fill <- c("#a13023")
names(v_colors_fill) <- v_names_colours

# Line
v_linetype <- c("solid")
names(v_linetype) <- v_names_colours

# Bracket position
y_end <- max(df_figS5_Rec_prop$ub) 


# For annotation_custom
npi_date <- grobTree(textGrob("Date of NPI implementation", 
                              x = 0.117, 
                              y = 0.97, 
                              hjust = 0,
                              gp = gpar(col = "#393F3E", 
                                        fontsize = annotation_size)))

calib_day <- grobTree(textGrob("Last day used for calibration", 
                               x = 0.744, 
                               y = 0.97, 
                               hjust = 0,
                               gp = gpar(col = "#393F3E", 
                                         fontsize = annotation_size)))

## Plot --------------------------------------------------------------

for(n_proj_type in c("SA")){
  
  ggplot(data = subset(df_figS5_Rec_prop, proj_type == n_proj_type)) +
    geom_line(aes(x = dates, 
                  y = value, 
                  colour = type),
              size = 0.9) +
    geom_ribbon(aes(x = dates,
                    y = value,
                    ymin = lb,
                    ymax = ub,
                    fill = type),
                alpha = 0.4,
                colour = NA) +
    scale_color_manual("", values = v_colors_colour) + 
    scale_fill_manual("", values =  v_colors_fill)  +
    geom_vline(xintercept = c(as.numeric(as.Date("2020-03-23")), as.numeric(n_date_last)),
               show.legend = TRUE,
               size = 0.5, 
               linetype = c("dashed", "twodash"),
               color = "#393F3E") +
    annotation_custom(npi_date) +
    # annotation_custom(calib_day) +
    scale_x_date(label_x_axis,
                 date_labels = "%b/%d", #Change the format of the date in the x axis
                 breaks = v_dates_breaks
    ) + 
    scale_y_continuous("",#outcome_eng, 
                       breaks = number_ticks(8),
                       labels = scales::percent) + # trans = "log"
    # labs(title    = "COVID-19 total recovered in MCMA, Mexico" #, caption  = gg_caption
    # ) +
    theme(text = element_text(size = text_size),
          # plot.caption = element_text(hjust = 0,
          #                             size = caption_size), 
          axis.text.x = element_text(angle = angle_x_axis, 
                                     hjust = 0.5,
                                     vjust = 0.5,
                                     colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.background = element_rect(fill = "white", 
                                          colour = "gray", 
                                          size = 0.15, 
                                          linetype = "solid"),
          panel.border = element_rect(colour = "black", 
                                      fill = NA, 
                                      size = 0.7), 
          strip.background = element_rect(fill   = "transparent",
                                          colour = "transparent"),
          legend.position = "none", 
          legend.justification = "center", 
          legend.direction = "vertical", 
          legend.title = element_blank(), 
          legend.key = element_rect(fill   = "transparent", 
                                    colour = "transparent")) 
  
  # Save plot
  ggsave(paste0("figs/figs_SA/figS5_Rec_prop_",n_proj_type,".jpg"), 
         width = 12, height = 9, dpi = 300)
}


# # Plot: facets ------------------------------------------------------------
# 
# ggplot(data = df_figS5_Rec_prop) +
#   geom_line(aes(x = dates, 
#                 y = value, 
#                 colour = type),
#             size = 0.9) +
#   geom_ribbon(aes(x = dates,
#                   y = value,
#                   ymin = lb,
#                   ymax = ub,
#                   fill = type),
#               alpha = 0.4,
#               colour = NA) +
#   facet_wrap(~ proj_type_label) +
#   scale_color_manual("", values = rep(v_colors_colour, 2)) + 
#   scale_fill_manual("", values =  rep(v_colors_fill, 2))  +
#   geom_vline(xintercept = c(as.numeric(as.Date("2020-03-23")), as.numeric(n_date_last)),
#              show.legend = TRUE,
#              size = 0.5, 
#              linetype = rep(c("dashed", "twodash"), 2),
#              color = "#393F3E") +
#   annotation_custom(npi_date) +
#   # annotation_custom(calib_day) +
#   scale_x_date(label_x_axis,
#                date_labels = "%b/%d", #Change the format of the date in the x axis
#                breaks = v_dates_breaks
#   ) + 
#   scale_y_continuous("",#outcome_eng, 
#                      breaks = number_ticks(8),
#                      labels = scales::percent) + # trans = "log"
#   # labs(title    = "COVID-19 total recovered in MCMA, Mexico" #, caption  = gg_caption
#   # ) +
#   theme(text = element_text(size = text_size),
#         # plot.caption = element_text(hjust = 0,
#         #                             size = caption_size), 
#         axis.text.x = element_text(angle = angle_x_axis, 
#                                    hjust = 0.5,
#                                    vjust = 0.5,
#                                    colour = "black"),
#         axis.text.y = element_text(colour = "black"),
#         panel.background = element_rect(fill = "white", 
#                                         colour = "gray", 
#                                         size = 0.15, 
#                                         linetype = "solid"),
#         panel.border = element_rect(colour = "black", 
#                                     fill = NA, 
#                                     size = 0.7), 
#         strip.background = element_rect(fill   = "transparent",
#                                         colour = "transparent"),
#         legend.position = "none", 
#         legend.justification = "center", 
#         legend.direction = "vertical", 
#         legend.title = element_blank(), 
#         legend.key = element_rect(fill   = "transparent", 
#                                   colour = "transparent")) 
# 
# # Save plot
# ggsave(paste0("figs/figs_SA/figS5_Rec_prop_SA_SQ.jpg"), 
#        width = 18, height = 9, dpi = 300)



# FIGURE S7: Percentage difference STATUS QUO VS HOLIDAYS -----------------

# Outcomes to be plotted
v_outcome     <-  c("Incident cases", 
                    "Incident deaths",
                    "Cumulative cases", 
                    "Cumulative deaths")

# Filter data and compute weekly data
df_CasesDeaths_daily <- df_out_total_prob_all %>%
  ungroup() %>%
  filter(county == state_i &
           Outcome %in% v_outcome &
           intervention_type == "BaseCase" &
           dates > n_date_last & 
           proj_type == "SA")

# Replicate status quo value
df_CasesDeaths_daily$BaseCase_value <- 0 
df_CasesDeaths_daily$BaseCase_value[df_CasesDeaths_daily$BaseCase_type == "Holidays"] <-  df_CasesDeaths_daily$value[df_CasesDeaths_daily$BaseCase_type == "StatusQuo"]
df_CasesDeaths_daily$BaseCase_value[df_CasesDeaths_daily$BaseCase_type == "StatusQuo"] <-  df_CasesDeaths_daily$value[df_CasesDeaths_daily$BaseCase_type == "StatusQuo"]

# Compute weekly percentage increases
df_figS7_PercDiff <- df_CasesDeaths_daily %>%
  mutate(diffPerc = (value - BaseCase_value)/BaseCase_value) %>%
  group_by(type, BaseCase_type, Intervention, intervention_type, Outcome, 
           dates, proj_type, proj_type_label) %>%
  summarise(value  = mean(diffPerc),
            median = quantile(diffPerc, probs = 0.5, names = FALSE),
            sd     = sd(diffPerc),
            lb     = quantile(diffPerc, probs = 0.025, names = FALSE),
            ub     = quantile(diffPerc, probs = 0.975, names = FALSE)) %>%
  filter(BaseCase_type == "Holidays")

# Save
save(df_figS7_PercDiff,
     file = paste0("figs/figs_SA/data_frames/df_figS7_PercDiff.RData"))

# Order outcomes
df_figS7_PercDiff$Outcome <- factor(df_figS7_PercDiff$Outcome,
                                    levels = v_outcome,
                                    ordered = TRUE)

# Colors
v_names_colours <- c(as.character(levels(df_figS7_PercDiff$Outcome)))

# Vector of colors to colour
v_colors_colour <- rep(c("#659583","#a13023"),2)
names(v_colors_colour) <- v_names_colours

# Vector of colors to fill
v_colors_fill <- v_colors_colour
names(v_colors_fill) <- v_names_colours

# Vector of colors to shape
v_colors_shape <- c(15, 16, 17,18)
names(v_colors_shape) <- v_names_colours

# For annotation_custom
start_holidays <- grobTree(textGrob("Start holidays", 
                                    x = 0.12, 
                                    y = 0.97, 
                                    hjust = 0,
                                    gp = gpar(col = "#393F3E", 
                                              fontsize = annotation_size)))

end_holidays <- grobTree(textGrob("End holidays", 
                                  x = 0.756, 
                                  y = 0.97, 
                                  hjust = 0,
                                  gp = gpar(col = "#393F3E", 
                                            fontsize = annotation_size)))

## Plot -------------------------------------------------------------------

for(n_proj_type in c("SA")){
  
  ggplot(data = subset(df_figS7_PercDiff, proj_type == n_proj_type),
         aes(x = dates,
             y = value,
             colour = Outcome,
             ymin = lb,
             ymax = ub,
             shape = Outcome)) +
    geom_point(#position = position_dodge(width = 5),
      size = 4) +
    geom_errorbar(#position = position_dodge(width = 5),
      size = 1.1) +
    facet_wrap(~ Outcome, ncol = 2, scales = "free_y") +
    # geom_ribbon(aes(x = dates,
    #                 y = mean_percent,
    #                 ymin = lb_percent,
    #                 ymax = ub_percent,
    #                 fill = Intervention),
    #             alpha = 0.4,
    #             colour = NA) +
    scale_color_manual(values = v_colors_colour) +
    scale_fill_manual(values =  v_colors_fill)  +
    scale_shape_manual(values = v_colors_shape) +
    guides(shape = guide_legend(title = ""),
           color = guide_legend(title = "",
                                ncol = 1,
                                override.aes = list(#linetype = v_linetype,
                                  color    = v_colors_colour,
                                  # fill     = v_colors_fill,
                                  shape    = v_colors_shape))) +
    geom_vline(xintercept = c(as.Date("2020-12-24"),
                              as.Date("2021-01-06")),
               show.legend = TRUE,
               size = 0.5,
               linetype = rep(c("dashed", "twodash"),4),
               color = "#393F3E") +
    # annotation_custom(start_holidays) +
    # annotation_custom(end_holidays) +
    scale_x_date(label_x_axis,
                 date_labels = "%b/%d", #Change the format of the date in the x axis
                 breaks = c(as.Date("2020-12-07"), as.Date("2020-12-14"),
                            as.Date("2020-12-21"), as.Date("2020-12-28"),
                            as.Date("2021-01-04"), as.Date("2021-01-11"),
                            as.Date("2021-01-18"), as.Date("2021-01-25"),
                            as.Date("2021-02-01"), as.Date("2021-02-08"),
                            as.Date("2021-02-15"), as.Date("2021-02-22"),
                            as.Date("2021-03-01"))
                 # breaks = "week"
    ) +
    scale_y_continuous("Difference from status quo (%)",
                       breaks = number_ticks(8),
                       labels = scales::percent
    ) + # trans = "log"
    #  labs(# title    = "COVID-19 incident cases in MCMA, Mexico",
    # caption  = gg_caption
    #       ) +
    theme(text = element_text(size = text_size),
          plot.caption = element_text(hjust = 0,
                                      size = caption_size),
          axis.text.x = element_text(angle = 90,
                                     hjust = 0.5,
                                     vjust = 0.5,
                                     colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.background = element_rect(fill = "white",
                                          colour = "gray",
                                          size = 0.15,
                                          linetype = "solid"),
          strip.background = element_rect(fill   = "transparent",
                                          colour = "transparent"),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      size = 0.7),
          legend.position = "none", # c(0.44,0.5),
          legend.text = element_text(size = 20),
          legend.justification = "center",
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.key = element_rect(fill   = "transparent",
                                    colour = "transparent"))
  
  # Save plot
  ggsave(paste0("figs/figs_SA/figS7_PercDiff_SQ_Holidays_",n_proj_type,".jpg"),
         width = 21, height = 15, dpi = 300)
}


# # Plot: facet -------------------------------------------------------------
# 
# # List of plots
# l_plots_figS7 <- vector(mode = "list", length = 4)
# i <- 1
# 
# for(n_outcome in v_outcome){
#   
#   l_plots_figS7[[i]] <- ggplot(data = subset(df_figS7_PercDiff, Outcome == n_outcome),
#                                aes(x = dates,
#                                    y = value,
#                                    colour = Outcome,
#                                    ymin = lb,
#                                    ymax = ub,
#                                    shape = Outcome)) +
#     geom_point(#position = position_dodge(width = 5),
#       size = 4) +
#     geom_errorbar(#position = position_dodge(width = 5),
#       size = 1.1) +
#     facet_wrap(~ proj_type_label, ncol = 2, scales = "fixed") +
#     # geom_ribbon(aes(x = dates,
#     #                 y = mean_percent,
#     #                 ymin = lb_percent,
#     #                 ymax = ub_percent,
#     #                 fill = Intervention),
#     #             alpha = 0.4,
#     #             colour = NA) +
#     scale_color_manual(values = rep(v_colors_colour, 2)) +
#     scale_fill_manual(values =  rep(v_colors_fill, 2))  +
#     scale_shape_manual(values = rep(v_colors_shape, 2)) +
#     guides(shape = guide_legend(title = ""),
#            color = guide_legend(title = "",
#                                 ncol = 1,
#                                 override.aes = list(#linetype = v_linetype,
#                                   color    = rep(v_colors_colour, 2),
#                                   # fill     = v_colors_fill,
#                                   shape    = rep(v_colors_shape, 2)))) +
#     geom_vline(xintercept = c(as.Date("2020-12-24"),
#                               as.Date("2021-01-06")),
#                show.legend = TRUE,
#                size = 0.5,
#                linetype = rep(c("dashed", "twodash"),2),
#                color = "#393F3E") +
#     # annotation_custom(start_holidays) +
#     # annotation_custom(end_holidays) +
#     scale_x_date(label_x_axis,
#                  date_labels = "%b/%d", #Change the format of the date in the x axis
#                  breaks = c(as.Date("2020-12-07"), as.Date("2020-12-14"),
#                             as.Date("2020-12-21"), as.Date("2020-12-28"),
#                             as.Date("2021-01-04"), as.Date("2021-01-11"),
#                             as.Date("2021-01-18"), as.Date("2021-01-25"),
#                             as.Date("2021-02-01"), as.Date("2021-02-08"),
#                             as.Date("2021-02-15"), as.Date("2021-02-22"),
#                             as.Date("2021-03-01"))
#                  # breaks = "week"
#     ) +
#     scale_y_continuous("",
#                        breaks = number_ticks(6),
#                        labels = scales::percent
#     ) + # trans = "log"
#     labs(subtitle = n_outcome
#     ) +
#     theme(text = element_text(size = text_size),
#           plot.caption = element_text(hjust = 0,
#                                       size = caption_size),
#           axis.text.x = element_text(angle = 90,
#                                      hjust = 0.5,
#                                      vjust = 0.5,
#                                      colour = "black"),
#           axis.text.y = element_text(colour = "black"),
#           panel.background = element_rect(fill = "white",
#                                           colour = "gray",
#                                           size = 0.15,
#                                           linetype = "solid"),
#           strip.background = element_rect(fill   = "transparent",
#                                           colour = "transparent"),
#           panel.border = element_rect(colour = "black",
#                                       fill = NA,
#                                       size = 0.7),
#           legend.position = "none", # c(0.44,0.5),
#           legend.text = element_text(size = 20),
#           legend.justification = "center",
#           legend.direction = "vertical",
#           legend.title = element_blank(),
#           legend.key = element_rect(fill   = "transparent",
#                                     colour = "transparent"))
#   
#   if(n_outcome == "Incident cases"){
#     l_plots_figS7[[i]] <- l_plots_figS7[[i]] +
#       labs(title = "Difference between higher social contacts and status quo")
#   }
#   
#   i <- i + 1
# }
# 
# # Combine graphs
# plot_S7 <- ggarrange(l_plots_figS7[[1]], 
#                      NULL,
#                      l_plots_figS7[[3]],
#                      NULL,
#                      l_plots_figS7[[2]], 
#                      NULL,
#                      l_plots_figS7[[4]],
#                      ncol = 1,
#                      heights = c(1, 0.05, 1, 0.05, 1, 0.05, 1),
#                      # labels = c("C","","D"), 
#                      hjust = 0,vjust = 1.5,
#                      font.label = list(size = 26, face = "bold"))
# 
# # Save plot
# ggsave(paste0("figs/figs_SA/figS7_PercDiff_SQ_Holidays_SA_SQ.jpg"),
#        plot = plot_S7, width = 17, height = 32, dpi = 300)
