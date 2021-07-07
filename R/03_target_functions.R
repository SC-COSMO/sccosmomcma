#' Generate calibration targets for Mexico
#'
#' \code{gen_targets} generates state-specific calibration targets for 
#' Mexico based on symptomatic cases. 
#'
#' @param v_states_calib States in Mexico used for calibration.
#' @param n_time_stamp   Date at which COVID19 series were released by Mexican SSA
#' @param n_date_ini     Initial date to generate targets
#' @param n_date_last    Last date to generate targets
#' @return 
#' A data.frame with state-specific calibration targets for Mexico.
#' @export
#'
gen_targets <-function(v_states_calib = "MCMA", 
                       n_time_stamp   = "2020-09-13",
                       n_date_ini     = NULL,
                       n_date_last    = "2020-08-31"){
  
  # Covid-19-mx-data by state and age_groups
  df_covid_ssa_state_age_groups <- fread(paste0("data-raw/covid_ssa_state_age_groups_",n_time_stamp,".csv"))

  # Covid-19-mx-data by state
  df_covid_mx <- fread(paste0("data-raw/covid_ssa_state_",n_time_stamp,".csv"))
  
  # Set ISO abbreviation by state
  v_acron_sts <- c("AGU","BCN","BCS","CAM","CHP","CHH","COA","COL","DUR","GUA",
                   "GRO","HID","JAL","CMX","MIC","MOR","NAY","NLE","OAX",
                   "PUE","QUE","ROO","SLP","SIN","SON","MEX","TAB","TAM","TLA",
                   "VER","YUC","ZAC","NTL","MCMA")
  
  df_covid_mx$state <- as.character(df_covid_mx$state)
  df_covid_mx$state[df_covid_mx$state == "Mexico"] <- "National"
  df_covid_mx$state <- as.factor(df_covid_mx$state)
  df_covid_mx$abbrev_state <- ""
  
  df_covid_ssa_state_age_groups$state <- as.character(df_covid_ssa_state_age_groups$state)
  df_covid_ssa_state_age_groups$state[df_covid_ssa_state_age_groups$state == "Mexico"] <- "National"
  df_covid_ssa_state_age_groups$state <- as.factor(df_covid_ssa_state_age_groups$state)
  df_covid_ssa_state_age_groups$abbrev_state <- ""
   
  v_states <- as.character(unique(df_covid_mx$state))
  names(v_states) <- v_acron_sts
  
  
  for(state_i in v_states){
    df_covid_mx$abbrev_state[df_covid_mx$state == state_i] <-  names(v_states)[which(v_states==state_i)]
    df_covid_ssa_state_age_groups$abbrev_state[df_covid_ssa_state_age_groups$state == state_i] <-  names(v_states)[which(v_states==state_i)]
  }
  
  ### Cumulative targets ###
  
  ### DIAGNOSED INFECTIONS
  df_targets_cases <- df_covid_mx %>%
    filter(state %in% v_states_calib,
           var_outcome == "Confirmed") %>% 
    ungroup() %>%
    group_by(country, state, county, population, var_outcome) %>% 
    # mutate(max_inc = max(incident_cases)) %>%
    # filter(Date <= Date[which(incident_cases == max_inc)]) %>% 
    rename(value = cum_cases) %>%
    mutate(Date = as.Date(date),
           Date0 = Date - Date[1],
           DateLast = n_date_last, 
           Target = "Cumulative confirmed infections",
           type = "Target",
           series = "Confirmed") %>%
    select(series, type, Target, population, Date, Date0, DateLast, value, 
           time_stamp, abbrev_state) %>%
    mutate(lb = epitools::pois.exact(x = value, pt = population)$lower*population,
           ub = epitools::pois.exact(x = value, pt = population)$upper*population,
           se = (ub-lb)/(1.96*2))
  
  ## Select targets from a starting and finishing date
  if (!is.null(n_date_ini) & !is.null(n_date_last)){
    df_targets_cases <- df_targets_cases %>% 
      filter(Date >= n_date_ini & Date <= n_date_last)
  }
  if(!is.null(n_date_last)){
    df_targets_cases <- df_targets_cases %>% 
      filter(Date <= n_date_last)
  }
  if(!is.null(n_date_ini)){
    df_targets_cases <- df_targets_cases %>% 
      filter(Date >= n_date_ini)
  }
  
  ### DIAGNOSED DEATHS
  df_targets_deaths <- df_covid_mx %>%
    filter(state %in% v_states_calib,
           var_outcome == "Deaths") %>% 
    ungroup() %>%
    group_by(country, state, county, population, var_outcome) %>% 
    # mutate(max_inc = max(incident_cases)) %>%
    # filter(Date <= Date[which(incident_cases == max_inc)]) %>% 
    rename(value = cum_cases) %>%
    mutate(Date = as.Date(date),
           Date0 = Date - Date[1],
           DateLast = n_date_last, 
           Target = "Cumulative COVID19 deaths infections",
           type = "Target",
           series = "Total DX COVID deaths") %>%
    select(series, type, Target, population, Date, Date0, DateLast, value, 
           time_stamp, abbrev_state) %>%
    mutate(lb = epitools::pois.exact(x = value, pt = population)$lower*population,
           ub = epitools::pois.exact(x = value, pt = population)$upper*population,
           se = (ub-lb)/(1.96*2))
  
  ## Select targets from a starting and finishing date
  if (!is.null(n_date_ini) & !is.null(n_date_last)){
    df_targets_deaths <- df_targets_deaths %>% 
      filter(Date >= n_date_ini & Date <= n_date_last)
  }
  if(!is.null(n_date_last)){
    df_targets_deaths <- df_targets_deaths %>% 
      filter(Date <= n_date_last)
  }
  if(!is.null(n_date_ini)){
    df_targets_deaths <- df_targets_deaths %>% 
      filter(Date >= n_date_ini)
  }
  
  #### Incident targets ####
  
  ### DIAGNOSED INFECTIONS
  df_targets_cases_inc <- df_covid_mx %>%
    filter(state %in% v_states_calib,
           var_outcome == "Confirmed") %>% 
    ungroup() %>%
    group_by(country, state, county, population, var_outcome) %>% 
    # mutate(max_inc = max(incident_cases)) %>%
    # filter(Date <= Date[which(incident_cases == max_inc)]) %>% 
    rename(value = new_cases) %>%
    mutate(Date = as.Date(date),
           Date0 = Date - Date[1],
           DateLast = n_date_last, 
           Target = "Incident confirmed infections",
           type = "Target",
           series = "Incident confirmed") %>%
    select(series, type, Target, population, Date, Date0, DateLast, value, 
           time_stamp, abbrev_state) %>%
    mutate(lb = epitools::pois.exact(x = value, pt = population)$lower*population,
           ub = epitools::pois.exact(x = value, pt = population)$upper*population,
           se = (ub-lb)/(1.96*2))
  
  ## Select targets from a starting and finishing date
  if (!is.null(n_date_ini) & !is.null(n_date_last)){
    df_targets_cases_inc <- df_targets_cases_inc %>% 
      filter(Date >= n_date_ini & Date <= n_date_last)
  }
  if(!is.null(n_date_last)){
    df_targets_cases_inc <- df_targets_cases_inc %>% 
      filter(Date <= n_date_last)
  }
  if(!is.null(n_date_ini)){
    df_targets_cases_inc <- df_targets_cases_inc %>% 
      filter(Date >= n_date_ini)
  }
  
  ### DIAGNOSED DEATHS
  df_targets_deaths_inc <- df_covid_mx %>%
    filter(state %in% v_states_calib,
           var_outcome == "Deaths") %>% 
    ungroup() %>%
    group_by(country, state, county, population, var_outcome) %>% 
    # mutate(max_inc = max(incident_cases)) %>%
    # filter(Date <= Date[which(incident_cases == max_inc)]) %>% 
    rename(value = new_cases) %>%
    mutate(Date = as.Date(date),
           Date0 = Date - Date[1],
           DateLast = n_date_last,
           Target = "Incident COVID19 deaths infections",
           type = "Target",
           series = "Total incident DX COVID deaths") %>%
    select(series, type, Target, population, Date, Date0, DateLast, value, 
           time_stamp, abbrev_state) %>%
    mutate(lb = epitools::pois.exact(x = value, pt = population)$lower*population,
           ub = epitools::pois.exact(x = value, pt = population)$upper*population,
           se = (ub-lb)/(1.96*2))
  
  ## Select targets from a starting and finishing date
  if (!is.null(n_date_ini) & !is.null(n_date_last)){
    df_targets_deaths_inc <- df_targets_deaths_inc %>% 
      filter(Date >= n_date_ini & Date <= n_date_last)
  }
  if(!is.null(n_date_last)){
    df_targets_deaths_inc <- df_targets_deaths_inc %>% 
      filter(Date <= n_date_last)
  }
  if(!is.null(n_date_ini)){
    df_targets_deaths_inc <- df_targets_deaths_inc %>% 
      filter(Date >= n_date_ini)
  }
  
  #### Covariance matrices ####
  ### Incident cases
  acf_cases_inc   <- acf(df_targets_cases_inc$value, lag.max = 1)
  v_acf_cases_inc <- acf_cases_inc$acf[2]^(0:(length(df_targets_cases_inc$value)-1))
  m_cov_cases_inc <- acor2cov(v_acor = v_acf_cases_inc, 
                              v_sd = df_targets_cases_inc$se)
  ### Incident deaths
  acf_deaths_inc   <- acf(df_targets_deaths_inc$value, lag.max = 1)
  v_acf_deaths_inc <- acf_deaths_inc$acf[2]^(0:(length(df_targets_deaths_inc$value)-1))
  m_cov_deaths_inc <- acor2cov(v_acor = v_acf_deaths_inc, 
                               v_sd = df_targets_deaths_inc$se)
  
  #### Cumulative targets by age groups ####
  ### Confirmed cases
  df_targets_cases_age <- df_covid_ssa_state_age_groups %>%
    filter(state %in% v_states_calib,
           var_outcome == "Confirmed") %>%
    ungroup() %>%
    group_by(country, state, county, population, var_outcome) %>%
    rename(value = cum_cases) %>%
    mutate(Date = as.Date(date),
           Date0 = Date - Date[1],
           DateLast = n_date_last,
           Target = "Cumulative confirmed infections by age group",
           type = "Target",
           series = "Confirmed") %>%
    select(series, type, Target, Date, Date0, DateLast, value,
           time_stamp, age_groups, abbrev_state) %>%
    spread(key = age_groups, value = value, fill = 0)

  ## Select targets from a starting and finishing date
  if (!is.null(n_date_ini) & !is.null(n_date_last)){
    df_targets_cases_age <- df_targets_cases_age %>%
      filter(Date >= n_date_ini & Date <= n_date_last)
  }
  if(!is.null(n_date_last)){
    df_targets_cases_age <- df_targets_cases_age %>%
      filter(Date <= n_date_last)
  }
  if(!is.null(n_date_ini)){
    df_targets_cases_age <- df_targets_cases_age %>%
      filter(Date >= n_date_ini)
  }

  ### Deaths
  df_targets_deaths_age <- df_covid_ssa_state_age_groups %>%
    filter(state %in% v_states_calib,
           var_outcome == "Deaths") %>%
    ungroup() %>%
    group_by(country, state, county, population, var_outcome) %>%
    rename(value = cum_cases) %>%
    mutate(Date = as.Date(date),
           Date0 = Date - Date[1],
           DateLast = n_date_last,
           Target = "Cumulative COVID19 deaths infections by age group",
           type = "Target",
           series = "Confirmed") %>%
    select(series, type, Target, Date, Date0, DateLast, value,
           time_stamp, age_groups, abbrev_state) %>%
    spread(key = age_groups, value = value, fill = 0)

  ## Select targets from a starting and finishing date
  if (!is.null(n_date_ini) & !is.null(n_date_last)){
    df_targets_deaths_age <- df_targets_deaths_age %>%
      filter(Date >= n_date_ini & Date <= n_date_last)
  }
  if(!is.null(n_date_last)){
    df_targets_deaths_age <- df_targets_deaths_age %>%
      filter(Date <= n_date_last)
  }
  if(!is.null(n_date_ini)){
    df_targets_deaths_age <- df_targets_deaths_age %>%
      filter(Date >= n_date_ini)
  }
  
  #### Indicent targets by age groups ####
  ### Confirmed cases
  df_targets_cases_inc_age <- df_covid_ssa_state_age_groups %>%
    filter(state %in% v_states_calib,
           var_outcome == "Confirmed") %>%
    ungroup() %>%
    group_by(country, state, county, population, var_outcome) %>%
    rename(value = new_cases) %>%
    mutate(Date = as.Date(date),
           Date0 = Date - Date[1],
           DateLast = n_date_last,
           Target = "Incident confirmed infections by age group",
           type = "Target",
           series = "Confirmed") %>%
    select(series, type, Target, Date, Date0, DateLast, value,
           time_stamp, age_groups, abbrev_state) %>%
    spread(key = age_groups, value = value, fill = 0)
  
  ## Select targets from a starting and finishing date
  if (!is.null(n_date_ini) & !is.null(n_date_last)){
    df_targets_cases_inc_age <- df_targets_cases_inc_age %>%
      filter(Date >= n_date_ini & Date <= n_date_last)
  }
  if(!is.null(n_date_last)){
    df_targets_cases_inc_age <- df_targets_cases_inc_age %>%
      filter(Date <= n_date_last)
  }
  if(!is.null(n_date_ini)){
    df_targets_cases_inc_age <- df_targets_cases_inc_age %>%
      filter(Date >= n_date_ini)
  }
  
  ### Deaths
  df_targets_deaths_inc_age <- df_covid_ssa_state_age_groups %>%
    filter(state %in% v_states_calib,
           var_outcome == "Deaths") %>%
    ungroup() %>%
    group_by(country, state, county, population, var_outcome) %>%
    rename(value = new_cases) %>%
    mutate(Date = as.Date(date),
           Date0 = Date - Date[1],
           DateLast = n_date_last,
           Target = "Incident COVID19 deaths infections by age group",
           type = "Target",
           series = "Confirmed") %>%
    select(series, type, Target, Date, Date0, DateLast, value,
           time_stamp, age_groups, abbrev_state) %>%
    spread(key = age_groups, value = value, fill = 0)
  
  ## Select targets from a starting and finishing date
  if (!is.null(n_date_ini) & !is.null(n_date_last)){
    df_targets_deaths_inc_age <- df_targets_deaths_inc_age %>%
      filter(Date >= n_date_ini & Date <= n_date_last)
  }
  if(!is.null(n_date_last)){
    df_targets_deaths_inc_age <- df_targets_deaths_inc_age %>%
      filter(Date <= n_date_last)
  }
  if(!is.null(n_date_ini)){
    df_targets_deaths_inc_age <- df_targets_deaths_inc_age %>%
      filter(Date >= n_date_ini)
  }
  
  #### Combine ALL targets ####
  df_all_targets <- bind_rows(df_targets_cases, 
                              df_targets_deaths,
                              df_targets_cases_inc,
                              df_targets_deaths_inc) %>%
    arrange(state, series, Date0)
  
  df_all_targets_age <- bind_rows(df_targets_cases_age,
                                  df_targets_deaths_age,
                                  df_targets_cases_inc_age,
                                  df_targets_deaths_inc_age)
  
  return(list(df_all_targets     = df_all_targets,
              cases              = df_targets_cases,
              deaths             = df_targets_deaths,
              cases_inc          = df_targets_cases_inc,
              deaths_inc         = df_targets_deaths_inc,
              cov_cases_inc      = m_cov_cases_inc,
              cov_deaths_inc     = m_cov_deaths_inc,
              df_all_targets_age = df_all_targets_age,
              cases_age          = df_targets_cases_age,
              deaths_age         = df_targets_deaths_age,
              cases_inc_age     = df_targets_cases_inc_age,
              deaths_inc_age     = df_targets_deaths_inc_age
              ))
              
}

#' Number of ticks for \code{ggplot2} plots
#'
#' Function for determining number of ticks on axis of \code{ggplot2} plots.
#' @param n integer giving the desired number of ticks on axis of
#' \code{ggplot2} plots. Non-integer values are rounded down.
#' @section Details:
#' Based on function \code{pretty}.
#' @export
number_ticks <- function(n) {
  function(limits) {
    pretty(limits, n + 1)
  }
}

#' Plot Targets
#'
#' \code{plot_targets} plots targets.
#'
#' @param l_targets List with calibration targets
#' @param print_plot Logical. Print plots
#' @param save_plot Logical. Save plots
#' @return
#' A ggplot2 object.
#' 
plot_targets <- function(l_targets, print_plot = TRUE, save_plot = TRUE){
  ### Obtain time stamp
  n_time_stamp <- max(l_targets$df_all_targets$DateLast, na.rm = TRUE)
  abbrev_state <- unique(l_targets$df_all_targets$abbrev_state)
  #abbrev_state <- unique(as.numeric(l_targets$df_all_targets$state))
  
  #### Plot aggregated targets ####
  gg_targets <- ggplot(l_targets$df_all_targets, aes(x = Date, y = value,
                                                     ymin = lb, ymax = ub,
                                                     color = type, shape = type)) +
    facet_wrap( ~ series , scales = "free_y") +
    geom_point(size = 3) + # for target: shape = 8
    geom_errorbar() +
    scale_shape_manual(values = c(1)) +
    scale_color_manual(values = c("black")) +
    scale_x_date("",
                 date_labels = "%m/%d/%y",
                 breaks = number_ticks(8)) +
    scale_y_continuous("",
                       breaks = number_ticks(6)) +
    # labs(title = “Cumulative total infectious cases of COVID-19 Mexico”,
    #      subtitle = “Days since first case in each state”,
    #      x = “Days since first case”, y = “Cumulative cases”) +
    theme_bw(base_size = 16) +
    theme(legend.position = "",
          axis.text.x = element_text(angle = 60, hjust = 1, size = 10))
  
  #### Plot age-specific targets ####
  # df_all_targets_age_long <- l_targets$df_all_targets_age %>% 
  #   gather(key = age_groups, value = cum_events, 
  #          -Target, -Date, -Date0, -DateLast, -time_stamp, -country, -state, -county, -population,
  #          -var_outcome, -series, -type, -abbrev_state)
  # df_all_targets_age_long$age_groups <- ordered(df_all_targets_age_long$age_groups, 
  #                                               unique(df_all_targets_age_long$age_groups))
  # 
  # gg_targets_ages_fill <- ggplot(df_all_targets_age_long, aes(x = Date, y = cum_events, 
  #                                                             fill = age_groups)) +
  #   geom_area(alpha = 0.7, position = "fill", stat = "identity") +
  #   facet_wrap(~ Target) + 
  #   scale_fill_ordinal("Age groups") + 
  #   xlab("Date") +
  #   ylab("Proportion") +
  #   theme_bw(base_size = 12)
  # 
  # gg_targets_ages <- ggplot(df_all_targets_age_long, aes(x = Date, y = cum_events, 
  #                                                        fill = age_groups)) +
  #   geom_area(alpha = 0.7) +
  #   facet_wrap(~ Target) + 
  #   scale_fill_ordinal("Age groups") + 
  #   xlab("Date") +
  #   ylab("Counts") +
  #   theme_bw(base_size = 12)
  
  if(print_plot){
    print(gg_targets)
    # print(gg_targets_ages_fill)
    # print(gg_targets_ages)
  }
  if(save_plot){
    ggsave(paste0("figs/03_targets_", abbrev_state, "_",n_time_stamp,".pdf"), 
           plot = gg_targets,
           width = 12, height = 8)
    # ggsave(paste0("figs/03_targets_ages_fill_", abbrev_state,"_",n_time_stamp,".pdf"), 
    #        plot = gg_targets_ages_fill,
    #        width = 10, height = 6)
    # ggsave(paste0("figs/03_targets_ages_", abbrev_state,"_", n_time_stamp,".pdf"), 
    #        plot = gg_targets_ages,
    #        width = 10, height = 6)
  }
}

#' Covariance matrix from autocorrelation coefficients and standard
#' deviations (or standard errors)
#'
#' \code{acor2cov} computes covariance matrix from autocorrelation coefficients 
#' and standard deviations (or standard errors) of means.
#'
#' @param v_acor Vector of autocorrelation coefficients indexed by time
#' @param v_sd Vector of standard deviations (or standard errors)
#' @return
#' A covariance matrix.
#'
acor2cov <- function(v_acor, v_sd){
  n_y <- length(v_acor)
  m_cor <- matrix(0, n_y, n_y)
  for (i in 0:(n_y-1)){
    m_cor[i+1, (i+1):n_y] <- v_acor[1:(n_y-i)]
  }
  m_cor <- m_cor + t(m_cor)
  diag(m_cor) <- 1
  m_cov <- MBESS::cor2cov(m_cor, v_sd)#cor2cov(covar, v_sd^2)
  return(m_cov)
}
