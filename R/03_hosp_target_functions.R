#' Generate hospitalization targets.
#'
#' \code{gen_hosp_targets} generates targets for hospitalization parameters. 
#'
#' @param n_time_stamp   Date at which hosp series were released
#' @param n_date_ini     Initial date to generate targets
#' @param n_date_last    Last date to generate targets
#' @return 
#' A data.frame with state-specific calibration targets for Mexico.
#' @export
#'
gen_hosp_targets <-function(n_time_stamp   = "2020-12-21",
                            n_date_ini     = NULL,
                            n_date_last    = "2020-12-07"){

  # Mexico City's hospitalization data. Source: ADIP
  load(paste0("C:/Data/vgracia/GitHub/Mexico-MCMA/data/df_hosp_MCMA_",n_time_stamp,".RData"))
  df_hosp_MCMA_raw <- df_hosp_MCMA

  df_hosp_MCMA <- df_hosp_MCMA[order(df_hosp_MCMA$dates),]
  df_hosp_MCMA$population <- unique(df_pop_state_cty_age$tot_pop)
  

# Total hospitalizations --------------------------------------------------

  df_targets_hosp <- df_hosp_MCMA %>%
    mutate(state = "MCMA",
           country = "Mexico",
           county = "MCMA",
           abbrev_state = "MCMA",
           var_outcome = "Hospitalized",
           Date = as.Date(dates),
           Date0 = Date - Date[1],
           DateLast = n_date_last, 
           Target = "Total hospitalizations",
           type = "Target",
           series = "Total hospitalizations",
           value = hosp_tot,
           time_stamp = n_time_stamp) %>%
    select(country, state, county, var_outcome, series, type, Target, population, 
           Date, Date0, DateLast, value, time_stamp, abbrev_state) %>%
    mutate(lb = epitools::pois.exact(x = value, pt = population)$lower*population,
           ub = epitools::pois.exact(x = value, pt = population)$upper*population,
           se = (ub-lb)/(1.96*2))
  
  # Select targets from a starting and finishing date
  if (!is.null(n_date_ini) & !is.null(n_date_last)){
    df_targets_hosp <- df_targets_hosp %>% 
      filter(Date >= n_date_ini & Date <= n_date_last)
  }
  if(!is.null(n_date_last)){
    df_targets_hosp <- df_targets_hosp %>% 
      filter(Date <= n_date_last)
  }
  if(!is.null(n_date_ini)){
    df_targets_hosp <- df_targets_hosp %>% 
      filter(Date >= n_date_ini)
  }
  

# Hosp with ventilator ----------------------------------------------------

  df_targets_vent <- df_hosp_MCMA %>%
    mutate(state = "MCMA",
           country = "Mexico",
           county = "MCMA",
           abbrev_state = "MCMA",
           var_outcome = "Hospitalized with ventilator",
           Date = as.Date(dates),
           Date0 = Date - Date[1],
           DateLast = n_date_last, 
           Target = "Hospitalizations with ventilator",
           type = "Target",
           series = "Hospitalizations with ventilator",
           value = bed_icu_tot,
           time_stamp = n_time_stamp) %>%
    select(country, state, county, var_outcome, series, type, Target, population, 
           Date, Date0, DateLast, value, time_stamp, abbrev_state) %>%
    mutate(lb = epitools::pois.exact(x = value, pt = population)$lower*population,
           ub = epitools::pois.exact(x = value, pt = population)$upper*population,
           se = (ub-lb)/(1.96*2))
  
  # Select targets from a starting and finishing date
  if (!is.null(n_date_ini) & !is.null(n_date_last)){
    df_targets_vent <- df_targets_vent %>% 
      filter(Date >= n_date_ini & Date <= n_date_last)
  }
  if(!is.null(n_date_last)){
    df_targets_vent <- df_targets_vent %>% 
      filter(Date <= n_date_last)
  }
  if(!is.null(n_date_ini)){
    df_targets_vent <- df_targets_vent %>% 
      filter(Date >= n_date_ini)
  }
  

# Hosp without ventilator -------------------------------------------------

  df_targets_novent <- df_hosp_MCMA %>%
    mutate(state = "MCMA",
           country = "Mexico",
           county = "MCMA",
           abbrev_state = "MCMA",
           var_outcome = "Hospitalized without ventilator",
           Date = as.Date(dates),
           Date0 = Date - Date[1],
           DateLast = n_date_last, 
           Target = "Hospitalizations without ventilator",
           type = "Target",
           series = "Hospitalizations without ventilator",
           value = bed_noicu_tot,
           time_stamp = n_time_stamp) %>%
    select(country, state, county, var_outcome, series, type, Target, population, 
           Date, Date0, DateLast, value, time_stamp, abbrev_state) %>%
    mutate(lb = epitools::pois.exact(x = value, pt = population)$lower*population,
           ub = epitools::pois.exact(x = value, pt = population)$upper*population,
           se = (ub-lb)/(1.96*2))
  
  # Select targets from a starting and finishing date
  if (!is.null(n_date_ini) & !is.null(n_date_last)){
    df_targets_novent <- df_targets_novent %>% 
      filter(Date >= n_date_ini & Date <= n_date_last)
  }
  if(!is.null(n_date_last)){
    df_targets_novent <- df_targets_novent %>% 
      filter(Date <= n_date_last)
  }
  if(!is.null(n_date_ini)){
    df_targets_novent <- df_targets_novent %>% 
      filter(Date >= n_date_ini)
  }
  

# Covariance matrices -----------------------------------------------------

  # Total hospitalizations
  acf_hosp   <- acf(df_targets_hosp$value, lag.max = 1)
  v_acf_hosp <- acf_hosp$acf[2]^(0:(length(df_targets_hosp$value)-1))
  m_cov_hosp <- hosp_acor2cov(v_acor = v_acf_hosp, 
                              v_sd = df_targets_hosp$se)
  
  # Ventilator
  acf_vent   <- acf(df_targets_vent$value, lag.max = 1)
  v_acf_vent <- acf_vent$acf[2]^(0:(length(df_targets_vent$value)-1))
  m_cov_vent <- hosp_acor2cov(v_acor = v_acf_vent, 
                               v_sd = df_targets_vent$se)
  
  # No ventilator
  acf_novent   <- acf(df_targets_novent$value, lag.max = 1)
  v_acf_novent <- acf_novent$acf[2]^(0:(length(df_targets_novent$value)-1))
  m_cov_novent <- hosp_acor2cov(v_acor = v_acf_novent, 
                         v_sd = df_targets_novent$se)


# Combine ALL targets -----------------------------------------------------

  df_all_targets <- bind_rows(df_targets_hosp, 
                              df_targets_vent,
                              df_targets_novent) %>%
    arrange(state, series, Date0)
  
  return(list(df_all_targets = df_all_targets,
              hosp           = df_targets_hosp,
              vent           = df_targets_vent,
              novent         = df_targets_novent,
              cov_hosp       = m_cov_hosp,
              cov_vent       = m_cov_vent,
              cov_novent     = m_cov_novent))
  
}

#' Plot hospitalization targets
#'
#' \code{plot_hosp_targets} plots hospitalization targets.
#'
#' @param l_hosp_targets List. Calibration targets.
#' @param print_plot Logical. Print plot if TRUE.
#' @param save_plot Logical. Save plot if TRUE.
#' @return
#' A ggplot2 object.
#' 
plot_hosp_targets <- function(l_hosp_targets, 
                              print_plot = TRUE, 
                              save_plot = TRUE){
  
  # Obtain time stamp and state abbreviation
  n_time_stamp <- max(l_hosp_targets$df_all_targets$DateLast, na.rm = TRUE)
  abbrev_state <- unique(l_hosp_targets$df_all_targets$abbrev_state)
  
  # Plot aggregated targets
  gg_targets <- ggplot(l_hosp_targets$df_all_targets, aes(x = Date, y = value,
                                                     ymin = lb, ymax = ub,
                                                     color = type, shape = type)) +
    facet_wrap( ~ series , scales = "free_y") +
    geom_point(size = 3) + # for target: shape = 8
    geom_errorbar() +
    scale_shape_manual(values = c(1)) +
    scale_color_manual(values = c("black")) +
    scale_x_date(breaks = number_ticks(8)) +
    scale_y_continuous(breaks = number_ticks(6)) +
    labs(title = "Hospitalizations of Covid-19 Mexico",
         # subtitle = “”,
         x = "",
         y = "") +
    theme_bw(base_size = 16) +
    theme(legend.position = "",
          axis.text.x = element_text(angle = 60, hjust = 1, size = 10))
  
  if(print_plot){
    print(gg_targets)
  }
  if(save_plot){
    ggsave(paste0("figs/05_hosp_targets_", abbrev_state, "_",n_time_stamp,".pdf"), 
           plot = gg_targets,
           width = 12, height = 8)
  }
}

#' Covariance matrix from autocorrelation coefficients and standard
#' deviations (or standard errors)
#'
#' \code{hosp_acor2cov} computes covariance matrix from autocorrelation coefficients 
#' and standard deviations (or standard errors) of means.
#'
#' @param v_acor Vector of autocorrelation coefficients indexed by time
#' @param v_sd Vector of standard deviations (or standard errors)
#' @return
#' A covariance matrix.
#'
hosp_acor2cov <- function(v_acor, v_sd){
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
