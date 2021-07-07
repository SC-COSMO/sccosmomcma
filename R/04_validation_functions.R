#' Generate outputs for validation
#'
#' \code{validate_out} produces output to validate the model.
#'
#' @param m_calib_post Matrix with calibrated parameters from posterior 
#' distribution
#' @param l_params_all List with all parameters of the decision model.
#' @param n_date_ini Initial date of calibration.
#' @param n_t_calib Number of calibration days.
#' @param n_lag_inf Lag in time series of infectious individuals.
#' @return
#' A list with outputs for validation.
#' 
validate_out <- function(m_calib_post, 
                         l_params_all,
                         n_date_ini, 
                         n_t_calib, 
                         n_lag_inf = 14,
                         l_dates_targets){
  
  ### Number of days in calibration
  v_targets_names <- names(l_dates_targets)
  l_v_dates  <- l_v_dates0 <- l_n_dates <- vector(mode = "list", length = length(l_dates_targets))
  names(l_v_dates)  <- names(l_v_dates0) <- names(l_n_dates) <- v_targets_names

  for(i in 1:length(l_dates_targets)){ # i = 1
    l_v_dates[[i]]  <- seq.Date(from = l_dates_targets[[i]][1],
                             to   = l_dates_targets[[i]][2],
                             by   ="day")
    
    l_v_dates0[[i]] <- 0:(length(l_v_dates[[i]])-1)
    l_n_dates[[i]]  <- length(l_v_dates0[[i]])
    
  }
  
  if(is.null(dim(m_calib_post))) { # If vector, change to matrix
    m_calib_post <- t(m_calib_post) 
  }
  
  ### Number of posterior samples
  n_samp <- nrow(m_calib_post)
  
  ### Print message
  print(paste0("Generating validation outcomes on ", 
               n_samp, 
               " parameter sets"))
  
  #### Initialize matrices to store model outputs ####
  ### Detected cases
  ## Incident
  m_out_infdx_inc <- matrix(NA, nrow = n_samp, ncol = l_n_dates[["cases_inc"]])
  colnames(m_out_infdx_inc) <- as.character(l_v_dates[["cases_inc"]])
  ## Cumulative
  m_out_infdx_cum <- matrix(NA, nrow = n_samp, ncol = l_n_dates[["cases"]])
  colnames(m_out_infdx_cum) <- as.character(l_v_dates[["cases"]])
  ### Detected COVID deaths
  ## Incident
  m_out_DXDcov_inc <- matrix(NA, nrow = n_samp, ncol = l_n_dates[["deaths_inc"]])
  colnames(m_out_DXDcov_inc) <- as.character(l_v_dates[["deaths_inc"]])
  ## Cumulative
  m_out_DXDcov_cum <- matrix(NA, nrow = n_samp, ncol = l_n_dates[["deaths"]])
  colnames(m_out_DXDcov_cum) <- as.character(l_v_dates[["deaths"]])
  
  #### Compute model-predicted outputs for each sample of posterior distribution ####
  ### Evaluate model at each posterior sample and store results
  # for(i in 1:3){ # i = 1
  #   l_out_post <- calibration_out(v_params_calib = m_calib_post[i, ], 
  #                                 l_params_all = l_params_all, 
  #                                 # v_states_calib = v_states_calib[i], 
  #                                 n_lag_inf = n_lag_inf, 
  #                                 n_lag_conf = 0)
  #   m_out_infdx_inc[i, ]  <- l_out_post$DXIncTot$value
  #   m_out_infdx_cum[i, ]  <- l_out_post$DXCumTot$value
  #   m_out_DXDcov_inc[i, ] <- l_out_post$DcovIncTot$value
  #   m_out_DXDcov_cum[i, ] <- l_out_post$DcovTot$value
  #   cat('\r', paste(round(i/n_samp * 100), "% done", sep = " ")) # display progress
  # }
  ### Evaluate model at each posterior sample and store results
  no_cores <- detectCores() - 1     # detect number of cores
  cl <- makeCluster(no_cores)       # initialize cluster object
  registerDoParallel(cl)
  opts <- list(attachExportEnv = TRUE)
  system.time(
    m_out_post_all <- foreach(i = 1:n_samp, .combine = rbind, .export = ls(globalenv()), # i = 1
                              .packages=c("sccosmomcma",
                                          "tidyverse",
                                          "dplyr",
                                          "lubridate",
                                          "epitools"),
                              .options.snow = opts) %dopar% { 
      # ## For progress bar.
      # if(!exists("pb")) pb <- tkProgressBar(title = "Parallel task for internal validation",
      #                                       min=1, max=n.coverage)
      # info <- sprintf("%s%% done", round(i/n_samp*100))
      # setTkProgressBar(pb, i, label = sprintf("Progress of simulations (%s)", info))
      ## Run model
      l_out_post <- calibration_out(v_params_calib = m_calib_post[i, ], 
                                    l_params_all = l_params_all, 
                                    #v_states_calib = v_states_calib[i], 
                                    n_lag_inf = n_lag_inf, 
                                    n_lag_conf = 0,
                                    l_dates_targets)
      v_out_post <- c(i, 
                      l_out_post$DXIncTot$value, 
                      l_out_post$DXCumTot$value, 
                      l_out_post$DcovIncTot$value, 
                      l_out_post$DcovTot$value)
      v_out_post
    }
  )
  stopCluster(cl)
  
  if(is.null(dim(m_out_post_all))) { # If vector, change to matrix
    m_out_post_all <- t(m_out_post_all) 
  }
  
  a_out_post_all      <- m_out_post_all[ , -1]
  if(is.null(dim(a_out_post_all))) { # If vector, change to matrix
    a_out_post_all <- t(a_out_post_all) 
  }
  
  v_cuts <- c(l_n_dates[["cases_inc"]],
              l_n_dates[["cases_inc"]]+l_n_dates[["cases"]],
              l_n_dates[["cases_inc"]]+l_n_dates[["cases"]] + l_n_dates[["deaths_inc"]])
  
  m_out_infdx_inc[, ]  <- a_out_post_all[, 1:v_cuts[1]]
  m_out_infdx_cum[, ]  <- a_out_post_all[, (v_cuts[1] + 1):v_cuts[2]]
  m_out_DXDcov_inc[, ] <- a_out_post_all[, (v_cuts[2] + 1):v_cuts[3]]
  m_out_DXDcov_cum[, ] <- a_out_post_all[, (v_cuts[3] + 1):dim(a_out_post_all)[2]]
  
  ### Create data frames with model predicted outputs
  df_out_infdx_inc <- data.frame(type = "Model", 
                                 Target = "Incident confirmed infections",
                                 m_out_infdx_inc, 
                                 check.names = FALSE)
  
  df_out_infdx_cum <- data.frame(type = "Model", 
                                 Target = "Cumulative confirmed infections",
                                 m_out_infdx_cum, 
                                 check.names = FALSE)
  
  df_out_DXDcov_inc <- data.frame(type = "Model", 
                                  Target = "Incident COVID19 deaths infections",
                                  m_out_DXDcov_inc, 
                                  check.names = FALSE)
  
  df_out_DXDcov_cum <- data.frame(type = "Model", 
                                  Target = "Cumulative COVID19 deaths infections",
                                  m_out_DXDcov_cum, 
                                  check.names = FALSE)
  
  ### Transform data frames to long format
  df_out_infdx_inc_lng <- reshape2::melt(df_out_infdx_inc, 
                                         id.vars = c("type", "Target"), 
                                         variable.name = "Date", 
                                         value.name = "Value")
  
  df_out_infdx_cum_lng <- reshape2::melt(df_out_infdx_cum, 
                                         id.vars = c("type", "Target"), 
                                         variable.name = "Date", 
                                         value.name = "Value")
  
  df_out_DXDcov_inc_lng <- reshape2::melt(df_out_DXDcov_inc, 
                                          id.vars = c("type", "Target"), 
                                          variable.name = "Date", 
                                          value.name = "Value")
  
  df_out_DXDcov_cum_lng <- reshape2::melt(df_out_DXDcov_cum, 
                                          id.vars = c("type", "Target"), 
                                          variable.name = "Date", 
                                          value.name = "Value")
  
  ### Compute posterior model-predicted 95% CI
  df_out_infdx_inc_summ <- summarise_data(df_out_infdx_inc_lng, 
                                          varname = "value",
                                          groupnames = c("type", "Target", "Date"))
  
  df_out_infdx_cum_summ <- summarise_data(df_out_infdx_cum_lng, 
                                          varname = "value",
                                          groupnames = c("type", "Target", "Date"))
  
  df_out_DXDcov_inc_summ <- summarise_data(df_out_DXDcov_inc_lng, 
                                           varname = "value",
                                           groupnames = c("type", "Target", "Date"))
  
  df_out_DXDcov_cum_summ <- summarise_data(df_out_DXDcov_cum_lng, 
                                           varname = "value",
                                           groupnames = c("type", "Target", "Date"))
  
  ### Generate data frame with time-series data
  df_out_all_summ <- bind_rows(df_out_infdx_inc_summ, 
                               df_out_infdx_cum_summ, 
                               df_out_DXDcov_inc_summ,
                               df_out_DXDcov_cum_summ)
  
  df_out_all_summ$Date <- as.Date(df_out_all_summ$Date)
  
  ### Generate list to return all data frames
  l_out_post_summ <- list(df_out_infdx_inc_summ  = df_out_infdx_inc_summ,
                          df_out_infdx_cum_summ  = df_out_infdx_cum_summ,
                          df_out_DXDcov_inc_summ = df_out_DXDcov_inc_summ,
                          df_out_DXDcov_cum_summ = df_out_DXDcov_cum_summ,
                          df_out_all_summ        = df_out_all_summ)
  
  return(l_out_post_summ)
}      

#' Summarize posterior output
#'
#' \code{summarise_data} is used to to calculate the mean, standard deviation and 
#' 95% credible interval.
#' @param data Data frame.
#' @param varname Name of a column containing the variable.
#' @param groupnames Vector of column names to be used as grouping variables.
#' @return 
#' A data frame containing the posterior output.
#' @export
summarise_data <- function(data, varname, groupnames){
  # summary_func <- function(x, col){
  #   c(mean = mean(x[[col]], na.rm = TRUE),
  #     median = quantile(x[[col]], probs = 0.5, na.rm = TRUE, names = FALSE),
  #     sd = sd(x[[col]], na.rm=TRUE),
  #     lb = quantile(x[[col]], probs = 0.025, na.rm = TRUE, names = FALSE),
  #     ub = quantile(x[[col]], probs = 0.975, na.rm = TRUE, names = FALSE))
  # }
  # data_sum <- plyr::ddply(data, groupnames, .fun = summary_func, 
  #                         varname)
  # data_sum <- plyr::rename(data_sum, c("mean" = varname))
  df_summ <- data %>% 
    group_by_at(vars(groupnames)) %>%
    summarise(mean = mean(Value, na.rm = TRUE),
              median = quantile(Value, probs = 0.5, na.rm = TRUE, names = FALSE),
              sd = sd(Value, na.rm = TRUE),
              lb = quantile(Value, probs = 0.025, na.rm = TRUE, names = FALSE),
              ub = quantile(Value, probs = 0.975, na.rm = TRUE, names = FALSE)) 
  colnames(df_summ)[colnames(df_summ)=="mean"] <- varname
  return(df_summ)
}

#' Plot Targets
#'
#' \code{plot_targets} plots targets.
#'
#' @param l_model_out List with decision model outputs
#' @param l_targets List calibration targets
#' @param print_plot Logical. Prints plot if TRUE
#' @param return_plot Logical. Returns plot if TRUE
#' @return
#' A ggplot2 object.
#' 
plot_model_out_vs_targets <- function(df_all_targets, 
                                      df_model_out,
                                      print_plot = TRUE, 
                                      return_plot = TRUE){
  gg_model_out_vs_targets <- ggplot(df_all_targets, 
                                    aes(x = Date, y = value,
                                        ymin = lb, ymax = ub,
                                        color = type, shape = type)) +
    geom_point(size = 3) + # for target: shape = 8
    geom_errorbar() +
    geom_line(data = df_model_out, 
              aes(x = Date, y = value), col = "red") +
    geom_ribbon(data = df_model_out,
                aes(ymin = lb, ymax = ub), alpha = 0.4, fill = "red") +
    facet_wrap(~ Target, scales = "free_y") +
    scale_shape_manual(values = c(1, 1)) +
    scale_color_manual(values = c("red", "black")) +
    scale_x_date(breaks = number_ticks(20)) +
    scale_y_continuous(breaks = number_ticks(6)) +
    # labs(title = “Cumulative total infectious cases of COVID-19 Mexico”,
    #      subtitle = “Days since first case in each state”,
    #      x = “Days since first case”, y = “Cumulative cases”) +
    ylab("Counts") + 
    theme_bw(base_size = 16) +
    theme(legend.position = "",
          axis.text.x = element_text(angle = 60, hjust = 1, size = 10))
  if(print_plot == T){
    print(gg_model_out_vs_targets)
  }
  if(return_plot){
    return(list(gg_model_out_vs_targets = gg_model_out_vs_targets))  
  }
}
