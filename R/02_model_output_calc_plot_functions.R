############### GENERIC/HELPER FUNCTIONS  #####################

#' Plot COSMO data
#' 
#' \code{plot_cosmo_data} a generic function called by various wrapper 
#' functions to plot a data frame returned by one of our \code{calc} methods by 
#' age group and/or in total over time.
#' 
#' @param df_plotting Data.frame of data to be plotted.
#' @param proportion Flag (default is FALSE) of whether to divide the outcome
#' of interest by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' the outcome of interest by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") data
#' should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @param label_y y-axis label.
#' @return 
#' A ggplot object.
#' @export
plot_cosmo_data <- function(df_plotting, 
                            proportion = FALSE, 
                            age_proportion = FALSE,
                            only_all = TRUE, 
                            print_plot = TRUE,
                            label_y)
  
{
  if(only_all & age_proportion){
    stop("You cannot compute age-specific proportions if wanting results for all the population")
  }
  
  if(proportion){
    if(age_proportion){
      df_popize_plot <- calc_popsize_totals(l_out_cosmo)
      df_plotting <- data.frame(time = df_plotting[,1], df_plotting[, -1]/df_popize_plot[, -1], 
                                   check.names = FALSE)
      label_y <- paste("Within-group", label_y, "(proportions)")
    } else{
      df_popize_plot <- calc_popsize_totals(l_out_cosmo)
      df_plotting <- data.frame(time = df_plotting[, 1], df_plotting[, -1]/df_popize_plot[, ncol(df_popize_plot)], 
                                   check.names = FALSE)
      label_y <- paste(label_y, "(proportions)")
    }
  }
  df_plotting_lng <- reshape2::melt(df_plotting, id.vars = "time", 
                                       variable.name = "Age group",
                                       value.name = "value")
  if(only_all){
    df_plotting_lng <- subset(df_plotting_lng, `Age group` == "All")
  } else{
    df_plotting_lng <- subset(df_plotting_lng, `Age group` != "All")
  }
  
  gg_cosmo_plot <- ggplot(data = df_plotting_lng, 
                          aes(x = time, y = value, 
                          color = `Age group`)) +
    geom_line() +
    ylab(label_y) +
    theme_bw(base_size = 16)
  if(print_plot) {print(gg_cosmo_plot)}
  return(gg_cosmo_plot)
  
}

############################### POPULATION SIZE #####################


#' Total population 
#' 
#' \code{calc_popsize_totals} calculate total population by age group and 
#' overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total population by age group and overall 
#' as columns over time.
#' @export
calc_popsize_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  i_popsize <- 1:(l_params_all$n_states_alive * l_params_all$n_ages) ### ADD variable of number of susceptibles and recovered
  
  v_names_seir <- l_params_all$v_names_states_alive# c("S", l_params_all$v_names_exp_states, l_params_all$v_names_inf_states, "R")
  
  m_out_cosmo  <- as.matrix(df_out_cosmo[, (1 + i_popsize)])
  
  a_pop_dist <- array(t(m_out_cosmo), # c(df_out_cosmo[, (1 + i_popsize)], recursive = TRUE)
                      dim = list(length(v_names_ages), length(v_names_seir), length(df_out_cosmo$time)), 
                      dimnames = list(v_names_ages,  v_names_seir, df_out_cosmo$time))
  
  m_popsize <- t(apply(a_pop_dist, 3, rowSums))
  
  df_popsize <- data.frame(time = df_out_cosmo$time, 
                           m_popsize,
                           N_ALL = rowSums(m_popsize), 
                           check.names = FALSE)
  colnames(df_popsize)[-1] <- c(levels(v_names_ages), "All")
  return(df_popsize)
}

#' Plot total population
#' 
#' \code{plot_popsize_totals} plot total population by age group and/or 
#' in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") population
#' should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_popsize_totals <- function(l_out_cosmo,
                                only_all = TRUE, 
                                print_plot = TRUE){
  
  df_popize_plot <- calc_popsize_totals(l_out_cosmo)
  df_popsize_plot_lng <- reshape2::melt(df_popize_plot, id.vars = "time", 
                                        variable.name = "Age group",
                                        value.name = "Population")
  if(only_all){
    df_popsize_plot_lng <- subset(df_popsize_plot_lng, `Age group` == "All")
  } else{
    df_popsize_plot_lng <- subset(df_popsize_plot_lng, `Age group` != "All")
  }
  
  gg_popsize <- ggplot(data = df_popsize_plot_lng, 
                       aes(x = time, y = Population, 
                           color = `Age group`)) +
    geom_line() +
    theme_bw(base_size = 16)
  if(print_plot) {print(gg_popsize)}
  return(gg_popsize)
}

############################### INFECTIONS #####################

#' Total cumulative infections
#' 
#' \code{calc_infcum_totals} calculate total number of cumulative infections by age group and 
#' overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters
#' @return 
#' A data.frame with the total number of cumulative infections by age group and overall 
#' as columns over time.
#' @export
calc_infcum_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  v_names_inf_tot <- paste("Itot", v_names_ages, sep = "_")
  df_Inftot <- data.frame(time = df_out_cosmo$time, 
                          df_out_cosmo[, v_names_inf_tot],
                          ItotALL = rowSums(df_out_cosmo[, v_names_inf_tot]), 
                          check.names = FALSE)
  colnames(df_Inftot)[-1] <- c(levels(v_names_ages), "All")
  return(df_Inftot)
}

#' Plot total cumulative infections
#' 
#' \code{plot_infcum_totals} plot total number of cumulative infections by 
#' age group and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param proportion Flag (default is FALSE) of whether to divide total number
#' of cumulative infections by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total number of cumulative infections by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") total 
#' number of cumulative infections should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_infcum_totals <- function(l_out_cosmo, 
                               proportion = FALSE, 
                               age_proportion = FALSE,
                               only_all = TRUE, 
                               print_plot = TRUE){
  
  df_Inftot_plot <- calc_infcum_totals(l_out_cosmo)
  label_y <- "Cumulative Infections"
  gg_Inftot <- plot_cosmo_data(df_plotting = df_Inftot_plot, 
                               proportion = proportion,
                               age_proportion = age_proportion,
                               only_all = only_all,
                               print_plot = print_plot,
                               label_y = label_y)
  return(gg_Inftot)
}

#' Total incident infections
#' 
#' \code{calc_infinc_totals} calculate total number of incident infections by 
#' age group and overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total number of incident infections by age group and overall 
#' as columns over time.
#' @export
calc_infinc_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  v_names_inf_tot <- paste("Itot", v_names_ages, sep = "_")
  m_ItotInc <- rbind(0,
                     matrixStats::colDiffs(as.matrix(df_out_cosmo[, v_names_inf_tot])))
  
  df_InfInctot <- data.frame(time = df_out_cosmo$time, 
                             m_ItotInc,
                             ItotIncALL = rowSums(m_ItotInc), 
                             check.names = FALSE)
  colnames(df_InfInctot)[-1] <- c(levels(v_names_ages), "All")
  return(df_InfInctot)
}

#' Plot total incident infections
#' 
#' \code{plot_infinc_totals} plot total number of incident infections by age group 
#' and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param proportion Flag (default is FALSE) of whether to divide total number of
#' incident infections by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total number of incident infections by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") total 
#' number of incident infections should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_infinc_totals <- function(l_out_cosmo, 
                               proportion = FALSE, 
                               age_proportion = FALSE,
                               only_all = TRUE, 
                               print_plot = TRUE){
  
  df_InfInctot_plot <- calc_infinc_totals(l_out_cosmo)
  label_y <- "Incident Infections"
  gg_InfInctot <- plot_cosmo_data(df_plotting = df_InfInctot_plot, 
                               proportion = proportion,
                               age_proportion = age_proportion,
                               only_all = only_all,
                               print_plot = print_plot,
                               label_y = label_y)
  return(gg_InfInctot)
}

#' Total cumulative diagnosed infectious infections
#' 
#' \code{calc_idxcum_totals} calculate total number of cumulative diagnosed 
#' infectious infections by age group and overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total number of cumulative diagnosed infectious infections by age group and overall 
#' as columns over time.
#' @export
calc_idxcum_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  v_names_idx_tot <- paste("IDXtot", v_names_ages, sep = "_")
  df_IDXtot <- data.frame(time = df_out_cosmo$time, 
                          df_out_cosmo[, v_names_idx_tot],
                          IDXtotALL = rowSums(df_out_cosmo[, v_names_idx_tot]), 
                          check.names = FALSE)
  colnames(df_IDXtot)[-1] <- c(levels(v_names_ages), "All")
  return(df_IDXtot)
}

#' Plot total cumulative diagnosed infectious infections
#' 
#' \code{plot_idxcum_totals} plot total number of cumulative diagnosed infectious 
#' infections by age group and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param proportion Flag (default is FALSE) of whether to divide total number of 
#' cumulative diagnosed infectious by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total number of cumulative diagnosed infectious by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") total 
#' total number of cumulative diagnosed infectious should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_idxcum_totals <- function(l_out_cosmo, 
                               proportion = FALSE, 
                               age_proportion = FALSE,
                               only_all = TRUE, 
                               print_plot = TRUE){

  df_IDXtot_plot <- calc_idxcum_totals(l_out_cosmo)
  label_y <- "Incident Diagnosed Infectious Infections"
  gg_IDXtot <- plot_cosmo_data(df_plotting = df_IDXtot_plot, 
                                  proportion = proportion,
                                  age_proportion = age_proportion,
                                  only_all = only_all,
                                  print_plot = print_plot,
                                  label_y = label_y)
  return(gg_IDXtot)
}

#' Total incident diagnosed infectious infections
#' 
#' \code{calc_idxinc_totals} calculate total number of incident diagnosed infectious 
#' infections by age group and overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total number of incident diagnosed infectious infections 
#' by age group and overall as columns over time.
#' @export
calc_idxinc_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  v_names_idx_tot <- paste("IDXtot", v_names_ages, sep = "_")
  m_IndIDXtotInc <- rbind(0,
                          matrixStats::colDiffs(as.matrix(df_out_cosmo[, v_names_idx_tot])))
  df_InfIDXInctot <- data.frame(time = df_out_cosmo$time, 
                                m_IndIDXtotInc,
                                IDXInctotALL = rowSums(m_IndIDXtotInc), 
                                check.names = FALSE)
  colnames(df_InfIDXInctot)[-1] <- c(levels(v_names_ages), "All")
  return(df_InfIDXInctot)
}

#' Plot total incident diagnosed infectious infections
#' 
#' \code{plot_idxinc_totals} plot total number of incident diagnosed 
#' infectious infections by age group and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param proportion Flag (default is FALSE) of whether to divide total number of 
#' incident diagnosed infectious by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total number of incident diagnosed infectious by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") total 
#' total number of incident diagnosed infectious should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_idxinc_totals <- function(l_out_cosmo, 
                               proportion = FALSE, age_proportion = FALSE,
                               only_all = TRUE, print_plot = TRUE){
  df_InfIDXInctot_plot <- calc_idxinc_totals(l_out_cosmo)
  label_y <- "Incident Diagnosed Infectious Infections"
  gg_InfIDXInctot <- plot_cosmo_data(df_plotting = df_InfIDXInctot_plot, 
                               proportion = proportion,
                               age_proportion = age_proportion,
                               only_all = only_all,
                               print_plot = print_plot,
                               label_y = label_y)
  return(gg_InfIDXInctot)
}

#' Total cumulative diagnosed infections (E and I)
#' 
#' \code{calc_dxcum_totals} calculate total number of cumulative diagnosed 
#' infections (E and I) by age group and overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total number of cumulative diagnosed infections 
#' (E and I) by age group and overall as columns over time.
#' @export
calc_dxcum_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  v_names_dx_tot <- paste("DXtot", v_names_ages, sep = "_")
  df_DXtot <- data.frame(time = df_out_cosmo$time, 
                         df_out_cosmo[, v_names_dx_tot],
                         DXtotALL = rowSums(df_out_cosmo[, v_names_dx_tot]), 
                         check.names = FALSE)
  colnames(df_DXtot)[-1] <- c(levels(v_names_ages), "All")
  return(df_DXtot)
}

#' Plot total cumulative diagnosed infections (E and I)
#' 
#' \code{plot_dxcum_totals} plot total number of cumulative diagnosed 
#' infections (E and I) by age group and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param proportion Flag (default is FALSE) of whether to divide total 
#' cumulative diagnosed infections (E and I) by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total cumulative diagnosed infections (E and I) by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") total 
#' total cumulative diagnosed infections (E and I) should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_dxcum_totals <- function(l_out_cosmo, 
                              proportion = FALSE, 
                              age_proportion = FALSE,
                              only_all = TRUE, 
                              print_plot = TRUE){
  
  df_DXtot_plot <- calc_dxcum_totals(l_out_cosmo)
  label_y <- "Cumulative diagnosed cases"
  gg_DXtot <- plot_cosmo_data(df_plotting = df_DXtot_plot, 
                                     proportion = proportion,
                                     age_proportion = age_proportion,
                                     only_all = only_all,
                                     print_plot = print_plot,
                                     label_y = label_y)
  return(gg_DXtot)
}

#' Total incident diagnosed infections (E and I)
#' 
#' \code{calc_dxinc_totals} calculate total number of incident diagnosed 
#' infections (E and I) by age group and overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total number of incident diagnosed infections (E and I) 
#' by age group and overall as columns over time.
#' @export
calc_dxinc_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  v_names_dx_tot <- paste("DXtot", v_names_ages, sep = "_")
  m_DXtotInc <- rbind(0,
                      matrixStats::colDiffs(as.matrix(df_out_cosmo[, v_names_dx_tot])))
  df_DXInctot <- data.frame(time = df_out_cosmo$time, 
                            m_DXtotInc,
                            DXInctotALL = rowSums(m_DXtotInc), 
                            check.names = FALSE)
  colnames(df_DXInctot)[-1] <- c(levels(v_names_ages), "All")
  return(df_DXInctot)
}

#' Plot total incident diagnosed cases (E and I)
#' 
#' \code{plot_dxinc_totals} plot total number of incident diagnosed 
#' cases (E and I) by age group and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param proportion Flag (default is FALSE) of whether to divide total incident
#'  diagnosed cases (E and I) by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total incident diagnosed cases (E and I) by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") total 
#' total incident diagnosed cases (E and I) should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_dxinc_totals <- function(l_out_cosmo, 
                              proportion = FALSE, age_proportion = FALSE,
                              only_all = TRUE, print_plot = TRUE){
  df_DXInctot_plot <- calc_dxinc_totals(l_out_cosmo)
  label_y <- "Incident diagnosed cases"
  gg_DXInctot <- plot_cosmo_data(df_plotting = df_DXInctot_plot, 
                              proportion = proportion,
                              age_proportion = age_proportion,
                              only_all = only_all,
                              print_plot = print_plot,
                              label_y = label_y)
  return(gg_DXInctot)
}

#' Total diagnosed infections (E and I)
#' 
#' \code{calc_dx_totals} calculate total number of diagnosed 
#' infections (E and I)by age group and overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total number of diagnosed infections (E and I) 
#' by age group and overall as columns over time.
#' @export
calc_dx_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  m_names_dx <- matrix(NA, 
                       nrow = length(v_names_ages), 
                       ncol = length(l_params_all$v_names_dx_states))
  df_DXtot_ages <- as.data.frame(matrix(NA, 
                                        ncol = length(v_names_ages),
                                        nrow = (l_params_all$n_t+1)))
  colnames(df_DXtot_ages) <- v_names_ages
  for(i in 1:length(v_names_ages)){
    m_names_dx[i, ] <- paste(l_params_all$v_names_dx_states, v_names_ages[i], sep = "_")
    df_DXtot_ages[, i] <- rowSums(df_out_cosmo[, m_names_dx[i, ]])
  }
  
  df_DXtot <- data.frame(time = df_out_cosmo$time, 
                         df_DXtot_ages,
                         DXtotALL = rowSums(df_DXtot_ages), 
                         check.names = FALSE)
  colnames(df_DXtot)[-1] <- c(levels(v_names_ages), "All")
  return(df_DXtot)
}

#' Plot total diagnosed infections (E and I)
#' 
#' \code{plot_dx_totals} plot total number of diagnosed 
#' infections (E and I) by age group and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param proportion Flag (default is FALSE) of whether to divide total 
#' diagnosed infections (E and I) by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total diagnosed infections (E and I) by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") total 
#' total diagnosed infections (E and I) should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_dx_totals <- function(l_out_cosmo, 
                           proportion = FALSE, 
                           age_proportion = FALSE,
                           only_all = TRUE,
                           print_plot = TRUE){
  
  df_DXtot_plot <- calc_dx_totals(l_out_cosmo)
  label_y <- "Diagnosed Infections"
  gg_DXtot <- plot_cosmo_data(df_plotting = df_DXtot_plot, 
                                 proportion = proportion,
                                 age_proportion = age_proportion,
                                 only_all = only_all,
                                 print_plot = print_plot,
                                 label_y = label_y)
  return(gg_DXtot)
}


#' Total infections
#' 
#' \code{calc_inf_totals} calculate total number of infections 
#' by age group and overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total number of infections by age group and overall 
#' as columns over time.
#' @export
calc_inf_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  m_names_inf <- matrix(NA, 
                        nrow = length(v_names_ages), 
                        ncol = length(l_params_all$v_names_inftot_states))
  df_Inftot_ages <- as.data.frame(matrix(NA, 
                                         ncol = length(v_names_ages),
                                         nrow = (l_params_all$n_t+1)))
  colnames(df_Inftot_ages) <- v_names_ages
  for(i in 1:length(v_names_ages)){
    m_names_inf[i, ] <- paste(l_params_all$v_names_inftot_states,
                              v_names_ages[i], sep = "_") 
    df_Inftot_ages[, i] <- rowSums(df_out_cosmo[, m_names_inf[i, ]])
  }
  
  df_Inftot <- data.frame(time = df_out_cosmo$time, 
                          df_Inftot_ages,
                          ItotALL = rowSums(df_Inftot_ages), 
                          check.names = FALSE)
  colnames(df_Inftot)[-1] <- c(levels(v_names_ages), "All")
  return(df_Inftot)
}

#' Plot total infections
#' 
#' \code{plot_inf_totals} plot total number of infections by 
#' age group and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param proportion Flag (default is FALSE) of whether to divide total
#' number of infections by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total number of infections by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") total
#' number of infections should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_inf_totals <- function(l_out_cosmo, 
                            proportion = FALSE, 
                            age_proportion = FALSE,
                            only_all = TRUE, 
                            print_plot = TRUE){
  
  df_Inftot_plot <- calc_inf_totals(l_out_cosmo)
  label_y <- "Infections"
  gg_Inftot <- plot_cosmo_data(df_plotting = df_Inftot_plot, 
                              proportion = proportion,
                              age_proportion = age_proportion,
                              only_all = only_all,
                              print_plot = print_plot,
                              label_y = label_y)
  return(gg_Inftot)
}

#' Total diagnosed infections
#' 
#' \code{calc_idx_totals} calculate total number of diagnosed 
#' infections by age group and overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total number of diagnosed infections by age group and overall 
#' as columns over time.
#' @export
calc_idx_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  m_names_inf <- matrix(NA, 
                        nrow = length(v_names_ages), 
                        ncol = length(l_params_all$v_names_infidx_states))
  df_IDXtot_ages <- as.data.frame(matrix(NA, 
                                         ncol = length(v_names_ages),
                                         nrow = (l_params_all$n_t+1)))
  colnames(df_IDXtot_ages) <- v_names_ages
  for(i in 1:length(v_names_ages)){
    m_names_inf[i, ] <- paste(l_params_all$v_names_infidx_states, 
                              v_names_ages[i], sep = "_")  
    df_IDXtot_ages[, i] <- rowSums(df_out_cosmo[, m_names_inf[i, ]])
  }
  
  df_IDXtot <- data.frame(time = df_out_cosmo$time, 
                          df_IDXtot_ages,
                          IDXtotALL = rowSums(df_IDXtot_ages), 
                          check.names = FALSE)
  colnames(df_IDXtot)[-1] <- c(levels(v_names_ages), "All")
  return(df_IDXtot)
}

#' Plot total diagnosed infections
#' 
#' \code{plot_idx_totals} plot total number of diagnosed 
#' infections by age group and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param proportion Flag (default is FALSE) of whether to divide total number of 
#' diagnosed infections (E and I) by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total number of diagnosed infections (E and I) by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") 
#' total number of diagnosed infections (E and I) should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_idx_totals <- function(l_out_cosmo, 
                            proportion = FALSE,
                            age_proportion = FALSE,
                            only_all = TRUE, 
                            print_plot = TRUE){
  
  df_IDXtot_plot <- calc_idx_totals(l_out_cosmo)
  label_y <- "Diagnosed Infections"
  gg_IDXtot <- plot_cosmo_data(df_plotting = df_IDXtot_plot, 
                               proportion = proportion,
                               age_proportion = age_proportion,
                               only_all = only_all,
                               print_plot = print_plot,
                               label_y = label_y)
  return(gg_IDXtot)
}

#' Total infections (E and I)
#' 
#' \code{calc_expinf_totals} calculate total number of  
#' infections (E and I) by age group and overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total number of infections (E and I) by age group and overall 
#' as columns over time.
#' @export
calc_expinf_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  v_names_expinf_states <- c(l_params_all$v_names_exp_states,
                             l_params_all$v_names_inf_states,
                             l_params_all$v_names_expidx_states,
                             l_params_all$v_names_infidx_states)
  
  m_names_dx <- matrix(NA, 
                       nrow = length(v_names_ages), 
                       ncol = length(v_names_expinf_states))
  df_ExpInftot_ages <- as.data.frame(matrix(NA, 
                                        ncol = length(v_names_ages),
                                        nrow = (l_params_all$n_t+1)))
  colnames(df_ExpInftot_ages) <- v_names_ages
  for(i in 1:length(v_names_ages)){
    m_names_dx[i, ] <- paste(v_names_expinf_states, v_names_ages[i], sep = "_")
    df_ExpInftot_ages[, i] <- rowSums(df_out_cosmo[, m_names_dx[i, ]])
  }
  
  df_ExpInftot <- data.frame(time = df_out_cosmo$time, 
                         df_ExpInftot_ages,
                         ExpInftotALL = rowSums(df_ExpInftot_ages), 
                         check.names = FALSE)
  colnames(df_ExpInftot)[-1] <- c(levels(v_names_ages), "All")
  return(df_ExpInftot)
}

#' Plot total infections (E and I)
#' 
#' \code{plot_expinf_totals} plot total number of  
#' infections (E and I) by age group and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters
#' @param proportion Flag (default is FALSE) of whether to divide total number of
#' infections (E and I) by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total number of infections (E and I) by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") total
#' number of infections (E and I) should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_expinf_totals <- function(l_out_cosmo, 
                           proportion = FALSE, age_proportion = FALSE,
                           only_all = TRUE, print_plot = TRUE){
  df_ExpInftot_plot <- calc_expinf_totals(l_out_cosmo)
  label_y <- "Infections (E and I)"
  gg_ExpInftot <- plot_cosmo_data(df_plotting = df_ExpInftot_plot, 
                              proportion = proportion,
                              age_proportion = age_proportion,
                              only_all = only_all,
                              print_plot = print_plot,
                              label_y = label_y)
  return(gg_ExpInftot)
}


############################### SYMPTOMATICS #####################

#' Total symptomatic infections (E and I)
#' 
#' \code{calc_sx_totals} calculate total number of symptomatic 
#' infections (E and I) by age group and overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total number of symptomatic infections (E and I) 
#' by age group and overall as columns over time.
#' @export
calc_sx_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  v_state_names <- c(l_params_all$v_names_exp_l2_states, 
                     l_params_all$v_names_expidx_l2_states,
                     l_params_all$v_names_inf_l2_states,
                     l_params_all$v_names_infidx_l2_states)
  
  m_names_sx <- matrix(NA, 
                       nrow = length(v_names_ages), 
                       ncol = length(v_state_names))
  df_SXtot_ages <- as.data.frame(matrix(NA, 
                                        ncol = length(v_names_ages),
                                        nrow = (l_params_all$n_t+1)))
  colnames(df_SXtot_ages) <- v_names_ages
  for(i in 1:length(v_names_ages)){
    m_names_sx[i, ] <- paste(v_state_names, v_names_ages[i], sep = "_")
    df_SXtot_ages[, i] <- rowSums(df_out_cosmo[, m_names_sx[i, ]])
  }
  
  df_SXtot <- data.frame(time = df_out_cosmo$time, 
                         df_SXtot_ages,
                         SXtotALL = rowSums(df_SXtot_ages), 
                         check.names = FALSE)
  colnames(df_SXtot)[-1] <- c(levels(v_names_ages), "All")
  return(df_SXtot)
}

#' Plot total symptomatic infections
#' 
#' \code{plot_sx_totals} plot total number of symptomatic (E and I) 
#' by age group and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param proportion Flag (default is FALSE) of whether to divide total
#'  number of symptomatic (E and I) by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total number of symptomatic (E and I)  by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") total 
#' number of symptomatic (E and I) should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_sx_totals <- function(l_out_cosmo, 
                           proportion = FALSE, 
                           age_proportion = FALSE,
                           only_all = TRUE, 
                           print_plot = TRUE){
  
  df_Sxtot_plot <- calc_sx_totals(l_out_cosmo)
  label_y <- "Symptomatic Infections"
  gg_Sxtot <- plot_cosmo_data(df_plotting = df_Sxtot_plot, 
                                 proportion = proportion,
                                 age_proportion = age_proportion,
                                 only_all = only_all,
                                 print_plot = print_plot,
                                 label_y = label_y)
  return(gg_Sxtot)
}


########################### DEATHS ############################


#' Total Covid-19 deaths
#' 
#' \code{calc_deaths_totals} calculate total number of Covid-19 
#' deaths by age group and overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total number of Covid-19 deaths by age group and overall 
#' as columns over time.
#' @export
calc_deaths_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  v_names_deaths_tot <- paste("Dcov", v_names_ages, sep = "_")
  df_Dcov <- data.frame(time = df_out_cosmo$time, 
                        df_out_cosmo[, v_names_deaths_tot],
                        DcovALL = rowSums(df_out_cosmo[, v_names_deaths_tot]), 
                        check.names = FALSE)
  colnames(df_Dcov)[-1] <- c(levels(v_names_ages), "All")
  return(df_Dcov)
}

#' Plot total Covid-19 deaths
#' 
#' \code{plot_deaths_totals} plot total number of deaths from Covid-19 
#' by age group and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param proportion Flag (default is FALSE) of whether to divide total deaths 
#' from Covid-19 by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total deaths from Covid-19 by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") 
#' total deaths from Covid-19 should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_deaths_totals <- function(l_out_cosmo, 
                               proportion = FALSE, 
                               age_proportion = FALSE,
                               only_all = TRUE, 
                               print_plot = TRUE){

  df_Dcov_plot <- calc_deaths_totals(l_out_cosmo)
  label_y <- "Covid-19 deaths"
  gg_Dcov <- plot_cosmo_data(df_plotting = df_Dcov_plot, 
                               proportion = proportion,
                               age_proportion = age_proportion,
                               only_all = only_all,
                               print_plot = print_plot,
                               label_y = label_y)
  return(gg_Dcov)
}

#' Total Covid-19 deaths from diagnosed infections
#' 
#' \code{calc_deathsdx_totals} calculate total number of Covid-19 deaths 
#' from diagnosed infections by age group and overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total number of Covid-19 deaths from diagnosed 
#' infections by age group and overall as columns over time.
#' @export
calc_deathsdx_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  v_names_deaths_tot <- paste("DcovDX", v_names_ages, sep = "_")
  df_DcovDX <- data.frame(time = df_out_cosmo$time, 
                          df_out_cosmo[, v_names_deaths_tot],
                          DcovDXALL = rowSums(df_out_cosmo[, v_names_deaths_tot]), 
                          check.names = FALSE)
  colnames(df_DcovDX)[-1] <- c(levels(v_names_ages), "All")
  return(df_DcovDX)
}

#' Plot total Covid-19 deaths
#' 
#' \code{plot_deathsdx_totals} plot total number of Covid-19 deaths 
#' from diagnosed infections by age group and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param proportion Flag (default is FALSE) of whether to divide total
#' number of Covid-19 deaths by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total number of Covid-19 deaths by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") 
#' total number of Covid-19 deaths should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_deathsdx_totals <- function(l_out_cosmo, 
                                 proportion = FALSE, 
                                 age_proportion = FALSE,
                                 only_all = TRUE, 
                                 print_plot = TRUE){
  df_DcovDX_plot <- calc_deathsdx_totals(l_out_cosmo)
  label_y <- "Diagnosed Covid-19 deaths"
  gg_DcovDX <- plot_cosmo_data(df_plotting = df_DcovDX_plot, 
                             proportion = proportion,
                             age_proportion = age_proportion,
                             only_all = only_all,
                             print_plot = print_plot,
                             label_y = label_y)
  return(gg_DcovDX)
}

#' Total Incident Covid-19 deaths from diagnosed infections
#' 
#' \code{calc_incdeathsdx_totals} calculate total incident number of 
#' Covid-19 deaths from diagnosed infections by age group and overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total incident number of Covid-19 deaths from diagnosed 
#' infections by age group and overall as columns over time.
#' @export
calc_incdeathsdx_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  v_names_deaths_tot <- paste("DcovDX", v_names_ages, sep = "_")
  m_DcovDXInc <- rbind(0,
                       matrixStats::colDiffs(as.matrix(df_out_cosmo[, v_names_deaths_tot])))
  
  df_DcovDXInc <- data.frame(time = df_out_cosmo$time,
                             m_DcovDXInc,
                             DcovDXIncALL = rowSums(m_DcovDXInc), 
                             check.names = FALSE)
  colnames(df_DcovDXInc)[-1] <- c(levels(v_names_ages), "All")
  return(df_DcovDXInc)
}

#' Plot total incident Covid-19 deaths from diagnosed infections
#' 
#' \code{plot_incdeathsdx_totals} plot total incident number of Covid-19 
#' deaths from diagnosed infections by age group and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters
#' @param proportion Flag (default is FALSE) of whether to divide total
#' incident number of Covid-19 deaths by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total incident number of Covid-19 deaths by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") total incident
#'  number of Covid-19 deaths should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @return 
#' A ggplot object.
#' @export
plot_incdeathsdx_totals <- function(l_out_cosmo, 
                                    proportion = FALSE, 
                                    age_proportion = FALSE,
                                    only_all = TRUE,
                                    print_plot = TRUE){
  
  df_IncDcovDX_plot <- calc_incdeathsdx_totals(l_out_cosmo)
  label_y <- "Diagnosed Covid-19 deaths"
  gg_IncDcovDX <- plot_cosmo_data(df_plotting = df_IncDcovDX_plot, 
                                  proportion = proportion,
                                  age_proportion = age_proportion,
                                  only_all = only_all,
                                  print_plot = print_plot,
                                  label_y = label_y)
  return(gg_IncDcovDX)
}

############################### RECOVERED/IMMUNE #####################

#' Total Recovered
#' 
#' \code{calc_rec_totals} calculate total number of recovered 
#' by age group and overall over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @return 
#' A data.frame with the total number of recovered by age group and overall 
#' as columns over time.
#' @export
calc_rec_totals <- function(l_out_cosmo){
  df_out_cosmo <- l_out_cosmo$df_out_cosmo
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  v_names_recovered_tot <- paste("R", v_names_ages, sep = "_")
  df_Rectot <- data.frame(time = df_out_cosmo$time, 
                          df_out_cosmo[, v_names_recovered_tot],
                          RecALL = rowSums(df_out_cosmo[, v_names_recovered_tot]), 
                          check.names = FALSE)
  colnames(df_Rectot)[-1] <- c(levels(v_names_ages), "All")
  return(df_Rectot)
}

#' Plot total recovered
#' 
#' \code{plot_rec_totals} plot total number of recovered 
#' by age group and/or in total over time.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param proportion Flag (default is FALSE) of whether to divide total 
#' number of recovered by total population.
#' @param age_proportion Flag (default is FALSE) of whether to divide 
#' total number of recovered by age-specific total population.
#' @param only_all Flag (default is TRUE) of whether only overall ("All") total 
#' number of recovered should be plotted.
#' @param print_plot Flag (default is TRUE) of whether the plot should be printed
#' in screen.
#' @param only_all logical. Should only overall total infections be plotted
#' @return 
#' A ggplot object.
#' @export
plot_rec_totals <- function(l_out_cosmo, 
                            proportion = FALSE, 
                            age_proportion = FALSE,
                            only_all = TRUE, 
                            print_plot = TRUE){
  df_Rectot_plot <- calc_rec_totals(l_out_cosmo)
  label_y <- "Recovered"
  gg_Rectot <- plot_cosmo_data(df_plotting = df_Rectot_plot, 
                                  proportion = proportion,
                                  age_proportion = age_proportion,
                                  only_all = only_all,
                                  print_plot = print_plot,
                                  label_y = label_y)
  return(gg_Rectot)
}



############################### REPRODUCTION NUMBERS #####################

#' Calculate basic reproduction number according to Keeling & Rohani
#' 
#' \code{calc_basic_reproduction_number_kr} calculates the basic reproduction 
#' number according to Keeling & Rohani.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param v_time Vector with times at which basic reproduction number is 
#' desired to be calculated.
#' @references 
#' \enumerate{
#' \item Emilia Vynnycky and Richard G White. An introduction to infectious 
#' disease modelling. Oxford, New York: Oxford University Press, 2010. 
#' \item Keeling, M. J., & Rohani, P. Modeling Infectious Diseases in 
#' Humans and Animals. Princeton University Press, 2008.
#' }
#' @return 
#' A scalar with the basic reproduction number.
#' @export
calc_basic_reproduction_number_kr <- function(l_out_cosmo, 
                                              v_time = 2:10){
  ## Obtain cumulative new infections
  df_infcum_totals <- calc_infcum_totals(l_out_cosmo = l_out_cosmo)
  ## Select days on which Rt is desired
  df_infcum_totals <- df_infcum_totals[v_time, ]
  ## Run log-linear regression to estimate exponential growth rate
  lm_Rt <- lm(log(All) ~ time, data = df_infcum_totals)
  ## Estimate of exponential growth rate
  lambda <- lm_Rt$coefficients[2]
  ## Obtain epidemic model structure parameters
  n_exp_states <- l_out_cosmo$l_params_all$n_exp_states
  n_inf_states <- l_out_cosmo$l_params_all$n_inf_states
  r_sigma <- l_out_cosmo$l_params_all$v_sigma[1]
  r_gamma <- l_out_cosmo$l_params_all$v_gamma[1]
  # avg_exp_dur <- n_exp_states/(r_sigma) # /parms$timestep
  # avg_inf_dur <- n_inf_states/(r_gamma) # /parms$timestep
  ## Estimate reproductive number using eq. 3.14 of Keeling and Rohani (2008)
  R0_est <- (lambda*((lambda/(r_sigma*n_exp_states)) + 1)^n_exp_states)/
    (r_gamma*(1 - ((lambda/(r_gamma*n_inf_states)) + 1)^(-n_inf_states)))
  names(R0_est) <- "R0_est"
  return(R0_est)
}

#' Estimate the time dependent reproduction number according to Wallinga & 
#' Teunis
#' 
#' \code{calc_reproduction_number_wt} calculates the time-dependent reproduction 
#' number according to Wallinga & Teunis.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param v_time Vector with times at which the time-dependent reproduction 
#' number is desired to be calculated.
#' @param nsim_chosen Number of simulations used by the WT algorithm to 
#' compute confidence intervals which are often not used. Hence the number
#' is small. The default for WT in the package is 10,000 which uses a very
#' large amount of memory and is not appropriate for parallelization because
#' we will run out of RAM. The default is set to a much smaller number.
#' @references 
#' \enumerate{
#' \item Wallinga, J., & Teunis, P. (2004). Different epidemic curves for severe 
#' acute respiratory syndrome reveal similar impacts of control measures. 
#' American Journal of Epidemiology, 160(6), 509–516. 
#' https://doi.org/10.1093/aje/kwh255
#' }
#' @return 
#' A data.frame with with the time-dependent reproduction number.
#' @export
calc_reproduction_number_wt <- function(l_out_cosmo, v_time = 2:20, nsim_chosen = 200){
  ## Obtain incident new infections
  df_infinc_totals <- calc_infinc_totals(l_out_cosmo = l_out_cosmo)
  ## Select days on which Rt is desired
  df_infinc_totals <- df_infinc_totals %>%
                            filter(time %in% v_time)
  ## Obtain epidemic model structure parameters
  # n_exp_states <- l_out_cosmo$l_params_all$n_exp_states
  # n_inf_states <- l_out_cosmo$l_params_all$n_inf_states
  r_sigma <- l_out_cosmo$l_params_all$v_sigma[1]
  r_gamma <- l_out_cosmo$l_params_all$v_gamma[1]
  ## Derive mean and SD parameters for gamma distribution 
  mu_gen_time <- 1/r_sigma + 0.5*(1/r_gamma) # 5.8
  sd_gen_time <- sqrt(1^2 + ((1/2)*2.1)^2) # 3
  
  gen_time_dist <- R0::generation.time("gamma", c(mu_gen_time, sd_gen_time)) # From Jason Adrews: this is the generation interval, defined as a gamma distribution with a mean of 5.8 days
  Rt_est_wt     <- R0::est.R0.TD(epid = df_infinc_totals$All, 
                                 t = df_infinc_totals$time,
                                 begin = 1,
                                 end = last(df_infinc_totals$time)-first(df_infinc_totals$time)+1,                                  
                                 GT = gen_time_dist,
                                 nsim = nsim_chosen)
  df_Rt_est_wt <- data.frame(time = df_infinc_totals$time[-length(Rt_est_wt$R)],
                             Rt = Rt_est_wt$R[-length(Rt_est_wt$R)],
                             lb = Rt_est_wt$conf.int[-nrow(Rt_est_wt$conf.int), 1],
                             ub = Rt_est_wt$conf.int[-nrow(Rt_est_wt$conf.int), 2]
  )
  return(df_Rt_est_wt)
}
