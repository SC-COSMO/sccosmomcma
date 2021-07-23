#------------------------------------------------------------------------------#
# This script plots the calibrated posterior distributions.                    #
#                                                                              #
# Author:                                                                      #
#    Fernando Alarid-Escudero, PhD, <fernando.alarid@cide.edu>                 #
#                                                                              #
#------------------------------------------------------------------------------#
rm(list = ls())

# Load libraries and functions ---------------------------------------------

library(ggplot2)
library(ggcorrplot)
library(GGally)      # To use ggpairs
library(dampack)
library(data.table)
library(dplyr)

# Functions
source("R/03_calibration_functions.R")

# Number of ticks for plots: dampack package
number_ticks <- function(n) {
  function(limits) {
    pretty(limits, n + 1)
  }
}

for(n_proj_type in c("SQ", "SA")){
  
  cat("Projection scenario:",n_proj_type,"\n")
  
  # Projection type directory
  if(n_proj_type == "SA"){
    n_dir_proj <- "figs_SA/"
  }else if(n_proj_type == "SQ"){
    n_dir_proj <- "figs_paper/"
  }
  
  # Calibrated parameters
  load(paste0("output/03_map_output_IMIS_MCMA_",n_proj_type,"_2020-12-13.RData"))
  
  
  # NPI effectiveness under end-of-year holiday period ----------------------
  
  m_calib_post_eff <- m_calib_post
  m_calib_post_eff[, 3:7] <- 1-m_calib_post_eff[, 3:7]
  
  ## Multiplier
  # Right after NPI implementation in March 23
  mean(m_calib_post[, "r_soc_dist_factor"])
  quantile(m_calib_post[, "r_soc_dist_factor"], probs = c(0.025, 0.975))
  
  # At the end of calibration period
  mean(m_calib_post[, "r_soc_dist_factor_5"])
  quantile(m_calib_post[, "r_soc_dist_factor_5"], probs = c(0.025, 0.975))
  
  # Stricter NPI
  stricter_npi <- matrixStats::rowMins(m_calib_post[, 3:7])
  mean(stricter_npi)
  quantile(stricter_npi, probs = c(0.025, 0.975))
  
  stricter_npi_sq <- (stricter_npi - m_calib_post[, "r_soc_dist_factor_5"])/m_calib_post[, "r_soc_dist_factor_5"]
  
  cat("Stricter NPI: ",
      round(mean(stricter_npi_sq)*100),"%", 
      " [",round(quantile(stricter_npi_sq, probs = 0.025)*100,2),", ",
      round(quantile(stricter_npi_sq, probs = 0.975)*100,2),"]", " physical distancing compared to status quo\n",
      sep = "")
  
  # Holiday bump
  holiday_bump <- m_calib_post[, "r_soc_dist_factor_5"] + 0.3
  mean(holiday_bump)
  quantile(holiday_bump, probs = c(0.025, 0.975))
  
  holiday_bump_sq <- (holiday_bump - m_calib_post[, "r_soc_dist_factor_5"])/m_calib_post[, "r_soc_dist_factor_5"]
  
  cat("Holiday bump: ", round(mean(holiday_bump_sq)*100,2), "%", 
      " [",round(quantile(holiday_bump_sq, probs = 0.025)*100,2),", ",
      round(quantile(holiday_bump_sq, probs = 0.975)*100,2),"]", " contacts compared to status quo\n",
      sep = "")
  
  ## Effectiveness
  # Right after NPI implementation in March 23
  mean(m_calib_post_eff[, "r_soc_dist_factor"])
  quantile(m_calib_post_eff[, "r_soc_dist_factor"], probs = c(0.025, 0.975))
  
  # At the end of calibration period
  mean(m_calib_post_eff[, "r_soc_dist_factor_5"])
  quantile(m_calib_post_eff[, "r_soc_dist_factor_5"], probs = c(0.025, 0.975))
  
  holiday_bump_eff <- m_calib_post_eff[, "r_soc_dist_factor_5"] - 0.3
  mean(holiday_bump_eff)
  quantile(holiday_bump_eff, probs = c(0.025, 0.975))

  cat("Holiday bump: ", round(mean(holiday_bump_eff)*100,2), "%", 
      " [",round(quantile(holiday_bump_eff, probs = 0.025)*100,2),", ",
      round(quantile(holiday_bump_eff, probs = 0.975)*100,2),"]", " physical distancing - pre-pandemic levels\n",
      sep = "")
  
  stricter_npi_eff <- matrixStats::rowMaxs(m_calib_post_eff[, 3:7])
  mean(stricter_npi_eff)
  quantile(stricter_npi_eff, probs = c(0.025, 0.975))
  
  
  # Posterior visualization -------------------------------------------------
  
  ## Posterior pairwise correlation -----------------------------------------
  
  gg_post_pairs_corr <- ggpairs(data.frame(m_calib_post_eff),
                                upper = list(continuous = wrap("cor",
                                                               color = "black",
                                                               size = 5)),
                                diag = list(continuous = wrap("barDiag",
                                                              alpha = 0.8)),
                                lower = list(continuous = wrap("points", 
                                                               alpha = 0.3,
                                                               size = 0.5)),
                                columnLabels = c("beta",
                                                 "tau",
                                                 "eta[1]",
                                                 "eta[2]",
                                                 "eta[3]",
                                                 "eta[4]",
                                                 "eta[5]",
                                                 "nu[lb]",
                                                 "nu[ub]",
                                                 "nu[rate]",
                                                 "nu[mid]"),
                                labeller = "label_parsed") +
    theme_bw(base_size = 18) +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_text(size=6),
          axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank())

  # Save plot
  ggsave(gg_post_pairs_corr,
         filename = paste0("figs/",n_dir_proj,"figS2_posterior-IMIS-correlation-1k_",n_proj_type,".jpg"),
         width = 12, height = 8, dpi = 300)
  
  
  ## Posterior vs prior comparison ------------------------------------------
  
  # Prior vs posterior distribution
  m_calib_post_eff <- m_calib_post
  m_calib_post_eff[, 3:7] <- 1 - m_calib_post[, 3:7]
  
  m_samp_prior <- sample.prior(1000)
  m_samp_prior_eff <- m_samp_prior
  m_samp_prior_eff[, 3:7] <- 1-m_samp_prior_eff[, 3:7]
  
  # Ordered prior vs posterior distribution in terms of NPI effectiveness
  df_samp_prior_eff <- melt(cbind(PDF = "Prior", 
                                  as.data.frame(m_samp_prior_eff)), 
                            variable.name = "Parameter")
  df_samp_post_imis_eff  <- melt(cbind(PDF = "Posterior IMIS", 
                                       as.data.frame(m_calib_post_eff)), 
                                 variable.name = "Parameter")
  df_samp_prior_post_eff <- bind_rows(df_samp_prior_eff, 
                                      df_samp_post_imis_eff)
  df_samp_prior_post_eff$PDF <- ordered(df_samp_prior_post_eff$PDF, 
                                        levels = c("Prior", "Posterior IMIS")) # "Posterior SIR", 
  
  df_samp_prior_post_eff$Parameter <- factor(df_samp_prior_post_eff$Parameter,
                                             levels = levels(df_samp_prior_post_eff$Parameter),
                                             ordered = TRUE, 
                                             labels = c(expression(beta),
                                                        expression(tau),
                                                        expression(eta[1]),
                                                        expression(eta[2]),
                                                        expression(eta[3]),
                                                        expression(eta[4]),
                                                        expression(eta[5]),
                                                        expression(nu[lb]),
                                                        expression(nu[ub]),
                                                        expression(nu[rate]),
                                                        expression(nu[mid])
                                                     ))
  
  gg_post_imis <- ggplot(df_samp_prior_post_eff, 
                         aes(x = value, y = ..density.., fill = PDF)) +
    facet_wrap(~Parameter, scales = "free", 
               ncol = 4,
               labeller = label_parsed) +
    scale_x_continuous(breaks = number_ticks(4)) +
    geom_density(alpha=0.5) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom",
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  # Save plot
  ggsave(gg_post_imis,
         filename = paste0("figs/",n_dir_proj,"figS3_posterior_v_prior-IMIS-1k_",n_proj_type,".jpg"),
         width = 10, height = 8, dpi = 300)
  
}

