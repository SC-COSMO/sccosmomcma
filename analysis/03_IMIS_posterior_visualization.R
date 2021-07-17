#------------------------------------------------------------------------------#
# This script plots the calibrated posterior distributions.                    #
#                                                                              #
# Author:                                                                      #
#    Fernando Alarid-Escudero, PhD, <fernando.alarid@cide.edu>                 #
#                                                                              #
#------------------------------------------------------------------------------#
rm(list = ls())

# Load libraries and data -------------------------------------------------

library(ggplot2)
library(ggcorrplot)
library(GGally)      # To use ggpairs
library(dampack)

# Calibrated parameters
load("output/03_map_output_IMIS_MCMA_SQ_2020-12-13.RData")


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

stricter_npi <- matrixStats::rowMins(m_calib_post[, 3:7])
mean(stricter_npi)
quantile(stricter_npi, probs = c(0.025, 0.975))

## Effectiveness
# Right after NPI implementation in March 23
mean(m_calib_post_eff[, "r_soc_dist_factor"])
quantile(m_calib_post_eff[, "r_soc_dist_factor"], probs = c(0.025, 0.975))
# At the end of calibration period
mean(m_calib_post_eff[, "r_soc_dist_factor_5"])
quantile(m_calib_post_eff[, "r_soc_dist_factor_5"], probs = c(0.025, 0.975))

holiday_bump <- m_calib_post_eff[, "r_soc_dist_factor_5"] - 0.3
mean(holiday_bump)
quantile(holiday_bump, probs = c(0.025, 0.975))

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
gg_post_pairs_corr

# Save plot
ggsave(filename = "figs/figs_paper/figS2_posterior-IMIS-correlation-1k.jpg",
       width = 12, height = 8, dpi = 300)


## Posterior vs prior comparison ------------------------------------------

df_samp_prior_post_ordered_eff$Parameter <- factor(df_samp_prior_post_ordered_eff$Parameter, 
                                                   levels = levels(df_samp_prior_post_ordered_eff$Parameter),
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

gg_post_imis <- ggplot(df_samp_prior_post_ordered_eff, 
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
       filename = "figs/figs_paper/figS3_posterior_v_prior-IMIS-1k.jpg",
       width = 10, height = 8, dpi = 300)
