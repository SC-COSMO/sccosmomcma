library(dampack)
# Grid of X-axis values
x <- seq(0, 20, 0.1)

mean_time <- 7

lambda_exp <- 1/mean_time
sx_gamma_params <- gamma_params(mu = mean_time, sigma = 2.9, scale = F)

pdf("figs/01_exponential_vs_gamma.pdf", width = 8, height = 6)
# Exponential
plot(x, pexp(x, lambda_exp), type = "l",
     ylab = "F(x)", lwd = 2, col = "#8C1515", 
     ylim = c(0, 1), 
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
# Gamma
lines(x, pgamma(x, shape = sx_gamma_params$shape, rate = sx_gamma_params$rate), 
      col = "#1D5949", lty = 1, lwd = 2)
# 1 assymptote
abline(h = 1, lty = 2)
# Adding a legend
legend("right", 
       c(as.expression(bquote(bold("Dwell time"))), 
         "Exponential (constant rate)", "Gamma (non-constant rate)"),
       lty = c(0, 1, 1), col = c(NA, "#8C1515", "#1D5949"), 
       box.lty = 0, lwd = 2, cex = 1.5)
dev.off()
