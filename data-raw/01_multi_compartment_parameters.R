library(dampack)
library(kdensity)

n_samples <- 10000

# Source: He et al. Nature paper
latent_params <- gamma_params(mu = 3, sigma = 1, scale = F)
v_latent_duration <- rgamma(n_samples, 
                            shape = latent_params$shape, 
                            rate  = latent_params$rate)

infectious_params <- gamma_params(mu = 3.1, sigma = 2.1, scale = F)
v_infectious_duration <- rgamma(n_samples, 
                                shape = infectious_params$shape, 
                                rate  = infectious_params$rate)

## Source: 
# Lauer, S. A., Grantz, K. H., Bi, Q., Jones, F. K., Zheng, Q., Meredith, H. R., 
# Azman, A. S., Reich, N. G., & Lessler, J. (2020). The Incubation Period of 
# Coronavirus Disease 2019 (COVID-19) From Publicly Reported Confirmed Cases: 
# Estimation and Application. Annals of Internal Medicine, 2019. https://doi.org/10.7326/m20-0504
time_to_symptoms_params <- gamma_params(mu = 5, sigma = 2.9, scale = F) 
v_time_to_symptoms <- rgamma(n_samples, 
                             shape = time_to_symptoms_params$shape, 
                             rate  = time_to_symptoms_params$rate)

mean(v_latent_duration)
sd(v_latent_duration)
quantile(v_latent_duration, 
         probs = c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99))

mean(v_infectious_duration)
sd(v_infectious_duration)
quantile(v_infectious_duration, 
         probs = c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99))

mean(v_time_to_symptoms)
sd(v_time_to_symptoms)
quantile(v_time_to_symptoms, 
         probs = c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99))

# psych::pairs.panels(cbind(v_latent_duration, v_time_to_symptoms))

### Validation to: 
## - Ashcroft, P., Huisman, J. S., Lehtinen, S., Bouman, J. A., Althaus, C. L., 
##   Regoes, R. R., & Bonhoeffer, S. (2020). COVID-19 infectivity profile 
##   correction. Swiss Medical Weekly, August, 3–5. 
##   https://doi.org/10.4414/smw.2020.20336
df_ashcroft <- read.csv(file = "data-raw/WebPlotDigitizer_Ashcroft_2020.csv")
plot(CDF/100~days, data = df_ashcroft)
v_time_symptoms_infectiousness <- v_latent_duration - v_time_to_symptoms
# plot(density(v_time_symptoms_infectiousness))
plot(ecdf(v_time_symptoms_infectiousness), add = TRUE, col = "red") 
abline(v = -1)
abline(v = 0)
abline(v = 0.5)
abline(v = 1.5)
abline(h = 0.5)
abline(v = -2.5)
abline(h = 0.12)

abline(h = 0.9)
abline(v = 1.9)
abline(v = 9.3)

abline(h = 0.05)
abline(h = 0.34)
abline(v = -2.0)


sum(v_time_symptoms_infectiousness < -2.5)/n_samples

### Determine number of compartments
v_infection_duration <- (v_latent_duration + v_infectious_duration)
v_infection_duration_gt_symptom_time_ind <- v_infection_duration > v_time_to_symptoms
v_infection_duration_gt_symptom_time <- v_time_to_symptoms[v_infection_duration_gt_symptom_time_ind]
plot(ecdf(v_infection_duration_gt_symptom_time), add = TRUE, col = "blue") 

mean(v_time_to_symptoms > (v_latent_duration + v_infectious_duration))

mean(v_time_to_symptoms > (v_latent_duration) & v_time_to_symptoms < (v_latent_duration + v_infectious_duration))
mean(v_time_to_symptoms < (v_latent_duration))


# take the mean from the rgamma vector: the duration
# our rate parameter is equal to 1/duration
# find number of compartments s.t. dgamma density matches the rgamma vector

number_compartments <- function(rate_of_exit, l_params) {
  
  max_number_compartments <- floor(1/rate_of_exit)
  v_times <- seq(0,30, length.out = 1000)
  v_sse <- vector(length = max_number_compartments)
  for (i in 1:max_number_compartments) {
    out_gamma <- dgamma(v_times, shape=i, rate = (rate_of_exit * i))
    target_gamma <- dgamma(v_times, shape=l_params$shape, rate=l_params$rate)
    sse <- sum((out_gamma - target_gamma)^2)  
    v_sse[i] <- sse
  }
  v_sse  
  return(which.min(v_sse))
}

number_compartments(1/mean(v_time_to_symptoms), time_to_symptoms_params)

n_compartments_E <- number_compartments(1/mean(v_latent_duration), latent_params)
n_compartments_I <- number_compartments(1/mean(v_infectious_duration), infectious_params)

# get the daily increase in the likelihood of symptoms for the first n_compartments_E +
# n_compartments_I days 
for (i in 1:(n_compartments_E + n_compartments_I)) {
  #  x = seq(i-1, i)
  #  plot(pgamma(q=x, shape=time_to_symptoms_params$shape, rate=time_to_symptoms_params$rate))
  print(i)
  prev_frac <- pgamma(q=i-1, shape=time_to_symptoms_params$shape, rate=time_to_symptoms_params$rate)
  end_frac <- pgamma(q=i, shape=time_to_symptoms_params$shape, rate=time_to_symptoms_params$rate)
  frac_detected <- (end_frac-prev_frac)/(1-prev_frac)
  print(prev_frac) 
  print(end_frac) 
  print(frac_detected) 
}

vec <- numeric((4*(n_compartments_E + n_compartments_I)))

for (i in 1:(4*(n_compartments_E + n_compartments_I))) {
  #  x = seq(i-1, i)
  #  plot(pgamma(q=x, shape=time_to_symptoms_params$shape, rate=time_to_symptoms_params$rate))
  print(i)
  prev_frac <- pgamma(q=i-1, shape=time_to_symptoms_params$shape, rate=time_to_symptoms_params$rate)
  end_frac <- pgamma(q=i, shape=time_to_symptoms_params$shape, rate=time_to_symptoms_params$rate)
  frac_detected <- (end_frac-prev_frac)/(1-prev_frac)
  vec[i] <- frac_detected
  print(prev_frac) 
  print(end_frac) 
  print(frac_detected) 
}
plot(vec[1:5], ylim = c(0, max(vec[1:5])))
points(vec5, col = "red")

rate_time_to_symptoms <- 1/mean(v_time_to_symptoms)
max_number_compartments <- floor(mean(v_time_to_symptoms))
v_times <- seq(0,30, length.out = 1000)
v_sse <- vector(length = max_number_compartments)
for (i in 1:max_number_compartments) {
  out_gamma <- dgamma(v_times, shape=i, rate = (rate_time_to_symptoms * i))
  symptoms_gamma <- dgamma(v_times, shape=time_to_symptoms_params$shape, rate=time_to_symptoms_params$rate)
  sse <- sum((out_gamma - symptoms_gamma)^2)  
  v_sse[i] <- sse
}
v_sse
which.min(v_sse)

rate_inf <- 1/mean(v_infectious_duration)
max_number_compartments <- floor(mean(v_infectious_duration))
v_times <- seq(0,30, length.out = 1000)
v_sse <- vector(length = max_number_compartments)
for (i in 1:max_number_compartments) {
  out_gamma <- dgamma(v_times, shape=i, rate = (rate_inf * i))
  inf_gamma <- dgamma(v_times, shape=infectious_params$shape, rate=infectious_params$rate)
  sse <- sum((out_gamma - inf_gamma)^2)  
  v_sse[i] <- sse
}
v_sse
which.min(v_sse)

plot(v_times, dgamma(v_times, shape=1, rate = (rate_inf * 1)))
points(v_times, dgamma(v_times, shape=infectious_params$shape, rate=infectious_params$rate), col="red")
