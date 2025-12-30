

library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(42)   # Ensures reproducibility

# -----------------------
# STEP 1: Basic Measurements & Histogram
# -----------------------
true_value <- 100
noise_sd <- 5
n_measurements <- 100

measurements <- rnorm(n_measurements, mean = true_value, sd = noise_sd)
data <- data.frame(measurements)

# Plot 1: Distribution of measurements
ggplot(data, aes(x = measurements)) +
  geom_histogram(binwidth = 2, fill = "steelblue", color = "black") +
  geom_vline(xintercept = true_value, color = "red", linewidth = 1) +
  labs(title = "Distribution of Measurements",
       x = "Measured Value",
       y = "Frequency")

# STEP 2: Mean, SD, Standard Error
mean_value <- mean(measurements)
std_dev <- sd(measurements)
standard_error <- std_dev / sqrt(n_measurements)

# STEP 3: Confidence Interval
confidence_level <- 0.95
alpha <- 1 - confidence_level
z <- qnorm(1 - alpha/2)
lower_ci <- mean_value - z * standard_error
upper_ci <- mean_value + z * standard_error

# Plot 2: CI for measurement
ci_data <- data.frame(
  mean = mean_value,
  lower = lower_ci,
  upper = upper_ci
)

ggplot(ci_data, aes(x = 1, y = mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(title = "Estimated Value with 95% Confidence Interval",
       x = "",
       y = "Measured Quantity") +
  theme_minimal()

# -----------------------
# STEP 4: Monte Carlo Uncertainty Propagation
# -----------------------
V_mean <- 10
V_sd   <- 0.2
I_mean <- 2
I_sd   <- 0.05
N <- 10000

V_samples <- rnorm(N, mean = V_mean, sd = V_sd)
I_samples <- rnorm(N, mean = I_mean, sd = I_sd)

P_samples <- V_samples * I_samples
P_mean <- mean(P_samples)
P_sd   <- sd(P_samples)
CI_95 <- quantile(P_samples, probs = c(0.025, 0.975))

# Plot 3: Monte Carlo histogram of power
hist(
  P_samples,
  breaks = 60,
  probability = TRUE,
  main = "Monte Carlo Uncertainty Propagation for Power",
  xlab = "Power (W)",
  col = "lightblue",
  border = "white"
)
abline(v = P_mean, col = "red", lwd = 2)
abline(v = CI_95[1], col = "darkgreen", lwd = 2, lty = 2)
abline(v = CI_95[2], col = "darkgreen", lwd = 2, lty = 2)
legend(
  "topright",
  legend = c("Mean", "95% CI"),
  col = c("red", "darkgreen"),
  lwd = 2,
  lty = c(1, 2)
)

# -----------------------
# STEP 5: Bayesian Estimation
# -----------------------
set.seed(1)
n <- length(measurements)
sample_mean <- mean(measurements)

# Prior
prior_mean <- 95
prior_sd <- 10

posterior_variance <- 1 / ( (1/prior_sd^2) + (n/noise_sd^2) )
posterior_mean <- posterior_variance * ( (prior_mean/prior_sd^2) + (n * sample_mean / noise_sd^2) )
posterior_sd <- sqrt(posterior_variance)

x <- seq(70, 130, length.out = 1000)
prior_density <- dnorm(x, prior_mean, prior_sd)
posterior_density <- dnorm(x, posterior_mean, posterior_sd)

# Plot 4: Bayesian prior + posterior + measurements
hist(measurements, probability = TRUE, col = "lightgray",
     main = "Bayesian Estimation of Physical Measurement",
     xlab = "Measured Quantity")
lines(x, prior_density, col = "blue", lwd = 2)
lines(x, posterior_density, col = "red", lwd = 2)
abline(v = true_value, col = "darkgreen", lwd = 2, lty = 2)
legend("topright",
       legend = c("Measurements", "Prior", "Posterior", "True Value"),
       col = c("lightgray", "blue", "red", "darkgreen"),
       lwd = c(5, 2, 2, 2))

# Plot 5: Frequentist vs Bayesian
freq_mean <- sample_mean
plot(x, posterior_density, type = "l", col = "red", lwd = 2,
     main = "Frequentist vs Bayesian Estimation",
     xlab = "Measured Quantity", ylab = "Density")
abline(v = freq_mean, col = "blue", lwd = 2)
abline(v = true_value, col = "darkgreen", lwd = 2, lty = 2)
legend("topright",
       legend = c("Bayesian Posterior", "Frequentist Mean", "True Value"),
       col = c("red", "blue", "darkgreen"),
       lwd = 2)

# -----------------------
# STEP 6: Advanced Bayesian analysis for Unknown Noise
# -----------------------
alpha_prior <- 3
beta_prior  <- 20

alpha_post <- alpha_prior + n/2
beta_post  <- beta_prior + sum((measurements - sample_mean)^2)/2

posterior_var_mean <- beta_post / (alpha_post - 1)
posterior_sd_est   <- sqrt(posterior_var_mean)

# Plot 6: Posterior distribution of noise variance
x_var <- seq(1, 100, length.out = 1000)
posterior_var_density <- dgamma(1/x_var, shape = alpha_post, rate = beta_post) * (1 / x_var^2)
plot(x_var, posterior_var_density, type = "l", col = "purple", lwd = 2,
     main = "Posterior Distribution of Unknown Noise Variance",
     xlab = "Variance", ylab = "Density")
abline(v = posterior_var_mean, col = "darkorange", lwd = 2, lty = 2)
legend("topright",
       legend = c("Posterior Density", "Posterior Mean"),
       col = c("purple", "darkorange"),
       lwd = 2, lty = c(1,2))

# Plot 7: Bayesian posterior with estimated noise
posterior_density_twist <- dnorm(x, posterior_mean, posterior_sd_est)
plot(x, posterior_density_twist, type = "l", col = "red", lwd = 2,
     main = "Bayesian Posterior with Unknown Noise Estimation",
     xlab = "Measured Quantity", ylab = "Density")
abline(v = true_value, col = "darkgreen", lwd = 2, lty = 2)
legend("topright",
       legend = c("Posterior (unknown noise)", "True Value"),
       col = c("red", "darkgreen"),
       lwd = 2)

# Plot 8: Comparison of Bayesian posteriors
plot(x, posterior_density, type = "l", col = "blue", lwd = 2,
     main = "Comparison of Bayesian Posteriors",
     xlab = "Measured Quantity", ylab = "Density")
lines(x, posterior_density_twist, col = "red", lwd = 2, lty = 2)
abline(v = true_value, col = "darkgreen", lwd = 2, lty = 2)
legend("topright",
       legend = c("Posterior (known noise)", "Posterior (unknown noise)", "True Value"),
       col = c("blue", "red", "darkgreen"),
       lwd = 2, lty = c(1,2,2))

# Final summary of statistics
posterior_mean
posterior_sd
posterior_sd_est
freq_mean
lower_ci
upper_ci
P_mean
P_sd
