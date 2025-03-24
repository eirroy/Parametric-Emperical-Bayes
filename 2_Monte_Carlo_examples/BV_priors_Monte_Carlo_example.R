##############################################################################
# PEB_monte_carlo_example.R
#
# Demonstration of a biological-variation (BV)–based PEB approach for
# personal reference intervals, using simulated data with known BV parameters.
# Thresholds begin at the population reference interval (popRI).
##############################################################################

# --- 0) Load/Install Required Packages --------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(ggplot2, patchwork, dplyr, purrr)

# For reproducibility
set.seed(42)

##############################################################################
# 1) Set simulation parameters for the BV-based PEB approach
##############################################################################

# Population (mean) parameter on the original scale
population_mean <- 100  # (pop_mu)

# Between-subject biological variation
cv_between <- 0.20      # (CV_g)

# Within-subject biological variation
cv_within  <- 0.10      # (CV_i)

# Analytical variation (half of cv_within)
cv_analytical <- 0.5 * cv_within  # (CV_a)

# Total variation (within-subject + analytical)
cv_total <- sqrt(cv_within^2 + cv_analytical^2)

# Number of individuals and measurements
number_of_patients   <- 5000
samples_per_individual <- 5

# Significance level for 95% intervals
alpha      <- 0.05
z_score_95 <- qnorm(1 - alpha / 2)

# Simulate each person's “homeostatic set point” from a normal distribution
individual_means <- rnorm(
  number_of_patients, 
  mean = population_mean, 
  sd   = population_mean * cv_between
)

# Create a nested dataset (patient × repeated measurements)
bv_data <- data.frame(
  patient_id      = rep(seq_len(number_of_patients), each = samples_per_individual),
  individual_mean = rep(individual_means, each = samples_per_individual),
  sample_number   = rep(seq_len(samples_per_individual), times = number_of_patients)
) %>%
  mutate(
    # Simulate measurements around each individual's mean, 
    # with variation ~ (individual_mean * cv_total).
    measurement = rnorm(
      n(), 
      mean = individual_mean, 
      sd   = individual_mean * cv_total
    )
  )

##############################################################################
# 2) Define priors using known coefficients of variation (log scale)
##############################################################################

# Convert CVs to standard deviations on the log scale
sigma_total_log   <- sqrt(log(1 + cv_total^2)) 
sigma_between_log <- sqrt(log(1 + cv_between^2))

# Population standard deviation on the log scale
sigma_population_log <- sqrt(sigma_total_log^2 + sigma_between_log^2)

# Population mean on the log scale
mu_population_log <- log(population_mean)

# Intraclass correlation on the log scale
b1 <- sigma_between_log^2 / (sigma_between_log^2 + sigma_total_log^2)

##############################################################################
# 3) Log-transform the simulated data
##############################################################################
bv_data <- bv_data %>%
  rename(
    patient_id_new   = patient_id, 
    concentration    = measurement
  ) %>%
  mutate(
    log_concentration = log(concentration)
  )

##############################################################################
# 4) Iterative PEB: dynamic prediction intervals
##############################################################################
# For each patient, compute the predicted homeostatic set point and 
# prediction intervals, iterating from 0 to n prior measurements.

get_intervals_for_patient <- function(df) {
  # Ensure chronological order
  df <- df[order(df$sample_number), ]
  x <- df$log_concentration
  n <- length(x)
  
  # n = 0: No prior samples, use population parameters (popRI)
  df0 <- data.frame(
    patient_id = df$patient_id_new[1],
    number_of_prior_samples = 0,
    log_measurement = x[1],
    log_predicted_set_point = mu_population_log,
    log_lower_interval = mu_population_log - z_score_95 * sigma_population_log,
    log_upper_interval = mu_population_log + z_score_95 * sigma_population_log
  )
  
  # For prior samples n >= 1: 
  ns <- 1:n
  bn <- (b1 * ns) / (b1 * ns + (1 - b1))
  sample_mean <- cumsum(x) / ns
  pred_set <- mu_population_log + (sample_mean - mu_population_log) * bn
  peb_threshold <- z_score_95 * sqrt(1 - (b1 * bn)) * sigma_population_log
  
  # Next measurement (for prediction) is the (n+1)th element (NA for last)
  next_measurement <- c(x[-1], NA)
  
  df1 <- data.frame(
    patient_id = df$patient_id_new[1],
    number_of_prior_samples = ns,
    log_measurement = next_measurement,
    log_predicted_set_point = pred_set,
    log_lower_interval = pred_set - peb_threshold,
    log_upper_interval = pred_set + peb_threshold
  )
  
  rbind(df0, df1)
}

# Apply the function to each patient
final_res <- bv_data %>%
  split(.$patient_id_new) %>%
  map_dfr(get_intervals_for_patient)

##############################################################################
# 5) Back-transform & identify samples exceeding the PEB thresholds
##############################################################################
population_lower <- exp(mu_population_log - z_score_95 * sigma_population_log)
population_upper <- exp(mu_population_log + z_score_95 * sigma_population_log)

final_res <- final_res %>%
  mutate(
    measurement_original = ifelse(!is.na(log_measurement), exp(log_measurement), NA),
    predicted_original   = exp(log_predicted_set_point),
    lower_original       = exp(log_lower_interval),
    upper_original       = exp(log_upper_interval),
    peb_threshold_flag = case_when(
      !is.na(measurement_original) & 
        (measurement_original < lower_original | measurement_original > upper_original)
      ~ "outside",
      !is.na(measurement_original) ~ "inside",
      TRUE                         ~ NA_character_
    )
  )

##############################################################################
# 6) Example plot for a subset of patients
##############################################################################
# Show the first 30 patients, first 5 dynamic intervals each.

p <- final_res %>%
  filter(patient_id <= 30) %>%            
  group_by(patient_id) %>%
  slice(1:5) %>%                          
  ggplot(aes(x = number_of_prior_samples)) +
  geom_ribbon(
    aes(ymin = lower_original, ymax = upper_original),
    fill = "lightgreen", alpha = 0.5
  ) +
  geom_line(aes(y = predicted_original), color = "darkgreen") +
  geom_point(aes(y = measurement_original, color = peb_threshold_flag)) +
  facet_wrap(~ patient_id, scales = "free_y") +
  scale_color_manual(values = c("inside" = "forestgreen", "outside" = "red")) +
  labs(
    x = "Number of Prior Measurements (n)", 
    y = "PEB thresholds"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )

print(p)
