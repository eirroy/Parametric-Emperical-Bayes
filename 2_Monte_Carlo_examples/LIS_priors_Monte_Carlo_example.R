#######################################################################################################################
# LIS_priors_Monte_Carlo_example.R
#
# Demonstration of a laboratory information system (LIS) –based PEB approach for personal reference intervals (prRI):
# Establishing priors from simulated nested BV data with known parameters.
######################################################################################################################

# --- 0) Load/Install Required Packages --------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(refineR, robustbase, ggplot2, patchwork, dplyr)

# For reproducibility
set.seed(42)

##############################################################################
# 1) Set simulation parameters for establishing LIS-priors
##############################################################################

# Population (mean) parameter on the original scale
mu_pop <- 100 

# Between-subject variation (CV_g)
cv_between <- 0.20

# Within-subject variation (CV_i)
cv_within <- 0.10

# Analytical variation (half of cv_within)
cv_analytical <- 0.5 * cv_within  # (CV_a)

# Total variation (within-subject + analytical)
cv_total <- sqrt(cv_within^2 + cv_analytical^2)

# Significance level for 99% intervals
alpha <- 0.05
z_score_95 <- qnorm(1 - alpha / 2)

# Number of individuals and measurements
number_of_patients     <- 5000
samples_per_individual <- 5

# Between-subject SD on the original scale
sd_between <- mu_pop * cv_between

# Simulate each person's “homeostatic set point”
individual_means <- rnorm(
  number_of_patients,
  mean = mu_pop,
  sd   = sd_between
)

# Create a nested dataset (patient × repeated measurements)
sim_data <- data.frame(
  patient_id      = rep(seq_len(number_of_patients), each = samples_per_individual),
  individual_mean = rep(individual_means,            each = samples_per_individual),
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
# 2) Estimate priors with refineR + robust regression
##############################################################################

# Use patients 31+ to define priors; keep 1–30 for validation
df_first <- sim_data %>%
  filter(patient_id > 30) %>%
  group_by(patient_id) %>%
  slice_min(sample_number, with_ties = FALSE) %>%
  ungroup()

# Use refineR to obtain Box–Cox parameters
refineR_obj <- findRI(df_first$measurement)
lambda_boxcox        <- refineR_obj$Lambda
mu_population_Box_Cox    <- refineR_obj$Mu
sigma_population_Box_Cox <- refineR_obj$Sigma

# Exclude outliers (99% refineR interval) before robust regression
ref_range <- getRI(refineR_obj, RIperc = c(0.005, 0.995))
ref_lower <- ref_range[1, 2]
ref_upper <- ref_range[2, 2]

sim_data_filtered <- sim_data %>%
  filter(
    measurement >= ref_lower,
    measurement <= ref_upper,
    patient_id > 30
  ) %>%
  mutate(
    # Box–Cox transform each measurement
    Box_Cox_measurement = BoxCox(measurement, lambda_boxcox)
  )

# Prepare data for robust regression
paired_data <- sim_data_filtered %>%
  group_by(patient_id) %>%
  arrange(sample_number) %>%
  mutate(Box_Cox_measurement_prev = lag(Box_Cox_measurement)) %>%
  filter(!is.na(Box_Cox_measurement_prev)) %>%
  ungroup()

# Fit a robust linear model to estimate b1
robust_fit <- lmrob(Box_Cox_measurement ~ Box_Cox_measurement_prev, data = paired_data)
b1 <- coef(robust_fit)[2] # Extract regression slope

##############################################################################
# 3) Iterative PEB: dynamic prediction intervals
##############################################################################
# For each patient, compute the predicted set point and prediction intervals,
# iterating from 0 to n prior measurements.

get_intervals_for_patient <- function(d) {
  d <- d[order(d$sample_number), ]
  x <- d$Box_Cox_measurement
  n <- length(x)
  
  # n=0: no prior samples, use population parameters
  df0 <- data.frame(
    patient_id              = d$patient_id[1],
    number_of_prior_samples = 0,
    Box_Cox_measurement         = x[1],
    Box_Cox_predicted_set_point = mu_population_Box_Cox,
    Box_Cox_lower_interval      = mu_population_Box_Cox - z_score_95 * sigma_population_Box_Cox,
    Box_Cox_upper_interval      = mu_population_Box_Cox + z_score_95 * sigma_population_Box_Cox
  )
  
  # For prior samples n >= 1
  ns      <- 1:n
  bn_vec  <- (b1 * ns) / (b1 * ns + (1 - b1))
  mean_bc <- cumsum(x) / ns
  
  pred_set     <- mu_population_Box_Cox + (mean_bc - mu_population_Box_Cox) * bn_vec
  peb_threshold <- z_score_95 * sqrt(1 - (b1 * bn_vec)) * sigma_population_Box_Cox
  
  # Next measurement is the (n+1)-th element (NA for last)
  next_measurement <- c(x[-1], NA)
  
  df1 <- data.frame(
    patient_id              = d$patient_id[1],
    number_of_prior_samples = ns,
    Box_Cox_measurement         = next_measurement,
    Box_Cox_predicted_set_point = pred_set,
    Box_Cox_lower_interval      = pred_set - peb_threshold,
    Box_Cox_upper_interval      = pred_set + peb_threshold
  )
  
  rbind(df0, df1)
}

##############################################################################
# 4) Apply the PEB and assemble final results
##############################################################################

final_res <- sim_data %>%
  mutate(Box_Cox_measurement = BoxCox(measurement, lambda_boxcox)) %>%
  split(.$patient_id) %>%
  lapply(get_intervals_for_patient) %>%
  bind_rows()

##############################################################################
# 5) Back-transform & identify out-of-range measurements
##############################################################################

final_res <- final_res %>%
  mutate(
    measurement_original = ifelse(
      !is.na(Box_Cox_measurement),
      invBoxCox(Box_Cox_measurement, lambda_boxcox),
      NA_real_
    ),
    predicted_original = invBoxCox(Box_Cox_predicted_set_point, lambda_boxcox),
    lower_original     = invBoxCox(Box_Cox_lower_interval,      lambda_boxcox),
    upper_original     = invBoxCox(Box_Cox_upper_interval,      lambda_boxcox),
    peb_threshold = case_when(
      !is.na(measurement_original) &
        (measurement_original < lower_original |
           measurement_original > upper_original) ~ "outside",
      !is.na(measurement_original) ~ "inside",
      TRUE ~ NA_character_
    )
  )

##############################################################################
# 6) Example plot: personal reference interval for first 30 patients
##############################################################################

p <- final_res %>%
  filter(patient_id <= 30) %>%
  group_by(patient_id) %>%
  slice(1:5) %>%
  ggplot(aes(x = number_of_prior_samples)) +
  geom_ribbon(
    aes(ymin = lower_original, ymax = upper_original),
    fill = "lightblue", alpha = 0.4
  ) +
  geom_line(aes(y = predicted_original), color = "steelblue4") +
  geom_point(aes(y = measurement_original, color = peb_threshold)) +
  facet_wrap(~ patient_id, scales = "free_y") +
  scale_color_manual(values = c("inside" = "steelblue1", "outside" = "red")) +
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

