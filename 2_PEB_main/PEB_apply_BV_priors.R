##############################################################################
# PEB_monte_carlo_BV_Tidy.R
#
# Demonstration of a PEB approach using known Biological Variation (BV) priors.
# - log transform
# - uses known CVI, CVA, CVG
##############################################################################

# 0) Load/Install Required Packages ------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(here, ggplot2, dplyr)

##############################################################################
# 1) BV Parameters & Setup
##############################################################################
CVI    <- 0.10       # within-subject CV
CVA    <- 0.5 * CVI  # analytical CV
CVG    <- 0.20       # between-subject CV
pop_mu <- 100        # population mean on the original scale

# Derived BV parameters (on log scale)
cv_total   <- sqrt(CVI^2 + CVA^2)
sigma_within_log  <- sqrt(log(1 + cv_total^2))
sigma_between_log <- sqrt(log(1 + CVG^2))
sigma_pop <- sqrt(sigma_within_log^2 + sigma_between_log^2)  # total pop SD on log scale
B1        <- sigma_between_log^2 / (sigma_between_log^2 + sigma_within_log^2)
mu_pop    <- log(pop_mu)  # population mean on log scale

# Significance level
alpha   <- 0.05
z_score <- qnorm(1 - alpha / 2)

##############################################################################
# 2) Load and Prepare Data
##############################################################################
data <- read.csv(here("2_PEB_main", "PEB_data.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    measurement = as.numeric(measurement)
  )

# Sort by patient & sample number, convert patient_id to factor
data <- data %>%
  arrange(patient_id, sample_number) %>%
  mutate(patient_id = factor(patient_id))

##############################################################################
# 3) Define Log Transform
##############################################################################
log_transform <- function(x) {
  log(x)
}

inv_log_transform <- function(y) {
  exp(y)
}

# Apply the log transform
data <- data %>%
  mutate(transformed = log_transform(measurement))

##############################################################################
# 4) Define Iterative PEB Function
##############################################################################
PEB_function <- function(patient_data) {
  
  # Ensure ascending order by sample_number
  patient_data <- patient_data[order(patient_data$sample_number), ]
  n_rows       <- nrow(patient_data)
  
  out_list <- vector("list", length = n_rows + 1)
  
  for (n_obs in seq_len(n_rows + 1) - 1) {
    
    # If n_obs == 0, use the population-level parameters
    if (n_obs == 0) {
      y_hat         <- mu_pop
      threshold     <- z_score * sigma_pop
      next_meas     <- if (n_obs < n_rows) patient_data$transformed[n_obs + 1] else NA
    } else {
      # For n_obs >= 1, apply Empirical Bayes update
      prev_meas  <- patient_data$transformed[seq_len(n_obs)]
      m_bar      <- mean(prev_meas)
      next_meas  <- if (n_obs < n_rows) patient_data$transformed[n_obs + 1] else NA
      
      # PEB formulas
      Bn    <- (B1 * n_obs) / (B1 * n_obs + (1 - B1))
      y_hat <- mu_pop + (m_bar - mu_pop) * Bn
      threshold <- z_score * sqrt(1 - B1 * Bn) * sigma_pop
    }
    
    out_list[[n_obs + 1]] <- data.frame(
      patient_id            = patient_data$patient_id[1],
      n_prior               = n_obs,
      transform_measurement = next_meas,
      transform_yhat        = y_hat,
      transform_lower       = y_hat - threshold,
      transform_upper       = y_hat + threshold
    )
  }
  
  bind_rows(out_list)
}

##############################################################################
# 5) Apply PEB to Each Patient (Generate Personal Reference Intervals)
##############################################################################
prRI <- data %>%
  split(.$patient_id) %>%
  lapply(PEB_function) %>%
  bind_rows()

##############################################################################
# 6) Back-Transform & Flag Out-of-Range
##############################################################################
prRI <- prRI %>%
  mutate(
    measurement = if_else(
      is.na(transform_measurement),
      NA_real_,
      inv_log_transform(transform_measurement)
    ),
    yhat  = inv_log_transform(transform_yhat),
    lower = inv_log_transform(transform_lower),
    upper = inv_log_transform(transform_upper),
    PEB_flag = case_when(
      !is.na(measurement) & (measurement < lower | measurement > upper) ~ "Outside",
      !is.na(measurement)                                              ~ "Inside",
      TRUE                                                             ~ NA_character_
    )
  )

##############################################################################
# 7) Personal Reference Interval Plot (n > 1)
##############################################################################
# Example colors; adjust as desired
fill_color   <- "#8cce91"
line_color   <- "darkgreen"
inside_color <- "forestgreen"

prRI_plot <- prRI %>%
  group_by(patient_id) %>%
  ggplot(aes(x = n_prior)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = fill_color, alpha = 0.3) +
  geom_line(aes(y = yhat), color = line_color) +
  geom_point(aes(y = measurement, color = PEB_flag), na.rm = TRUE) +
  facet_wrap(~ patient_id, scales = "free_y") +
  scale_color_manual(values = c("Inside" = inside_color, "Outside" = "red")) +
  labs(
    x = "Number of Prior Measurements (n)",
    y = "Back-Transformed Values",
    color = "PEB Flag"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

print(prRI_plot)

##############################################################################
# 8) Calculate Result Pairs for Reference Change Values
##############################################################################
result_pairs <- data %>%
  group_by(patient_id) %>%
  arrange(sample_number) %>%
  mutate(
    Prev_Sample    = lag(measurement),
    Current_Sample = measurement
  ) %>%
  filter(!is.na(Prev_Sample)) %>%
  ungroup()

##############################################################################
# 9) Reference Change Value Threshold Plot (n = 1)
##############################################################################
result_pairs <- result_pairs %>%
  mutate(
    log_prev  = log_transform(Prev_Sample),
    log_yhat  = mu_pop + (log_prev - mu_pop) * B1,  # Single-step
    halfw     = z_score * sqrt(1 - B1^2) * sigma_pop,
    RCV_lower = inv_log_transform(log_yhat - halfw),
    RCV_upper = inv_log_transform(log_yhat + halfw),
    RCV_flag  = if_else(
      Current_Sample < RCV_lower | Current_Sample > RCV_upper,
      "Outside",
      "Inside"
    )
  )

RCV_plot <- ggplot(result_pairs, aes(x = Prev_Sample, y = Current_Sample)) +
  geom_ribbon(aes(ymin = RCV_lower, ymax = RCV_upper), 
              fill = "#8cce91", alpha = 0.2) +
  geom_point(aes(color = RCV_flag, shape = RCV_flag), size = 2.5) +
  scale_color_manual(name = "Threshold",
                     values = c("Inside" = "forestgreen", "Outside" = "red")) +
  scale_shape_manual(name = "Threshold",
                     values = c("Inside" = 16, "Outside" = 17)) +
  labs(
    title = "Reference Change Value Threshold: BV Priors",
    x = "Previous Result",
    y = "Current Sample"
  ) +
  theme_minimal()

print(RCV_plot)
