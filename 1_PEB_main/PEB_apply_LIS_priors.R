##############################################################################
# PEB_apply_LIS_priors
#
# Application of PEB approach using priors established for LIS data.
##############################################################################

# 0) Load/Install Required Packages ------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(here, ggplot2, dplyr)

##############################################################################
# 1) Specify LIS-priors 
##############################################################################
lambda    <- 0.655447  # Box-Cox transformation parameter
sigma_pop <- 4.653029  # Population SD on the Box–Cox scale
mu_pop    <- 29.42723  # Population mean on the Box–Cox scale
B1        <- 0.7533967 # Intraclass correlation on the Box–Cox scale

# Set significance level
alpha   <- 0.05
z_score <- qnorm(1 - alpha / 2)

##############################################################################
# 2) Load and Prepare Data
##############################################################################
data <- read.csv(here("1_PEB_main", "PEB_apply_data.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    measurement = as.numeric(measurement)
  )

# Sort by patient & sample number, convert patient_id to factor
data <- data %>%
  arrange(patient_id, sample_number) %>%
  mutate(patient_id = factor(patient_id))

##############################################################################
# 3) Define Box–Cox Transform
##############################################################################
box_cox <- function(x, lambda) {
  if (lambda == 0) {
    log(x)
  } else {
    (x^lambda - 1) / lambda
  }
}

inv_box_cox <- function(y, lambda) {
  if (lambda == 0) {
    exp(y)
  } else {
    ((lambda * y) + 1)^(1 / lambda)
  }
}

# Apply the Box–Cox transform
data <- data %>%
  mutate(transformed = box_cox(measurement, lambda))

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
      Bn         <- (B1 * n_obs) / (B1 * n_obs + (1 - B1))
      y_hat      <- mu_pop + (m_bar - mu_pop) * Bn
      threshold  <- z_score * sqrt(1 - B1 * Bn) * sigma_pop
    }
    
    out_list[[n_obs + 1]] <- data.frame(
      patient_id             = patient_data$patient_id[1],
      n_prior                = n_obs,
      transform_measurement  = next_meas,
      transform_yhat         = y_hat,
      transform_lower        = y_hat - threshold,
      transform_upper        = y_hat + threshold
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
      inv_box_cox(transform_measurement, lambda)
    ),
    yhat   = inv_box_cox(transform_yhat, lambda),
    lower  = inv_box_cox(transform_lower, lambda),
    upper  = inv_box_cox(transform_upper, lambda),
    PEB_flag = case_when(
      !is.na(measurement) & (measurement < lower | measurement > upper) ~ "Outside",
      !is.na(measurement)                                              ~ "Inside",
      TRUE                                                             ~ NA_character_
    )
  )

##############################################################################
# 7) Personal Reference Interval Plot (n > 1)
##############################################################################
fill_color   <- "skyblue3"
line_color   <- "steelblue4"
inside_color <- "deepskyblue4"

prRI_plot <- prRI %>%
  # filter(patient_id %in% unique(patient_id)[1:20]) %>%  # Uncomment to restrict plot to a selected range of patient_ids
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
    bc_prev  = box_cox(Prev_Sample, lambda),
    bc_yhat  = mu_pop + (bc_prev - mu_pop) * B1,
    halfw    = z_score * sqrt(1 - B1^2) * sigma_pop,
    RCV_lower = inv_box_cox(bc_yhat - halfw, lambda),
    RCV_upper = inv_box_cox(bc_yhat + halfw, lambda),
    RCV_flag  = if_else(
      Current_Sample < RCV_lower | Current_Sample > RCV_upper,
      "Outside",
      "Inside"
    )
  )

RCV_plot <- ggplot(result_pairs, aes(x = Prev_Sample, y = Current_Sample)) +
  geom_ribbon(aes(ymin = RCV_lower, ymax = RCV_upper), 
              fill = "skyblue3", alpha = 0.2) +
  geom_point(aes(color = RCV_flag, shape = RCV_flag), size = 2.5) +
  scale_color_manual(name = "Threshold",
                     values = c("Inside" = "deepskyblue4", "Outside" = "red")) +
  scale_shape_manual(name = "Threshold",
                     values = c("Inside" = 16, "Outside" = 17)) +
  labs(
    title = "Reference Change Value Threshold: LIS Priors",
    x = "Previous Result",
    y = "Current Sample"
  ) +
  theme_minimal()

print(RCV_plot)