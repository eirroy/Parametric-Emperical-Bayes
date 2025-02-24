##############################################################################
# PEB_LIS_Monte_Carlo_example.R
# Demonstration of a "LIS-based" Parametric Empirical Bayes approach
# for personal reference intervals using simulated data.
#
# Dependencies: refineR, robustbase, ggplot2, patchwork, dplyr
##############################################################################

library(dplyr)
library(refineR)
library(robustbase)
library(ggplot2)
library(patchwork)

set.seed(123)  # For reproducibility

##############################################################################
# 1) Simulate nested data
##############################################################################

# Population-level mean
pop_mu <- 100

# Between-subject CV and within-subject CV
CV_g <- 0.20  # 20% between-subject variation
CV_i <- 0.05  # 5% within-subject variation

# Derived SD for the population distribution
sd_g <- pop_mu * CV_g

n_patients   <- 5000  # Simulate 5,000 individuals (recommended for fitting refineR)
samples_each <- 5     # Each has 5 serial measurements

# Draw individual "true" homeostatic means from N(pop_mu, sd_g^2)
individual_means <- rnorm(
  n    = n_patients,
  mean = pop_mu,
  sd   = sd_g
)

# For each individual, draw 'samples_each' measurements
sim_data <- data.frame()
for (i in seq_len(n_patients)) {
  true_mean_i    <- individual_means[i]
  measurements_i <- rnorm(
    n    = samples_each,
    mean = true_mean_i,
    sd   = true_mean_i * CV_i
  )
  tmp_df <- data.frame(
    patient_id    = i,
    true_mean     = true_mean_i,
    sample_number = seq_len(samples_each),
    measurement   = measurements_i
  )
  sim_data <- rbind(sim_data, tmp_df)
}

##############################################################################
# 2) Use refineR on each patient's first measurement to identify population
#    reference distribution and estimate Box-Cox parameters.
##############################################################################
df_first <- sim_data %>%
  group_by(patient_id) %>%
  slice_min(sample_number, with_ties = FALSE) %>%
  ungroup()

# Apply refineR to first measurements
refineR_obj <- findRI(df_first$measurement)

# Extract 0.5%–99.5% cutoffs (adjust as you see fit)
ref_range <- getRI(refineR_obj, RIperc = c(0.005, 0.995))
ref_lower <- ref_range[1, 2]
ref_upper <- ref_range[2, 2]

# Extract transform parameters
lambda    <- refineR_obj$Lambda
sigma_pop <- refineR_obj$Sigma   # population SD on the Box-Cox scale
mu_LIS_bc <- refineR_obj$Mu      # population mean on the Box-Cox scale

##############################################################################
# 3) Exclude measurements outside the 99% refineR-based interval
##############################################################################
sim_data_filtered <- sim_data %>%
  filter(measurement >= ref_lower, measurement <= ref_upper)

##############################################################################
# 4) Box-Cox transform the remaining measurements
##############################################################################
sim_data_filtered <- sim_data_filtered %>%
  mutate(measurement_bc = BoxCox(measurement, lambda))

##############################################################################
# 5) Form "previous" vs "new" pairs for robust regression
##############################################################################
paired_data <- sim_data_filtered %>%
  group_by(patient_id) %>%
  arrange(sample_number) %>%
  mutate(bc_prev = lag(measurement_bc, 1)) %>%
  filter(!is.na(bc_prev)) %>%
  ungroup()

##############################################################################
# 6) Robust regression to estimate B1 (intraclass correlation)
##############################################################################
robust_fit <- lmrob(measurement_bc ~ bc_prev, data = paired_data)
B1_LIS     <- coef(robust_fit)[2]

##############################################################################
# 7) Iterative PEB:
#    For each subject and each measurement i, define "nobs = i-1",
#    then:
#       Bn       = (B1 * nobs) / [ B1 * nobs + (1 - B1) ]
#       predicted_bc = mu_LIS_bc + (mean_previous_bc - mu_LIS_bc)*Bn
#       threshold_bc = z * sqrt(1 - B1*Bn)*sigma_pop
##############################################################################
alpha   <- 0.05
z_score <- qnorm(1 - alpha/2)

results_list_LIS <- list()
patient_ids      <- unique(sim_data_filtered$patient_id)

for (p in patient_ids) {
  df_p  <- sim_data_filtered %>%
    filter(patient_id == p) %>%
    arrange(sample_number)
  n_obs <- nrow(df_p)
  
  # Container for each subject
  subject_results <- data.frame(
    patient_id       = integer(),
    observation      = integer(),
    measurement_bc   = numeric(),
    predicted_bc     = numeric(),
    lower_bc         = numeric(),
    upper_bc         = numeric()
  )
  
  for (i in seq_len(n_obs)) {
    nobs        <- i - 1  # how many prior measurements
    measure_i_bc <- df_p$measurement_bc[i]
    
    # Mean of prior data on BC scale
    if (nobs == 0) {
      mean_bc <- NA
    } else {
      mean_bc <- mean(df_p$measurement_bc[1:nobs])
    }
    
    # Bn => shrinks from population to the individual's mean
    Bn_LIS <- if (nobs > 0) {
      (B1_LIS * nobs) / (B1_LIS * nobs + (1 - B1_LIS))
    } else 0
    
    # predicted_bc => the best guess of this individual's set point
    pred_bc <- mu_LIS_bc + if (!is.na(mean_bc)) (mean_bc - mu_LIS_bc)*Bn_LIS else 0
    
    # threshold_bc => half-width of the predictive interval
    measure_thresh_bc <- z_score * sqrt(1 - (B1_LIS * Bn_LIS)) * sigma_pop
    
    subject_results <- rbind(
      subject_results,
      data.frame(
        patient_id     = p,
        observation    = i,
        measurement_bc = measure_i_bc,
        predicted_bc   = pred_bc,
        lower_bc       = pred_bc - measure_thresh_bc,
        upper_bc       = pred_bc + measure_thresh_bc
      )
    )
  }
  results_list_LIS[[p]] <- subject_results
}

final_res_LIS <- do.call(rbind, results_list_LIS)

##############################################################################
# 8) Back-transform + define outliers
##############################################################################
# We'll also define a 95% refineR popRI for comparison
ref95        <- getRI(refineR_obj, RIperc = c(0.025, 0.975))
popRI_lower  <- ref95[1, 2]
popRI_upper  <- ref95[2, 2]

final_res_LIS <- final_res_LIS %>%
  mutate(
    measurement_original = invBoxCox(measurement_bc, lambda),
    predicted_original   = invBoxCox(predicted_bc,   lambda),
    lower_original       = invBoxCox(lower_bc,       lambda),
    upper_original       = invBoxCox(upper_bc,       lambda),
    
    # "dynamic" outlier vs. PEB intervals
    Outlier_dynamic = case_when(
      measurement_original < lower_original | measurement_original > upper_original ~ "Outside",
      TRUE ~ "Inside"
    ),
    # "static" outlier vs. refineR 95% popRI
    Outlier_static = case_when(
      measurement_original < popRI_lower | measurement_original > popRI_upper ~ "Outside",
      TRUE ~ "Inside"
    )
  )

##############################################################################
# 9) Calculate false-positive rates
##############################################################################
num_out_dynamic <- sum(final_res_LIS$Outlier_dynamic == "Outside")
num_in_dynamic  <- sum(final_res_LIS$Outlier_dynamic == "Inside")
FPR_dynamic_LIS <- num_out_dynamic / (num_out_dynamic + num_in_dynamic)

num_out_static  <- sum(final_res_LIS$Outlier_static == "Outside")
num_in_static   <- sum(final_res_LIS$Outlier_static == "Inside")
FPR_static_LIS  <- num_out_static / (num_out_static + num_in_static)

cat("\n[LIS-based PEB] Dynamic-interval FPR:", round(FPR_dynamic_LIS, 4),
    "\n[LIS-based PEB] Static popRI FPR:     ", round(FPR_static_LIS, 4), "\n")

##############################################################################
# 10) Example plot for selected patients
##############################################################################
plot_list <- list()
for (p in 1:1) {
  df_p <- final_res_LIS %>% filter(patient_id == p)
  
  y_min <- min(df_p$lower_original, na.rm = TRUE)
  y_max <- max(df_p$upper_original, na.rm = TRUE)
  
  gg_tmp <- ggplot(df_p, aes(x = observation)) +
    # Dynamic intervals
    geom_ribbon(aes(ymin = lower_original, ymax = upper_original),
                fill = "lightblue", alpha = 0.3) +
    geom_line(aes(y = predicted_original), color = "steelblue3", linewidth = 1) +
    # Points (use size here, as linewidth is not applicable)
    geom_point(aes(y = measurement_original, color = Outlier_dynamic),
               size = 2) +
    scale_color_manual(values = c("Inside" = "black", "Outside" = "red")) +
    # Static popRI lines
    geom_hline(yintercept = popRI_lower, color = "black", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = popRI_upper, color = "black", linetype = "dashed", linewidth = 1) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(
      title = paste("Patient", p),
      x = "Number of Prior Results (n - 1)",
      y = "Concentration"
    ) +
    scale_x_continuous(breaks = seq(1, max(df_p$observation)),
                       labels = function(x) x - 1) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  plot_list[[p]] <- gg_tmp
}
wrap_plots(plot_list, ncol = 1)
