################################################################################
#   Script Name: Establish_LIS_priors_main.R
#   Description: Loads and filters LIS data, establishes a reference interval 
#                using each patient's first observation, excludes samples taken 
#                within 24 hours (if activated), and then forms pairs (using 
#                successive observations) for robust regression.
################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                1. PACKAGES & ENVIRONMENT SETUP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
required_packages <- c(
  "dplyr", "hms", "knitr", 
  "refineR", "robustbase", "here", "ggplot2", "patchwork", "ggtext"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# For reproducibility
set.seed(123)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                2. USER INPUTS & SETTINGS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data reading parameters:
date_format <- "%m/%d/%Y"      # Specify the date format in your data (e.g., "05/27/2022")
time_format <- "%H:%M:%S"      # Specify the time format (e.g., "18:09:43")
delimiter   <- ";"             # Specify the delimiter in your CSV file (e.g., ";" or ",")

# Filter activation flags:
activate_age_filter         <- TRUE  # Activate the age filter?
activate_sex_filter         <- TRUE   # Activate the sex filter?
activate_24hr_filter        <- TRUE   # Activate the 24-hour filter for pairing?
activate_ref_interval_filter<- TRUE   # Activate the reference interval filter?

# Basic filter parameters:
filter_age_min <- 18           # Minimum age
filter_sex     <- c("F")       # Included Sex (adjust as needed)

# Reference interval filter (here: central 99%)
ref_percetile_lower <- 0.005
ref_percetile_upper <- 0.995

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                3. LOAD & PREPARE DATA 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- here("PEB_main","LIS_priors.csv")  # Adjust the path if needed

# Read the CSV using the specified delimiter.
if(delimiter == ";"){
  data <- read.csv2(file_path, stringsAsFactors = FALSE)
} else {
  data <- read.csv(file_path, stringsAsFactors = FALSE)
}

# Convert Sample_date and Sample_time, and ensure Result is numeric.
data <- data %>%
  mutate(
    Sample_date      = as.Date(Sample_date, format = date_format),
    Sample_time      = as_hms(Sample_time),  # as_hms() works with standard time format
    Result           = as.numeric(Result)
  ) %>%
  # Create a combined datetime for accurate time differences.
  mutate(
    Sample_datetime = as.POSIXct(paste(Sample_date, Sample_time),
                                 format = "%Y-%m-%d %H:%M:%S")
  )

# Apply basic filters on age and sex if activated.
if(activate_age_filter){
  data <- data %>% filter(Age >= filter_age_min)
}
if(activate_sex_filter){
  data <- data %>% filter(Sex %in% filter_sex)
}

# Verify a few rows
head(data)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                4. ESTABLISH REFERENCE INTERVAL FROM BASELINE OBSERVATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For each patient, select the first (earliest) observation as baseline.
baseline_data <- data %>%
  group_by(ID) %>%
  slice_min(Sample_datetime, with_ties = FALSE) %>%
  ungroup()

# Compute the reference interval using the baseline Results.
refineR_obj <- findRI(baseline_data$Result)
# Optional: plot(refineR_obj)
lambda        <- refineR_obj$Lambda
sigma_pop_LIS <- refineR_obj$Sigma
mu_pop_LIS    <- refineR_obj$Mu

# 99% reference interval for filtering data (using 0.5% and 99.5% percentiles)
ref <- getRI(refineR_obj, RIperc = c(ref_percetile_lower, ref_percetile_upper))
ref_lower <- ref[1, 2]
ref_upper <- ref[2, 2]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                5. FORM PAIRS OF OBSERVATIONS (b	%24 hours apart)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For each patient, order observations and create pairs of successive measurements.
# If the 24-hour filter is activated, only pairs with a gap b	%24 hours are retained.
paired_data <- data %>%
  group_by(ID) %>%
  arrange(Sample_datetime, .by_group = TRUE) %>%
  mutate(
    lag_Result   = lag(Result),
    lag_datetime = lag(Sample_datetime),
    Time_diff_hours = as.numeric(difftime(Sample_datetime, lag_datetime, units = "hours"))
  ) %>%
  { if(activate_24hr_filter) filter(., !is.na(lag_Result), Time_diff_hours >= 24) 
    else filter(., !is.na(lag_Result)) } %>%  # Use all successive pairs if not activated.
  ungroup() %>%
  rename(
    Baseline_Result   = lag_Result,
    New_Result        = Result,
    Baseline_datetime = lag_datetime
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                6. APPLY REFERENCE INTERVAL FILTER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# If the reference interval filter is activated, exclude pairs where either 
# measurement falls outside the 99% reference interval.
if(activate_ref_interval_filter){
  filtered_data <- paired_data %>%
    filter(
      Baseline_Result >= ref_lower, Baseline_Result <= ref_upper,
      New_Result      >= ref_lower, New_Result      <= ref_upper
    )
} else {
  filtered_data <- paired_data
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                7. BOX-COX TRANSFORMATION & ROBUST REGRESSION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Apply Box-Cox transformation to both measurements.
filtered_data <- filtered_data %>%
  mutate(
    bc_Baseline = BoxCox(Baseline_Result, lambda),
    bc_New      = BoxCox(New_Result, lambda)
  )

# Perform robust linear regression: new measurement as a function of baseline.
robust_fit <- lmrob(bc_New ~ bc_Baseline, data = filtered_data)
B1_LIS    <- coef(robust_fit)[2]  # slope

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                8. PLOTTING THE REGRESSION & RESIDUAL DENSITY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Limit dataset to a maximum of 10,000 randomly sampled points.
set.seed(123)
max_points <- 10000
n_points <- min(nrow(filtered_data), max_points)
plot_data <- filtered_data %>% slice_sample(n = n_points)

# Compute predicted values and residuals.
plot_data <- plot_data %>%
  mutate(
    predicted_bc = predict(robust_fit, newdata = .),
    residual_bc  = bc_New - predicted_bc
  )

# Main regression plot.
p_main <- ggplot(plot_data, aes(x = bc_Baseline, y = bc_New)) +
  geom_point(color = "deepskyblue3", alpha = 0.1) +
  geom_line(aes(y = predicted_bc), size = 1, color = "#36454f") +
  labs(
    x = "Baseline (Box-Cox transformed)",
    y = "New (Box-Cox transformed)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "darkgrey", fill = NA, size = 0.5),
    text = element_text(face = "bold")
  )

# Residual density plot.
p_resid <- ggplot(plot_data, aes(x = residual_bc)) +
  geom_density(color = "#36454f", fill = "deepskyblue3", alpha = 0.45) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#36454f", size = 0.5) +
  labs(x = "Residual distribution", y = "Density") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "darkgrey", fill = NA, size = 0.5),
    text = element_text(face = "bold")
  )

# Create a caption string 
caption_text <- sprintf("**B1 = %.3f** at the Box-Cox scale **N; = %.3f**", B1_LIS, lambda)

# Combine the two plots 
combined_plot <- p_main + p_resid + 
  plot_annotation(caption = caption_text) &
  theme(plot.caption = element_markdown(size = 14, hjust = 0.5))

# Print the combined plot
print(combined_plot)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                9. OUTPUT THE LIS PRIORS AS A SUMMARY TABLE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary_table <- data.frame(
  LIS_priors = c("Box-Cox lambda", 
                "Population SD (Box-Cox scale)", 
                "Population Mean (Box-Cox scale)", 
                "Intraclass Correlation - B1 (Box-Cox scale)"),
  Value = c(lambda, sigma_pop_LIS, mu_pop_LIS, B1_LIS)
)

knitr::kable(summary_table, caption = "Summary of LIS Priors")

