## Simulate data from fitted NACC models
##
## This script reads the fitted Part 1 (GLM) and Part 2 (LTRC) models
## and generates simulated data with covariate distributions matching
## the original NACC dataset.

library(dplyr)
library(readr)
library(purrr)

# Set seed for reproducibility
set.seed(20240101)

# Read original data for covariate distribution reference
original_data <- read_csv("inst/extdata/NACC_mod_dat.csv")
n_target <- nrow(original_data)  # Target sample size: 5036
n_sim <- ceiling(n_target * 1.5)  # Oversample by 50% to account for truncation filtering

cat("Target sample size:", n_target, "\n")
cat("Oversampling to:", n_sim, "observations\n\n")

# Read fitted models
part_1_mod <- readRDS("inst/extdata/nacc_pt_1_mod.rds")
part_2_mod <- readRDS("inst/extdata/nacc_pt_2_mod.rds")

# Extract Part 2 model residuals (from those with dementia)
residuals_part2 <- part_2_mod$data$residuals

# Get unique ADC codes from original data
adc_codes <- unique(original_data$adc) %>% sort()

# ============================================================================
# SECTION 1: Generate covariate distributions
# ============================================================================

cat("Generating covariates from marginal distributions...\n")

# Continuous covariates ----
# Fit normal distributions to continuous variables
years_educ_dist <- original_data %>%
  filter(!is.na(years_education)) %>%
  summarize(mean = mean(years_education), sd = sd(years_education))

mod_l_dist <- original_data %>%
  filter(!is.na(mod_l), is.finite(mod_l)) %>%
  summarize(mean = mean(mod_l), sd = sd(mod_l))

# Binary/Categorical covariates ----
# Estimate proportions
is_female_prob <- mean(original_data$is_female, na.rm = TRUE)
is_married_prob <- mean(original_data$is_married, na.rm = TRUE)
comorbidity_prob <- mean(original_data$comorbidity, na.rm = TRUE)
is_race_black_prob <- mean(original_data$is_race_black, na.rm = TRUE)
is_race_other_prob <- mean(original_data$is_race_other, na.rm = TRUE)

# num_e4 distribution (0, 1, or 2)
num_e4_dist <- original_data %>%
  filter(!is.na(num_e4)) %>%
  group_by(num_e4) %>%
  summarize(prob = n() / nrow(original_data), .groups = "drop") %>%
  arrange(num_e4)

# t distribution (age at death - 65)
t_dist <- original_data %>%
  filter(!is.na(t)) %>%
  summarize(mean = mean(t), sd = sd(t))

# ============================================================================
# SECTION 2: Sample covariates
# ============================================================================

cat("Sampling", n_sim, "observations using joint demographic distribution...\n")

# Create joint distribution for key demographics to ensure all combinations exist
# This prevents rare combinations (e.g., race=Other + num_e4=2 + male) from being missed
demographic_combos <- expand.grid(
  num_e4 = 0:2,
  is_race_black = 0:1,
  is_race_other = 0:1,
  is_female = 0:1
) %>%
  # Remove impossible combinations (can't be both black and other)
  filter(!(is_race_black == 1 & is_race_other == 1)) %>%
  # Calculate observed proportions from original data
  mutate(
    combo_id = row_number(),
    obs_count = map_int(combo_id, function(i) {
      combo <- slice(., i)
      original_data %>%
        filter(
          num_e4 == combo$num_e4,
          is_race_black == combo$is_race_black,
          is_race_other == combo$is_race_other,
          is_female == combo$is_female
        ) %>%
        nrow()
    })
  ) %>%
  mutate(
    obs_prop = obs_count / sum(obs_count),
    # Ensure minimum representation of 0.5% for rare groups
    adj_prop = pmax(obs_prop, 0.005),
    # Re-normalize to sum to 1
    final_prop = adj_prop / sum(adj_prop)
  )

cat("Found", nrow(demographic_combos), "unique demographic combinations\n")
cat("Rare combinations (< 1%):", sum(demographic_combos$obs_prop < 0.01), "\n\n")

# Sample from joint demographic distribution
sampled_combos <- sample(
  x = demographic_combos$combo_id,
  size = n_sim,
  replace = TRUE,
  prob = demographic_combos$final_prop
)

# Create base data frame with joint demographic structure
sim_data <- tibble(
  combo_id = sampled_combos
) %>%
  left_join(
    demographic_combos %>% select(combo_id, num_e4, is_race_black, is_race_other, is_female),
    by = "combo_id"
  ) %>%
  select(-combo_id) %>%
  mutate(
    # ADC: sample from observed centers
    adc = sample(adc_codes, size = n(), replace = TRUE),

    # Other binary covariates (sampled independently)
    is_married = rbinom(n(), size = 1, prob = is_married_prob),
    comorbidity = rbinom(n(), size = 1, prob = comorbidity_prob),

    # Continuous covariates
    years_education = rnorm(n(), mean = years_educ_dist$mean, sd = years_educ_dist$sd),
    t = rnorm(n(), mean = t_dist$mean, sd = t_dist$sd),
    mod_l = rnorm(n(), mean = mod_l_dist$mean, sd = mod_l_dist$sd),

    # Time period
    death_year = sample(2005:2024, n(), replace = TRUE),

    # Clean up education
    years_education = floor(years_education)
  )

# Ensure reasonable bounds
sim_data <- sim_data %>%
  mutate(
    years_education = pmax(pmin(years_education, 20), 0),  # Bound to [0, 20]
    t = pmax(t, 0),  # t must be positive
    # Back-transform mod_l to l
    l = (exp(mod_l) / (1 + exp(mod_l))) * t
  )

# ============================================================================
# SECTION 3: Generate Part 1 model outcomes (dementia probability)
# ============================================================================

cat("Generating Part 1 outcomes (dementia status)...\n")

# Create dummy variables for factors if needed
sim_data <- sim_data %>%
  mutate(
    num_e4_1 = ifelse(num_e4 == 1, 1, 0),
    num_e4_2 = ifelse(num_e4 == 2, 1, 0),
    is_race_white = 1 - is_race_black - is_race_other
  )

# Use Part 1 model to get linear predictor
# Formula: got_dementia ~ t + years_education + is_female + is_married +
#          comorbidity + as.factor(num_e4) + is_race_black + is_race_other + l

# Extract Part 1 model coefficients
coef_part1 <- coef(part_1_mod)
intercept_part1 <- coef_part1[1]

# Compute linear predictor (excluding intercept)
# Order: t, years_education, is_female, is_married, comorbidity, num_e4_1, num_e4_2,
#        is_race_black, is_race_other, l
sim_data <- sim_data %>%
  mutate(
    lp_part1 = intercept_part1 +
      coef_part1["t"] * t +
      coef_part1["years_education"] * years_education +
      coef_part1["is_female"] * is_female +
      coef_part1["is_married"] * is_married +
      coef_part1["comorbidity"] * comorbidity +
      coef_part1["as.factor(num_e4)1"] * num_e4_1 +
      coef_part1["as.factor(num_e4)2"] * num_e4_2 +
      coef_part1["is_race_black"] * is_race_black +
      coef_part1["is_race_other"] * is_race_other +
      coef_part1["l"] * l,
    mu = exp(lp_part1) / (1 + exp(lp_part1)),
    got_dementia = rbinom(n_sim, size = 1, prob = mu)
  )

cat("Simulated", sum(sim_data$got_dementia), "dementia cases out of", n_sim, "\n")

# ============================================================================
# SECTION 4: Generate Part 2 model outcomes (for dementia cases only)
# ============================================================================

cat("Generating Part 2 outcomes (mod_y for dementia cases)...\n")

# Extract Part 2 model coefficients
coef_part2 <- part_2_mod$parameters$beta

# For those with dementia, generate mod_y from the model + residuals
sim_data <- sim_data %>%
  mutate(
    # Linear predictor from Part 2 model
    # Formula: Surv(mod_y, got_dementia) ~ t + years_education + is_female +
    #          is_married + comorbidity + num_e4_1 + num_e4_2 + is_race_black + is_race_other
    lp_part2 = coef_part2[1] * t +
      coef_part2[2] * years_education +
      coef_part2[3] * is_female +
      coef_part2[4] * is_married +
      coef_part2[5] * comorbidity +
      coef_part2[6] * num_e4_1 +
      coef_part2[7] * num_e4_2 +
      coef_part2[8] * is_race_black +
      coef_part2[9] * is_race_other
  )

# Sample residuals for dementia cases
n_dementia <- sum(sim_data$got_dementia)
sampled_residuals <- sample(residuals_part2, size = n_dementia, replace = TRUE)

# Assign residuals to dementia cases
dementia_idx <- which(sim_data$got_dementia == 1)
sim_data$residual <- NA
sim_data$residual[dementia_idx] <- sampled_residuals

# Generate mod_y as linear predictor + residual (for dementia cases)
# For non-dementia: use the upper bound mod_y = log(1 / (1 - 1)) = Inf, or just 1
sim_data <- sim_data %>%
  mutate(
    mod_y = ifelse(got_dementia == 1,
                   lp_part2 + residual,
                   Inf),  # Non-dementia cases: no event

    # Back-transform mod_y to y
    # y = exp(mod_y) / (1 + exp(mod_y)) for dementia cases
    # y = 1 for non-dementia cases
    y = ifelse(got_dementia == 1 & is.finite(mod_y),
               exp(mod_y) / (1 + exp(mod_y)),
               1)
  )

# ============================================================================
# SECTION 5: Apply truncation constraint and sample down
# ============================================================================

cat("\nApplying left-truncation constraint...\n")

# For dementia cases, we must have mod_y >= mod_l (event time after entry)
# This is the defining feature of left-truncated data
sim_data_valid <- sim_data %>%
  filter(
    # Keep all non-dementia cases
    got_dementia == 0 |
    # For dementia cases, keep only those where mod_y >= mod_l
    (got_dementia == 1 & mod_y >= mod_l)
  )

cat("After truncation filtering:", nrow(sim_data_valid), "valid observations\n")
cat("  Non-dementia cases:", sum(sim_data_valid$got_dementia == 0), "\n")
cat("  Dementia cases (mod_y >= mod_l):", sum(sim_data_valid$got_dementia == 1), "\n")

# Sample down to target size if we have more than needed
if (nrow(sim_data_valid) >= n_target) {
  set.seed(20240102)  # Different seed for sampling
  sim_data_sampled <- sim_data_valid %>%
    slice_sample(n = n_target)
  cat("\nSampled down to target size:", n_target, "\n")
} else {
  warning("Fewer valid observations (", nrow(sim_data_valid),
          ") than target (", n_target, "). Using all valid observations.")
  sim_data_sampled <- sim_data_valid
}

# Use the sampled data for final processing
sim_data <- sim_data_sampled

# ============================================================================
# SECTION 6: Generate timing variables
# ============================================================================

cat("Generating dementia and death timing...\n")

sim_data <- sim_data %>%
  mutate(
    # s: age at dementia onset (for those who got dementia)
    s = ifelse(got_dementia == 1,
               y * t,  # y is the fraction of lifespan from entry to dementia
               NA),

    # age_at_study_entry (back from l = age_at_entry - 65)
    age_at_study_entry = l + 65,

    # age_at_death (from t = age_at_death - 65)
    age_at_death = t + 65,

    # age_at_dementia
    age_at_dementia = ifelse(got_dementia == 1,
                             s + 65,
                             NA),

    # scaled_l
    scaled_l = l / t
  )

# ============================================================================
# SECTION 7: Clean up and output
# ============================================================================

cat("Finalizing simulated dataset...\n")

# Select and order columns to match original dataset structure
sim_data_final <- sim_data  %>%
  mutate(
    id = paste0("SIM_", seq_len(n()))
  ) %>%
  select(
    id,
    adc,
    death_year,
    got_dementia,
    age_at_study_entry,
    age_at_dementia,
    age_at_death,
    s,
    t,
    y,
    mod_y,
    l,
    scaled_l,
    mod_l,
    num_e4,
    is_female,
    is_married,
    comorbidity,
    is_race_white,
    is_race_black,
    is_race_other,
    years_education
  ) %>%
  select(id, adc, everything())

# Write output
output_file <- "inst/extdata/NACC_mod_dat_simulated.csv"
write_csv(sim_data_final, output_file)

cat("\n✓ Simulation complete!\n")
cat("Output file:", output_file, "\n")
cat("Sample size:", nrow(sim_data_final), "\n")
cat("Dementia cases:", sum(sim_data_final$got_dementia), "\n")
cat("Non-dementia cases:", sum(sim_data_final$got_dementia == 0), "\n")
cat("\nTruncation check for dementia cases:\n")
dementia_cases <- sim_data_final %>% filter(got_dementia == 1)
cat("  All mod_y >= mod_l:", all(dementia_cases$mod_y >= dementia_cases$mod_l), "\n")
cat("  Min(mod_y - mod_l):", round(min(dementia_cases$mod_y - dementia_cases$mod_l), 4), "\n")

cat("\nDemographic combination check:\n")
demo_check <- sim_data_final %>%
  group_by(num_e4, is_race_black, is_race_other, is_female) %>%
  summarize(n = n(), .groups = "drop") %>%
  arrange(num_e4, is_race_black, is_race_other, is_female)

cat("  Total unique combinations:", nrow(demo_check), "\n")
cat("  Expected combinations:", nrow(demographic_combos), "\n")
cat("  All combinations represented:", nrow(demo_check) == nrow(demographic_combos), "\n")

if (nrow(demo_check) < nrow(demographic_combos)) {
  missing <- anti_join(
    demographic_combos %>% select(num_e4, is_race_black, is_race_other, is_female),
    demo_check %>% select(num_e4, is_race_black, is_race_other, is_female),
    by = c("num_e4", "is_race_black", "is_race_other", "is_female")
  )
  cat("\n⚠ Missing combinations:\n")
  print(missing)
}

cat("\nDemographic distribution:\n")
print(demo_check)

cat("\nSummary statistics:\n")
print(summary(sim_data_final))
