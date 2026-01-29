# Weighted Analyses ----

library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(stringr)
library(patchwork)
library(fs)
library(lme4)
library(truncreg)
library(truncnorm)
library(statmod)
library(TMB)
library(ltrc)
library(sandwich)
library(lmtest)

mod_dat <- read_csv("inst/extdata/NACC_mod_dat.csv")


## Weighted Analysis 1 ----

# Weighting for APOE status

mod_dat <- read_csv("inst/extdata/NACC_mod_dat.csv")

weights_1 <- data.frame(
  num_e4 = c(0, 1, 2),
  obs_props = mod_dat %>% group_by(num_e4) %>% summarize(pct = n() / nrow(.)) %>% pull(pct),
  pop_props = c(0.85, 0.13, 0.02)
) %>%
  mutate(
    weight = pop_props / obs_props
  ) %>%
  select(-obs_props, -pop_props)

mod_dat <- mod_dat %>%
  left_join(weights_1, by = "num_e4")

part_1_mod <- glm(
  got_dementia ~ t + years_education + is_female + is_married + comorbidity + is_race_black + is_race_other + as.factor(num_e4) + l,
  data = mod_dat,
  family = "binomial",
  weights = mod_dat$weight
)
summary(part_1_mod)
# note that standard errors are not correct here, need to use sandwich estimator if you want to do inference

dem_mod_dat <- mod_dat %>%
  filter(got_dementia == 1) %>%
  mutate(
    num_e4_1 = ifelse(num_e4 == 1, 1, 0),
    num_e4_2 = ifelse(num_e4 == 2, 1, 0)
  )

tictoc::tic()
part_2_mod <- ltrc(
  survival::Surv(mod_y, got_dementia) ~ t + years_education + is_female + is_married + comorbidity + is_race_black + is_race_other + num_e4_1 + num_e4_2,
  data = dem_mod_dat, trunc_time = dem_mod_dat$mod_l, n_start = 10, int_knots = 1, weights = dem_mod_dat$weight
)
tictoc::toc()

part_2_mod_clean <- get_clean_model(part_2_mod)
part_2_mod_clean$parameters$beta

## Weighted Analysis 2 ----
# Weighting for Sex and Marital Status


mod_dat <- read_csv("inst/extdata/NACC_mod_dat.csv")

weights_2 <- data.frame(
  is_female = c(0, 0, 1, 1),
  is_married = c(0, 1, 0, 1),
  obs_props = mod_dat %>% group_by(is_female, is_married) %>% summarize(pct = n() / nrow(.)) %>% pull(pct),
  pop_props = c(0.091, 0.364, 0.127, 0.418)
) %>%
  mutate(
    weight = pop_props / obs_props
  ) %>%
  select(-obs_props, -pop_props)

mod_dat <- mod_dat %>%
  left_join(weights_2, by = c("is_female", "is_married"))

part_1_mod <- glm(
  got_dementia ~ t + years_education + is_female + is_married + comorbidity + is_race_black + is_race_other + as.factor(num_e4) + l,
  data = mod_dat,
  family = "binomial",
  weights = mod_dat$weight
)
summary(part_1_mod)
# note that standard errors are not correct here, need to use sandwich estimator if you want to do inference

dem_mod_dat <- mod_dat %>%
  filter(got_dementia == 1) %>%
  mutate(
    num_e4_1 = ifelse(num_e4 == 1, 1, 0),
    num_e4_2 = ifelse(num_e4 == 2, 1, 0)
  )

tictoc::tic()
part_2_mod <- ltrc(
  survival::Surv(mod_y, got_dementia) ~ t + years_education + is_female + is_married + comorbidity + is_race_black + is_race_other + num_e4_1 + num_e4_2,
  data = dem_mod_dat, trunc_time = dem_mod_dat$mod_l, n_start = 10, int_knots = 1, weights = dem_mod_dat$weight
)
tictoc::toc()

part_2_mod_clean <- get_clean_model(part_2_mod)
part_2_mod_clean$parameters$beta


## Weighted Analysis 3 ----
# Weighting for Race, APOE Status, and Gender
# Assumptions: use values in Lumsen et al paper (2015) and then assume those hold across races
# use US census data for race profile of Americans: 74.8% white, 13.7% black, 11.5% other


mod_dat <- read_csv("inst/extdata/NACC_mod_dat.csv")

weights_3 <- list(
  num_e4 = 0:2,
  is_race_black = 0:1,
  is_race_other = 0:1,
  is_female = 0:1
) %>%
  expand.grid() %>%
  filter(!(is_race_black == 1 & is_race_other == 1)) %>%
  arrange(num_e4, is_race_black, is_race_other, is_female) %>%
  mutate(
    obs_props = mod_dat %>% group_by(num_e4, is_race_black, is_race_other, is_female) %>% summarize(n = n(), pct = n() / nrow(.)) %>% pull(pct),
    pop_props = c(
      0.267784, # num_e4: 0 | race: white | sex: male
      0.267784, # num_e4: 0 | race: white | sex: female
      0.04117, # num_e4: 0 | race: other | sex: male
      0.04117, # num_e4: 0 | race: other | sex: female
      0.049046, # num_e4: 0 | race: black | sex: male
      0.049046, # num_e4: 0 | race: black | sex: female
      0.097614, # num_e4: 1 | race: white | sex: male
      0.097614, # num_e4: 1 | race: white | sex: female
      0.0150075, # num_e4: 1 | race: other | sex: male
      0.0150075, # num_e4: 1 | race: other | sex: female
      0.0178785, # num_e4: 1 | race: black | sex: male
      0.0178785, # num_e4: 1 | race: black | sex: female
      0.008602, # num_e4: 2 | race: white | sex: male
      0.008602, # num_e4: 2 | race: white | sex: female
      0.0013225, # num_e4: 2 | race: other | sex: male
      0.0013225, # num_e4: 2 | race: other | sex: female
      0.0015755, # num_e4: 2 | race: black | sex: male
      0.0015755  # num_e4: 2 | race: black | sex: female
    ),
    weight = pop_props / obs_props
  ) %>%
  select(-obs_props, -pop_props)


mod_dat <- mod_dat %>%
  left_join(weights_3, by = c("num_e4", "is_race_black", "is_race_other", "is_female"))

part_1_mod <- glm(
  got_dementia ~ t + years_education + is_female + is_married + comorbidity + is_race_black + is_race_other + as.factor(num_e4) + l,
  data = mod_dat,
  family = "binomial",
  weights = mod_dat$weight
)
summary(part_1_mod)
# note that standard errors are not correct here, need to use sandwich estimator if you want to do inference
vcov_robust <- vcovHC(part_1_mod, type = "HC0")
coeftest(part_1_mod, vcov = vcov_robust)

# CI
beta <- coef(part_1_mod)
se <- sqrt(diag(vcov_robust))
cbind(
  lower = beta - 1.96 * se,
  upper = beta + 1.96 * se
)

dem_mod_dat <- mod_dat %>%
  filter(got_dementia == 1) %>%
  mutate(
    num_e4_1 = ifelse(num_e4 == 1, 1, 0),
    num_e4_2 = ifelse(num_e4 == 2, 1, 0)
  )

tictoc::tic()
part_2_mod <- ltrc(
  survival::Surv(mod_y, got_dementia) ~ t + years_education + is_female + is_married + comorbidity + is_race_black + is_race_other + num_e4_1 + num_e4_2,
  data = dem_mod_dat, trunc_time = dem_mod_dat$mod_l, n_start = 10, int_knots = 1, weights = dem_mod_dat$weight
)
tictoc::toc()

part_2_mod_clean <- get_clean_model(part_2_mod)
part_2_mod_clean$parameters$beta
vcov(part_2_mod_clean) %>% diag() %>% sqrt()

readr::write_rds(part_2_mod_clean, "inst/extdata/weighted_pt_2_mod.rds")


test_mod <- readr::read_rds("inst/extdata/weighted_pt_2_mod.rds")
test_mod$parameters$beta
beta <- test_mod$parameters$beta
se <- sqrt(diag(vcov(test_mod)))
cbind(
  est = beta,
  lower = beta - 1.96 * se,
  upper = beta + 1.96 * se
)

actual_mod <- readr::read_rds("inst/extdata/nacc_pt_2_mod.rds")
actual_mod$parameters$beta
actual_mod$metrics$log_likelihood
