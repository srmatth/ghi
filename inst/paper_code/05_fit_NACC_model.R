## Setup ----

library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(stringr)
library(patchwork)
library(fs)
library(data.table)
library(ltrc)

A0 <- 65

DATA_DIR <- "inst/extdata" # Replace this with your data directory

if (fs::file_exists(fs::path(DATA_DIR, "NACC_mod_dat.csv"))) {
  mod_dat <- read_csv(fs::path(DATA_DIR, "NACC_mod_dat.csv"))
} else {
  mod_dat <- read_csv(fs::path(DATA_DIR, "NACC_mod_dat_simulated.csv"))
}

## Do the modeling ----

part_1_mod <- glm(
  got_dementia ~ t + years_education + is_female + is_married + comorbidity + as.factor(num_e4) + is_race_black + is_race_other + l,
  data = mod_dat,
  family = "binomial"
)
summary(part_1_mod)

dem_mod_dat <- mod_dat %>%
  filter(got_dementia == 1) %>%
  mutate(
    num_e4_1 = ifelse(num_e4 == 1, 1, 0),
    num_e4_2 = ifelse(num_e4 == 2, 1, 0)
  )

tictoc::tic()
part_2_mod <- ltrc(
  survival::Surv(mod_y, got_dementia) ~ t + years_education + is_female + is_married + comorbidity + num_e4_1 + num_e4_2 + is_race_black + is_race_other,
  data = dem_mod_dat, trunc_time = dem_mod_dat$mod_l, n_start = 10, int_knots = 1
)
tictoc::toc()

part_2_mod_clean <- get_clean_model(part_2_mod)
write_rds(part_2_mod_clean, fs::path(DATA_DIR, "nacc_pt_2_mod.rds"))
write_rds(part_1_mod, fs::path(DATA_DIR, "nacc_pt_1_mod.rds"))

# part_2_mod_clean$parameters$beta
# vcov(part_2_mod_clean) %>% diag() %>% sqrt()
# part_2_mod_clean$metrics$log_likelihood
