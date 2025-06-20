## Additional Explorations ----

## 1. Changes in Marital Status ----

## Setup ----

library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(stringr)
library(patchwork)
library(fs)

A0 <- 65

# Subjects who have died and for whom there is APOE info
nacc_sub <- read_csv("~/Documents/Data/NACC/2025-01-22_healthspan_analysis_subset.csv")%>%
  filter(NACCDIED == 1)

over_65 <- nacc_sub %>%
  filter(NACCAGE >= A0)

dementia_free_at_baseline <- over_65 %>%
  group_by(NACCID) %>%
  arrange(NACCVNUM) %>%
  slice_head(n = 1) %>%
  filter(DEMENTED == 0) %>%
  pull(NACCID)

over_65_dementia_free <- over_65 %>%
  filter(NACCID %in% dementia_free_at_baseline, NACCAPOE != 9, !is.na(EDUC)) # These are the visits that we care about

over_65_dementia_free %>%
  ggplot() +
  aes(x = NACCVNUM, y = MARISTAT, group = NACCID) +
  geom_line()

maristat <- over_65_dementia_free %>%
  group_by(NACCID) %>%
  arrange(NACCVNUM) %>%
  summarize(
    baseline_status = first(MARISTAT),
    ending_status = last(MARISTAT)
  ) %>%
  ungroup()

table(maristat$baseline_status, maristat$ending_status)

mean(maristat$baseline_status == maristat$ending_status)

## 2. Survival Models ----

library(survival)
library(readr)
library(dplyr)
library(ggplot2)
library(ltrc)

## Marginal Death Model ----

nacc_sub <- read_csv("~/Documents/Data/NACC/2025-01-22_healthspan_analysis_subset.csv")

over_65 <- nacc_sub %>%
  filter(NACCAGE >= 65)

dementia_free_at_baseline <- over_65 %>%
  group_by(NACCID) %>%
  arrange(NACCVNUM) %>%
  slice_head(n = 1) %>%
  filter(DEMENTED == 0) %>%
  pull(NACCID)

over_65_dementia_free <- over_65 %>%
  filter(NACCID %in% dementia_free_at_baseline, NACCAPOE != 9)

mod_dat <- over_65_dementia_free %>%
  mutate(
    birth_date = str_c(BIRTHYR, BIRTHMO, "01", sep = "-") %>% ymd(),
    visit_date = str_c(VISITYR, VISITMO, VISITDAY, sep = "-") %>% ymd(),
    death_date = str_c(NACCYOD, NACCMOD, "01", sep = "-") %>% ymd(),
    visit_age = as.numeric(difftime(visit_date, birth_date, units = "days")) / 365.25,
    death_age = as.numeric(difftime(death_date, birth_date, units = "days")) / 365.25,
    TBI = ifelse(TBI %in% 1:2, 1, 0),
    CBSTROKE = ifelse(CBSTROKE %in% 1:2, 1, 0),
    CVHATT = ifelse(CVHATT %in% 1:2, 1, 0),
    DIABETES = ifelse(DIABETES %in% 1:2, 1, 0),
    ARTHRIT = ifelse(ARTHRIT %in% 1:2, 1, 0)
  ) %>%
  group_by(NACCID, NACCADC) %>%
  arrange(NACCVNUM) %>%
  summarize(
    age_at_study_entry = min(visit_age),
    age_at_dementia = min(visit_age[DEMENTED == 1]),
    age_at_death = unique(death_age),
    age_at_last_visit = last(visit_age),
    apoe_gene = unique(NACCAPOE),
    num_e4 = unique(NACCNE4S),
    sex = unique(SEX),
    stroke = first(CBSTROKE),
    tbi = first(TBI),
    hatt = first(CVHATT),
    diab = first(DIABETES),
    arthritis = first(ARTHRIT),
    hispanic = first(HISPANIC),
    education = first(EDUC),
    marital_status = first(MARISTAT),
    race = first(RACE)
  ) %>%
  ungroup() %>%
  mutate(
    delta = ifelse(is.na(age_at_death), 0, 1),
    t = ifelse(delta == 1, age_at_death - A0, age_at_last_visit - A0),
    l = age_at_study_entry - A0,
    education = ifelse(education == 99, NA, education),
    is_married = ifelse(marital_status == 1, 1, 0),
    is_female = ifelse(sex == 1, 0, 1),
    comorbidity = ifelse(stroke == 1 | tbi == 1 | hatt == 1 | diab == 1 | arthritis == 1, 1, 0),
    is_race_white = ifelse(race == 1, 1, 0),
    is_race_black = ifelse(race == 2, 1, 0),
    is_race_other = ifelse(!(race %in% 1:2), 1, 0)
  ) %>%
  filter(!is.na(education)) %>%
  select(
    id = NACCID,
    adc = NACCADC,
    delta,
    age_at_study_entry,
    age_at_dementia,
    age_at_death,
    t,
    l,
    num_e4,
    is_female,
    is_married,
    comorbidity,
    is_race_white,
    is_race_black,
    is_race_other,
    years_education = education
  ) %>%
  filter(t > l) %>%
  mutate(
    num_e4_1 = ifelse(num_e4 == 1, 1, 0),
    num_e4_2 = ifelse(num_e4 == 2, 1, 0)
  )


cox_mod <- coxph(
  Surv(l, t, delta) ~ years_education + is_female + is_married + comorbidity + num_e4_1 + num_e4_2 + is_race_black + is_race_other,
  data = mod_dat
)

summary(cox_mod)


## Marginal AD Model ----

mod_dat <- read_csv("inst/extdata/NACC_mod_dat.csv") %>%
  mutate(
    num_e4_1 = ifelse(num_e4 == 1, 1, 0),
    num_e4_2 = ifelse(num_e4 == 2, 1, 0),
    mod_s = log(s),
    mod_l = log(l)
  )

tictoc::tic()
ad_mod <- ltrc(
  Surv(s, got_dementia) ~ years_education + is_female + is_married + comorbidity + is_race_black + is_race_other + num_e4_1 + num_e4_2,
  data = mod_dat,
  trunc_time = mod_dat$l,
  int_knots = 2
)
tictoc::toc()

ad_mod_clean <- get_clean_model(ad_mod)
summary(ad_mod_clean)
vcov(ad_mod_clean)

cox_mod <- coxph(
  Surv(l, s, got_dementia) ~ years_education + is_female + is_married + comorbidity + is_race_black + is_race_other + num_e4_1 + num_e4_2,
  data = mod_dat
)
summary(cox_mod)

# 3. Sensitivity Analysis for value of A* and B*

# This was done on the server by adding values to `l` in the prediction grid.


