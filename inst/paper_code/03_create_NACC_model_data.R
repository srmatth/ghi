## NACC Data Subsetting ----

## Setup ----

library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(stringr)
library(patchwork)
library(fs)
library(data.table)

A0 <- 65

## Subset the Data ----

nacc <- fread("~/Documents/Data/NACC/2025-01-22_UDS_genotype.csv")

colnames(nacc)

nacc_sub <- nacc %>%
  select(
    NACCID,
    NACCADC,
    VISITMO,
    VISITDAY,
    VISITYR,
    NACCVNUM,
    NACCDAYS,
    NACCFDYS,
    BIRTHMO,
    BIRTHYR,
    NACCDIED,
    NACCMOD,
    NACCYOD,
    SEX,
    HISPANIC,
    RACE,
    EDUC,
    MARISTAT,
    HANDED,
    NACCAGE,
    NACCAGEB,
    SMOKYRS,
    PACKSPER,
    ALCOCCAS,
    ALCFREQ,
    CVHATT,
    HATTYEAR,
    CBSTROKE,
    NACCSTYR,
    TBI,
    TBIYEAR,
    ARTHRIT,
    DIABETES,
    ALCOHOL,
    ABUSOTHR,
    HEIGHT,
    WEIGHT,
    NACCBMI,
    BPSYS,
    BPDIAS,
    HRATE,
    NORMCOG,
    DEMENTED,
    NACCUDSD,
    NACCAPOE,
    NACCNE4S
  )
write_csv(nacc_sub, "~/Documents/Data/NACC/2025-01-22_healthspan_analysis_subset.csv")

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
  filter(NACCID %in% dementia_free_at_baseline, NACCAPOE != 9) # These are the visits that we care about

## Create the data for the model ----

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
    got_dementia = max(DEMENTED),
    age_at_study_entry = min(visit_age),
    age_at_dementia = min(visit_age[DEMENTED == 1]),
    age_at_death = unique(death_age),
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
    got_dementia = ifelse(age_at_dementia >= age_at_death, 0, got_dementia),
    s = ifelse(got_dementia, age_at_dementia - A0, age_at_death - A0),
    t = age_at_death - A0,
    y = ifelse(got_dementia, s/t, 1),
    l = age_at_study_entry - A0,
    scaled_l = l / t,
    mod_y = log(y / (1 - y)),
    mod_l = log(scaled_l / (1 - scaled_l)),
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
    years_education = education
  )

write_csv(mod_dat, "inst/extdata/NACC_mod_dat.csv")

