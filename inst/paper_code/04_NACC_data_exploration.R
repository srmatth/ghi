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

mod_dat <- read_csv("inst/extdata/NACC_mod_dat.csv")

## Data Exploration ----

centers <- mod_dat %>%
  group_by(adc) %>%
  summarize(
    n = n(),
    pct_female = mean(is_female),
    pct_married = mean(is_married),
    study_entry = dataanalysis::continuous_summary(age_at_study_entry),
    death = dataanalysis::continuous_summary(age_at_death),
    dementia = mean(got_dementia),
    pct_black = mean(is_race_black)
  )
centers %>%
  ggplot() +
  aes(x = pct_female, y = pct_married) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

mod_dat %>%
  group_by(is_female, is_married) %>%
  summarize(n = n())


mod_dat %>%
  group_by(is_race_black, is_female) %>%
  summarize(n = n())

mod_dat %>%
  group_by(is_race_black, is_married) %>%
  summarize(n = n())

mod_dat %>%
  group_by(is_female, is_race_black, is_married) %>%
  summarize(n = n())

## Descriptive Statistics ----

mean(mod_dat$age_at_study_entry)
sd(mod_dat$age_at_study_entry)

mean(mod_dat$got_dementia)
sum(mod_dat$got_dementia)
mean(mod_dat$got_dementia==0)
sum(mod_dat$got_dementia==0)

mean(mod_dat$age_at_death)
sd(mod_dat$age_at_death)

sum(mod_dat$is_female)
mean(mod_dat$is_female)
sum(mod_dat$is_female==0)
mean(mod_dat$is_female==0)

mean(mod_dat$years_education)
sd(mod_dat$years_education)

sum(mod_dat$is_married)
mean(mod_dat$is_married)
sum(mod_dat$is_married==0)
mean(mod_dat$is_married==0)

sum(mod_dat$comorbidity)
mean(mod_dat$comorbidity)
sum(mod_dat$comorbidity==0)
mean(mod_dat$comorbidity==0)

sum(mod_dat$is_race_white)
mean(mod_dat$is_race_white)
sum(mod_dat$is_race_black)
mean(mod_dat$is_race_black)
sum(mod_dat$is_race_other)
mean(mod_dat$is_race_other)

sum(mod_dat$num_e4 == 0)
mean(mod_dat$num_e4 == 0)
sum(mod_dat$num_e4 == 1)
mean(mod_dat$num_e4 == 1)
sum(mod_dat$num_e4 == 2)
mean(mod_dat$num_e4 == 2)

mean(mod_dat$age_at_dementia[mod_dat$got_dementia == 1])
sd(mod_dat$age_at_dementia[mod_dat$got_dementia == 1])
mean(mod_dat$y[mod_dat$got_dementia == 1])
sd(mod_dat$y[mod_dat$got_dementia == 1])

## Descriptive Statistics Stratified by APOE status ----
mean(mod_dat$age_at_study_entry[mod_dat$num_e4 == 0])
sd(mod_dat$age_at_study_entry[mod_dat$num_e4 == 0])
mean(mod_dat$age_at_study_entry[mod_dat$num_e4 == 1])
sd(mod_dat$age_at_study_entry[mod_dat$num_e4 == 1])
mean(mod_dat$age_at_study_entry[mod_dat$num_e4 == 2])
sd(mod_dat$age_at_study_entry[mod_dat$num_e4 == 2])

mean(mod_dat$got_dementia[mod_dat$num_e4 == 0])
sum(mod_dat$got_dementia[mod_dat$num_e4 == 0])
mean(mod_dat$got_dementia[mod_dat$num_e4 == 1])
sum(mod_dat$got_dementia[mod_dat$num_e4 == 1])
mean(mod_dat$got_dementia[mod_dat$num_e4 == 2])
sum(mod_dat$got_dementia[mod_dat$num_e4 == 2])

mean(mod_dat$got_dementia[mod_dat$num_e4 == 0]==0)
sum(mod_dat$got_dementia[mod_dat$num_e4 == 0]==0)
mean(mod_dat$got_dementia[mod_dat$num_e4 == 1]==0)
sum(mod_dat$got_dementia[mod_dat$num_e4 == 1]==0)
mean(mod_dat$got_dementia[mod_dat$num_e4 == 2]==0)
sum(mod_dat$got_dementia[mod_dat$num_e4 == 2]==0)

mean(mod_dat$age_at_death[mod_dat$num_e4 == 0])
sd(mod_dat$age_at_death[mod_dat$num_e4 == 0])
mean(mod_dat$age_at_death[mod_dat$num_e4 == 1])
sd(mod_dat$age_at_death[mod_dat$num_e4 == 1])
mean(mod_dat$age_at_death[mod_dat$num_e4 == 2])
sd(mod_dat$age_at_death[mod_dat$num_e4 == 2])

sum(mod_dat$is_female[mod_dat$num_e4 == 0])
mean(mod_dat$is_female[mod_dat$num_e4 == 0])
sum(mod_dat$is_female[mod_dat$num_e4 == 1])
mean(mod_dat$is_female[mod_dat$num_e4 == 1])
sum(mod_dat$is_female[mod_dat$num_e4 == 2])
mean(mod_dat$is_female[mod_dat$num_e4 == 2])

sum(mod_dat$is_female[mod_dat$num_e4 == 0]==0)
mean(mod_dat$is_female[mod_dat$num_e4 == 0]==0)
sum(mod_dat$is_female[mod_dat$num_e4 == 1]==0)
mean(mod_dat$is_female[mod_dat$num_e4 == 1]==0)
sum(mod_dat$is_female[mod_dat$num_e4 == 2]==0)
mean(mod_dat$is_female[mod_dat$num_e4 == 2]==0)

mean(mod_dat$years_education[mod_dat$num_e4 == 0])
sd(mod_dat$years_education[mod_dat$num_e4 == 0])
mean(mod_dat$years_education[mod_dat$num_e4 == 1])
sd(mod_dat$years_education[mod_dat$num_e4 == 1])
mean(mod_dat$years_education[mod_dat$num_e4 == 2])
sd(mod_dat$years_education[mod_dat$num_e4 == 2])

sum(mod_dat$is_married[mod_dat$num_e4 == 0])
mean(mod_dat$is_married[mod_dat$num_e4 == 0])
sum(mod_dat$is_married[mod_dat$num_e4 == 1])
mean(mod_dat$is_married[mod_dat$num_e4 == 1])
sum(mod_dat$is_married[mod_dat$num_e4 == 2])
mean(mod_dat$is_married[mod_dat$num_e4 == 2])

sum(mod_dat$is_married[mod_dat$num_e4 == 0]==0)
mean(mod_dat$is_married[mod_dat$num_e4 == 0]==0)
sum(mod_dat$is_married[mod_dat$num_e4 == 1]==0)
mean(mod_dat$is_married[mod_dat$num_e4 == 1]==0)
sum(mod_dat$is_married[mod_dat$num_e4 == 2]==0)
mean(mod_dat$is_married[mod_dat$num_e4 == 2]==0)

sum(mod_dat$comorbidity[mod_dat$num_e4 == 0])
mean(mod_dat$comorbidity[mod_dat$num_e4 == 0])
sum(mod_dat$comorbidity[mod_dat$num_e4 == 1])
mean(mod_dat$comorbidity[mod_dat$num_e4 == 1])
sum(mod_dat$comorbidity[mod_dat$num_e4 == 2])
mean(mod_dat$comorbidity[mod_dat$num_e4 == 2])

sum(mod_dat$comorbidity[mod_dat$num_e4 == 0]==0)
mean(mod_dat$comorbidity[mod_dat$num_e4 == 0]==0)
sum(mod_dat$comorbidity[mod_dat$num_e4 == 1]==0)
mean(mod_dat$comorbidity[mod_dat$num_e4 == 1]==0)
sum(mod_dat$comorbidity[mod_dat$num_e4 == 2]==0)
mean(mod_dat$comorbidity[mod_dat$num_e4 == 2]==0)

sum(mod_dat$is_race_white[mod_dat$num_e4 == 0])
mean(mod_dat$is_race_white[mod_dat$num_e4 == 0])
sum(mod_dat$is_race_white[mod_dat$num_e4 == 1])
mean(mod_dat$is_race_white[mod_dat$num_e4 == 1])
sum(mod_dat$is_race_white[mod_dat$num_e4 == 2])
mean(mod_dat$is_race_white[mod_dat$num_e4 == 2])

sum(mod_dat$is_race_black[mod_dat$num_e4 == 0])
mean(mod_dat$is_race_black[mod_dat$num_e4 == 0])
sum(mod_dat$is_race_black[mod_dat$num_e4 == 1])
mean(mod_dat$is_race_black[mod_dat$num_e4 == 1])
sum(mod_dat$is_race_black[mod_dat$num_e4 == 2])
mean(mod_dat$is_race_black[mod_dat$num_e4 == 2])

sum(mod_dat$is_race_other[mod_dat$num_e4 == 0])
mean(mod_dat$is_race_other[mod_dat$num_e4 == 0])
sum(mod_dat$is_race_other[mod_dat$num_e4 == 1])
mean(mod_dat$is_race_other[mod_dat$num_e4 == 1])
sum(mod_dat$is_race_other[mod_dat$num_e4 == 2])
mean(mod_dat$is_race_other[mod_dat$num_e4 == 2])

mean(mod_dat$age_at_dementia[mod_dat$got_dementia == 1])
sd(mod_dat$age_at_dementia[mod_dat$got_dementia == 1])
mean(mod_dat$y[mod_dat$got_dementia == 1])
sd(mod_dat$y[mod_dat$got_dementia == 1])
