## Explorations of NACC for Reviewers ----

## Setup ----

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

A0 <- 65

# get the data
nacc_sub <- read_csv("~/Documents/Data/NACC/2025-01-22_healthspan_analysis_subset.csv")

dataset <- nacc_sub %>%
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
    last_visit_age = max(visit_age),
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
    is_race_other = ifelse(!(race %in% 1:2), 1, 0),
    follow_up = last_visit_age - age_at_study_entry
  ) %>%
  select(
    id = NACCID,
    adc = NACCADC,
    got_dementia,
    age_at_study_entry,
    age_at_dementia,
    age_at_death,
    last_visit_age,
    follow_up,
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


# Subjects who are dementia-free at age A0
over_65 <- nacc_sub %>%
  filter(NACCAGE >= A0)

dementia_free_at_baseline <- over_65 %>%
  group_by(NACCID) %>%
  arrange(NACCVNUM) %>%
  slice_head(n = 1) %>%
  filter(DEMENTED == 0) %>%
  pull(NACCID)

over_65_dementia_free <- over_65 %>%
  filter(NACCID %in% dementia_free_at_baseline) # These are the visits that we care about
length(unique(over_65_dementia_free$NACCID))


#3 Included subjects:
included <- read_csv("inst/extdata/NACC_mod_dat.csv")


## Excluded subjects:
population <- over_65_dementia_free %>%
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
    last_visit_age = max(visit_age),
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
    is_race_other = ifelse(!(race %in% 1:2), 1, 0),
    follow_up = last_visit_age - age_at_study_entry
  ) %>%
  select(
    id = NACCID,
    adc = NACCADC,
    got_dementia,
    age_at_study_entry,
    age_at_dementia,
    age_at_death,
    last_visit_age,
    follow_up,
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

excluded <- population %>%
  filter(!(id %in% included$id))
length(unique(excluded$id))

population %>%
  filter(!is.na(t)) %>%
  nrow()

## Summary Table (Reviewer 2 Question 2) ----

dataanalysis::continuous_summary(included$age_at_study_entry)
dataanalysis::continuous_summary(excluded$age_at_study_entry)

dataanalysis::continuous_summary(included$age_at_death - included$age_at_study_entry)
dataanalysis::continuous_summary(excluded$follow_up)

dataanalysis::categorical_summary(included$got_dementia, level = 0)
dataanalysis::categorical_summary(excluded$got_dementia, level = 0)
dataanalysis::categorical_summary(included$got_dementia, level = 1)
dataanalysis::categorical_summary(excluded$got_dementia, level = 1)


dataanalysis::categorical_summary(included$is_female, level = 0)
dataanalysis::categorical_summary(excluded$is_female, level = 0)
dataanalysis::categorical_summary(included$is_female, level = 1)
dataanalysis::categorical_summary(excluded$is_female, level = 1)

dataanalysis::continuous_summary(included$years_education)
dataanalysis::continuous_summary(excluded$years_education)

dataanalysis::categorical_summary(is.na(included$years_education) %>% as.numeric(), level = 1)
dataanalysis::categorical_summary(is.na(excluded$years_education) %>% as.numeric(), level = 1)


dataanalysis::categorical_summary(included$is_married, level = 0)
dataanalysis::categorical_summary(excluded$is_married, level = 0)
dataanalysis::categorical_summary(included$is_married, level = 1)
dataanalysis::categorical_summary(excluded$is_married, level = 1)

dataanalysis::categorical_summary(included$comorbidity, level = 0)
dataanalysis::categorical_summary(excluded$comorbidity, level = 0)
dataanalysis::categorical_summary(included$comorbidity, level = 1)
dataanalysis::categorical_summary(excluded$comorbidity, level = 1)


dataanalysis::categorical_summary(included$is_race_white, level = 1)
dataanalysis::categorical_summary(excluded$is_race_white, level = 1)
dataanalysis::categorical_summary(included$is_race_black, level = 1)
dataanalysis::categorical_summary(excluded$is_race_black, level = 1)
dataanalysis::categorical_summary(included$is_race_other, level = 1)
dataanalysis::categorical_summary(excluded$is_race_other, level = 1)

dataanalysis::categorical_summary(is.na(excluded$is_race_white), level = TRUE)

dataanalysis::categorical_summary(included$num_e4, level = 0)
dataanalysis::categorical_summary(excluded$num_e4, level = 0)
dataanalysis::categorical_summary(included$num_e4, level = 1)
dataanalysis::categorical_summary(excluded$num_e4, level = 1)
dataanalysis::categorical_summary(included$num_e4, level = 2)
dataanalysis::categorical_summary(excluded$num_e4, level = 2)


dataanalysis::categorical_summary(excluded$num_e4, level = 9)

## Age at enrollment plot ----

custom_sum <- function(x) {
  df <- data.frame(
    xintercept = mean(x)
  )
  return(df)
}

included %>%
  mutate(`e4 Alleles` = num_e4) %>%
  ggplot() +
  aes(x = age_at_study_entry) +
  geom_histogram(bins = 10) +
  facet_wrap(~`e4 Alleles`, labeller = "label_both") +
  theme_bw() +
  theme(

  ) +
  xlab("Age at Study Entry") +
  ylab("Number of Patients") +
  scale_x_continuous(
    breaks = c(65, 75, 85, 95, 105)
  )



## NACC Center Follow-Up (Reviewer 2 Question 3) ----
centers <- dataset %>%
  group_by(adc) %>%
  summarize(
    n = n(),
    pct_female = mean(is_female),
    pct_married = mean(is_married),
    study_entry_age = dataanalysis::continuous_summary(age_at_study_entry),
    follow_up_time = dataanalysis::continuous_summary(follow_up),
    last_visit_age = dataanalysis::continuous_summary(last_visit_age),
    dementia = mean(got_dementia),
    pct_black = mean(is_race_black)
  )
centers %>%
  ggplot() +
  aes(x = pct_female, y = pct_married) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

mean(population$age_at_study_entry)

## NACC Center Descriptive Table ----

library(xtable)
mod_dat %>%
  group_by(adc) %>%
  summarize(
    n = n(),
    study_entry = format(round(mean(age_at_study_entry), 1), nsmall = 1),
    dementia = format(round(mean(got_dementia) * 100, 1), nsmall = 1),
    pct_female = format(round(mean(is_female) * 100, 1), nsmall = 1),
    pct_married = format(round(mean(is_married) * 100, 1), nsmall = 1),
    pct_black = format(round(mean(is_race_black) * 100, 1), nsmall = 1),
    apoe0 = format(round(mean(num_e4 == 0) * 100, 1), nsmall = 1),
    apoe1 = format(round(mean(num_e4 == 1) * 100, 1), nsmall = 1),
    apoe2 = format(round(mean(num_e4 == 2) * 100, 1), nsmall = 1)
  ) %>%
  mutate(
    across(c("dementia", "pct_female", "pct_married", "pct_black", "apoe0", "apoe1", "apoe2"),
           ~paste0(.x, "%"))
  ) %>%
  xtable(
    digits = c(0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1)
  ) %>%
  print(include.rownames = FALSE)

population %>%
  group_by(adc) %>%
  summarize(
    n = n(),
    study_entry = format(round(mean(age_at_study_entry), 1), nsmall = 1),
    dementia = format(round(mean(got_dementia) * 100, 1), nsmall = 1),
    pct_female = format(round(mean(is_female) * 100, 1), nsmall = 1),
    pct_married = format(round(mean(is_married) * 100, 1), nsmall = 1),
    pct_black = format(round(mean(is_race_black) * 100, 1), nsmall = 1),
    apoe0 = format(round(mean(num_e4 == 0) * 100, 1), nsmall = 1),
    apoe1 = format(round(mean(num_e4 == 1) * 100, 1), nsmall = 1),
    apoe2 = format(round(mean(num_e4 == 2) * 100, 1), nsmall = 1)
  ) %>%
  mutate(
    across(c("dementia", "pct_female", "pct_married", "pct_black", "apoe0", "apoe1", "apoe2"),
           ~paste0(.x, "%"))
  ) %>%
  xtable(
    digits = c(0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1)
  ) %>%
  print(include.rownames = FALSE)



## Random effects part 1 model (by center) (Reviewer 2 Question 3) ----

A0 <- 65

mod_dat <- read_csv("inst/extdata/NACC_mod_dat.csv")

## Do the modeling ----

part_1_mod_me <- glmer(
  got_dementia ~ t + years_education + is_female + is_married + comorbidity +
    as.factor(num_e4) + is_race_black + is_race_other + l + (1|adc),
  data = mod_dat,
  family = "binomial"
)
summary(part_1_mod_me)

part_1_mod_me_more <- glmer(
  got_dementia ~ t + years_education + is_female + is_married + comorbidity +
    as.factor(num_e4) + is_race_black + is_race_other + l +
    (is_female + is_married + comorbidity + is_race_black + is_race_other + as.factor(num_e4) + 1 |adc),
  data = mod_dat,
  family = "binomial"
)
summary(part_1_mod_me_more)

## Random effects part 2 model (by center) (Reviewer 2 Question 3) ----

# This uses a normal distribution assumption, otherwise the theory is not there

dem_mod_dat <- mod_dat %>%
  filter(got_dementia == 1) %>%
  mutate(
    num_e4_1 = ifelse(num_e4 == 1, 1, 0),
    num_e4_2 = ifelse(num_e4 == 2, 1, 0)
  )

fit_basic <- truncreg(
  mod_y ~ t + years_education + is_female + is_married + comorbidity + num_e4_1 + num_e4_2 + is_race_black + is_race_other,
  point = dem_mod_dat$mod_l,
  data = dem_mod_dat,
  direction = "left"
)

summary(fit_basic)
resids <- residuals(fit_basic)

# Now do the lmm version (coded from scratch)
# log-likelihood for truncated mixed model (random intercept)
loglik_trunc_lmm <- function(par, y, X, L, id, gh_points = 20) {

  # unpack parameters
  p <- ncol(X)
  beta     <- par[1:p]
  sigma    <- exp(par[p+1])
  sigma_b  <- exp(par[p+2])

  # Gauss–Hermite nodes & weights
  gh <- gauss.quad(gh_points, kind="hermite")
  nodes   <- gh$nodes
  weights <- gh$weights

  logLik <- 0

  for (i in unique(id)) {
    idx <- id == i
    y_i <- y[idx]
    X_i <- X[idx,,drop=FALSE]
    L_i <- L[idx]

    # integrate over random effect b_i via GH quadrature:
    # b_i = sqrt(2)*sigma_b*node
    integrand_vals <- numeric(gh_points)

    for (k in 1:gh_points) {
      b_i <- sqrt(2) * sigma_b * nodes[k]
      mu_i <- as.vector(X_i %*% beta + b_i)

      # density of truncated normal for each measurement
      dens <- dtruncnorm(
        y_i,
        a = L_i,
        b = Inf,
        mean = mu_i,
        sd = sigma
      )

      integrand_vals[k] <- prod(dens) * weights[k]
    }

    # random effects normalizing constant
    Li <- sum(integrand_vals) / sqrt(pi)

    logLik <- logLik + log(Li)
  }

  return(-logLik)  # optim minimizes
}

# Construct X matrix
X <- model.matrix(~ t + years_education + is_female + is_married + comorbidity + num_e4_1 + num_e4_2 + is_race_black + is_race_other,
                  data = dem_mod_dat)

# Starting values
start <- c(
  rep(0, ncol(X)),       # beta
  log(sd(dem_mod_dat$mod_y)),         # log(sigma)
  log(sd(dem_mod_dat$mod_y) / 2)      # log(sigma_b)
)

fit <- optim(
  start,
  loglik_trunc_lmm,
  y  = dem_mod_dat$mod_y,
  X  = X,
  L  = dem_mod_dat$mod_l,
  id = dem_mod_dat$adc,
  method = "BFGS",
  control = list(maxit = 500),
  hessian = TRUE
)

fit$par
vcov_mat <- solve(fit$hessian)

# Standard errors
se <- sqrt(diag(vcov_mat))

vcov_mat
se

fit$par - 1.96 * se
fit$par + 1.96 * se

# Truncated lmm with random slopes now
TMB::compile("inst/paper_code/model/truncated_lmm.cpp")
dyn.load("inst/paper_code/model/truncated_lmm.so")

# Construct X and Z
X <- model.matrix(~ t + years_education + is_female + is_married +
                    comorbidity + num_e4_1 + num_e4_2 +
                    is_race_black + is_race_other, data = dem_mod_dat)

Z <- model.matrix(~ is_female + is_married + comorbidity +
                    num_e4_1 + num_e4_2 + is_race_black + is_race_other,
                  data = dem_mod_dat)

# Add random intercept as first column
Z <- cbind(1, Z)

# Map subjects to 0-based index for TMB
id <- as.integer(as.factor(dem_mod_dat$adc)) - 1

data <- list(
  y = dem_mod_dat$mod_y,
  X = X,
  Z = Z,
  id = id,
  L = dem_mod_dat$mod_l
)

# Starting values
p <- ncol(X)
q <- ncol(Z)
n_subj <- length(unique(id))

parameters <- list(
  beta = rep(0, p),
  log_sigma = 0,
  b_raw = rep(0, n_subj * q),
  log_sd_b = rep(0, q)
)

obj <- MakeADFun(data, parameters, random = "b_raw", DLL = "truncated_lmm")

opt <- nlminb(obj$par, obj$fn, obj$gr)

# standard errors
rep <- sdreport(obj)
summary(rep)


#### Sensitivity Analysis for Excluding Patients ----

# the main idea here is that we can do the analysis as if we were 10 years in the
# past to see if there is any difference in the final model estimates





A0 <- 65

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
  filter(visit_date <= ymd("2015-01-01"), death_date <= ymd("2015-01-01")) %>%
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
library(ltrc)
tictoc::tic()
part_2_mod <- ltrc(
  survival::Surv(mod_y, got_dementia) ~ t + years_education + is_female + is_married + comorbidity + num_e4_1 + num_e4_2 + is_race_black + is_race_other,
  data = dem_mod_dat, trunc_time = dem_mod_dat$mod_l, n_start = 10, int_knots = 1
)
tictoc::toc()

part_2_mod_clean <- get_clean_model(part_2_mod)

summary(part_2_mod_clean)
sds <- vcov(part_2_mod_clean) %>% diag() %>% sqrt()

part_2_mod_clean$parameters$beta -1.96*sds

part_2_mod_clean$parameters$beta +1.96*sds

pnorm(-abs(part_2_mod_clean$parameters$beta / sds))


## Check Truncated Normal Distribution Fit

library(truncnorm)

x <- seq(-3, 3, by = 0.01)
y <- dtruncnorm(x, a = -1.5)
plot(x, y, type = "l")

samp <- rtruncnorm(10000, a = -1.5)
mean(samp)
sd(samp)

part_2_mod <- readRDS("inst/extdata/nacc_pt_2_mod.rds")
mod_dat <- read_csv("inst/extdata/NACC_mod_dat.csv")
dem_mod_dat <- mod_dat %>%
  filter(got_dementia == 1) %>%
  mutate(
    num_e4_1 = ifelse(num_e4 == 1, 1, 0),
    num_e4_2 = ifelse(num_e4 == 2, 1, 0)
  )

resids <- as.numeric(part_2_mod$data$residuals)
trunc_times <- dem_mod_dat$mod_l - part_2_mod$data$predictors %*% part_2_mod$parameters$beta

lklhd <- function(pars, e, l) {
  mu <- pars[1]
  sigma <- pars[2]
  # dens <- dnorm(e, mu, sigma)
  # prob <- pnorm(l, mean = mu, sd = sigma)
  contribs <- log(dnorm(e, mean = mu, sd = sigma)) - log(1 - pnorm(l, mean = mu, sd = sigma))
  cat("Mu:", mu, " | Sigma:", sigma, " | Likelihood:", sum(contribs), "\n")
  return(-sum(contribs))
}

res <- optim(
  par = c(0, 2),
  fn = lklhd,
  e = resids,
  l = trunc_times,
  hessian = TRUE
) # mu = -1.2582, sigma = 1.731
# Oh wait, this won't work because the support changes for each observation

# Naive Approach: Use mean and SD of residuals
est_mu <- -1.258
est_sigma <- 1.731

# QQ-Plot
observed <- (pnorm(resids, est_mu, est_sigma) - pnorm(trunc_times, est_mu, est_sigma)) / (1 - pnorm(trunc_times, est_mu, est_sigma))

qqplot(
  qbeta(ppoints(length(observed)),1,1), observed
)
abline(0, 1, col = "red")

obs_vals <- qnorm(observed, est_mu, est_sigma)
theoretical_vals <- qtruncnorm(qbeta(ppoints(length(observed)),1,1), a = trunc_times, mean = mean(resids), sd = est_sigma)

plot(obs_vals, theoretical_vals)

qqplot(resids, theoretical_vals)

qqplot(
  theoretical_vals, obs_vals
)
abline(0, 1, col = "red")

ggplot() +
  aes(x = sort(resids), y = sort(theoretical_vals)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope=1, intercept=0, color="red") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 18)
  ) +
  xlab("Observed Residual") +
  ylab("Theoretical Value")

data.frame(
  x = seq(-2.5, 7.5, by = 0.01),
  y = dtruncnorm(seq(-2.5, 7.5, by = 0.01), a = trunc_times, mean = est_mu, sd = est_sigma)
)
ggplot() +
  geom_density(aes(x = resids)) +
  geom_line(aes(x = seq(-2.5, 7.5, by = 0.01), y = dtruncnorm(seq(-2.5, 7.5, by = 0.01), a = trunc_times, mean = est_mu, sd = est_sigma)))

## Then let's check it against the nuisance parameter
x <- seq(-6, 6, by = 0.01)
obs_log_haz <- ltrc:::predict_log_hazard(
  x = x,
  gamma = part_2_mod$parameters$gamma,
  knots = part_2_mod$parameters$interior_knots,
  boundary_knots = part_2_mod$parameters$boundary_knots
)

theoretical_log_haz <- log(dnorm(x, est_mu, est_sigma) / pnorm(x, est_mu, est_sigma, lower.tail = FALSE))

ggplot() +
  geom_line(
    aes(x = x, y = obs_log_haz)
  ) +
  geom_line(
    aes(x = x, y = theoretical_log_haz),
    color = "red"
  ) +
  geom_vline(xintercept = -4.65, color = "red", lty = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2.13, color = "red", lty = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -1.44, color = "gray", lty = "dashed") +
  geom_vline(xintercept = 3.12, color = "gray", lty = "dashed") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 18)
  ) +
  xlab("Residual Value") +
  ylab("Log Hazard")


## Model checking ----

mod_dat <- read_csv("inst/extdata/NACC_mod_dat.csv")

part_1_mod <- glm(
  got_dementia ~ t + years_education + is_female + is_married + comorbidity + as.factor(num_e4) + is_race_black + is_race_other + l,
  data = mod_dat,
  family = "binomial"
)
summary(part_1_mod)

library(sure)
sure::resids(part_1_mod)

table(mod_dat$got_dementia)
prop.table(table(mod_dat$got_dementia))

colSums(is.na(mod_dat))

library(car)
vif(part_1_mod)

## Partial residuals

# Extract model matrix and coefficients
X <- model.matrix(part_1_mod)
beta <- coef(part_1_mod)

# Linear predictor
eta_hat <- X %*% beta

# Partial residual for t
partial_resid_t <- eta_hat + residuals(part_1_mod, type = "working") -
  X[, "t"] * beta["t"]

plot_dat <- data.frame(
  t = mod_dat$t,
  partial_resid = partial_resid_t
)

ggplot(plot_dat, aes(x = t + 65, y = partial_resid)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
  labs(
    x = "Age at Death",
    y = "Partial residual (logit scale)"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  )

# Partial residual for l
partial_resid_l <- eta_hat + residuals(part_1_mod, type = "working") -
  X[, "l"] * beta["l"]

plot_dat <- data.frame(
  l = mod_dat$l,
  partial_resid = partial_resid_l
)

ggplot(plot_dat, aes(x = l + 65, y = partial_resid)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
  labs(
    x = "Age at Study Entry",
    y = "Partial residual (logit scale)"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 18)
  )
# Partial residual for education
partial_resid_years_education <- eta_hat + residuals(part_1_mod, type = "working") -
  X[, "years_education"] * beta["years_education"]

plot_dat <- data.frame(
  years_education = mod_dat$years_education,
  partial_resid = partial_resid_years_education
)

ggplot(plot_dat, aes(x = years_education, y = partial_resid)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
  labs(
    x = "Years of Education",
    y = "Partial residual (logit scale)"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 18)
  )

library(visreg)

visreg(part_1_mod, "t", scale = "response")
visreg(part_1_mod, "years_education", scale = "response")
visreg(part_1_mod, "l", scale = "response")


library(splines)

part_1_mod_spline <- glm(
  got_dementia ~ ns(t, 3) + ns(years_education, 3) +
    is_female + is_married + comorbidity +
    as.factor(num_e4) + is_race_black + is_race_other +
    ns(l, 3),
  data = mod_dat,
  family = binomial
)
anova(part_1_mod, part_1_mod_spline, test = "Chisq")

infl <- influence.measures(part_1_mod)
summary(infl)

cd_dat <- data.frame(
  index = seq_len(nrow(mod_dat)),
  cooks = cooks.distance(part_1_mod)
)

ggplot(cd_dat, aes(x = index, y = cooks)) +
  geom_point(alpha = 0.4) +
  geom_hline(
    yintercept = 4 / nrow(mod_dat),
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    y = "Cook's Distance",
    x = "Observation index"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 18)
  )

lev_dat <- data.frame(
  index = seq_len(nrow(mod_dat)),
  leverage = hatvalues(part_1_mod),
  cooks = cooks.distance(part_1_mod)
)

ggplot(lev_dat, aes(x = index, y = leverage)) +
  geom_point(alpha = 0.4) +
  geom_hline(
    yintercept = 2 * mean(lev_dat$leverage),
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    y = "Leverage",
    x = "Observation index"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 18)
  )

ggplot(lev_dat, aes(x = cooks, y = leverage)) +
  geom_point(alpha = 0.4) +
  geom_hline(
    yintercept = 2 * mean(lev_dat$leverage),
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(
    xintercept = 4 / nrow(mod_dat),
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    y = "Leverage",
    x = "Cook's Distance"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 18)
  )

# fit the model without influential points
keep_indx <- which(lev_dat$leverage < 2 * mean(lev_dat$leverage) & lev_dat$cooks < 4/nrow(mod_dat))
part_1_mod_no_influential <- glm(
  got_dementia ~ t + years_education + is_female + is_married + comorbidity + as.factor(num_e4) + is_race_black + is_race_other + l,
  data = mod_dat,
  family = "binomial",
  subset = lev_dat$leverage < 2 * mean(lev_dat$leverage) & lev_dat$cooks < 4/nrow(mod_dat)
)
summary(part_1_mod_no_influential)
part_1_mod_no_influential <- glm(
  got_dementia ~ t + years_education + is_female + is_married + comorbidity + as.factor(num_e4) + is_race_black + is_race_other + l,
  data = mod_dat,
  family = "binomial",
  subset = lev_dat$cooks < 4/nrow(mod_dat)
)
summary(part_1_mod_no_influential)
part_1_mod_no_influential <- glm(
  got_dementia ~ t + years_education + is_female + is_married + comorbidity + as.factor(num_e4) + is_race_black + is_race_other + l,
  data = mod_dat,
  family = "binomial",
  subset = !(lev_dat$leverage > 2 * mean(lev_dat$leverage) & lev_dat$cooks > 4/nrow(mod_dat))
)
summary(part_1_mod_no_influential)



AIC(part_1_mod)
BIC(part_1_mod)

library(pROC)

roc_obj <- roc(mod_dat$got_dementia, fitted(part_1_mod))
plot(roc_obj)
auc(roc_obj)

plot(
  fitted(part_1_mod),
  residuals(part_1_mod, type = "deviance"),
  xlab = "Fitted values",
  ylab = "Deviance residuals"
)
abline(h = 0, lty = 2)


# Use e4 as a continuous measure
part_1_mod_num <- glm(
  got_dementia ~ t + years_education + is_female + is_married +
    comorbidity + num_e4 + is_race_black + is_race_other + l,
  data = mod_dat,
  family = binomial
)
summary(part_1_mod_num)
anova(part_1_mod, part_1_mod_num, test = "Chisq")


# Part 2 model
part_2_mod_clean <- read_rds("inst/extdata/nacc_pt_2_mod.rds")

# Partial residual plots

X <- model.matrix(~ . - 1, data = as.data.frame(part_2_mod_clean$data$predictors))
beta <- part_2_mod_clean$parameters$beta
resid <- as.numeric(part_2_mod_clean$data$residuals)

partial_resid_plot <- function(j, var_name) {

  if (var_name == "Age at Death") add_val <- 65
  else add_val <- 0
  pr <- resid + X[, j] * beta[j]

  plot_dat <- data.frame(
    x = X[, j],
    pr = pr
  )

  ggplot(plot_dat, aes(x = x + add_val, y = pr)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
    labs(
      x = var_name,
      y = "Partial residual"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      text = element_text(family = "Times"),
      strip.background = element_blank(),
      strip.text = element_text(size = 15),
      axis.text = element_text(size = 15),
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title = element_text(size = 18)
    )
}

# Example
partial_resid_plot(j = 1, var_name = "Age at Death")
partial_resid_plot(j = 2, var_name = "Years of Education")

# Leverage
X <- part_2_mod_clean$data$predictors
H <- X %*% solve(t(X) %*% X) %*% t(X)
lev <- diag(H)

lev_dat <- data.frame(
  index = seq_along(lev),
  leverage = lev,
  race = case_when(
    X[,8] == 1 ~ "Black",
    X[,9] == 1 ~ "Other",
    TRUE ~ "White"
  )
)

ggplot(lev_dat, aes(x = index, y = leverage, color = race)) +
  geom_point(alpha = 0.4) +
  geom_hline(
    yintercept = 2 * mean(lev),
    linetype = "dashed",
    color = "red"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 18)
  ) +
  labs(
    color = "Race",
    y = "Leverage",
    x = "Observation index"
  )

high_lev_obs <- lev_dat$index[lev_dat$lev > 0.02]
X[high_lev_obs,] %>% View()
X[high_lev_obs,] %>% colMeans()

# Residuals

res_dat <- data.frame(
  fitted = as.numeric(part_2_mod_clean$data$fitted_response),
  resid = resid,
  martingale_resid = as.numeric(part_2_mod_clean$data$martingale_residuals)
)

ggplot(res_dat, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE, color = "blue") +
  geom_hline(yintercept = mean(res_dat$resid), linetype = "dashed", color="red") +
  labs(
    x = "Fitted values",
    y = "Residuals"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 18)
  )

ggplot(res_dat, aes(x = fitted, y = martingale_resid)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE, color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color="red") +
  labs(
    x = "Fitted value",
    y = "Martingale Residual"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 18)
  )

ggplot(res_dat, aes(x = martingale_resid)) +
  geom_histogram(bins=25) +
  geom_vline(xintercept = mean(res_dat$martingale_resid), linetype = "dashed", color="blue") +
  labs(
    x = "Martingale Residual",
    y = "Count"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 18)
  )

ggplot(
  data.frame(lev = lev, resid = resid, race = lev_dat$race),
  aes(x = lev, y = abs(resid), color = race)
) +
  geom_point(alpha = 0.4) +
  labs(
    color = "Race",
    x = "Leverage",
    y = "Absolute Residual"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 18)
  )


## Fixed ADRC Effects ----


A0 <- 65

mod_dat <- read_csv("inst/extdata/NACC_mod_dat.csv")

## Do the modeling ----

part_1_mod <- glm(
  got_dementia ~ t + years_education + is_female + is_married + comorbidity + as.factor(num_e4) + is_race_black + is_race_other + l + as.factor(adc),
  data = mod_dat,
  family = "binomial"
)
summary(part_1_mod)

part_1_mod <- glm(
  got_dementia ~ t + years_education + is_female + is_married + comorbidity + is_race_black + is_race_other + l + as.factor(adc)*t + as.factor(num_e4) ,
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

library(ltrc)
tictoc::tic()
part_2_mod <- ltrc(
  survival::Surv(mod_y, got_dementia) ~ t + years_education + is_female + is_married + comorbidity + num_e4_1 + num_e4_2 + is_race_black + is_race_other + as.factor(adc),
  data = dem_mod_dat, trunc_time = dem_mod_dat$mod_l, n_start = 10, int_knots = 1
)
tictoc::toc()

part_2_mod_clean <- get_clean_model(part_2_mod)


## Simulation Study for Upsampling based on Covariate ----
# Create the data (from the simulation study)
alpha <- c(1.7, 0.05, 0.05)
beta <- c(-0.1, -0.5, 2, 0.2)
omega <- c(0, -0.1, 1, -1)
N_SIM <- 1000
sample_size <- 800
error_fn <- rnorm

res_diffs <- matrix(NA, nrow = N_SIM, ncol = length(beta) + 1)
sub_ests <- matrix(NA, nrow = N_SIM, ncol = length(beta) + 1)

for (i in 1:N_SIM) {
  if (i %% 10 == 0) print(paste("Working on Interation", i))
  true_dat <- data.frame(
    l = runif(4*sample_size, 0, 5),
    x_1 = runif(4*sample_size, -3, 3),
    x_2 = rbinom(4*sample_size, 1, 0.4)
  ) %>%
    mutate(
      t = exp(alpha[1] + alpha[2] * x_1 + alpha[3] * x_2 + rnorm(4*sample_size, 0, .3)), #runif(3*n, l, 25),
      eps = error_fn(4*sample_size),
      y = exp(omega[1] + omega[2] * t + omega[3] * x_1 + omega[4] * x_2 + eps) / (1 + exp(omega[1] + omega[2] * t + omega[3] * x_1 + omega[4] * x_2 + eps)),
      s = y * t
    )
  obs_dat <- true_dat %>%
    filter(
      s > l, t > l
    ) %>%
    mutate(
      prob_disease = exp(beta[1] * t + beta[2] * x_1 + beta[3] * x_2 + beta[4] * l) / (1 + exp(beta[1] * t + beta[2] * x_1 + beta[3] * x_2 + beta[4] * l)),
      delta_s = rbinom(nrow(.), 1, prob_disease),
      y = ifelse(delta_s, y, 1)
    )
  sim_dat <- obs_dat %>%
    sample_n(size = sample_size)

  # Now downsample when x_2 = 0

  sim_dat_1 <- sim_dat %>%
    filter(x_2 == 1)
  sim_dat_0 <- sim_dat %>%
    filter(x_2 == 0) %>%
    sample_n(size = 100)

  sim_dat_sub <- rbind(sim_dat_1, sim_dat_0)

  # Logistic regression model

  part_1_mod <- glm(delta_s ~ t + x_1 + x_2 + l, data = sim_dat, family = binomial)
  part_1_mod_sub <- glm(delta_s ~ t + x_1 + x_2 + l, data = sim_dat_sub, family = binomial)

  sub_ests[i,] <- part_1_mod_sub$coefficients
  res_diffs[i,] <- part_1_mod_sub$coefficients - part_1_mod$coefficients

}

colMeans(sub_ests)
colMeans(res_diffs)






got_s <- sim_dat %>%
  filter(delta_s == 1) %>%
  mutate(
    mod_y = log(y / (1 - y)),
    scaled_l = l / t,
    mod_l = log(scaled_l / (1 - scaled_l))
  )

part_2_mod <- tryCatch({
  ltrc(survival::Surv(mod_y, delta_s) ~ t + x_1 + x_2, data = got_s, trunc_time = got_s$mod_l, int_knots = 2)
},
error = function(e) {
  print(e)
  return(NULL)
})
if (is.null(part_2_mod)) return(NULL)








