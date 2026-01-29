## Function to run the analysis as if it were any point in the past


library(ltrc)

run_analysis_in_past <- function(date) {
  cat("\n\n===========================================")
  cat("\n RUNNING ANALYSIS FOR:", date)
  cat("\n===========================================\n\n")
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
    filter(visit_date <= ymd(date), death_date <= ymd(date))
  N <- nrow(mod_dat)
  mod_dat <- mod_dat%>%
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
  # print(summary(part_1_mod))

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

  # print(summary(part_2_mod_clean))
  sds <- vcov(part_2_mod_clean) %>% diag() %>% sqrt()

  part_1_lowers <- part_1_mod$coefficients - 1.96 * as.numeric(summary(part_1_mod)$coef[,2])
  part_1_uppers <- part_1_mod$coefficients + 1.96 * as.numeric(summary(part_1_mod)$coef[,2])
  part_2_lowers <- part_2_mod_clean$parameters$beta -1.96*sds
  part_2_uppers <- part_2_mod_clean$parameters$beta +1.96*sds


  # Save the parameters for the e4 alleles
  sub_res <- data.frame(
    cutoff_date = date,
    pop_size = N,
    sample_size = nrow(mod_dat),
    num_dementia = nrow(dem_mod_dat),
    model = c("Part 1", "Part 1", "Part 2", "Part 2"),
    label = c("Num e4:1", "Num e4:2"),
    estimate = c(part_1_mod$coefficients[7:8], part_2_mod_clean$parameter$beta[6:7]),
    ci_lower = c(part_1_lowers[7:8], part_2_lowers[6:7]),
    ci_upper = c(part_1_uppers[7:8], part_2_uppers[6:7])
  )

  return(sub_res)

}

dates_to_analyze <- paste0("20", 10:25, "-01-01")

big_res <- data.frame()
for (d in dates_to_analyze) {
  tryCatch({
    tmp_res <- run_analysis_in_past(d)
    big_res <- rbind(big_res, tmp_res)
  }, error = function(e) {
    print(e)
  })
}


readr::write_csv(big_res, "inst/extdata/model_stationarity.csv")
run_analysis_in_past("2010-01-01")

library(ggplot2)

big_res %>%
  mutate(cutoff_date = lubridate::ymd(cutoff_date)) %>%
  ggplot() +
  geom_ribbon(
    aes(x = cutoff_date, ymin = ci_lower, ymax = ci_upper),
    fill = "gray"
  ) +
  geom_line(
    aes(x = cutoff_date, y = estimate)
  ) +
  facet_grid(label~model, scales = "free") +
  theme_bw()  +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12)
  ) +
  xlab("Year Analysis was Run") +
  ylab("Estimated Effect")

big_res %>%
  distinct(cutoff_date, pop_size, sample_size, num_dementia) %>%
  mutate(pct_of_sample = sample_size / pop_size)



