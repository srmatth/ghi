## NACC model bootstrap on the server

## Setup ----
library(ltrc)
library(dplyr)
library(readr)
library(purrr)

args <- commandArgs(trailingOnly = TRUE)
FILE_NUM <- as.integer(args[1])
cat("FILE NUMBER:", FILE_NUM, "\n")

mod_dat <- read_csv("NACC_mod_dat.csv")

predict_ghi <- function(mod1, mod2, t, educ, is_female, is_married, has_comorbidity, e4_alleles_1, e4_alleles_2, is_race_black, is_race_other, l) {
  lp_mu <- sum(mod1$coefficients[-1] *
                 c(t, educ, is_female, is_married, has_comorbidity, e4_alleles_1, e4_alleles_2, is_race_black, is_race_other, l))
  mu <- exp(lp_mu) / (1 + exp(lp_mu))

  lp_y <- sum(mod2$parameters$beta *
                c(t, educ, is_female, is_married, has_comorbidity, e4_alleles_1, e4_alleles_2, is_race_black, is_race_other))
  mod_l <- log((l/t) / (1 - (l/t)))
  lower_bound <- mod_l - lp_y
  tau <- dem_mod_dat %>%
    mutate(
      l_resid_scale = mod_l - lp_y
    ) %>%
    pull(l_resid_scale)
  km_df <- data.frame(epsilon = mod2$data$residuals, tau = tau, delta = 1) %>%
    mutate(
      g_of_e = exp(lp_y + epsilon) / (1 + exp(lp_y + epsilon)),
      g_of_t = exp(lp_y + tau) / (1 + exp(lp_y + tau))
    ) %>%
    filter(epsilon > lower_bound, g_of_e > g_of_t)
  if (nrow(km_df) < 3) return(data.frame())
  km_fit <- tryCatch({
    survfit(Surv(time = g_of_t, time2 = g_of_e, event = delta) ~ 1, data = km_df)
  }, error = function(e) {
    return(NULL)
  })
  if (is.null(km_fit)) return(data.frame())
  km_times <- c(0, km_fit$time)
  km_surv <- c(1, km_fit$surv)
  # left reimann sum
  km_estimate <- sum(km_surv[-length(km_surv)] * diff(km_times))

  area_expectation <- (1 - mu) + (mu) * km_estimate

  data.frame(
    t = t,
    educ = educ,
    is_female = is_female,
    is_married = is_married, has_comorbidity = has_comorbidity,
    e4_alleles_1 = e4_alleles_1,
    e4_alleles_2 = e4_alleles_2,
    is_race_black = is_race_black,
    is_race_other = is_race_other,
    l = l,
    lp_mu = lp_mu,
    mu = mu,
    lp_y = lp_y,
    expected_y = km_estimate,
    ghi_estimate = area_expectation
  )

}

## Part 1 Model ----

# (for prediction)
part_1_mod <- glm(
  got_dementia ~ t + years_education + is_female + is_married + comorbidity + as.factor(num_e4) + is_race_black + is_race_other + l,
  data = mod_dat,
  family = "binomial"
)

dem_mod_dat <- mod_dat %>%
  filter(got_dementia == 1) %>%
  mutate(
    num_e4_1 = ifelse(num_e4 == 1, 1, 0),
    num_e4_2 = ifelse(num_e4 == 2, 1, 0)
  )

pred_grid <- list(
  t = 5:25,
  educ = c(16),
  is_female = c(0,1),
  is_married = c(0,1),
  has_comorbidity = c(0,1),
  e4_alleles_1 = c(0,1),
  e4_alleles_2 = c(0,1),
  is_race_black = c(0,1),
  is_race_other = c(0),
  l = 1
) %>%
  expand.grid()

## Bootstrap iterations -----

all_bs_preds <- data.frame()
for (k in 1:6) {
  cat("Iteration:", k, "\n")
  dem_mod_dat_sub <- dem_mod_dat %>%
    slice_sample(n = nrow(.), replace = TRUE)

  part_2_mod_sub <- ltrc(
    survival::Surv(mod_y, got_dementia) ~ t + years_education + is_female + is_married + comorbidity + num_e4_1 + num_e4_2 + is_race_black + is_race_other,
    data = dem_mod_dat_sub, trunc_time = dem_mod_dat_sub$mod_l, n_start = 10, int_knots = 1
  )
  part_2_mod_sub_clean <- get_clean_model(part_2_mod_sub)

  sub_res <- purrr::pmap_dfr(
    .l = pred_grid,
    .f = predict_ghi,
    mod1 = part_1_mod,
    mod2 = part_2_mod_sub_clean
  )

  all_bs_preds <- all_bs_preds %>%
    rbind(sub_res)
}


out_file <- paste0("NACC_mod_bs/file_", FILE_NUM, ".csv")
write_csv(all_bs_preds, out_file)

