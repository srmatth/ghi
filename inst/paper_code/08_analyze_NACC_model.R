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
library(xtable)

A0 <- 65

mod_dat <- read_csv("inst/extdata/NACC_mod_dat.csv")
part_2_mod_clean <- read_rds("inst/extdata/nacc_pt_2_mod.rds")
part_1_mod <- read_rds("inst/extdata/nacc_pt_1_mod.rds")

## Table of Coefficients ----

omegas <- part_2_mod_clean$parameters$beta
sds <- vcov(part_2_mod_clean) %>% diag() %>% sqrt()

betas <- part_1_mod$coefficients[-1]
sds_b <- vcov(part_1_mod) %>% diag() %>% sqrt() %>% `[`(-1)

data.frame(
  estimate_m = round(betas,3),
  ci_m = paste0("(", round(betas - 1.96 * sds_b, 3), ", ", round(betas + 1.96 * sds_b,3), ")"),
  p_val_m = round(2*pnorm(-abs(betas / sds_b)), 4),
  estimate = c(round(omegas,3), "_"),
  ci = c(paste0("(", round(omegas - 1.96 * sds, 3), ", ", round(omegas + 1.96 * sds,3), ")"), "_"),
  p_val = c(round(2*pnorm(-abs(omegas / sds)), 4), "_")
) %>%
  xtable(digits = c(0, 3, 3, 4, 3, 3, 4),row.names = FALSE) %>%
  print(type = "latex")

## Get Predicted Values and with variability
predict_ghi <- function(mod1, mod2, t, educ, is_female, is_married, has_comorbidity, e4_alleles_1, e4_alleles_2, is_race_black, is_race_other, l) {
  lp_mu <- sum(mod1$coefficients[-1] *
                 c(t, educ, is_female, is_married, has_comorbidity, e4_alleles_1, e4_alleles_2, is_race_black, is_race_other, l))
  mu <- exp(lp_mu) / (1 + exp(lp_mu))

  var_mu <- (mu * (1 - mu))^2 * c(1, t, educ, is_female, is_married, has_comorbidity, e4_alleles_1, e4_alleles_2, is_race_black, is_race_other, l) %*% vcov(mod1) %*% c(1, t, educ, is_female, is_married, has_comorbidity, e4_alleles_1, e4_alleles_2, is_race_black, is_race_other, l)


  lp_y <- sum(mod2$parameters$beta *
                c(t, educ, is_female, is_married, has_comorbidity, e4_alleles_1, e4_alleles_2, is_race_black, is_race_other))

  sd_lp_y <- c(t, educ, is_female, is_married, has_comorbidity, e4_alleles_1, e4_alleles_2, is_race_black, is_race_other) %*% vcov(mod2) %*% c(t, educ, is_female, is_married, has_comorbidity, e4_alleles_1, e4_alleles_2, is_race_black, is_race_other)

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
    var_mu = var_mu,
    lp_y = lp_y,
    sd_lp_y = sd_lp_y,
    expected_y = km_estimate,
    ghi_estimate = area_expectation
  )

}

pred_grid <- list(
  t = c(0.5, 1:25),
  educ = c(16),
  is_female = c(0,1),
  is_married = c(0,1),
  has_comorbidity = c(0,1),
  e4_alleles_1 = c(0,1),
  e4_alleles_2 = c(0,1),
  is_race_black = c(0, 1),
  is_race_other = c(0, 1),
  l = 0.25
) %>%
  expand.grid()

res <- purrr::pmap_dfr(
  .l = pred_grid,
  .f = predict_ghi,
  mod1 = part_1_mod,
  mod2 = part_2_mod_clean
)

## Get Bootstrap iterations
bs_res <- data.frame()
for (f in dir_ls("inst/sim_res/NACC_mod_bs_with_race")) {
  tmp <- read_csv(f)
  bs_res <- bs_res %>%
    bind_rows(tmp)
}


bs_sds <- bs_res %>%
  left_join(
    res %>%
      mutate(lower_lp = lp_y - 5*sd_lp_y, upper_lp = lp_y + 5*sd_lp_y) %>%
      select(t, educ, is_female, is_married, has_comorbidity, e4_alleles_1, e4_alleles_2, is_race_black, is_race_other, l, lower_lp, upper_lp)
  ) %>%
  # filter(lp_y > lower_lp, lp_y < upper_lp) %>%
  group_by(t, educ, is_female, is_married, has_comorbidity, e4_alleles_1, e4_alleles_2, is_race_black, is_race_other, l) %>%
  summarize(
    var_y = var(expected_y)
  )

res_w_var <- res %>%
  left_join(bs_sds) %>%
  mutate(
    var_ghi = (1 - expected_y)^2 * var_mu + (1 - mu)^2 * var_y,
    sd_ghi = sqrt(var_ghi)
  )

#### Big Plots ----

res_w_var %>%
  mutate(
    Race = case_when(
      is_race_black == 1 ~ "Race: Black",
      is_race_other == 1 ~ "Race: Other",
      TRUE ~ "Race: White"
    ),
    Race = factor(Race, levels = c("Race: White", "Race: Black", "Race: Other")),
    married_str = ifelse(is_married == 1, "Married", "Unmarried"),
    sex_str = ifelse(is_female == 1, "Female", "Male"),
    comorbidity_str = ifelse(has_comorbidity == 1, "with Comorbidity", "without Comorbidity"),
    Group = paste0(married_str, " ", sex_str, "\n", comorbidity_str)
  ) %>%
  filter(is_female == 0) %>%
  filter(educ == 16) %>%
  filter(!(e4_alleles_2 == 1 & e4_alleles_1 == 1)) %>%
  filter(!(is_race_black == 1 & is_race_other == 1)) %>%
  mutate(
    num_e4 = case_when(
      e4_alleles_2 == 1 ~ 2,
      e4_alleles_1 == 1 ~ 1,
      TRUE ~ 0
    ),
    Comorbidity = ifelse(has_comorbidity == 1, "Yes", "No"),
    Sex = ifelse(is_female == 1, "Female", "Male"),
    lower = ghi_estimate - 2 * sd_ghi,
    upper = ghi_estimate + 2 * sd_ghi
  ) %>%
  ggplot() +
  aes(x = t + 65, y = ghi_estimate, color = as.factor(num_e4), fill = as.factor(num_e4)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_line(aes(lty = as.factor(num_e4)), size = 1) +
  facet_grid(Group ~ Race, labeller = "label_value") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Age at Death", y = "Predicted GHI", color = "Number of e4 Alleles", fill = "Number of e4 Alleles", lty = "Number of e4 Alleles") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  xlim(c(65, 90))

res_w_var %>%
  mutate(
    Race = case_when(
      is_race_black == 1 ~ "Race: Black",
      is_race_other == 1 ~ "Race: Other",
      TRUE ~ "Race: White"
    ),
    Race = factor(Race, levels = c("Race: White", "Race: Black", "Race: Other")),
    married_str = ifelse(is_married == 1, "Married", "Unmarried"),
    sex_str = ifelse(is_female == 1, "Female", "Male"),
    comorbidity_str = ifelse(has_comorbidity == 1, "with Comorbidity", "without Comorbidity"),
    Group = paste0(married_str, " ", sex_str, "\n", comorbidity_str)
  ) %>%
  filter(is_female == 1) %>%
  filter(educ == 16) %>%
  filter(!(e4_alleles_2 == 1 & e4_alleles_1 == 1)) %>%
  filter(!(is_race_black == 1 & is_race_other == 1)) %>%
  mutate(
    num_e4 = case_when(
      e4_alleles_2 == 1 ~ 2,
      e4_alleles_1 == 1 ~ 1,
      TRUE ~ 0
    ),
    Comorbidity = ifelse(has_comorbidity == 1, "Yes", "No"),
    Sex = ifelse(is_female == 1, "Female", "Male"),
    lower = ghi_estimate - 2 * sd_ghi,
    upper = ghi_estimate + 2 * sd_ghi
  ) %>%
  ggplot() +
  aes(x = t + 65, y = ghi_estimate, color = as.factor(num_e4), fill = as.factor(num_e4)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_line(aes(lty = as.factor(num_e4)), size = 1) +
  facet_grid(Group ~ Race, labeller = "label_value") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Age at Death", y = "Predicted GHI", color = "Number of e4 Alleles", fill = "Number of e4 Alleles", lty = "Number of e4 Alleles") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  xlim(c(65, 90))

## Prediction Table ----

# The table was just made by hand, looking at the `res_w_var` data frame

