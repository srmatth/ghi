# More Model Diagnostics ----

library(ggplot2)
library(ResourceSelection)
library(dplyr)
library(ltrc)

## Part 1 Model ----

mod_dat <- read_csv("inst/extdata/NACC_mod_dat.csv")

part_1_mod <- glm(
  got_dementia ~ t + years_education + is_female + is_married + comorbidity + as.factor(num_e4) + is_race_black + is_race_other + l,
  data = mod_dat,
  family = "binomial",
  y = TRUE
)
summary(part_1_mod)


### Surrogate Residuals ----
surrogate_residuals_logistic <- function(model, nsim = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Extract observed outcomes and fitted probabilities
  y <- model$y
  if (is.null(y)) stop("Model must be fit with y = TRUE")

  p <- fitted(model)

  n <- length(y)
  res <- matrix(NA_real_, nrow = n, ncol = nsim)

  for (s in seq_len(nsim)) {
    U <- numeric(n)

    # Bernoulli CDF intervals
    for (i in seq_len(n)) {
      if (y[i] == 0) {
        U[i] <- runif(1, min = 0, max = 1 - p[i])
      } else {
        U[i] <- runif(1, min = 1 - p[i], max = 1)
      }
    }

    res[, s] <- qnorm(U)
  }

  if (nsim == 1) {
    drop(res)
  } else {
    res
  }
}

r <- surrogate_residuals_logistic(part_1_mod, seed = 26339)
R <- surrogate_residuals_logistic(part_1_mod, nsim = 100, seed = 26339)

qq_data <- apply(R, 2, sort)
qq_mean <- rowMeans(qq_data)
r_mean <- rowMeans(R)

env_bottom <- apply(qq_data, 1, min)
env_top <- apply(qq_data, 1, max)

#### QQ plot for surrogate residuals ----

ggplot() +
  aes(sample = r) +
  geom_ribbon(
    aes(x = qnorm(ppoints(nrow(R))), ymin = env_bottom, ymax = env_top),
    fill = "blue",
    alpha = 0.25
  ) +
  geom_qq(alpha = 0.25) +
  geom_qq_line(color = "red", lty = "dashed") +
  xlab("Theoretical Quantiles") +
  ylab("Sample Quantiles")  +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  )

#### Surrogate residuals vs. linear predictor ----
eta <- predict(part_1_mod, type = "link")

ggplot() +
  aes(x = eta, y = r) +
  geom_point(alpha = 0.25) +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  geom_smooth(method = "loess", se = FALSE, alpha = 0.5) +
  xlab("Linear Predictor") +
  ylab("Surrogate Residual")  +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  )

#### Surrogate residuals vs. covariates ----

ggplot() +
  aes(x = mod_dat$years_education, y = r) +
  geom_point(alpha = 0.25) +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  geom_smooth(method = "loess", se = FALSE, alpha = 0.5) +
  xlab("Years of Education") +
  ylab("Surrogate Residual")  +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  )

ggplot() +
  aes(x = mod_dat$t + 65, y = r) +
  geom_point(alpha = 0.25) +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  geom_smooth(method = "loess", se = FALSE, alpha = 0.5) +
  xlab("Age at Death") +
  ylab("Surrogate Residual")  +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  )

ggplot() +
  aes(x = mod_dat$l + 65, y = r) +
  geom_point(alpha = 0.25) +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  geom_smooth(method = "loess", se = FALSE, alpha = 0.5) +
  xlab("Age at Study Entry") +
  ylab("Surrogate Residual")  +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  )

### Hosmer-Lemeshow Test ----

hl <- hoslem.test(
  x = mod_dat$got_dementia,
  y = fitted(part_1_mod),
  g = 10
)
hl

hl_data <- data.frame(
  observed = hl$observed[, 2],
  expected = hl$expected[, 2],
  group = seq_len(nrow(hl$observed))
)

# hl_data %>%
#   ggplot() +
#   aes(x = expected, y = observed) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0, color = "red", lty = "dashed") +
#   xlab("Expected Events") +
#   ylab("Observed Events") +
#   theme_bw() +
#   theme(
#     legend.position = "bottom",
#     text = element_text(family = "Times"),
#     strip.background = element_blank(),
#     strip.text = element_text(size = 15),
#     axis.text = element_text(size = 15),
#     axis.title = element_text(size = 18)
#   )

calib <- mod_dat %>%
  mutate(
    p_hat = fitted(part_1_mod),
    group = ntile(p_hat, 10)
  ) %>%
  group_by(group) %>%
  summarize(
    mean_p = mean(p_hat),
    obs_rate = mean(got_dementia),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    se = sqrt(obs_rate * (1 - obs_rate) / n),
    lower = obs_rate - 1.96 * se,
    upper = obs_rate + 1.96 * se
  )

calib %>%
  ggplot() +
  geom_point(
    aes(x = mean_p, y = obs_rate),
    size = 2
  ) +
  geom_segment(
    aes(x = mean_p, xend = mean_p, y = lower, yend = upper)
  ) +
  geom_segment(
    aes(x = mean_p - 0.005, xend = mean_p + 0.005, y = lower, yend = lower)
  ) +
  geom_segment(
    aes(x = mean_p - 0.005, xend = mean_p + 0.005, y = upper, yend = upper)
  ) +
  geom_abline(slope = 1, intercept = 0, color = "red", lty = "dashed") +
  xlab("Mean Predicted Probability") +
  ylab("Observed Event Rate") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  )

# Loess-smoothing approach for the calibration plot
predicted <- fitted(part_1_mod)
lo <- loess(got_dementia ~ predicted, data = mod_dat)
x <- seq(0, 1, length.out = 200)
y_hat <- predict(lo, newdata = data.frame(predicted = x))

ggplot() +
  aes(x = x, y = y_hat) +
  geom_line(size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "red", lty = "dashed") +
  xlab("Predicted Probability") +
  ylab("(Smoothed) Observed Probability")  +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  ) +
  xlim(c(0,1)) +
  ylim(c(0,1))

## Part 2 Model ----

dem_mod_dat <- mod_dat %>%
  filter(got_dementia == 1) %>%
  mutate(
    num_e4_1 = ifelse(num_e4 == 1, 1, 0),
    num_e4_2 = ifelse(num_e4 == 2, 1, 0)
  )

tictoc::tic()
part_2_mod <- ltrc(
  survival::Surv(mod_y, got_dementia) ~ I(t^3) + I(t^2) + t + years_education + is_female + is_married + comorbidity + is_race_black + is_race_other + num_e4_1 + num_e4_2,
  data = dem_mod_dat, trunc_time = dem_mod_dat$mod_l, n_start = 10, int_knots = 1
)
tictoc::toc()

part_2_mod_clean <- get_clean_model(part_2_mod)
readr::write_rds(part_2_mod_clean, "inst/extdata/nacc_pt2_mod_with_martingale_resids.rds")
part_2_mod_clean <- readr::read_rds("inst/extdata/nacc_pt2_mod_with_martingale_resids.rds")
part_2_mod_clean$parameters$beta

### Martingale Residuals ----

martingale_resids <- part_2_mod_clean$data$martingale_residuals

#### Cumulative Martingale Residuals ----

cum_mart_resid_plot <- function(mart, x, x_label = "Age at Death") {
  n <- length(mart)

  # jitter if the variable is discrete
  if (length(x) != length(unique(x))) {
    x <- x + rnorm(length(x), 0, 0.1)
  }

  ord <- order(x)
  x_ord <- x[ord]
  m_ord <- mart[ord]

  R_obs <- cumsum(m_ord) / sqrt(n)

  B <- 1000
  R_sim <- matrix(NA_real_, nrow = n, ncol = B)

  for (b in seq_len(B)) {
    G <- rnorm(n)
    M_star <- m_ord * G
    R_sim[, b] <- cumsum(M_star) / sqrt(n)
  }

  lower <- apply(R_sim, 1, quantile, probs = 0.025)
  upper <- apply(R_sim, 1, quantile, probs = 0.975)

  if (x_label == "Age at Death") plot_x <- x_ord + 65
  else plot_x <- x_ord

  p <- ggplot() +
    geom_ribbon(
      aes(x = plot_x, ymin = lower, ymax = upper),
      fill = "blue",
      alpha = 0.25
    ) +
    geom_line(
      aes(x = plot_x, y = R_obs)
    ) +
    geom_hline(yintercept = 0, color = "red", lty = "dashed") +
    xlab(x_label) +
    ylab("Cumulative Martingale Residual") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      text = element_text(family = "Times"),
      strip.background = element_blank(),
      strip.text = element_text(size = 15),
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 18)
    )

  return(p)
}

cum_mart_resid_plot(martingale_resids, dem_mod_dat$t, x_label = "Age at Death")
cum_mart_resid_plot(martingale_resids, dem_mod_dat$years_education, x_label = "Years of Education")

# What about the age at study entry? What does this graph tell us?
cum_mart_resid_plot(martingale_resids, dem_mod_dat$l, x_label = "Age at Study Entry")

#### Martingale Residuals vs. Covariates ----

ggplot() +
  aes(x = dem_mod_dat$years_education, y = martingale_resids) +
  geom_point(alpha = 0.25) +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  geom_smooth(method = "loess", se = FALSE, alpha = 0.5) +
  xlab("Years of Education") +
  ylab("Martingale Residual")  +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  )

ggplot() +
  aes(x = dem_mod_dat$t + 65, y = martingale_resids) +
  geom_point(alpha = 0.25) +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  geom_smooth(method = "loess", se = FALSE, alpha = 0.5) +
  xlab("Age at Death") +
  ylab("Martingale Residual")  +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  )

# This chart doesn't really mean anything
ggplot() +
  aes(x = dem_mod_dat$l + 65, y = martingale_resids) +
  geom_point(alpha = 0.25) +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  geom_smooth(method = "loess", se = FALSE, alpha = 0.5) +
  xlab("Age at Study Entry") +
  ylab("Martingale Residual")  +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  )

