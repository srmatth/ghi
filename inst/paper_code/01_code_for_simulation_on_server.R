####### SIMULATION RUN ON SERVER ----


######################################
######### SETUP ######################
######################################

# Load Libraries
library(dplyr)
library(readr)
library(ltrc)
library(survival)

args <- commandArgs(trailingOnly = TRUE)
FILE_NUM <- as.integer(args[1])
cat("FILE NUMBER:", FILE_NUM, "\n")

# Set Parameter Values
alpha <- c(1.7, 0.05, 0.05)
beta <- c(-0.1, -0.5, 2, 0.2)
omega <- c(0, -0.1, 1, -1)
N_SIM <- 5
SAMPLE_SIZE <- 400

# Set prediction grid
covariate_grid <- list(
  t = 7:12,
  x_1 = seq(-1.5, 1.5, by = 0.25),
  x_2 = c(0,1),
  l = 1
) %>%
  expand.grid() %>%
  mutate(
    l_over_t = l / t,
    eps_lower = log(l_over_t / (1 - l_over_t)) - omega[2] * t - omega[3] * x_1 - omega[4] * x_2
  )

# Define simulation run function
ghi_bs_sim_run <- function(i, sample_size, error_fn, error_name, covariate_grid, n_bootstrap = 300,
                           alpha = c(1.7, 0.05, 0.05), beta = c(-0.1, -0.5, 2, 0.2), omega = c(0, -0.1, 1, -1)) {
  cat("Iteration", i, "\n")
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

  # Get actual estimates

  part_1_mod <- glm(delta_s ~ t + x_1 + x_2 + l, data = sim_dat, family = binomial)

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

  part_2_mod_clean <- get_clean_model(part_2_mod)
  beta_hat <- part_1_mod$coefficients[-1]
  omega_hat <- part_2_mod_clean$parameters$beta

  part_2_residuals <- c(part_2_mod_clean$data$residuals)
  part_2_mod_clean$data$fitted_response

  pred_df <- data.frame()
  lps <- covariate_grid$t * omega_hat[1] + covariate_grid$x_1 * omega_hat[2] + covariate_grid$x_2 * omega_hat[3]

  for (j in 1:nrow(covariate_grid)) {

    ## Compute Estimate
    mu_0 <- 1 - predict(part_1_mod, newdata = covariate_grid[j,1:4], type = "response")

    sub_lp <- lps[j]
    cutoff <- covariate_grid[j,6]
    tau <- got_s %>%
      mutate(
        l_resid_scale = mod_l - omega_hat[1] * t - omega_hat[2] * x_1 - omega_hat[3] * x_2
      ) %>%
      pull(l_resid_scale)
    km_df <- data.frame(epsilon = part_2_residuals, tau = tau, delta = 1) %>%
      mutate(
        g_of_e = exp(sub_lp + epsilon) / (1 + exp(sub_lp + epsilon)),
        g_of_t = exp(sub_lp + tau) / (1 + exp(sub_lp + tau))
      ) %>%
      filter(epsilon > cutoff)
    km_fit <- tryCatch({
      survfit(Surv(time = g_of_t, time2 = g_of_e, event = delta) ~ 1, data = km_df)
    }, error = function(e) {
      return(NULL)
    })
    if (is.null(km_fit)) next
    km_times <- c(0, km_fit$time)
    km_surv <- c(1, km_fit$surv)
    # left reimann sum
    km_estimate <- sum(km_surv[-length(km_surv)] * diff(km_times))

    area_expectation <- mu_0 + (1 - mu_0) * km_estimate

    ## Compute Variance
    g_n <- mu_0
    f_n <- km_estimate

    x_mat <- matrix(c(1, as.numeric(covariate_grid[j,1:4])), ncol = 1)
    sigma_11 <- (g_n * (1 - g_n))^2 * t(x_mat) %*% vcov(part_1_mod) %*% x_mat

    ## Save the results
    tmp_preds <- data.frame(
      t = covariate_grid[j,1],
      x_1 = covariate_grid[j,2],
      x_2 = covariate_grid[j,3],
      l = covariate_grid[j,4],
      lower_bound = covariate_grid[j,6],
      mu_0 = mu_0,
      sub_lp = sub_lp,
      expected_y = km_estimate,
      ghi = area_expectation,
      sigma_11 = sigma_11
    )

    pred_df <- rbind(pred_df, tmp_preds)

  }

  sub_res <- list(
    iter_num = i,
    n = sample_size,
    pct_truncation = 1 - nrow(obs_dat) / nrow(true_dat),
    error_name = error_name,
    beta = beta,
    omega = omega,
    covariate_grid = covariate_grid,
    data = sim_dat,
    part_1_model = part_1_mod,
    part_2_model = part_2_mod_clean,
    part_2_residuals = part_2_residuals,
    part_2_modified_residuals = km_df$epsilon,
    prediction_df = pred_df
  )

  ## bootstrap variance estimate
  bs_data <- list()
  k <- 0
  cat("Bootstrap iterations: ")
  while (k < n_bootstrap) {
    cat(paste0(k, ","))
    samp_indx <- sample(1:nrow(got_s), size = nrow(got_s), replace = TRUE)
    got_s_bs <- got_s %>%
      slice(samp_indx)

    part_2_mod <- tryCatch({
      ltrc(survival::Surv(mod_y, delta_s) ~ t + x_1 + x_2, data = got_s_bs, trunc_time = got_s_bs$mod_l, int_knots = 2)
    },
    error = function(e) {
      print(e)
      i <- i - 1
      return(NULL)
    })
    if (is.null(part_2_mod)) next

    k <- k + 1

    part_2_mod_clean <- get_clean_model(part_2_mod)
    omega_hat <- part_2_mod_clean$parameters$beta

    part_2_residuals <- c(part_2_mod_clean$data$residuals)

    pred_df <- data.frame()
    lps <- covariate_grid$t * omega_hat[1] + covariate_grid$x_1 * omega_hat[2] + covariate_grid$x_2 * omega_hat[3]

    for (j in 1:nrow(covariate_grid)) {

      sub_lp <- lps[j]
      cutoff <- covariate_grid[j,6]
      tau <- got_s %>%
        mutate(
          l_resid_scale = mod_l - omega_hat[1] * t - omega_hat[2] * x_1 - omega_hat[3] * x_2
        ) %>%
        pull(l_resid_scale)
      km_df <- data.frame(epsilon = part_2_residuals, tau = tau, delta = 1) %>%
        mutate(
          g_of_e = exp(sub_lp + epsilon) / (1 + exp(sub_lp + epsilon)),
          g_of_t = exp(sub_lp + tau) / (1 + exp(sub_lp + tau))
        ) %>%
        filter(epsilon > cutoff)
      km_fit <- tryCatch({
        survfit(Surv(time = g_of_t, time2 = g_of_e, event = delta) ~ 1, data = km_df)
      }, error = function(e) {
        return(NULL)
      })
      if (is.null(km_fit)) next
      km_times <- c(0, km_fit$time)
      km_surv <- c(1, km_fit$surv)
      # left reimann sum
      km_estimate <- sum(km_surv[-length(km_surv)] * diff(km_times))

      ## Save the results
      tmp_preds <- data.frame(
        t = covariate_grid[j,1],
        x_1 = covariate_grid[j,2],
        x_2 = covariate_grid[j,3],
        l = covariate_grid[j,4],
        lower_bound = covariate_grid[j,6],
        sub_lp = sub_lp,
        expected_y = km_estimate
      )

      pred_df <- rbind(pred_df, tmp_preds)

    }
    bs_data <- append(bs_data, list(list(
      iter_num = k,
      idx = samp_indx,
      prediction_df = pred_df
    )))
  }
  cat("\n")

  ## Return the results
  sub_res <- append(sub_res, list(bs_data = bs_data))

  sub_res

}

big_result <- list()
for (m in 1:N_SIM) {
  test <- ghi_bs_sim_run(i = m, sample_size = SAMPLE_SIZE, error_fn = rnorm, error_name = "Standard Normal", covariate_grid = covariate_grid)
  big_result <- append(big_result, list(test))
}

out_file <- paste0("ghi_sim_res/error_a_", SAMPLE_SIZE, "_", FILE_NUM, ".rds")
write_rds(big_result, out_file)



