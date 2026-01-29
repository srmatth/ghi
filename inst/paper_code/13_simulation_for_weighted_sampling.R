## Simulation for testing the weighted sampling

## Simulation Study for Upsampling based on Covariate ----

# Load Libraries
library(dplyr)
library(readr)
library(ltrc)
library(survival)

args <- commandArgs(trailingOnly = TRUE)
FILE_NUM <- as.integer(args[1])
cat("FILE NUMBER:", FILE_NUM, "\n")

# Create the data (from the simulation study)
alpha <- c(1.7, 0.05, 0.05)
beta <- c(-0.1, -0.5, 2, 0.2)
omega <- c(0, -0.1, 1, -1)
N_SIM <- 20
sample_size <- 800
error_fn <- rnorm

res <- data.frame()

for (i in 1:N_SIM) {
  print(paste("Working on Interation", i))
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

  got_s <- sim_dat %>%
    filter(delta_s == 1) %>%
    mutate(
      mod_y = log(y / (1 - y)),
      scaled_l = l / t,
      mod_l = log(scaled_l / (1 - scaled_l))
    )
  got_s_sub <- sim_dat_sub %>%
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
  part_2_mod_sub <- tryCatch({
    ltrc(survival::Surv(mod_y, delta_s) ~ t + x_1 + x_2, data = got_s_sub, trunc_time = got_s_sub$mod_l, int_knots = 2)
  },
  error = function(e) {
    print(e)
    return(NULL)
  })
  if (is.null(part_2_mod) | is.null(part_2_mod_sub)) next

  part_2_mod_clean <- get_clean_model(part_2_mod)
  part_2_mod_sub_clean <- get_clean_model(part_2_mod_sub)

  beta_hat <- part_1_mod$coefficients[-1]
  beta_hat_sub <- part_1_mod_sub$coefficients[-1]
  omega_hat <- part_2_mod_clean$parameters$beta
  omega_hat_sub <- part_2_mod_sub_clean$parameters$beta

  sub_res <- data.frame(
    file_num = FILE_NUM,
    iter = i,
    parameters = c(paste0("Beta ", 1:4), paste0("Omega ", 1:3)),
    full_est = c(beta_hat, omega_hat),
    sub_est = c(beta_hat_sub, omega_hat_sub),
    truth = c(beta, omega[2:4])
  )

  res <- rbind(res, sub_res)

}

out_file <- paste0("subsample_sim/file_", FILE_NUM, ".rds")
write_rds(big_result, out_file)

q("no")





##### ANALYZING SIMULATION DONT COPY TO SERVER ----

files <- fs::dir_ls("inst/sim_res/subsample_sim")

all_res <- data.frame()
for (f in files) {
  tmp <- readr::read_csv(f)
  all_res <- rbind(all_res, tmp)
}

all_res <- all_res %>%
  mutate(diff = full_est - sub_est)

all_res %>%
  filter(sub_est < 3, sub_est > -3) %>%
  filter(full_est < 3, full_est > -3) %>%
  ggplot() +
  aes(x = sub_est) +
  geom_histogram() +
  facet_wrap(~parameters, scales = "free_x")

all_res %>%
  filter(sub_est < 3, sub_est > -3) %>%
  filter(full_est < 3, full_est > -3) %>%
  group_by(parameters) %>%
  summarize(
    n = n(),
    avg_full = mean(full_est),
    avg_sub = mean(sub_est),
    truth = mean(truth),
    avg_diff = mean(diff)
  )





