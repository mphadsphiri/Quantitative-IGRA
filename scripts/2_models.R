# Title: Quantitative IGRA model analyses
# Author: Mphatso Phiri
# Date: 2025-05-10 (last edited 2026-04-02)

## Correlation TB1 & TB2 ---------

# Multivariate normal model

# overall
m0_tb1tb2_corr_overall <- brm(
  mvbind(tb1_nil, tb2_nil) ~ 1,
  data = d,
  family = gaussian(),
  prior = c(
    prior(normal(0, 5), class = "Intercept", resp = "tb1nil"),
    prior(normal(0, 5), class = "Intercept", resp = "tb2nil"),
    prior(exponential(1), class = "sigma", resp = "tb1nil"),
    prior(exponential(1), class = "sigma", resp = "tb2nil"),
    prior(lkj_corr_cholesky(2), class = "Lrescor")
  ),
  iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
  control = CONTROL,
  seed = SEED)

# age & sex adjusted
m0_tb1tb2_corr_overall_agesex <- brm(
  mvbind(tb1_nil, tb2_nil) ~ s(age, by=sex),
  data = d,
  family = gaussian(),
  prior = c(
    prior(normal(0, 5), class = "Intercept", resp = "tb1nil"),
    prior(normal(0, 5), class = "Intercept", resp = "tb2nil"),
    prior(exponential(1), class = "sigma", resp = "tb1nil"),
    prior(exponential(1), class = "sigma", resp = "tb2nil"),
    prior(lkj_corr_cholesky(2), class = "Lrescor")
  ),
  iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
  control = CONTROL,
  seed = SEED)


# age 10-19


m0_tb1tb2_corr_1019 <- brm(
  mvbind(tb1_nil, tb2_nil) ~ 1,
  data = d |> filter(age_10y %in% c("10-19")),
  family = gaussian(),
  prior = c(
    prior(normal(0, 5), class = "Intercept", resp = "tb1nil"),
    prior(normal(0, 5), class = "Intercept", resp = "tb2nil"),
    prior(exponential(1), class = "sigma", resp = "tb1nil"),
    prior(exponential(1), class = "sigma", resp = "tb2nil"),
    prior(lkj_corr_cholesky(2), class = "Lrescor")
  ),
  iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
  control = CONTROL,
  seed = SEED)

# age 20-29

m0_tb1tb2_corr_2029 <- brm(
  mvbind(tb1_nil, tb2_nil) ~ 1,
  data = d |> filter(age_10y %in% c("20-29")),
  family = gaussian(),
  prior = c(
    prior(normal(0, 5), class = "Intercept", resp = "tb1nil"),
    prior(normal(0, 5), class = "Intercept", resp = "tb2nil"),
    prior(exponential(1), class = "sigma", resp = "tb1nil"),
    prior(exponential(1), class = "sigma", resp = "tb2nil"),
    prior(lkj_corr_cholesky(2), class = "Lrescor")
  ),
  iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
  control = CONTROL,
  seed = SEED)

# age 30-40 only

m0_tb1tb2_corr_3040 <- brm(
  mvbind(tb1_nil, tb2_nil) ~ 1,
  data = d |> filter(age_10y %in% c("30-40")),
  family = gaussian(),
  prior = c(
    prior(normal(0, 5), class = "Intercept", resp = "tb1nil"),
    prior(normal(0, 5), class = "Intercept", resp = "tb2nil"),
    prior(exponential(1), class = "sigma", resp = "tb1nil"),
    prior(exponential(1), class = "sigma", resp = "tb2nil"),
    prior(lkj_corr_cholesky(2), class = "Lrescor")
  ),
  iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
  control = CONTROL,
  seed = SEED)


# sex

# male 

m0_tb1tb2_corr_sexMale <- brm(
  mvbind(tb1_nil, tb2_nil) ~ 1,
  data = d |> filter(sex %in% c("Male")),
  family = gaussian(),
  prior = c(
    prior(normal(0, 5), class = "Intercept", resp = "tb1nil"),
    prior(normal(0, 5), class = "Intercept", resp = "tb2nil"),
    prior(exponential(1), class = "sigma", resp = "tb1nil"),
    prior(exponential(1), class = "sigma", resp = "tb2nil"),
    prior(lkj_corr_cholesky(2), class = "Lrescor")
  ),
  iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
  control = CONTROL,
  seed = SEED)


# female

m0_tb1tb2_corr_sexFemale <- brm(
  mvbind(tb1_nil, tb2_nil) ~ 1,
  data = d |> filter(sex %in% c("Female")),
  family = gaussian(),
  prior = c(
    prior(normal(0, 5), class = "Intercept", resp = "tb1nil"),
    prior(normal(0, 5), class = "Intercept", resp = "tb2nil"),
    prior(exponential(1), class = "sigma", resp = "tb1nil"),
    prior(exponential(1), class = "sigma", resp = "tb2nil"),
    prior(lkj_corr_cholesky(2), class = "Lrescor")
  ),
  iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
  control = CONTROL,
  seed = SEED)


# TB1 & TB2 responses among QFT-Plus positives ---------------------------------


#### TB1 -------

# age & sex spline interaction

m0_tb1_nil_qia_hurdle_age_sex_sp <- brm(bf(tb1_nil_qia | trunc(ub = 10) ~ s(age, by=sex, k=3),
                                           hu ~ s(age, by=sex, k=3)),
                                        data = d,
                                        family = hurdle_lognormal(),
                                        prior = c(
                                          prior(normal(logit(0.8),2), class = "Intercept", dpar = "hu"),
                                          prior(normal(0,1), class = "b", dpar = "hu"),
                                          prior(student_t(2,0,2), class = "sds", coef = "s(age, by = sex, k = 3)", dpar = "hu"),
                                          # prior(normal(0,5), class = "Intercept"),
                                          prior(normal(0,1), class = "b"),
                                          prior(exponential(1), class="sigma")
                                        ),
                                        iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
                                        control = CONTROL, seed = SEED)


# HIV

m0_tb1_nil_qia_hurdle_hiv <- brm(bf(tb1_nil_qia | trunc(ub = 10) ~ hiv_status,
                                    hu ~ s(age, by=sex, k=3)),
                                 data = d,
                                 family = hurdle_lognormal(),
                                 prior = c(
                                   prior(normal(logit(0.8),1), class = "Intercept", dpar = "hu"),
                                   prior(normal(0,1), class = "b", dpar = "hu"),
                                   prior(student_t(2,0,2), class = "sds", coef = "s(age, by = sex, k = 3)", dpar = "hu"),
                                   prior(normal(0,5), class = "Intercept"),
                                   prior(normal(0,1), class = "b"),
                                   prior(exponential(1), class="sigma")
                                 ),
                                 iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
                                 control = CONTROL, seed = SEED)



# Previous TB treatment

m0_tb1_nil_qia_hurdle_tbrx <- brm(bf(tb1_nil_qia | trunc(ub = 10) ~ tbrxhx2,
                                     hu ~ s(age, by=sex, k=3)),
                                  data = d,
                                  family = hurdle_lognormal(),
                                  prior = c(
                                    prior(normal(logit(0.8),1), class = "Intercept", dpar = "hu"),
                                    prior(normal(0,1), class = "b", dpar = "hu"),
                                    prior(student_t(2,0,2), class = "sds", coef = "s(age, by = sex, k = 3)", dpar = "hu"),
                                    prior(normal(0,5), class = "Intercept"),
                                    prior(normal(0,1), class = "b"),
                                    prior(exponential(1), class="sigma")
                                  ),
                                  iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
                                  control = CONTROL, seed = SEED)

# Household TB contact

m0_tb1_nil_qia_hurdle_hhtb <- brm(bf(tb1_nil_qia | trunc(ub = 10) ~ hh_tbrx2,
                                     hu ~ s(age, by=sex, k=3)),
                                  data = d,
                                  family = hurdle_lognormal(),
                                  prior = c(
                                    prior(normal(logit(0.8),1), class = "Intercept", dpar = "hu"),
                                    prior(normal(0,1), class = "b", dpar = "hu"),
                                    prior(student_t(2,0,2), class = "sds", coef = "s(age, by = sex, k = 3)", dpar = "hu"),
                                    prior(normal(0,5), class = "Intercept"),
                                    prior(normal(0,1), class = "b"),
                                    prior(exponential(1), class="sigma")
                                  ),
                                  iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
                                  control = CONTROL, seed = SEED)


# full model: age, sex & HIV, previous TB treatment

m0_tb1_nil_qia_hurdle_full <- brm(bf(tb1_nil_qia | trunc(ub = 10) ~ s(age, by=sex, k=3) + hiv_status + tbrxhx2 + hh_tbrx2,
                                     hu ~ s(age, by=sex, k=3)),
                                  data = d,
                                  family = hurdle_lognormal(),
                                  prior = c(
                                    prior(normal(logit(0.8),1), class = "Intercept", dpar = "hu"),
                                    prior(normal(0,1), class = "b", dpar = "hu"),
                                    prior(student_t(2,0,2), class = "sds", coef = "s(age, by = sex, k = 3)", dpar = "hu"),
                                    prior(normal(0,5), class = "Intercept"),
                                    prior(normal(0,1), class = "b"),
                                    prior(exponential(1), class="sigma")
                                  ),
                                  iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
                                  control = CONTROL, seed = SEED)

#### TB2 -------------

# age & sex spline interaction

m0_tb2_nil_qia_hurdle_age_sex_sp <- brm(bf(tb2_nil_qia | trunc(ub = 10) ~ s(age, by=sex, k=3),
                                           hu ~ s(age, by=sex, k=3)),
                                        data = d,
                                        family = hurdle_lognormal(),
                                        prior = c(
                                          prior(normal(logit(0.8),1), class = "Intercept", dpar = "hu"),
                                          prior(normal(0,1), class = "b", dpar = "hu"),
                                          prior(student_t(2,0,2), class = "sds", coef = "s(age, by = sex, k = 3)", dpar = "hu"),
                                          prior(normal(0,5), class = "Intercept"),
                                          prior(normal(0,1), class = "b"),
                                          prior(exponential(1), class="sigma")
                                        ),
                                        iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
                                        control = CONTROL, seed = SEED)

# HIV

m0_tb2_nil_qia_hurdle_hiv <- brm(bf(tb2_nil_qia | trunc(ub = 10) ~ hiv_status,
                                    hu ~ s(age, by=sex, k=3)),
                                 data = d,
                                 family = hurdle_lognormal(),
                                 prior = c(
                                   prior(normal(logit(0.8),1), class = "Intercept", dpar = "hu"),
                                   prior(normal(0,1), class = "b", dpar = "hu"),
                                   prior(student_t(2,0,2), class = "sds", coef = "s(age, by = sex, k = 3)", dpar = "hu"),
                                   prior(normal(0,5), class = "Intercept"),
                                   prior(normal(0,1), class = "b"),
                                   prior(exponential(1), class="sigma")
                                 ),
                                 iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
                                 control = CONTROL, seed = SEED)

# Previous TB treatment

m0_tb2_nil_qia_hurdle_tbrx <- brm(bf(tb2_nil_qia | trunc(ub = 10) ~ tbrxhx2,
                                     hu ~ s(age, by=sex, k=3)),
                                  data = d,
                                  family = hurdle_lognormal(),
                                  prior = c(
                                    prior(normal(logit(0.8),1), class = "Intercept", dpar = "hu"),
                                    prior(normal(0,1), class = "b", dpar = "hu"),
                                    prior(student_t(2,0,2), class = "sds", coef = "s(age, by = sex, k = 3)", dpar = "hu"),
                                    prior(normal(0,5), class = "Intercept"),
                                    prior(normal(0,1), class = "b"),
                                    prior(exponential(1), class="sigma")
                                  ),
                                  iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
                                  control = CONTROL, seed = SEED)

# Household TB contact

m0_tb2_nil_qia_hurdle_hhtb <- brm(bf(tb2_nil_qia | trunc(ub = 10) ~ hh_tbrx2,
                                     hu ~ s(age, by=sex, k=3)),
                                  data = d,
                                  family = hurdle_lognormal(),
                                  prior = c(
                                    prior(normal(logit(0.8),1), class = "Intercept", dpar = "hu"),
                                    prior(normal(0,1), class = "b", dpar = "hu"),
                                    prior(student_t(2,0,2), class = "sds", coef = "s(age, by = sex, k = 3)", dpar = "hu"),
                                    prior(normal(0,5), class = "Intercept"),
                                    prior(normal(0,1), class = "b"),
                                    prior(exponential(1), class="sigma")
                                  ),
                                  iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
                                  control = CONTROL, seed = SEED)


# full model: age, sex & HIV, previous TB treatment

m0_tb2_nil_qia_hurdle_full <- brm(bf(tb2_nil_qia | trunc(ub = 10) ~ s(age, by=sex, k=3) + hiv_status + tbrxhx2 + hh_tbrx2,
                                     hu ~ s(age, by=sex, k=3)),
                                  data = d,
                                  family = hurdle_lognormal(),
                                  prior = c(
                                    prior(normal(logit(0.8),1), class = "Intercept", dpar = "hu"),
                                    prior(normal(0,1), class = "b", dpar = "hu"),
                                    prior(student_t(2,0,2), class = "sds", coef = "s(age, by = sex, k = 3)", dpar = "hu"),
                                    prior(normal(0,5), class = "Intercept"),
                                    prior(normal(0,1), class = "b"),
                                    prior(exponential(1), class="sigma")
                                  ),
                                  iter = ITER, warmup = WARMUP, chains = CHAINS, cores = CORES,
                                  control = CONTROL, seed = SEED)







# TB2-TB1 differential (continuous) ---------------------------------------------------


# age & sex

m0_tb2tb1_diff_agesex_sp <- brm(tb2tb1_diff ~ s(age, by=sex, k=3),
                                data = d,
                                prior = c(
                                  prior(normal(0,1.5), class="Intercept"),
                                  prior(student_t(2,0,2), class = "sds"),
                                  prior(normal(0,1), class = "b")
                                ),
                                iter = ITER, chains = CHAINS, cores = CORES,
                                control = CONTROL, seed = SEED)



# HIV

m0_tb2tb1_diff_hiv <- brm(tb2tb1_diff ~ hiv_status,
                          data = d,
                          prior = c(
                            prior(normal(0,1.5), class="Intercept"),
                            prior(normal(0,1), class="b")
                          ),
                          iter = ITER, chains = CHAINS, cores = CORES,
                          control = CONTROL, seed = SEED)

# TB treatment

m0_tb2tb1_diff_tbrx <- brm(tb2tb1_diff ~ tbrxhx2,
                           data = d,
                           prior = c(
                             prior(normal(0,1.5), class="Intercept"),
                             prior(normal(0,1), class="b")
                           ),
                           iter = ITER, chains = CHAINS, cores = CORES,
                           control = CONTROL, seed = SEED)


# Household TB contact

m0_tb2tb1_diff_hhtb <- brm(tb2tb1_diff ~ hh_tbrx2,
                           data = d,
                           prior = c(
                             prior(normal(0,1.5), class="Intercept"),
                             prior(normal(0,1), class="b")
                           ),
                           iter = ITER, chains = CHAINS, cores = CORES,
                           control = CONTROL, seed = SEED)


# Full
m0_tb2tb1_diff_full <- brm(tb2tb1_diff ~ s(age, by=sex, k=3) + hiv_status + tbrxhx2 + hh_tbrx2,
                           data = d,
                           prior = c(
                             prior(normal(0,1.5), class="Intercept"),
                             prior(normal(0,1), class="b")
                           ),
                           iter = ITER, chains = CHAINS, cores = CORES,
                           control = CONTROL, seed = SEED)



# TB2-TB1 > 0.6 IU/ml (categorical) --------------------------------------------


# age & sex

m0_tb2tb1_06_agesex_sp <- brm(tb2tb1_06 ~ s(age, by=sex, k=3),
                              data = d,
                              family = bernoulli(link="logit"),
                              prior = c(
                                prior(normal(logit(0.2),1.5), class="Intercept"),
                                prior(student_t(2,0,2), class = "sds"),
                                prior(normal(0,1), class = "b")
                              ),
                              iter = ITER, chains = CHAINS, cores = CORES,
                              control = CONTROL, seed = SEED)



# HIV

m0_tb2tb1_06_hiv <- brm(tb2tb1_06 ~ hiv_status,
                        data = d,
                        family = bernoulli(link="logit"),
                        prior = c(
                          prior(normal(logit(0.2),1.5), class="Intercept"),
                          prior(normal(0,1), class="b")
                        ),
                        iter = ITER, chains = CHAINS, cores = CORES,
                        control = CONTROL, seed = SEED)



# TB treatment

m0_tb2tb1_06_tbrx <- brm(tb2tb1_06 ~ tbrxhx2,
                         data = d,
                         family = bernoulli(link="logit"),
                         prior = c(
                           prior(normal(logit(0.2),1.5), class="Intercept"),
                           prior(normal(0,1), class="b")
                         ),
                         iter = ITER, chains = CHAINS, cores = CORES,
                         control = CONTROL, seed = SEED)


# Household TB contact

m0_tb2tb1_06_hhtb <- brm(tb2tb1_06 ~ hh_tbrx2,
                         data = d,
                         family = bernoulli(link="logit"),
                         prior = c(
                           prior(normal(logit(0.2),1.5), class="Intercept"),
                           prior(normal(0,1), class="b")
                         ),
                         iter = ITER, chains = CHAINS, cores = CORES,
                         control = CONTROL, seed = SEED)


# Full
m0_tb2tb1_06_full <- brm(tb2tb1_06 ~ s(age, by=sex, k=3) + hiv_status + tbrxhx2 + hh_tbrx2,
                         data = d,
                         family = bernoulli(link="logit"),
                         prior = c(
                           prior(normal(logit(0.2),1.5), class="Intercept"),
                           prior(normal(0,1), class="b")
                         ),
                         iter = ITER, chains = CHAINS, cores = CORES,
                         control = CONTROL, seed = SEED)

# QFT-Plus thresholds ---------------------------------

# define thresholds
thresholds <- c(0.1, 0.2, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

head(d)

# re-define result for each threshold
data_thresh <- map(thresholds, function(thresh) {
  d |>
    select(age, sex, tb1, tb2, nil, mitogen) |>
    mutate(
      threshold = thresh,
      QFT_reclass = QFT_thresh_fn(tb1, tb2, nil, mitogen, threshold = thresh),
      # calculate nil-subtracted antigen response
      tb1_nil = tb1 - nil,
      tb2_nil = tb2 - nil,
      mitogen_nil = mitogen - nil
    )
})

# Combine all results
data_thresh <- bind_rows(data_thresh)


# prepare model data - exclude indeterminates

data_thresh <- data_thresh |>
  filter(QFT_reclass %in% c("Positive", "Negative")) |>
  mutate(
    IGRA = as.integer(QFT_reclass == "Positive"),
    # sex_binary = as.integer(sex == "Male"),  # Adjust based on your sex coding
    # age_scaled = scale(age)[,1],  # Standardize age for better sampling
    # age_centered = age - mean(age, na.rm = TRUE)
  )

# Fit models for each threshold
models <- list()

for (thresh in thresholds) {
  cat(sprintf("\nFitting model for threshold %.2f...\n", thresh))

  # filter data for threshold
  df1 <- data_thresh |>
    filter(threshold == thresh)
  
  # same priors as for QFT-Plus threshold
  model <- brm(
    IGRA ~ s(age, by = sex, k = 3),
    data = df1,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(logit(0.2), 2.5), class = "Intercept"), # ~20% prevalence at ≥0.35
      prior(normal(0, 1), class = "b")
    ),
    chains = CHAINS, iter = ITER, warmup = WARMUP, seed = SEED,
    silent = 2, refresh = 0,  # suppress iteration messages
  )
  models[[as.character(thresh)]] <- model
}

