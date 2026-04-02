# Title: Quantitative IGRA model analyses
# Author: Mphatso Phiri
# Date: 2025-05-10 (last edited 2026-04-02)

# summarise correlation parameter
summarise_corr_fn <- function(x) {
  model_name <- deparse(substitute(x))
  as_draws_df(x, variable = "rescor__tb1nil__tb2nil") |>
    summarise_draws(mean, sd, rhat, ess_bulk, \(x) posterior::quantile2(x, probs = c(0.025, 0.5, 0.975))) |> 
    mutate(model = model_name, .before = everything()) |> 
    mutate(`Correlation (95% CrI)` = sprintf("%.2f (%.2f-%.2f)", round(mean, 2), round(`q2.5`, 2), round(`q97.5`, 2))) |> 
    select(model, `Correlation (95% CrI)`, mean, sd, rhat, ess_bulk) |> 
    rename(Group = model, Mean = mean, SD = sd, Rhat = rhat, ESS_bulk = ess_bulk) |> 
    mutate(Group = case_when(
      Group == "m0_tb1tb2_corr_overall" ~ "Overall",
      Group == "m0_tb1tb2_corr_1019"    ~ "Age: 10-19",
      Group == "m0_tb1tb2_corr_2029"    ~ "Age: 20-29",
      Group == "m0_tb1tb2_corr_3040"    ~ "Age: 30-40",
      Group == "m0_tb1tb2_corr_sexMale"  ~ "Sex: Male",
      Group == "m0_tb1tb2_corr_sexFemale"  ~ "Sex: Female",
      TRUE                             ~ Group
    ))
}

# summarise ORs

# function to summarise effects
summarise_or <- function(fit, coef_keep=coefs, rename_map=coef_rename_map, probs = c(0.025, 0.975)) {
  model_name <- rlang::as_label(rlang::enquo(fit))   # extracts the variable name
  
  fit |>
    as_draws_df() |>
    summarise_draws(
      mean, sd, mad, rhat,
      q = \(x) posterior::quantile2(x, probs = probs),
      ess_bulk
    ) |>
    dplyr::filter(variable %in% coef_keep) |>
    dplyr::mutate(
      # decide on the ORIGINAL brms parameter name
      .hu_plogis = stringr::str_detect(variable, "Intercept_hu"),
      .sigma = stringr::str_detect(variable, "^sigma"),
      
      # transform
      across(c(mean, `q2.5`, `q97.5`), \(x) ifelse(.hu_plogis, plogis(x),
                                                   ifelse(.sigma, (x), exp(x)))),
      
      # now rename for presentation
      variable = dplyr::recode(variable, !!!rename_map),
      
      across(where(is.numeric), \(x) round(x, 2)),
      or = mean,
      or_ci = sprintf("%.2f (%.2f-%.2f)", mean, `q2.5`, `q97.5`)
    ) |>
    dplyr::transmute(
      model = model_name,
      variable, or, q2.5, q97.5, sd, rhat, ess_bulk, or_ci
    )
  
}

# coef names 
coefs <- c(
  "b_hiv_statusPositive",
  "b_tbrxhx2Yes",
  "b_hh_tbrx2Yes",
  "Intercept",
  "Intercept_hu",
  "sigma",
  
  
  # Joint
  "b_tb1nilqia_hiv_statusPositive",
  "b_tb1nilqia_tbrxhx2Yes",
  "b_tb1nilqia_hh_tbrx2Yes",
  "b_tb1nilqia_Intercept",
  "b_tb2nilqia_Intercept",
  "Intercept_tb1nilqia",
  "Intercept_tb2nilqia",
  "Intercept_hu_tb1nilqia",
  "Intercept_hu_tb2nilqia",
  "sigma_tb1nilqia",
  "sigma_tb2nilqia"
)

# new coef names
coef_rename_map <- c(
  # new name = old name
  
  # univariate models
  "Intercept" = "Intercept",
  "HIV positive" = "b_hiv_statusPositive",
  "Past TB treatment" = "b_tbrxhx2Yes",
  "Household TB contact" = "b_hh_tbrx2Yes",
  "Hurdle intercept" = "Intercept_hu",
  "Sigma" = "sigma",
  
  # joint
  "HIV Positive" = "b_tb1nilqia_hiv_statusPositive",
  "Past TB treatment" = "b_tb1nilqia_tbrxhx2Yes",
  "Household TB contact" = "b_tb1nilqia_hh_tbrx2Yes",
  
  "TB1: Intercept" = "Intercept_tb1nilqia",
  "TB2: Intercept" = "Intercept_tb2nilqia",
  "TB1: Hurdle Intercept" = "Intercept_hu_tb1nilqia",
  "TB2 Hurdle Intercept" = "Intercept_hu_tb2nilqia",
  "TB1: Sigma" = "sigma_tb1nilqia",
  "TB2: Sigma" = "sigma_tb2nilqia"
  
)

# coef order
coef_order <- c(
  # Univariate models
  "HIV positive",
  "Past TB treatment",
  "Household TB contact",
  "Intercept",
  "Sigma",
  "Hurdle intercept",
  
  # Joint
  "HIV Positive",
  "Past TB treatment",
  "Household TB contact",
  
  "TB1: Intercept",
  "TB2: Intercept",
  "TB1: Hurdle Intercept",
  "TB2 Hurdle Intercept",
  "TB1: Sigma",
  "TB2: Sigma"
)


# re-calculate QFT-Plus results based on hypothetical thresholds
QFT_thresh_fn <- function(tb1, tb2, nil, mitogen, threshold = 0.35) {
  tb1_nil <- tb1 - nil
  tb2_nil <- tb2 - nil
  mitogen_nil <- mitogen - nil
  
  indet_nil <- nil > 8.0
  
  tb1_positive <- (tb1_nil >= threshold) & (tb1_nil >= 0.25 * nil)
  tb2_positive <- (tb2_nil >= threshold) & (tb2_nil >= 0.25 * nil)
  antigen_positive <- tb1_positive | tb2_positive
  
  mitogen_ok <- mitogen_nil >= 0.5
  
  case_when(
    indet_nil ~ "Indeterminate",
    antigen_positive ~ "Positive",
    !antigen_positive & mitogen_ok ~ "Negative",
    TRUE ~ "Indeterminate"
  )
}

