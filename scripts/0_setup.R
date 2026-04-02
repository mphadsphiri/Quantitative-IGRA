# Title: Functions
# Author: Mphatso
# Date: 2025-05-10 (last edited 2026-04-02)

# Objective: 

# Set up --------------
## packages --------------

# packages
library(here)
library(dplyr)
library(tidyr)
library(forcats)
library(purrr)
library(tidyverse)
library(kableExtra)
# library(gtsummary)
# library(gt)
library(ggplot2)
library(posterior)
library(rstan)
library(cmdstanr)
library(brms)
library(patchwork)
library(tidybayes)
library(scales)
library(ggbreak)

# Cite ggbreak
# S Xu, M Chen, T Feng, L Zhan, L Zhou, G Yu. Use ggbreak to effectively utilize plotting space to deal with large datasets and
# outliers. Frontiers in Genetics. 2021, 12:774846. doi: 10.3389/fgene.2021.774846 

# settings
options(mc.cores = parallel::detectCores())
options(file_refit = "always")
options(brms.backend = "cmdstanr")
rstan_options(auto_write = FALSE)

# model fitting settings
ITER <- 2000; WARMUP <- ITER * 0.5; CHAINS <- CORES <- 4; 
CONTROL <- list(adapt_delta = 0.999, max_treedepth = 15)

SEED <- 1257
NDRAWS <- 1000

theme_set(theme_bw(base_size = 11))

# data

d <- readRDS(here("data/data.rds"))

d <- d |> 
  select(age, age_5y, age_10y, sex, hiv_status, artstatus, fever, cough_2wks, tbrxhx2, hh_tbrx2,
         tb1_raw, tb2_raw, nil_raw, mitogen_raw,
         tb1, tb2, nil, mitogen, tb1_nil, tb2_nil, mitogen_nil,
         QFTP_result_final)

# prepare data
d <- d |> 
  
  # redefine nil-subtracted TB1 & TB2 as ≥ 0 based on QIAGEN threshold (0.35 IU/ml)
  mutate(tb1_nil_qia = case_when(tb1_nil < 0.35 ~ 0, TRUE ~ tb1_nil),
         tb2_nil_qia = case_when(tb2_nil < 0.35 ~ 0, TRUE ~ tb2_nil)) |> 
  
  # TB2-TB1 differential
  mutate(tb2tb1_diff = tb2_nil - tb1_nil) |> 
  
  # TB2-TB1 differential > 0.6 IU/ml
  mutate(tb2tb1_06 = ifelse(tb2tb1_diff > 0.6, 1, 0)) |> 
  mutate(tb2tb1_06_cat = ifelse(tb2tb1_06 == 1, "Yes", "No")) |> 
  mutate(tb2tb1_06_cat = factor(tb2tb1_06_cat, levels = c("Yes", "No"))) |> 

  # tube-specific result
  
  mutate(tb1_result = case_when(nil > 8 ~ "Indeterminate",
                                tb1_nil >= 0.35 & tb1_nil >= 0.25 * nil ~ "Positive",
                                mitogen-nil < 0.5 & (tb1_nil < 0.35 | (tb1_nil >= 0.35 & tb1_nil < 0.25 * nil)) ~ "Indeterminate",
                                mitogen-nil >= 0.5 & (tb1_nil < 0.35 | (tb1_nil >= 0.35 & tb1_nil < 0.25 * nil)) ~ "Negative")) |>
  
  mutate(tb2_result = case_when(nil > 8 ~ "Indeterminate",
                                tb2_nil >= 0.35 & tb2_nil >= 0.25 * nil ~ "Positive",
                                mitogen-nil < 0.5 & (tb2_nil < 0.35 | (tb2_nil >= 0.35 & tb2_nil < 0.25 * nil)) ~ "Indeterminate",
                                mitogen-nil >= 0.5 & (tb2_nil < 0.35 | (tb2_nil >= 0.35 & tb2_nil < 0.25 * nil)) ~ "Negative")) |>
  
  mutate(reason_indet = case_when(tb1_result == "Indeterminate" & nil > 8 ~ "High Nil",
                                  tb1_result == "Indeterminate" & mitogen - nil < 0.5 ~ "Low Mitogen",
                                  tb2_result == "Indeterminate" & nil > 8 ~ "High Nil",
                                  tb2_result == "Indeterminate" & mitogen - nil < 0.5 ~ "Low Mitogen",
                                  tb1_result != "Indeterminate" &  tb1_result != "Indeterminate"~ NA_character_)) |>
  
  mutate(
    tb1_neg = tb1_nil < 0.35 | (tb1_nil >= 0.35 & tb1_nil < 0.25 * nil),
    tb2_neg = tb2_nil < 0.35 | (tb2_nil >= 0.35 & tb2_nil < 0.25 * nil),
    
    QFT_result_recalc = case_when(
      # Nil very high
      nil > 8 ~ "Indeterminate",
      
      # Positive if TB1 or TB2 positive
      (tb1_nil >= 0.35 & tb1_nil >= 0.25 * nil) |
        (tb2_nil >= 0.35 & tb2_nil >= 0.25 * nil) ~ "Positive",
      
      # Both TB1 & TB2 negative + low mitogen
      tb1_neg & tb2_neg & mitogen-nil < 0.5 ~ "Indeterminate",
      
      # Both TB1 & TB2 negative + adequate mitogen
      tb1_neg & tb2_neg & mitogen-nil >= 0.5 ~ "Negative",
      
      TRUE ~ NA_character_
    )
  ) |> 
  
  mutate(QFTP_result_check = case_when(tb1_result == "Positive" | tb2_result == "Positive" ~ "Positive",
                                       tb1_result == "Negative" & tb2_result == "Negative" ~ "Negative",
                                       TRUE ~ "Indeterminate")) |> 
  
  mutate(tube_positive = case_when(tb1_result == "Positive" & tb2_result == "Positive" ~ "Both TB1 & TB2",
                                   tb1_result == "Positive" & tb2_result != "Positive" ~ "TB1 only",
                                   tb1_result != "Positive" & tb2_result == "Positive" ~ "TB2 only",
                                   tb1_result != "Positive" & tb2_result != "Positive" ~ "None")) |> 
  
  mutate(tube_positive = factor(tube_positive, levels = c("Both TB1 & TB2", "TB1 only", "TB2 only", "None"))) |> 
  
  # tube result combined
  mutate(tube_result_concord = case_when(tb1_result == "Positive" & tb2_result == "Positive" ~ "Both TB1 & TB2 positive",
                                         tb1_result == "Negative" & tb2_result == "Negative" ~ "Both TB1 & TB2 negative",
                                         tb1_result == "Indeterminate" & tb2_result == "Indeterminate" ~ "Both TB1 & TB2 indeterminate",
                                         
                                         tb1_result == "Positive" & tb2_result == "Negative" ~ "TB1 positive & TB2 negative",
                                         tb1_result == "Positive" & tb2_result == "Indeterminate" ~ "TB1 positive & TB2 indeterminate",
                                         
                                         tb1_result == "Negative" & tb2_result == "Positive" ~ "TB1 negative & TB2 positive",
                                         tb1_result == "Negative" & tb2_result == "Indeterminate" ~ "TB1 negative & TB2 indeterminate",
                                         
                                         tb1_result == "Indeterminate" & tb2_result == "Positive" ~ "TB1 indeterminate & TB2 positive",
                                         tb1_result == "Indeterminate" & tb2_result == "Negative" ~ "TB1 indeterminate & TB2 negative"
                                         
  )) |>
  
  mutate(tube_result_concord = factor(tube_result_concord, levels = c("Both TB1 & TB2 positive",
                                                                      "Both TB1 & TB2 negative",
                                                                      "Both TB1 & TB2 indeterminate",
                                                                      "TB1 positive & TB2 negative",
                                                                      "TB1 negative & TB2 positive",
                                                                      "TB1 positive & TB2 indeterminate",
                                                                      "TB1 indeterminate & TB2 positive")))



# Session Info
sessionInfo()
