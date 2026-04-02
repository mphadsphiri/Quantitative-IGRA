# Title: Quantitative IGRA - descriptive
# Author: Mphatso Phiri
# Date: 2025-05-10 (last edited 2026-04-02)

# Objective: descriptive analysis

## data --------------

# prepare data for desc analyses

df <- d
# recoded values "> 10¶" or "> 10#" to 10.00.
# '¶' result outside ELISA range
# '#' non-numerical OD value entered therefore assigned OD value 4

# original values >10 recoded
sum(df$tb1_raw == "> 10¶" | df$tb1_raw == "> 10#") # 48
mean(df$tb1_raw =="> 10¶" | df$tb1_raw == "> 10#") # 1.7%

sum(df$tb2_raw == "> 10¶" | df$tb2_raw == "> 10#") # 13
mean(df$tb2_raw == "> 10¶" | df$tb2_raw =="> 10#") # 0.4%

# QFT-Plus positive, negative or indeterminate
table(df$QFTP_result_final)

df <- df |> 
  filter(QFTP_result_final %in% c("Positive", "Negative", "Indeterminate"))

table(df$QFTP_result_final)

# recode NA to explicing "Missing" & relevel sex for descriptive tables 
desc_vars <- c("fever", "cough_2wks", "hiv_status", "artstatus",
               "tbrxhx2", "hh_tbrx2")
df <- df|>
  mutate(across(all_of(desc_vars),
                ~ fct_explicit_na(as.factor(.), na_level = "Missing")))  # recode NA

levels(df$hiv_status) # confirm missing

# df <- df |> mutate(sex = fct_relevel(sex, c("Male", "Female"))) # sex relevel
# 
# levels(df$sex) # confirm levels


# Participant characteristics -----------------------------------------
## overall

tb1_overall <- df |> 
  summarise(
    
    # n
    n = n(),
    
    # age
    mean_age = sprintf("%.1f", mean(age)),
    sd_age = sprintf("%.1f", sd(age)),
    median_age = sprintf("%.1f", median(age)),
    p25_age = sprintf("%.1f", quantile(age, 0.25)),
    p75_age = sprintf("%.1f", quantile(age, 0.75)),
    
    # sex
    male = sum(sex %in% "Male"),
    female = sum(sex %in% "Female"),
    sex_missing = sum(sex %in% "Missing"),
    
    #  HIV
    hiv_positive = sum(hiv_status %in% "Positive"),
    hiv_negative = sum(hiv_status %in% "Negative"),
    hiv_missing = sum(hiv_status %in% "Missing"),
    
    #  ART
    art_yes = sum(artstatus %in% "Yes" & hiv_status %in% "Positive"),
    art_no= sum(artstatus %in% "No" & hiv_status %in% "Positive"),
    art_missing = sum(artstatus %in% c("Missing") & hiv_status %in% "Positive"),
    
    # cough ≥ 2 weeks
    cough_2wks_yes = sum(cough_2wks %in% "Yes"),
    cough_2wks_no = sum(cough_2wks %in% "No"),
    cough_2wks_missing = sum(cough_2wks %in% "Missing"),
    
    # fever
    fever_yes = sum(fever %in% "Yes"),
    fever_no = sum(fever %in% "No"),
    fever_missing = sum(fever %in% "Missing"),
    
    # TB treatment
    tbrxhx2_yes = sum(tbrxhx2 %in% "Yes"),
    tbrxhx2_no = sum(tbrxhx2 %in% "No"),
    tbrxhx2_missing = sum(tbrxhx2 %in% "Missing"),
    
    # Household TB contact
    hh_tbrx2_yes = sum(hh_tbrx2 %in% "Yes"),
    hh_tbrx2_no = sum(hh_tbrx2 %in% "No"),
    hh_tbrx2_missing = sum(hh_tbrx2 %in% "Missing"),
    

    # QFT-Plus
    QFT_positive = sum(QFTP_result_final %in% "Positive"),
    QFT_negative = sum(QFTP_result_final %in% "Negative"),
    QFT_indet = sum(QFTP_result_final %in% "Indeterminate"),
    
    
  ) |> 
  
  # summaries
  mutate(
    # n
    n_prop = round(100 * n / n, 1),
    n_prop2 = paste0(sprintf("%.1f", n_prop), "%", " (", n, ")"),
    
    # age
    mean_age_sd = paste0(mean_age, " (", sd_age, ")"), 
    median_age_IQR = paste0(median_age, " (", p25_age, ", ", p75_age, ")"), 
    
    # sex
    male_prop = round(100 * male / (male + female + sex_missing), 1), 
    male_prop2 = paste0(male_prop, "%", " ", "(", male, "/", male + female, ")"), 
    
    female_prop = round(100 * female / (male + female + sex_missing), 1), 
    female_prop2 = paste0(female_prop, "%", " ", "(", male, "/", male + female, ")"), 
    
    sex_missing_prop = round(100 * sex_missing / (male + female + sex_missing), 1),
    sex_missing_prop2 = paste0(sex_missing_prop, "%", " ", "(", sex_missing, "/", male + female + sex_missing, ")"), 
    
    # HIV 
    
    hiv_prop = round(100 * hiv_positive / (hiv_positive + hiv_negative), 1), 
    hiv_prop2 = paste0(hiv_prop, "%", " ", "(", hiv_positive, "/", hiv_positive + hiv_negative, ")"), 
    
    hiv_neg_prop = round(100 * hiv_negative / (hiv_positive + hiv_negative), 1), 
    hiv_neg_prop2 = paste0(hiv_neg_prop, "%", " ", "(", hiv_negative, "/", hiv_positive + hiv_negative, ")"), 
    
    hiv_missing_prop = round(100 * hiv_missing / (hiv_positive + hiv_negative + hiv_missing), 1),
    hiv_missing_prop2 = paste0(hiv_missing_prop, "%", " ", "(", hiv_missing, "/", hiv_positive + hiv_negative + hiv_missing, ")"), 
    
    # ART
    art_yes_prop = round(100 * art_yes / (art_yes + art_no), 1), 
    art_yes_prop2 = paste0(art_yes_prop, "%", " ", "(", art_yes, "/", art_yes + art_no, ")"), 
    
    art_missing_prop = round(100 * art_missing / (art_yes + art_no + art_missing), 1),
    art_missing_prop2 = paste0(art_missing_prop, "%", " ", "(", art_missing, "/", art_yes + art_no + art_missing, ")"), 
    
    # cough ≥ 2 weeks
    
    cough_2wks_prop = round(100 * cough_2wks_yes / (cough_2wks_yes + cough_2wks_no), 1),
    cough_2wks_prop2 = paste0(cough_2wks_prop, "%", " ", "(", cough_2wks_yes, "/", cough_2wks_yes + cough_2wks_no, ")"),
    
    cough_2wks_missing_prop = round(100 * cough_2wks_missing / (cough_2wks_yes + cough_2wks_no + cough_2wks_missing), 1),
    cough_2wks_missing_prop2 = paste0(cough_2wks_missing_prop, "%", " ", "(", cough_2wks_missing, "/", cough_2wks_yes + cough_2wks_no + cough_2wks_missing, ")"), 
    
    # fever
    fever_prop = round(100 * fever_yes / (fever_yes + fever_no), 1), 
    fever_prop2 = paste0(fever_prop, "%", " ", "(", fever_yes, "/", fever_yes + fever_no, ")"), 
    
    fever_missing_prop = round(100 * fever_missing / (fever_yes + fever_no + fever_missing), 1),
    fever_missing_prop2 = paste0(fever_missing_prop, "%", " ", "(", fever_missing, "/", fever_yes + fever_no + fever_missing, ")"), 
    
    # TB treatment
    
    tbrxhx2_prop = round(100 * tbrxhx2_yes / (tbrxhx2_yes + tbrxhx2_no), 1),
    tbrxhx2_prop2 = paste0(tbrxhx2_prop, "%", " ", "(", tbrxhx2_yes, "/", tbrxhx2_yes + tbrxhx2_no, ")"), 
    
    tbrxhx2_missing_prop = round(100 * tbrxhx2_missing / (tbrxhx2_yes + tbrxhx2_no + tbrxhx2_missing), 1),
    tbrxhx2_missing_prop2 = paste0(tbrxhx2_missing_prop, "%", " ", "(", tbrxhx2_missing, "/", tbrxhx2_yes + tbrxhx2_no + tbrxhx2_missing, ")"), 
    
    # Household TB contact
    hh_tbrx2_prop = round(100 * hh_tbrx2_yes / (hh_tbrx2_yes + hh_tbrx2_no), 1), 
    hh_tbrx2_prop2 = paste0(hh_tbrx2_prop, "%", " ", "(", hh_tbrx2_yes, "/", hh_tbrx2_no, ")"), 
    
    hh_tbrx2_missing_prop = round(100 * hh_tbrx2_missing / (hh_tbrx2_yes + hh_tbrx2_no + hh_tbrx2_missing), 1),
    hh_tbrx2_missing_prop2 = paste0(hh_tbrx2_missing_prop, "%", " ", "(", hh_tbrx2_missing, "/", hh_tbrx2_yes + hh_tbrx2_no + hh_tbrx2_missing, ")"), 
    
    # QFT-Plus
    
    QFT_pos_prop = round(100 * QFT_positive / n, 1),
    QFT_pos_prop2 = paste0(sprintf("%.1f", QFT_pos_prop), "%", " ", "(", QFT_positive, "/", n, ")"),
    
    QFT_neg_prop = round(100 * QFT_negative / n, 1),
    QFT_neg_prop2 = paste0(sprintf("%.1f", QFT_neg_prop), "%", " ", "(", QFT_negative, "/", n, ")"),
    
    QFT_indet_prop = round(100 * QFT_indet / n, 1),
    QFT_indet_prop2 = paste0(sprintf("%.1f", QFT_indet_prop), "%", " ", "(", QFT_indet, "/", n, ")")
    
  )


# overall wide

tb1_overall_long <- tb1_overall |> 
  mutate(across(everything(), as.character)) |>
  mutate(age_10y = "Overall") |>
  pivot_longer(cols = -c(age_10y), names_to = "Characteristic", values_to = "Overall") 

# overall long

tb1_overall_long

## age (10-year) group

tb1_age_10y <- df |> 
  group_by(age_10y) |> 
  
  summarise(
    # n
    n = n(),
    
    # age
    mean_age = sprintf("%.1f", mean(age)),
    sd_age = sprintf("%.1f",sd(age)),
    median_age = sprintf("%.1f", median(age)), 
    p25_age = sprintf("%.1f", quantile(age, 0.25)), 
    p75_age = sprintf("%.1f", quantile(age, 0.75)),
    
    # sex
    male = sum(sex %in% "Male"),
    female = sum(sex %in% "Female"),
    sex_missing = sum(sex %in% "Missing"),
    
    #  HIV
    hiv_positive = sum(hiv_status %in% "Positive"),
    hiv_negative = sum(hiv_status %in% "Negative"),
    hiv_missing = sum(hiv_status %in% "Missing"),
    
    #  ART
    art_yes = sum(artstatus %in% "Yes" & hiv_status %in% "Positive"),
    art_no= sum(artstatus %in% "No" & hiv_status %in% "Positive"),
    art_missing = sum(artstatus %in% c("Missing") & hiv_status %in% "Positive"),
    
    # fever
    fever_yes = sum(fever %in% "Yes"),
    fever_no = sum(fever %in% "No"),
    fever_missing = sum(fever %in% "Missing"),
    
    # cough ≥ 2 weeks
    cough_2wks_yes = sum(cough_2wks %in% "Yes"),
    cough_2wks_no = sum(cough_2wks %in% "No"),
    cough_2wks_missing = sum(cough_2wks %in% "Missing"),
    
    # TB treatment
    tbrxhx2_yes = sum(tbrxhx2 %in% "Yes"),
    tbrxhx2_no = sum(tbrxhx2 %in% "No"),
    tbrxhx2_missing = sum(tbrxhx2 %in% "Missing"),
    
    # Household TB contact
    hh_tbrx2_yes = sum(hh_tbrx2 %in% "Yes"),
    hh_tbrx2_no = sum(hh_tbrx2 %in% "No"),
    hh_tbrx2_missing = sum(hh_tbrx2 %in% "Missing"),
    
    # QFT-Plus
    QFT_positive = sum(QFTP_result_final %in% "Positive"),
    QFT_negative = sum(QFTP_result_final %in% "Negative"),
    QFT_indet = sum(QFTP_result_final %in% "Indeterminate"),
    
  ) |> 
  
  # summaries
  mutate(
    # n
    n_prop = round(100 * n / sum(n), 1),
    n_prop2 = paste0(sprintf("%.1f", n_prop), "%", " (", n, ")"),
    
    # age
    mean_age_sd = paste0(mean_age, " (", sd_age, ")"), 
    median_age_IQR = paste0(median_age, " (", p25_age, ", ", p75_age, ")"), 
    
    # sex
    male_prop = round(100 * male / (male + female + sex_missing), 1), 
    male_prop2 = paste0(male_prop, "%", " ", "(", male, "/", male + female, ")"), 
    
    female_prop = round(100 * female / (male + female + sex_missing), 1), 
    female_prop2 = paste0(female_prop, "%", " ", "(", male, "/", male + female, ")"), 
    
    sex_missing_prop = round(100 * sex_missing / (male + female + sex_missing), 1),
    sex_missing_prop2 = paste0(sex_missing_prop, "%", " ", "(", sex_missing, "/", male + female + sex_missing, ")"), 
    
    
    # HIV 
    
    hiv_prop = round(100 * hiv_positive / (hiv_positive + hiv_negative), 1), 
    hiv_prop2 = paste0(hiv_prop, "%", " ", "(", hiv_positive, "/", hiv_positive + hiv_negative, ")"), 
    
    hiv_neg_prop = round(100 * hiv_negative / (hiv_positive + hiv_negative), 1), 
    hiv_neg_prop2 = paste0(hiv_neg_prop, "%", " ", "(", hiv_negative, "/", hiv_positive + hiv_negative, ")"), 
    
    hiv_missing_prop = round(100 * hiv_missing / (hiv_positive + hiv_negative + hiv_missing), 1),
    hiv_missing_prop2 = paste0(hiv_missing_prop, "%", " ", "(", hiv_missing, "/", hiv_positive + hiv_negative + hiv_missing, ")"), 
    
    # ART
    art_yes_prop = round(100 * art_yes / (art_yes + art_no), 1), 
    art_yes_prop2 = paste0(art_yes_prop, "%", " ", "(", art_yes, "/", art_yes + art_no, ")"), 
    
    art_missing_prop = round(100 * art_missing / (art_yes + art_no + art_missing), 1),
    art_missing_prop2 = paste0(art_missing_prop, "%", " ", "(", art_missing, "/", art_yes + art_no + art_missing, ")"), 
    
    # cough ≥ 2 weeks
    
    cough_2wks_prop = round(100 * cough_2wks_yes / (cough_2wks_yes + cough_2wks_no), 1),
    cough_2wks_prop2 = paste0(cough_2wks_prop, "%", " ", "(", cough_2wks_yes, "/", cough_2wks_yes + cough_2wks_no, ")"),
    
    cough_2wks_missing_prop = round(100 * cough_2wks_missing / (cough_2wks_yes + cough_2wks_no + cough_2wks_missing), 1),
    cough_2wks_missing_prop2 = paste0(cough_2wks_missing_prop, "%", " ", "(", cough_2wks_missing, "/", cough_2wks_yes + cough_2wks_no + cough_2wks_missing, ")"), 
    
    # fever
    fever_prop = round(100 * fever_yes / (fever_yes + fever_no), 1), 
    fever_prop2 = paste0(fever_prop, "%", " ", "(", fever_yes, "/", fever_yes + fever_no, ")"), 
    
    fever_missing_prop = round(100 * fever_missing / (fever_yes + fever_no + fever_missing), 1),
    fever_missing_prop2 = paste0(fever_missing_prop, "%", " ", "(", fever_missing, "/", fever_yes + fever_no + fever_missing, ")"), 
    
    # TB treatment
    
    tbrxhx2_prop = round(100 * tbrxhx2_yes / (tbrxhx2_yes + tbrxhx2_no), 1),
    tbrxhx2_prop2 = paste0(tbrxhx2_prop, "%", " ", "(", tbrxhx2_yes, "/", tbrxhx2_yes + tbrxhx2_no, ")"), 
    
    tbrxhx2_missing_prop = round(100 * tbrxhx2_missing / (tbrxhx2_yes + tbrxhx2_no + tbrxhx2_missing), 1),
    tbrxhx2_missing_prop2 = paste0(tbrxhx2_missing_prop, "%", " ", "(", tbrxhx2_missing, "/", tbrxhx2_yes + tbrxhx2_no + tbrxhx2_missing, ")"), 
    
    # Household TB contact
    hh_tbrx2_prop = round(100 * hh_tbrx2_yes / (hh_tbrx2_yes + hh_tbrx2_no), 1), 
    hh_tbrx2_prop2 = paste0(hh_tbrx2_prop, "%", " ", "(", hh_tbrx2_yes, "/", hh_tbrx2_no, ")"), 
    
    hh_tbrx2_missing_prop = round(100 * hh_tbrx2_missing / (hh_tbrx2_yes + hh_tbrx2_no + hh_tbrx2_missing), 1),
    hh_tbrx2_missing_prop2 = paste0(hh_tbrx2_missing_prop, "%", " ", "(", hh_tbrx2_missing, "/", hh_tbrx2_yes + hh_tbrx2_no + hh_tbrx2_missing, ")"), 
    
    # QFT-Plus
    
    QFT_pos_prop = round(100 * QFT_positive / n, 1),
    QFT_pos_prop2 = paste0(sprintf("%.1f", QFT_pos_prop), "%", " ", "(", QFT_positive, "/", n, ")"),
    
    QFT_neg_prop = round(100 * QFT_negative / n, 1),
    QFT_neg_prop2 = paste0(sprintf("%.1f", QFT_neg_prop), "%", " ", "(", QFT_negative, "/", n, ")"),
    
    QFT_indet_prop = round(100 * QFT_indet / n, 1),
    QFT_indet_prop2 = paste0(sprintf("%.1f", QFT_indet_prop), "%", " ", "(", QFT_indet, "/", n, ")")
    
  )

tb1_age_10y


# age wide
tb1_age_10y_long <- tb1_age_10y |> 
  mutate(across(everything(), as.character)) |>
  
  pivot_longer(-c(age_10y), names_to = "Characteristic", values_to = "value") |> 
  
  ungroup() |> 
  
  pivot_wider(id_cols = Characteristic, names_from = c(age_10y), values_from = value)


# Combined overall & age group 

tb1_combined <- left_join(tb1_overall_long |>
                            select(-age_10y),
                          tb1_age_10y_long)

tb1_combined


# select variables for Table 1
tb1 <- tb1_combined |> 
  filter(Characteristic %in% c(
    "n_prop2",
    "median_age_IQR",
    "male_prop2",
    # "female_prop2",
    "hiv_prop2",
    "art_yes_prop2",
    "fever_prop2",
    "diarrhea_prop2",
    "cough_prop2",
    "cough_2wks_prop2",
    "any_tb_sympts_prop2",
    "tbrxhx2_prop2",
    "hh_tbrx2_prop2",
    "QFT_pos_prop2",
    "QFT_neg_prop2",
    "QFT_indet_prop2")
    
  ) |> 
  
  mutate(Characteristic = case_match(Characteristic,
                                     "n_prop2" ~ "n (%)", # old ~ new
                                     "median_age_IQR" ~ "Age years, median (IQR)", # old ~ new
                                     "male_prop2" ~ "Male",
                                     # "female_prop2" ~ "Female",
                                     "hiv_prop2" ~ "HIV positive",
                                     "art_yes_prop2" ~ "HIV positive,n ART",
                                     "fever_prop2" ~ "Fever",
                                     "diarrhea_prop2" ~ "Diarrhea",
                                     "cough_prop2" ~ "Any cough",
                                     "cough_2wks_prop2" ~ "Cough ≥ 2 weeks",
                                     "any_tb_sympts_prop2" ~ "Any TB symptoms",
                                     "tbrxhx2_prop2" ~ "Previous TB treatment",
                                     "hh_tbrx2_prop2" ~ "Household TB contact",
                                     "QFT_pos_prop2" ~ "QFT-Plus positive",
                                     "QFT_neg_prop2" ~ "QFT-Plus negative",
                                     "QFT_indet_prop2" ~ "QFT-Plus indeterminate")
         
  )



tb1


# Figure 1: TB1 & TB2 distributions --------------------------------------------

# box plot
f1 <- df |> 
  mutate(sex = fct_rev(sex)) |> 
  select(age_5y, sex, tb1_nil, tb2_nil) |> 
  pivot_longer(cols = c(tb1_nil, tb2_nil), names_to = "tube", values_to = "value") |> 
  ggplot() +
  # geom_point(aes(age_5y, value, color=sex), alpha=0.3) + # , position = position_jitter()
  geom_boxplot(aes(age_5y, value, color = sex, fill=sex), alpha=0.3, position = position_dodge()) +
  # geom_smooth(aes(factor(age), value, group = sex, color = sex), method = "loess") +
  facet_wrap(~ifelse(tube=="tb1_nil", "TB1","TB2")) +
  labs(x = "Age (years)",
       y = "IFN-γ (ΙU/mL)"
       # title = "Distribution of QuantiFERON-TB Gold Plus Nil-subtracted TB1 & TB2 values by age & sex"
  ) +
  
  scale_y_continuous(limits = c(-5,10)) +
  # scale_y_log10() +
  
  scale_color_brewer(name = "Sex", type = "qual", palette = "Set1") +
  scale_fill_brewer(name = "Sex", type = "qual", palette = "Set1") +
  # scale_x_discrete(expand = c(rep(0.0004,4))) +
  
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  theme_bw() +
  theme(#legend.position = c(0.8,0.2),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 12),
    plot.subtitle = element_text(size = 11),
    legend.position = "bottom")


f1


# Distribution of continuous TB1 & TB2 distributions -----------------

# predictions


# new data

df_new <- expand_grid(
  age = seq(10, 40, by = 0.01),
  sex = factor(c("Female", "Male"), levels = c("Female", "Male"))
)


# linear predictor

# TB1
df_pred_mu_tb1 <- add_linpred_draws(
  newdata = df_new,
  object = m0_tb1_nil_qia_hurdle_age_sex_sp,
  ndraws = NDRAWS) |> 
  mutate(tube = "TB1")

# TB2
df_pred_mu_tb2 <- add_linpred_draws(
  newdata = df_new,
  object = m0_tb2_nil_qia_hurdle_age_sex_sp,
  ndraws = NDRAWS) |> 
  mutate(tube = "TB2")

# join
df_pred_mu <- df_pred_mu_tb1 |> 
  select(tube, .draw, .row, age, sex, .linpred) |> 
  rename(.linpred1 = .linpred) |>
  left_join(df_pred_mu_tb2 |> 
              select(tube, .draw, .row, age, sex, .linpred) |> 
              rename(.linpred2 = .linpred),
            by = c(".draw", ".row", "age", "sex")) |> 
  pivot_longer(cols = c(.linpred1, .linpred2), names_to = "tube", values_to = ".linpred") |> 
  mutate(mu = exp(.linpred))


# summarise

df_pred_mu <- df_pred_mu |> 
  group_by(tube, age, sex) |> 
  mean_qi(mu, .width = 0.95) |>
  ungroup()

# plot

f2 <- df_pred_mu |> 
  mutate(sex = fct_rev(sex)) |> 
  ggplot(aes(age, mu)) +
  geom_line(aes(color = sex), linewidth = 0.8) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = sex), alpha = 0.25) +
  coord_cartesian(ylim = c(0, 10)) +
  scale_y_continuous(
    breaks = seq(0, 10, by = 2)) +
  scale_color_brewer(name = "Sex", palette = "Set1") +
  scale_fill_brewer(name = "Sex", palette = "Set1") +
  facet_wrap(~ tube, labeller = as_labeller(c(.linpred1 = "TB1", .linpred2 = "TB2"))) +
  labs(
    x = "Age (years)",
    y = "Predicted mean response (IU/mL)",
    caption = "IU/mL: International units per millilitre."
    # title = "Predicted QuantiFERON-TB Gold Plus responses among participants with >0.35 IU/mL by age and sex"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 11),
    text = element_text(size = 11),
    legend.position = "bottom")


f2


# Factors associated with the magnitude of TB1 and TB2 ----------------


# TB1

# univariate
m0_tb1_nil_qia_hurdle_hiv_OR <- summarise_or(m0_tb1_nil_qia_hurdle_hiv)
m0_tb1_nil_qia_hurdle_tbrx_OR <- summarise_or(m0_tb1_nil_qia_hurdle_tbrx)
m0_tb1_nil_qia_hurdle_hhtb_OR <- summarise_or(m0_tb1_nil_qia_hurdle_hhtb)

# combine univariate
tb1_hurdle_OR <- bind_rows(
  m0_tb1_nil_qia_hurdle_hiv_OR,
  m0_tb1_nil_qia_hurdle_tbrx_OR,
  m0_tb1_nil_qia_hurdle_hhtb_OR,
)

# full

tb1_hurdle_full_OR <- summarise_or(m0_tb1_nil_qia_hurdle_full)

# combine univariate & full model
tb1_combined <- right_join(
  tb1_hurdle_OR |>
    rename(OR_CI = or_ci) |> 
    select(variable, OR_CI) |> 
    filter(!variable %in% c("Intercept", "sigma", "Intercept_hu")),
  tb1_hurdle_full_OR |>
    rename(AOR_CI = or_ci) |> 
    select(variable, AOR_CI),
  by = "variable")
  

tb1_combined

tbl2_tb1 <- tb1_combined |> 
  mutate(
    variable = case_when(
      variable == "b_hiv_statusPositive"            ~ "HIV positive",
      variable == "b_tbrxhx2Yes"                    ~ "Previous TB treatment",
      variable == "b_hh_tbrx2Yes"                   ~ "Household TB contact",
      TRUE                                          ~ variable,
    )
  ) |> 
  rename("Covariate" = variable,
         "OR (95% CrI)" = OR_CI,
         "AOR (95% CrI)" = AOR_CI)

tbl2_tb1


# TB2

# univariate
m0_tb2_nil_qia_hurdle_hiv_OR <- summarise_or(m0_tb2_nil_qia_hurdle_hiv)
m0_tb2_nil_qia_hurdle_tbrx_OR <- summarise_or(m0_tb2_nil_qia_hurdle_tbrx)
m0_tb2_nil_qia_hurdle_hhtb_OR <- summarise_or(m0_tb2_nil_qia_hurdle_hhtb)

# combined univariate
tb2_hurdle_OR <- bind_rows(
  m0_tb2_nil_qia_hurdle_hiv_OR,
  m0_tb2_nil_qia_hurdle_tbrx_OR,
  m0_tb2_nil_qia_hurdle_hhtb_OR,
)

# full

tb2_hurdle_full_OR <- summarise_or(m0_tb2_nil_qia_hurdle_full)

# combine univariate & full model
tb2_combined <- right_join(
  tb2_hurdle_OR |>
    rename(OR_CI = or_ci) |> 
    select(variable, OR_CI) |> 
    filter(!variable %in% c("Intercept", "sigma", "Intercept_hu")),
  tb2_hurdle_full_OR |>
    rename(AOR_CI = or_ci) |> 
    select(variable, AOR_CI),
  by = "variable")

tb2_combined

tbl2_tb2 <- tb2_combined |> 
  mutate(
    variable = case_when(
      variable == "b_hiv_statusPositive"            ~ "HIV positive",
      variable == "b_tbrxhx2Yes"                    ~ "Previous TB treatment",
      variable == "b_hh_tbrx2Yes"                   ~ "Household TB contact",
      TRUE                                          ~ variable,
    )
  ) |> 
  rename("Covariate" = variable,
         "OR (95% CrI)" = OR_CI,
         "AOR (95% CrI)" = AOR_CI)

tbl2_tb2


# combine TB1 & TB2 tables

tbl2 <- right_join(
  
  tbl2_tb1 |>
    rename(`TB1: OR (95% CrI` = `OR (95% CrI)`,
         `TB1: AOR (95% CrI)` = `AOR (95% CrI)`), 
  
  tbl2_tb2 |>
      rename(`TB2: OR (95% CrI)` = `OR (95% CrI)`, 
             `TB2: AOR (95% CrI)` = `AOR (95% CrI)`), 
    by = "Covariate"
  )


tbl2



# Factors associated with TB1-TB2>0.6  -----------------------------------------

# Figure 3: age & sex

## predictions

# new data
df_new <- expand_grid(
  age = seq(10, 40, by = 1),
  sex = factor(c("Female", "Male"), levels = c("Male", "Female"))
)

# predict
tb2tb1_df_pred <- m0_tb2tb1_06_agesex_sp |> 
  add_epred_draws(newdata = df_new, ndraws = NDRAWS)


## plot
f3 <- tb2tb1_df_pred |> 
  mean_qi() |>
  ggplot() + 
  
  geom_ribbon(aes(age, .epred, ymin = .lower, ymax = .upper, fill = sex), alpha = 0.3) +
  geom_line(aes(age, .epred, color = sex), linewidth = 1) +
  labs(x = "Age (years)", y = "Probability of TB2 minus TB1 ≥ 0.6 IU/mL (%)",
       fill = "Credible interval",
       caption = "IU/mL: International Units per millilitre") + 
  
  scale_fill_brewer(name = "Sex", type = "qual", palette = "Set1") +
  scale_color_brewer(name = "Sex", type = "qual", palette = "Set1") +

  scale_x_continuous(breaks = seq(10, 40, by = 5)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix = NULL)) +
  
  theme_bw()+
  theme(legend.position = c(0.12,0.85),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        plot.subtitle = element_text(size = 11))

f3


# Table 3: Factors associated with TB2-TB1>0.6 


## Table 3
# univariate
m0_tb2tb1_06_hiv_OR <- summarise_or(m0_tb2tb1_06_hiv)
m0_tb2tb1_06_tbrx_OR <- summarise_or(m0_tb2tb1_06_tbrx)
m0_tb2tb1_06_hhtb_OR <- summarise_or(m0_tb2tb1_06_hhtb)

# combined univariate
tb2tb1_06_OR <- bind_rows(
  m0_tb2tb1_06_hiv_OR,
  m0_tb2tb1_06_tbrx_OR,
  m0_tb2tb1_06_hhtb_OR
)

# full

tb2tb1_06_full_OR <- summarise_or(m0_tb2tb1_06_full)

# combine univariate & full model
tb2tb1_06_combined <- right_join(
  tb2tb1_06_OR |>
    rename(OR_CI = or_ci) |> 
    select(variable, OR_CI) |> 
    filter(!variable %in% c("Intercept", "sigma", "Intercept_hu")),
  tb2tb1_06_full_OR |> select(variable, or_ci) |> rename(AOR_CI = or_ci),
  by = "variable"
)


tbl3 <- tb2tb1_06_combined |> 
  mutate(
    variable = case_when(
      variable == "b_hiv_statusPositive"            ~ "HIV positive",
      variable == "b_tbrxhx2Yes"                    ~ "Previous TB treatment",
      variable == "b_hh_tbrx2Yes"                   ~ "Household TB contact",
      TRUE                                          ~ variable,
    )
  ) |> 
  rename("Covariate" = variable,
         "OR (95% CrI)" = OR_CI,
         "AOR (95% CrI)" = AOR_CI)

tbl3


# Tube-specific response QuantiFERON-TB Gold Plus results -------------
## overall

tb2_QFT_overall <- df |> 
  
  summarise(
    
    # n
    n = n(),
    
    tube_result_concord_both_pos = sum(tube_result_concord %in% "Both TB1 & TB2 positive"),
    tube_result_concord_both_neg = sum(tube_result_concord %in% "Both TB1 & TB2 negative"),
    
    tube_result_concord_tb1ptb2n = sum(tube_result_concord %in% "TB1 positive & TB2 negative"),
    tube_result_concord_tb1ntb2p = sum(tube_result_concord %in% "TB1 negative & TB2 positive"),
    
    tube_result_concord_both_ind = sum(tube_result_concord %in% "Both TB1 & TB2 indeterminate"),
    tube_result_concord_tb1ptb2i = sum(tube_result_concord %in% "TB1 positive & TB2 indeterminate"),
    tube_result_concord_tb1itb2p = sum(tube_result_concord %in% "TB1 indeterminate & TB2 positive"),
    
    .groups = "drop"
    
  ) |> 
  
  # summaries
  
  mutate(
    
    tube_result_concord_both_pos_prop = round(100 * tube_result_concord_both_pos / n, 1),
    tube_result_concord_both_pos_prop2 = paste0(sprintf("%.1f", tube_result_concord_both_pos_prop), "%", " ", "(", tube_result_concord_both_pos, "/", n, ")"),
    
    tube_result_concord_both_neg_prop = round(100 * tube_result_concord_both_neg / n, 1),
    tube_result_concord_both_neg_prop2 = paste0(sprintf("%.1f", tube_result_concord_both_neg_prop), "%", " ", "(", tube_result_concord_both_neg, "/", n, ")"),
    
    tube_result_concord_both_ind_prop = round(100 * tube_result_concord_both_ind / n, 1),
    tube_result_concord_both_ind_prop2 = paste0(sprintf("%.1f", tube_result_concord_both_ind_prop), "%", " ", "(", tube_result_concord_both_ind, "/", n, ")"),
    
    tube_result_concord_tb1ptb2n_prop = round(100 * tube_result_concord_tb1ptb2n / n, 1),
    tube_result_concord_tb1ptb2n_prop2 = paste0(sprintf("%.1f", tube_result_concord_tb1ptb2n_prop), "%", " ", "(", tube_result_concord_tb1ptb2n, "/", n, ")"),
    
    tube_result_concord_tb1ntb2p_prop = round(100 * tube_result_concord_tb1ntb2p / n, 1),
    tube_result_concord_tb1ntb2p_prop2 = paste0(sprintf("%.1f", tube_result_concord_tb1ntb2p_prop), "%", " ", "(", tube_result_concord_tb1ntb2p, "/", n, ")"),
    
    tube_result_concord_tb1ptb2i_prop = round(100 * tube_result_concord_tb1ptb2i / n, 1),
    tube_result_concord_tb1ptb2i_prop2 = paste0(sprintf("%.1f", tube_result_concord_tb1ptb2i_prop), "%", " ", "(", tube_result_concord_tb1ptb2i, "/", n, ")"),
    
    tube_result_concord_tb1itb2p_prop = round(100 * tube_result_concord_tb1itb2p / n, 1),
    tube_result_concord_tb1itb2p_prop2 = paste0(sprintf("%.1f", tube_result_concord_tb1itb2p_prop), "%", " ", "(", tube_result_concord_tb1itb2p, "/", n, ")"),
    
  )

tb2_QFT_overall <- tb2_QFT_overall |>  # add age group column for binding rows
  mutate(age_10y = "Overall", .before = 1)


## age group 

tb2_QFT_age <- df |> 
  group_by(age_10y) |>
  summarise(
    
    # n
    n = n(),
    
    tube_result_concord_both_pos = sum(tube_result_concord %in% "Both TB1 & TB2 positive"),
    tube_result_concord_both_neg = sum(tube_result_concord %in% "Both TB1 & TB2 negative"),
    tube_result_concord_both_ind = sum(tube_result_concord %in% "Both TB1 & TB2 indeterminate"),
    
    tube_result_concord_tb1ptb2n = sum(tube_result_concord %in% "TB1 positive & TB2 negative"),
    tube_result_concord_tb1ptb2i = sum(tube_result_concord %in% "TB1 positive & TB2 indeterminate"),
    
    tube_result_concord_tb1ntb2p = sum(tube_result_concord %in% "TB1 negative & TB2 positive"),
    tube_result_concord_tb1ntb2i = sum(tube_result_concord %in% "TB1 negative & TB2 indeterminate"),
    
    tube_result_concord_tb1itb2p = sum(tube_result_concord %in% "TB1 indeterminate & TB2 positive"),
    tube_result_concord_tb1itb2n = sum(tube_result_concord %in% "TB1 indeterminate & TB2 negative"),
    
    .groups = "drop"
    
  ) |> 
  
  # summaries
  
  mutate( 
    
    # QFT-Plus
    
    tube_result_concord_both_pos_prop = round(100 * tube_result_concord_both_pos / n, 1),
    tube_result_concord_both_pos_prop2 = paste0(sprintf("%.1f", tube_result_concord_both_pos_prop), "%", " ", "(", tube_result_concord_both_pos, "/", n, ")"),
    
    tube_result_concord_both_neg_prop = round(100 * tube_result_concord_both_neg / n, 1),
    tube_result_concord_both_neg_prop2 = paste0(sprintf("%.1f", tube_result_concord_both_neg_prop), "%", " ", "(", tube_result_concord_both_neg, "/", n, ")"),
    
    tube_result_concord_both_ind_prop = round(100 * tube_result_concord_both_ind / n, 1),
    tube_result_concord_both_ind_prop2 = paste0(sprintf("%.1f", tube_result_concord_both_ind_prop), "%", " ", "(", tube_result_concord_both_ind, "/", n, ")"),
    
    tube_result_concord_tb1ptb2n_prop = round(100 * tube_result_concord_tb1ptb2n / n, 1),
    tube_result_concord_tb1ptb2n_prop2 = paste0(sprintf("%.1f", tube_result_concord_tb1ptb2n_prop), "%", " ", "(", tube_result_concord_tb1ptb2n, "/", n, ")"),
    
    tube_result_concord_tb1ptb2i_prop = round(100 * tube_result_concord_tb1ptb2i / n, 1),
    tube_result_concord_tb1ptb2i_prop2 = paste0(sprintf("%.1f", tube_result_concord_tb1ptb2i_prop), "%", " ", "(", tube_result_concord_tb1ptb2i, "/", n, ")"),
    
    tube_result_concord_tb1ntb2p_prop = round(100 * tube_result_concord_tb1ntb2p / n, 1),
    tube_result_concord_tb1ntb2p_prop2 = paste0(sprintf("%.1f", tube_result_concord_tb1ntb2p_prop), "%", " ", "(", tube_result_concord_tb1ntb2p, "/", n, ")"),
    
    tube_result_concord_tb1ntb2i_prop = round(100 * tube_result_concord_tb1ntb2i / n, 1),
    tube_result_concord_tb1ntb2i_prop2 = paste0(sprintf("%.1f", tube_result_concord_tb1ntb2i_prop), "%", " ", "(", tube_result_concord_tb1ntb2i, "/", n, ")"),
    
    tube_result_concord_tb1itb2p_prop = round(100 * tube_result_concord_tb1itb2p / n, 1),
    tube_result_concord_tb1itb2p_prop2 = paste0(sprintf("%.1f", tube_result_concord_tb1itb2p_prop), "%", " ", "(", tube_result_concord_tb1itb2p, "/", n, ")"),
    
    tube_result_concord_tb1itb2n_prop = round(100 * tube_result_concord_tb1itb2n / n, 1),
    tube_result_concord_tb1itb2n_prop2 = paste0(sprintf("%.1f", tube_result_concord_tb1itb2n_prop), "%", " ", "(", tube_result_concord_tb1itb2n, "/", n, ")")
    
  )


tb2_QFT_age


# sex

tb2_QFT_sex <- df |>
  mutate(sex = fct_rev(sex)) |> 
  group_by(sex) |>
  summarise(
    # n
    n = n(),
    # tube-specific QFT-Plus result
    tube_result_concord_both_pos = sum(tube_result_concord %in% "Both TB1 & TB2 positive"),
    tube_result_concord_both_neg = sum(tube_result_concord %in% "Both TB1 & TB2 negative"),
    tube_result_concord_both_ind = sum(tube_result_concord %in% "Both TB1 & TB2 indeterminate"),
    
    tube_result_concord_tb1ptb2n = sum(tube_result_concord %in% "TB1 positive & TB2 negative"),
    tube_result_concord_tb1ptb2i = sum(tube_result_concord %in% "TB1 positive & TB2 indeterminate"),
    
    tube_result_concord_tb1ntb2p = sum(tube_result_concord %in% "TB1 negative & TB2 positive"),
    tube_result_concord_tb1ntb2i = sum(tube_result_concord %in% "TB1 negative & TB2 indeterminate"),
    
    tube_result_concord_tb1itb2p = sum(tube_result_concord %in% "TB1 indeterminate & TB2 positive"),
    tube_result_concord_tb1itb2n = sum(tube_result_concord %in% "TB1 indeterminate & TB2 negative"),
    
    .groups = "drop"
    
  ) |> 
  
  # summaries
  
  mutate( 
    # QFT-Plus
    tube_result_concord_both_pos_prop = round(100 * tube_result_concord_both_pos / n, 1),
    tube_result_concord_both_pos_prop2 = paste0(sprintf("%.1f", tube_result_concord_both_pos_prop), "%", " ", "(", tube_result_concord_both_pos, "/", n, ")"),
    
    tube_result_concord_both_neg_prop = round(100 * tube_result_concord_both_neg / n, 1),
    tube_result_concord_both_neg_prop2 = paste0(sprintf("%.1f", tube_result_concord_both_neg_prop), "%", " ", "(", tube_result_concord_both_neg, "/", n, ")"),
    
    tube_result_concord_tb1ptb2n_prop = round(100 * tube_result_concord_tb1ptb2n / n, 1),
    tube_result_concord_tb1ptb2n_prop2 = paste0(sprintf("%.1f", tube_result_concord_tb1ptb2n_prop), "%", " ", "(", tube_result_concord_tb1ptb2n, "/", n, ")"),
    
    tube_result_concord_tb1ntb2p_prop = round(100 * tube_result_concord_tb1ntb2p / n, 1),
    tube_result_concord_tb1ntb2p_prop2 = paste0(sprintf("%.1f", tube_result_concord_tb1ntb2p_prop), "%", " ", "(", tube_result_concord_tb1ntb2p, "/", n, ")"),
    
    tube_result_concord_both_ind_prop = round(100 * tube_result_concord_both_ind / n, 1),
    tube_result_concord_both_ind_prop2 = paste0(sprintf("%.1f", tube_result_concord_both_ind_prop), "%", " ", "(", tube_result_concord_both_ind, "/", n, ")"),
    
    tube_result_concord_tb1ptb2i_prop = round(100 * tube_result_concord_tb1ptb2i / n, 1),
    tube_result_concord_tb1ptb2i_prop2 = paste0(sprintf("%.1f", tube_result_concord_tb1ptb2i_prop), "%", " ", "(", tube_result_concord_tb1ptb2i, "/", n, ")"),
    
    tube_result_concord_tb1itb2p_prop = round(100 * tube_result_concord_tb1itb2p / n, 1),
    tube_result_concord_tb1itb2p_prop2 = paste0(sprintf("%.1f", tube_result_concord_tb1itb2p_prop), "%", " ", "(", tube_result_concord_tb1itb2p, "/", n, ")"),
    
  )


tb2_QFT_sex


## combine covariate-specific outcome summaries

tb2_outcomes_covs_summ <- bind_rows(
  tb2_QFT_overall |> 
    mutate(characteristic = "overall",
           levels = "all"),
  tb2_QFT_age |> 
    mutate(characteristic = "age", .before = age_10y) |> 
    rename(levels = age_10y),
  tb2_QFT_sex |> 
    mutate(characteristic = "sex", .before = sex) |> 
    rename(levels = sex),
  
)

tbl4 <- tb2_outcomes_covs_summ |> 
  select(levels,
         tube_result_concord_both_pos_prop2,
         tube_result_concord_both_neg_prop2,
         tube_result_concord_both_ind_prop2,
         tube_result_concord_tb1ptb2n_prop2,
         tube_result_concord_tb1ptb2i_prop2,
         tube_result_concord_tb1ntb2p_prop2,
         tube_result_concord_tb1ntb2i_prop2,
         tube_result_concord_tb1itb2p_prop2,
         tube_result_concord_tb1itb2n_prop2
  ) |> 
  
  mutate(Characteristic = case_when(levels == "all" ~ "Overall",
                                    levels == "10-19" ~ "Age: 10-19",
                                    levels == "20-29" ~ "Age: 20-29",
                                    levels == "30-40" ~ "Age: 30-40",
                                    levels == "Male" ~ "Sex: Male",
                                    levels == "Female" ~ "Sex: Female"),
         .before = levels) |> 
  
  select(Characteristic,
         `TB1 & TB2 positive, % (n/N)` = tube_result_concord_both_pos_prop2,
         `TB1 & TB2 negative, % (n/N)` = tube_result_concord_both_neg_prop2,
         `TB1 positive & TB2 negative, % (n/N)` = tube_result_concord_tb1ptb2n_prop2,
         `TB1 negative & TB2 positive, % (n/N)` = tube_result_concord_tb1ntb2p_prop2,
         `TB1 & TB2 indeterminate, % (n/N)` = tube_result_concord_both_ind_prop2,
         `TB1 positive & TB2 indeterminate, % (n/N)` = tube_result_concord_tb1ptb2i_prop2,
         ` TB1 indeterminate & TB2 positive,  % (n/N)` = tube_result_concord_tb1itb2p_prop2
  )


tbl4


# Thresholds -------------------------------------------------------------------

# summarise QFT-Plus results across thresholds

thresh_results <- data_thresh |> 
  group_by(threshold, QFT_reclass) |>
  summarise(n = n()) |>
  arrange(threshold, QFT_reclass)

thresh_results

thresh_results_age_sex <- data_thresh |> 
  group_by(threshold, age, sex, QFT_reclass) |>
  summarise(n = n()) |>
  arrange(threshold, age, sex, QFT_reclass) |>
  ungroup() |> 
  filter(threshold %in% c(0.35)) |>
  group_by(age, sex, QFT_reclass) |> 
  mutate(n = sum(n))

thresh_results_age_sex


# Prevalence, ARTI & M:F ratios

# age & sex predictions

df_pred_list <- map_dfr(names(models), function(thresh) {
  model <- models[[thresh]]
  
  # data range for threshold
  data_thresh <- data_thresh |> filter(threshold == as.numeric(thresh))
  
  df_pred <- expand_grid(
    age = unique(data_thresh$age),
    sex = unique(data_thresh$sex)
  ) |> 
    mutate(threshold = as.numeric(thresh))
  
  # posterior predictions
  set.seed(SEED)
  add_epred_draws(model, newdata = df_pred, ndraws = NDRAWS) 
  
})

# predict older age for ARTI calculation
df_pred_list_older <- map_dfr(names(models), function(thresh) {
  model <- models[[thresh]]
  
  # data range for threshold
  data_thresh <- data_thresh |> filter(threshold == as.numeric(thresh))
  
  df_pred <- expand_grid(
    age = unique(data_thresh$age) + 1,
    sex = unique(data_thresh$sex)
  ) |> 
    mutate(threshold = as.numeric(thresh))
  
  # posterior predictions with full draws
  set.seed(SEED)
  add_epred_draws(model, newdata = df_pred, ndraws = NDRAWS) |>
    
    mutate(age=age-1) |>
    
    rename(.epred_older = .epred)
})

# join

df_pred_list <- df_pred_list |>
  select(.row, .draw, age, sex, threshold, .epred) |> 
  left_join(df_pred_list_older |> 
              select(.row, .draw, age, sex, threshold, .epred_older),
            by = c(".row", ".draw", "age", "sex", "threshold"))


# select sample predictions
df_pred_list |>
  ungroup() |>
  filter(age %in% c(10, 20, 30, 40, 50), threshold %in% c(0.1, 0.2, 0.3, 0.35, 0.4, 0.5)) |>
  group_by(threshold, sex, age) |>
  mean_qi(.epred) |>
  mutate(across(c(.epred, .lower, .upper), ~percent(., accuracy = 0.1)))

# calculate ARTI & M:F ratios

## ARTI
df_pred_list_arti <- df_pred_list |> 
  mutate(arti = (.epred_older - .epred) / (1 - .epred)) |> 
  select(-.epred_older, arti)

# M:F ratios
df_pred_list_mfr <- df_pred_list |> 
  select(-.epred_older) |>
  pivot_wider(id_cols = c(age, threshold, .draw), names_from = sex, values_from = .epred) |> 
  rename(m = Male, f = Female) |> 
  mutate(mfr = m / f)

# plot

# prevalence

f4 <- df_pred_list |> 
  mutate(sex = fct_rev(sex)) |> 
  filter(!threshold > 0.5) |>
  group_by(threshold, age, sex) |> 
  mean_qi(.epred, .width = 0.95) |> 
  ggplot(aes(x = age, y = .epred, color = sex)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = sex), alpha = 0.3, linewidth = 0) +
  facet_wrap(~threshold, ncol = 2,
             labeller = labeller(threshold = function(x) paste("Threshold =", x, "IU/ml"))) +
  labs(x = "Age (years)",
       y = "Predicted Probability of Positive QFT-Plus (%)") +
  
  scale_color_brewer(name = "Sex", type = "qual", palette = "Set1") +
  scale_fill_brewer(name = "Sex", type = "qual", palette = "Set1") +
  
  scale_x_continuous(breaks = seq(10, 40, by = 5)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  
  theme_bw(base_size = 11)+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        plot.subtitle = element_text(size = 11))

f4

# ARTI

f5 <- df_pred_list_arti |> 
  mutate(sex = fct_rev(sex)) |> 
  filter(!threshold > 0.5) |>
  group_by(threshold, age, sex) |> 
  mean_qi(arti, .width = 0.95) |> 
  ggplot(aes(x = age, y = arti, color = sex)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = sex), alpha = 0.3, linewidth = 0) +
  facet_wrap(~threshold, ncol = 2, labeller = labeller(threshold = function(x) paste("Threshold =", x, "IU/ml"))) +
  labs(x = "Age (years)",
       y = "Risk of conversion (%)") +
  
  scale_color_brewer(name = "Sex", type = "qual", palette = "Set1") +
  scale_fill_brewer(name = "Sex", type = "qual", palette = "Set1") +
  
  scale_x_continuous(breaks = seq(10, 40, by = 5)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  
  theme_bw(base_size = 11)+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        plot.subtitle = element_text(size = 11))

f5

# M:F ratio

f6 <- df_pred_list_mfr |> 
  filter(!threshold > 0.5) |>
  group_by(threshold, age) |> 
  mean_qi(mfr, .width = 0.95) |> 
  ggplot(aes(x = age, y = mfr, color = as.factor(threshold))) +
  geom_line(linewidth = 0.6) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.15, linewidth = 0) +
  facet_wrap(~threshold, ncol = 2, labeller = labeller(threshold = function(x) paste("Threshold =", x, "IU/ml"))) +
  labs(x = "Age (years)",
       y = "Male-to-Female prevalence ratio",
       caption = "IU/mL: International units per millilitre")  +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 11),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

f6



# SUPPLEMENTARY material ----------------------------------------------------------------

# Figure 1: TB1 & TB2 distributions

supp_f1 <- df |> 
  select(age_10y, sex, tb1_nil, tb2_nil) |> 
  pivot_longer(cols = c(tb1_nil, tb2_nil), names_to = "tube", values_to = "value") |> 
  ggplot() +
  geom_histogram(aes(x = value, fill = sex), binwidth = 0.4, alpha = 0.6) +
  facet_grid(~ifelse(tube=="tb1_nil", "TB1-Nil","TB2-Nil")) +
  
  labs(x = "IFN-γ (ΙU/mL",
       y = "Count",
       title = "Distribution of QuantiFERON-TB Gold Plus TB1 & TB2 values by age & sex") +
  
  scale_y_break(c(75,2250)) +
  
  scale_fill_brewer(name = "Sex", type = "qual", palette = "Set1") +
  
  theme_bw() +
  theme(legend.position = "bottom")

supp_f1 

# Supp Table 1 TB1 & TB2 distributions across characteristics --------------------------

tb1tb2_dist_summ <- df |>
  summarise(
    n = n(),
    `TB1<0` = paste0(sprintf("%.1f", round(100 * sum(tb1_nil < 0) / n(), 1)), "%", " (", sum(tb1_nil < 0), "/", n(), ")"),
    `TB2<0` = paste0(sprintf("%.1f", round(100 * sum(tb2_nil < 0) / n(), 1)), "%", " (", sum(tb2_nil < 0), "/", n(), ")"),
    
    `TB1==0` = paste0(sprintf("%.1f", round(100 * sum(tb1_nil == 0) / n(), 1)), "%", " (", sum(tb1_nil == 0), "/", n(), ")"),
    `TB2==0` = paste0(sprintf("%.1f", round(100 * sum(tb2_nil == 0) / n(), 1)), "%", " (", sum(tb2_nil == 0), "/", n(), ")"),
    
    `TB1<0.35` = paste0(sprintf("%.1f", round(100 * sum(tb1_nil < 0.35) / n(), 1)), "%", " (", sum(tb1_nil < 0.35), "/", n(), ")"),    
    `TB2<0.35` = paste0(sprintf("%.1f", round(100 * sum(tb2_nil < 0.35) / n(), 1)), "%", " (", sum(tb2_nil < 0.35), "/", n(), ")"),
    
    `TB1≥0.35` = paste0(sprintf("%.1f", round(100 * sum(tb1_nil >= 0.35) / n(), 1)), "%", " (", sum(tb1_nil >= 0.35), "/", n(), ")"),    
    `TB2≥0.35` = paste0(sprintf("%.1f", round(100 * sum(tb2_nil >= 0.35) / n(), 1)), "%", " (", sum(tb2_nil >= 0.35), "/", n(), ")"),
    
    `TB1>10` = paste0(sprintf("%.1f", round(100 * sum(tb1_raw == "> 10¶" | tb1_raw == "> 10#") / n(), 1)), "%", " (", sum(tb1_raw == "> 10¶" | tb1_raw == "> 10#"), "/", n(), ")"),
    `TB2>10` = paste0(sprintf("%.1f", round(100 * sum(tb2_raw == "> 10¶" | tb2_raw == "> 10#") / n(), 1)), "%", " (", sum(tb2_raw == "> 10¶" | tb2_raw == "> 10#"), "/", n(), ")"),
    
    `TB1<0.2` = paste0(sprintf("%.1f", round(100 * sum(tb1_nil < 0.2) / n(), 1)), "%", " (", sum(tb1_nil < 0.2), "/", n(), ")"),
    `TB2<0.2` = paste0(sprintf("%.1f", round(100 * sum(tb2_nil < 0.2) / n(), 1)), "%", " (", sum(tb2_nil < 0.2), "/", n(), ")"),
    
    `TB1 0.2-0.7` = paste0(sprintf("%.1f", round(100 * sum(tb1_nil >= 0.2 & tb1_nil <= 0.7) / n(), 1)), "%", " (", sum(tb1_nil >= 0.2 & tb1_nil <= 0.7), "/", n(), ")"),
    `TB2 0.2-0.7` = paste0(sprintf("%.1f", round(100 * sum(tb2_nil >= 0.2 & tb2_nil <= 0.7) / n(), 1)), "%", " (", sum(tb2_nil >= 0.2 & tb2_nil <= 0.7), "/", n(), ")"),
    
    `TB1>0.7` = paste0(sprintf("%.1f", round(100 * sum(tb1_nil > 0.7) / n(), 1)), "%", " (", sum(tb1_nil > 0.7), "/", n(), ")"),
    `TB2>0.7` = paste0(sprintf("%.1f", round(100 * sum(tb2_nil > 0.7) / n(), 1)), "%", " (", sum(tb2_nil > 0.7), "/", n(), ")")
  ) |> 
  
  pivot_longer(cols = -n, names_to = "Summary", values_to = "% (n/N)") 

tb1tb2_dist_summ


# Supp Table 2: Correlation TB1 & TB2 --------------------------------------------------------

corr_tb1tb2_overall <- summarise_corr_fn(m0_tb1tb2_corr_overall)
corr_tb1tb2_1019 <- summarise_corr_fn(m0_tb1tb2_corr_1019)
corr_tb1tb2_2029 <- summarise_corr_fn(m0_tb1tb2_corr_2029)
corr_tb1tb2_3040 <- summarise_corr_fn(m0_tb1tb2_corr_3040)

corr_tb1tb2_sexMale <- summarise_corr_fn(m0_tb1tb2_corr_sexMale)
corr_tb1tb2_sexFemale <- summarise_corr_fn(m0_tb1tb2_corr_sexFemale)

# combine
corr_tb1tb2 <- bind_rows(corr_tb1tb2_overall, corr_tb1tb2_1019, corr_tb1tb2_2029,
                         corr_tb1tb2_3040, corr_tb1tb2_sexMale, corr_tb1tb2_sexFemale) |> 
  
  select(Group, `Correlation (95% CrI)`)

corr_tb1tb2


## Supp Figure 2: PPC

# PPC

# observed data

d_pp <- d |>
  transmute(
    age = age,
    sex = sex,
    tb1_nil_qia = tb1_nil_qia,
    tb2_nil_qia = tb2_nil_qia
  )

# prop zero
obs_zero <- d_pp |>
  rename("TB1" = tb1_nil_qia,
         "TB2" = tb2_nil_qia)|>
  pivot_longer(cols = c("TB1", "TB2"),
               names_to = "tube")|>
  select(tube, age, sex, value) |> 
  mutate(is_zero = (value == 0)) |> 
  group_by(tube, age, sex)|>
  summarise(
    obs_prop = mean(is_zero),
    n = n(),
    .groups = "drop"
  )

# cont
obs_cont <- d_pp |> 
  filter(tb1_nil_qia > 0) |> 
  group_by(age, sex) |> 
  summarise(mean_tb1 = mean(tb1_nil_qia),
            med_tb1 = median(tb1_nil_qia))


# Hurdle

nd <- expand.grid(age = 10:40,
                  sex = factor(c("Male", "Female"), levels = c("Male", "Female")))

levels(nd$sex)

# 1) posterior

# hu

# TB1
hu_draws_tb1 <- m0_tb1_nil_qia_hurdle_age_sex_sp |> 
  add_linpred_draws(
    newdata = nd,
    dpar = "hu",
    ndraws = 1000
  ) |>
  mutate(p0 = plogis(.linpred)) |> 
  mutate(tube = "TB1")

# TB2
hu_draws_tb2 <- m0_tb2_nil_qia_hurdle_age_sex_sp |> 
  add_linpred_draws(
    newdata = nd,
    dpar = "hu",
    ndraws = 1000
  ) |>
  mutate(p0 = plogis(.linpred)) |> 
  mutate(tube = "TB2")

# join
hu_draws <- bind_rows(hu_draws_tb1, hu_draws_tb2) |> 
  mutate(tube = factor(tube, levels = c("TB1", "TB2")))


# 1) posterior predicted TB1

# predictions

# TB1
yrep_draws_tb1 <- predicted_draws(
  newdata = d_pp,
  object = m0_tb1_nil_qia_hurdle_age_sex_sp,
  ndraws = 1000
)|>
  as.data.frame() |>
  select(.draw, sex, age, .prediction) |> 
  mutate(tube = "TB1")

# TB2
yrep_draws_tb2 <- predicted_draws(
  newdata = d_pp,
  object = m0_tb2_nil_qia_hurdle_age_sex_sp,
  ndraws = 1000
)|>
  as.data.frame() |>
  select(.draw, sex, age, .prediction) |> 
  mutate(tube = "TB2")

# joint
yrep_draws <- bind_rows(yrep_draws_tb1, yrep_draws_tb2) |> 
  mutate(is_zero = (.prediction == 0))|>
  mutate(tube = factor(tube, levels = c("TB1", "TB2"))) |> 
  mutate(.prediction_trunc = case_when(.prediction > 10 ~ 10, TRUE ~ .prediction), .after = .prediction)


# posterior predicted zeros
yrep_prop_zero <- yrep_draws|>
  group_by(tube, .draw, sex, age)|>
  summarise(mean_is_zero = mean(is_zero), .groups = "drop")|>
  group_by(tube, sex, age)|>
  summarise(
    pred_mean = mean(mean_is_zero),
    pred_l95  = quantile(mean_is_zero, 0.025),
    pred_u95  = quantile(mean_is_zero, 0.975),
    .groups = "drop"
  )


# posterior predicted > 0 positives
yrep_pos <- yrep_draws|>
  mutate(.prediction = case_when(.prediction > 10 ~ 10, TRUE ~ .prediction)) |>
  filter(.prediction >= 0.35) |> 
  group_by(.draw, sex, age)|>
  summarise(mean_tb1 = mean(.prediction), .groups = "drop")|>
  group_by(sex, age)|>
  summarise(
    pred_mean = mean(mean_tb1),
    pred_l95  = quantile(mean_tb1, 0.025),
    pred_u95  = quantile(mean_tb1, 0.975),
    .groups = "drop"
  )


# 4) plot

# obs vs predicted zeros (same panel)

supp_f2 <- obs_zero |>
  mutate(sex = fct_rev(sex)) |> 
  left_join(yrep_prop_zero, by = c("tube", "age", "sex"))|>
  ggplot(aes(x = age, y = 1-obs_prop, color = sex)) +
  # observed
  geom_point(aes(shape = sex, fill = sex), size = 1.5, position = position_dodge(0.6)
  ) +
  # predictive intervals
  geom_errorbar(aes(ymin = 1-pred_l95, ymax = 1-pred_u95, linetype = sex),
                width = 0.4, position = position_dodge(0.6)
  ) +
  
  # predictive mean (open)
  geom_point(aes(y = 1-pred_mean, shape = sex), size = 1.8, fill = NA,
             position = position_dodge(0.6)
  ) +
  
  coord_cartesian(ylim = c(0, .75)) +
  scale_y_continuous(
    breaks = c(0, 0.75, seq(0, 1, by = 0.25)),
    labels = scales::percent_format(accuracy = 1)
  ) +
  
  scale_color_brewer(name = "Sex", palette = "Set1") +
  scale_fill_brewer(name = "Sex", palette = "Set1") +
  scale_shape_manual(name = "Sex", values = c("Male" = 24, "Female" = 21)) +
  scale_linetype_manual(name = "Sex", values = c("Male" = "solid", "Female" = "dashed")) +
  
  facet_wrap(~ tube) +
  
  labs(
    x = "Age (years)",
    y = "Proportion with TB1 ≥ 0.35 IU/mL (%)",
    caption = "Symbols: filled are observed values; open are posterior predictive means. \nBars are 95% credible intervals. Solid are for females, dashed are for males."
  ) +
  
  guides(
    fill = "none",  # suppress duplicate
    linetype = "none",  # suppress duplicate
    color = guide_legend(
      override.aes = list(
        shape = c(24, 21),
        fill = c("white", "white"),
        linetype = "solid"))
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 11),
    text = element_text(size = 11),
    legend.position = "bottom"
  )

supp_f2


## Supp Figure 3: TB2-TB1 differential continuous distributions -------------

# histogram
supp_f3 <- df |> 
  mutate(sex = fct_rev(sex)) |> 
  select(age_10y, sex, tb2tb1_diff) |> 
  ggplot() +
  geom_histogram(aes(x = tb2tb1_diff), fill = "black", binwidth = 0.4, alpha = 1) +
  labs(x = "TB2-TB1 differential (IFN-γ ΙU/mL)",
       y = "Count",
       # title = "Distribution of QuantiFERON-TB Gold Plus Nil-subtracted TB1 & TB2 values by age & sex",
       # caption = "The y-axis is discontinuous to improve visualisation outside the large mass at zero") +
  )+
  # scale_y_break(c(75,2400)) +
  
  scale_fill_brewer(name = "Sex", type = "qual", palette = "Set1") +
  theme_bw() +
  theme(legend.position = "bottom")

supp_f3 

