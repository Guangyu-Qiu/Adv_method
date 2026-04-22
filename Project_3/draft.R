################################################################################
# Framingham Heart Study - Stroke Risk Analysis
# Complete Analysis Script
################################################################################

library(readr)
library(ggplot2)
library(tidyverse)
library(gtsummary)
library(tibble)
library(survival)
library(survminer)
library(gridExtra)

################################################################################
# SECTION 1: DATA CLEANING
################################################################################

data <- read_csv("frmgham2.csv")

# Verify baseline structure
baseline_check <- data %>%
  group_by(RANDID) %>%
  summarise(has_baseline = any(TIME == 0)) %>%
  ungroup()
table(baseline_check$has_baseline)

baseline_count <- data %>%
  filter(TIME == 0) %>%
  count(RANDID)
table(baseline_count$n)

# Generate baseline analysis dataset
analysis_baseline <- data %>%
  arrange(RANDID, TIME) %>%
  group_by(RANDID) %>%
  filter(TIME == 0) %>%
  ungroup() %>%
  transmute(
    id       = RANDID,
    sex      = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")),
    age      = AGE,
    chol     = TOTCHOL,
    bp       = SYSBP,
    smoke    = factor(CURSMOKE, levels = c(0, 1), labels = c("No", "Yes")),
    bmi      = BMI,           # continuous only (per revision)
    diabetes = factor(DIABETES, levels = c(0, 1), labels = c("No", "Yes")),
    hd       = factor(PREVCHD, levels = c(0, 1), labels = c("No", "Yes")),
    med      = factor(BPMEDS,  levels = c(0, 1), labels = c("No", "Yes")),
    # 10-year censored outcome
    time_raw   = TIMESTRK,
    stroke_raw = STROKE,
    time   = pmin(TIMESTRK, 3652),
    stroke = ifelse(TIMESTRK > 3652, 0, STROKE)
  ) %>%
  # IMPORTANT: exclude individuals with time == 0 (stroke at baseline)
  filter(time_raw > 0)

# Missingness check
missing_summary <- analysis_baseline %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "N_Missing") %>%
  filter(N_Missing > 0)
print(missing_summary)

# Range checks
summary(analysis_baseline %>% dplyr::select(age, chol, bp, bmi))

################################################################################
# SECTION 2: TABLE 1 — BASELINE CHARACTERISTICS
################################################################################

tbl1 <- analysis_baseline %>%
  dplyr::select(sex, age, chol, bp, smoke, bmi, diabetes, hd, med) %>%
  tbl_summary(
    by = sex,
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      age      = "Age (years)",
      chol     = "Total Cholesterol (mg/dL)",
      bp       = "Systolic Blood Pressure (mmHg)",
      smoke    = "Current Smoker",
      bmi      = "BMI (kg/m²)",
      diabetes = "Diabetes",
      hd       = "Prevalent Coronary Heart Disease",
      med      = "Anti-hypertensive Medications"
    ),
    missing = "ifany",
    digits  = all_continuous() ~ 2
  ) %>%
  add_overall(last = TRUE) %>%
  bold_labels()

tbl1

################################################################################
# SECTION 3: DESCRIPTIVE — RISK FACTOR CHANGES OVER 3 PERIODS
# (uses all periods, not just baseline)
################################################################################

# BP and diabetes change over periods
longitudinal_summary <- data %>%
  filter(PERIOD %in% c(1, 2, 3)) %>%
  mutate(Period = factor(PERIOD, labels = c("Period 1 (~1956)", "Period 2 (~1962)", "Period 3 (~1968)"))) %>%
  group_by(Period) %>%
  summarise(
    N               = n(),
    BP_mean         = round(mean(SYSBP, na.rm = TRUE), 2),
    BP_sd           = round(sd(SYSBP, na.rm = TRUE), 2),
    Diabetes_pct    = round(mean(DIABETES, na.rm = TRUE) * 100, 1),
    Cholesterol_mean = round(mean(TOTCHOL, na.rm = TRUE), 2),
    Cholesterol_sd  = round(sd(TOTCHOL, na.rm = TRUE), 2),
    BMI_mean        = round(mean(BMI, na.rm = TRUE), 2),
    Smoke_pct       = round(mean(CURSMOKE, na.rm = TRUE) * 100, 1)
  )

print(longitudinal_summary)

# Nicer gtsummary version of longitudinal table
tbl_longitudinal <- data %>%
  filter(PERIOD %in% c(1, 2, 3)) %>%
  mutate(
    Period   = factor(PERIOD, labels = c("Period 1", "Period 2", "Period 3")),
    diabetes = factor(DIABETES, levels = c(0,1), labels = c("No","Yes")),
    smoke    = factor(CURSMOKE,  levels = c(0,1), labels = c("No","Yes"))
  ) %>%
  dplyr::select(Period, SYSBP, TOTCHOL, BMI, diabetes, smoke) %>%
  tbl_summary(
    by = Period,
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      SYSBP    = "Systolic BP (mmHg)",
      TOTCHOL  = "Total Cholesterol (mg/dL)",
      BMI      = "BMI (kg/m²)",
      diabetes = "Diabetes",
      smoke    = "Current Smoker"
    ),
    missing = "ifany",
    digits  = all_continuous() ~ 2
  ) %>%
  bold_labels()

tbl_longitudinal

################################################################################
# SECTION 4: EXPLORATORY — KM CURVES + LOG-RANK TESTS
################################################################################

# Helper: plot KM by a categorical variable, stratified by sex
km_plot <- function(data, var, var_label, filename) {
  
  surv_obj <- Surv(data$time, data$stroke)
  
  # Combined (both sexes)
  formula_combined <- as.formula(paste0("Surv(time, stroke) ~ ", var))
  fit_combined <- survfit(formula_combined, data = data)
  
  p_combined <- ggsurvplot(
    fit_combined,
    data       = data,
    pval       = TRUE,
    pval.method = TRUE,
    conf.int   = TRUE,
    risk.table = TRUE,
    title      = paste0("KM Curve by ", var_label, " (All)"),
    xlab       = "Time (days)",
    ylab       = "Stroke-free Survival Probability",
    legend.title = var_label,
    palette    = "jco",
    ggtheme    = theme_bw(base_size = 13)
  )
  
  # Stratified by sex
  formula_sex <- as.formula(paste0("Surv(time, stroke) ~ ", var))
  
  male_data   <- filter(data, sex == "Male")
  female_data <- filter(data, sex == "Female")
  
  fit_male   <- survfit(formula_sex, data = male_data)
  fit_female <- survfit(formula_sex, data = female_data)
  
  p_male <- ggsurvplot(
    fit_male,
    data       = male_data,
    pval       = TRUE,
    pval.method = TRUE,
    conf.int   = TRUE,
    risk.table = TRUE,
    title      = paste0("KM Curve by ", var_label, " (Male)"),
    xlab       = "Time (days)",
    ylab       = "Stroke-free Survival Probability",
    legend.title = var_label,
    palette    = "jco",
    ggtheme    = theme_bw(base_size = 13)
  )
  
  p_female <- ggsurvplot(
    fit_female,
    data       = female_data,
    pval       = TRUE,
    pval.method = TRUE,
    conf.int   = TRUE,
    risk.table = TRUE,
    title      = paste0("KM Curve by ", var_label, " (Female)"),
    xlab       = "Time (days)",
    ylab       = "Stroke-free Survival Probability",
    legend.title = var_label,
    palette    = "jco",
    ggtheme    = theme_bw(base_size = 13)
  )
  
  # Save
  ggsave(paste0(filename, "_all.png"),
         plot   = print(p_combined),
         width  = 10, height = 8, dpi = 300, bg = "white")
  ggsave(paste0(filename, "_male.png"),
         plot   = print(p_male),
         width  = 10, height = 8, dpi = 300, bg = "white")
  ggsave(paste0(filename, "_female.png"),
         plot   = print(p_female),
         width  = 10, height = 8, dpi = 300, bg = "white")
  
  invisible(list(combined = p_combined, male = p_male, female = p_female))
}

# KM for each categorical variable
km_plot(analysis_baseline, "smoke",    "Smoking Status",              "km_smoke")
km_plot(analysis_baseline, "diabetes", "Diabetes",                    "km_diabetes")
km_plot(analysis_baseline, "hd",       "Prevalent CHD",               "km_hd")
km_plot(analysis_baseline, "med",      "Anti-hypertensive Medication","km_med")

################################################################################
# SECTION 5: VARIABLE SELECTION VIA BACKWARD AIC
# Strategy:
#   1. Univariable Cox for each of the 8 candidate variables
#   2. Include those with p < 0.20 in multivariable model
#   3. Backward selection by AIC (step())
#   4. Force-add the 3 PI-specified variables if dropped
################################################################################

# Core variables always retained: age, bp, diabetes
core_vars      <- c("age", "bp", "diabetes")
candidate_vars <- c("age", "bp", "diabetes", "chol", "smoke", "bmi", "hd", "med")

run_selection <- function(sex_group) {
  
  df <- analysis_baseline %>% filter(sex == sex_group)
  
  # --- Step 1: Univariable screening ---
  uni_results <- map_dfr(candidate_vars, function(v) {
    form <- as.formula(paste0("Surv(time, stroke) ~ ", v))
    fit  <- coxph(form, data = df)
    s    <- summary(fit)$coefficients
    data.frame(
      variable = v,
      coef     = round(s[1, "coef"], 4),
      HR       = round(exp(s[1, "coef"]), 4),
      p_value  = round(s[1, "Pr(>|z|)"], 4)
    )
  })
  
  cat("\n====", sex_group, "— Univariable Results ====\n")
  print(uni_results)
  
  # Variables passing p < 0.20
  selected_vars <- uni_results %>%
    filter(p_value < 0.20) %>%
    pull(variable)
  
  # Force-include core variables
  selected_vars <- union(core_vars, selected_vars)
  
  cat("\nVariables entering backward selection:", paste(selected_vars, collapse = ", "), "\n")
  
  # --- Step 2: Backward AIC ---
  full_formula <- as.formula(
    paste0("Surv(time, stroke) ~ ", paste(selected_vars, collapse = " + "))
  )
  full_model <- coxph(full_formula, data = df)
  
  # step() uses extractAIC for coxph; direction = "backward"
  final_model_aic <- step(full_model, direction = "backward", trace = 1)
  
  # --- Step 3: Force core variables back if removed ---
  final_vars <- names(coef(final_model_aic))
  # For factor variables, strip the level suffix to get base name
  final_base <- unique(gsub("(Yes|No|Male|Female|\\d+)", "", final_vars))
  final_base  <- trimws(final_base)
  
  missing_core <- setdiff(core_vars, names(coef(final_model_aic)))
  # More robust check using grepl
  core_in_model <- sapply(core_vars, function(cv) {
    any(grepl(cv, names(coef(final_model_aic)), ignore.case = TRUE))
  })
  missing_core <- core_vars[!core_in_model]
  
  if (length(missing_core) > 0) {
    cat("\nForce-adding core variables dropped by AIC:", paste(missing_core, collapse = ", "), "\n")
    all_final_vars <- union(
      names(coef(final_model_aic)) %>% gsub("(Yes|No|Male|Female)", "", .) %>% trimws(),
      missing_core
    )
    # Re-fit with forced core variables
    final_formula <- as.formula(
      paste0("Surv(time, stroke) ~ ",
             paste(union(selected_vars[selected_vars %in% names(coef(final_model_aic)) |
                                         sapply(selected_vars, function(x) any(grepl(x, names(coef(final_model_aic)))))],
                         missing_core),
                   collapse = " + "))
    )
    final_model_aic <- coxph(final_formula, data = df)
  }
  
  cat("\n====", sex_group, "— Final Model ====\n")
  print(summary(final_model_aic))
  
  return(final_model_aic)
}

cox_male   <- run_selection("Male")
cox_female <- run_selection("Female")

################################################################################
# SECTION 6: PH ASSUMPTION — SCHOENFELD RESIDUALS
################################################################################

ph_plots <- function(model, sex_label) {
  
  ph_test <- cox.zph(model)
  cat("\n====", sex_label, "— Schoenfeld PH Test ====\n")
  print(ph_test)
  
  # Save global table
  ph_df <- as.data.frame(ph_test$table) %>%
    rownames_to_column("Variable") %>%
    mutate(across(where(is.numeric), ~ round(., 4)))
  
  write_csv(ph_df, paste0("ph_test_", tolower(sex_label), ".csv"))
  
  # Plot all Schoenfeld residuals
  n_vars <- nrow(ph_test$table) - 1  # exclude GLOBAL
  
  png(paste0("schoenfeld_", tolower(sex_label), ".png"),
      width = 6 * ceiling(sqrt(n_vars + 1)),
      height = 5 * ceiling((n_vars + 1) / ceiling(sqrt(n_vars + 1))),
      res = 300)
  
  ggcoxzph(ph_test,
           point.col  = "steelblue",
           point.size = 0.8,
           point.alpha = 0.5,
           ggtheme    = theme_bw(base_size = 12)) %>%
    print()
  
  dev.off()
  
  invisible(ph_test)
}

ph_male   <- ph_plots(cox_male,   "Male")
ph_female <- ph_plots(cox_female, "Female")

################################################################################
# SECTION 7: 10-YEAR STROKE PROBABILITY TABLE
#
# Profile grid:
#   age      : 40, 50, 60
#   diabetes : Yes / No
#   hypertension (bp > 160): Yes (bp = 170) / No (bp = 120)
#   med      : Yes / No  (PI-chosen extra variable)
#   All other variables fixed at sample mean (continuous) or mode (categorical)
################################################################################

make_prob_table <- function(model, df, sex_label) {
  
  # Determine which variables are in the final model
  model_vars <- attr(model$terms, "term.labels")
  cat("\nFinal model variables for", sex_label, ":", paste(model_vars, collapse = ", "), "\n")
  
  # Reference values: mean for continuous, mode for categorical
  mode_val <- function(x) {
    ux <- na.omit(unique(x))
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  ref <- list(
    age      = mean(df$age,  na.rm = TRUE),
    chol     = mean(df$chol, na.rm = TRUE),
    bp       = mean(df$bp,   na.rm = TRUE),
    bmi      = mean(df$bmi,  na.rm = TRUE),
    smoke    = mode_val(df$smoke),
    diabetes = mode_val(df$diabetes),
    hd       = mode_val(df$hd),
    med      = mode_val(df$med)
  )
  
  # Build profile grid
  profiles <- expand.grid(
    age_group  = c(40, 50, 60),
    diab       = c("Yes", "No"),
    htn        = c("Yes", "No"),   # hypertension: Yes = bp 170, No = bp 120
    med_status = c("Yes", "No")
  ) %>%
    mutate(
      age      = age_group,
      bp       = ifelse(htn == "Yes", 170, 120),
      diabetes = factor(diab,       levels = c("No","Yes")),
      med      = factor(med_status, levels = c("No","Yes"))
    )
  
  # Fill remaining variables from reference
  for (v in setdiff(model_vars, c("age", "bp", "diabetes", "med"))) {
    profiles[[v]] <- ref[[v]]
  }
  
  # Predict 10-year survival probability (S(t=3652))
  # survfit with newdata gives survival curve per row
  probs <- map_dbl(seq_len(nrow(profiles)), function(i) {
    nd  <- profiles[i, , drop = FALSE]
    sf  <- survfit(model, newdata = nd)
    # Interpolate at day 3652 (10 years)
    t_idx <- which.min(abs(sf$time - 3652))
    surv_at_10 <- sf$surv[t_idx]
    round((1 - surv_at_10) * 100, 2)   # return as % stroke probability
  })
  
  result <- profiles %>%
    dplyr::select(age_group, diab, htn, med_status) %>%
    mutate(
      `10-yr Stroke Prob (%)` = probs
    ) %>%
    rename(
      Age             = age_group,
      Diabetes        = diab,
      `Hypertension (BP>160)` = htn,
      `BP Medication` = med_status
    )
  
  cat("\n====", sex_label, "— 10-Year Stroke Probability Table ====\n")
  print(result)
  
  write_csv(result, paste0("stroke_prob_", tolower(sex_label), ".csv"))
  
  return(result)
}

prob_male   <- make_prob_table(cox_male,   filter(analysis_baseline, sex == "Male"),   "Male")
prob_female <- make_prob_table(cox_female, filter(analysis_baseline, sex == "Female"), "Female")

################################################################################
# SECTION 8: COMBINED PROBABILITY TABLE (Male + Female side by side)
################################################################################

combined_prob <- prob_male %>%
  rename(`Male Prob (%)` = `10-yr Stroke Prob (%)`) %>%
  left_join(
    prob_female %>% rename(`Female Prob (%)` = `10-yr Stroke Prob (%)`),
    by = c("Age", "Diabetes", "Hypertension (BP>160)", "BP Medication")
  )

print(combined_prob)
write_csv(combined_prob, "stroke_prob_combined.csv")

# Formatted gtsummary-style table
library(gt)

combined_gt <- combined_prob %>%
  gt() %>%
  tab_header(
    title    = "10-Year Probability of Stroke by Risk Profile",
    subtitle = "Framingham Heart Study — Stratified by Sex"
  ) %>%
  tab_spanner(
    label   = "Risk Profile",
    columns = c(Age, Diabetes, `Hypertension (BP>160)`, `BP Medication`)
  ) %>%
  tab_spanner(
    label   = "10-Year Stroke Probability",
    columns = c(`Male Prob (%)`, `Female Prob (%)`)
  ) %>%
  fmt_number(
    columns  = c(`Male Prob (%)`, `Female Prob (%)`),
    decimals = 2
  ) %>%
  cols_label(
    `Male Prob (%)`   = "Male (%)",
    `Female Prob (%)` = "Female (%)"
  ) %>%
  tab_style(
    style     = cell_fill(color = "#f0f4ff"),
    locations = cells_body(rows = seq(1, nrow(combined_prob), 2))
  ) %>%
  tab_footnote(
    footnote = "Other covariates fixed at sample mean (continuous) or mode (categorical). Hypertension defined as SBP > 160 mmHg (Yes = 170, No = 120 used for prediction).",
    locations = cells_title()
  )

gtsave(combined_gt, "stroke_prob_combined.html")

################################################################################
# SECTION 9: FOREST PLOT — HAZARD RATIOS FROM FINAL MODELS
################################################################################

extract_hr <- function(model, sex_label) {
  s <- summary(model)$coefficients
  data.frame(
    sex      = sex_label,
    variable = rownames(s),
    HR       = exp(s[, "coef"]),
    lower    = exp(s[, "coef"] - 1.96 * s[, "se(coef)"]),
    upper    = exp(s[, "coef"] + 1.96 * s[, "se(coef)"]),
    p_value  = s[, "Pr(>|z|)"]
  )
}

hr_df <- bind_rows(
  extract_hr(cox_male,   "Male"),
  extract_hr(cox_female, "Female")
) %>%
  mutate(
    sig   = ifelse(p_value < 0.05, "p < 0.05", "p ≥ 0.05"),
    label = sprintf("%.2f (%.2f–%.2f)", HR, lower, upper)
  )

p_forest <- ggplot(hr_df, aes(x = HR, y = variable, color = sex, shape = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Male" = "#2166ac", "Female" = "#d6604d")) +
  scale_shape_manual(values = c("p < 0.05" = 16, "p ≥ 0.05" = 1)) +
  scale_x_log10() +
  labs(
    title    = "Hazard Ratios for 10-Year Stroke Risk",
    subtitle = "Cox Proportional Hazards Model — Stratified by Sex",
    x        = "Hazard Ratio (log scale)",
    y        = NULL,
    color    = "Sex",
    shape    = "Significance"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(face = "bold")
  )

ggsave("forest_plot_hr.png", plot = p_forest,
       width = 10, height = 6, dpi = 300, bg = "white")

cat("\n\n=== ANALYSIS COMPLETE ===\n")
cat("Output files generated:\n")
cat("  - stroke_prob_combined.csv / .html\n")
cat("  - forest_plot_hr.png\n")
cat("  - schoenfeld_male.png / schoenfeld_female.png\n")
cat("  - ph_test_male.csv / ph_test_female.csv\n")
cat("  - km_*.png (KM curves per variable)\n")