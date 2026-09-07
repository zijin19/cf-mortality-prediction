# Shared analysis settings. Values can be overridden with environment variables.

config <- list(
  seed = 2024L,
  split_seed = 123L,
  n_patients = as.integer(Sys.getenv("CF_N_PATIENTS", "1000")),
  study_years = 2011:2019,
  prediction_horizon = 2,
  landmark_ages = 6:50,
  holdout_fraction = 0.10,
  outer_folds = 5L,
  inner_folds = 5L,
  lasso_nlambda = as.integer(Sys.getenv("CF_LASSO_NLAMBDA", "30")),
  decision_threshold = 0.50,
  min_landmark_at_risk = 40L,
  min_landmark_events = 5L
)

predictors <- c(
  "age_at_diagnosis", "current_age", "sex", "df508",
  "pancreatic_insufficiency", "cfrd", "bmi_percentile",
  "fev1_pct_predicted", "fvc_pct_predicted", "iv_antibiotic_cat",
  "home_oxygen", "b_cepacia", "p_aeruginosa", "s_aureus",
  "corticosteroids", "depression"
)

# The final thesis landmark models excluded age variables because landmark age
# already supplies the time scale.
landmark_predictors <- setdiff(
  predictors, c("age_at_diagnosis", "current_age")
)

factor_levels <- list(
  sex = c("Female", "Male"),
  df508 = c("Other", "Heterozygous", "Homozygous"),
  iv_antibiotic_cat = c("0", "1", "2", "3+")
)
