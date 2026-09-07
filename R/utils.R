clip <- function(x, lower, upper) pmin(pmax(x, lower), upper)

expit <- function(x) 1 / (1 + exp(-x))

coerce_predictor_types <- function(data) {
  data$sex <- factor(data$sex, levels = factor_levels$sex)
  data$df508 <- factor(data$df508, levels = factor_levels$df508)
  data$iv_antibiotic_cat <- factor(
    data$iv_antibiotic_cat,
    levels = factor_levels$iv_antibiotic_cat,
    ordered = TRUE
  )
  data
}

make_stratified_patient_folds <- function(data, k, seed, outcome = "event_2yr") {
  patient_outcome <- data |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(stratum = max(.data[[outcome]], na.rm = TRUE), .groups = "drop")

  set.seed(seed)
  patient_outcome <- patient_outcome |>
    dplyr::group_by(stratum) |>
    dplyr::mutate(fold = sample(rep(seq_len(k), length.out = dplyr::n()))) |>
    dplyr::ungroup()

  patient_outcome$fold[match(data$patient_id, patient_outcome$patient_id)]
}

safe_auc <- function(observed, predicted) {
  if (length(unique(observed)) < 2L) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(observed, predicted, quiet = TRUE)))
}

calibration_stats <- function(observed, predicted) {
  if (length(unique(observed)) < 2L) {
    return(c(slope = NA_real_, intercept = NA_real_))
  }
  eps <- 1e-6
  logit_prediction <- qlogis(clip(predicted, eps, 1 - eps))
  fit <- tryCatch(
    glm(observed ~ logit_prediction, family = binomial()),
    error = function(e) NULL
  )
  if (is.null(fit)) return(c(slope = NA_real_, intercept = NA_real_))
  c(slope = unname(coef(fit)[2]), intercept = unname(coef(fit)[1]))
}

net_benefit <- function(observed, predicted, threshold) {
  positive <- predicted >= threshold
  tp <- sum(positive & observed == 1)
  fp <- sum(positive & observed == 0)
  (tp - fp * threshold / (1 - threshold)) / length(observed)
}

evaluate_predictions <- function(observed, predicted, threshold) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]
  predicted <- clip(predicted[keep], 0, 1)
  cal <- calibration_stats(observed, predicted)

  tibble::tibble(
    n = length(observed),
    events = sum(observed),
    c_statistic = safe_auc(observed, predicted),
    calibration_slope = unname(cal["slope"]),
    calibration_intercept = unname(cal["intercept"]),
    brier_score = mean((predicted - observed)^2),
    net_benefit = net_benefit(observed, predicted, threshold)
  )
}

decision_curve <- function(observed, predicted, thresholds) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]
  predicted <- predicted[keep]
  prevalence <- mean(observed)

  purrr::map_dfr(thresholds, function(threshold) {
    tibble::tibble(
      threshold = threshold,
      net_benefit_model = net_benefit(observed, predicted, threshold),
      net_benefit_treat_all = prevalence -
        (1 - prevalence) * threshold / (1 - threshold),
      net_benefit_treat_none = 0
    )
  })
}

model_label <- function(framework, selection) {
  paste(
    dplyr::recode(framework,
      logistic = "Logistic", cox = "Cox", landmark = "Landmark Cox"
    ),
    dplyr::recode(selection,
      forward = "Forward", backward = "Backward", group_lasso = "Group lasso"
    ),
    sep = " - "
  )
}
