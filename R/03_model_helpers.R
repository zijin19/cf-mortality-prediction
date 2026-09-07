# Model construction and prediction helpers for the full 3 x 3 comparison.

source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))

build_landmark_data <- function(data, ages = config$landmark_ages,
                                horizon = config$prediction_horizon) {
  data <- coerce_predictor_types(data)

  purrr::map_dfr(ages, function(landmark_age) {
    at_risk <- data |>
      dplyr::filter(
        current_age <= landmark_age,
        followup_age > landmark_age,
        is.na(event_age) | event_age > landmark_age
      ) |>
      dplyr::group_by(patient_id) |>
      dplyr::slice_max(current_age, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        landmark_age = landmark_age,
        event_2yr = as.integer(
          event_indicator == 1 & !is.na(event_age) &
            event_age <= landmark_age + horizon
        ),
        has_complete_window = event_2yr == 1 | followup_age >= landmark_age + horizon,
        followup_time = dplyr::if_else(
          event_2yr == 1, event_age - landmark_age, as.numeric(horizon)
        )
      ) |>
      dplyr::filter(has_complete_window, followup_time > 0)

    at_risk
  }) |>
    dplyr::arrange(landmark_age, patient_id)
}

eligible_landmark_ages <- function(data) {
  data |>
    dplyr::group_by(landmark_age) |>
    dplyr::summarise(
      n_at_risk = dplyr::n(), events = sum(event_2yr), .groups = "drop"
    ) |>
    dplyr::filter(
      n_at_risk >= config$min_landmark_at_risk,
      events >= config$min_landmark_events
    ) |>
    dplyr::pull(landmark_age)
}

make_stepwise_formula <- function(framework, terms, null = FALSE) {
  rhs <- if (null) character(0) else terms
  if (framework == "logistic") {
    return(stats::reformulate(rhs, response = "event_2yr"))
  }

  survival_response <- if (framework == "cox") {
    "survival::Surv(current_age, age_exit, event_2yr)"
  } else {
    "survival::Surv(followup_time, event_2yr)"
  }
  # Keep cluster outside the formula. coxph() promotes cluster() from a formula
  # to a separate argument, while stepAIC() still sees the original scope term.
  # That mismatch makes forward/backward scope checks fail.
  specials <- if (framework == "landmark") "strata(landmark_age)" else character(0)
  rhs_all <- c(rhs, specials)
  if (length(rhs_all) == 0L) rhs_all <- "1"
  stats::as.formula(paste(
    survival_response, "~",
    paste(rhs_all, collapse = " + ")
  ))
}

fit_stepwise <- function(data, framework, direction) {
  terms <- if (framework == "landmark") landmark_predictors else predictors
  full_formula <- make_stepwise_formula(framework, terms, null = FALSE)
  null_formula <- make_stepwise_formula(framework, terms, null = TRUE)
  # stepAIC() refits models from their stored calls. Bind both formulas and the
  # data to this function frame so the symbol `model_data` cannot fall through
  # to utils::data() during those updates.
  environment(full_formula) <- environment()
  environment(null_formula) <- environment()
  model_data <- as.data.frame(data)
  penalty <- log(dplyr::n_distinct(model_data$patient_id))

  if (framework == "logistic") {
    full_model <- glm(full_formula, data = model_data, family = binomial())
    if (direction == "backward") {
      return(MASS::stepAIC(
        full_model,
        scope = list(lower = null_formula, upper = full_formula),
        direction = "backward", trace = 0, k = penalty
      ))
    }
    null_model <- glm(null_formula, data = model_data, family = binomial())
    return(MASS::stepAIC(
      null_model,
      scope = list(lower = null_formula, upper = full_formula),
      direction = "forward", trace = 0, k = penalty
    ))
  }

  # MASS::stepAIC() delegates scope handling to factor.scope(), which cannot
  # reliably process the special strata() term used by the landmark model.
  # Use an explicit grouped BIC search for both Cox frameworks instead. Each
  # candidate here is a complete predictor term, so factors enter or leave as
  # a group and strata(landmark_age) is retained in every landmark candidate.
  selected_terms <- select_cox_terms_bic(
    model_data, framework, terms, direction, penalty
  )

  # Clustering changes the variance estimate, not the Cox partial likelihood
  # used for selection, so add patient clustering only in the final refit.
  final_model <- fit_cox_candidate(
    model_data, framework, selected_terms, clustered = TRUE
  )
  if (has_cox_convergence_warning(final_model)) {
    stop(
      "The selected ", framework,
      " model triggered a convergence warning during the clustered refit."
    )
  }
  final_model
}

fit_cox_candidate <- function(data, framework, selected_terms,
                              clustered = FALSE) {
  formula <- make_stepwise_formula(
    framework, selected_terms, null = length(selected_terms) == 0L
  )
  environment(formula) <- environment()

  arguments <- list(
    formula = formula,
    data = data,
    ties = "breslow",
    model = TRUE,
    x = TRUE
  )
  if (clustered) arguments$cluster <- data$patient_id

  # Sparse event/predictor combinations can make a Cox coefficient diverge.
  # Capture only convergence-related warnings so the search can reject those
  # candidates explicitly; unrelated warnings remain visible to the user.
  convergence_warnings <- character(0)
  fit <- withCallingHandlers(
    do.call(survival::coxph, arguments),
    warning = function(warning) {
      message <- conditionMessage(warning)
      convergence_pattern <- paste(
        "coefficient may be infinite",
        "did not converge",
        "ran out of iterations",
        sep = "|"
      )
      if (grepl(convergence_pattern, message, ignore.case = TRUE)) {
        convergence_warnings <<- c(convergence_warnings, message)
        invokeRestart("muffleWarning")
      }
    }
  )
  attr(fit, "convergence_warnings") <- unique(convergence_warnings)
  fit
}

has_cox_convergence_warning <- function(model) {
  length(attr(model, "convergence_warnings", exact = TRUE)) > 0L
}

cox_partial_bic <- function(model, penalty) {
  log_likelihood <- as.numeric(stats::logLik(model))
  -2 * log_likelihood + penalty * length(stats::coef(model))
}

select_cox_terms_bic <- function(data, framework, candidate_terms,
                                 direction, penalty) {
  selected_terms <- if (direction == "forward") character(0) else candidate_terms
  current_model <- fit_cox_candidate(
    data, framework, selected_terms, clustered = FALSE
  )
  current_bic <- cox_partial_bic(current_model, penalty)
  current_valid <- !has_cox_convergence_warning(current_model)

  repeat {
    terms_to_change <- if (direction == "forward") {
      setdiff(candidate_terms, selected_terms)
    } else {
      selected_terms
    }
    if (length(terms_to_change) == 0L) break

    proposed_term_sets <- lapply(terms_to_change, function(term) {
      if (direction == "forward") {
        c(selected_terms, term)
      } else {
        setdiff(selected_terms, term)
      }
    })
    proposed_models <- lapply(proposed_term_sets, function(term_set) {
      tryCatch(
        fit_cox_candidate(data, framework, term_set, clustered = FALSE),
        error = function(error) NULL
      )
    })
    proposed_bic <- vapply(proposed_models, function(model) {
      if (is.null(model)) Inf else cox_partial_bic(model, penalty)
    }, numeric(1))
    proposed_valid <- vapply(proposed_models, function(model) {
      !is.null(model) && !has_cox_convergence_warning(model)
    }, logical(1))

    if (!current_valid) {
      # A saturated backward starting model can be unstable. Prefer a stable
      # one-term deletion immediately; if none exists yet, keep simplifying
      # through the best finite candidate until a stable fit is reached.
      eligible <- which(proposed_valid)
      if (length(eligible) == 0L) eligible <- which(is.finite(proposed_bic))
      if (length(eligible) == 0L) {
        stop("No finite Cox candidate was available during ", direction, " selection.")
      }
      best <- eligible[which.min(proposed_bic[eligible])]
    } else {
      valid_bic <- proposed_bic
      valid_bic[!proposed_valid] <- Inf
      best <- which.min(valid_bic)
      if (!is.finite(valid_bic[best]) || valid_bic[best] >= current_bic - 1e-8) {
        break
      }
    }
    selected_terms <- proposed_term_sets[[best]]
    current_model <- proposed_models[[best]]
    current_bic <- proposed_bic[best]
    current_valid <- proposed_valid[best]
  }

  if (!current_valid) {
    stop(
      "Cox ", direction,
      " selection could not find a model without a convergence warning."
    )
  }
  selected_terms
}

fit_landmark_stepwise <- function(data, direction) {
  ages <- eligible_landmark_ages(data)
  if (length(ages) == 0L) stop("No landmark ages have enough training events.")
  fit_data <- dplyr::filter(data, landmark_age %in% ages)
  list(
    framework = "landmark",
    selection = direction,
    eligible_ages = ages,
    model = fit_stepwise(fit_data, "landmark", direction)
  )
}

make_design_spec <- function(data, terms) {
  data <- coerce_predictor_types(data)
  design_formula <- stats::reformulate(terms)
  design_terms <- stats::terms(design_formula, data = data)
  matrix_full <- stats::model.matrix(design_terms, data = data)
  assignment <- attr(matrix_full, "assign")
  keep <- colnames(matrix_full) != "(Intercept)"
  matrix <- matrix_full[, keep, drop = FALSE]
  assignment <- assignment[keep]
  term_labels <- attr(design_terms, "term.labels")
  group_names <- term_labels[assignment]
  groups <- match(group_names, unique(group_names))

  center <- colMeans(matrix)
  scale <- apply(matrix, 2, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1

  list(
    terms = design_terms,
    columns = colnames(matrix),
    center = center,
    scale = scale,
    group = groups,
    group_names = group_names
  )
}

apply_design_spec <- function(spec, data) {
  data <- coerce_predictor_types(data)
  matrix_full <- stats::model.matrix(spec$terms, data = data)
  matrix <- matrix_full[, colnames(matrix_full) != "(Intercept)", drop = FALSE]

  missing_columns <- setdiff(spec$columns, colnames(matrix))
  if (length(missing_columns) > 0L) {
    missing_matrix <- matrix(0, nrow(matrix), length(missing_columns))
    colnames(missing_matrix) <- missing_columns
    matrix <- cbind(matrix, missing_matrix)
  }
  matrix <- matrix[, spec$columns, drop = FALSE]
  matrix <- sweep(matrix, 2, spec$center, "-")
  sweep(matrix, 2, spec$scale, "/")
}

fit_group_logistic <- function(data, seed) {
  spec <- make_design_spec(data, predictors)
  x <- apply_design_spec(spec, data)
  y <- ifelse(data$event_2yr == 1, 1, -1)
  fold_id <- make_stratified_patient_folds(
    data, config$inner_folds, seed, outcome = "event_2yr"
  )

  cv <- gglasso::cv.gglasso(
    x = x, y = y, group = spec$group, loss = "logit", pred.loss = "loss",
    foldid = fold_id, nfolds = config$inner_folds,
    nlambda = config$lasso_nlambda
  )

  list(
    framework = "logistic", selection = "group_lasso", spec = spec,
    cv = cv, lambda = cv$lambda.1se
  )
}

fit_group_cox <- function(data, terms, seed, framework = "cox") {
  spec <- make_design_spec(data, terms)
  x <- apply_design_spec(spec, data)
  y <- cbind(data$followup_time, data$event_2yr)
  fold_id <- make_stratified_patient_folds(
    data, config$inner_folds, seed, outcome = "event_2yr"
  )

  cv <- grpreg::cv.grpsurv(
    X = x, y = y, group = spec$group, penalty = "grLasso",
    nfolds = config$inner_folds, fold = fold_id, seed = seed,
    nlambda = config$lasso_nlambda
  )

  list(
    framework = framework, selection = "group_lasso", spec = spec,
    cv = cv, lambda = cv$lambda.min
  )
}

fit_landmark_group_lasso <- function(data, seed) {
  ages <- eligible_landmark_ages(data)
  models <- purrr::map(
    ages,
    function(age) {
      age_data <- dplyr::filter(data, landmark_age == age)
      fit_group_cox(
        age_data, landmark_predictors, seed + as.integer(age),
        framework = "landmark"
      )
    }
  )
  names(models) <- as.character(ages)
  list(
    framework = "landmark", selection = "group_lasso",
    eligible_ages = ages, models = models
  )
}

predict_stepwise_risk <- function(model, data, framework,
                                  horizon = config$prediction_horizon) {
  if (framework == "logistic") {
    return(as.numeric(stats::predict(model, newdata = data, type = "response")))
  }

  # A forward BIC search can legitimately select a null Cox model, especially
  # in an outer-CV training fold with few events. coxph stores NULL rather than
  # numeric(0) coefficients for that model, but predict.coxph(type="expected")
  # performs a matrix product and requires a numeric zero-length vector.
  # Normalising the stored value preserves the null model and lets coxph use
  # its fitted baseline hazard (including landmark-specific strata).
  prediction_model <- model
  if (length(stats::coef(prediction_model)) == 0L) {
    prediction_model$coefficients <- numeric(0)
  }
  prediction_data <- data
  prediction_data$followup_time <- horizon
  if (framework == "cox") {
    prediction_data$age_exit <- prediction_data$current_age + horizon
  }
  expected_events <- as.numeric(stats::predict(
    prediction_model, newdata = prediction_data, type = "expected"
  ))
  clip(1 - exp(-expected_events), 0, 1)
}

predict_group_logistic <- function(fit, data) {
  x <- apply_design_spec(fit$spec, data)
  linear_predictor <- drop(stats::predict(
    fit$cv$gglasso.fit, x, s = fit$lambda, type = "link"
  ))
  clip(expit(linear_predictor), 0, 1)
}

evaluate_survival_functions <- function(functions, horizon, n_expected) {
  if (is.function(functions)) return(1 - functions(horizon))
  if (length(functions) != n_expected) {
    stop("Unexpected number of survival functions returned by grpreg.")
  }
  vapply(functions, function(f) 1 - f(horizon), numeric(1))
}

predict_group_cox <- function(fit, data, horizon = config$prediction_horizon) {
  x <- apply_design_spec(fit$spec, data)
  survival_functions <- stats::predict(
    fit$cv$fit, X = x, type = "survival", lambda = fit$lambda
  )
  clip(
    evaluate_survival_functions(survival_functions, horizon, nrow(x)),
    0, 1
  )
}

fit_nine_models <- function(interval_data, landmark_data, seed = config$seed) {
  message("Fitting Logistic models...")
  logistic_forward <- fit_stepwise(interval_data, "logistic", "forward")
  logistic_backward <- fit_stepwise(interval_data, "logistic", "backward")
  logistic_group_lasso <- fit_group_logistic(interval_data, seed)

  message("Fitting Cox models...")
  cox_forward <- fit_stepwise(interval_data, "cox", "forward")
  cox_backward <- fit_stepwise(interval_data, "cox", "backward")
  cox_group_lasso <- fit_group_cox(interval_data, predictors, seed + 100L)

  message("Fitting Landmark Cox models...")
  landmark_forward <- fit_landmark_stepwise(landmark_data, "forward")
  landmark_backward <- fit_landmark_stepwise(landmark_data, "backward")
  landmark_group_lasso <- fit_landmark_group_lasso(landmark_data, seed + 200L)

  list(
    logistic_forward = logistic_forward,
    logistic_backward = logistic_backward,
    logistic_group_lasso = logistic_group_lasso,
    cox_forward = cox_forward,
    cox_backward = cox_backward,
    cox_group_lasso = cox_group_lasso,
    landmark_forward = landmark_forward,
    landmark_backward = landmark_backward,
    landmark_group_lasso = landmark_group_lasso
  )
}

prediction_frame <- function(data, predicted, framework, selection) {
  tibble::tibble(
    patient_id = data$patient_id,
    landmark_age = if ("landmark_age" %in% names(data)) data$landmark_age else NA_real_,
    observed = data$event_2yr,
    predicted = as.numeric(predicted),
    framework = framework,
    selection = selection,
    model = model_label(framework, selection)
  )
}

predict_landmark_stepwise <- function(fit, data) {
  evaluation <- dplyr::filter(data, landmark_age %in% fit$eligible_ages)
  prediction_frame(
    evaluation,
    predict_stepwise_risk(fit$model, evaluation, "landmark"),
    "landmark", fit$selection
  )
}

predict_landmark_group_lasso <- function(fit, data) {
  purrr::map_dfr(fit$eligible_ages, function(age) {
    evaluation <- dplyr::filter(data, landmark_age == age)
    if (nrow(evaluation) == 0L) return(tibble::tibble())
    prediction_frame(
      evaluation,
      predict_group_cox(fit$models[[as.character(age)]], evaluation),
      "landmark", "group_lasso"
    )
  })
}

predict_nine_models <- function(models, interval_data, landmark_data) {
  dplyr::bind_rows(
    prediction_frame(
      interval_data,
      predict_stepwise_risk(models$logistic_forward, interval_data, "logistic"),
      "logistic", "forward"
    ),
    prediction_frame(
      interval_data,
      predict_stepwise_risk(models$logistic_backward, interval_data, "logistic"),
      "logistic", "backward"
    ),
    prediction_frame(
      interval_data, predict_group_logistic(models$logistic_group_lasso, interval_data),
      "logistic", "group_lasso"
    ),
    prediction_frame(
      interval_data,
      predict_stepwise_risk(models$cox_forward, interval_data, "cox"),
      "cox", "forward"
    ),
    prediction_frame(
      interval_data,
      predict_stepwise_risk(models$cox_backward, interval_data, "cox"),
      "cox", "backward"
    ),
    prediction_frame(
      interval_data, predict_group_cox(models$cox_group_lasso, interval_data),
      "cox", "group_lasso"
    ),
    predict_landmark_stepwise(models$landmark_forward, landmark_data),
    predict_landmark_stepwise(models$landmark_backward, landmark_data),
    predict_landmark_group_lasso(models$landmark_group_lasso, landmark_data)
  )
}
