# Patient-level outer CV plus final evaluation on untouched holdout patients.

source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "03_model_helpers.R"))

interval_split <- readr::read_csv(
  here::here("data", "processed", "interval_split.csv"), show_col_types = FALSE
) |>
  coerce_predictor_types()
analysis_split <- readr::read_csv(
  here::here("data", "processed", "analysis_split.csv"), show_col_types = FALSE
) |>
  coerce_predictor_types()
final_models <- readRDS(here::here("outputs", "models", "nine_models.rds"))

summarise_predictions <- function(predictions) {
  predictions |>
    dplyr::group_by(framework, selection, model) |>
    dplyr::group_modify(~ evaluate_predictions(
      .x$observed, .x$predicted, config$decision_threshold
    )) |>
    dplyr::ungroup()
}

cv_predictions <- purrr::map_dfr(seq_len(config$outer_folds), function(fold) {
  message("Outer CV fold ", fold, " of ", config$outer_folds)

  interval_train <- interval_split |>
    dplyr::filter(split == "training", outer_fold != fold)
  interval_validation <- interval_split |>
    dplyr::filter(split == "training", outer_fold == fold)
  longitudinal_train <- analysis_split |>
    dplyr::filter(split == "training", outer_fold != fold)
  longitudinal_validation <- analysis_split |>
    dplyr::filter(split == "training", outer_fold == fold)

  landmark_train <- build_landmark_data(longitudinal_train)
  landmark_validation <- build_landmark_data(longitudinal_validation)
  models <- fit_nine_models(
    interval_train, landmark_train, seed = config$seed + 1000L * fold
  )

  predict_nine_models(models, interval_validation, landmark_validation) |>
    dplyr::mutate(set = "Cross-validation", fold = fold)
})

holdout_intervals <- dplyr::filter(interval_split, split == "holdout")
holdout_longitudinal <- dplyr::filter(analysis_split, split == "holdout")
holdout_landmarks <- build_landmark_data(holdout_longitudinal)

holdout_predictions <- predict_nine_models(
  final_models, holdout_intervals, holdout_landmarks
) |>
  dplyr::mutate(set = "Holdout", fold = NA_integer_)

cv_fold_detail <- cv_predictions |>
  dplyr::group_by(fold) |>
  dplyr::group_modify(~ summarise_predictions(.x)) |>
  dplyr::ungroup() |>
  dplyr::mutate(set = "CV fold")

cv_mean <- cv_fold_detail |>
  dplyr::group_by(framework, selection, model) |>
  dplyr::summarise(
    dplyr::across(
      c(n, events, c_statistic, calibration_slope, calibration_intercept,
        brier_score, net_benefit),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(set = "CV mean", fold = NA_integer_)

holdout_summary <- summarise_predictions(holdout_predictions) |>
  dplyr::mutate(set = "Holdout", fold = NA_integer_)

performance_summary <- dplyr::bind_rows(cv_mean, holdout_summary) |>
  dplyr::select(
    model, framework, selection, set, fold, n, events,
    c_statistic, calibration_slope, calibration_intercept,
    brier_score, net_benefit
  ) |>
  dplyr::arrange(framework, selection, set) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 4)))

all_predictions <- dplyr::bind_rows(cv_predictions, holdout_predictions)

decision_curve_results <- holdout_predictions |>
  dplyr::group_by(framework, selection, model) |>
  dplyr::group_modify(~ decision_curve(
    .x$observed, .x$predicted, seq(0.01, 0.20, by = 0.01)
  )) |>
  dplyr::ungroup()

dir.create(here::here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(performance_summary, here::here("outputs", "tables", "performance_summary.csv"))
readr::write_csv(cv_fold_detail, here::here("outputs", "tables", "cv_fold_detail.csv"))
readr::write_csv(all_predictions, here::here("outputs", "tables", "all_predictions.csv"))
readr::write_csv(
  decision_curve_results,
  here::here("outputs", "tables", "decision_curve_analysis.csv")
)

message("Evaluation complete for all nine model combinations.")
