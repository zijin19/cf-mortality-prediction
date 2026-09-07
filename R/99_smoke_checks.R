# Fail loudly if core reproducibility or leakage safeguards are violated.

split_map <- readr::read_csv(
  here::here("data", "processed", "patient_split.csv"), show_col_types = FALSE
)
performance <- readr::read_csv(
  here::here("outputs", "tables", "performance_summary.csv"), show_col_types = FALSE
)
predictions <- readr::read_csv(
  here::here("outputs", "tables", "all_predictions.csv"), show_col_types = FALSE
)
final_models <- readRDS(
  here::here("outputs", "models", "nine_models.rds")
)

stepwise_cox_models <- list(
  final_models$cox_forward,
  final_models$cox_backward,
  final_models$landmark_forward$model,
  final_models$landmark_backward$model
)

stopifnot(
  nrow(split_map) == dplyr::n_distinct(split_map$patient_id),
  !any(is.na(split_map$split)),
  dplyr::n_distinct(performance$model) == 9L,
  dplyr::n_distinct(dplyr::filter(predictions, set == "Holdout")$model) == 9L,
  all(predictions$predicted >= 0 & predictions$predicted <= 1, na.rm = TRUE),
  all(sort(unique(stats::na.omit(predictions$fold))) == seq_len(config$outer_folds)),
  all(vapply(stepwise_cox_models, function(model) {
    length(attr(model, "convergence_warnings", exact = TRUE)) == 0L
  }, logical(1)))
)

training_ids <- split_map$patient_id[split_map$split == "training"]
holdout_ids <- split_map$patient_id[split_map$split == "holdout"]
stopifnot(length(intersect(training_ids, holdout_ids)) == 0L)

message(paste(
  "Smoke checks passed: nine models, stable Cox fits, valid risks,",
  "and no patient split overlap."
))
