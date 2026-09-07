# Fit all nine models to the complete training partition, leaving holdout
# patients untouched until final evaluation.

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

training_intervals <- dplyr::filter(interval_split, split == "training")
training_longitudinal <- dplyr::filter(analysis_split, split == "training")
training_landmarks <- build_landmark_data(training_longitudinal)

final_models <- fit_nine_models(
  training_intervals, training_landmarks, seed = config$seed
)

dir.create(here::here("outputs", "models"), recursive = TRUE, showWarnings = FALSE)
saveRDS(final_models, here::here("outputs", "models", "nine_models.rds"))

landmark_counts <- training_landmarks |>
  dplyr::group_by(landmark_age) |>
  dplyr::summarise(n_at_risk = dplyr::n(), events = sum(event_2yr), .groups = "drop") |>
  dplyr::mutate(
    used = landmark_age %in% final_models$landmark_forward$eligible_ages
  )
dir.create(here::here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  landmark_counts, here::here("outputs", "tables", "landmark_training_counts.csv")
)

message("Saved all nine final fitted models.")
