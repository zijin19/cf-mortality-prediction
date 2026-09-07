# Create an untouched 10% holdout and five outer CV folds at patient level.

source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))

interval_data <- readr::read_csv(
  here::here("data", "processed", "interval_data.csv"), show_col_types = FALSE
) |>
  coerce_predictor_types()
analysis_data <- readr::read_csv(
  here::here("data", "processed", "analysis_data.csv"), show_col_types = FALSE
) |>
  coerce_predictor_types()

patient_strata <- interval_data |>
  dplyr::group_by(patient_id) |>
  dplyr::summarise(ever_event = max(event_indicator), .groups = "drop")

set.seed(config$split_seed)
split_map <- patient_strata |>
  dplyr::group_by(ever_event) |>
  dplyr::mutate(
    shuffled = sample.int(dplyr::n()),
    holdout = shuffled <= pmax(1L, round(dplyr::n() * config$holdout_fraction))
  ) |>
  dplyr::ungroup() |>
  dplyr::select(-shuffled)

training_patients <- split_map |>
  dplyr::filter(!holdout) |>
  dplyr::group_by(ever_event) |>
  dplyr::mutate(outer_fold = sample(rep(seq_len(config$outer_folds), length.out = dplyr::n()))) |>
  dplyr::ungroup()

split_map <- split_map |>
  dplyr::left_join(
    training_patients |> dplyr::select(patient_id, outer_fold),
    by = "patient_id"
  ) |>
  dplyr::mutate(split = ifelse(holdout, "holdout", "training")) |>
  dplyr::select(patient_id, ever_event, split, outer_fold)

interval_split <- interval_data |>
  dplyr::inner_join(split_map, by = "patient_id")
analysis_split <- analysis_data |>
  dplyr::inner_join(split_map, by = "patient_id")

readr::write_csv(split_map, here::here("data", "processed", "patient_split.csv"))
readr::write_csv(interval_split, here::here("data", "processed", "interval_split.csv"))
readr::write_csv(analysis_split, here::here("data", "processed", "analysis_split.csv"))

message(
  "Patient split: ", sum(split_map$split == "training"), " training; ",
  sum(split_map$split == "holdout"), " holdout."
)
