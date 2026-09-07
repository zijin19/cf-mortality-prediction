# Clean the synthetic source tables and construct two-year interval outcomes.

source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))

patient <- readr::read_csv(here::here("data", "raw", "patient.csv"), show_col_types = FALSE)
encounter <- readr::read_csv(
  here::here("data", "raw", "annual_encounter.csv"), show_col_types = FALSE
)

encounter_clean <- encounter |>
  dplyr::filter(
    !is.na(fev1_pct_predicted), !is.na(fvc_pct_predicted), !is.na(bmi_percentile),
    dplyr::between(fev1_pct_predicted, 8, 150),
    dplyr::between(fvc_pct_predicted, 8, 150),
    dplyr::between(bmi_percentile, 1, 100)
  ) |>
  dplyr::mutate(
    iv_antibiotic_cat = dplyr::case_when(
      iv_antibiotic_times == 0 ~ "0",
      iv_antibiotic_times == 1 ~ "1",
      iv_antibiotic_times == 2 ~ "2",
      TRUE ~ "3+"
    )
  )

analysis_data <- encounter_clean |>
  dplyr::inner_join(patient, by = "patient_id") |>
  coerce_predictor_types() |>
  dplyr::arrange(patient_id, current_age)

horizon <- config$prediction_horizon

interval_data <- analysis_data |>
  dplyr::mutate(
    event_2yr = as.integer(
      event_indicator == 1 & !is.na(event_age) &
        event_age > current_age & event_age <= current_age + horizon
    ),
    has_complete_window = event_2yr == 1 | followup_age >= current_age + horizon,
    followup_time = dplyr::if_else(
      event_2yr == 1, event_age - current_age, as.numeric(horizon)
    ),
    age_exit = current_age + followup_time
  ) |>
  dplyr::filter(has_complete_window, followup_time > 0) |>
  dplyr::select(
    patient_id, report_year, current_age, age_exit, followup_time, event_2yr,
    dplyr::all_of(predictors), event_indicator, event_type, event_age, followup_age
  )

dir.create(here::here("data", "processed"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(analysis_data, here::here("data", "processed", "analysis_data.csv"))
readr::write_csv(interval_data, here::here("data", "processed", "interval_data.csv"))

message(
  "Analysis data: ", nrow(analysis_data), " annual records. Two-year data: ",
  nrow(interval_data), " intervals; event rate ",
  scales::percent(mean(interval_data$event_2yr), accuracy = 0.01), "."
)
