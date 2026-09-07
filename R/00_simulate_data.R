# Generate fully synthetic longitudinal registry data. The values are informed
# only by aggregate summaries in the thesis; no patient record is reconstructed.

source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))

set.seed(config$seed)

n_patients <- config$n_patients
study_years <- config$study_years

patient_base <- tibble::tibble(
  patient_id = seq_len(n_patients),
  age_in_2011 = round(clip(rnorm(n_patients, mean = 25, sd = 13), 6, 60)),
  sex = sample(c("Female", "Male"), n_patients, TRUE, c(0.463, 0.537)),
  age_at_diagnosis = round(clip(
    ifelse(
      runif(n_patients) < 0.65,
      runif(n_patients, 0.05, 1),
      rgamma(n_patients, shape = 2, scale = 3)
    ),
    0, 18
  ), 2),
  df508 = sample(
    c("Homozygous", "Heterozygous", "Other"), n_patients, TRUE,
    c(0.456, 0.408, 0.136)
  ),
  pancreatic_insufficiency = rbinom(n_patients, 1, 0.827),
  patient_effect = rnorm(n_patients)
) |>
  dplyr::mutate(
    birth_year = 2011L - age_in_2011,
    latent_severity = 0.45 * patient_effect +
      0.25 * pancreatic_insufficiency +
      0.15 * (df508 == "Homozygous")
  )

simulate_patient <- function(i) {
  p <- patient_base[i, ]
  years <- study_years
  ages <- p$age_in_2011 + seq_along(years) - 1
  time_index <- seq_along(years) - 1

  burden <- as.numeric(p$latent_severity) + 0.055 * time_index +
    cumsum(rnorm(length(years), 0, 0.08))
  fev1 <- clip(79 - 12 * burden - 0.33 * pmax(ages - 25, 0) + rnorm(length(years), 0, 7), 8, 150)
  fvc <- clip(fev1 + 13 + rnorm(length(years), 0, 6), 8, 150)
  bmi <- clip(55 - 8 * burden - 0.12 * pmax(ages - 30, 0) + rnorm(length(years), 0, 13), 1, 99)

  cfrd <- integer(length(years))
  cfrd[1] <- rbinom(1, 1, expit(-2.0 + 0.04 * ages[1] + 0.45 * burden[1]))
  p_aeruginosa <- integer(length(years))
  p_aeruginosa[1] <- rbinom(1, 1, expit(-0.15 + 0.35 * burden[1]))
  if (length(years) > 1L) {
    for (j in 2:length(years)) {
      cfrd[j] <- if (cfrd[j - 1] == 1) 1L else
        rbinom(1, 1, expit(-4.2 + 0.03 * ages[j] + 0.45 * burden[j]))
      p_aeruginosa[j] <- rbinom(
        1, 1,
        if (p_aeruginosa[j - 1] == 1) 0.78 else expit(-1.35 + 0.45 * burden[j])
      )
    }
  }
  b_cepacia <- rbinom(length(years), 1, expit(-3.35 + 0.55 * burden))
  s_aureus <- rbinom(length(years), 1, expit(0.20 - 0.20 * burden))
  iv_courses <- pmin(rpois(length(years), exp(-1.0 + 0.55 * burden)), 5L)
  home_oxygen <- rbinom(length(years), 1, expit(-4.2 + 1.25 * burden + 0.03 * ages))
  corticosteroids <- rbinom(length(years), 1, expit(-3.7 + 1.05 * burden))
  depression <- rbinom(length(years), 1, expit(-2.8 + 0.018 * ages + 0.20 * burden))

  annual_event_probability <- expit(
    -5.15 + 0.035 * pmax(ages - 20, 0) + 0.035 * (70 - fev1) +
      0.55 * home_oxygen + 0.32 * cfrd + 0.22 * iv_courses + 0.42 * b_cepacia
  )
  event_draw <- runif(length(years)) < annual_event_probability
  event_index <- if (any(event_draw)) which(event_draw)[1] else NA_integer_

  if (is.na(event_index)) {
    last_index <- length(years)
    event_indicator <- 0L
    event_type <- "Censored"
    event_age <- NA_real_
    followup_age <- ages[last_index] + 1
  } else {
    last_index <- event_index
    event_indicator <- 1L
    event_type <- sample(c("Death", "Lung transplant"), 1, prob = c(0.68, 0.32))
    event_age <- ages[event_index] + runif(1, 0.1, 0.95)
    followup_age <- event_age
  }

  encounter <- tibble::tibble(
    patient_id = p$patient_id,
    report_year = years[seq_len(last_index)],
    current_age = ages[seq_len(last_index)],
    fev1_pct_predicted = round(fev1[seq_len(last_index)], 1),
    fvc_pct_predicted = round(fvc[seq_len(last_index)], 1),
    bmi_percentile = round(bmi[seq_len(last_index)], 1),
    cfrd = cfrd[seq_len(last_index)],
    iv_antibiotic_times = iv_courses[seq_len(last_index)],
    home_oxygen = home_oxygen[seq_len(last_index)],
    b_cepacia = b_cepacia[seq_len(last_index)],
    p_aeruginosa = p_aeruginosa[seq_len(last_index)],
    s_aureus = s_aureus[seq_len(last_index)],
    corticosteroids = corticosteroids[seq_len(last_index)],
    depression = depression[seq_len(last_index)]
  )

  summary <- tibble::tibble(
    patient_id = p$patient_id,
    event_indicator = event_indicator,
    event_type = event_type,
    event_age = event_age,
    followup_age = followup_age
  )
  list(encounter = encounter, summary = summary)
}

simulated <- purrr::map(seq_len(n_patients), simulate_patient)
encounter <- purrr::map_dfr(simulated, "encounter")
event_summary <- purrr::map_dfr(simulated, "summary")

# Add modest missingness to exercise the same complete-record rules used for
# the three continuous clinical predictors in the thesis.
set.seed(config$seed + 1L)
for (variable in c("fev1_pct_predicted", "fvc_pct_predicted", "bmi_percentile")) {
  missing_rows <- sample.int(nrow(encounter), size = floor(0.04 * nrow(encounter)))
  encounter[[variable]][missing_rows] <- NA_real_
}

patient <- patient_base |>
  dplyr::select(-patient_effect, -latent_severity) |>
  dplyr::left_join(event_summary, by = "patient_id")

dir.create(here::here("data", "raw"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(patient, here::here("data", "raw", "patient.csv"))
readr::write_csv(encounter, here::here("data", "raw", "annual_encounter.csv"))

message(
  "Synthetic data: ", nrow(patient), " patients, ", nrow(encounter),
  " annual records, ", sum(patient$event_indicator), " composite events."
)
