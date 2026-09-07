# Run the complete public demonstration from the project root:
#   Sys.setenv(CF_N_PATIENTS = 1000, CF_LASSO_NLAMBDA = 30)
#   source("main.R")

required_packages <- c(
  "broom", "dplyr", "gglasso", "ggplot2", "grpreg", "here", "MASS",
  "pROC", "purrr", "readr", "rlang", "scales", "survival", "tibble", "tidyr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the missing packages before running the pipeline: ",
    paste(missing_packages, collapse = ", "),
    "\nRun: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

steps <- c(
  "R/00_simulate_data.R",
  "R/01_prepare_data.R",
  "R/02_split_data.R",
  "R/03_model_helpers.R",
  "R/04_fit_final_models.R",
  "R/05_evaluate_models.R",
  "R/06_create_figures.R",
  "R/99_smoke_checks.R"
)

for (step in steps) {
  message("\n===== Running ", step, " =====")
  source(here::here(step), local = FALSE)
}

message("\nPipeline complete. See outputs/tables/ and outputs/figures/.")
