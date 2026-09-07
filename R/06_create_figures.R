# Reporting figures from synthetic holdout predictions.

source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))

predictions <- readr::read_csv(
  here::here("outputs", "tables", "all_predictions.csv"), show_col_types = FALSE
) |>
  dplyr::filter(set == "Holdout", is.finite(predicted))
performance <- readr::read_csv(
  here::here("outputs", "tables", "performance_summary.csv"), show_col_types = FALSE
)
dca <- readr::read_csv(
  here::here("outputs", "tables", "decision_curve_analysis.csv"), show_col_types = FALSE
)

model_order <- c(
  "Logistic - Forward", "Logistic - Backward", "Logistic - Group lasso",
  "Cox - Forward", "Cox - Backward", "Cox - Group lasso",
  "Landmark Cox - Forward", "Landmark Cox - Backward", "Landmark Cox - Group lasso"
)

calibration_points <- predictions |>
  dplyr::group_by(model) |>
  # Quintiles are more stable than octiles in the small, low-event holdout set.
  dplyr::mutate(bin = dplyr::ntile(predicted, 5)) |>
  dplyr::group_by(model, bin) |>
  dplyr::summarise(
    predicted = mean(predicted), observed = mean(observed), n = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(model = factor(model, levels = model_order))

plot_limit <- max(
  0.10,
  calibration_points$predicted,
  calibration_points$observed,
  na.rm = TRUE
)
# Add padding, then round upward to the next 0.05. Basing the limit on grouped
# calibration points prevents a few raw high-risk predictions from compressing
# every panel into the lower-left corner.
plot_limit <- min(1, ceiling(plot_limit * 1.1 * 20) / 20)

calibration_plot <- ggplot2::ggplot(
  calibration_points, ggplot2::aes(predicted, observed)
) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#B84A39") +
  ggplot2::geom_line(colour = "#245B78", linewidth = 0.7) +
  ggplot2::geom_point(ggplot2::aes(size = n), colour = "#245B78", alpha = 0.85) +
  ggplot2::facet_wrap(~ model, ncol = 3) +
  ggplot2::coord_equal(xlim = c(0, plot_limit), ylim = c(0, plot_limit)) +
  ggplot2::scale_size_continuous(range = c(1.5, 4), guide = "none") +
  ggplot2::labs(
    title = "Calibration on synthetic holdout data",
    subtitle = "Points summarize quintiles of predicted two-year risk; dashed line is ideal calibration",
    x = "Mean predicted risk", y = "Observed event proportion",
    caption = paste(
      "Shared axes reflect the grouped calibration range;",
      "individual predictions may extend beyond the display."
    )
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold")
  )

performance_long <- performance |>
  dplyr::filter(set %in% c("CV mean", "Holdout")) |>
  tidyr::pivot_longer(
    c(c_statistic, brier_score), names_to = "metric", values_to = "value"
  ) |>
  dplyr::mutate(
    model = factor(model, levels = rev(model_order)),
    metric = dplyr::recode(
      metric, c_statistic = "C-statistic", brier_score = "Brier score"
    )
  )

performance_plot <- ggplot2::ggplot(
  performance_long,
  ggplot2::aes(value, model, colour = set, shape = set)
) +
  ggplot2::geom_point(size = 2.5, position = ggplot2::position_dodge(width = 0.5)) +
  ggplot2::facet_wrap(~ metric, scales = "free_x") +
  ggplot2::scale_colour_manual(values = c("CV mean" = "#245B78", "Holdout" = "#D17B28")) +
  ggplot2::labs(
    title = "Synthetic-data model performance",
    subtitle = "These values demonstrate the workflow and do not reproduce thesis estimates",
    x = NULL, y = NULL, colour = NULL, shape = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

dca_long <- dca |>
  tidyr::pivot_longer(
    c(net_benefit_model, net_benefit_treat_all, net_benefit_treat_none),
    names_to = "strategy", values_to = "net_benefit"
  ) |>
  dplyr::mutate(
    model_name = factor(.data$model, levels = model_order),
    strategy = dplyr::recode(
      strategy,
      net_benefit_model = "Model",
      net_benefit_treat_all = "Treat all",
      net_benefit_treat_none = "Treat none"
    )
  )

decision_curve_plot <- ggplot2::ggplot(
  dca_long,
  ggplot2::aes(threshold, net_benefit, colour = strategy, linetype = strategy)
) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::facet_wrap(~ model_name, ncol = 3, scales = "free_y") +
  ggplot2::scale_colour_manual(values = c(
    "Model" = "#245B78", "Treat all" = "#777777", "Treat none" = "#222222"
  )) +
  ggplot2::scale_linetype_manual(values = c(
    "Model" = "solid", "Treat all" = "dashed", "Treat none" = "dotted"
  )) +
  ggplot2::labs(
    title = "Decision curves on synthetic holdout data",
    subtitle = "Exploratory extension; decision-curve analysis was not part of the original thesis",
    x = "Threshold probability", y = "Net benefit", colour = NULL, linetype = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold")
  )

dir.create(here::here("outputs", "figures"), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
  here::here("outputs", "figures", "calibration_holdout.png"),
  calibration_plot, width = 10, height = 9, dpi = 180, bg = "white"
)
ggplot2::ggsave(
  here::here("outputs", "figures", "model_performance.png"),
  performance_plot, width = 10, height = 6.5, dpi = 180, bg = "white"
)
ggplot2::ggsave(
  here::here("outputs", "figures", "decision_curves.png"),
  decision_curve_plot, width = 10, height = 9, dpi = 180, bg = "white"
)

message("Created three synthetic-data figures.")
