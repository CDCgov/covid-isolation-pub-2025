### This script creates a multi-panel figure on sensitivity analyses
### It is based on the "overall" weighted results, restricted to categories 1-3
### The plotted means are means weighted by transmission potential
### In addition to the figure, the script also outputs a CSV file with key
### summary statistics that correspond to what is plotted in the figure.
### Some of these statistics (specifically, about what happens if we assume
### post-isolation precautions have no effect) are included in the manuscript
### text. Therefore, they need to be output numerically, not just plotted.

## ===============================================#
## Setup -----------------
## ===============================================#
library(dplyr)
library(readr)
library(purrr)
library(yaml)

dotenv::load_dot_env(file = ".env")
args <- commandArgs(trailingOnly = TRUE)

mode_dependent_output_subdirectory <- dplyr::if_else( # nolint: object_length_linter, line_length_linter.
  length(args) > 0,
  args[1],
  "development"
) # implying development is the default

products_config <- yaml::read_yaml("products_config.yml")

results_dir <- file.path(
  Sys.getenv("MAIN_OUTPUT_DIR"),
  mode_dependent_output_subdirectory
)

model_scenarios_dir <- file.path(
  results_dir,
  products_config$scenarios$directory
)


# create paths for input to this script

overall_compare_guidance_df_path <- file.path(
  model_scenarios_dir,
  products_config$scenarios$overall_compare_guidance_df
)

# create paths for the output (the figure file)
fig_path <- file.path(
  results_dir,
  products_config$products$directory
)
sens_analysis_fig_path_pdf <- file.path(
  fig_path,
  "figure_manuscript_sensitivity_analysis.pdf"
)
sens_analysis_fig_path_png <- file.path(
  fig_path,
  "figure_manuscript_sensitivity_analysis.png"
)
if (!dir.exists(fig_path)) {
  dir.create(fig_path, recursive = TRUE)
}

# create paths for the output (the csv file)

sens_analysis_table_path <- file.path(
  fig_path,
  "table_manuscript_sensitivity_analysis.csv"
)

## ===============================================#
## Data -----------------
## ===============================================#

overall_compare_guidance_df <- readr::read_csv(
  file = overall_compare_guidance_df_path
)

# drop category 4

overall_compare_guidance_df <- overall_compare_guidance_df |>
  filter(symp_type_cat != 4)

## ===============================================#
## Extract data, name levels, calculate means
## ===============================================#

iso_eff_df <- overall_compare_guidance_df |>
  dplyr::ungroup() |>
  dplyr::filter(
    label %in% c("main", "iso_eff_1", "iso_eff_0.7")
  ) |>
  dplyr::select(label, diff_averted, intrinsic_I_t) |>
  dplyr::mutate(label = dplyr::case_match(
    label,
    "main" ~ "90%*",
    "iso_eff_1" ~ "100%",
    "iso_eff_0.7" ~ "70%"
  )) |>
  dplyr::mutate(label = factor(label,
    levels = c(
      "70%",
      "90%*",
      "100%"
    ), ordered = T
  ))

iso_eff_means <- iso_eff_df |>
  dplyr::group_by(label) |>
  summarise(wt_mean_diff_averted = weighted.mean(
    x = diff_averted,
    w = intrinsic_I_t, na.rm = T
  ))

post_iso_eff_df <- overall_compare_guidance_df |>
  dplyr::ungroup() |>
  dplyr::filter(
    label %in% c("main", "post_iso_eff_0", "post_iso_eff_0.25")
  ) |>
  dplyr::select(label, diff_averted, intrinsic_I_t) |>
  dplyr::mutate(label = dplyr::case_match(
    label,
    "main" ~ "50%*",
    "post_iso_eff_0" ~ "0%",
    "post_iso_eff_0.25" ~ "25%"
  )) |>
  dplyr::mutate(label = factor(label,
    levels = c(
      "0%",
      "25%",
      "50%*"
    ), ordered = T
  ))

post_iso_eff_means <- post_iso_eff_df |>
  dplyr::group_by(label) |>
  summarise(wt_mean_diff_averted = weighted.mean(
    x = diff_averted,
    w = intrinsic_I_t, na.rm = T
  ))

test_count_df <- overall_compare_guidance_df |>
  dplyr::ungroup() |>
  dplyr::filter(
    label %in% c("main", "test_count_1")
  ) |>
  dplyr::select(label, diff_averted, intrinsic_I_t) |>
  dplyr::mutate(label = dplyr::case_match(
    label,
    "main" ~ "2 tests*",
    "test_count_1" ~ "1 test"
  )) |>
  dplyr::mutate(label = factor(label,
    levels = c(
      "1 test",
      "2 tests*"
    ), ordered = T
  ))

test_count_means <- test_count_df |>
  dplyr::group_by(label) |>
  summarise(wt_mean_diff_averted = weighted.mean(
    x = diff_averted,
    w = intrinsic_I_t, na.rm = T
  ))

iso_lag_df <- overall_compare_guidance_df |>
  dplyr::ungroup() |>
  dplyr::filter(
    label %in% c("main", "iso_lag_2")
  ) |>
  dplyr::select(label, diff_averted, intrinsic_I_t) |>
  dplyr::mutate(label = dplyr::case_match(
    label,
    "main" ~ "0 days*",
    "iso_lag_2" ~ "2 days"
  )) |>
  dplyr::mutate(label = factor(label,
    levels = c(
      "0 days*",
      "2 days"
    ), ordered = T
  ))

iso_lag_means <- iso_lag_df |>
  dplyr::group_by(label) |>
  summarise(wt_mean_diff_averted = weighted.mean(
    x = diff_averted,
    w = intrinsic_I_t, na.rm = T
  ))

infect_transform_df <- overall_compare_guidance_df |>
  dplyr::ungroup() |>
  dplyr::filter(
    label %in% c(
      "main",
      "infect_transform_antigen",
      "infect_transform_culture",
      "infect_transform_logVL_x_antigen",
      "infect_transform_logVL_x_culture"
    )
  ) |>
  dplyr::select(label, diff_averted, intrinsic_I_t) |>
  dplyr::mutate(label = dplyr::case_match(
    label,
    "main" ~ "log(viral load)*",
    "infect_transform_antigen" ~ "Pr(antigen+)",
    "infect_transform_culture" ~ "Pr(culture+)",
    "infect_transform_logVL_x_antigen" ~ "log(VL) x Pr(antigen+)",
    "infect_transform_logVL_x_culture" ~ "log(VL) x Pr(culture+)"
  )) |>
  dplyr::mutate(label = factor(label,
    levels = rev(c(
      "log(viral load)*",
      "Pr(antigen+)",
      "Pr(culture+)",
      "log(VL) x Pr(antigen+)",
      "log(VL) x Pr(culture+)"
    )), ordered = T
  ))

infect_transform_means <- infect_transform_df |>
  dplyr::group_by(label) |>
  summarise(wt_mean_diff_averted = weighted.mean(
    x = diff_averted,
    w = intrinsic_I_t, na.rm = T
  ))

infec_threshold_df <- overall_compare_guidance_df |>
  dplyr::ungroup() |>
  dplyr::filter(
    label %in% c("main", "infec_threshold_0")
  ) |>
  dplyr::select(label, diff_averted, intrinsic_I_t) |>
  dplyr::mutate(label = dplyr::case_match(
    label,
    "main" ~ "1.5 log10(IU)/mL*",
    "infec_threshold_0" ~ "0 log10(IU)/mL"
  )) |>
  dplyr::mutate(label = factor(label,
    levels = rev(c(
      "1.5 log10(IU)/mL*",
      "0 log10(IU)/mL"
    )), ordered = T
  ))

infec_threshold_means <- infec_threshold_df |>
  dplyr::group_by(label) |>
  summarise(wt_mean_diff_averted = weighted.mean(
    x = diff_averted,
    w = intrinsic_I_t, na.rm = T
  ))

## ===============================================#
## Create and export summary stats CSV
## ===============================================#

summary_fx <- function(df) {
  df |>
    group_by(label) |>
    summarise(
      min = min(diff_averted),
      Q1 = quantile(diff_averted, 0.25, na.rm = TRUE),
      median = median(diff_averted, na.rm = TRUE),
      Q2 = quantile(diff_averted, 0.75, na.rm = TRUE),
      max = max(diff_averted),
      weighted_mean = weighted.mean(
        x = diff_averted,
        w = intrinsic_I_t, na.rm = T
      )
    ) |>
    mutate(label = as.character(label))
}

sens_dfs <- list(
  iso_eff = iso_eff_df,
  post_iso_eff = post_iso_eff_df,
  test_count = test_count_df,
  iso_lag = iso_lag_df,
  infect_transform = infect_transform_df,
  infec_threshold = infec_threshold_df
)

sens_summary_stats <- map(sens_dfs, summary_fx) |> bind_rows(.id = "data_frame")

write.csv(
  x = sens_summary_stats, file = sens_analysis_table_path,
  row.names = FALSE
)

## ===============================================#
## Figure -----------------
## ===============================================#

## multi-panel figure of sensitivity analysis results

sens_analysis_fig_fx <- function() {
  par(mfrow = c(6, 1), mar = c(1, 11, 4, 1), oma = c(6, 0, 0, 0), xpd = NA)

  x_min <- -70
  x_max <- 70

  boxplot(
    formula = iso_eff_df$diff_averted ~ iso_eff_df$label,
    horizontal = TRUE, outcex = 0.1, outcol = "grey",
    ylab = NA, xlab = NA, las = 1, ylim = c(x_min, x_max)
  )
  points(
    x = iso_eff_means$wt_mean_diff_averted,
    y = iso_eff_means$label,
    pch = 17
  )
  title(
    main = "A. Proportion of transmission potential averted when in isolation",
    adj = 0
  )

  boxplot(
    formula = post_iso_eff_df$diff_averted ~ post_iso_eff_df$label,
    horizontal = TRUE, outcex = 0.1, outcol = "grey",
    ylab = NA, xlab = NA, las = 1, ylim = c(x_min, x_max)
  )
  points(
    x = post_iso_eff_means$wt_mean_diff_averted,
    y = post_iso_eff_means$label,
    pch = 17
  )
  title(
    main = "B. Proportion of transmission potential averted when in \npost-isolation precautions", # nolint: line_length_linter
    adj = 0
  )

  boxplot(
    formula = test_count_df$diff_averted ~ test_count_df$label,
    horizontal = TRUE, outcex = 0.1, outcol = "grey",
    ylab = NA, xlab = NA, las = 1, ylim = c(x_min, x_max)
  )
  points(
    x = test_count_means$wt_mean_diff_averted,
    y = test_count_means$label,
    pch = 17
  )
  title(main = "C. Number of tests (under previous guidance)", adj = 0)

  boxplot(
    formula = iso_lag_df$diff_averted ~ iso_lag_df$label,
    horizontal = TRUE, outcex = 0.1, outcol = "grey",
    ylab = NA, xlab = NA, las = 1, ylim = c(x_min, x_max)
  )
  points(
    x = iso_lag_means$wt_mean_diff_averted,
    y = iso_lag_means$label,
    pch = 17
  )
  title(main = "D. Days from symptom onset to isolation", adj = 0)

  boxplot(
    formula = infect_transform_df$diff_averted ~ infect_transform_df$label,
    horizontal = TRUE, outcex = 0.1, outcol = "grey",
    ylab = NA, xlab = NA, las = 1, ylim = c(x_min, x_max)
  )
  points(
    x = infect_transform_means$wt_mean_diff_averted,
    y = infect_transform_means$label,
    pch = 17
  )
  title(
    main = "E. Assumption relating viral load to transmission potential",
    adj = 0
  )

  boxplot(
    formula = infec_threshold_df$diff_averted ~ infec_threshold_df$label,
    horizontal = TRUE, outcex = 0.1, outcol = "grey",
    ylab = NA,
    xlab = "Percentange point difference in transmission potential (updated guidance - previous guidance)", # nolint: line_length_linter
    las = 1, ylim = c(x_min, x_max)
  )
  points(
    x = infec_threshold_means$wt_mean_diff_averted,
    y = infec_threshold_means$label,
    pch = 17
  )
  title(
    main = "F. Assumed minimum viral load threshold for transmission",
    adj = 0
  )

  mtext(
    text = expression(paste(italic("Previous guidance more protective"))),
    side = 1, line = 4.5, adj = 0, cex = 0.7
  )
  mtext(
    text = expression(paste(italic("Updated guidance more protective"))),
    side = 1, line = 4.5, adj = 1, cex = 0.7
  )

  par(xpd = NA)
  arrows(x0 = -30, x1 = -50, y0 = -2.8, length = 0.1)
  arrows(x0 = 30, x1 = 50, y0 = -2.8, length = 0.1)
}

pdf(file = sens_analysis_fig_path_pdf, width = 6.5, height = 8)

sens_analysis_fig_fx()

dev.off()

png(
  filename = sens_analysis_fig_path_png,
  width = 6.5, height = 8, units = "in", res = 300
)

sens_analysis_fig_fx()

dev.off()
