### This script plots residual risk of transmission
### as a function of time since symptom onset or improvement
### under various infectiousness assumptions

## ===============================================#
## Setup -----------------
## ===============================================#
library(dplyr)
library(readr)
library(ggplot2)
library(tidyselect)
library(yaml)
library(dotenv)
library(devtools)

devtools::load_all("isolation")

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

parameters_config <- yaml::read_yaml(file.path(
  "analysis",
  "parameters.yml"
))

# set seed to be able to replicate randomness in sub-sample for plotting
# (used in development mode; production mode does not sub-sample)

set.seed(parameters_config$overall_parameters$seed)

# create paths for input to this script (output from run_ind_scenarios.R)

overall_compare_guidance_df_path <- file.path(
  model_scenarios_dir,
  products_config$scenarios$overall_compare_guidance_df
)

# create path to bring in underlying trajectory parameters

posterior_path <- file.path(
  results_dir,
  products_config$vl_model$stan_posterior_directory,
  products_config$vl_model$posterior_sampled
)

colors <- c(
  "Log VL" = "#a6bddb",
  "Antigen" = "#fb9a99",
  "Culture" = "#b2df8a",
  "Log VL x Antigen" = "#918adf",
  "Log VL x Culture" = "#dfdc8a"
)

## ===============================================#
## Data -----------------
## ===============================================#

overall_compare_guidance_df <- readr::read_csv(
  file = overall_compare_guidance_df_path
)

posterior_df <- readr::read_csv(
  file = posterior_path
)

# before taking any further steps, we filter category 4 out of the data

overall_compare_guidance_df <- overall_compare_guidance_df |>
  dplyr::filter(symp_type_cat %in% c(1:3))

# restrict to "main" results and the only needed columns

main_symp_cat_df <- overall_compare_guidance_df |>
  dplyr::ungroup() |>
  dplyr::filter(label == "main") |>
  dplyr::select(
    posterior_id, symp_type_cat,
    si_time
  )

# merge in data on the underlying parameters

main_symp_cat_df <- merge(x = main_symp_cat_df, y = posterior_df, all.x = TRUE)

# we want to make a dataframe that is long format
# for each infectiousness transform, we have a time
# since symp onset, and for that time we have the columns
# of mean, median, q10, q25, q75, q90.

calculate_infectiousness_over_time <- function(
    main_symp_cat_df,
    min_day,
    max_day,
    infectious_dt,
    min_day_display,
    max_day_display) {
  time_vector <- seq(
    min_day,
    max_day,
    infectious_dt
  )

  # time_vector is the time range over which infectiousness
  # is calculated, the same as in the main intervention simulations
  # in contrast, display_time_vector is simply the time points
  # for which residual infectiousness is plotted in this figure
  display_time_vector <- seq(min_day_display, max_day_display, 1)

  I_t_summarized <- data.frame()

  for (infect_transform in c(
    "logVL", "antigen", "culture",
    "logVL_x_antigen", "logVL_x_culture"
  )) {
    print(infect_transform)

    I_t_samples <- matrix(
      nrow = length(display_time_vector),
      ncol = nrow(main_symp_cat_df)
    )

    for (traj in seq_len(nrow(main_symp_cat_df))) {
      traj_vars <- main_symp_cat_df[traj, ]
      vl_curve <- isolation::triangle_vl(
        t = time_vector,
        dp = traj_vars$dp,
        tp = traj_vars$tp,
        wp = traj_vars$wp,
        wr = traj_vars$wr
      )

      I_t <- isolation::simulate_intrinsic_infectiousness(
        t = time_vector,
        vl_curve = vl_curve,
        prob_antigen_positive = isolation::prob_antigen_positive(
          logVL = vl_curve,
          antigen_50 = traj_vars$antigen_50,
          sigma_antigen = traj_vars$sigma_antigen
        ),
        prob_culture_positive = isolation::prob_culture_positive(
          t = time_vector,
          tp = traj_vars$tp,
          logVL = vl_curve,
          culture_50 = traj_vars$culture_50,
          culture_beta = traj_vars$culture_beta,
          sigma_culture = traj_vars$sigma_culture
        ),
        infec_threshold = parameters_config$scenario_parameters$infec_threshold,
        transform_fun = infect_transform
      )

      I_t_temp <- data.frame(time = time_vector, I_t = I_t) |>
        # the use of an offset term lets us shift the time series to be
        # wrt to symp improvement as well
        dplyr::mutate(day = ceiling(time - traj_vars$t_offset)) |>
        dplyr::group_by(day) |>
        dplyr::summarise(daily_risk = sum(I_t)) |>
        dplyr::mutate(
          total_risk = sum(daily_risk),
          normalized_daily_risk = daily_risk / total_risk * 100,
          residual_normalized_daily_risk = 100 - cumsum(normalized_daily_risk)
        ) |>
        # Restrict to days in specified range
        dplyr::filter(day >= min_day_display & day <= max_day_display)

      I_t_temp <- I_t_temp |>
        dplyr::select(day, residual_normalized_daily_risk)

      # the following if statements "fill out" values for
      # residual_normalized_daily_risk where day is missing
      # even though we will set min_day_display > min_day
      # and set max_day_display < max_day,
      # we still need this step because t_offset can lead to
      # missing values of day at the extremes
      if (min(I_t_temp$day) > min_day_display) {
        I_t_temp <- dplyr::bind_rows(data.frame(
          "day" = min_day_display:(min(I_t_temp$day) - 1),
          "residual_normalized_daily_risk" = 100
        ), I_t_temp)
      }

      if (max(I_t_temp$day) < max_day_display) {
        I_t_temp <- dplyr::bind_rows(I_t_temp, data.frame(
          "day" = (max(I_t_temp$day) + 1):max_day_display,
          "residual_normalized_daily_risk" = 0
        ))
      }

      stopifnot(identical(display_time_vector, I_t_temp$day))
      stopifnot(I_t_temp$residual_normalized_daily_risk >= -1e-10)

      I_t_samples[, traj] <- I_t_temp$residual_normalized_daily_risk
    }

    I_t_summarized <- as.data.frame(I_t_samples) |>
      # every column is a particular sample of residual_normalized_daily_risk
      # over time (each row is a new timepoint)
      # we want to calculate the box plot stats for each time
      # so hence the rowwise and then the summarise
      dplyr::rowwise() |>
      dplyr::summarise(
        mean = mean(dplyr::c_across(tidyselect::where(is.numeric)),
          na.rm = TRUE
        ),
        middle = quantile(dplyr::c_across(tidyselect::where(is.numeric)),
          probs = .5, na.rm = TRUE
        ),
        q75 = quantile(dplyr::c_across(tidyselect::where(is.numeric)),
          probs = .75, na.rm = TRUE
        ),
        q25 = quantile(dplyr::c_across(tidyselect::where(is.numeric)),
          probs = .25, na.rm = TRUE
        ),
        q90 = quantile(dplyr::c_across(tidyselect::where(is.numeric)),
          probs = .9, na.rm = TRUE
        ),
        q10 = quantile(dplyr::c_across(tidyselect::where(is.numeric)),
          probs = .1, na.rm = TRUE
        )
      ) |>
      dplyr::mutate(
        infect_transform = infect_transform,
        day = display_time_vector
      ) |>
      dplyr::bind_rows(I_t_summarized)
  }

  # Rename for display and reorder so that main analysis (logVL) is first
  I_t_summarized <- I_t_summarized |>
    dplyr::mutate(infect_transform = dplyr::case_match(
      infect_transform,
      "logVL" ~ "Log VL",
      "antigen" ~ "Antigen",
      "culture" ~ "Culture",
      "logVL_x_antigen" ~ "Log VL x Antigen",
      "logVL_x_culture" ~ "Log VL x Culture"
    ))

  I_t_summarized$infect_transform <- factor(I_t_summarized$infect_transform,
    levels = names(colors),
    ordered = TRUE
  )

  return(I_t_summarized)
}

# Plot
make_residual_risk_boxplot <- function(I_t_summarized, xlabel, img_path) {
  boxplot <- ggplot2::ggplot(
    I_t_summarized,
    ggplot2::aes(x = day, y = middle, fill = infect_transform)
  ) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = q10, ymax = q90),
      position = ggplot2::position_dodge(), width = 0.81
    ) +
    ggplot2::geom_crossbar(ggplot2::aes(ymin = q25, ymax = q75),
      position = ggplot2::position_dodge(), width = 0.81
    ) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(
      values = colors,
      name = "VL to infectiousness\ntransform"
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 15)
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(-10, 15, by = 5),
      expand = c(0.01, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 100, by = 20),
      expand = c(0.01, 0)
    ) +
    ggplot2::labs(
      x = xlabel,
      y = "Residual risk of transmission\n(% of transmission remaining)"
    )

  # Save image
  ggplot2::ggsave(img_path, width = 12, height = 7, dpi = 360)
}

base_image_path <- file.path(
  results_dir,
  products_config$products$directory
)

if (!dir.exists(base_image_path)) {
  dir.create(base_image_path, recursive = TRUE)
}

# subsample for development mode, for production use all posterior_ids
if (mode_dependent_output_subdirectory != "production") {
  subsampled_ids <- sample(main_symp_cat_df$posterior_id, size = 100)
} else {
  subsampled_ids <- main_symp_cat_df$posterior_id
}

main_symp_cat_df <- main_symp_cat_df |>
  dplyr::filter(posterior_id %in% subsampled_ids)

# calculate infectiousness wrt to symp onset
main_symp_cat_df$t_offset <- -0.5

I_t_summarized <- calculate_infectiousness_over_time(
  main_symp_cat_df,
  min_day = parameters_config$scenario_parameters$minimum_days_before_symptom_onset_considered, # nolint: line_length_linter.
  max_day = parameters_config$scenario_parameters$max_si +
    parameters_config$scenario_parameters$min_post_iso_days,
  infectious_dt = parameters_config$scenario_parameters$infectiousness_integration_dt, # nolint: line_length_linter.
  min_day_display = -5,
  max_day_display = 15
)

make_residual_risk_boxplot(
  I_t_summarized,
  "Time since symptom onset (days)",
  img_path = file.path(
    base_image_path,
    "figure_supplemental_residual_risk_symptom_onset.png"
  )
)

# calculate infectiousness wrt to symp improvement
# subtract 1.5 days because 1.5 days after observed symptom improvement
# time is when we model symptom improvement as actually happening
# (see, e.g., documentation for simulate_isolation_proposed())
main_symp_cat_df$t_offset <- main_symp_cat_df$si_time - 1.5

I_t_summarized <- calculate_infectiousness_over_time(
  main_symp_cat_df,
  min_day = parameters_config$scenario_parameters$minimum_days_before_symptom_onset_considered, # nolint: line_length_linter.
  max_day = parameters_config$scenario_parameters$max_si +
    parameters_config$scenario_parameters$min_post_iso_days,
  infectious_dt = parameters_config$scenario_parameters$infectiousness_integration_dt, # nolint: line_length_linter.
  min_day_display = -7,
  max_day_display = 8
)

make_residual_risk_boxplot(
  I_t_summarized,
  "Time since symptom improvement (days)",
  file.path(
    base_image_path,
    "figure_supplemental_residual_risk_symptom_improvement.png"
  )
)
