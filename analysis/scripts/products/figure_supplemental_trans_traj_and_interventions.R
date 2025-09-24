### This script plots a multi-panel figure that illustrates infectiousness
### trajectories and transmission potential averted by previous guidance
### and by updated guidance
### the sampling is restricted to symptom categories 1-3.

## ===============================================#
## Setup -----------------
## ===============================================#
library(dplyr)
library(readr)
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

if (!dir.exists(file.path(
  results_dir, products_config$products$directory
))) {
  dir.create(
    file.path(
      results_dir, products_config$products$directory
    ),
    recursive = TRUE
  )
}

# get and set seed; randomness in this script determines which simulated
# cases are selected to be plotted, but doesn't affect analysis/calculations

parameters_config <- yaml::read_yaml(file.path(
  "analysis",
  "parameters.yml"
))

set.seed(parameters_config$overall_parameters$seed)

# create paths for input to this script (output from 05_run_ind_scenarios.R)

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

# create paths for the output (the figure file)

trans_traj_and_interventions_plot_path <- file.path(
  results_dir,
  products_config$products$directory,
  "figure_supplemental_trans_traj_and_interventions.png"
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

# restrict to "main" results and only the needed columns

main_symp_cat_df <- overall_compare_guidance_df |>
  dplyr::ungroup() |>
  dplyr::filter(label == "main") |>
  dplyr::select(
    posterior_id, symp_type_cat,
    si_time, antigen_detect_prob, diff_averted, intrinsic_I_t
  )

# merge in data on the underlying parameters

main_symp_cat_df <- merge(x = main_symp_cat_df, y = posterior_df, all.x = TRUE)

# calculate deciles of transmission potential averted and sample by decile

main_symp_cat_df <- main_symp_cat_df |>
  dplyr::mutate(diff_averted_decile = dplyr::ntile(diff_averted, n = 10))

sample_traj <- main_symp_cat_df |>
  dplyr::filter(!is.na(diff_averted_decile)) |>
  dplyr::group_by(diff_averted_decile) |>
  dplyr::sample_n(size = 5)

sample_traj <- sample_traj[order(sample_traj$diff_averted), ]

png(
  filename = trans_traj_and_interventions_plot_path,
  width = 6.5, height = 8, units = "in", res = 300
)

par(
  mfrow = c(10, 5), mar = c(0.1, 0.1, 0.1, 0.1),
  oma = c(3, 3, 4, 3), yaxs = "i"
)

time_vector <- seq(
  from = parameters_config$scenario_parameters$minimum_days_before_symptom_onset_considered, # nolint: line_length_linter.
  to = parameters_config$scenario_parameters$max_si +
    parameters_config$scenario_parameters$min_post_iso_days,
  by = parameters_config$scenario_parameters$infectiousness_integration_dt
)

y_lim_max <-
  parameters_config$scenario_parameters$infectiousness_integration_dt * 10

for (traj in seq_len(nrow(sample_traj))) {
  traj_vars <- sample_traj[traj, ]
  vl_curve <- isolation::triangle_vl(
    t = time_vector,
    dp = traj_vars$dp,
    tp = traj_vars$tp,
    wp = traj_vars$wp,
    wr = traj_vars$wr
  )

  prob_ant_pos <- isolation::prob_antigen_positive(
    logVL = vl_curve,
    antigen_50 = traj_vars$antigen_50,
    sigma_antigen = traj_vars$sigma_antigen
  )

  I_t <- isolation::simulate_intrinsic_infectiousness(
    t = time_vector,
    vl_curve = vl_curve,
    infec_threshold = parameters_config$scenario_parameters$infec_threshold,
    transform_fun = "logVL"
  )

  curr_list <- isolation::simulate_isolation_current(
    It = I_t,
    t = time_vector,
    prob_antigen_positive = prob_ant_pos,
    symp_type_cat = traj_vars$symp_type_cat,
    symp_improv = traj_vars$si_time,
    iso_lag = parameters_config$scenario_parameters$iso_lag,
    test_count = parameters_config$scenario_parameters$test_count,
    test_interval = parameters_config$scenario_parameters$test_interval,
    min_iso_days_mild = parameters_config$scenario_parameters$min_iso_days_mild,
    min_iso_days_moderate =
      parameters_config$scenario_parameters$min_iso_days_moderate,
    min_total_precaution_days =
      parameters_config$scenario_parameters$min_total_precaution_days,
    iso_eff = parameters_config$scenario_parameters$iso_eff,
    post_iso_eff = parameters_config$scenario_parameters$post_iso_eff,
    symp_improv_adjust = 1.5,
    iso_lag_adjust = 0.5
  )

  prop_list <- isolation::simulate_isolation_proposed(
    It = I_t,
    t = time_vector,
    symp_improv_in = traj_vars$si_time,
    symp_improv_adjust = 1.5,
    iso_lag = parameters_config$scenario_parameters$iso_lag,
    iso_lag_adjust = 0.5,
    min_post_iso_days = parameters_config$scenario_parameters$min_post_iso_days,
    iso_eff = parameters_config$scenario_parameters$iso_eff,
    post_iso_eff = parameters_config$scenario_parameters$post_iso_eff
  )

  plot(time_vector, I_t,
    ylim = c(0, 1) * y_lim_max,
    type = "l", xlab = NA, ylab = NA, yaxt = "n", xaxt = "n"
  )

  curr_averted_vec <- I_t - curr_list$I_t
  curr_averted_times <- time_vector[curr_averted_vec > 0]
  curr_averted_pos <- curr_averted_vec[curr_averted_vec > 0]

  polygon(
    x = c(curr_averted_times, rev(curr_averted_times)),
    y = c(curr_averted_pos, rep(0, times = length(curr_averted_times))),
    col = "green", border = NA
  )

  prop_averted_vec <- I_t - prop_list$I_t
  prop_averted_times <- time_vector[prop_averted_vec > 0]
  prop_averted_pos <- prop_averted_vec[prop_averted_vec > 0]

  polygon(
    x = c(prop_averted_times, rev(prop_averted_times)),
    y = c(prop_averted_pos, rep(0, times = length(prop_averted_times))),
    col = "orange", border = NA
  )

  either_averted_vec <- pmin(curr_averted_vec, prop_averted_vec)
  either_averted_times <- time_vector[either_averted_vec > 0]
  either_averted_pos <- either_averted_vec[either_averted_vec > 0]

  polygon(
    x = c(either_averted_times, rev(either_averted_times)),
    y = c(
      either_averted_pos,
      rep(0, times = length(either_averted_times))
    ),
    col = "lightgrey", border = NA
  )

  recalc_diff_averted <-
    ((sum(prop_averted_vec) - sum(curr_averted_vec)) / sum(I_t)) * 100

  stopifnot(all.equal(traj_vars$diff_averted, recalc_diff_averted))

  arrows(
    x0 = -0.5, y0 = 0.85 * y_lim_max, x1 = traj_vars$si_time - 1.5 + 0.1,
    angle = 90, length = 0.03, lwd = 2, code = 3
  )

  diff_averted_str <- paste(
    sprintf("%.1f", round(traj_vars$diff_averted, digits = 1)), "pp"
  )

  text(
    x = 33, y = 0.8 * y_lim_max, labels = bquote(bold(.(diff_averted_str))),
    adj = 1, cex = 0.7
  )

  ant_det_prob_str <- paste0(
    sprintf(
      "%.1f",
      round(traj_vars$antigen_detect_prob * 100, digits = 1)
    ), "%"
  )

  text(
    x = 33, y = 0.6 * y_lim_max, labels = bquote(italic(.(ant_det_prob_str))),
    adj = 1, cex = 0.7
  )

  text(
    x = 33, y = 0.4 * y_lim_max,
    labels = paste("cat.", traj_vars$symp_type_cat),
    adj = 1, cex = 0.7
  )

  if (traj %in% 46:50) {
    axis(side = 1)
  }

  if (traj == 3) {
    legend("top",
      legend = c(
        "not averted",
        "averted under previous and updated guidance",
        "averted under previous guidance only",
        "averted under updated guidance only"
      ),
      fill = c(NA, "lightgrey", "green", "orange"),
      border = c("black", NA, NA, NA),
      xpd = NA, bty = "n", inset = -0.7, ncol = 2,
      title = "status of transmission potential"
    )
  }
}

mtext(
  text = "time (days from end of day of symptom onset)",
  side = 1, line = 2, cex = 0.7, outer = TRUE
)

mtext(
  text = "transmission potential",
  side = 2, line = 1, cex = 0.7, outer = TRUE
)

mtext(
  text = expression(paste(italic("previous guidance more protective"))),
  side = 4, line = 1, adj = 1, cex = 0.7, outer = TRUE
)
mtext(
  text = expression(paste(italic("updated guidance more protective"))),
  side = 4, line = 1, adj = 0, cex = 0.7, outer = TRUE
)

dev.off()
