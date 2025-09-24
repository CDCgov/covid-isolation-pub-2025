## ===============================================#
## Setup -----------------
## ===============================================#
library(dplyr)
library(readr)
library(tibble)
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

parameters_config <- yaml::read_yaml(file.path(
  "analysis",
  "parameters.yml"
))

products_config <- yaml::read_yaml("products_config.yml")

# *for the purposes of this script* the "MAIN_OUTPUT_DIR" serves as both the
# source of the input (the processed output from the stan model)
# and the home for the output (the results of running the intervention model)
results_output_dir <- file.path(
  Sys.getenv("MAIN_OUTPUT_DIR"),
  mode_dependent_output_subdirectory
)

# this script takes as input draws from the posterior
# each row is a set of parameters that, collectively, are sufficient to define
# a viral load trajectory; this script is agnostic to how they were generated
# (e.g., whether they reflect individual variability or variability in the mean)
# this input is the output from subsample_extracted_parameters.R

posterior_path <- file.path(
  results_output_dir,
  products_config$vl_model$stan_posterior_directory,
  products_config$vl_model$posterior_sampled
)

# there are two types of outputs from this script:
# 1) a dataframe of the intervention parameters that are used
# 2) a summary dataframe with differences between guidance
# all of these outputs will live in a directory called model_scenarios

model_scenarios_dir <- file.path(
  results_output_dir,
  products_config$scenarios$directory
)

if (!dir.exists(model_scenarios_dir)) {
  dir.create(model_scenarios_dir, recursive = TRUE)
}

# also need archive directory within model_scenarios
# and a directory within that

archive_dir <- file.path(
  model_scenarios_dir,
  "archive",
  paste(products_config$scenarios$run_name,
    strftime(Sys.Date(), "%y%m%d"),
    sep = "_"
  )
)

if (!dir.exists(archive_dir)) {
  dir.create(archive_dir, recursive = TRUE)
}

# create paths for output

intervention_params_path <- file.path(
  model_scenarios_dir,
  products_config$scenarios$intervention_params
)

compare_guidance_df_path <- file.path(
  model_scenarios_dir,
  products_config$scenarios$compare_guidance_df
)

intervention_params_archive_path <- file.path(
  archive_dir,
  paste(
    products_config$scenarios$run_name,
    strftime(Sys.Date(), "%y%m%d"),
    products_config$scenarios$intervention_params,
    sep = "_"
  )
)

compare_guidance_df_archive_path <- file.path(
  archive_dir,
  paste(
    products_config$scenarios$run_name,
    strftime(Sys.Date(), "%y%m%d"),
    products_config$scenarios$compare_guidance_df,
    sep = "_"
  )
)

## ===============================================#
## Functions -----------------
## ===============================================#
estimate_isolation_impact <- function(params,
                                      posterior_vars,
                                      time_vector,
                                      vl_curve,
                                      prob_antigen_positive,
                                      prob_culture_positive) {
  ## 1) Convert to infectiousness
  nat_hist_ts <- isolation::simulate_intrinsic_infectiousness(
    t = time_vector,
    vl_curve = vl_curve,
    prob_antigen_positive = prob_antigen_positive,
    prob_culture_positive = prob_culture_positive,
    infec_threshold = params$infec_threshold,
    transform_fun = params$infect_transform
  )

  ## 2.1) Apply proposed guidance
  prop_list <- isolation::simulate_isolation_proposed(
    It = nat_hist_ts,
    t = time_vector,
    symp_improv_in = posterior_vars$si_time,
    symp_improv_adjust = 1.5,
    iso_lag = params$iso_lag,
    iso_lag_adjust = 0.5,
    min_post_iso_days = params$min_post_iso_days,
    iso_eff = params$iso_eff,
    post_iso_eff = params$post_iso_eff
  )

  ## 2.2) Apply current guidance
  curr_list <- isolation::simulate_isolation_current(
    It = nat_hist_ts,
    t = time_vector,
    prob_antigen_positive = prob_antigen_positive,
    symp_type_cat = posterior_vars$symp_type_cat,
    symp_improv = posterior_vars$si_time,
    iso_lag = params$iso_lag,
    test_count = params$test_count,
    test_interval = params$test_interval,
    min_iso_days_mild = params$min_iso_days_mild,
    min_iso_days_moderate = params$min_iso_days_moderate,
    min_total_precaution_days = params$min_total_precaution_days,
    iso_eff = params$iso_eff,
    post_iso_eff = params$post_iso_eff,
    symp_improv_adjust = 1.5,
    iso_lag_adjust = 0.5
  )

  delta_t <- params$infectiousness_integration_dt
  p_iso_prop <- sum(prop_list$p_iso, na.rm = TRUE) * delta_t
  p_iso_curr <- sum(curr_list$p_iso, na.rm = TRUE) * delta_t
  p_post_iso_prop <- sum(prop_list$p_post_iso, na.rm = TRUE) * delta_t
  p_post_iso_curr <- sum(curr_list$p_post_iso, na.rm = TRUE) * delta_t

  ## 3) Store properties and compute differences in efficacy
  return(list(
    antigen_detect_prob = curr_list$antigen_detect_prob,
    prop_averted = prop_list$averted,
    curr_averted = curr_list$averted,
    p_iso_prop = p_iso_prop,
    p_iso_curr = p_iso_curr,
    p_post_iso_prop = p_post_iso_prop,
    p_post_iso_curr = p_post_iso_curr,
    diff_averted = prop_list$averted - curr_list$averted,
    diff_iso = p_iso_prop - p_iso_curr,
    diff_post_iso = p_post_iso_prop - p_post_iso_curr,
    t = time_vector,
    nat_hist_It = nat_hist_ts,
    curr_It = curr_list$I_t,
    prop_It = prop_list$I_t
  ))
}

## ===============================================#
## Data -----------------
## ===============================================#

posterior_df <- readr::read_csv(file = posterior_path)

## ===============================================#
## Parameters -----------------
## ===============================================#
# start a dataframe of parameters, the first row of which will be
# the "baseline" conditions for the main analysis
# note that the first two parameters are not about defining the interventions
# per se but are part of defining disease natural history
# ("max_si" and "infec_threshold")

params_df <- tibble::as_tibble(parameters_config[["scenario_parameters"]])

sens_list <- parameters_config$sens_list

for (i in seq_along(sens_list)) {
  if (length(sens_list[[i]]) > 0) {
    params_tmp <- params_df[rep(1, length(sens_list[[i]])), ]
    params_tmp[, names(sens_list[i])] <- sens_list[[i]]
    params_tmp$analysis <- names(sens_list[i])
    params_tmp$label <- sprintf(
      "%s_%s",
      names(sens_list[i]),
      sens_list[[i]]
    )
    params_df <- dplyr::bind_rows(params_df, params_tmp)
  }
}

# write params_df to a csv file
readr::write_csv(
  x = params_df,
  file = intervention_params_path
)
readr::write_csv(
  x = params_df,
  file = intervention_params_archive_path
)

## ===============================================#
## Simulate guidance per person -----------------
## ===============================================#

## Define vector of times t over which to simulate infectiousness
## and interventions; t is relative to symptom onset (more precisely, t = 0
## is the end of the 24-hour period when symptoms started)

## The time vector depends on minimum_days_before_symptom_onset_considered,
## max_si, min_post_iso_days, and infectiousness_integration_dt.
## the code does not support doing sensitivity analyses on these values

stopifnot(
  parameters_config$scenario_parameters$minimum_days_before_symptom_onset_considered == # nolint: line_length_linter.
    unique(params_df$minimum_days_before_symptom_onset_considered),
  parameters_config$scenario_parameters$max_si ==
    unique(params_df$max_si),
  parameters_config$scenario_parameters$min_post_iso_days ==
    unique(params_df$min_post_iso_days),
  parameters_config$scenario_parameters$infectiousness_integration_dt ==
    unique(params_df$infectiousness_integration_dt)
)

time_vector <- seq(
  from = parameters_config$scenario_parameters$minimum_days_before_symptom_onset_considered, # nolint: line_length_linter.
  to = parameters_config$scenario_parameters$max_si +
    parameters_config$scenario_parameters$min_post_iso_days,
  by = parameters_config$scenario_parameters$infectiousness_integration_dt
)

## Create an empty list with length equal to the number of rows of input
## Each row of the input will yield a dataframe with one row of results
## for the main analysis and followed by a row for each sensitivity analysis.
## The motivation for looping over parameters and then looping over analyses
## within a set of parameters is that much of the computational burden
## (e.g., making the time-series of viral load, culture, and antigen curves)
## can be readily shared across sensitivity analyses.
## In theory, it would be possible to share other parts of the computational
## burden (e.g., time in isolation) across some--but not other--sensitivity
## analyses; as a practical matter, the loss in clarity and flexibility is
## likely not worth the gain in efficiency.

compare_guidance_dfs_list <- vector(
  mode = "list",
  length = nrow(posterior_df)
)

for (ss in seq_len(nrow(posterior_df))) {
  posterior_row <- posterior_df[ss, ]

  scenario_df <- posterior_row |>
    dplyr::select(
      posterior_id,
      symp_type_cat,
      si_time
    ) |>
    dplyr::mutate(
      prop_averted = NA, curr_averted = NA,
      p_iso_prop = NA, p_iso_curr = NA,
      p_post_iso_prop = NA, p_post_iso_curr = NA,
      diff_averted = NA,
      diff_iso = NA, diff_post_iso = NA,
      intrinsic_I_t = NA,
      antigen_detect_prob = NA
    )

  scenario_df <- params_df |>
    dplyr::select(label, analysis) |>
    dplyr::bind_cols(scenario_df)

  # calculate logVL trajectory

  vl_timeseries <- isolation::triangle_vl(
    t = time_vector,
    dp = posterior_row$dp,
    tp = posterior_row$tp,
    wp = posterior_row$wp,
    wr = posterior_row$wr
  )

  # calculate antigen positivity trajectory

  prob_antigen_positive <- plogis(
    q = vl_timeseries,
    location = posterior_row$antigen_50,
    scale = posterior_row$sigma_antigen
  )

  # calculate culture positivity trajectory

  prob_culture_positive <- plogis(
    q = vl_timeseries,
    location = posterior_row$culture_50 +
      posterior_row$culture_beta * (time_vector - posterior_row$tp),
    scale = posterior_row$sigma_culture
  )

  ## Loop through each scenario (sensitivity analysis)
  ## Estimate infectiousness, simulate guidance, store results in df
  for (ii in seq_len(nrow(scenario_df))) {
    impact_list <- estimate_isolation_impact(
      params = params_df[ii, ],
      posterior_vars = scenario_df[ii, ],
      time_vector = time_vector,
      vl_curve = vl_timeseries,
      prob_antigen_positive = prob_antigen_positive,
      prob_culture_positive = prob_culture_positive
    )

    scenario_df$antigen_detect_prob[ii] <- impact_list$antigen_detect_prob
    scenario_df$intrinsic_I_t[ii] <- sum(impact_list$nat_hist_It, na.rm = TRUE)
    scenario_df$prop_averted[ii] <- impact_list$prop_averted
    scenario_df$curr_averted[ii] <- impact_list$curr_averted
    scenario_df$p_iso_prop[ii] <- impact_list$p_iso_prop
    scenario_df$p_iso_curr[ii] <- impact_list$p_iso_curr
    scenario_df$p_post_iso_prop[ii] <- impact_list$p_post_iso_prop
    scenario_df$p_post_iso_curr[ii] <- impact_list$p_post_iso_curr
    scenario_df$diff_averted[ii] <- impact_list$diff_averted
    scenario_df$diff_iso[ii] <- impact_list$diff_iso
    scenario_df$diff_post_iso[ii] <- impact_list$diff_post_iso
  }

  compare_guidance_dfs_list[[ss]] <- scenario_df
}

## ===============================================#
## Write outputs -----------------
## ===============================================#
compare_guidance_df <- dplyr::bind_rows(compare_guidance_dfs_list)

## Export to csv
readr::write_csv(
  x = compare_guidance_df,
  file = compare_guidance_df_path
)

readr::write_csv(
  x = compare_guidance_df,
  file = compare_guidance_df_archive_path
)
