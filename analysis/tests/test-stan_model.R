# This test simulates data according to our modeled data
# generating process for the viral load using fixed parameters.
# We fit the Stan model to simulated data using those fixed
# parameters to see if they are recovered.

# In totality, this consists of (a) picking viral load pooled parameters,
# individual-level variability, antigen and culture positivity levels,
# and the strength of the relationship between the viral load clearance
# time and time to symptom improvement, (b) simulating individual-level
# data consistent with these parameters, (c) and evaluating the model fit
# on this data to see if the parameters are recovered.

devtools::load_all("isolation")
dotenv::load_dot_env(file = ".env")

mode_dependent_output_subdirectory <- "test"

parameters_config <- yaml::read_yaml(file.path(
  "analysis",
  "parameters.yml"
))

products_config <- yaml::read_yaml(
  "products_config.yml"
)

data_input_dir <- Sys.getenv("DATA_INPUT_DIR")
results_output_dir <- file.path(
  Sys.getenv("MAIN_OUTPUT_DIR"),
  mode_dependent_output_subdirectory
)

set.seed(parameters_config$overall_parameters$seed)

save_stan_fit <- function(
    vl_stan_fit,
    main_name,
    stan_fit_directory) {
  if (!dir.exists(stan_fit_directory)) {
    dir.create(stan_fit_directory, recursive = TRUE)
  }

  archive_dir <- file.path(stan_fit_directory, "archive")
  if (!dir.exists(archive_dir)) {
    dir.create(archive_dir, recursive = TRUE)
  }

  archive_stan_fit_path <- file.path(
    archive_dir,
    paste0(
      main_name, "_",
      strftime(Sys.Date(), "%y%m%d"), ".rds"
    )
  )

  main_stan_fit_path <- file.path(
    stan_fit_directory,
    paste0(main_name, ".rds")
  )

  saveRDS(
    vl_stan_fit,
    main_stan_fit_path
  )

  if (file.exists(archive_stan_fit_path)) {
    file.remove(archive_stan_fit_path)
  }

  file.copy(main_stan_fit_path, archive_stan_fit_path)
}

print("Creating simulation dataset.")
# We can set this parameters as command-line arguments since a user may
# want to adjust them without editing the code or directly
# from the Makefile.
parser <- optparse::OptionParser()
parser <- optparse::add_option(parser, c("-p", "--pairs"),
  type = "integer", default = 5,
  help = "Number of pairs of people",
  metavar = "number",
  dest = "pairs"
)
parser <- optparse::add_option(parser, c("-s", "--samples"),
  type = "integer", default = 20,
  help = "Number of viral load samples per person",
  metavar = "number",
  dest = "samples"
)
parser <- optparse::add_option(parser, c("-d", "--start-day"),
  type = "integer", default = -7,
  help = "Days with respect to symptom onset at which viral load values start",
  metavar = "number",
  dest = "start_day"
)
args <- optparse::parse_args(parser, commandArgs(trailingOnly = TRUE),
  positional_arguments = 0
)$options
# First, we set the number of VL samples per person.
samples <- args$samples
# Next, we set the number of people. We consider pairs of people
# because we want the individual-level parameters to be balanced
# around the pooled parameters. To do this, for every one person
# we set to have a positive deviation from the pooled parameter,
# we set one person to have a negative deviation.
# For testing purposes, we have set this to the smallest number
# that still gets us a large enough simulated dataset to have
# trustworthy samples.
pairs <- args$pairs
# Finally, set how many days pre symptom onset do we start the
# VL trajectory
start_day <- args$start_day

# Pooled parameters we are trying to recreate
tp <- 1.5
wp <- 5
wr <- 10
dp <- 6.5
sigma <- 0.5
antigen_50 <- 5.5
sigma_antigen <- 0.8
culture_50 <- 4
sigma_culture <- 0.2
culture_beta <- 0.1
si_shape <- 1.7
si_beta_0 <- 1.4
si_beta_wr <- 0.8
# The deviation from the pooled parameters.
spacer <- 0.1
# Parameters for censoring of the viral load
lod <- 1.75
loq_lower <- 3.5
loq_upper <- 10

example_input_data <- dplyr::tibble(
  # pair 1, 2, 3, etc.
  pair_num = rep(seq(1, pairs), each = 2 * samples),
  # whether the person in the pair is a + deviation or a - deviation
  pair_idx = rep(c(rep(1, samples), rep(-1, samples)), times = pairs),
  contiguous_id = rep(seq(1, 2 * pairs),
    each = samples
  ),
  # We pretend that we have VL data on the same days for all people.
  # We deal with censoring (i.e., a day stochastically dropping out) below.
  days_since_symp_onset = rep(
    seq(
      start_day,
      start_day + samples - 1
    ),
    times = 2 * pairs
  )
)

example_input_data <- example_input_data |>
  # Make individual-level VL parameters.
  dplyr::mutate(
    tp_indiv = tp + spacer * pair_num * pair_idx,
    dp_indiv = dp * exp(spacer * pair_num * pair_idx),
    wp_indiv = wp * exp(spacer * pair_num * pair_idx),
    wr_indiv = wr * exp(spacer * pair_num * pair_idx)
  ) |>
  # turn params into VL estimate using triangle + add noise to data
  dplyr::mutate(logVL = isolation::triangle_vl(
    days_since_symp_onset,
    dp_indiv,
    tp_indiv,
    wp_indiv,
    wr_indiv
  ) + rnorm(n = nrow(example_input_data), mean = 0, sd = sigma)) |>
  # Calculate probability of antigen positivity to determine whether
  # antigen test is positive or negative.
  # Technically, antigen positivity is a stochastic function of
  # underlying/triangle VL, but by using logVL (which experiences
  # normally distributed noise) we get the same effect without needing
  # to calculate a binomial variable.
  dplyr::mutate(
    antigen = as.integer(plogis(
      q = logVL,
      location = antigen_50,
      scale = sigma_antigen
    ) > 0.5),
    culture = as.integer(plogis(
      q = logVL,
      # Our model assumes culture positivity is a function of time too.
      location = culture_50 + culture_beta * (days_since_symp_onset - tp_indiv),
      scale = sigma_culture
    ) > 0.5)
  ) |>
  # To improve numerical stability in the Stan model,
  # we use the individual-level variability in clearance time
  # as the predictor of symptom improvement time, rather than the
  # full clearance time. The mean clearance time is just incorporated
  # into the Weibull distribution scale parameters, so this choice
  # doesn't change what we are fitting.
  dplyr::mutate(si_scale = exp(
    si_beta_0 + si_beta_wr * spacer * pair_num * pair_idx
  ))

example_input_data <- example_input_data |>
  # For each participant
  dplyr::group_by(contiguous_id) |>
  # Get their symptom improvement time randomly drawn from the Weibull
  dplyr::mutate(
    symp_duration = isolation::estimate_symptom_improvement(
      # There are multiple rows with the same parameter values for this person
      max(si_scale),
      si_shape,
      max_si = parameters_config$scenario_parameters$max_si
    )
  ) |>
  dplyr::mutate(symp_duration_censored = 0)

# Censoring
example_input_data <- example_input_data |>
  dplyr::mutate(lod = lod) |>
  dplyr::mutate(loq_lower = loq_lower) |>
  # Make an artificial loq_upper for this dataset (just so we have the variable)
  dplyr::mutate(loq_upper = loq_upper) |>
  # VL below 0 is below LOD (negative test), btwn 0 and 3 is below LOQ_lower
  dplyr::mutate(logVL = dplyr::case_when(
    logVL <= lod ~ parameters_config$global_data_markers$lod,
    logVL <= loq_lower & logVL > lod ~ parameters_config$global_data_markers$loq_lower,
    logVL >= loq_upper ~ parameters_config$global_data_markers$loq_upper,
    TRUE ~ logVL
  )) |>
  # We don't have logVL every day for each participant in the real data, so
  # we drop some values out. While we drop out values deterministically to
  # make the simulating easier, it doesn't make a difference for verifying
  # simulation results from Stan since data structure is the same.
  dplyr::mutate(logVL = dplyr::case_when(
    days_since_symp_onset %% 4 == 0 ~
      parameters_config$global_data_markers$skipped_test,
    TRUE ~ logVL
  )) |>
  # Do the same kind of censoring for antigen and culture results
  dplyr::mutate(antigen = dplyr::case_when(
    days_since_symp_onset %% 2 == 0 ~
      parameters_config$global_data_markers$skipped_test,
    TRUE ~ antigen
  )) |>
  dplyr::mutate(culture = dplyr::case_when(
    days_since_symp_onset %% 3 == 0 ~
      parameters_config$global_data_markers$skipped_test,
    TRUE ~ culture
  )) |>
  # We do a similar process for symptom improvement times, but
  # we don't want to drop out any symptom improvement times that
  # are equal to 1 day because then we wouldn't have even known
  # whether the person had symptoms or not.
  dplyr::mutate(symp_duration_censored = dplyr::case_when(
    (contiguous_id %% 8 == 0) & (symp_duration != 1) ~ 1,
    TRUE ~ 0
  )) |>
  dplyr::mutate(symp_duration = dplyr::case_when(
    symp_duration_censored == 1 ~ symp_duration - 1,
    TRUE ~ symp_duration
  ))

# based on the simulated data, create a minimum working example csv file
# that can be made public and used to demonstrate fitting the stan model
# (this output file is not used in this particular test, but can be used
# to demonstrate parts of the analytic pipeline)

example_input_data_mwe <- example_input_data |> dplyr::select(
  days_since_symp_onset, logVL, antigen, culture, symp_duration,
  symp_duration_censored, lod, loq_lower, loq_upper, contiguous_id
)

example_input_data_mwe <- dplyr::bind_rows(
  example_input_data_mwe |> dplyr::mutate(symp_type_cat = 1),
  example_input_data_mwe |> dplyr::mutate(symp_type_cat = 2),
  example_input_data_mwe |> dplyr::mutate(symp_type_cat = 3),
  example_input_data_mwe |> dplyr::mutate(symp_type_cat = 4),
)

new_contiguous_id <- dplyr::consecutive_id(example_input_data_mwe$contiguous_id)

example_input_data_mwe$contiguous_id <- new_contiguous_id

preprocessed_output_dir <- file.path(
  results_output_dir, products_config$preprocessed_data$directory
)

if (!dir.exists(preprocessed_output_dir)) {
  dir.create(preprocessed_output_dir, recursive = TRUE)
}

readr::write_csv(
  example_input_data_mwe,
  file.path(
    preprocessed_output_dir,
    products_config$vl_model$stan_ready_preprocessed_data_output
  )
)

example_input_list <- isolation::generate_stan_input_list(
  data = example_input_data,
  parameters_config = parameters_config
)

symp_model_path <- system.file("stan", "triangle_vl_symp_hazard.stan",
  package = "isolation", mustWork = TRUE
)

print("Fitting Stan model.")
fit_joint_symp <- rstan::stan(
  file = symp_model_path,
  data = example_input_list,
  cores = 4, chains = 4,
  iter = parameters_config$stan_fit_parameters$MCMC_iterations,
  control = list(
    adapt_delta = parameters_config$stan_fit_parameters$adapt_delta,
    max_treedepth = parameters_config$stan_fit_parameters$max_treedepth
  ),
  seed = parameters_config$overall_parameters$seed,
  init_r = 2
)

print("Check HMC sampler diagnostics.")
testthat::test_that("HMC diagnostics are good", {
  testthat::expect_equal(rstan::get_num_divergent(fit_joint_symp), 0)
  testthat::expect_equal(rstan::get_num_max_treedepth(fit_joint_symp), 0)
  testthat::expect_equal(length(rstan::get_low_bfmi_chains(fit_joint_symp)), 0)
})

main_name <- products_config$vl_model$test_name

save_stan_fit(
  vl_stan_fit = fit_joint_symp,
  main_name = main_name,
  stan_fit_directory = file.path(
    results_output_dir,
    products_config$vl_model$stan_fit_directory
  )
)

print("Extracting statistics on posterior samples from Stan model.")
vl_posterior <- isolation::extract_posteriors_to_df(
  vl_stan_fit = fit_joint_symp,
  parameters_config = parameters_config,
  extract_symp_improv_params = TRUE
)

stan_posterior_directory <- file.path(
  results_output_dir,
  products_config$vl_model$stan_posterior_directory
)

if (!dir.exists(stan_posterior_directory)) {
  dir.create(stan_posterior_directory, recursive = TRUE)
}

archive_dir <- file.path(stan_posterior_directory, "archive")

if (!dir.exists(archive_dir)) {
  dir.create(archive_dir, recursive = TRUE)
}

vl_samples_path <- file.path(
  stan_posterior_directory,
  paste0(main_name, ".csv")
)

archive_vl_samples_path <- file.path(
  archive_dir,
  paste0(
    main_name,
    strftime(Sys.Date(), "%y%m%d"), ".csv"
  )
)

if (file.exists(archive_vl_samples_path)) {
  file.remove(archive_vl_samples_path)
}

readr::write_csv(vl_posterior$individual_level_params, archive_vl_samples_path)
readr::write_csv(vl_posterior$individual_level_params, vl_samples_path)

# Do we recover the parameters we used to simulate the
# data we put into the Stan model?
# Print out the percent deviation between the estimated and actual parameters
percent_deviation <- function(estimated, actual) {
  return((estimated - actual) / actual * 100)
}

testthat::test_that("Are pooled triangle viral load parameters approximately recovered?", {
  print("Checking if pooled triangle viral load parameters are approximately recovered.")
  print(paste(
    "peak VL magnitude estimated - actual % deviation:",
    percent_deviation(mean(vl_posterior$pooled_params$dp_mean), dp)
  ))
  testthat::expect_equal(
    mean(vl_posterior$pooled_params$dp_mean), dp,
    # By default, tolerance is a relative fractional deviation unless the values are small
    tol = 0.1
  )
  print(paste(
    "peak VL timing estimated - actual % deviation:",
    percent_deviation(mean(vl_posterior$pooled_params$tp_mean), tp)
  ))
  testthat::expect_equal(
    mean(vl_posterior$pooled_params$tp_mean), tp,
    tol = 0.1
  )
  print(paste(
    "proliferation time estimated - actual % deviation:",
    percent_deviation(mean(vl_posterior$pooled_params$wp_mean), wp)
  ))
  testthat::expect_equal(
    mean(vl_posterior$pooled_params$wp_mean), wp,
    tol = 0.1
  )
  print(paste(
    "clearance time estimated - actual % deviation:",
    percent_deviation(mean(vl_posterior$pooled_params$wr_mean), wr)
  ))
  testthat::expect_equal(
    mean(vl_posterior$pooled_params$wr_mean), wr,
    tol = 0.1
  )
  print(paste(
    "RT-PCR logVL measurement error estimated - actual % deviation:",
    percent_deviation(mean(vl_posterior$pooled_params$sigma), sigma)
  ))
  testthat::expect_equal(
    mean(vl_posterior$pooled_params$sigma), sigma,
    tol = 0.1
  )
})

testthat::test_that("Are antigen and culture positivity parameters approximately recovered?", {
  print("Checking if antigen and culture positivity 50% positivity logVLs are approximately recovered.")
  print(paste(
    "logVL for 50% antigen positivity estimated - actual % deviation:",
    percent_deviation(mean(vl_posterior$pooled_params$antigen_50), antigen_50)
  ))
  testthat::expect_equal(
    mean(vl_posterior$pooled_params$antigen_50), antigen_50,
    tol = 0.1
  )
  print(paste(
    "logVL for 50% culture positivity estimated - actual % deviation:",
    percent_deviation(mean(vl_posterior$pooled_params$culture_50), culture_50)
  ))
  testthat::expect_equal(
    mean(vl_posterior$pooled_params$culture_50), culture_50,
    tol = 0.1
  )
  print("Checking if variance in 50% logVL for antigen/culture positivity has the right sign.")
  # Because our estimates of these quantities absorb uncertainty from dropping data in and out,
  # we can only reliably test that the sign is right rather than trying to put some bound on the
  # estimate accuracy. Moreover, what's important for predicting an individual's antigen positivity
  # is not so much the exact value of this parameter but that we have the overall directionality
  # correct.
  print(paste(
    "estimated variance in antigen 50% log VL:", mean(vl_posterior$pooled_params$sigma_antigen),
    "actual:", sigma_antigen
  ))
  testthat::expect_identical(
    sign(mean(vl_posterior$pooled_params$sigma_antigen)), sign(sigma_antigen)
  )
  print(paste(
    "estimated variance in culture 50% log VL:", mean(vl_posterior$pooled_params$sigma_culture),
    "actual:", sigma_culture
  ))
  testthat::expect_identical(
    sign(mean(vl_posterior$pooled_params$sigma_culture)), sign(sigma_culture)
  )
  print(paste(
    "estimated culture beta (time dependence in culture positivity):",
    mean(vl_posterior$pooled_params$culture_beta), "actual:", culture_beta
  ))
  testthat::expect_identical(
    sign(mean(vl_posterior$pooled_params$culture_beta)), sign(culture_beta)
  )
})

testthat::test_that("Are Weibull symptom improvement time parameters approximately recovered?", {
  print(paste(
    "estimated symp improv Weibull shape:", mean(vl_posterior$pooled_params$si_shape),
    "actual:", si_shape
  ))
  testthat::expect_lte(
    abs(mean(vl_posterior$pooled_params$si_shape) - si_shape), 0.5
  )
  print(paste(
    "estimated base Weibull scale for the symptom improvement distribution:", mean(vl_posterior$pooled_params$si_beta_0_exponentiated),
    "actual:", exp(si_beta_0)
  ))
  testthat::expect_lte(
    abs(mean(vl_posterior$pooled_params$si_beta_0_exponentiated) - exp(si_beta_0)), 0.5
  )
  print(paste(
    "estimated association between clearance time and symptom improvement:",
    mean(vl_posterior$pooled_params$si_beta_wr), "actual:", (si_beta_wr * spacer)
  ))
  testthat::expect_lte(
    abs(mean(vl_posterior$pooled_params$si_beta_wr) - (si_beta_wr * spacer)), 0.5
  )
})
