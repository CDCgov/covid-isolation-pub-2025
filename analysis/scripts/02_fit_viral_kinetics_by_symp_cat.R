## ===============================#
## Prep data from INHERENT/RVTN for
## Stan and run with Stan VL/symp model

## ===============================#
## Setup --------------
## ===============================#
library(dplyr)
library(readr)
library(yaml)
library(dotenv)
library(devtools)
library(rstan)

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

data_input_dir <- file.path(
  Sys.getenv("MAIN_OUTPUT_DIR"),
  mode_dependent_output_subdirectory
)

results_output_dir <- file.path(
  Sys.getenv("MAIN_OUTPUT_DIR"),
  mode_dependent_output_subdirectory
)

set.seed(parameters_config$overall_parameters$seed)

default_iter <- parameters_config$stan_fit_parameters$MCMC_iterations
default_warmup <- parameters_config$stan_fit_parameters$warmup
default_control <- list(
  adapt_delta = parameters_config$stan_fit_parameters$adapt_delta,
  max_treedepth = parameters_config$stan_fit_parameters$max_treedepth
)

## ===============================#
## Convenience functions for saving stan fits
## ===============================#

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

## ===============================#
## Function for fitting symp + VL model
## ===============================#

## guarantee model starts sampling from some point on the prior
## remember that init points in MCMC do not impact stationary distribution
## because a markov chain is formed
## instead, these changes just help the chain to converge faster
## note that by initing in different points, we still allow for
## the chains to explore different space in case the posterior
## is multimodal (bad) or Rhat wants to diverge (bad)
init_function <- function(stan_list) {
  std_shrinker <- 0.1
  init_vals <- list(
    # tp
    tp_mean = rnorm(
      1,
      stan_list$tp_prior[1],
      std_shrinker * sqrt(stan_list$tp_prior[2])
    ),
    tp_std = abs(rnorm(
      1,
      stan_list$tp_std_prior[1],
      std_shrinker * sqrt(stan_list$tp_std_prior[2])
    )),
    tp_raw = rnorm(stan_list$n_id, 0, std_shrinker),
    # dp
    log_dp_mean = rnorm(1, 0, std_shrinker * sqrt(stan_list$prior_sd)),
    log_dp_sd = abs(rnorm(1, 0, std_shrinker * sqrt(stan_list$prior_sd))),
    dp_raw = rnorm(stan_list$n_id, 0, std_shrinker),
    ## the next two params are times on a log scale
    ## so we can have a diffuse prior, but we should
    ## initialize over a closer range to ensure measurable logp
    # wp
    log_wp_mean = rnorm(1, 0, std_shrinker * sqrt(stan_list$prior_sd)),
    log_wp_sd = abs(rnorm(1, 1, std_shrinker * sqrt(stan_list$prior_sd))),
    wp_raw = rnorm(stan_list$n_id, 0, std_shrinker),
    # wr
    log_wr_mean = rnorm(1, 0, std_shrinker * sqrt(stan_list$prior_sd)),
    log_wr_sd = abs(rnorm(1, 1, std_shrinker * sqrt(stan_list$prior_sd))),
    wr_raw = rnorm(stan_list$n_id, 0, std_shrinker),
    ## for the rest of these params, we have put a diffuse prior on them
    ## but sometimes that prior is too diffuse to provide measurable logp
    ## over the whole support
    ## init from different but limited values known to have density
    sigma = abs(rnorm(
      1,
      3,
      std_shrinker * sqrt(stan_list$sigma_prior[2])
    )),
    antigen_50 = rnorm(
      1,
      stan_list$antigen_50_prior[1],
      std_shrinker * sqrt(stan_list$antigen_50_prior[2])
    ),
    sigma_antigen = abs(rnorm(
      1,
      2,
      std_shrinker * sqrt(stan_list$sigma_prior[2])
    )),
    culture_50 = rnorm(
      1,
      stan_list$culture_50_prior[1],
      std_shrinker * sqrt(stan_list$culture_50_prior[2])
    ),
    sigma_culture = abs(rnorm(
      1,
      2,
      std_shrinker * sqrt(stan_list$sigma_prior[2])
    )),
    culture_beta = rnorm(
      1,
      stan_list$culture_beta_prior[1],
      std_shrinker * sqrt(stan_list$culture_beta_prior[2])
    ),
    si_shape = 1,
    si_beta_0 = rnorm(
      1,
      stan_list$si_beta_0_prior[1],
      0
    ),
    si_beta_wr = rnorm(
      1,
      stan_list$si_beta_wr_prior[1],
      0
    )
  )

  return(init_vals)
}

# this function runs the stan viral kinetics + symptom duration model
# to each symptom category to fit the parameters and then export those fits
# for later postprocessing
run_vl_symp_duration_model <- function(
    preprocessed_data,
    main_name,
    parameters_config,
    data_input_dir,
    stan_fit_directory,
    qc_figure_directory,
    stan_posterior_directory,
    iter = default_iter,
    warmup = default_warmup,
    control = default_control,
    seed = parameters_config$overall_parameters$seed) {
  # make stan ready lists for joint data
  if (iter <= warmup) {
    stop("MCMC iterations is less than or equal to warmup period.")
  }
  stan_list <- isolation::generate_stan_input_list(
    data = preprocessed_data,
    parameters_config = parameters_config
  )

  print(stan_list$n_id)

  symp_model_path <- system.file("stan", "triangle_vl_symp_hazard.stan",
    package = "isolation", mustWork = TRUE
  )

  fit_joint_symp <- rstan::stan(
    file = symp_model_path,
    data = stan_list,
    cores = 4, chains = 4,
    iter = default_iter,
    warmup = default_warmup,
    control = default_control,
    seed = seed,
    init = lapply(
      1:4, # 4 chains, never run with fewer
      function(id) init_function(stan_list)
    )
  )

  save_stan_fit(
    vl_stan_fit = fit_joint_symp,
    main_name = main_name,
    stan_fit_directory = stan_fit_directory
  )

  if (!dir.exists(qc_figure_directory)) {
    dir.create(qc_figure_directory, recursive = TRUE)
  }

  if (!dir.exists(stan_posterior_directory)) {
    dir.create(stan_posterior_directory, recursive = TRUE)
  }

  isolation::vl_fits_to_pdf(
    vl_stan_fit = fit_joint_symp,
    data_list = stan_list,
    parameters_config = parameters_config,
    main_name = main_name,
    figure_dir = qc_figure_directory,
    data_output_dir = stan_posterior_directory
  )

  return(fit_joint_symp)
}

## ===============================#
## Run model --------------
## ===============================#

# read from CSV
joint_dataset <- readr::read_csv(
  file.path(
    data_input_dir,
    products_config$preprocessed_data$directory,
    products_config$vl_model$stan_ready_preprocessed_data_output
  )
)

# now we fit the model by symptom category to the data
main_name <- products_config$vl_model$run_name

# provide more/less limited initialization values
# to ensure the sampler can start and speed convergence
# the point of an MCMC warmup is that the init values
# do not change the downstream results/samples

for (symp_cat in parameters_config$global_data_markers$symptom_categories) {
  joint_dataset_specific_cat <- joint_dataset |>
    dplyr::filter(symp_type_cat == symp_cat) |>
    dplyr::mutate(contiguous_id = dplyr::consecutive_id(contiguous_id))

  if (mode_dependent_output_subdirectory != "production") {
    joint_dataset_specific_cat <- joint_dataset_specific_cat |>
      dplyr::filter(contiguous_id %in% sample(
        1:max(contiguous_id),
        1 + ceiling(
          parameters_config$overall_parameters$development_dataset_subset *
            max(contiguous_id)
        )
      )) |>
      dplyr::mutate(contiguous_id = dplyr::consecutive_id(contiguous_id))
  }

  fit_joint_symp <- run_vl_symp_duration_model(
    preprocessed_data = joint_dataset_specific_cat,
    main_name = paste0(main_name, symp_cat),
    parameters_config = parameters_config,
    data_input_dir = data_input_dir,
    stan_fit_directory = file.path(
      results_output_dir,
      products_config$vl_model$stan_fit_directory
    ),
    qc_figure_directory = file.path(
      "quality_control"
    ),
    stan_posterior_directory = file.path(
      results_output_dir,
      products_config$vl_model$stan_posterior_directory
    )
  )
}
