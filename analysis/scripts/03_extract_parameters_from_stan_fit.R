## ===============================#
## Setup --------------
## ===============================#
library(dplyr)
library(readr)
library(purrr)
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

results_output_dir <- file.path(
  Sys.getenv("MAIN_OUTPUT_DIR"),
  mode_dependent_output_subdirectory
)

set.seed(parameters_config$overall_parameters$seed)

## ===============================#
## Functions---------------
## ===============================#
run_param_extraction_and_save <- function(
    stan_fit_name,
    symp_cat,
    parameters_config,
    extract_symp_improv_params,
    stan_fit_directory,
    stan_posterior_directory) {
  if (!dir.exists(stan_posterior_directory)) {
    dir.create(stan_posterior_directory, recursive = TRUE)
  }

  archive_dir <- file.path(stan_posterior_directory, "archive")

  if (!dir.exists(archive_dir)) {
    dir.create(archive_dir, recursive = TRUE)
  }

  # read stan fit
  stan_path <- file.path(
    stan_fit_directory,
    paste0(stan_fit_name, ".rds")
  )
  vl_stan_fit <- readRDS(stan_path)
  rstan::check_hmc_diagnostics(vl_stan_fit)

  # get output df of VL posteriors from Stan fit
  vl_posterior <- isolation::extract_posteriors_to_df(
    vl_stan_fit = vl_stan_fit,
    parameters_config = parameters_config,
    extract_symp_improv_params = extract_symp_improv_params
  )

  vl_posterior$individual_level_params$symp_type_cat <- symp_cat
  vl_posterior$pooled_params$symp_type_cat <- symp_cat

  vl_samples_path <- file.path(
    stan_posterior_directory,
    paste0(stan_fit_name, ".csv")
  )

  archive_vl_samples_path <- file.path(
    archive_dir,
    paste0(
      stan_fit_name,
      strftime(Sys.Date(), "%y%m%d"), ".csv"
    )
  )

  if (file.exists(archive_vl_samples_path)) {
    file.remove(archive_vl_samples_path)
  }

  readr::write_csv(
    vl_posterior$individual_level_params,
    archive_vl_samples_path
  )
  readr::write_csv(vl_posterior$individual_level_params, vl_samples_path)
  return(vl_posterior)
}

## ===============================#
## export CSV----------------------
## ===============================#

base_name <- products_config$vl_model$run_name

stan_fit_directory <- file.path(
  results_output_dir,
  products_config$vl_model$stan_fit_directory
)

stan_posterior_directory <- file.path(
  results_output_dir,
  products_config$vl_model$stan_posterior_directory
)

posteriors <- purrr::map(
  parameters_config$global_data_markers$symptom_categories,
  function(symp_cat) {
    run_param_extraction_and_save(
      stan_fit_name = paste0(base_name, symp_cat),
      symp_cat = symp_cat,
      parameters_config = parameters_config,
      extract_symp_improv_params = TRUE,
      stan_fit_directory = stan_fit_directory,
      stan_posterior_directory = stan_posterior_directory
    )
  }
)

vl_individual_params <- purrr::map_dfr(
  posteriors,
  function(x) x$individual_level_params
)

vl_pooled_params <- purrr::map_dfr(
  posteriors,
  function(x) x$pooled_params
)

# we exported the extracted parameters from each category's fit to a CSV
# but, for ease for the rest of the pipeline, it is helpful to join
# each of the extracted dataframes from the symptom categories together
# and export to one file that has a constant name and does not change location
# that is what we do here

posterior_individual_level_path <- file.path( # nolint: object_length_linter.
  stan_posterior_directory,
  products_config$vl_model$posterior_individual
)

posterior_pooled_path <- file.path(
  stan_posterior_directory,
  products_config$vl_model$posterior_pooled
)

readr::write_csv(
  vl_individual_params,
  posterior_individual_level_path
)

readr::write_csv(
  vl_pooled_params,
  posterior_pooled_path
)

# check that we wrote the right file to the right place!
posterior_pooled_check <- readr::read_csv(posterior_pooled_path)
posterior_individual_check <- readr::read_csv(posterior_individual_level_path)

# do the dims of the two files not match (one should be wrt individual params)
stopifnot(dim(posterior_individual_check)[1] != dim(posterior_pooled_check)[1])
stopifnot(dim(
  posterior_pooled_check
)[1] == (
  # mcmc chain samples * 4 chains * number of symptom categories
  parameters_config$stan_fit_parameters$MCMC_iterations -
    parameters_config$stan_fit_parameters$warmup) * 4 *
  length(parameters_config$global_data_markers$symptom_categories))
