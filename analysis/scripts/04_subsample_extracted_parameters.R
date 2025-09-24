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

parameters_config <- yaml::read_yaml(file.path(
  "analysis",
  "parameters.yml"
))

products_config <- yaml::read_yaml("products_config.yml")

# *for the purposes of this script* the "MAIN_OUTPUT_DIR" serves as both the
# source of the input (the processed output from the stan model)
# and the home for the output (the sub-sampled parameters)
results_output_dir <- file.path(
  Sys.getenv("MAIN_OUTPUT_DIR"),
  mode_dependent_output_subdirectory
)

# the input here is the "full" extracted output from the Stan model
# at the moment, this is set to be the individual-level extractions

posterior_full_path <- file.path(
  results_output_dir,
  products_config$vl_model$stan_posterior_directory,
  products_config$vl_model$posterior_individual
)

# the output from this script is a dataframe sub-sampled from the Stan extract

posterior_sampled_path <- file.path(
  results_output_dir,
  products_config$vl_model$stan_posterior_directory,
  products_config$vl_model$posterior_sampled
)

posterior_sampled_archive_path <- file.path(
  results_output_dir,
  products_config$vl_model$stan_posterior_directory,
  "archive",
  paste0(
    "VL_symp_cat_sampled_",
    strftime(Sys.Date(), "%y%m%d"), ".csv"
  )
)

set.seed(parameters_config$overall_parameters$seed)

## ===============================================#
## Data -----------------
## ===============================================#

posterior_df_full <- readr::read_csv(file = posterior_full_path)

## ===============================================#
## Sample a subset of the parameter sets ----------
## ===============================================#

# create a unique posterior_id for each row (set of parameters) from the
# extract from stan
# this id is not essential, but could be used to trace back the source of a
# given set of sub-sampled parameters
# one benefit to creating this id is that subsequent scripts in the pipeline
# can ignore the particular id columns that arise in a specific stan extract
# (e.g., there are different id columns for individualized versus pooled)

posterior_df_full <- posterior_df_full |>
  dplyr::mutate(posterior_id = seq_len(nrow(posterior_df_full)))

# select variables that are needed for the intervention simulations
# this ensures that the required variables--and only the required variables--
# are in the output

posterior_df_full <- posterior_df_full |>
  dplyr::select(
    posterior_id, symp_type_cat,
    antigen_50, sigma_antigen,
    culture_50, sigma_culture, culture_beta,
    tp, dp, wp, wr,
    si_time
  )

# now sample from the posterior, the same number of draws per symp_cat

samples_per_cat <- parameters_config$sub_sample_parameters$samples_per_symp_cat

symp_type_cats <- parameters_config$global_data_markers$symptom_categories

# check that symptom categories are as expected
stopifnot(identical(
  symp_type_cats,
  as.integer(unique(posterior_df_full$symp_type_cat))
))

# check that enough samples are available to draw
min_available_samples_per_cat <- min(table(posterior_df_full$symp_type_cat))
stopifnot(min_available_samples_per_cat >= samples_per_cat)

posterior_df_list <- vector(mode = "list", length = length(symp_type_cats))
for (i in seq_along(posterior_df_list)) {
  posterior_df_cat <- posterior_df_full |>
    dplyr::filter(symp_type_cat == symp_type_cats[i])
  posterior_df_list[[i]] <- posterior_df_cat[
    sample(x = seq_len(nrow(posterior_df_cat)), size = samples_per_cat),
  ]
}

## ===============================================#
## Write output -----------------
## ===============================================#
posterior_df <- dplyr::bind_rows(posterior_df_list)

## Export to csv
readr::write_csv(
  x = posterior_df,
  file = posterior_sampled_path
)

readr::write_csv(
  x = posterior_df,
  file = posterior_sampled_archive_path
)
