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
# source of the input (the results of running the intervention model)
# and the source of the output (the results weighted by symptom category)
results_output_dir <- file.path(
  Sys.getenv("MAIN_OUTPUT_DIR"),
  mode_dependent_output_subdirectory
)

# this script takes as input the output from 05_run_ind_scenarios.R
# the number of rows of input is the size of the sub-sample from the posterior
# times the number of analyses conducted per parameter set
# for example, if there are 10,000 posterior draws sub-sampled from each of 4
# symptom categories and 1 main analysis plus 11 sensitivity analyses run
# for each parameter set, then the input for this script will have 480,000 rows

model_scenarios_dir <- file.path(
  results_output_dir,
  products_config$scenarios$directory
)

compare_guidance_path <- file.path(
  model_scenarios_dir,
  products_config$scenarios$compare_guidance_df
)

# another input for this script is weights for symptom categories

symp_cat_wts_path <- file.path(
  results_output_dir,
  products_config$preprocessed_data$directory,
  products_config$scenarios$symp_cat_weights
)

# there is one output from this script:
# a dataframe that is structurally like the output of 05_run_ind_scenarios.R,
# but weighted to create an "overall" category
# the output contains one row per analysis per draw; the columns contain
# data on the amount of infectiousness averted under the each guidance,
# expected time spent in isolation and post-isolation precautions under each
# guidance, and differences thereof
# (there are also columns identifying the symptom category, posterior id,
# the analysis, and the symptom improvement time)

# create paths for output

overall_compare_guidance_path <- file.path(
  model_scenarios_dir,
  products_config$scenarios$overall_compare_guidance_df
)

archive_dir <- file.path(
  model_scenarios_dir,
  "archive",
  paste(products_config$scenarios$run_name,
    strftime(Sys.Date(), "%y%m%d"),
    sep = "_"
  )
)

overall_compare_guidance_df_archive_path <- file.path(
  archive_dir,
  paste(
    products_config$scenarios$run_name,
    strftime(Sys.Date(), "%y%m%d"),
    products_config$scenarios$overall_compare_guidance_df,
    sep = "_"
  )
)

set.seed(parameters_config$overall_parameters$seed)

## ===============================================#
## Data -----------------
## ===============================================#

compare_guidance_df <- readr::read_csv(file = compare_guidance_path)

symp_cat_wts <- readr::read_csv(symp_cat_wts_path)

## ===============================================#
## Apply the sampling -----------------
## ===============================================#

# reduce compare_guidance_df to be a df solely of unique values of
# posterior_id and the corresponding symp_type_cat
# (the de-duplication is necessary because compare_guidance_df is expanded such
# that each posterior_id appears for each (sensitivity) analysis)

id_symp_cat_df <- compare_guidance_df |>
  dplyr::select(posterior_id, symp_type_cat) |>
  dplyr::distinct(posterior_id, .keep_all = TRUE)

# check that there is only one symp_type_cat per posterior_id

stopifnot(identical(
  id_symp_cat_df,
  compare_guidance_df |>
    dplyr::distinct(posterior_id, symp_type_cat)
))

# determine which posterior_ids should be in the sample

sampled_ids_df <- isolation::sample_ids_overall_symp_cat(
  symp_id_df = id_symp_cat_df,
  symp_cat_wts = symp_cat_wts,
  n_overall_sample =
    parameters_config$sub_sample_parameters$samples_per_symp_cat
)

# apply the sample

overall_compare_guidance_df <- compare_guidance_df |>
  dplyr::filter(posterior_id %in% sampled_ids_df$posterior_id)

stopifnot(identical(
  sort(sampled_ids_df$posterior_id),
  sort(unique(overall_compare_guidance_df$posterior_id))
))

## ===============================================#
## Write outputs -----------------
## ===============================================#

## Export to csv
readr::write_csv(
  x = overall_compare_guidance_df,
  file = overall_compare_guidance_path
)

readr::write_csv(
  x = overall_compare_guidance_df,
  file = overall_compare_guidance_df_archive_path
)
