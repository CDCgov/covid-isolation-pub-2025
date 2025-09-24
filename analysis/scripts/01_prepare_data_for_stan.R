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

## ===============================#
## Convenience functions for loading data
## ===============================#

# this function reads each dataset and extracts the symptom categories
# of interest for the given fit from it before preprocessing the data for stan
read_and_preprocess_data <- function(
    dataset_path,
    parameters_config,
    symp_cats_to_include,
    dataset_name,
    at_least_one_test,
    mode_dependent_output_subdirectory_ = mode_dependent_output_subdirectory,
    products_config_ = products_config) {
  raw_dataset <- readr::read_csv(dataset_path)
  raw_dataset_category <- raw_dataset |>
    dplyr::filter(symp_type_cat %in% symp_cats_to_include)
  preprocessed_cat <- isolation::stan_preprocess(
    data_to_prep = raw_dataset_category,
    parameters_config = parameters_config,
    dataset_name = dataset_name,
    at_least_one_test = at_least_one_test
  )
  mermaid_output <- file(file.path(
    data_input_dir,
    products_config_$preprocessed_data$directory,
    paste0(dataset_name, products_config_$vl_model$consort_diagram)
  ))
  writeLines(preprocessed_cat$mermaid_diagram, mermaid_output)
  close(mermaid_output)
  return(preprocessed_cat$data_to_prep)
}

# this function does the actual reading in of datasets and then joining
read_preprocess_join_data <- function(
    data_input_dir,
    parameters_config,
    symp_cats_to_include,
    inherent_path,
    rvtn_path,
    stan_trimmed) {
  inherent <- read_and_preprocess_data(
    dataset_path = inherent_path,
    parameters_config = parameters_config,
    symp_cats_to_include = symp_cats_to_include,
    dataset_name = "INHERENT",
    at_least_one_test = stan_trimmed
  )

  rvtn <- read_and_preprocess_data(
    dataset_path = rvtn_path,
    parameters_config = parameters_config,
    symp_cats_to_include = symp_cats_to_include,
    dataset_name = "RVTN",
    at_least_one_test = stan_trimmed
  )

  # join the two data sets
  joint_dataset <- dplyr::bind_rows(inherent, rvtn) |>
    # need to renumber the ids to be contiguous across the two data sets
    dplyr::mutate(contiguous_id = dplyr::consecutive_id(id))
  return(joint_dataset)
}

# first we need to read in the files from each study
# and these data are useful later, so we need to store them

inherent_path <- file.path(
  data_input_dir,
  products_config$preprocessed_data$directory,
  products_config$preprocessed_data$INHERENT
)

rvtn_path <- file.path(
  data_input_dir,
  products_config$preprocessed_data$directory,
  products_config$preprocessed_data$RVTN
)

joint_dataset_stan_trimmed <- read_preprocess_join_data(
  data_input_dir = data_input_dir,
  parameters_config = parameters_config,
  symp_cats_to_include = parameters_config$global_data_markers$symptom_categories, # nolint: line_length_linter.
  inherent_path = inherent_path,
  rvtn_path = rvtn_path,
  stan_trimmed = TRUE
)

# write to CSV so that next script can use it
joint_dataset_stan_trimmed <- readr::write_csv(
  joint_dataset_stan_trimmed,
  file.path(
    data_input_dir,
    products_config$preprocessed_data$directory,
    products_config$vl_model$stan_ready_preprocessed_data_output
  )
)

joint_dataset_full <- read_preprocess_join_data(
  data_input_dir = data_input_dir,
  parameters_config = parameters_config,
  symp_cats_to_include = parameters_config$global_data_markers$symptom_categories, # nolint: line_length_linter.
  inherent_path = inherent_path,
  rvtn_path = rvtn_path,
  stan_trimmed = FALSE
)

# write to CSV so that next script can use it
joint_dataset_full <- readr::write_csv(
  joint_dataset_full,
  file.path(
    data_input_dir,
    products_config$preprocessed_data$directory,
    products_config$vl_model$full_preprocessed_data_output
  )
)

# generate symptom category weights based on contents of joint_dataset
symp_cat_wts <- isolation::calculate_symp_cat_weights(
  joint_dataset = joint_dataset_stan_trimmed
)

# write symptom category weights to csv
readr::write_csv(
  x = symp_cat_wts,
  file = file.path(
    data_input_dir,
    products_config$preprocessed_data$directory,
    products_config$scenarios$symp_cat_weights
  )
)
