### This script generates tables on probability of being antigen or culture
### positive, stratified by VL and day from symptom onset

## ===============================================#
## Setup -----------------
## ===============================================#
library(dplyr)
library(tidyr)
library(readr)
library(yaml)
library(openxlsx)

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
output_dir <- file.path(
  Sys.getenv("MAIN_OUTPUT_DIR"),
  mode_dependent_output_subdirectory
)

table_path_out <- file.path(
  output_dir,
  products_config$products$directory
)

if (!dir.exists(table_path_out)) {
  dir.create(table_path_out, recursive = TRUE)
}

# read the preprocessed data file from CSV
joint_dataset <- readr::read_csv(
  file.path(
    data_input_dir,
    products_config$preprocessed_data$directory,
    products_config$vl_model$full_preprocessed_data_output
  )
)

# create path for the output (the table files)

culture_by_day_by_VL_table_path <- file.path(
  table_path_out,
  "table_supplemental_culture_by_day_by_VL.xlsx"
)

antigen_by_day_by_VL_table_path <- file.path(
  table_path_out,
  "table_supplemental_antigen_by_day_by_VL.xlsx"
)

###########################

# restrict to main analysis symptom categories (exclude category 4)
joint_dataset <- joint_dataset |>
  filter(symp_type_cat != 4)

# exclude participant-days from which no PCR / logVL result is available
logVL_dataset <- joint_dataset |>
  filter(logVL != parameters_config$global_data_markers$skipped_test)

# code categories of days since symptom onset
logVL_dataset <- logVL_dataset |>
  mutate(
    days_since_symp_onset_level =
      cut(
        x = days_since_symp_onset,
        breaks = c(-10, 0, 2, 4, 6, 8, 100),
        labels = c("<=0", "1-2", "3-4", "5-6", "7-8", ">=9")
      )
  )

# code categories of PCR result / log viral load
logVL_dataset <- logVL_dataset |>
  mutate(
    log_VL_level =
      cut(
        x = logVL,
        breaks = c(
          parameters_config$global_data_markers$lod,
          parameters_config$global_data_markers$loq_lower,
          3, 4, 5, 6, 100
        ), right = F,
        labels = c(
          "PCR-", "PCR+, VL < LOQ", "[3, 4)", "[4, 5)", "[5, 6)",
          ">=6"
        )
      )
  )

# function that makes table for "assay_type" either "culture" or "antigen"

tabulate_by_assay_result_fx <- function(assay_type) {
  # make table of antigen or culture results stratified by
  # days symptom from and viral load / PCR result

  table_df <- logVL_dataset |>
    filter({{ assay_type }} !=
      parameters_config$global_data_markers$skipped_test) |>
    group_by(days_since_symp_onset_level, log_VL_level) |>
    summarise(
      n_pos = sum({{ assay_type }} == 1),
      n_tested = n()
    ) |>
    mutate(prop_pos = n_pos / n_tested) |>
    mutate(pos_str = paste0(
      round(prop_pos * 100, digits = 0), "%", " (",
      n_pos, "/", n_tested, ")"
    )) |>
    select(days_since_symp_onset_level, log_VL_level, pos_str) |>
    pivot_wider(names_from = log_VL_level, values_from = pos_str)

  # create a table as a formatted XLSX file

  table_wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(table_wb, "Sheet1")
  openxlsx::mergeCells(table_wb,
    sheet = "Sheet1",
    cols = 2:7, rows = 1
  )
  openxlsx::writeData(table_wb,
    sheet = "Sheet1",
    x = t(c(NA, "Viral load category")),
    startRow = 1, colNames = FALSE
  )
  openxlsx::writeData(table_wb,
    sheet = "Sheet1",
    x = t(c(
      "Days since symptom onset",
      colnames(table_df)[2:7]
    )),
    startRow = 2, colNames = FALSE
  )
  openxlsx::writeData(table_wb,
    sheet = "Sheet1",
    x = table_df,
    startRow = 3, colNames = FALSE
  )

  # add borders and bold style to header cells

  openxlsx::addStyle(
    wb = table_wb, sheet = "Sheet1",
    style = openxlsx::createStyle(
      border = "bottom",
      textDecoration = "bold"
    ),
    rows = 1, cols = 2:7
  )

  openxlsx::addStyle(
    wb = table_wb, sheet = "Sheet1",
    style = openxlsx::createStyle(
      border = "bottom",
      textDecoration = "bold"
    ),
    rows = 2, cols = 1:7
  )

  return(table_wb)
}

# run function to generate tables

culture_wb <- tabulate_by_assay_result_fx(assay_type = culture)

antigen_wb <- tabulate_by_assay_result_fx(assay_type = antigen)

# save workbooks

openxlsx::saveWorkbook(
  wb = culture_wb, file = culture_by_day_by_VL_table_path, overwrite = TRUE
)

openxlsx::saveWorkbook(
  wb = antigen_wb, file = antigen_by_day_by_VL_table_path, overwrite = TRUE
)
