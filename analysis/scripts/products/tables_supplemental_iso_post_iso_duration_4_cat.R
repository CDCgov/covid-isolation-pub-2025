### This script generates tables on duration of isolation and post-isolation
### precautions. Category 4 is included in this script.
### Both tables (one based on means, the other based on medians)
### could be used in the supplement.

## ===============================================#
## Setup -----------------
## ===============================================#
library(dplyr)
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

# create paths for input to this script (output from 05_run_ind_scenarios.R)

compare_guidance_df_path <- file.path(
  model_scenarios_dir,
  products_config$scenarios$compare_guidance_df
)

overall_compare_guidance_df_path <- file.path(
  model_scenarios_dir,
  products_config$scenarios$overall_compare_guidance_df
)

# create path for the output (the table file)

duration_iso_post_iso_table_path <- file.path(
  results_dir,
  products_config$products$directory,
  "table_supplemental_iso_post_iso_duration_4_cat_median.xlsx"
)

mean_duration_iso_post_iso_table_path <- file.path(
  results_dir,
  products_config$products$directory,
  "table_supplemental_iso_post_iso_duration_4_cat_mean.xlsx"
)

## ===============================================#
## Functions -----------------
## ===============================================#

# function that creates a string of the median and the IQR
median_iqr_str_fx <- function(x) {
  paste0(
    sprintf(
      "%.1f",
      round(median(x, na.rm = TRUE), digits = 1)
    ),
    " (",
    sprintf(
      "%.1f",
      round(quantile(x, probs = c(0.25), na.rm = TRUE), digits = 1)
    ),
    ", ",
    sprintf(
      "%.1f",
      round(quantile(x, probs = c(0.75), na.rm = TRUE), digits = 1)
    ),
    ")"
  )
}

# function that creates a string of the mean
mean_str_fx <- function(x) {
  paste0(
    sprintf(
      "%.1f",
      round(mean(x, na.rm = TRUE), digits = 1)
    )
  )
}

## ===============================================#
## Data -----------------
## ===============================================#

compare_guidance_df <- readr::read_csv(
  file = compare_guidance_df_path
)

overall_compare_guidance_df <- readr::read_csv(
  file = overall_compare_guidance_df_path
)

# combine the main output dataframe with the "overall" sampled df

overall_compare_guidance_df <- overall_compare_guidance_df |>
  mutate(symp_type_cat = 5)

compare_guidance_df <- dplyr::bind_rows(
  compare_guidance_df, overall_compare_guidance_df
)

# calculate time spent in any precautions

compare_guidance_df <- compare_guidance_df |>
  mutate(
    p_iso_and_post_iso_prop = p_iso_prop + p_post_iso_prop,
    p_iso_and_post_iso_curr = p_iso_curr + p_post_iso_curr,
    diff_iso_and_post_iso = p_iso_and_post_iso_prop - p_iso_and_post_iso_curr
  )

stopifnot(
  all.equal(
    compare_guidance_df$diff_iso_and_post_iso,
    compare_guidance_df$diff_iso + compare_guidance_df$diff_post_iso
  )
)

# change antigen detection probabilities to percentages

compare_guidance_df <- compare_guidance_df |>
  mutate(
    antigen_detect_prob = antigen_detect_prob * 100
  )

# calculate summary statistics

summary_table <- compare_guidance_df |>
  dplyr::group_by(label, symp_type_cat) |>
  dplyr::summarise(
    n = n(),
    across(
      c(
        antigen_detect_prob,
        p_iso_prop, p_iso_curr, diff_iso,
        p_post_iso_prop, p_post_iso_curr, diff_post_iso,
        p_iso_and_post_iso_prop, p_iso_and_post_iso_curr, diff_iso_and_post_iso
      ),
      list(
        "median_iqr_str" = median_iqr_str_fx,
        "mean_str" = mean_str_fx
      )
    )
  )

summary_table <- summary_table |> dplyr::mutate(
  symp_type_cat = dplyr::case_match(
    symp_type_cat,
    1 ~ "1. Shortness of breath",
    2 ~ "2. Fever or body aches",
    3 ~ "3. Mild respiratory symptoms",
    4 ~ "4. Other non-specific symptoms",
    5 ~ "Overall"
  )
)

#####################################
### table based on medians
#####################################

# duration of isolation, post-isolation precautions, and overall table
# this is conceptually three tables, but will be printed with shared headers

dur_iso_table <- summary_table |>
  dplyr::ungroup() |>
  dplyr::filter(label == "main") |>
  dplyr::select(
    symp_type_cat,
    antigen_detect_prob_median_iqr_str,
    p_iso_curr_median_iqr_str,
    p_iso_prop_median_iqr_str,
    diff_iso_median_iqr_str
  )

dur_post_iso_table <- summary_table |>
  dplyr::ungroup() |>
  dplyr::filter(label == "main") |>
  dplyr::select(
    symp_type_cat,
    antigen_detect_prob_median_iqr_str,
    p_post_iso_curr_median_iqr_str,
    p_post_iso_prop_median_iqr_str,
    diff_post_iso_median_iqr_str
  )

dur_iso_and_post_iso_table <- summary_table |>
  dplyr::ungroup() |>
  dplyr::filter(label == "main") |>
  dplyr::select(
    symp_type_cat,
    antigen_detect_prob_median_iqr_str,
    p_iso_and_post_iso_curr_median_iqr_str,
    p_iso_and_post_iso_prop_median_iqr_str,
    diff_iso_and_post_iso_median_iqr_str
  )

# now write to XLSX file
# set up workbook and headers
duration_iso_wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(duration_iso_wb, "Sheet1")
openxlsx::writeData(duration_iso_wb,
  sheet = "Sheet1",
  x =
    t(c(
      "Symptom category",
      "Antigen detection probability %, median (IQR)",
      "Previous guidance",
      "Updated guidance",
      "Difference (updated – previous)"
    )),
  startRow = 1, colNames = FALSE
)

# add portion of table on duration of isolation
openxlsx::mergeCells(duration_iso_wb,
  sheet = "Sheet1",
  cols = 3:5, rows = 2
)
openxlsx::writeData(duration_iso_wb,
  sheet = "Sheet1",
  x =
    t(c(NA, NA, "Expected days in isolation, median (IQR)")),
  startRow = 2, colNames = FALSE
)
openxlsx::writeData(duration_iso_wb,
  sheet = "Sheet1",
  x = dur_iso_table,
  startRow = 3, colNames = FALSE
)

# add portion of table on duration of post-isolation precautions
openxlsx::mergeCells(duration_iso_wb,
  sheet = "Sheet1",
  cols = 3:5, rows = 8
)
openxlsx::writeData(duration_iso_wb,
  sheet = "Sheet1",
  x =
    t(c(NA, NA, "Expected days in post-isolation precautions, median (IQR)")),
  startRow = 8, colNames = FALSE
)
openxlsx::writeData(duration_iso_wb,
  sheet = "Sheet1",
  x = dur_post_iso_table,
  startRow = 9, colNames = FALSE
)

# add portion of table on isolation and post-isolation combined
openxlsx::mergeCells(duration_iso_wb,
  sheet = "Sheet1",
  cols = 3:5, rows = 14
)
openxlsx::writeData(duration_iso_wb,
  sheet = "Sheet1",
  x =
    t(c(NA, NA, "Expected days in isolation and post-isolation precautions combined, median (IQR)")), # nolint: line_length_linter
  startRow = 14, colNames = FALSE
)
openxlsx::writeData(duration_iso_wb,
  sheet = "Sheet1",
  x = dur_iso_and_post_iso_table,
  startRow = 15, colNames = FALSE
)

# add borders to certain cells

openxlsx::addStyle(
  wb = duration_iso_wb, sheet = "Sheet1",
  style = openxlsx::createStyle(border = c("bottom")),
  rows = c(2, 8, 14), cols = 3
)

openxlsx::addStyle(
  wb = duration_iso_wb, sheet = "Sheet1",
  style = openxlsx::createStyle(border = c("bottom")),
  rows = c(1, 7, 13), cols = 1:5, gridExpand = TRUE
)

openxlsx::addStyle(
  wb = duration_iso_wb, sheet = "Sheet1",
  style = openxlsx::createStyle(textDecoration = "bold"),
  rows = c(1, 2, 8, 14), cols = 1:5,
  gridExpand = TRUE, stack = TRUE
)

# save workbook
openxlsx::saveWorkbook(
  wb = duration_iso_wb,
  file = duration_iso_post_iso_table_path,
  overwrite = TRUE
)

#####################################
### table based on means
#####################################

# duration of isolation, post-isolation precautions, and overall table
# this is conceptually three tables, but will be printed with shared headers

mean_dur_iso_table <- summary_table |>
  dplyr::ungroup() |>
  dplyr::filter(label == "main") |>
  dplyr::select(
    symp_type_cat,
    antigen_detect_prob_mean_str,
    p_iso_curr_mean_str,
    p_iso_prop_mean_str,
    diff_iso_mean_str
  )

mean_dur_post_iso_table <- summary_table |>
  dplyr::ungroup() |>
  dplyr::filter(label == "main") |>
  dplyr::select(
    symp_type_cat,
    antigen_detect_prob_mean_str,
    p_post_iso_curr_mean_str,
    p_post_iso_prop_mean_str,
    diff_post_iso_mean_str
  )

mean_dur_iso_and_post_iso_table <- summary_table |>
  dplyr::ungroup() |>
  dplyr::filter(label == "main") |>
  dplyr::select(
    symp_type_cat,
    antigen_detect_prob_mean_str,
    p_iso_and_post_iso_curr_mean_str,
    p_iso_and_post_iso_prop_mean_str,
    diff_iso_and_post_iso_mean_str
  )

# now write to XLSX file
# set up workbook and headers
mean_duration_iso_wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(mean_duration_iso_wb, "Sheet1")
openxlsx::writeData(mean_duration_iso_wb,
  sheet = "Sheet1",
  x =
    t(c(
      "Symptom category",
      "Antigen detection probability %, mean",
      "Previous guidance",
      "Updated guidance",
      "Difference (updated – previous)"
    )),
  startRow = 1, colNames = FALSE
)

# add portion of table on duration of isolation
openxlsx::mergeCells(mean_duration_iso_wb,
  sheet = "Sheet1",
  cols = 3:5, rows = 2
)
openxlsx::writeData(mean_duration_iso_wb,
  sheet = "Sheet1",
  x =
    t(c(NA, NA, "Expected days in isolation, mean")),
  startRow = 2, colNames = FALSE
)
openxlsx::writeData(mean_duration_iso_wb,
  sheet = "Sheet1",
  x = mean_dur_iso_table,
  startRow = 3, colNames = FALSE
)

# add portion of table on duration of post-isolation precautions
openxlsx::mergeCells(mean_duration_iso_wb,
  sheet = "Sheet1",
  cols = 3:5, rows = 8
)
openxlsx::writeData(mean_duration_iso_wb,
  sheet = "Sheet1",
  x =
    t(c(NA, NA, "Expected days in post-isolation precautions, mean")), # nolint: line_length_linter
  startRow = 8, colNames = FALSE
)
openxlsx::writeData(mean_duration_iso_wb,
  sheet = "Sheet1",
  x = mean_dur_post_iso_table,
  startRow = 9, colNames = FALSE
)

# add portion of table on isolation and post-isolation combined
openxlsx::mergeCells(mean_duration_iso_wb,
  sheet = "Sheet1",
  cols = 3:5, rows = 14
)
openxlsx::writeData(mean_duration_iso_wb,
  sheet = "Sheet1",
  x =
    t(c(NA, NA, "Expected days in isolation and post-isolation precautions combined, mean")), # nolint: line_length_linter
  startRow = 14, colNames = FALSE
)
openxlsx::writeData(mean_duration_iso_wb,
  sheet = "Sheet1",
  x = mean_dur_iso_and_post_iso_table,
  startRow = 15, colNames = FALSE
)

# add borders to certain cells

openxlsx::addStyle(
  wb = mean_duration_iso_wb, sheet = "Sheet1",
  style = openxlsx::createStyle(border = c("bottom")),
  rows = c(2, 8, 14), cols = 3
)

openxlsx::addStyle(
  wb = mean_duration_iso_wb, sheet = "Sheet1",
  style = openxlsx::createStyle(border = c("bottom")),
  rows = c(1, 7, 13), cols = 1:5, gridExpand = TRUE
)

openxlsx::addStyle(
  wb = mean_duration_iso_wb, sheet = "Sheet1",
  style = openxlsx::createStyle(textDecoration = "bold"),
  rows = c(1, 2, 8, 14), cols = 1:5,
  gridExpand = TRUE, stack = TRUE
)

# save workbook
openxlsx::saveWorkbook(
  wb = mean_duration_iso_wb,
  file = mean_duration_iso_post_iso_table_path,
  overwrite = TRUE
)
