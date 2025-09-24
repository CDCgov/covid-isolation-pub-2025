## ===============================#
## Setup --------------
## ===============================#
library(dplyr)
library(readr)
library(yaml)
library(dotenv)

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

# helper functions to write mean and median to strings
# for the programmatic generation of the table
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

percentage_str_fx <- function(num, denom) {
  paste0(
    num,
    " (",
    sprintf(
      "%.1f",
      round(num / denom * 100, digits = 1)
    ),
    ")"
  )
}

# let us first focus on just table 1
# we need the number of people by symptom category and overall
# we need the time frame for when infections were gathered
# we need the ages -- percent above 50 and percent below 18
# we need the distribution of when people occured
# we need the percent who had received a vaccine dose within last year

extract_variables_per_patient <- function(df, supplemental) {
  if (!supplemental) {
    patient_level <- df |>
      dplyr::filter(symp_type_cat %in% c(1:3))
  } else {
    patient_level <- df
  }
  patient_level <- patient_level |>
    dplyr::mutate(
      months_since_last_covid_dose = dplyr::if_else(
        months_since_last_covid_dose == -1,
        NA,
        months_since_last_covid_dose
      )
    ) |>
    dplyr::group_by(contiguous_id, symp_type_cat) |>
    dplyr::summarise(
      age_cat = dplyr::first(age_cat),
      symp_duration = dplyr::first(symp_duration),
      cdc_period = as.factor(dplyr::first(cdc_period)),
      months_since_last_covid_dose = dplyr::first(months_since_last_covid_dose),
      symp_count = max(symp_count),
      symp_duration_censored = dplyr::first(symp_duration_censored)
    ) |>
    dplyr::ungroup()
  return(patient_level)
}

participant_demographic_stats <- function(
    grouped_by_patient_level) {
  summarized <- grouped_by_patient_level |>
    dplyr::summarise(
      `Participants; no.` = dplyr::n(),
      `Category prevalence; %` = sprintf(
        "%.1f",
        round(100 * dplyr::n() / nrow(grouped_by_patient_level),
          digits = 1
        )
      ),
      # 0 is the label for people below 18
      `Aged 0-17 at enrollment; no. (%)` = percentage_str_fx(sum(age_cat == 0), dplyr::n()),
      # 1 is the label for people 18-50
      `Aged 18-49 at enrollment; no. (%)` = percentage_str_fx(sum(age_cat == 1), dplyr::n()),
      # 2 is the label for people  50-64 -- there are no 65+ people in the study
      `Aged 50-64 at enrollment; no. (%)` = percentage_str_fx(sum(age_cat == 2), dplyr::n()),
      # prop in each time period category
      `Delta: Sep-Dec 17, 2021; no. (%)` = percentage_str_fx(sum(cdc_period == "Delta: Sep-Dec 17, '21"), dplyr::n()),
      `Omicron BA.1/BA.2: Dec 18, 2021-Jun 17, 2022; no. (%)` = percentage_str_fx(sum(cdc_period == "Omicron BA1/BA2: Dec 18, '21-Jun 17, '22"), dplyr::n()),
      `Omicron BA.4/BA.5: Jun 18, 2022-Jan 14, 2023; no. (%)` = percentage_str_fx(sum(cdc_period == "Omicron BA4/5: Jun 18, '22-Jan 14, '23"), dplyr::n()),
      `Omicron XBB and descendants: Jan 15-Oct 31, 2023; no. (%)` = percentage_str_fx(sum(cdc_period == "XBB etc: Jan 15-Oct 31, '23"), dplyr::n()),
      # months since last vaccine dose
      `Months since last vaccine dose; median (IQR)` = median_iqr_str_fx(
        months_since_last_covid_dose
      ),
      `Max COVID symptoms; median (IQR)` = median_iqr_str_fx(
        symp_count[symp_duration_censored == 0]
      ),
      `Symptom duration (days to improvement); median (IQR)` = median_iqr_str_fx(
        symp_duration[symp_duration_censored == 0]
      ),
    )
  return(summarized)
}

# grab the variables we need for table 1
patient_level <- extract_variables_per_patient(joint_dataset,
  supplemental = FALSE
)

patient_level_supplemental <- extract_variables_per_patient(joint_dataset,
  supplemental = TRUE
)

# make the table by category
by_category <- patient_level |> dplyr::group_by(symp_type_cat)
by_category <- participant_demographic_stats(by_category)

by_category_supplemental <- patient_level_supplemental |> dplyr::group_by(symp_type_cat)
by_category_supplemental <- participant_demographic_stats(by_category_supplemental)

# make the table overall
overall <- patient_level |> dplyr::group_by()
overall <- participant_demographic_stats(overall)

overall_supplemental <- patient_level_supplemental |> dplyr::group_by()
overall_supplemental <- participant_demographic_stats(overall_supplemental)

participant_demographics <- as.data.frame(t(by_category[
  ,
  !(colnames(by_category) %in% c("symp_type_cat"))
]))
colnames(participant_demographics) <- c(
  "1. Shortness of breath",
  "2. Fever or body aches",
  "3. Mild respiratory symptoms"
)
overall <- as.data.frame(t(overall))
colnames(overall) <- c("Overall")

participant_demographics <- dplyr::bind_cols(
  participant_demographics,
  overall
)

participant_demographics_supplemental <- as.data.frame(t(by_category_supplemental[
  ,
  !(colnames(by_category) %in% c("symp_type_cat"))
]))
colnames(participant_demographics_supplemental) <- c(
  "1. Shortness of breath",
  "2. Fever or body aches",
  "3. Mild respiratory symptoms",
  "4. Other non-specific symptoms"
)
overall_supplemental <- as.data.frame(t(overall_supplemental))
colnames(overall_supplemental) <- c("Overall")

participant_demographics_supplemental <- dplyr::bind_cols(
  participant_demographics_supplemental,
  overall_supplemental
)

write.csv(participant_demographics,
  file = file.path(
    table_path_out,
    "table_supplemental_participant_demographics_3_cat.csv"
  ),
  row.names = TRUE
)

write.csv(participant_demographics_supplemental,
  file = file.path(
    table_path_out,
    "table_supplemental_participant_demographics_4_cat.csv"
  ),
  row.names = TRUE
)
