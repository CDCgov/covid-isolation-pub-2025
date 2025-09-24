library(testthat)
library(readr)

dotenv::load_dot_env(file = ".env")
args <- commandArgs(trailingOnly = TRUE)

mode_dependent_output_subdirectory <- dplyr::if_else( # nolint: object_length_linter, line_length_linter.
  length(args) > 0,
  args[1],
  "development"
) # development is the default

products_config <- yaml::read_yaml("products_config.yml")

# Define folder path and file names
folder_path <- file.path(
  Sys.getenv("MAIN_OUTPUT_DIR"),
  mode_dependent_output_subdirectory,
  products_config$preprocessed_data$directory
)
file_names <- c(
  products_config$preprocessed_data$RVTN,
  products_config$preprocessed_data$INHERENT
)

# Create the full file paths
data_paths <- file.path(folder_path, file_names)

for (data_path in data_paths) {
  test_that(sprintf("All columns in %s
  are numeric and have no missing values", data_path), {
    message("Checking file: ", data_path) # nolint: indentation_linter
    data <- read_csv(data_path)

    for (col_name in names(data)) { # nolint: indentation_linter
      # Check if the column is numeric
      expect_true(is.numeric(data[[col_name]]),
        info = paste(col_name, "is not numeric.")
      )

      # Check for missing values #nolint: indentation_linter
      expect_true(all(!is.na(data[[col_name]])),
        info = paste(col_name, "contains missing values.")
      )
    }
  })
}
