## ===============================================#
## Setup -----------------
## ===============================================#
library(dplyr)
library(yaml)
library(dotenv)
library(DiagrammeR)
library(webshot)
library(htmltools)

dotenv::load_dot_env(file = ".env")
args <- commandArgs(trailingOnly = TRUE)

mode_dependent_output_subdirectory <- dplyr::if_else( # nolint: object_length_linter, line_length_linter.
  length(args) > 0,
  args[1],
  "development"
) # implying development is the default

data_input_dir <- file.path(
  Sys.getenv("MAIN_OUTPUT_DIR"),
  mode_dependent_output_subdirectory
)
output_dir <- path.expand(
  file.path(
    Sys.getenv("MAIN_OUTPUT_DIR"),
    mode_dependent_output_subdirectory
  )
)

# for the programmatic generation of mermaid diagrams,
# need to specify where to get an html reader
webshot::install_phantomjs(force = TRUE)
# the .env file must also have the following:
# export OPENSSL_CONF=/dev/null

products_config <- yaml::read_yaml("products_config.yml")

read_and_render_mermaid <- function(dataset_name) {
  mermaid_diagram <- paste(readLines(file.path(
    data_input_dir,
    products_config$preprocessed_data$directory,
    paste0(dataset_name, products_config$vl_model$consort_diagram)
  )), collapse = "\n")
  print("NOTE: some warning messages may be returned when running this script.
  These are just warnings of not having an open plotting window
  because mermaid diagrams are automatically printed to the plotting window!")
  DiagrammeR::mermaid(mermaid_diagram) |>
    htmltools::html_print() |>
    webshot::webshot(
      file = file.path(
        output_dir,
        products_config$products$directory,
        paste0("figure_supplemental_consort_diagram_", dataset_name, ".png")
      ),
      cliprect = "viewport", zoom = 2,
      vheight = 400,
      expand = c(-120, 0, 0, 0)
    )
}

if (!dir.exists(file.path(
  output_dir, products_config$products$directory
))) {
  dir.create(
    file.path(
      output_dir, products_config$products$directory
    ),
    recursive = TRUE
  )
}

read_and_render_mermaid("RVTN")
read_and_render_mermaid("INHERENT")
