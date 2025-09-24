### Script that plots association between time to symptom improvement
### and time to viral load clearance

## ===============================================#
## Setup -----------------
## ===============================================#
library(dplyr)
library(readr)
library(ggplot2)
library(yaml)

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

vl_model_dir <- file.path(
  results_dir,
  products_config$vl_model$stan_posterior_directory
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

# create paths for input to this script (output in VL_symp_cat_sampled.csv)

VL_symp_cat_sampled_df_path <- file.path(
  vl_model_dir,
  products_config$vl_model$posterior_sampled
)

# create paths for the output (the figure file)

symp_improv_VL_association_plot_path <- file.path(
  results_dir,
  products_config$products$directory,
  "figure_supplemental_symp_improv_VL_association.png"
)

## ===============================================#
## Data -----------------
## ===============================================#

VL_symp_cat_sampled_df <- readr::read_csv(
  file = VL_symp_cat_sampled_df_path
)

## code labels and generate plot

facet_labels <- c(
  "1" = "A. Shortness of breath",
  "2" = "B. Fever or body aches",
  "3" = "C. Mild respiratory symptoms"
)

png(
  filename = symp_improv_VL_association_plot_path,
  width = 6.5, height = 4, units = "in", res = 300
)

VL_symp_cat_sampled_df |>
  filter(symp_type_cat %in% 1:3) |>
  ggplot2::ggplot(ggplot2::aes(x = factor(si_time), y = wr)) +
  ggplot2::geom_violin(trim = TRUE, scale = "count", fill = "grey") +
  ggplot2::facet_wrap(~symp_type_cat,
    labeller = ggplot2::as_labeller(facet_labels)
  ) +
  ggplot2::labs(
    x = "days from onset to observed symptom improvement",
    y = "days from VL peak to clearance"
  ) +
  ggplot2::scale_x_discrete(
    limits = factor(1:15),
    breaks = seq(1, 15, by = 2)
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 30)) +
  ggplot2::stat_summary(fun = median, cex = 0.3) +
  ggplot2::theme_minimal()

dev.off()
