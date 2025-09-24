### Script that creates figures about difference in transmission potential
### averted and a corresponding summary table.
### The figure and table in this script are supplementary, because category 4
### is included.

## ===============================================#
## Setup -----------------
## ===============================================#
library(dplyr)
library(readr)
library(yaml)
library(openxlsx)
library(colorspace)

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

# get and set seed; randomness in this script only relevant to visualization
# and does not affect any calculations

parameters_config <- yaml::read_yaml(file.path(
  "analysis",
  "parameters.yml"
))

set.seed(parameters_config$overall_parameters$seed)

# create paths for input to this script (output from 05_run_ind_scenarios.R)

compare_guidance_df_path <- file.path(
  model_scenarios_dir,
  products_config$scenarios$compare_guidance_df
)

overall_compare_guidance_df_path <- file.path(
  model_scenarios_dir,
  products_config$scenarios$overall_compare_guidance_df
)

# create paths for the various output (the figure files)

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
supp_symp_cat_w_cat_4_fig_path <- file.path(
  results_dir,
  products_config$products$directory,
  "figure_supplemental_results_diff_averted_4_cat.png"
)

diff_averted_summary_table_w_cat_4_path <- file.path(
  results_dir,
  products_config$products$directory,
  "table_supplemental_results_diff_averted_4_cat.xlsx"
)

## ===============================================#
## Data -----------------
## ===============================================#

compare_guidance_df <- readr::read_csv(
  file = compare_guidance_df_path
)

overall_compare_guidance_df <- readr::read_csv(
  file = overall_compare_guidance_df_path
)

# modify the "overall" df to label the rows as category "5"
# this is solely a convenience that allows us to create a single df for plotting
# also, randomly permute the order of rows so as to reduce visual artifacts
# of plotting order (the rows were ordered by symptom category)

overall_compare_guidance_df <- overall_compare_guidance_df |>
  dplyr::mutate(symp_type_cat_orig = symp_type_cat, symp_type_cat = 5)

overall_compare_guidance_df <-
  overall_compare_guidance_df[
    sample(seq_len(nrow(overall_compare_guidance_df)), replace = FALSE),
  ]

compare_guidance_df <- compare_guidance_df |>
  dplyr::mutate(symp_type_cat_orig = symp_type_cat)

compare_guidance_df <- dplyr::bind_rows(
  compare_guidance_df, overall_compare_guidance_df
)

# process data to restrict to main analysis and set levels

main_symp_cat_df <- compare_guidance_df |>
  dplyr::ungroup() |>
  dplyr::filter(label == "main") |>
  dplyr::select(
    symp_type_cat, symp_type_cat_orig,
    si_time, antigen_detect_prob, diff_averted, intrinsic_I_t
  ) |>
  dplyr::mutate(
    symp_type_cat = dplyr::case_match(
      symp_type_cat,
      1 ~ "A. Shortness of breath",
      2 ~ "B. Fever or body aches",
      3 ~ "C. Mild respiratory symptoms",
      4 ~ "D. Other non-specific symptoms",
      5 ~ "E. Overall"
    )
  ) |>
  dplyr::mutate(symp_type_cat = factor(symp_type_cat,
    levels = c(
      "A. Shortness of breath",
      "B. Fever or body aches",
      "C. Mild respiratory symptoms",
      "D. Other non-specific symptoms",
      "E. Overall"
    ),
    ordered = TRUE
  ))

# calculate mean difference averted by symptom category (and overall)

mean_diff_by_symp_cat <- main_symp_cat_df |>
  dplyr::group_by(symp_type_cat) |>
  dplyr::summarise(mean_diff_averted = mean(diff_averted, na.rm = TRUE))

mean_diff_by_symp_cat_wt_I_t <- main_symp_cat_df |>
  dplyr::group_by(symp_type_cat) |>
  dplyr::summarise(wt_mean_diff_averted = weighted.mean(
    x = diff_averted,
    w = intrinsic_I_t, na.rm = TRUE
  ))


## ===============================================#
## Plot figures -----------------
## ===============================================#

# set up colors for antigen detection probability

grey_blue_solid <- sequential_hcl(
  n = 100, h = 270,
  c1 = 0, c2 = 150, cmax = 150,
  l1 = 50, l2 = 50, alpha = 1, power = 1
)
grey_blue <- sequential_hcl(
  n = 100, h = 270,
  c1 = 0, c2 = 150, cmax = 150,
  l1 = 50, l2 = 50, alpha = 0.07, power = 1
)

main_symp_cat_df <- main_symp_cat_df |>
  dplyr::mutate(ant_det_col_level = floor(antigen_detect_prob * 100) + 1)

### Plot supplemental figure similar to Figure 3 but with cat 4 ###
# differences in transmission potential averted
# points are colored by antigen detection probability
# plot includes histograms above each dot plot
# plotted means are means weighted by total infectiousness

x_min <- -50
x_max <- 50

y_min <- -0.3
y_max <- 4.8

dot_plt_width <- 0.3

png(
  filename = supp_symp_cat_w_cat_4_fig_path,
  width = 6.5, height = 8, units = "in", res = 300
)

par(oma = c(5, 1, 1, 3), xpd = NA, yaxs = "i")

layout(
  mat = matrix(data = c(2, 1), nrow = 1, byrow = TRUE),
  widths = c(18, 1)
)

par(mar = c(10, 1, 9, 0))

image(
  x = 1, y = 0:100,
  z = matrix(data = 1:100, nrow = 1, ncol = 100, byrow = T),
  col = grey_blue_solid[1:100],
  ylab = NA, xlab = NA, axes = F
)
axis(side = 4, cex.axis = 0.7, padj = -1.5)
mtext(
  text = "Probability of detection by antigen testing (%)",
  side = 4, line = 1.5, cex = 0.7
)

par(mar = c(1, 0, 0, 0))

plot(
  x = main_symp_cat_df$diff_averted,
  y = 5 - as.numeric(main_symp_cat_df$symp_type_cat) +
    runif(
      n = length(main_symp_cat_df$symp_type_cat),
      min = -dot_plt_width / 2, max = dot_plt_width / 2
    ),
  col = grey_blue[main_symp_cat_df$ant_det_col_level],
  yaxt = "n", ylab = NA, pch = 16, cex = 0.5, ylim = c(y_min, y_max),
  xlim = c(x_min, x_max),
  xlab = NA
)

main_symp_cat_boxplot_list <- boxplot(
  formula =
    main_symp_cat_df$diff_averted ~ main_symp_cat_df$symp_type_cat,
  at = 4:0, whisklty = 0, staplelty = 0,
  horizontal = TRUE, outline = F, col = NA, xaxt = "n", yaxt = "n",
  ylab = NA, xlab = NA, add = T, boxwex = 0.3
)
points(x = mean_diff_by_symp_cat_wt_I_t$wt_mean_diff_averted, y = 4:0, pch = 17)

text(
  x = x_min, y = 4:0 + 0.6,
  labels = levels(main_symp_cat_df$symp_type_cat), pos = 4, cex = 0.7
)

mtext(
  text = expression(paste(italic("Previous guidance more protective"))),
  side = 1, line = 3.5, adj = 0, cex = 0.7
)
mtext(
  text = expression(paste(italic("Updated guidance more protective"))),
  side = 1, line = 3.5, adj = 1, cex = 0.7
)
mtext(
  text = "Percentage point difference in transmission potential (updated guidance - previous guidance)", # nolint: line_length_linter
  side = 1, line = 2, cex = 0.7
)

par(xpd = NA)
arrows(x0 = -20, x1 = -40, y0 = -1.1, length = 0.1)
arrows(x0 = 20, x1 = 40, y0 = -1.1, length = 0.1)

# before actually adding any histograms, create the counts for each
# so we can get the max count in any bin in any histogram
# so we can appropriately apply the same scale across all plots

hist_list <- vector("list", 5)
max_counts <- vector("integer", 5)

for (cat in 1:5) {
  hist_list[[cat]] <- hist(
    x = main_symp_cat_df$diff_averted[
      as.numeric(main_symp_cat_df$symp_type_cat) == cat
    ],
    breaks = seq(from = -54.5, to = 54.5, by = 1), plot = F
  )

  max_counts[cat] <- max(hist_list[[cat]]$counts)
}

# make a slightly adjusted dot plot width so histogram is slightly above plot

dot_plt_width <- 0.3 + 0.02

plot_ratio <- (max(max_counts) * 1.2) / (1 - dot_plt_width)

for (cat in 1:5) {
  par(new = TRUE)

  y_coord_bottom <- (5 - cat) + dot_plt_width / 2
  y_coord_top <- (5 - cat) + 1 - dot_plt_width / 2

  plot(hist_list[[cat]],
    col = "grey", border = NA, freq = TRUE,
    ylim = c(y_min - y_coord_bottom, y_max - y_coord_bottom) *
      plot_ratio,
    main = NA, yaxt = "n", xaxt = "n",
    xlim = c(x_min, x_max), xlab = NA
  )
}

dev.off()

### Supplemental table on difference averted ###
# make table that reports out the summary statistics plotted in the
# "by symptom cat" figure (IQR, median, mean), as well as the number of
# simulated persons not included
# category 4 is included

median_iqr_diff_str <- paste0(
  sprintf(
    "%.1f",
    round(main_symp_cat_boxplot_list$stats[3, ], digits = 1)
  ),
  " (",
  sprintf(
    "%.1f",
    round(main_symp_cat_boxplot_list$stats[2, ], digits = 1)
  ),
  ", ",
  sprintf(
    "%.1f",
    round(main_symp_cat_boxplot_list$stats[4, ], digits = 1)
  ),
  ")"
)

mean_diff_str <- sprintf(
  "%.1f",
  round(mean_diff_by_symp_cat$mean_diff_averted, digits = 1)
)

wt_mean_diff_str <- sprintf(
  "%.1f",
  round(mean_diff_by_symp_cat_wt_I_t$wt_mean_diff_averted, digits = 1)
)

diff_averted_summary_df <-
  data.frame(
    "symp_type_cat" = mean_diff_by_symp_cat$symp_type_cat,
    "n" = main_symp_cat_boxplot_list$n,
    "mean_diff_averted" = mean_diff_str,
    "wt_mean_diff_averted" = wt_mean_diff_str,
    "median_iqr_diff_averted" = median_iqr_diff_str
  )

diff_averted_summary_wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(diff_averted_summary_wb, "Sheet1")
openxlsx::mergeCells(diff_averted_summary_wb,
  sheet = "Sheet1",
  cols = 3:5, rows = 1
)
openxlsx::writeData(diff_averted_summary_wb,
  sheet = "Sheet1",
  x =
    t(c(
      "Symptom category",
      "n",
      "Difference in transmission potential averted (updated – previous)" # nolint: line_length_linter
    )),
  startRow = 1, colNames = FALSE
)

openxlsx::writeData(diff_averted_summary_wb,
  sheet = "Sheet1",
  x =
    t(c(NA, NA, "mean", "weighted mean", "median (IQR)")),
  startRow = 2, colNames = FALSE
)

openxlsx::writeData(diff_averted_summary_wb,
  sheet = "Sheet1",
  x = diff_averted_summary_df,
  startRow = 3, colNames = FALSE
)

openxlsx::saveWorkbook(
  wb = diff_averted_summary_wb,
  file = diff_averted_summary_table_w_cat_4_path,
  overwrite = TRUE
)
