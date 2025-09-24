#' Plot viral load fits
#'
#' Plot individual-level VL fits to a PDF and save underlying data; these plots
#' are used to check model performance by visual inspection
#'
#' @param vl_stan_fit fitted stan fit object
#' @param data_list underlying data used for `vl_stan_fit`
#' @param parameters_config config file
#' @param main_name name under which output data will be saved
#' @param figure_dir directory where PDF fit file will be saved
#' @param data_output_dir directory where data output will be saved
#'
#' @return saves PDF plots and data frames of observed data and corresponding
#' fitted viral load trajectories
#' @export
vl_fits_to_pdf <- function(
    vl_stan_fit,
    data_list,
    parameters_config,
    main_name,
    figure_dir,
    data_output_dir) {
  # create fig dir
  if (!dir.exists(figure_dir)) {
    dir.create(figure_dir)
  }
  # setup
  grey_col <- rgb(red = 0, green = 0, blue = 0, alpha = 0.1)
  blue_col <- rgb(red = 0, green = 0, blue = 1, alpha = 0.4)
  ts <- seq(-10, 21, 0.1)
  pdf(file.path(figure_dir, paste0(main_name, ".pdf")),
    width = 9, height = 9
  )
  layout(mat = matrix(1:25, ncol = 5, nrow = 5))
  par(mar = c(3, 2, 1, 1), oma = c(1, 2, 1, 2))

  # need to export these data for better plotting in Python
  data_df <- data.frame(
    id = integer(),
    time = double(),
    logVL = double(),
    antigen = integer(),
    culture = integer()
  )

  simulated_df <- data.frame(
    id = integer(),
    time = double(),
    logVL = double(),
    sample_num = integer()
  )
  global_data_markers <- parameters_config$global_data_markers
  for (i in 1:data_list$n_id) {
    # plot fitted parameter estimates for dp[i], tp[i]
    # wp[i], and wr[i] (not pooled)
    # first plot the individual's data
    logVL <- data_list$y[data_list$id == i] # nolint: object_name_linter.
    antigen <- data_list$antigen_result[data_list$id == i]
    days_since_symp_onset <- data_list$t[data_list$id == i]
    culture <- data_list$culture_result[data_list$id == i]
    plot(
      x = days_since_symp_onset,
      # if logVL is missing, plot as -1
      y = ifelse(logVL == global_data_markers$skipped_test, -1,
        ifelse(
          logVL == global_data_markers$lod |
            logVL == global_data_markers$loq_lower,
          # if logVL is negative or below LOQ, plot as 0
          0, logVL # nolint: object_name_linter.
        )
      ),
      pch = ifelse(antigen == 1, 3, # + for positive, x for negative, o for NA
        ifelse(antigen == 0, 4, 1)
      ), frame = FALSE,
      xlim = c(-10, 21),
      ylim = c(-2, 12),
      cex = 0.7
    )
    temp_data <- data.frame(
      id = rep(i, length(antigen)),
      time = days_since_symp_onset,
      logVL = ifelse(logVL == global_data_markers$skipped_test, -1,
        ifelse(
          logVL == global_data_markers$lod,
          0, ifelse(logVL == global_data_markers$loq_lower,
            1, logVL
          )
        )
      ),
      antigen = antigen,
      culture = culture
    )
    data_df <- rbind(data_df, temp_data)
    mtext("days since symptom onset", 1, cex = 0.5, line = 2)
    mtext("-log10 VL", 2, cex = 0.5, line = 2)
    abline(h = 0, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.4))
    abline(
      h = -1, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.4),
      lty = 2
    )
    counter <- 1
    for (chain in 1:4) {
      chain_samples <- vl_stan_fit@sim$samples[[chain]]
      size_template <- chain_samples[[paste0("dp[", i, "]")]]
      for (index in seq(0.75 * length(size_template),
        length(size_template),
        length.out = 15
      )) {
        # plot the VL samples as lines over the points
        dp <- (chain_samples[[paste0("dp[", i, "]")]][index])
        tp <- (chain_samples[[paste0("tp[", i, "]")]][index])
        wp <- (chain_samples[[paste0("wp[", i, "]")]][index])
        wr <- (chain_samples[[paste0("wr[", i, "]")]][index])
        lines(
          x = ts,
          y = sapply(
            X = ts,
            FUN = isolation::triangle_vl,
            dp = dp,
            tp = tp,
            wp = wp,
            wr = wr
          ),
          col = blue_col
        )
        temp_simulated <- data.frame(
          id = rep(i, length(ts)),
          time = ts,
          logVL = sapply(
            X = ts,
            FUN = isolation::triangle_vl,
            dp = dp,
            tp = tp,
            wp = wp,
            wr = wr
          ),
          sample_num = rep(counter, length(ts))
        )
        simulated_df <- rbind(temp_simulated, simulated_df)
        counter <- counter + 1
      }
    }
    ## plot observed viral load values on top to make sure that they are visible
    points(
      x = days_since_symp_onset,
      y = ifelse(logVL == global_data_markers$skipped_test, -1,
        ifelse(
          logVL == global_data_markers$lod |
            logVL == global_data_markers$loq_lower,
          0, logVL
        )
      ),
      pch = ifelse(antigen == 1, 3,
        ifelse(antigen == 0, 4, 1)
      )
    )
  }
  dev.off()
  readr::write_csv(
    data_df,
    file.path(
      data_output_dir,
      paste0("data_", main_name, ".csv")
    )
  )
  readr::write_csv(
    simulated_df,
    file.path(
      data_output_dir,
      paste0("samples_", main_name, ".csv")
    )
  )
}
