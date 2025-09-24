#' Extract Stan output posteriors to data frames
#'
#' extract the pooled posterior samples from the stan fit
#' and output to a dataframe
#'
#' @param vl_stan_fit viral load stan fit
#' @param parameters_config parameters config object
#' @param extract_symp_improv_params whether to extract
#' symptom parameters from stan fit object
#' (i.e., weibull proportional hazards coefs)
#'
#' @return dataframe of individual-level parameters +
#' dataframe of pooled parameters
#' @export
extract_posteriors_to_df <- function(
    vl_stan_fit,
    parameters_config,
    extract_symp_improv_params = TRUE) {
  midpoints <- parameters_config$prior_midpoint_values

  individual_level_params <- data.frame()
  # the number of individuals in the symptom category
  # there are some variables that are individual-specific
  # for instance "wr_raw" -- use the number of times that
  # occurs to automatically figure out number of individuals
  # 1 for chain 1 (just for arguments sake to extract what we need)
  n_indiv <- vl_stan_fit@sim$dims_oi$wr_raw
  for (i in 1:n_indiv) {
    specific_individual_params <- data.frame(
      tp = c(rstan::extract(vl_stan_fit,
        pars = paste0("tp[", i, "]"),
        permuted = FALSE
      )),
      dp = c(rstan::extract(vl_stan_fit,
        pars = paste0("dp[", i, "]"),
        permuted = FALSE
      )),
      wp = c(rstan::extract(vl_stan_fit,
        pars = paste0("wp[", i, "]"),
        permuted = FALSE
      )),
      wr = c(rstan::extract(vl_stan_fit,
        pars = paste0("wr[", i, "]"),
        permuted = FALSE
      ))
    ) |>
      dplyr::mutate(
        study_participant_number_in_category = i,
        # can only assign posterior draw id the row number because
        # permuted is FALSE
        posterior_draw_id = dplyr::row_number()
      )

    if (extract_symp_improv_params) {
      specific_individual_params <- specific_individual_params |>
        dplyr::mutate(si_time = c(rstan::extract(vl_stan_fit,
          pars = paste0("si_simulated[", i, "]"),
          permuted = FALSE
        )))
    }

    stopifnot(sum(is.na(individual_level_params)) == 0)

    individual_level_params <- individual_level_params |>
      dplyr::bind_rows(specific_individual_params)
  }

  pooled_params <- data.frame(
    # need the four hinge parameters: tp, dp, wp, wr
    # also need their stds
    # need antigen_50 and sigma and culture_50, beta, and sigma

    # tp needs no math transformation
    tp_mean = c(rstan::extract(vl_stan_fit,
      pars = "tp_mean",
      permuted = FALSE
    )),
    dp_mean = exp(c(rstan::extract(
      vl_stan_fit,
      pars = "log_dp_mean",
      permuted = FALSE
    ))) * midpoints$dp,
    wp_mean = exp(c(rstan::extract(
      vl_stan_fit,
      pars = "log_wp_mean",
      permuted = FALSE
    ))) * midpoints$wp,
    wr_mean = exp(c(rstan::extract(
      vl_stan_fit,
      pars = "log_wr_mean",
      permuted = FALSE
    ))) * midpoints$wr,
    sigma = c(rstan::extract(
      vl_stan_fit,
      pars = "sigma",
      permuted = FALSE
    )),

    # now the stds
    tp_std = c(rstan::extract(vl_stan_fit,
      pars = "tp_std",
      permuted = FALSE
    )),
    dp_std = exp(c(rstan::extract(
      vl_stan_fit,
      pars = "log_dp_sd",
      permuted = FALSE
    ))),
    wp_std = exp(c(rstan::extract(
      vl_stan_fit,
      pars = "log_wp_sd",
      permuted = FALSE
    ))),
    wr_std = exp(c(rstan::extract(
      vl_stan_fit,
      pars = "log_wr_sd",
      permuted = FALSE
    ))),

    # antigen_50
    antigen_50 = c(rstan::extract(vl_stan_fit,
      pars = "antigen_50",
      permuted = FALSE
    )),
    sigma_antigen = c(rstan::extract(vl_stan_fit,
      pars = "sigma_antigen",
      permuted = FALSE
    )),

    # culture values
    culture_50 = c(rstan::extract(vl_stan_fit,
      pars = "culture_50",
      permuted = FALSE
    )),
    culture_beta = c(rstan::extract(vl_stan_fit,
      pars = "culture_beta",
      permuted = FALSE
    )),
    sigma_culture = c(rstan::extract(vl_stan_fit,
      pars = "sigma_culture",
      permuted = FALSE
    ))
  ) |>
    dplyr::mutate(posterior_draw_id = dplyr::row_number())

  stopifnot(sum(is.na(pooled_params)) == 0)

  if (extract_symp_improv_params) {
    pooled_symp_params <- data.frame(
      si_shape = c(rstan::extract(vl_stan_fit,
        pars = "si_shape",
        permuted = FALSE
      )),
      si_beta_0_exponentiated = exp(c(rstan::extract(
        vl_stan_fit,
        pars = "si_beta_0",
        permuted = FALSE
      ))),
      si_beta_wr = c(rstan::extract(
        vl_stan_fit,
        pars = "si_beta_wr",
        permuted = FALSE
      ))
    ) |>
      dplyr::mutate(posterior_draw_id = dplyr::row_number())

    stopifnot(sum(is.na(pooled_symp_params)) == 0)
    stopifnot(dim(pooled_symp_params)[1] == dim(pooled_params)[1])

    # join overall VL pooled params with symp params
    pooled_params <- pooled_params |>
      dplyr::left_join(pooled_symp_params, by = c("posterior_draw_id"))
  }

  individual_level_params <- pooled_params |>
    dplyr::select(
      posterior_draw_id,
      antigen_50, sigma_antigen,
      culture_50, sigma_culture, culture_beta
    ) |>
    dplyr::right_join(individual_level_params, by = c("posterior_draw_id"))

  stopifnot(dim(individual_level_params)[1] / dim(pooled_params)[1] == n_indiv)
  stopifnot(dim(pooled_params)[1] == vl_stan_fit@sim$chains * (parameters_config$stan_fit_parameters$MCMC_iterations - parameters_config$stan_fit_parameters$warmup)) # nolint: line_length_linter.

  # each id should show up n_indiv number of times because ids are numbers
  # for a particular posterior sample
  stopifnot(all(sapply(
    unique(individual_level_params$posterior_draw_id),
    function(x) sum(individual_level_params$posterior_draw_id == x) == n_indiv
  )))

  return(list(
    individual_level_params = individual_level_params,
    pooled_params = pooled_params
  ))
}
