#' Make input list for Stan
#'
#' stan expects a list of variables -- turn a dataframe
#' that is ready for stan into a list which also has all the
#' parameters config information
#'
#' @param data output of stan_preprocess
#' @param parameters_config parameters_config.yml
#'
#' @return data list ready for stan input
#' @export
generate_stan_input_list <- function(data, parameters_config) {
  # priors and global data markers
  priors <- parameters_config$prior_distribution_parameters
  prior_midpoint_values <- parameters_config$prior_midpoint_values

  data_full <- isolation::create_stan_variables(
    data,
    parameters_config
  )
  symptom_data <- data_full$symptom_data
  data <- data_full$by_day_data # `data` shorthand since used so often

  vl_list_for_stan <- list(
    # base items about data
    N = nrow(data),
    n_id = max(data$contiguous_id),
    id = data$contiguous_id, # code expects contiguous IDs
    # data
    t = data$days_since_symp_onset,
    y = data$logVL,
    pcr_avail = data$pcr_avail,
    # this is the by data set LOD/LOQ value used in Stan likelihood function
    lod = data$lod,
    below_lod = data$below_lod,
    loq_lower = data$loq_lower,
    below_loq_lower = data$below_loq_lower,
    loq_upper = data$loq_upper,
    VL_within_detectable_range = data$VL_within_detectable_range,
    above_loq_upper = data$above_loq_upper,
    antigen_result = data$antigen,
    antigen_avail = data$antigen_avail,
    culture_result = data$culture,
    culture_avail = data$culture_avail,
    si_time = symptom_data$symp_duration,
    cens_duration = symptom_data$symp_duration_censored,
    max_si = parameters_config$scenario_parameters$max_si,
    # everything is relative to midpoint values
    dp_midpoint = prior_midpoint_values$dp,
    wp_midpoint = prior_midpoint_values$wp,
    wr_midpoint = prior_midpoint_values$wr,
    # priors
    tp_prior = priors$tp_prior, # tp_prior is a vector from the yaml
    tp_std_prior = priors$tp_std_prior,
    sigma_prior = priors$sigma_prior, # sigma_prior is a vector from the yaml
    prior_sd = priors$prior_sd,
    antigen_50_prior = priors$antigen_50_prior, # vector
    culture_50_prior = priors$culture_50_prior, # vector
    culture_beta_prior = priors$culture_beta_prior, # vector
    si_beta_0_prior = priors$si_beta_0_prior, # vector
    si_beta_wr_prior = priors$si_beta_wr_prior, # vector
    si_shape_prior = priors$si_shape_prior # vector
  )
  return(vl_list_for_stan)
}
