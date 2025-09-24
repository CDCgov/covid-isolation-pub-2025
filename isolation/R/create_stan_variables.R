#' Create Stan model indicator variables
#'
#' for easy running of the stan model, indicator variables
#' and such data manipulations are done outside the model
#'
#' @param data output of stan_preprocess
#' @param parameters_config parameters_config.yml
#'
#' @return list of dataframes with additional variables needed for stan
#' @export
create_stan_variables <- function(data, parameters_config) {
  # global data markers
  global_data_markers <- parameters_config$global_data_markers

  # the data dataframe has multiple entries for each person
  # for each VL measurement
  # but we have added variables that are individual-specific
  # so they are repeated, such as the symptom duration
  # we want to extract the symptom duration and whether it is censored
  # by taking the first value per person
  # but, we also want to check that there aren't multiple values
  # of the symptom duration for an individual
  symptom_data <- data |>
    dplyr::select(symp_duration, symp_duration_censored, contiguous_id) |>
    dplyr::group_by(contiguous_id) |>
    dplyr::summarise(
      max_symp_duration = max(symp_duration),
      min_symp_duration = min(symp_duration),
      max_symp_duration_censored = max(symp_duration_censored),
      min_symp_duration_censored = min(symp_duration_censored),
      symp_duration = dplyr::first(symp_duration),
      symp_duration_censored = dplyr::first(symp_duration_censored)
    ) |>
    dplyr::ungroup()
  # we expect each individual to have the same value of symp duration
  # and symp duration censored for all their occurences in the dataframe
  stopifnot(symptom_data$min_symp_duration == symptom_data$max_symp_duration)
  stopifnot(
    symptom_data$min_symp_duration_censored ==
      symptom_data$max_symp_duration_censored
  )

  # add the indicator variables for whether
  # test is available and above or below LOQ/LOD
  data <- data |>
    dplyr::mutate(
      pcr_avail = as.integer(logVL != global_data_markers$skipped_test),
      below_lod = as.integer(logVL == global_data_markers$lod),
      below_loq_lower = as.integer(logVL == global_data_markers$loq_lower),
      above_loq_upper = as.integer(logVL == global_data_markers$loq_upper),
      antigen_avail = as.integer(antigen != global_data_markers$skipped_test),
      culture_avail = as.integer(culture != global_data_markers$skipped_test)
    ) |>
    dplyr::mutate(VL_within_detectable_range = as.integer(
      !below_lod & !below_loq_lower & !above_loq_upper & pcr_avail
    ))

  # above_loq_upper is annoying in that RVTN does not have
  # an loq_upper, so we have labeled their loq_upper as skipped_test
  # in `stan_preprocess`.
  # now, we want to check that all values for which loq_upper is skipped test
  # have above_loq_upper as 0 -- meaning none of them are above
  # this loq_upper we have created just to make the data analysis easier
  stopifnot(data$above_loq_upper[data$loq_upper == global_data_markers$skipped_test] == 0) # nolint: line_length_linter.

  if (
    global_data_markers$skipped_test > global_data_markers$loq_upper &&
      global_data_markers$lod < global_data_markers$loq_lower
  ) {
    # whenever above LOQ, logVL \in [loq, skipped_test]
    stopifnot(
      data$logVL[data$VL_within_detectable_range == 1] > global_data_markers$loq_lower & # nolint: line_length_linter.
        # technically this second condition is not guaranteed
        # because you could have chosen loq_upper < logVL
        # but that would be a not great choice around data markers
        # so if this condition is wrong, something is FUNKY with the data!
        data$logVL[data$VL_within_detectable_range == 1] < global_data_markers$loq_upper # nolint: line_length_linter.
    )
    # because VL_within_detectable_range does not include LOD/LOQ values
    stopifnot(sum(data$VL_within_detectable_range) <= sum(data$pcr_avail))
  }
  # below_lod and below_loq_lower are subsets of not missing (i.e., available)
  stopifnot(
    sum(data$below_lod) +
      sum(data$below_loq_lower) +
      sum(data$above_loq_upper) <= sum(data$pcr_avail)
  )
  # all conditions should together account for all the data
  stopifnot(
    sum(data$below_lod) +
      sum(data$below_loq_lower) +
      sum(data$above_loq_upper) +
      sum(!data$pcr_avail) + # missing values
      # actually VL measurements that can be interpreted as such
      sum(data$VL_within_detectable_range) == nrow(data)
  )
  # neither should be one at the same time
  stopifnot(
    sum(
      as.logical(data$below_lod) &
        as.logical(data$below_loq_lower) &
        as.logical(data$above_loq_upper) &
        as.logical(!data$pcr_avail) &
        as.logical(data$VL_within_detectable_range)
    ) == 0
  )

  return(list(
    symptom_data = symptom_data,
    by_day_data = data
  ))
}
