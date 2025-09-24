#' Simulate previous guidance for isolation and post-isolation precautions
#'
#' Simulate how much transmission potential may be averted by following
#' CDC's guidance for duration of isolation and post-isolation precautions for
#' persons with COVID-19, based on the guidance that existed prior to
#' 2024-03-01. (The reference to "current" in the function name is because the
#' initial version of this code was prepared prior to 2024-03-01.)
#'
#' @param It Vector of time-varying infectiousness obtained from
#'   `simulate_intrinsic_infectiousness()`
#' @param t Vector of time values
#' @param prob_antigen_positive Vector of probabilities of being antigen
#'   positive
#' @param symp_type_cat An integer indicating person's symptom category
#' @param symp_improv Observed time of symptom improvement
#' @param iso_lag Lag from symptom onset to isolation
#' @param test_count Number of antigen tests taken (if remaining symptomatic)
#' @param test_interval Time between antigen tests
#' @param min_iso_days_moderate Minimum days in isolation from symptom onset
#'   for persons with moderate symptoms
#' @param min_iso_days_mild Minimum days in isolation from symptom onset
#'   for persons with mild symptoms
#' @param min_total_precaution_days Minimum total days spent in isolation and
#'   post-isolation precautions combined
#' @param iso_eff Efficacy of isolation in reducing transmission potential
#' @param post_iso_eff Efficacy of post-isolation precautions in reducing
#'   transmission potential
#' @param symp_improv_adjust Presumed time between actual and observed symptom
#'   improvement time
#' @param iso_lag_adjust Presumed time between actual symptom onset and end of
#'   "Day 0", which is t=0
#'
#' @return A list consisting of the following:
#'  * antigen_detect_prob: the probability the individual's infection is
#' detected by antigen testing
#'  * p_iso: vector of time-varying probability that a person is in isolation
#'  * p_post_iso: vector of time-varying probability that a person is taking
#' post-isolation precautions
#'  * I_t: time-varying transmission potential (area under the curve) when
#' following the guidance
#'  * averted : percentage of transmission potential averted by following the
#' guidance
#' @export
simulate_isolation_current <- function(It, t,
                                       prob_antigen_positive,
                                       symp_type_cat,
                                       symp_improv,
                                       iso_lag,
                                       test_count,
                                       test_interval,
                                       min_iso_days_moderate,
                                       min_iso_days_mild,
                                       min_total_precaution_days,
                                       iso_eff,
                                       post_iso_eff,
                                       symp_improv_adjust,
                                       iso_lag_adjust) {
  ## times when a person decides to test
  test_t <- seq(
    from = iso_lag,
    to = iso_lag + (test_count - 1) * test_interval,
    by = test_interval
  )

  ## probability of antigen detection based on times of testing
  ant_det <- get_antigen_prob_detection(
    test_t = test_t,
    t = t,
    antigen_sens = prob_antigen_positive,
    symp_improv = symp_improv,
    symp_improv_adjust = symp_improv_adjust,
    iso_lag_adjust = iso_lag_adjust
  )

  ## Calculate times during isolation and post-isolation
  if (symp_type_cat == 1) {
    min_iso_days <- min_iso_days_moderate
  } else if (symp_type_cat > 0 && symp_type_cat < 5) {
    min_iso_days <- min_iso_days_mild
  } else {
    stop("Only symptom categories allowed: 1-4")
  }

  # adjust "iso_lag" by 0.5 to account for partial days
  # an "iso_lag" of 0, logically, corresponds to immediate isolation
  # but t = 0 actually corresponds to the end of "Day 0"; we shift the start
  # of isolation back by 0.5 days to account for the fact that symptoms began
  # at some point during "Day 0".
  iso_lag <- iso_lag - iso_lag_adjust

  # adjust "symp_improv" by 1.5 to account for partial days and the data
  # generating mechanism
  # Suppose someone has an observed symptom improvement time of 4. This means
  # that improvement actually happened on Day 3 (e.g., if a person reported 2
  # symptoms in the 24 hour period ending at t = 4, and reported 3 symptoms
  # during the 24 hour period ending at t = 3, then the decrease in symptoms
  # must have occurred between t = 2 and t = 3).
  symp_improv <- symp_improv - symp_improv_adjust

  # make indicator for isolation period, assuming antigen detected
  iso_end_ant_pos <- max(min_iso_days, symp_improv)
  iso_times_ant_pos <- as.numeric(t > iso_lag & t < iso_end_ant_pos)
  if (iso_end_ant_pos > iso_lag) {
    iso_times_ant_pos[t %in% c(iso_lag, iso_end_ant_pos)] <- 0.5
  }

  # make indicator for post-isolation period, assuming antigen detected
  post_iso_end_ant_pos <- max(min_total_precaution_days, symp_improv)
  post_iso_times_ant_pos <- as.numeric(
    t > iso_end_ant_pos & t < post_iso_end_ant_pos
  )

  if (post_iso_end_ant_pos > iso_end_ant_pos) {
    post_iso_times_ant_pos[t %in% c(
      iso_end_ant_pos,
      post_iso_end_ant_pos
    )] <- 0.5
  }

  # make indicator for isolation period, assuming antigen not detected
  # make isolation end 0.5 days earlier to account for partial days
  iso_end_ant_neg <- symp_improv
  iso_times_ant_neg <- as.numeric(t > iso_lag & t < iso_end_ant_neg)
  if (iso_end_ant_neg > iso_lag) {
    iso_times_ant_neg[t %in% c(iso_lag, iso_end_ant_neg)] <- 0.5
  }

  # create indicator weighted by whether antigen detected
  iso_times <- iso_times_ant_pos * ant_det + iso_times_ant_neg * (1 - ant_det)
  post_iso_times <- post_iso_times_ant_pos * ant_det

  I_t_intervention <- It * (
    1 - (iso_times * iso_eff + post_iso_times * post_iso_eff)
  )
  inf_averted <- get_infectiousness_averted(
    I_guid = I_t_intervention,
    I_null = It
  )

  return(list(
    antigen_detect_prob = ant_det,
    p_iso = iso_times,
    p_post_iso = post_iso_times,
    I_t = I_t_intervention,
    averted = inf_averted
  ))
}
