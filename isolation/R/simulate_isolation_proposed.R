#' Simulate updated guidance for isolation and post-isolation precautions
#'
#' Simulate how much transmission potential may be averted by following
#' CDC's guidance for duration of isolation and post-isolation precautions for
#' persons with COVID-19, based on the guidance the updated guidance from
#' 2024-03-01. (The reference to "proposed" in the function name is because the
#' initial version of this code was prepared prior to 2024-03-01.)
#'
#' @param It Vector of time-varying infectiousness obtained from
#'   `simulate_intrinsic_infectiousness()`
#' @param t Vector of time values
#' @param symp_improv_in Observed time of symptom improvement
#' @param symp_improv_adjust Presumed time between actual and observed symptom
#'   improvement time
#' @param iso_lag Lag from symptom onset to isolation
#' @param iso_lag_adjust Presumed time between actual symptom onset and end of
#'   "Day 0", which is t=0
#' @param min_post_iso_days Minimum days in post-isolation precautions
#' @param iso_eff Efficacy of isolation in reducing transmission potential
#' @param post_iso_eff Efficacy of post-isolation precautions in reducing
#'   transmission potential
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
simulate_isolation_proposed <- function(It, t, symp_improv_in,
                                        symp_improv_adjust,
                                        iso_lag, iso_lag_adjust,
                                        min_post_iso_days,
                                        iso_eff, post_iso_eff) {
  iso_lag <- iso_lag - iso_lag_adjust

  # adjust observed symptom improvement time by 1.5 to account for partial days
  # and the data generating mechanism.
  # Suppose someone has an observed symptom improvement time of 4. This means
  # that improvement actually happened on Day 3 (e.g., if a person reported 2
  # symptoms in the 24 hour period ending at t = 4, and reported 3 symptoms
  # during the 24 hour period ending at t = 3, then the decrease in symptoms
  # must have occurred between t = 2 and t = 3).
  symp_improv <- symp_improv_in - symp_improv_adjust

  # make indicator for isolation period
  iso_end <- symp_improv
  iso_times <- as.numeric(t > iso_lag & t < iso_end)
  if (iso_end > iso_lag) {
    iso_times[t %in% c(iso_lag, iso_end)] <- 0.5
  }

  # make indicator for post-isolation period
  # include condition that if symptoms improve before iso_lag,
  # the person does not engage in post isolation precautions
  post_iso_end <- symp_improv + min_post_iso_days
  post_iso_times <- as.integer(
    t > iso_end & t < post_iso_end & iso_end >= iso_lag
  )

  if (iso_end >= iso_lag) {
    post_iso_times[t %in% c(iso_end, post_iso_end)] <- 0.5
  }

  I_t_intervention <- It * (
    1 - (iso_times * iso_eff + post_iso_times * post_iso_eff)
  )

  inf_averted <- isolation::get_infectiousness_averted(
    I_guid = I_t_intervention,
    I_null = It
  )


  return(list(
    p_iso = iso_times,
    p_post_iso = post_iso_times,
    I_t = I_t_intervention,
    averted = inf_averted
  ))
}
