#' Calculate probability of detection by antigen testing
#'
#' The probability that a case will be detected through antigen
#' testing is calculated as 1 minus the product of the probabilities
#' that each test will be negative. This is one probability per case;
#' it is not time-varying. A test is not taken if the patient's symptoms
#' have already improved.
#'
#' @param test_t Times at which testing is indicated
#' @param t Vector of times at which antigen sensitivity is known
#' @param antigen_sens  Vector of antigen test sensitivities at time t
#' @param symp_improv  Observed time of symptom improvement
#' @param symp_improv_adjust Presumed time between actual and observed symptom
#'   improvement time
#' @param iso_lag_adjust Presumed time between actual symptom onset and end of
#'   "Day 0", which is t=0
#'
#' @return probability of antigen detection
#' @export
get_antigen_prob_detection <- function(test_t, t, antigen_sens,
                                       symp_improv, symp_improv_adjust,
                                       iso_lag_adjust) {
  # subtract 0.5 from test times to account for partial days
  # this parallels the approach taken for iso_lag
  # so a test time of 0 means at symptom onset rather than end of Day 0
  test_t <- test_t - iso_lag_adjust

  # adjust "symp_improv" by 1.5 to account for partial days and the data
  # generating mechanism
  # Suppose someone has an observed symptom improvement time of 4. This means
  # that improvement actually happened on Day 3 (e.g., if a person reported 2
  # symptoms in the 24 hour period ending at t = 4, and reported 3 symptoms
  # during the 24 hour period ending at t = 3, then the decrease in symptoms
  # must have occurred between t = 2 and t = 3).

  symp_improv <- symp_improv - symp_improv_adjust

  # assume that test are not taken if symptom improvement has already occurred
  # (but we do allow for tests at the same time as symptom improvement)
  return(1 - prod(1 - antigen_sens[t %in% test_t[test_t <= symp_improv]]))
}
