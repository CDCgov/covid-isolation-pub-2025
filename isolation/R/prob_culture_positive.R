#' Calculate probability a culture assay is positive
#'
#' Probability a culture assay is positive depends on log viral load and time
#' relative to peak viral load, plus parameters that relate these inputs to
#' sensitivity.
#'
#' @param logVL underlying log viral load
#' @param t time at which logVL measurements were taken
#' @param culture_50 logVL at which culture positive
#' probability is 50% at t = tp
#' @param sigma_culture "steepness"/variability in sigmoid curve
#' as culture probability increases with increasing logVL
#' @param culture_beta change in log-odds of culture positive as time
#' since peak VL increases
#' @param tp peak viral load time
#'
#' @return probability of testing culture positive based on VL and time
#'   relative to peak viral load
#' @export
prob_culture_positive <- function(
    logVL, # nolint: object_name_linter.
    t, culture_50, sigma_culture, culture_beta, tp) {
  stopifnot(length(logVL) == length(t))

  prob_culture_positive <- plogis(
    q = logVL,
    location = culture_50 + culture_beta * (t - tp),
    scale = sigma_culture
  )

  return(prob_culture_positive)
}
