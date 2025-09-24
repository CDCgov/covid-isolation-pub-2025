#' Calculate probability an antigen test is positive
#'
#' Probability an antigen test is positive depends on log viral load and
#' parameters that relate viral load to test sensitivity
#'
#' @param logVL underlying log viral load
#' @param antigen_50 logVL at which antigen positive probability is 50%
#' @param sigma_antigen "steepness"/variability in sigmoid curve
#' as antigen probability increases with increasing logVL
#'
#' @return probability antigen positive for each given logVL
#' @export
prob_antigen_positive <- function(
    logVL, # nolint: object_name_linter.
    antigen_50, sigma_antigen) {
  prob_antigen_positive <- plogis(
    q = logVL,
    location = antigen_50,
    scale = sigma_antigen
  )

  return(prob_antigen_positive)
}
