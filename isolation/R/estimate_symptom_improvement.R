#' Simulate symptom improvement time
#'
#' Simulates symptom improvment using a discreteWeibull
#'
#' @param si_scale escale parameter for weibull
#' @param si_shape shape parameter for weibull
#' @param max_si maximum time to symptom improvement
#'
#' @return estimated symptom improvement
#' @export
estimate_symptom_improvement <- function(si_scale, si_shape, max_si) {
  # unlike in the stan model,
  # here we use 1 / si_scale because in the stan extraction we transform
  # to calculate si_scale as exp(-something) which is the proper convention
  q_weibull <- exp(-(1 / si_scale)^si_shape)
  # initial sim_si value at first to ensure we sample a value
  sim_si <- max_si + 1
  while (sim_si > max_si) {
    # sample a value below or equal to max_si
    sim_si <- DiscreteWeibull::rdweibull(
      n = 1,
      q = q_weibull,
      beta = si_shape, zero = FALSE
    )
  }

  return(sim_si)
}
