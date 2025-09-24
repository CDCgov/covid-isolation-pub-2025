#' Calculate dt accounting for mid-point approximation
#'
#' Make a vector of the length of time intervals (i.e., dt), but with the
#' first and last values halved.
#'
#' @param t A vector of time points of viral load observations.
#'
#' @return A vector of dt values the same length as t,
#' except with start and end values set equal to 1/2 * dt.
#' @export
#'
#' @examples
#' # Here dt is 0.1; this should return vector of length 31
#' # where values are 0.1, except first and last values are 0.05
#' dt_mid_approx_fx(seq(from = -1, to = 2, by = 0.1))
get_dt_mid_approx <- function(t) {
  dt <- diff(t)
  if (length(dt) > 0) {
    if (diff(range(dt)) > .Machine$double.eps^0.5) {
      warning("dts are inconsistent")
    }
  }
  dt <- c(dt[1] / 2, dt)
  dt[length(dt)] <- dt[length(dt)] / 2
  return(dt)
}
