#' Get viral load values from triangle parameters
#'
#' Calulate (log) viral load at given times based on parameters for
#' "triangular" VL function
#'
#' @param t times at which viral load should be estimated
#' @param dp peak value of viral load -- assumes that peak VL > LOD
#' @param tp time of peak vira load
#' @param wp proliferation time (time from logVL = 0 up to peak)
#' @param wr clearance time (time from peak back down to logVL = 0)
#'
#' @return triangular viral load at each time t
#' @export
triangle_vl <- function(t, dp, tp, wp, wr) {
  log_vl <- rep(NA, length(t))

  log_vl <- dplyr::if_else(t <= tp,
    (dp / wp) * (t - (tp - wp)),
    dp - (dp / wr) * (t - tp)
  )

  return(log_vl)
}
