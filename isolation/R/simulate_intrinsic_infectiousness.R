#' Simulate intrinsic time-varying infectiousness
#'
#' Create a time series of an individual's intrinsic infectiousness (i.e.,
#' their potential for infecting others in the absence of interventions).
#' Multiple functional forms are supported to relate log viral load, antigen or
#' culture positivity, and combinations thereof to infectiousness.
#'
#' @param t time vector
#' @param vl_curve time-varying log viral load values for one person
#' @param prob_antigen_positive probabilities of antigen assay testing positive
#' @param prob_culture_positive probabilities of having culturable virus
#' @param infec_threshold logVL to use as base of infectiousness triangle
#' @param transform_fun viral load to infectiousness transformation (logVL,
#' antigen, culture, logVL_x_antigen, or logVL_x_culture)
#'
#' @return a time series of infectiousness, I(t), multiplied by time step, dt;
#' summing these values is numerical integration of area under infectiousness
#' curve
#' @export
simulate_intrinsic_infectiousness <- function(t, vl_curve,
                                              prob_antigen_positive,
                                              prob_culture_positive,
                                              infec_threshold,
                                              transform_fun) {
  dt <- get_dt_mid_approx(t = t)
  if (!(transform_fun %in% c(
    "logVL", "antigen", "culture",
    "logVL_x_antigen", "logVL_x_culture"
  ))) {
    stop("Simulate nat hist: transform_fun not valid")
  }
  if (transform_fun == "logVL") {
    I_t_dt <- pmax(0, vl_curve - infec_threshold) * dt
    return(I_t_dt)
  } else if (transform_fun == "antigen") {
    I_t_dt <- prob_antigen_positive * dt
    return(I_t_dt)
  } else if (transform_fun == "culture") {
    I_t_dt <- prob_culture_positive * dt
    return(I_t_dt)
  } else if (transform_fun == "logVL_x_antigen") {
    I_t_dt <- pmax(0, vl_curve - infec_threshold) * prob_antigen_positive * dt
    return(I_t_dt)
  } else if (transform_fun == "logVL_x_culture") {
    I_t_dt <- pmax(0, vl_curve - infec_threshold) * prob_culture_positive * dt
    return(I_t_dt)
  }
}
