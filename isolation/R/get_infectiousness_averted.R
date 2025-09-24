#' Percentage of infectiousness averted
#'
#' Calculate the percentage reduction in onward transmission potential when
#' following guidance versus in the absence of any intervention
#'
#' @param I_guid infectiousness timeseries when following guidance
#' @param I_null infectiousness timeseries when no intervention
#'
#' @return value of the percentage of infectiousness averted
#' @export
get_infectiousness_averted <- function(I_guid, I_null) {
  I_null_sum <- sum(I_null, na.rm = TRUE)
  I_guid_sum <- sum(I_guid, na.rm = TRUE)
  return(((I_null_sum - I_guid_sum) / I_null_sum) * 100)
}
