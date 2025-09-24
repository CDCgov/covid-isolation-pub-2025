#' Calculate symptom category weights
#'
#' Get weights for each symptom category based on observed proportion
#'
#' @param joint_dataset a data frame that includes symp_type_cat
#'
#' @return A data frame with one row for each symp_type_cat and the
#' corresponding weight.
#' @export
calculate_symp_cat_weights <- function(joint_dataset) {
  if (!("symp_type_cat" %in% colnames(joint_dataset)) |
    !("contiguous_id" %in% colnames(joint_dataset))) {
    stop("required column names are missing")
  }

  symp_type_cat_df <- joint_dataset |>
    dplyr::distinct(symp_type_cat, contiguous_id)

  if (length(unique(symp_type_cat_df$contiguous_id)) !=
    nrow(symp_type_cat_df)) {
    stop("check for multiple categories in one contiguous id")
  }

  symp_cat_counts <- table(symp_type_cat_df$symp_type_cat)
  symp_cat_weights <- symp_cat_counts / sum(symp_cat_counts)
  symp_cat_weights <- data.frame(symp_cat_weights)

  if (ncol(symp_cat_weights) != 2) {
    stop("number of columns in symp_cat_weights incorrect")
  }

  colnames(symp_cat_weights) <- c("symp_type_cat", "wt")

  return(symp_cat_weights)
}
