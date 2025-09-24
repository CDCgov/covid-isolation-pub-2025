#' Weighted sample to reflect frequency of symptom categories
#'
#' @param symp_id_df a data frame relating posterior IDs and symptom categories
#' @param symp_cat_wts a data frame with weights for each symptom category
#' @param n_overall_sample number of samples in the overall "category"
#'
#' @return a dataframe of posterior IDs that have been sampled
#' @export
#'
sample_ids_overall_symp_cat <- function(symp_id_df, symp_cat_wts,
                                        n_overall_sample) {
  stopifnot(sum(symp_cat_wts$wt) == 1)

  # calculate the number of samples that should be drawn from each category

  symp_cat_wts <- symp_cat_wts |>
    dplyr::mutate(expected_count = wt * n_overall_sample) |>
    dplyr::mutate(integer_count = round(x = expected_count, digits = 0))

  # ensure that number to draw per category sums to the total target size

  discrepancy <- sum(symp_cat_wts$integer_count) - n_overall_sample

  if (discrepancy > 0) {
    # decrease by 1 those categories that were rounded up the most
    symp_cat_wts <- symp_cat_wts |>
      dplyr::mutate(round_diff = integer_count - expected_count) |>
      dplyr::arrange(desc(round_diff))

    symp_cat_wts$integer_count[1:discrepancy] <-
      symp_cat_wts$integer_count[1:discrepancy] - 1
  }

  if (discrepancy < 0) {
    # increase by 1 those categories that were rounded down the most
    symp_cat_wts <- symp_cat_wts |>
      dplyr::mutate(round_diff = integer_count - expected_count) |>
      dplyr::arrange(round_diff)

    symp_cat_wts$integer_count[1:abs(discrepancy)] <-
      symp_cat_wts$integer_count[1:abs(discrepancy)] + 1
  }

  stopifnot(sum(symp_cat_wts$integer_count) == n_overall_sample)

  symp_cat_wts <- symp_cat_wts |> dplyr::arrange(symp_type_cat)

  stratified_samples_df_list <- purrr::map2(
    .x = symp_cat_wts$symp_type_cat,
    .y = symp_cat_wts$integer_count,
    .f = ~ {
      symp_id_df_subset <- dplyr::filter(symp_id_df, symp_type_cat == .x)
      sampled_rows <- sample(x = 1:nrow(symp_id_df_subset), size = .y, replace = FALSE) # nolint
      return(symp_id_df_subset[sampled_rows, ])
    }
  )

  for (i in seq_along(stratified_samples_df_list)) {
    symp_id_df_subset <- symp_id_df |>
      dplyr::filter(symp_type_cat == symp_cat_wts$symp_type_cat[i])
    stratified_samples_df_list[[i]] <- symp_id_df_subset[
      sample(
        x = seq_len(nrow(symp_id_df_subset)),
        size = symp_cat_wts$integer_count[i],
        replace = FALSE
      ),
    ]
  }

  sampled_ids_df <- dplyr::bind_rows(stratified_samples_df_list)

  stopifnot(nrow(sampled_ids_df) == n_overall_sample)

  counts_by_cat_in_sample <- data.frame(table(sampled_ids_df$symp_type_cat))

  stopifnot(counts_by_cat_in_sample$Var1 == symp_cat_wts$symp_type_cat &
    counts_by_cat_in_sample$Freq == symp_cat_wts$integer_count) # nolint

  return(sampled_ids_df)
}
