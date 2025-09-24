library(testthat)
library(dplyr)

n_samples_per_cat <- 100

example_input_data <- dplyr::tibble(
  posterior_id = 1:(4 * n_samples_per_cat),
  symp_type_cat = sort(rep(x = c(1, 2, 3, 4), times = n_samples_per_cat))
)

wts_data <- dplyr::tibble(
  symp_type_cat = c(1, 2, 3, 4),
  wt = c(0.4, 0.3, 0.2, 0.1)
)

test_that("Has the right columns", {
  sampled_ids_df <- isolation::sample_ids_overall_symp_cat(
    symp_id_df = example_input_data,
    symp_cat_wts = wts_data,
    n_overall_sample = n_samples_per_cat
  )
  actual_names <- colnames(sampled_ids_df)
  expected_names <- c("symp_type_cat", "posterior_id")
  expect_setequal(actual_names, expected_names)
})

test_that("Category proportions match intended weights", {
  sampled_ids_df <- isolation::sample_ids_overall_symp_cat(
    symp_id_df = example_input_data,
    symp_cat_wts = wts_data,
    n_overall_sample = n_samples_per_cat
  )
  weights_by_cat <- sampled_ids_df %>%
    dplyr::group_by(symp_type_cat) %>%
    dplyr::summarize(N = n()) %>%
    dplyr::mutate(samp_proportion = N / sum(N)) %>%
    dplyr::left_join(wts_data, by = c("symp_type_cat"))
  expect_equal(weights_by_cat$samp_proportion, weights_by_cat$wt,
    tolerance = 1 / n_samples_per_cat
  )
})

###

wts_data_round_down <- dplyr::tibble(
  symp_type_cat = c(1, 2, 3, 4),
  wt = c(0.406, 0.307, 0.2, 0.087)
)

# naively, counts per category would be 41, 31, 20, and 9, which is >100
stopifnot(sum(wts_data_round_down$wt) == 1)
stopifnot(sum(round(wts_data_round_down$wt * 100)) > 100)

test_that("Extra rounding down works as expected", {
  sampled_ids_df <- isolation::sample_ids_overall_symp_cat(
    symp_id_df = example_input_data,
    symp_cat_wts = wts_data_round_down,
    n_overall_sample = 100
  )

  counts_by_cat_in_sample <- data.frame(table(sampled_ids_df$symp_type_cat))

  expect_equal(counts_by_cat_in_sample$Freq, c(40, 31, 20, 9))
})

###

wts_data_round_up <- dplyr::tibble(
  symp_type_cat = c(1, 2, 3, 4),
  wt = c(0.404, 0.303, 0.2, 0.093)
)

# naively, counts per category would be 40, 30, 20, and 9, which is <100
stopifnot(sum(wts_data_round_up$wt) == 1)
stopifnot(sum(round(wts_data_round_up$wt * 100)) < 100)

test_that("Extra rounding up works as expected", {
  sampled_ids_df <- isolation::sample_ids_overall_symp_cat(
    symp_id_df = example_input_data,
    symp_cat_wts = wts_data_round_up,
    n_overall_sample = 100
  )

  counts_by_cat_in_sample <- data.frame(table(sampled_ids_df$symp_type_cat))

  expect_equal(counts_by_cat_in_sample$Freq, c(41, 30, 20, 9))
})
