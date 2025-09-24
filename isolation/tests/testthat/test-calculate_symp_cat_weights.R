# testing the function calculate_symp_cat_weights

# tiny toy dataset
example_data <- data.frame(
  "symp_type_cat" = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4),
  "contiguous_id" = c(1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 8)
)

# randomize the order of rows
example_data <- example_data[sample(1:nrow(example_data), replace = F), ]

# hand-calculated correct result for toy dataset
example_output <- data.frame(
  "symp_type_cat" = as.factor(c(1, 2, 3, 4)),
  "wt" = c(3 / 8, 2 / 8, 1 / 8, 2 / 8)
)

# create misnamed columns in dataset (missing first "p" in "symp_type_cat")
misnamed_data <- data.frame(
  "sym_type_cat" = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4),
  "contiguous_id" = c(1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 8)
)

# create dataset with multiple categories per contiguous id
duplicate_contig_id_data <- data.frame(
  "sym_type_cat" = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4),
  "contiguous_id" = c(1, 1, 2, 3, 3, 4, 5, 6, 6, 7, 8)
)

test_that("quantitative results from calculate_symp_cat_weights are as expected", {
  expect_equal(
    calculate_symp_cat_weights(joint_dataset = example_data),
    example_output
  )
})

test_that("returns error when incorrect column names", {
  expect_error(calculate_symp_cat_weights(joint_dataset = misnamed_data))
})

test_that("returns error when multiple categories per id", {
  expect_error(calculate_symp_cat_weights(joint_dataset = duplicate_contig_id_data))
})
