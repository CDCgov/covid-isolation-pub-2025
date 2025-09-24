# test get_dt_mid_approx function

example_input <- seq(from = -1, to = 2, by = 0.1)

# this is based on example input has overall length 31 and the "ends"
# should be 1/2 of dt
expected_output <- c(0.05, rep(x = 0.1, times = 29), 0.05)

test_that("get_dt_mid_approx behaves as expected", {
  expect_equal(
    get_dt_mid_approx(example_input),
    expected_output
  )
  expect_true(is.na(get_dt_mid_approx(t = 1)))
  expect_warning(get_dt_mid_approx(c(1, 2, 4)))
})
