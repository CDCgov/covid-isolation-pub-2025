# test the isolation functions `prob_antigen_positive`
# and `prob_culture_positive`

# first let's test prob antigen positive
# let's test a base case sigmoid curve
# let's test the 50% value
# let's test that the derivative varies with logVL

synthetic_logVL_like_x_vals <- seq(-5, 5, 0.5) # nolint: object_name_linter.

sigmoid_y_base <- 1 / (1 + exp(-synthetic_logVL_like_x_vals))

shift <- 3
sigmoid_y_shifted <- 1 / (1 + exp(-(synthetic_logVL_like_x_vals - shift)))

scale <- 2
sigmoid_y_scaled <- 1 / (1 + exp(-(1 / scale) * synthetic_logVL_like_x_vals))

sigmoid_y_shifted_and_scaled <- 1 / (1 + exp(-((1 / scale) *
  (synthetic_logVL_like_x_vals - shift))))

testthat::test_that("prob antigen positive behaves
like a sigmoid as expected", {
  # test the base sigmoid
  testthat::expect_equal(
    sigmoid_y_base,
    prob_antigen_positive(synthetic_logVL_like_x_vals, 0, 1)
  )
  # test the shifted sigmoid
  testthat::expect_equal(
    sigmoid_y_shifted,
    prob_antigen_positive(synthetic_logVL_like_x_vals, shift, 1)
  )
  # test the scaled sigmoid
  testthat::expect_equal(
    sigmoid_y_scaled,
    prob_antigen_positive(synthetic_logVL_like_x_vals, 0, scale)
  )
  # test the shifted and scaled sigmoid
  testthat::expect_equal(
    sigmoid_y_shifted_and_scaled,
    prob_antigen_positive(synthetic_logVL_like_x_vals, shift, scale)
  )
})

# now test the 50% value
testthat::test_that("prob antigen 50% happens at expected values", {
  # test that the 50% value is 0 for a non shifted sigmoid
  testthat::expect_equal(
    0.5,
    prob_antigen_positive(0, 0, 1)
  )
  # 50% value shifts as we expect
  testthat::expect_equal(
    0.5,
    prob_antigen_positive(shift, shift, 1)
  )
  # 50% value is not impacted by scale
  testthat::expect_equal(
    0.5,
    prob_antigen_positive(shift, shift, scale)
  )
  # less than the shift should be smaller than 50%
  testthat::expect_gt(
    0.5,
    prob_antigen_positive(shift - 0.1, shift, scale)
  )
  # and vice versa
  testthat::expect_lt(
    0.5,
    prob_antigen_positive(shift + 0.1, shift, scale)
  )
})

# now test with some real logVL data
logVL <- triangle_vl(
  t = synthetic_logVL_like_x_vals,
  tp = 1.5,
  dp = 6,
  wp = 5,
  wr = 9
)
testthat::test_that("prob antigen positive varies with logVL", {
  # test that logVL and prob antigen positive covary
  testthat::expect_equal(
    sign(diff(logVL)),
    sign(diff(prob_antigen_positive(logVL, 0, 1)))
  )
  testthat::expect_equal(
    sign(diff(logVL)),
    sign(diff(prob_antigen_positive(logVL, shift, scale)))
  )
  testthat::expect_equal(
    which.max(logVL),
    which.max(prob_antigen_positive(logVL, shift, scale))
  )
})

# now perform similar tests for the culture function
# now we have to be careful about this additional term
# that describes how culture positivity also varies
# with time since the peak
sigmoid_y_time_varying_50 <- 1 / (1 + exp(-((1 / scale) *
  (synthetic_logVL_like_x_vals - (shift + scale * (synthetic_logVL_like_x_vals - shift))))))
testthat::test_that("prob culture positive behaves
like a sigmoid as expected", {
  # test the base sigmoid
  # base sigmoid under no time-varying change to the 50% value
  testthat::expect_equal(
    sigmoid_y_base,
    prob_culture_positive(synthetic_logVL_like_x_vals,
      synthetic_logVL_like_x_vals,
      culture_50 = 0, sigma_culture = 1,
      culture_beta = 0, tp = 0
    )
  )
  testthat::expect_equal(
    sigmoid_y_base,
    prob_culture_positive(synthetic_logVL_like_x_vals,
      synthetic_logVL_like_x_vals,
      culture_50 = 0, sigma_culture = 1,
      culture_beta = 0, tp = shift
    )
  )
  # are time-varying changes to the 50% value having the desired impact?
  testthat::expect_equal(
    sigmoid_y_time_varying_50,
    prob_culture_positive(synthetic_logVL_like_x_vals,
      synthetic_logVL_like_x_vals,
      culture_50 = shift, sigma_culture = scale,
      culture_beta = scale, tp = shift
    )
  )
})

testthat::test_that("prob culture positive 50% behaves as expected", {
  # test the 50% value under base sigmoid
  testthat::expect_equal(
    0.5,
    prob_culture_positive(shift,
      shift,
      culture_50 = shift, sigma_culture = scale,
      culture_beta = 0, tp = 0
    )
  )
  # tp should have no impact when culture_beta = 0
  testthat::expect_equal(
    0.5,
    prob_culture_positive(shift,
      shift,
      culture_50 = shift, sigma_culture = scale,
      culture_beta = 0, tp = shift
    )
  )
  # are time-varying changes to the 50% value having the desired impact?
  testthat::expect_equal(
    0.5,
    prob_culture_positive(shift - scale * shift,
      0,
      culture_50 = shift, sigma_culture = scale,
      culture_beta = scale, tp = shift
    )
  )
  testthat::expect_equal(
    0.5,
    prob_culture_positive(shift,
      shift,
      culture_50 = shift, sigma_culture = scale,
      culture_beta = scale, tp = shift
    )
  )
  # the 50% value changes with time -- we can test that this change
  # is working correctly by making sure that the 50% value drops
  testthat::expect_true(all(rep(0.5, length(synthetic_logVL_like_x_vals)) -
    prob_culture_positive(rep(shift, length(synthetic_logVL_like_x_vals)),
      synthetic_logVL_like_x_vals,
      culture_50 = shift, sigma_culture = scale,
      culture_beta = scale, tp = min(synthetic_logVL_like_x_vals)
    ) >= 0))
})

testthat::test_that("prob culture positive returns an error", {
  # test that the `stopifnot` error condition works correctly
  testthat::expect_error(prob_culture_positive(c(),
    c(0),
    culture_50 = shift, sigma_culture = scale,
    culture_beta = -scale,
    tp = shift
  ))
  testthat::expect_error(prob_culture_positive(c(0),
    c(),
    culture_50 = shift, sigma_culture = scale,
    culture_beta = -scale,
    tp = shift
  ))
  testthat::expect_error(prob_culture_positive(c(0),
    c(1, 2),
    culture_50 = shift, sigma_culture = scale,
    culture_beta = -scale,
    tp = shift
  ))
})
