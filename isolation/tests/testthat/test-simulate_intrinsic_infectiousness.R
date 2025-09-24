# testing the function simulate_intrinsic_infectiousness

example_t <- seq(from = -2, to = 10, by = 0.1)

example_tp <- 1.5

# the specific numbers used in this test are arbitrary,
# but broadly similar to the observed values

example_vl_curve <- isolation::triangle_vl(
  t = example_t,
  dp = 7,
  tp = example_tp,
  wp = 3,
  wr = 7
)

example_prob_antigen_positive <- plogis(
  q = example_vl_curve,
  location = 5.5,
  scale = 0.8
)

example_prob_culture_positive <- plogis(
  q = example_vl_curve,
  location = 5.5 + 0.5 * (example_t - example_tp),
  scale = 0.8
)

# the following test checks the area under the infectiousness curve
# for three different values for the infec_threshold (0, -0.5, and 0.5)
# this test calculates the area under the curve using the function,
# and then compares to the same value calculated manually.

test_that("Area under infectiousness curve is as expected for logVL", {
  testthat::expect_equal(
    object = sum(isolation::simulate_intrinsic_infectiousness(
      t = example_t,
      vl_curve = example_vl_curve,
      prob_antigen_positive = example_prob_antigen_positive,
      prob_culture_positive = example_prob_culture_positive,
      infec_threshold = 0,
      transform_fun = "logVL"
    )),
    expected = sum(
      example_vl_curve[example_vl_curve > 0] *
        isolation::get_dt_mid_approx(example_t)[example_vl_curve > 0]
    )
  )
  testthat::expect_equal(
    object = sum(isolation::simulate_intrinsic_infectiousness(
      t = example_t,
      vl_curve = example_vl_curve,
      prob_antigen_positive = example_prob_antigen_positive,
      prob_culture_positive = example_prob_culture_positive,
      infec_threshold = 0.5,
      transform_fun = "logVL"
    )),
    expected = sum(
      (example_vl_curve - 0.5)[(example_vl_curve - 0.5) > 0] *
        isolation::get_dt_mid_approx(example_t)[(example_vl_curve - 0.5) > 0]
    )
  )
  testthat::expect_equal(
    object = sum(isolation::simulate_intrinsic_infectiousness(
      t = example_t,
      vl_curve = example_vl_curve,
      prob_antigen_positive = example_prob_antigen_positive,
      prob_culture_positive = example_prob_culture_positive,
      infec_threshold = -0.5,
      transform_fun = "logVL"
    )),
    expected = sum(
      (example_vl_curve + 0.5)[(example_vl_curve + 0.5) > 0] *
        isolation::get_dt_mid_approx(example_t)[(example_vl_curve + 0.5) > 0]
    )
  )
})

# the following four tests check additional infectiousness assumptions
# these tests are similar to the above test, but the difference is that
# we compare the actual time series of infectiousness, rather than the sum

test_that("Infectiousness curve is as expected for antigen", {
  testthat::expect_equal(
    object = isolation::simulate_intrinsic_infectiousness(
      t = example_t,
      vl_curve = example_vl_curve,
      prob_antigen_positive = example_prob_antigen_positive,
      prob_culture_positive = example_prob_culture_positive,
      infec_threshold = 0.5,
      transform_fun = "antigen"
    ),
    expected = example_prob_antigen_positive * get_dt_mid_approx(t = example_t)
  )
})

test_that("Infectiousness curve is as expected for culture", {
  testthat::expect_equal(
    object = isolation::simulate_intrinsic_infectiousness(
      t = example_t,
      vl_curve = example_vl_curve,
      prob_antigen_positive = example_prob_antigen_positive,
      prob_culture_positive = example_prob_culture_positive,
      infec_threshold = 0.5,
      transform_fun = "culture"
    ),
    expected = example_prob_culture_positive * get_dt_mid_approx(t = example_t)
  )
})

test_that("Infectiousness curve is as expected for logVL_x_antigen", {
  testthat::expect_equal(
    object = sum(isolation::simulate_intrinsic_infectiousness(
      t = example_t,
      vl_curve = example_vl_curve,
      prob_antigen_positive = example_prob_antigen_positive,
      prob_culture_positive = example_prob_culture_positive,
      infec_threshold = 0.5,
      transform_fun = "logVL_x_antigen"
    )),
    expected = sum(((example_vl_curve - 0.5) *
      example_prob_antigen_positive)[(example_vl_curve - 0.5) > 0] *
      isolation::get_dt_mid_approx(example_t)[(example_vl_curve - 0.5) > 0])
  )
})

test_that("Infectiousness curve is as expected for logVL_x_culture", {
  testthat::expect_equal(
    object = sum(isolation::simulate_intrinsic_infectiousness(
      t = example_t,
      vl_curve = example_vl_curve,
      prob_antigen_positive = example_prob_antigen_positive,
      prob_culture_positive = example_prob_culture_positive,
      infec_threshold = 0.5,
      transform_fun = "logVL_x_culture"
    )),
    expected = sum(((example_vl_curve - 0.5) *
      example_prob_culture_positive)[(example_vl_curve - 0.5) > 0] *
      isolation::get_dt_mid_approx(example_t)[(example_vl_curve - 0.5) > 0])
  )
})

test_that("Returns error if wrong transformation", {
  testthat::expect_error(isolation::simulate_intrinsic_infectiousness(
    t = example_t,
    vl_curve = example_vl_curve,
    prob_antigen_positive = example_prob_antigen_positive,
    prob_culture_positive = example_prob_culture_positive,
    infec_threshold = 0.1,
    transform_fun = "test"
  ))
})
