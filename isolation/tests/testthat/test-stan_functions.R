# to reduce numeric instabilities, in stan define si_scale_inv
# as 1 / si_scale so that stan isn't continuously
# inverting something

# this means that we need a function to convert between q
# (what R DiscreteWeibull takes) and si_scale_inv since the custom
# function in stan is based off just si_scale_inv
q_to_si_scale_inv <- function(q, beta) {
  # recall q = exp(-si_scale_inv ^ beta)
  # with our new definition of si_scale_inv
  return((-log(q))^(1 / beta))
}

# max symp improv time
max_si <- 28
seed <- 9

set.seed(seed)

# get the functions form the stan model
symp_model_path <- system.file("stan", "triangle_vl_symp_hazard.stan",
  package = "isolation", mustWork = TRUE
)
rstan::expose_stan_functions(symp_model_path)

# first check for a range of q, beta, and si times -- do we get the
# same log cdf from Stan as we do from DiscreteWeibull?
# do cdf first because the other functions in the truncation process
# usd the CDF
testthat::test_that(
  "Discrete Weibull CDF in Stan implemented same
  as in DiscreteWeibull R package?",
  {
    # q is constrainted from 0 to to 1
    for (q in seq(0.1, 0.9, 0.1)) {
      # recognizing that beta = shape > 2 easily
      # causes numerical errors for big day
      # also -- only need to test three values -- < 1, = 1, > 1
      # these are the three distribution shapes
      for (beta in seq(0.5, 1.5, 0.5)) {
        for (day in seq(
          1,
          max_si,
          1
        )) {
          testthat::expect_equal(
            discrete_weibull_lcdf(
              day,
              # need to convert from q to si_scale_inv for stan fcn
              q_to_si_scale_inv(q, beta),
              beta
            ),
            # recall that discrete weibull package has no right truncation
            # right truncated pdf = unnormalized pdf / truncated CDF
            # or in log space, log(unnormalized pdf) - log(CDF)
            log(DiscreteWeibull::pdweibull(day,
              q,
              beta,
              zero = FALSE
            ))
          )
        }
      }
    }
  }
)

# first check for a range of q, beta, and si times -- do we get the
# same log pmf from Stan as we do from DiscreteWeibull?
testthat::test_that(
  "Discrete Weibull PDF in Stan implemented same
  as in DiscreteWeibull R package?",
  {
    # q is constrainted from 0 to to 1
    for (q in seq(0.1, 0.9, 0.1)) {
      # recognizing that beta = shape > 2 easily
      # causes numerical errors for big day
      # also -- only need to test three values -- < 1, = 1, > 1
      # these are the three distribution shapes
      for (beta in seq(0.5, 1.5, 0.5)) {
        for (day in seq(
          1,
          max_si,
          1
        )) {
          testthat::expect_equal(
            discrete_weibull_truncated_lpmf(
              day,
              # need to convert from q to si_scale_inv for stan fcn
              q_to_si_scale_inv(q, beta),
              beta,
              max_si # stan probability is conditioned on max_si
            ),
            # recall that discrete weibull package has no right truncation
            # right truncated pdf = unnormalized pdf / truncated CDF
            # or in log space, log(unnormalized pdf) - log(CDF)
            log(DiscreteWeibull::ddweibull(day,
              q,
              beta,
              zero = FALSE
            )) - log(DiscreteWeibull::pdweibull(max_si,
              q,
              beta,
              zero = FALSE
            ))
          )
        }
      }
    }
  }
)

# same check for the survival function
testthat::test_that(
  "Discrete Weibull survival function in Stan implemented
  same as in DiscreteWeibull R package?",
  {
    for (q in c(0.3, 0.6)) {
      for (beta in seq(0.5, 1.5, 0.5)) {
        for (day in seq(
          1,
          # recognizing that there are A LOT of numerical
          # instability problems for large day values
          # because in stan we directly calculate the log cCDF
          # but in R we are calculating the log of 1 - CDF manually
          max_si - 19,
          1
        )) {
          # like above, recall that stan discrete weibull
          # is conditioned on max_si
          # but DiscreteWeibull package is not
          # cdf of truncated dist is CDF(day) / CDF(max_si)
          # so cCDF is 1 - CDF(day) / CDF(max_si)
          testthat::expect_equal(
            discrete_weibull_truncated_lccdf(
              day,
              q_to_si_scale_inv(q, beta),
              beta,
              max_si
            ),
            log1p(-DiscreteWeibull::pdweibull(
              day,
              q,
              beta,
              zero = FALSE
            ) / DiscreteWeibull::pdweibull(
              max_si,
              q,
              beta,
              zero = FALSE
            )),
            tolerance = 0.01
          )
        }
      }
    }
  }
)

# in reality, we simulate symp improv times by drawing from a
# random variable discrete weibull distribution
# check that the random variable discrete weibull follows
# the same distribution that is being modeled in stan
# by empirically calculating the probability of observing
# values and checking it against calculated PMF from Stan
# need to be careful to account for truncating the distribution at max_si
testthat::test_that(
  "Discrete Weibull in Stan produces probabilities
    consistent with symp improv simulation times?",
  {
    for (rep in 1:10) {
      # draw some random q, beta
      rand_q <- runif(n = 1, min = 0, max = 1)
      rand_beta <- runif(n = 1, min = 0.5, max = 1.5)
      # draw some random si times
      # (conditioned on it being less than max_si)
      r_rand_sis <- sapply(1:10000, function(x) {
        isolation::estimate_symptom_improvement(
          # need the 1 / si_scale_inv because this part
          # uses the conventional def of si_scale
          1 / q_to_si_scale_inv(rand_q, rand_beta),
          rand_beta,
          max_si
        )
      })
      stan_rand_sis <- sapply(1:10000, function(x) {
        discrete_weibull_truncated_rng(
          q_to_si_scale_inv(rand_q, rand_beta),
          rand_beta,
          max_si
          # set random number generator
          # this line -- despite the dynamic seed --
          # fixes the output of the function. this seems to be an rstan problem?
          # base_rng__ = rstan::get_rng(seed = seed + x) # nolint: commented_code_linter, line_length_linter.
        )
      })
      # empirical calculation of probabilities of occurence
      r_empirical_probabilities <- data.frame(prop.table(table(r_rand_sis)))
      stan_empirical_probabilities <- data.frame(
        prop.table(table(stan_rand_sis))
      )

      empirical_probabilities <- dplyr::full_join(r_empirical_probabilities,
        stan_empirical_probabilities,
        by = c("r_rand_sis" = "stan_rand_sis")
      ) |>
        # there may be some si values observed in one way of
        # generating but not the other
        dplyr::mutate_if(is.numeric, dplyr::coalesce, 0) |>
        dplyr::select(
          si_time = r_rand_sis,
          r_empirical_probability = Freq.x,
          stan_empirical_probability = Freq.y
        )

      # for each si time that shows up, check that the empirical probability
      # and stan probabilities match
      for (row_num in seq_len(nrow(empirical_probabilities))) {
        si_time <- empirical_probabilities$si_time[row_num]
        r_empirical_probability <- empirical_probabilities$r_empirical_probability[row_num] # nolint: line_length_linter.
        stan_empirical_probability <- empirical_probabilities$stan_empirical_probability[row_num] # nolint: line_length_linter.
        stan_analytical_probability <- exp(
          discrete_weibull_truncated_lpmf(
            si_time,
            q_to_si_scale_inv(rand_q, rand_beta),
            rand_beta,
            max_si
          )
        )
        testthat::expect_lte(
          abs(
            r_empirical_probability - stan_analytical_probability
          ), # nolint: indentation_linter.
          0.02
        )
        testthat::expect_lte(
          abs(
            stan_empirical_probability - stan_analytical_probability
          ), # nolint: indentation_linter.
          0.02
        )
      }
    }
  }
)

# since we have stan functions loaded,
# also test the triangle viral load function
# and the bernoulli_logit_lpmf
ts <- seq(-7, 15, 0.5)
tp <- 2
dp <- 8
wp <- 5
wr <- 9
shift <- 3
testthat::test_that("LogVL in Stan same as in R for scenarios", {
  # tri_vl is not vectorized and expects a singular real as t
  # so we need to use sapply to calculate the function over multiple ts
  testthat::expect_equal(
    triangle_vl(t = ts, dp = dp, tp = tp, wp = wp, wr = wr),
    sapply(ts, function(x) {
      tri_vl(t = x, tp = tp, wp = wp, wr = wr, dp = dp)
    })
  )
  # test that the triangle VL functions have the expected behavior
  # if we decrease the proliferation/clearance times,
  # we expect the curve to again always be lower
  testthat::expect_true(all(triangle_vl(
    t = ts,
    dp = dp, tp = tp, wp = wp,
    wr = wr
  ) - triangle_vl(
    t = ts,
    dp = dp, tp = tp, wp = wp - 1,
    wr = wr - 1
  ) >= 0))
  testthat::expect_true(all(sapply(
    ts,
    function(x) {
      tri_vl(t = x, tp = tp, wp = wp, wr = wr, dp = dp)
    }
  ) - sapply(
    ts,
    function(x) {
      tri_vl(t = x, tp = tp, wp = wp - 1, wr = wr - 1, dp = dp)
    }
  ) >= 0))
  # check that shifting the peak timing actually just shifts the curve
  testthat::expect_equal(
    triangle_vl(t = ts + shift, dp = dp, tp = tp, wp = wp, wr = wr),
    triangle_vl(t = ts, dp = dp, tp = tp - shift, wp = wp, wr = wr)
  )
  testthat::expect_equal(
    sapply(
      ts + shift,
      function(x) {
        tri_vl(t = x, tp = tp, wp = wp, wr = wr, dp = dp)
      }
    ),
    sapply(
      ts,
      function(x) {
        tri_vl(t = x, tp = tp - shift, wp = wp, wr = wr, dp = dp)
      }
    )
  )
})

testthat::test_that("Bernoulli logit lpmf implemented in
Stan same as scenarios", {
  k50 <- 4.5
  scale <- 0.4
  testthat::expect_equal(
    log(prob_antigen_positive(
      triangle_vl(t = ts, dp = dp, tp = tp, wp = wp, wr = wr),
      k50, scale
    )),
    sapply(
      triangle_vl(t = ts, dp = dp, tp = tp, wp = wp, wr = wr),
      function(x) {
        # wrapped bernoulli logit lpmf just wraps stan's bernoulli
        # logit lpmf function which calculates the probability of
        # observing a 0 or 1 (bernoulli) via a logistic function
        # wrapping the function makes it accessible to R, so we
        # can test that our implementation of the function in R is correct
        # we are interested in comparing our calculation of the probability
        # *positive*, so we always want to calculate the probability
        # of observing a 1 (first argument) as a function of some underlying
        # x value based on some k50 and scale
        wrapped_bernoulli_logit_lpmf(1, x, k50, 1 / scale)
      }
    )
  )
})
