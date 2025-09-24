# test the function simulate_isolation_proposed.R
# this test also de facto tests get_infectiousness_averted.R

# set levels of natural history parameters

tp <- c(-1, 8) # times when peak occurs
wp <- c(1, 5) # proliferation time
wr <- c(1, 8) # clearance time
dp <- c(1, 9) # peak logVL
symp_improv_observed <- 1:16 # systematically checking across si times is key
antigen_50 <- 5.5 # logVL with 50% chance antigen positive
sigma_antigen <- 0.8 # scale parameter for logistic fx for antigen positive
transform_fun <- c("logVL", "antigen") # options for infectiousness model
infec_threshold <- 1.5 # threshold for base of logVL infectiousness "triangle"

# set levels of intervention parameters

iso_lag <- c(0, 2)
min_post_iso_days <- 5
iso_eff <- 0.9
post_iso_eff <- 0.5
symp_improv_adjust <- 1.5
iso_lag_adjust <- 0.5

# create all combinations of natural history parameters
param_combos <- tidyr::expand_grid(
  tp, wp, wr, dp,
  symp_improv_observed,
  antigen_50, sigma_antigen,
  transform_fun, infec_threshold,
  iso_lag, min_post_iso_days,
  iso_eff, post_iso_eff,
  symp_improv_adjust, iso_lag_adjust
)

# define the time series that will be used across all simulated persons
times <- seq(from = -10, to = 21, by = 0.1)

# set up vectors where summed I_t will be stored

summed_I_t_from_fx <- rep(NA, times = nrow(param_combos))
summed_I_t_from_test <- rep(NA, times = nrow(param_combos))

# set up vectors where averted infectiousness measures will be stored

infec_averted_from_fx <- rep(NA, times = nrow(param_combos))
infec_averted_from_test <- rep(NA, times = nrow(param_combos))

# set up list to store vectors replicating p_iso and p_post_iso

p_iso_vecs_from_fx <- vector("list", length = nrow(param_combos))
p_iso_vecs_from_test <- vector("list", length = nrow(param_combos))

p_post_iso_vecs_from_fx <- vector("list", length = nrow(param_combos))
p_post_iso_vecs_from_test <- vector("list", length = nrow(param_combos))

for (i in seq_len(nrow(param_combos))) {
  ## Get viral load curve time series from parameters that define triangle
  vl_curve <- isolation::triangle_vl(
    t = times,
    dp = param_combos$dp[i],
    tp = param_combos$tp[i],
    wp = param_combos$wp[i],
    wr = param_combos$wr[i]
  )

  ## Get time series for prob_antigen_positive
  ## there isn't a specific function for this, but can match the code in
  ## extract_posteriors_to_df.R

  prob_antigen_positive <- plogis(
    q = vl_curve,
    location = param_combos$antigen_50[i],
    scale = param_combos$sigma_antigen[i]
  )

  ## Get natural history (infectiousness)
  ## for the purpose of testing simulate_isolation_proposed, we don't test
  ## simulate_intrinsic_infectiousness; however it would be possible to test
  ## from the natural history parameters that form the triangle

  nat_hist_ts <- isolation::simulate_intrinsic_infectiousness(
    t = times,
    vl_curve = vl_curve,
    prob_antigen_positive = prob_antigen_positive,
    prob_culture_positive = prob_culture_positive,
    infec_threshold = param_combos$infec_threshold[i],
    transform_fun = param_combos$transform_fun[i]
  )

  ## Apply proposed guidance
  ## this generates the result we'll compare to in our test
  prop_list <- isolation::simulate_isolation_proposed(
    It = nat_hist_ts,
    t = times,
    symp_improv_in = param_combos$symp_improv_observed[i],
    symp_improv_adjust = param_combos$symp_improv_adjust[i],
    iso_lag = param_combos$iso_lag[i],
    iso_lag_adjust = param_combos$iso_lag_adjust[i],
    min_post_iso_days = param_combos$min_post_iso_days[i],
    iso_eff = param_combos$iso_eff[i],
    post_iso_eff = param_combos$post_iso_eff[i]
  )

  # store p_iso and p_post_iso

  p_iso_vecs_from_fx[[i]] <- prop_list$p_iso
  p_post_iso_vecs_from_fx[[i]] <- prop_list$p_post_iso

  # store summed I_t

  summed_I_t_from_fx[i] <- sum(prop_list$I_t)

  # store the resulting amount of infectiousness averted

  infec_averted_from_fx[i] <- prop_list$averted

  # the way "simulate_isolation_proposed" works is to calculate
  # probabilities that a person is isolating
  # or taking post-isolation precautions for each time t.

  # but alternatively and equivalently, we can think about epochs of relative
  # infectiousness (e.g., isolation period and post isolation period)
  # and calculate the amount of infectiousness averted during each

  # initialize a vector to store the amount of infectiousness averted

  infec_averted_pieces <- c()

  # initialize vectors of zeros for p_iso and p_post_iso

  p_iso <- rep(0, times = length(times))
  p_post_iso <- rep(0, times = length(times))

  # need to define when the person starts isolating relative to t = 0
  # where t = 0 is the end of the 24 hour period when symptoms started
  # if we assume that people with symptoms start isolating immediately
  # (i.e., iso_lag = 0), we model isolation starting at iso_lag - iso_lag_adjust
  # where iso_lag_adjust will typically be 0.5; we deal with the fact that
  # symptom timing observations are interval censored by taking a simple
  # mid-point approximation.

  iso_start <- param_combos$iso_lag[i] - param_combos$iso_lag_adjust[i]

  # similarly need to adjust the modeled observed symptom improvement time
  # for the data-generating mechanism of symptom data
  # suppose a person has an observed symptom improvement time of 4. This means
  # that improvement actually happened on Day 3 (e.g., if a person reported 2
  # symptoms in the 24 hour period ending at t = 4, and reported 3 symptoms
  # during the 24 hour period ending at t = 3, then the decrease in symptoms
  # must have occurred between t = 2 and t = 3,
  # and so is approximated as t = 2.5)

  symp_improv_actual <- param_combos$symp_improv_observed[i] -
    param_combos$symp_improv_adjust[i]

  # set the end point of isolation to be the symptom improvement time

  iso_end <- symp_improv_actual

  # for some people iso_start could be as long or longer than
  # symptom improvement; these people would never isolate, so
  # make logical indicator for whether isolation occurs
  # depends on iso_end being strictly greater than iso_start

  ever_isolate <- iso_end > iso_start

  if (ever_isolate) {
    which_iso_times <- times > iso_start & times < iso_end

    iso_infec_averted <- sum(nat_hist_ts[which_iso_times]) *
      param_combos$iso_eff[i]

    p_iso[which_iso_times] <- p_iso[which_iso_times] + 1

    which_iso_edges_times <- times %in% c(iso_start, iso_end)

    iso_infec_averted_edges <- sum(nat_hist_ts[which_iso_edges_times]) *
      0.5 * param_combos$iso_eff[i]

    p_iso[which_iso_edges_times] <- p_iso[which_iso_edges_times] + 0.5

    infec_averted_pieces <- c(
      infec_averted_pieces,
      iso_infec_averted,
      iso_infec_averted_edges
    )
  }

  # make logical indicator for post-isolation precautions
  # this differs from ever_isolate in that here it is sufficient for
  # iso_end to be greater than *or equal* to iso_start;
  # the logic here is that someone might realize that they have respiratory
  # virus symptoms *at the same time* those symptoms
  # improve and therefore might go straight to taking post-isolation
  # precautions; however, if symptom improvement time is short enough
  # and/or iso_lag is long enough that symptoms improve *before* iso_lag has
  # passed, then assume that no precautions are taken

  ever_post_isolate <- iso_end >= iso_start

  if (ever_post_isolate) {
    which_post_iso_times <- times > iso_end &
      times < iso_end + param_combos$min_post_iso_days[i]

    post_iso_infec_averted <- sum(nat_hist_ts[which_post_iso_times]) *
      param_combos$post_iso_eff[i]

    p_post_iso[which_post_iso_times] <- p_post_iso[which_post_iso_times] + 1

    which_post_iso_edges_times <-
      times %in% c(iso_end, iso_end + param_combos$min_post_iso_days[i])

    post_iso_infec_averted_edges <-
      sum(nat_hist_ts[which_post_iso_edges_times]) *
        0.5 * param_combos$post_iso_eff[i]

    p_post_iso[which_post_iso_edges_times] <-
      p_post_iso[which_post_iso_edges_times] + 0.5

    infec_averted_pieces <- c(
      infec_averted_pieces,
      post_iso_infec_averted,
      post_iso_infec_averted_edges
    )
  }

  # store p_iso and p_post_iso

  p_iso_vecs_from_test[[i]] <- p_iso
  p_post_iso_vecs_from_test[[i]] <- p_post_iso

  # calculate and store summed I_t

  summed_I_t_from_test[i] <- sum(nat_hist_ts) - sum(infec_averted_pieces)

  # sum up the amount of infectiousness averted under isolation
  # and post-isolation precautions, relative to infectiousness potential
  # in absence of intervention (natural history)

  cumulative_infec_averted <- sum(infec_averted_pieces) / sum(nat_hist_ts) * 100

  # store the result

  infec_averted_from_test[i] <- cumulative_infec_averted
}

test_that("simulate_isolation_proposed yields p_iso as expected ", {
  expect_equal(p_iso_vecs_from_fx, p_iso_vecs_from_test)
})

test_that("simulate_isolation_proposed yields p_post_iso as expected ", {
  expect_equal(p_post_iso_vecs_from_fx, p_post_iso_vecs_from_test)
})

test_that("simulate_isolation_proposed summed I_t is as expected ", {
  expect_equal(summed_I_t_from_fx, summed_I_t_from_test)
})

test_that("simulate_isolation_proposed averts expected infectiousness", {
  expect_equal(infec_averted_from_fx, infec_averted_from_test)
})
