# test simulate_isolation_current.R

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
symp_type_cat <- c(1, 3) # include "moderate" and "mild" categories

# set levels of intervention parameters

iso_lag <- c(0, 2)
test_count <- 2
test_interval <- 2
min_iso_days_moderate <- 10
min_iso_days_mild <- 5
min_total_precaution_days <- 10
iso_eff <- 0.9
post_iso_eff <- 0.5
symp_improv_adjust <- 1.5
iso_lag_adjust <- 0.5

# create all combinations of parameters
param_combos <- tidyr::expand_grid(
  tp, wp, wr, dp,
  symp_improv_observed,
  antigen_50, sigma_antigen,
  transform_fun,
  symp_type_cat,
  infec_threshold, iso_lag, test_count, test_interval,
  min_iso_days_moderate, min_iso_days_mild,
  min_total_precaution_days,
  iso_eff, post_iso_eff,
  symp_improv_adjust, iso_lag_adjust
)

# define the time series that will be used across all simulated persons

times <- seq(from = -10, to = 21, by = 0.1)

# set up vector to store probability antigen detection

ant_det_vec <- rep(NA, times = nrow(param_combos))

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
  ## for the purpose of testing simulate_isolation_current, we don't test
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

  ## Apply current guidance
  ## this generates the result we'll compare to in our test
  curr_list <- isolation::simulate_isolation_current(
    It = nat_hist_ts,
    t = times,
    prob_antigen_positive = prob_antigen_positive,
    symp_type_cat = param_combos$symp_type_cat[i],
    symp_improv = param_combos$symp_improv_observed[i],
    iso_lag = param_combos$iso_lag[i],
    test_count = param_combos$test_count[i],
    test_interval = param_combos$test_interval[i],
    min_iso_days_moderate = param_combos$min_iso_days_moderate[i],
    min_iso_days_mild = param_combos$min_iso_days_mild[i],
    min_total_precaution_days = param_combos$min_total_precaution_days[i],
    iso_eff = param_combos$iso_eff[i],
    post_iso_eff = param_combos$post_iso_eff[i],
    symp_improv_adjust = param_combos$symp_improv_adjust[i],
    iso_lag_adjust = param_combos$iso_lag_adjust[i]
  )

  # store p_iso and p_post_iso

  p_iso_vecs_from_fx[[i]] <- curr_list$p_iso
  p_post_iso_vecs_from_fx[[i]] <- curr_list$p_post_iso

  # store summed I_t

  summed_I_t_from_fx[i] <- sum(curr_list$I_t)

  # store the resulting amount of infectiousness averted

  infec_averted_from_fx[i] <- curr_list$averted

  # the way "simulate_isolation_current" works is to calculate
  # probabilities that a person is isolating
  # or taking post-isolation precautions for each time t.

  # but alternatively and equivalently, we can think about epochs of relative
  # infectiousness (e.g., isolation period and post isolation period)
  # and calculate the amount of infectiousness averted during each
  # this is more complicated for "simulate_isolation_current" than for
  # "simulate_isolation_proposed" because in the former, we average across
  # the probability of detection by antigen testing.

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

  # first consider an epoch of isolation that begins with symptom onset
  # (plus potentially a non-zero lag to enter isolation)
  # and ends with symptom improvement
  # this epoch is not dependent on test results
  # however, it does not happen for everyone; some people have such short
  # symptom duration or such long lag to enter isolation that
  # symp_improv_actual is not strictly greater than iso_start; make a logical
  # indicator to capture this

  ever_symp_isolate <- symp_improv_actual > iso_start

  if (ever_symp_isolate) {
    which_iso_symp_times <-
      which(times > iso_start & times < symp_improv_actual)

    iso_symp_infec_averted <-
      sum(nat_hist_ts[which_iso_symp_times]) * param_combos$iso_eff[i]

    p_iso[which_iso_symp_times] <- p_iso[which_iso_symp_times] + 1

    which_iso_symp_edges_times <-
      which(times %in% c(iso_start, symp_improv_actual))

    iso_symp_infec_averted_edges <-
      sum(nat_hist_ts[which_iso_symp_edges_times]) *
        0.5 * param_combos$iso_eff[i]

    p_iso[which_iso_symp_edges_times] <- p_iso[which_iso_symp_edges_times] + 0.5

    infec_averted_pieces <- c(
      infec_averted_pieces,
      iso_symp_infec_averted,
      iso_symp_infec_averted_edges
    )
  }

  # now consider epochs that are dependent on a positive test result
  # need to calculate probability of getting at least one positive result

  ## times when a person decides to test
  test_t <- seq(
    from = param_combos$iso_lag[i],
    to = param_combos$iso_lag[i] +
      (param_combos$test_count[i] - 1) * param_combos$test_interval[i],
    by = param_combos$test_interval[i]
  )

  ## probability of antigen detection based on times of testing
  ant_det <- get_antigen_prob_detection(
    test_t = test_t,
    t = times,
    antigen_sens = prob_antigen_positive,
    symp_improv = param_combos$symp_improv_observed[i],
    symp_improv_adjust = param_combos$symp_improv_adjust[i],
    iso_lag_adjust = param_combos$iso_lag_adjust[i]
  )

  ant_det_vec[i] <- ant_det

  # depending on at least one positive antigen test, there may be a period of
  # extended isolation precautions because of the minimum duration of isolation
  # component of the guidance

  # the minimum duration of isolation depends on the symptom category

  if (param_combos$symp_type_cat[i] == 1) {
    min_iso_days <- param_combos$min_iso_days_moderate[i]
  } else {
    min_iso_days <- param_combos$min_iso_days_mild[i]
  }

  # get indicator for extended isolation; depends on min_iso_days being longer
  # than when symptoms improve and being later than the lag to start isolation
  # also, iso_start must be before or at the same time as symptom improvement
  # (which corresponds to when testing may happen relative to symptoms)

  ever_post_symp_isolate <- min_iso_days > symp_improv_actual &
    min_iso_days > iso_start &
    iso_start <= symp_improv_actual

  if (ever_post_symp_isolate) {
    which_iso_post_symp_times <-
      which(times > symp_improv_actual & times < min_iso_days)

    iso_post_symp_infec_averted <-
      sum(nat_hist_ts[which_iso_post_symp_times]) *
        param_combos$iso_eff[i] * ant_det

    p_iso[which_iso_post_symp_times] <-
      p_iso[which_iso_post_symp_times] + ant_det

    which_iso_post_symp_edges_times <-
      which(times %in% c(symp_improv_actual, min_iso_days))

    iso_post_symp_infec_averted_edges <-
      sum(nat_hist_ts[which_iso_post_symp_edges_times]) *
        0.5 * param_combos$iso_eff[i] * ant_det

    p_iso[which_iso_post_symp_edges_times] <-
      p_iso[which_iso_post_symp_edges_times] + 0.5 * ant_det

    infec_averted_pieces <- c(
      infec_averted_pieces,
      iso_post_symp_infec_averted,
      iso_post_symp_infec_averted_edges
    )
  }

  # make logical indicator for post-isolation precautions
  # for post-isolation precautions to occur, isolation itself must end before
  # minimum total precaution days is reached
  # also, iso_start has to be less than or equal symptom improvement
  # (or else won't test)

  ever_post_isolate <-
    param_combos$min_total_precaution_days[i] > symp_improv_actual &
      param_combos$min_total_precaution_days[i] > min_iso_days &
      iso_start <= symp_improv_actual

  if (ever_post_isolate) {
    which_post_iso_times <-
      times > max(symp_improv_actual, min_iso_days) &
        times < param_combos$min_total_precaution_days[i]

    post_iso_infec_averted <-
      sum(nat_hist_ts[which_post_iso_times]) *
        param_combos$post_iso_eff[i] * ant_det

    p_post_iso[which_post_iso_times] <-
      p_post_iso[which_post_iso_times] + ant_det

    which_post_iso_edges_times <-
      which(times %in% c(
        max(symp_improv_actual, min_iso_days),
        param_combos$min_total_precaution_days[i]
      ))

    post_iso_infec_averted_edges <-
      sum(nat_hist_ts[which_post_iso_edges_times]) *
        0.5 * param_combos$post_iso_eff[i] * ant_det

    p_post_iso[which_post_iso_edges_times] <-
      p_post_iso[which_post_iso_edges_times] + 0.5 * ant_det

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

test_that("simulate_isolation_current yields p_iso as expected ", {
  expect_equal(p_iso_vecs_from_fx, p_iso_vecs_from_test)
})

test_that("simulate_isolation_current yields p_post_iso as expected ", {
  expect_equal(p_post_iso_vecs_from_fx, p_post_iso_vecs_from_test)
})

test_that("simulate_isolation_current summed I_t is as expected ", {
  expect_equal(summed_I_t_from_fx, summed_I_t_from_test)
})

test_that("simulate_isolation_current averts expected infectiousness", {
  expect_equal(infec_averted_from_fx, infec_averted_from_test)
})

# the following dataframe may be helpful for troubleshooting

combos_with_results <- param_combos |>
  dplyr::mutate(ant_det_vec,
    summed_I_t_from_fx, summed_I_t_from_test,
    "rel_diff_summed_I_f" =
      (summed_I_t_from_test - summed_I_t_from_fx) /
        summed_I_t_from_fx,
    infec_averted_from_fx, infec_averted_from_test,
    "rel_diff_infec_averted" =
      (infec_averted_from_test - infec_averted_from_fx) /
        infec_averted_from_fx
  )

# test for errors when symp_type_cat is not in 1 - 4.

test_that(
  "Errors for categories not in 1 - 4",
  expect_error(simulate_isolation_current(
    It = nat_hist_ts,
    t = times,
    prob_antigen_positive = prob_antigen_positive,
    symp_type_cat = 5,
    symp_improv = param_combos$symp_improv_observed[i],
    iso_lag = param_combos$iso_lag[i],
    test_count = param_combos$test_count[i],
    test_interval = param_combos$test_interval[i],
    min_iso_days_moderate = param_combos$min_iso_days_moderate[i],
    min_iso_days_mild = param_combos$min_iso_days_mild[i],
    min_total_precaution_days = param_combos$min_total_precaution_days[i],
    iso_eff = param_combos$iso_eff[i],
    post_iso_eff = param_combos$post_iso_eff[i],
    symp_improv_adjust = param_combos$symp_improv_adjust[i],
    iso_lag_adjust = param_combos$iso_lag_adjust[i]
  ))
)
