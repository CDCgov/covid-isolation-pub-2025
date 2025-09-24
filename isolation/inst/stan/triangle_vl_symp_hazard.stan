// Model for simple triangle viral load
// Derived from Kissler et al. code
// https://github.com/skissler/Ct_SequentialInfections/blob/main/code/fit_posteriors.stan
// Simplified to not have covariates
// Addended with a model of symptom improvement deriving from clearance time
// and antigen positivity and culture positivity models

// the following function implements the triangle-shaped viral load curve
// and calculates viral load as a function of time
// tp is the time when the peak occurs
// wp is elapsed time from onset to the peak
// wr is elapsed time from peak to clearance
// dp is the magnitude of peak viral load
functions {
  // tri_vl does not need to be vectorized because it is the only call
  // (and needs to be called) at the individual-level bc of pooling
  real tri_vl(real t, real tp, real wp, real wr, real dp) {
    // model the "upslope" of viral proliferations
    if (t <= tp)
      return ((dp / wp) * (t - (tp - wp)));
    // model the "downslope" of viral clearance
    else
      return(dp - (dp / wr) * (t - tp));
  }

  // PDF of discrete weibull as written in R `DiscreteWeibull` package
  // for having the minimum support on days be days >= 1 (so no days = 0)
  // q ^ (days - 1) ^ beta - q ^ days ^ beta
  // we solve for q and beta in the traditional Weibull framework
  // to find that beta = si_shape
  // q = exp(-(1 / si_scale) ^ beta)
  // in practice, because we need to output log pmf, we need to calculate log(q)
  // log(q) = -(1 / si_scale) ^ beta
  // since it is only ever 1 / si_scale that appears, set si_scale_inv to be
  // 1 / si_scale and just define si_scale_inv as the parameter of interest
  // to prevent having to invert an exp(-something) expression
  // this is different from when fitting the overdispersion neg bin param
  // or exponential rate where the parameter should be of the form 1 / param
  // to reduce numerical instability challenges
  // here, si_scale is defined elsewhere and we are just transforming it to q

  // first we model the untruncated CDF because that is called in truncated function
  real discrete_weibull_lcdf(int max_si, real si_scale_inv, real beta) {
    // CDF = 1 - q ^ days ^ beta
    // lCDF = log1m_exp(log(q ^ days ^ beta))
    return log1m_exp(-max_si ^ beta * si_scale_inv ^ beta);
  }

  // truncated PDF
  real discrete_weibull_truncated_lpmf(int days, real si_scale_inv, real beta, int max_si) {
    // truncated PDF = PDF / CDF(max_si)
    // log(q ^ (days - 1) ^ beta - q ^ days ^ beta) - log(CDF(max_si))
    // see this first term as log of a difference --> use log_diff_exp stan function
    // log_diff_exp(log(q ^ (days - 1) ^ beta), log(q ^ days ^ beta))
    // log_diff_exp(log(q) * (days - 1) ^ beta, log(q) * days ^ beta)
    return log_diff_exp(-(days - 1) ^ beta * si_scale_inv ^ beta,
        -days ^ beta * si_scale_inv ^ beta) - discrete_weibull_lcdf(max_si | si_scale_inv, beta);
  }

  // for right-censoring of data, we need a complementary CDF (survival function)
  real discrete_weibull_truncated_lccdf(int days, real si_scale_inv, real beta, int max_si) {
    // truncated CDF = (1 - q ^ days ^ beta) / discreteweibull_cdf(max_si)
    // truncated CDF = (1 - q ^ days ^ beta) / (1 - q ^ max_si ^ beta)
    // cCDF = 1 - truncated CDF
    // log1m_exp(log(truncated CDF))
    return log1m_exp(discrete_weibull_lcdf(days | si_scale_inv, beta) - discrete_weibull_lcdf(max_si | si_scale_inv, beta));
  }

  // want to generate the discrete weibull symp improv times directly in stan
  real discrete_weibull_truncated_rng(real si_scale_inv, real beta, int max_si) {
    // use the inverse CDF theorem to get the discrete Weibull value
    // basically, pass a uniform random number through the discrete Weibull inverse CDF
    // take the ceiling of that output -- that is our simulated symptom improvement time
    // ceiling because the CDF tells us probability *up to and including* a given value
    real random_cdf_value = uniform_rng(0, 1); // called "u" throughout
    // discrete weibull CDF is 1 - q ^ days ^ beta
    // with truncation that becomes (1 - q ^ days ^ beta) / (1 - q ^ max_si ^ beta)
    // can easily see that this function hits 1 at max_si days
    // invert the CDF: u * (1 - q ^ max_si ^ beta) = 1 - q ^ days ^ beta
    // 1 - u * (1 - q ^ max_si ^ beta) = q ^ days ^ beta
    // log(1 - u * (1 - q ^ max_si ^ beta)) = days ^ beta * log(q)
    // days = (log(1 - u * (1 - q ^ max_si ^ beta)) / (log(q))) ^ (1 / beta)
    // plug in q = exp(-(1 / si_scale) ^ beta)
    real log_q = -(si_scale_inv ^ beta);
    real days = (log1m_exp(log(random_cdf_value) + discrete_weibull_lcdf(max_si | si_scale_inv, beta)) / log_q) ^ (1 / beta);
    return ceil(days);
  }

  // want to be able to test that our custom functions of prob_antigen_positive and prob_culture_positive
  // match with what we are modeling in stan
  real wrapped_bernoulli_logit_lpmf(int result, real x, real k_50, real scale) {
    return bernoulli_logit_lpmf(result | scale * (x - k_50));
  }
}

// the input data
data {
  // base facts about data
  int<lower=0> N;  // total number of observations (across all participants)
  int<lower=0> n_id; // number of study participants
  int<lower=0> id[N]; // vector of participant IDs for each observation
  // VL data
  real t[N]; // vector of time points of observations (wrt to symp onset time per person)
  real y[N]; // vector of observed log viral loads
  int<lower=0, upper=1> pcr_avail[N]; // whether a PCR test is available at that time for that person (not missing)
  real lod[N]; // by-sample (i.e., study-specific) limit of detection on log VL scale
  int<lower=0, upper=1> below_lod[N]; // whether the PCR test at that time is negative (below LOD)
  real loq_lower[N]; // by-sample (i.e., study-specific) lower bound limit of quantification on log VL scale
  real loq_upper[N]; //by-sample (i.e., study-specific) upper bound limit of quantification on log VL scale
  int<lower=0, upper=1> below_loq_lower[N]; // whether given PCR test is below lower bound LOQ (but above LOD)
  int<lower=0, upper=1> VL_within_detectable_range[N]; // whether given PCR test is above LOQ (and above LOD and not missing)
  int<lower=0, upper=1> above_loq_upper[N]; // whether given PCR test is above upper bound LOQ
  int<lower=0> antigen_result[N]; // results of antigen test
  int<lower=0, upper=1> antigen_avail[N]; // whether an antigen test is available at that time for that person (not missing)
  int<lower=0> culture_result[N]; // results of culture assay
  int<lower=0, upper=1> culture_avail[N]; // whether a culture test is available (not missing)
  int<lower=1> si_time[n_id]; // vector of observed symptom improvement times
  int<lower=0, upper=1> cens_duration[n_id]; // whether an individual's symptom improvement time is right censored
  int<lower=1> max_si; // maximum allowed time to symptom improvement

  // midpoints values to which everything is relative (mean prior for pooled param)
  real dp_midpoint; // prior mean peak viral load
  real<lower=0> wp_midpoint; // prior mean time from onset to peak (duration, so > 0)
  real<lower=0> wr_midpoint; // prior mean time from peak to clearance (duration, so > 0)

  // priors
  // all two dimensional priors are the mean and standard deviation of a normal
  // distribution on the parameter (potentially in log space)
  real tp_prior[2]; // prior of peak time wrt to symp onset
  real<lower=0> tp_std_prior[2]; // prior on individual-level variability in peak time
  real<lower=0> sigma_prior[2]; // prior observation noise (Cauchy scale)
  real<lower=0> prior_sd; // prior of individual-level variability in pooled params -- meshes with mean values above
  real antigen_50_prior[2]; // prior antigen 50% positive logVL (Cauchy scale)
  real culture_50_prior[2]; // prior antigen 50% positive logVL (at t = tp) (Cauchy scale)
  real culture_beta_prior[2]; // prior on coef on t for culture lod
  real<lower=0> si_shape_prior[2]; // prior Weibull shape param (recall si_shape = 1 is exponential)
  real si_beta_0_prior[2]; // prior base Weibull scale / Weibull hazard rate
  real si_beta_wr_prior[2]; // prior coef on clearance time (wr) towards symp improv
}

// variable declarations for parameters of the model
parameters {
  real tp_mean; // hierarchical mean of the peak time (pooled)
  real<lower=0> tp_std; // individual-level variability in tp
  real tp_raw[n_id]; // individual-level deviation in peak viral load time

  real log_dp_mean; // mean of peak viral load on the log scale (pooled)
  real<lower=0> log_dp_sd; // st dev of peak viral load on log scale across individuals
  real dp_raw[n_id]; // we want to fit peak viral load for each individual

  real log_wp_mean; // mean of proliferation time on log scale (pooled)
  real<lower=0> log_wp_sd; // st dev of proliferation time on log scale across individuals
  real wp_raw[n_id]; // we want to fit proliferation time for each individual

  real log_wr_mean; // mean of clearance time on log scale (pooled)
  real<lower=0> log_wr_sd; // st dev of clearance time on log scale across individuals
  real wr_raw[n_id]; // we want to fit clearance time for each individual

  real<lower=0> sigma; // observation noise (observed vs. expected viral load)

  real antigen_50; // logVL where antigen prob detection 50%
  real<lower=0> sigma_antigen; // slope/variability in antigen positivity about VL with 50% positivity

  real culture_50; // logVL value where culture prob detection 50% at t = tp (individual-level)
  real<lower=0> sigma_culture; // slope/variability in culture positivity about VL with 50% positivity
  real culture_beta; // coefficient on t for change in culture_50 by day

  real<lower=0> si_shape; // Weibull shape (recall si_shape = 1 --> exponential)
  real si_beta_0; // base Weibull scale (Weibull hazard rate)
  real si_beta_wr; // coefficient on individual-level deviation in clearance time (hazard function on wr for si)

}

// variable declarations and statements for transformed parameters
transformed parameters {

  // we need to make individual level params from the pooled params + std normal individual-level deviations
  real dp[n_id]; // individual-level peak log VL
  real wp[n_id]; // individual-level proliferation time
  real wr[n_id]; // individual-level clearance time
  real tp[n_id]; // we want to fit peak time for each individual -- different from Kissler
  real mu[N]; // "true" value of logVL at a given time according to triangle model
  real si_scale_inv[n_id]; // Weibull scale/e^hazard functions

  for (i in 1:n_id) { // loop through each participant
    // and construct their dp, wp, wr, and tp values from individual-level std normal deviations
    dp[i] = exp(log_dp_mean + log_dp_sd * dp_raw[i]) * dp_midpoint;
    wp[i] = exp(log_wp_mean + log_wp_sd * wp_raw[i]) * wp_midpoint;
    wr[i] = exp(log_wr_mean + log_wr_sd * wr_raw[i]) * wr_midpoint;
    tp[i] = tp_mean + tp_std * tp_raw[i]; // tp goes positive and negative, so not exponentiated

    // symptom improvement Weibull scale term
    // as a proportional hazards function from clearance time
    si_scale_inv[i] = exp(-(si_beta_0 + si_beta_wr * wr_raw[i]));
  }

  for (i in 1:N) {
    // expected logVL according to triangle model
    // takes time as a variable input to determine where along triangle we are
    // and individual-level VL params define the shape of the triangle
    mu[i] = tri_vl(t[i], tp[id[i]], wp[id[i]], wr[id[i]], dp[id[i]]);
  }
}

// all the terms for which the log probability should be accumulated
model {

  tp_mean ~ normal(tp_prior[1], tp_prior[2]); // hierarchical mean
  tp_std ~ normal(tp_std_prior[1], tp_std_prior[2]) T[0,]; // sd of individual-level draws
  tp_raw ~ std_normal(); // draw an individual-level guess

  log_dp_mean ~ normal(0, prior_sd); // hierarchical mean
  log_dp_sd ~ normal(0, prior_sd) T[0,]; // sd of individual-level draws
  dp_raw ~ std_normal(); // draw an individual-level guess

  log_wp_mean ~ normal(0, prior_sd); // hierarchical mean
  log_wp_sd ~ normal(0, prior_sd) T[0,]; // sd of individual-level draws
  wp_raw ~ std_normal(); // draw an individual-level guess

  log_wr_mean ~ normal(0, prior_sd); // hierarchical mean
  log_wr_sd ~ normal(0, prior_sd) T[0,]; // sd of individual-level draws
  wr_raw ~ std_normal(); // draw an individual-level guess

  sigma ~ normal(sigma_prior[1], sigma_prior[2]) T[0,]; // measurement noise in RT-PCR log VL

  // antigen logistic regression params
  antigen_50 ~ normal(antigen_50_prior[1], antigen_50_prior[2]);
  sigma_antigen ~ normal(sigma_prior[1], sigma_prior[2]) T[0,]; // truncate at 0 to enforce VL incr --> prob antigen positive incr

  // culture logistic regression params
  culture_50 ~ normal(culture_50_prior[1], culture_50_prior[2]);
  sigma_culture ~ normal(sigma_prior[1], sigma_prior[2]) T[0,]; // same as for antigen
  culture_beta ~ normal(culture_beta_prior[1], culture_beta_prior[2]); // daily (could be + or -) rate at which culture positivity log odds drops

  si_shape ~ normal(si_shape_prior[1], si_shape_prior[2]) T[0,]; // Weibull shape param > 0 always
  // these next two params get exponentially transformed to make Weibull scale
  si_beta_0 ~ normal(si_beta_0_prior[1], si_beta_0_prior[2]); // baseline Weibull hazard/failure "rate", can be + or -
  si_beta_wr ~ normal(si_beta_wr_prior[1], si_beta_wr_prior[2]); // clearance time hazard rate, can be + or -

  for (i in 1:N) {
    // only accumulate probability for logVL if the test result is available
    if (pcr_avail[i] == 1)

    // in all cases, PCR measurements have some measurement noise, sigma, from mu
      if (below_lod[i] == 1)
        // COVID not detected means logVL must be below LOD -- interval censoring
        target += normal_lcdf(lod[i] | mu[i], sigma);
      else if (below_loq_lower[i] == 1)
        // means that logVL must be between LOD and loq -- interval censoring
        // phi(loq) - phi(lod)
        target += log_diff_exp(normal_lcdf(loq_lower[i] | mu[i], sigma), normal_lcdf(lod[i] | mu[i], sigma));
      else if (VL_within_detectable_range[i] == 1)
        // logVL measurement is subject to measurement noise sigma wrt to true value mu
        target += normal_lupdf(y[i] | mu[i], sigma);
      else if (above_loq_upper[i] == 1)
        // logVL is above the upper bound LOQ (above the standard curve's quantification range)
        // so we know LOQ is above some given value, but don't know where because the standard breaks down
        target += normal_lccdf(loq_upper[i] | mu[i], sigma);

    // both antigen and culture results are either negative or positive
    // basically encode this as a logistic regression as a function of VL: log-odds of testing positive
    // scale linearly with logVL, and for culture data time since peak
    if (antigen_avail[i] == 1)
      target += wrapped_bernoulli_logit_lpmf(antigen_result[i] | mu[i], antigen_50, 1 / sigma_antigen);

    if (culture_avail[i] == 1)
      target += wrapped_bernoulli_logit_lpmf(culture_result[i] | mu[i], culture_50 + culture_beta * (t[i] - tp[id[i]]), 1 / sigma_culture);
  }

  for (i in 1:n_id) { // symptom improvement is at the participant level
    // recall we calculated si_scale_inv based on wr_raw in transformed parameters
    if (cens_duration[i] == 1) // if the person has not improved by the time for which we have their last data
      // we know improvement must have happened some time after the last observed symptom data time
      target += discrete_weibull_truncated_lccdf(si_time[i] | si_scale_inv[i], si_shape, max_si);
    else
      // otherwise, probability of observing that improvement time based on participant scale and category shape
      target += discrete_weibull_truncated_lpmf(si_time[i] | si_scale_inv[i], si_shape, max_si);
  }
}

generated quantities {
  real si_simulated[n_id];
  for (i in 1:n_id) {
    si_simulated[i] = discrete_weibull_truncated_rng(si_scale_inv[i], si_shape, max_si);
  }
}
