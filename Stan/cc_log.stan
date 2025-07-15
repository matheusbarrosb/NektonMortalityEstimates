data {
  int<lower=1> K;                  // number of species
  int<lower=1> N_total;            // total number of counts (sum over all species)
  array[K] int<lower=0> A;         // number of age bins per species
  array[N_total] int<lower=0> counts;    // concatenated observed counts
  array[K] int<lower=1> start_idx;       // start index for each species in counts/rel_age
  array[K] int<lower=1> end_idx;         // end index for each species in counts/rel_age
  vector[N_total] rel_age;               // concatenated relative ages for each count
  // mortality priors
  vector<lower=0>[K] M_prior_mean;
  vector<lower=0>[K] M_prior_sd;
  // selectivity priors
  vector<lower=0>[K] a50_prior_mean; 
  vector<lower=0>[K] a50_prior_sd; 
}
parameters {
  vector<lower=0>[K] k;    // selectivity steepness
  vector<lower=0>[K] a50;  // age at 50% selectivity
  vector<lower=0>[K] M;    // mortality
}
model {
  // Priors --------------------------------------------------------------------
  k   ~ normal(0, 10);
  a50 ~ normal(a50_prior_mean, a50_prior_sd);
  M   ~ normal(M_prior_mean, M_prior_sd);

  // Likelihood ----------------------------------------------------------------
  for (s in 1:K) {
    int n = A[s];
    vector[n] Sa;
    vector[n] log_Na;
    vector[n] log_pred;
    vector[n] pred;
    vector[n] p;

    // Selectivity (still in normal scale for interpretability)
    for (a in 1:n) {
      int idx = start_idx[s] + a - 1;
      Sa[a] = 1 / (1 + exp(-k[s] * (rel_age[idx] - a50[s])));
    }
    
    // Log numbers-at-age
    log_Na[1] = 0; // log(1)
    for (a in 2:n) {
      log_Na[a] = log_Na[a-1] - M[s];
    }
    // Log predicted catch
    for (a in 1:n) {
      log_pred[a] = log(Sa[a]) + log_Na[a];
    }
    // Predicted catch (exponentiate back)
    for (a in 1:n) {
      pred[a] = exp(log_pred[a]);
    }
    // Relative catch
    p = pred / sum(pred);

    counts[start_idx[s]:end_idx[s]] ~ multinomial(p);

  }
}