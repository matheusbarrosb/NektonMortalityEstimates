data {
  int<lower=1> K;
  int<lower=1> N_total;
  array[K] int<lower=1> A;
  array[N_total] int<lower=0> counts;
  array[K] int<lower=1> start_idx;
  array[K] int<lower=1> end_idx;
  vector[N_total] rel_age;
  vector<lower=0>[K] M_prior_mean;
  vector<lower=0>[K] M_prior_sd;
  vector<lower=0>[K] a50_prior_mean;
  vector<lower=0>[K] a50_prior_sd;
}
parameters {
  vector<lower=0>[K] k;
  vector<lower=0>[K] a50;
  vector<lower=0>[K] M;
  vector<lower=0>[K] alpha;
}
model {
  // Priors --------------------------------------------------------------------
  k     ~ normal(0, 2);
  a50   ~ normal(a50_prior_mean, a50_prior_sd);
  M     ~ normal(M_prior_mean, M_prior_sd);
  alpha ~ exponential(1);
  
  // Likelihood ----------------------------------------------------------------
  for (s in 1:K) {
    int n = A[s];
    vector[n] Sa;
    vector[n] Na;
    vector[n] pred;
    vector[n] p;
    int idx_start = start_idx[s];
    int idx_end = end_idx[s];
    //int counts_sub[n];
    array[n] int counts_sub;

    // Selectivity
    for (a in 1:n) {
      int idx = idx_start + a - 1;
      Sa[a] = 1 / (1 + exp(-k[s] * (rel_age[idx] - a50[s])));
    }
    // Relative numbers-at-age
    Na[1] = 1; 
    for (a in 2:n) {
      Na[a] = Na[a-1] * exp(-M[s]);
    }
    // Predicted catch
    for (a in 1:n) {
      pred[a] = Sa[a] * Na[a];
    }
    // Relative catch
    p = pred / sum(pred);

    // Build counts_sub array for this species
    for (a in 1:n) {
      counts_sub[a] = counts[idx_start + a - 1];
    }

    counts_sub ~ dirichlet_multinomial(alpha[s] * p);
  }
}

