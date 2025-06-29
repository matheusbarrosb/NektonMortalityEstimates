data {
  int<lower=1> K;                  // number of species
  int<lower=1> N_total;            // total number of counts (sum over all species)
  int<lower=1> A[K];               // number of age bins per species
  int<lower=0> counts[N_total];    // concatenated observed counts
  int<lower=1> start_idx[K];       // start index for each species in counts
  int<lower=1> end_idx[K];         // end index for each species in counts
  vector[N_total] rel_age;         // concatenated relative ages for each count
  vector<lower=0>[K] M_prior_mean; // mortality priors
  vector<lower=0>[K] a50_prior_mean; // selectivity priors
}
parameters {
  vector<lower=0>[K] k;
  vector<lower=0>[K] a50;
  vector<lower=0>[K] M;
}
model {
  // Priors --------------------------------------------------------------------
  k   ~ normal(0, 10);
  a50 ~ normal(a50_prior_mean, 10);
  M   ~ normal(M_prior_mean, 0.01);

  // Likelihood ----------------------------------------------------------------
  for (s in 1:K) {
    int n = A[s];
    vector[n] Sa;
    vector[n] Na;
    vector[n] pred;
    vector[n] p;

    // Selectivity
    for (a in 1:n) {
      int idx = start_idx[s] + a - 1;
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
    p = pred / sum(pred);

    counts[start_idx[s]:end_idx[s]] ~ multinomial(p);
  }
}
