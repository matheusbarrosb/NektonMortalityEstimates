data {
  int<lower=1> K;                  // number of species
  int<lower=1> N_total;            // total number of counts (sum over all species)
  array[K] int<lower=0> A;
  //int<lower=1> A[K];               // number of age bins per species
  array[N_total] int<lower=0> counts;    // concatenated observed counts
  array[K] int<lower=1> start_idx;
  array[K] int<lower=1> end_idx;
  vector[N_total] rel_age;         // concatenated relative ages for each count
  // mortality priors
  vector<lower=0>[K] M_prior_mean;
  vector<lower=0>[K] M_prior_sd;
  // selectivity priors
  vector<lower=0>[K] a50_prior_mean; 
  vector<lower=0>[K] a50_prior_sd; 
}
parameters {
  vector<lower=0>[K] k;
  vector<lower=0>[K] a50;
  vector<lower=0>[K] M;
}
model {
  // Priors --------------------------------------------------------------------
  k   ~ exponential(1);
  a50 ~ normal(a50_prior_mean, a50_prior_sd);
  M   ~ normal(M_prior_mean, M_prior_sd);


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
      int idx = start_idx[s] + a - 1;
      Na[a] = Na[a-1] * exp(-M[s]);
    }
    // Predicted catch
    for (a in 1:n) {
      pred[a] = Sa[a] * Na[a];
    }
    // Relative catch
    p = pred / sum(pred);

    counts[start_idx[s]:end_idx[s]] ~ multinomial(p);
  }
}
