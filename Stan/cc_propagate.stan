data {
  int<lower=1> K;                        // number of species
  int<lower=1> N_total;                  // total number of fish across all species
  array[K] int<lower=1> A;               // number of age bins per species
  array[N_total] int<lower=0> counts;    // concatenated observed counts (by age bin)
  array[K] int<lower=1> start_idx;       // start index for each species in counts (1-based)
  array[K] int<lower=1> end_idx;         // end index for each species in counts (1-based)
  vector[N_total] length_obs;            // observed lengths for each fish 
  array[K] int<lower=1, upper=2> growth_model; // 1=VB, 2=linear
  // Priors for mortality
  vector<lower=0>[K] M_prior_mean;
  vector<lower=0>[K] M_prior_sd;
  // Priors for selectivity
  vector<lower=0>[K] a50_prior_mean; 
  vector<lower=0>[K] a50_prior_sd; 
  // Priors for VB growth
  vector<lower=0>[K] Linf_prior_mean;
  vector<lower=0>[K] Linf_prior_sd;
  vector<lower=0>[K] k_vb_prior_mean;
  vector<lower=0>[K] k_vb_prior_sd;
  vector<lower=0>[K] Linf_min;
  // Priors for linear growth
  vector<lower=0>[K] g_lin_prior_mean;
  vector<lower=0>[K] g_lin_prior_sd;
}
parameters {
  // Selectivity params
  vector<lower=0>[K] k_sel;    // selectivity steepness
  vector<lower=0>[K] a50;      // selectivity a50
  // Mortality
  vector<lower=0>[K] M;
  // VB Growth
  vector<lower=Linf_min>[K] Linf;     // von Bertalanffy Linf
  vector<lower=0>[K] k_vb;     // von Bertalanffy k
  vector[K] t0;                // von Bertalanffy t0
  // Linear Growth
  vector<lower=0>[K] g_lin;    // linear growth rate
  vector[K] L0_lin;            // linear length at age 0
}
model {
  // Priors
  k_sel   ~ normal(0, 2);
  a50     ~ normal(a50_prior_mean, a50_prior_sd);
  M       ~ normal(M_prior_mean, M_prior_sd);
  Linf    ~ normal(Linf_prior_mean, Linf_prior_sd);
  k_vb    ~ normal(k_vb_prior_mean, k_vb_prior_sd);
  g_lin   ~ normal(g_lin_prior_mean, g_lin_prior_sd);

  // Likelihood
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
    vector[n] rel_age_s;

    // Compute relative ages for this species with appropriate growth model
    for (a in 1:n) {
      int idx = idx_start + a - 1;
      if (growth_model[s] == 1) {
        // VB
        rel_age_s[a] = - (1/k_vb[s]) * log1m(length_obs[idx]/Linf[s]); //  -(1/k_list[i]) * log(1 - (dfs[[i]]$mids / linf_list[i]))
      } else {
        // Linear 
        rel_age_s[a] = length_obs[idx]/g_lin[s];
      }
    }

    // Selectivity-at-age
    for (a in 1:n) {
      Sa[a] = 1 / (1 + exp(-k_sel[s] * (rel_age_s[a] - a50[s])));
    }
    // Numbers-at-age (relative)
    Na[1] = 1;
    for (a in 2:n) {
      Na[a] = Na[a-1] * exp(-M[s]);
    }
    // Predicted catch
    for (a in 1:n) {
      pred[a] = Sa[a] * Na[a];
    }
    p = pred / sum(pred);

    // Extract counts for this species
    for (a in 1:n) {
      counts_sub[a] = counts[idx_start + a - 1];
    }
    counts_sub ~ multinomial(p);
  }
}