#### Packages ####
list_of_packages = c("TropFishR", "dplyr", "Matrix", "here", "tools", "ggplot2", "patchwork")

new_packages = list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  install.packages(new_packages)
}

require(TropFishR)
require(dplyr)
require(Matrix)
require(here)
require(tools)
require(ggplot2)
require(patchwork)

## Source functions ------------------------------------------------------------
function_dir = here::here("R")
source_files = list.files(function_dir, pattern = "\\.R$", full.names = TRUE)
for (file in source_files) source(file)

### Load data ------------------------------------------------------------------
data_dir  = here::here("data")
files     = list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
data_list = lapply(files, read.csv) 
names(data_list) = tools::file_path_sans_ext(basename(files))

## Parameter info --------------------------------------------------------------
k_list    = c(0.25/365, 0.11/30, 0.513, 0.325/365, 0.61, 2.43/365, 1.35, 0.815)
linf_list = c(410, 145, NA, 336.85, NA, 87.27, NA, NA)
absolute  = ifelse(is.na(linf_list), TRUE, FALSE)

## Process data ----------------------------------------------------------------
# include dates
for (i in 1:length(data_list)) {
  data_list[[i]]$samp_date = as.Date(data_list[[i]]$samp_date, format = "%d.%m.%Y")
}

# cutting underrepresented lengths
data_list$CALSAP_PaP = data_list$CALSAP_PaP[data_list$CALSAP_PaP$Length <= 40, ]
data_list$CYNNEB_PaP = data_list$CYNNEB_PaP[data_list$CYNNEB_PaP$Length <= 100, ]
data_list$BAICHR_PaP = data_list$BAICHR_PaP[data_list$BAICHR_PaP$Length <= 70, ]

# create length frequency objects
lfqs = list()
bin_sizes = c(2,2,1,1,1,1,1,0.5)
for (i in 1:length(data_list)) {
  lfqs[[i]] = lfqCreate(data     = data_list[[i]],
                        Lname    = "Length",
                        Dname    = "samp_date",
                        bin_size = bin_sizes[i])
};names(lfqs) = names(data_list)

# create catch composition dataframes
catches = list()
mids    = list()

for (i in 1:length(lfqs)) {
  catches[[i]] = rowSums(lfqs[[i]]$catch)
  mids[[i]]    = lfqs[[i]]$midLengths
}; names(catches) = names(data_list); names(mids) = names(data_list)

dfs = list()
for (i in 1:length(catches)) {
  dfs[[i]] = data.frame(catch = catches[[i]], mids  = mids[[i]])
}; names(dfs) = names(data_list)

#dfs$LITSET_PaP$mids = dfs$LITSET_PaP$mids * 10 # convert to mm for LITSET

# pass data to Stan
N = rep(NA, length(dfs))
for (i in 1:length(dfs)) {
  N[i] = nrow(dfs[[i]])
}

sp_names = c(rep("ARIFEL", N[1]),
             rep("BAICHR", N[2]),
             rep("CALSAP", N[3]),
             rep("CYNARE", N[4]),
             rep("CYNNEB", N[5]),
             rep("FUNGRA", N[6]),
             rep("LAGRHO", N[7]),
             rep("LITSET", N[8]))

sp_index = as.numeric(factor(sp_names, levels = c("ARIFEL", "BAICHR", "CALSAP", "CYNARE", "CYNNEB", "FUNGRA", 
"LAGRHO", "LITSET")))

counts_concat = unlist(catches)
zero_index = which(counts_concat == 0) # remove zero counts
counts_concat = counts_concat[-zero_index] # All counts in one long vector
A = sapply(catches, length)     # Number of age bins per species
K = length(A)
start_idx = cumsum(c(1, head(A, -1)))
end_idx = cumsum(A)

# add relative age to each dataframe
rel_ages = list()
for (i in 1:length(dfs)) {
  if (absolute[i]) {
    rel_ages[[i]] = dfs[[i]]$mids / k_list[i]
  } else {
    rel_ages[[i]] = -(1/k_list[i]) * log(1 - (dfs[[i]]$mids / linf_list[i]))
  }
};names(rel_ages) = names(dfs)

# After removing zeros:
counts_concat = unlist(catches)
zero_index = which(counts_concat == 0)
counts_concat = counts_concat[-zero_index]

# same for rel_age
rel_age = as.vector(c(do.call(c, rel_ages)))
rel_age = rel_age[-zero_index]

# Now, recalculate A, start_idx, end_idx
A = sapply(catches, function(x) sum(x != 0)) # recalculated after zero removal!
K = length(A)
start_idx = cumsum(c(1, head(A, -1)))
end_idx = cumsum(A)

# temporary fix to test the model
rel_age[is.na(rel_age) | is.nan(rel_age)] = 470


# # create data list for Stan
data_list = list(N_total = length(counts_concat),
                 K = K,
                 A = A,
                 counts = counts_concat,
                 rel_age = rel_age,
                 start_idx = start_idx,
                 end_idx = end_idx)

# add priors for M and a50
data_list$M_prior_mean = c(catch_curve_summary_table$M_mean[1:8])
data_list$M_prior_sd = c(catch_curve_summary_table$M_sd[1:8] * 4) # multiply by 4 to widen the prior
data_list$a50_prior_mean = c(400, 80, 15, 120, 60, 90, 40, 10)
data_list$a50_prior_sd   = c(100, 10,   3,   5,  2, 90, 40,  2)

# fit model --------------------------------------------------------------------
rstan::stanc("stan/cc.stan")
stan_model = rstan::stan_model("stan/cc.stan") # compile model
n_chains = 1
fit = rstan::sampling(stan_model, 
                      data    = data_list,
                      iter    = 1000,
                      warmup  = 300,
                      chains  = n_chains,
                      seed    = 2027)

# verify convergence
rstan::check_hmc_diagnostics(fit)

# Plotting ---------------------------------------------------------------------
# 1. Selectivity curves
post <- as.data.frame(rstan::extract(fit, pars = c("k", "a50", "M")))

K <- length(grep("^k\\[", colnames(post))) # number of species

rel_age_df <- purrr::imap_dfr(rel_ages, function(rel, sp) {
  tibble(species = sp, rel_age = rel)
})

k_meds   <- apply(post[ , grep("^k\\.", names(post))], 2, median)
a50_meds <- apply(post[ , grep("^a50\\.", names(post))], 2, median)

rel_age_df <- rel_age_df %>%
  mutate(
    k = k_meds[species],
    a50 = a50_meds[species],
    selectivity = 1 / (1 + exp(-k * (rel_age - a50)))
  )

Sa = rep(NA, length(rel_age_df$rel_age))
for (i in 1:data_list$N_total) {
  
  for(k in 1:length(k_meds)) {
    
    Sa[i] = 1 / (1 + exp(-k_meds[k] * (rel_age_df$rel_age[i] - a50_meds[k])))
    
  }
}

rel_age_df$selectivity <- Sa

library(ggplot2)
ggplot(rel_age_df, aes(x = rel_age, y = selectivity, color = as.factor(species))) +
  geom_line(linewidth = 1) +
  labs(x = "Relative age", y = "Selectivity Sₐ", color = "Species") +
  custom_theme() +
  facet_wrap(~ species, scales = "free_x")

# 2. Observed vs. predicted catch compositions
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

K <- length(A)
n_draws <- 500  # Number of posterior draws to use for predictive intervals
draw_ids <- sample(nrow(post), n_draws)

pred_list <- vector("list", K)

for (s in 1:K) {
  idxs <- start_idx[s]:end_idx[s]
  n_bin <- length(idxs)
  rel_age_s <- rel_age[idxs]
  obs_s <- counts_concat[idxs]
  pred_mat <- matrix(NA, nrow = n_draws, ncol = n_bin)
  for (d in 1:n_draws) {
    k <- post[draw_ids[d], paste0("k.", s)]
    a50 <- post[draw_ids[d], paste0("a50.", s)]
    M <- post[draw_ids[d], paste0("M.", s)]
    Sa <- 1 / (1 + exp(-k * (rel_age_s - a50)))
    Na <- cumprod(c(1, rep(exp(-M), n_bin - 1)))
    pred <- Sa * Na
    pred_mat[d, ] <- pred / sum(pred)
  }
  # summarize posterior predictive mean and interval
  pred_df <- data.frame(
    bin = seq_along(idxs),
    rel_age = rel_age_s,
    obs_prop = obs_s / sum(obs_s),
    pred_mean = colMeans(pred_mat),
    pred_lower = apply(pred_mat, 2, quantile, 0.025),
    pred_upper = apply(pred_mat, 2, quantile, 0.975),
    species = s
  )
  pred_list[[s]] <- pred_df
}

# exclude species 1 first
pred_list = pred_list[-1]  # Exclude ARIFEL

all_pred_df = bind_rows(pred_list)

sp_labels = c("Silver perch", "Blue crab", "White trout", "Spotted seatrout",
  "Gulf killifish", "Pinfish", "White shrimp")

all_pred_df$species = factor(all_pred_df$species, 
                              levels = 2:8, 
                              labels = sp_labels)

ggplot(all_pred_df, aes(x = rel_age)) +
  geom_point(aes(y = obs_prop), color = "black", size = 2, pch = 20) +
  geom_line(aes(y = pred_mean), size = 1, color = "red") +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper), alpha = 0.15) +
  facet_wrap(~ species, scales = "free", nrow = 4) +
  labs(x = "Relative age (days)", y = "Proportion in catch") +
  custom_theme()

fig_dir = here::here("res", "figures")
ggsave(
  filename = file.path(fig_dir, "catch_comps.pdf"),
  plot = last_plot(),
  width = 4.5, height = 5.5, dpi = 300
)

# summary of M posterior distributions (means, sds, quantiles) for each species
catch_comp_summary_table = data.frame(
  sp = sp_labels,
  M_mean = sapply(2:8, function(s) mean(post[[paste0("M.", s)]])),
  M_sd = sapply(2:8, function(s) sd(post[[paste0("M.", s)]])),
  M_2.5 = sapply(2:8, function(s) quantile(post[[paste0("M.", s)]], probs = 0.025)),
  M_97.5 = sapply(2:8, function(s) quantile(post[[paste0("M.", s)]], probs = 0.975))
);rownames(catch_curve_summary_table) = NULL
print(catch_comp_summary_table)
