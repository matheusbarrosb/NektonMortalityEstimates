#### Packages ####
list_of_packages = c("TropFishR", "dplyr", "Matrix", "here", "tools", "ggplot2", "patchwork", "cmdstanr")
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
require(cmdstanr) 
require(purrr)
require(cmdstanr)

## Source functions ------------------------------------------------------------
function_dir = here::here("R")
source_files = list.files(function_dir, pattern = "\\.R$", full.names = TRUE)
for (file in source_files) source(file)

### Load data ------------------------------------------------------------------
data_dir  = here::here("data")
files     = list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
# exclude 'ShrimpRegression.csv' file before processing
files = files[!grepl("ShrimpRegression.csv", files)]
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

# convert from CL to TL for white shrimp
data_list$LITSET_PaP$Length = data_list$LITSET_PaP$Length * 4.781 + 5.847
scal = .24
# cutting underrepresented lengths
data_list$CALSAP_PaP = data_list$CALSAP_PaP[data_list$CALSAP_PaP$Length <= 40, ]
data_list$CYNNEB_PaP = data_list$CYNNEB_PaP[data_list$CYNNEB_PaP$Length <= 100, ]
data_list$BAICHR_PaP = data_list$BAICHR_PaP[data_list$BAICHR_PaP$Length <= 70, ]
#data_list$LITSET_PaP = data_list$LITSET_PaP[data_list$LITSET_PaP$Length <= 100, ]

# create length frequency objects
lfqs = list()
bin_sizes = c(2,2,1,1,1,1,1,1)
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
zero_index    = which(counts_concat == 0) # remove zero counts
counts_concat = counts_concat[-zero_index] # All counts in one long vector
A             = sapply(catches, length)     # Number of age bins per species
K             = length(A)
start_idx     = cumsum(c(1, head(A, -1)))
end_idx       = cumsum(A)

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

rel_age[is.na(rel_age) | is.nan(rel_age)] = 470

# also concatenate mids like counts
length_obs = unlist(mids)
length_obs = length_obs[-zero_index] # remove zero counts

# # create data list for Stan
data_list = list(N_total = length(counts_concat),
                 length_obs = length_obs,
                 growth_model = c(1, 1, 2, 1, 2 ,1, 2, 2), # 1 VB, 2 linear 
                 K = K,
                 A = A,
                 counts = counts_concat,
                 rel_age = rel_age,
                 start_idx = start_idx,
                 end_idx = end_idx)

# add priors and hyperpriors
data_list$M_prior_mean = c(catch_curve_summary_table$M_mean[1:8])
data_list$M_prior_sd = c(catch_curve_summary_table$M_sd[1:8] * 4) # multiply by 4 to widen the prior
# data_list$M_prior_sd[6] = 0.1
# data_list$M_prior_sd[8] = 0.0005
data_list$a50_prior_mean = c(400, 80, 15, 120, 60, 80, 40, 60)
data_list$a50_prior_sd   = c(100, 1,   3,   5,  2, 5, 40,  100)

# fit model --------------------------------------------------------------------
stan_model_1 = cmdstanr::cmdstan_model("stan/cc.stan")
stan_model_2 = cmdstanr::cmdstan_model("stan/ccdm.stan")
stan_model_3 = cmdstanr::cmdstan_model("stan/cc_log.stan")

n_chains = 4

fit = stan_model_1$sample(
  data          = data_list,
  iter_sampling = 10000,
  iter_warmup   = 2000,
  chains        = n_chains,
  seed          = 2027
)


if (all(fit$summary(variables = c("a50", "M"))$rhat < 1.01)) {
  message("All parameters converged successfully.")
} else {
  message("Some parameters did not converge. Check the Rhat values.")
}

fit$summary(variables = c("a50", "M")) %>%
  dplyr::select(variable, mean, sd, rhat, ess_bulk, ess_tail)

# Produce convergence plots for each M, a50, and k from cmdstanr fit with ggplot2
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)

# Extract draws for M, a50, and k
post = fit$draws(variables = c("M", "a50", "k"), format = "df")

# Plotting MCMC traceplots -----------------------------------------------------
post %>%
  pivot_longer(cols = everything(), names_to = c(".value", "species"), names_pattern = "(.*)\\[(\\d+)\\]") %>%
  mutate(species = as.numeric(species)) %>% mutate(iteration = row_number() %% nrow(post) + 1) %>%
  # add chain number for every 10K iterations for each species, there were only four chains so it has to be 1 - 4
  mutate(chain = (iteration - 1) %/% 10000 + 1) %>%
  filter(!is.na(M) & !is.na(a50) & !is.na(k)) %>% 
  filter(iteration %% 10 == 0) %>%
  filter(species != 1) %>%
  pivot_longer(cols = c(M), names_to = "parameter") %>%
  mutate(species = factor(species, 
                          levels = 2:8, 
                          labels = c("Silver perch", "Blue crab", "White trout", "Spotted seatrout",
                                     "Gulf killifish", "Pinfish", "White shrimp"))) %>%
  mutate(value = ifelse(species == "White shrimp", value * scal, value)) %>%
  ggplot(aes(x = iteration, y = value, color = as.factor(chain))) +
  geom_line() +
  facet_wrap(~species, scales = "free_y", nrow = 4) +
  labs(x = "Iteration",
       y = expression(M~"(day"^{-1}*")"),
       color = "Chain") +
  custom_theme() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = PNWColors::pnw_palette("Bay", n = 4))

ggsave(
  filename = file.path(here::here("res", "figures"), "M_traceplots.pdf"),
  plot = last_plot(),
  width = 6, height = 5, dpi = 300
)

post %>%
  pivot_longer(cols = everything(), names_to = c(".value", "species"), names_pattern = "(.*)\\[(\\d+)\\]") %>%
  mutate(species = as.numeric(species)) %>% mutate(iteration = row_number() %% nrow(post) + 1) %>%
  mutate(species = factor(species, 
                          levels = 2:8, 
                          labels = c("Silver perch", "Blue crab", "White trout", "Spotted seatrout",
                                     "Gulf killifish", "Pinfish", "White shrimp"))) %>%
  # add chain number for every 10K iterations for each species, there were only four chains so it has to be 1 - 4
  mutate(chain = (iteration - 1) %/% 10000 + 1) %>%
  filter(!is.na(a50)) %>%
  filter(iteration %% 20 == 0) %>%
  filter(species != 1) %>%
  
  ggplot(aes(x = iteration, y = a50, color = as.factor(chain))) +
  geom_line() +
  facet_wrap(~species, scales = "free_y", nrow = 4) +
  labs(x = "Iteration",
       y = expression(a[50]~"(days)"),
       color = "Chain") +
  custom_theme() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = PNWColors::pnw_palette("Bay", n = 4))

ggsave(
  filename = file.path(here::here("res", "figures"), "a50_traceplots.pdf"),
  plot = last_plot(),
  width = 6, height = 5, dpi = 300
)

post %>%
  pivot_longer(cols = everything(), names_to = c(".value", "species"), names_pattern = "(.*)\\[(\\d+)\\]") %>%
  mutate(species = as.numeric(species)) %>% mutate(iteration = row_number() %% nrow(post) + 1) %>%
  mutate(species = factor(species, 
                          levels = 2:8, 
                          labels = c("Silver perch", "Blue crab", "White trout", "Spotted seatrout",
                                     "Gulf killifish", "Pinfish", "White shrimp"))) %>%
  # add chain number for every 10K iterations for each species, there were only four chains so it has to be 1 - 4
  mutate(chain = (iteration - 1) %/% 10000 + 1) %>%
  filter(!is.na(k)) %>%
  filter(iteration %% 2 == 0) %>%
  filter(species != 1) %>%
  
  ggplot(aes(x = iteration, y = k, color = as.factor(chain))) +
  geom_line() +
  facet_wrap(~species, scales = "free_y", nrow = 4) +
  labs(x = "Iteration",
       y = expression(k~"(day"^{-1}*")"),
       color = "Chain") +
  custom_theme() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = PNWColors::pnw_palette("Bay", n = 4))

ggsave(
  filename = file.path(here::here("res", "figures"), "k_traceplots.pdf"),
  plot = last_plot(),
  width = 6, height = 5, dpi = 300
)

# Plotting potential scale reduction factor
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)
# Extract the Rhat values for M, a50, and k
fit$summary(variables = c("M", "a50", "k")) %>%
  select(variable, rhat) %>%
  mutate(species = as.numeric(gsub(".*\\[(\\d+)\\]", "\\1", variable))) %>%
  filter(species != 1) %>%
  mutate(species = factor(species, 
                          levels = 2:8, 
                          labels = c("Silver perch", "Blue crab", "White trout", "Spotted seatrout",
                                     "Gulf killifish", "Pinfish", "White shrimp"))) %>%
  mutate(parameter = gsub("\\[.*\\]", "", variable)) %>%
  mutate(parameter = ifelse(parameter == "k", "s", parameter)) %>%
  filter(!is.na(rhat)) %>%
  ggplot(aes(x = species, y = rhat, color = parameter)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "Species", y = expression(R^{"^"})) +
  custom_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(), legend.position = "right") +
  ylim(0.9995, 1.001)

ggsave(
  filename = file.path(here::here("res", "figures"), "rhat.pdf"),
  plot = last_plot(),
  width = 4, height = 3, dpi = 300
)

# Plotting ---------------------------------------------------------------------
# 1. Selectivity curves
post = fit$draws(variables = c("k", "a50", "M"), format = "df")

K = data_list$K

rel_age_df = purrr::imap_dfr(rel_ages, function(rel, sp) {
  tibble(species = sp, rel_age = rel)
})

k_cols = grep("^k\\[", names(post), value = TRUE)
a50_cols = grep("^a50\\[", names(post), value = TRUE)

k_meds   = apply(post[, k_cols], 2, median)
k_sds   = apply(post[, k_cols], 2, sd)
a50_meds = apply(post[ , a50_cols], 2, median)
a50_sds   = apply(post[ , a50_cols], 2, sd)

k_df   = tibble(species = names(dfs), k = as.numeric(k_meds), k_sd = as.numeric(k_sds))
a50_df = tibble(species = names(dfs), a50 = as.numeric(a50_meds), a50_sd = as.numeric(a50_sds))

rel_age_df = rel_age_df %>%
  filter(species != "ARIFEL_PaP") %>% 
  left_join(k_df, by = "species") %>%
  left_join(a50_df, by = "species") %>%
  mutate(mu_sa = 1 / (1 + exp(-k * (rel_age - a50))),
         low_sa = 1 / (1 + exp(-(k - k_sd) * (rel_age - (a50 + a50_sd)))),
         high_sa = 1 / (1 + exp(-(k + k_sd) * (rel_age - (a50 - a50_sd)))))

# use species common names
sp_labels = c("Silver perch", "Blue crab", "White trout", "Spotted seatrout",
  "Gulf killifish", "Pinfish", "White shrimp")
rel_age_df$species = as.numeric(as.factor(rel_age_df$species)) # convert to numeric for factor levels
rel_age_df$species = factor(rel_age_df$species, 
                            levels = 1:7, 
                            labels = sp_labels)

ggplot(rel_age_df, aes(x = rel_age, y = mu_sa)) +
  geom_line(linewidth = 0.5) +
  #geom_ribbon(aes(ymin = low_sa, ymax = high_sa), alpha = 0.1) +
  labs(x = "Relative age (days)", y = "Selectivity", color = "Species") +
  custom_theme() +
  facet_wrap(~ species, scales = "free_x", nrow = 4)

# save the selectivity plot
fig_dir = here::here("res", "figures")
ggsave(
  filename = file.path(fig_dir, "selectivity_curves.pdf"),
  plot = last_plot(),
  width = 4.5, height = 5.5, dpi = 300
)

# 2. Observed vs. predicted catch compositions
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

K <- length(A)
n_draws <- 500  # Number of posterior draws to use for predictive intervals
draw_ids <- sample(nrow(post), n_draws)
pred_list <- vector("list", K)
k_names   <- paste0("k[",   1:K, "]")
a50_names <- paste0("a50[", 1:K, "]")
M_names   <- paste0("M[",   1:K, "]")

for (s in 1:K) {
  idxs <- start_idx[s]:end_idx[s]
  n_bin <- length(idxs)
  rel_age_s <- rel_age[idxs]
  obs_s <- counts_concat[idxs]
  pred_mat <- matrix(NA_real_, nrow = n_draws, ncol = n_bin)
  for (d in 1:n_draws) {
    k   <- as.numeric(post[draw_ids[d], k_names[s]])
    a50 <- as.numeric(post[draw_ids[d], a50_names[s]])
    M   <- as.numeric(post[draw_ids[d], M_names[s]])
    Sa <- 1 / (1 + exp(-k * (rel_age_s - a50)))
    Na <- cumprod(c(1, rep(exp(-M), n_bin - 1)))
    pred <- Sa * Na
    pred_vec <- as.numeric(pred / sum(pred))
    pred_mat[d, ] <- pred_vec # pred_vec should be length n_bin
  }
  pred_df <- data.frame(
    bin = seq_along(idxs),
    rel_age = rel_age_s,
    obs_prop = obs_s / sum(obs_s),
    pred_mean = colMeans(pred_mat, na.rm = TRUE),
    pred_lower = apply(pred_mat, 2, quantile, 0.025, na.rm = TRUE),
    pred_upper = apply(pred_mat, 2, quantile, 0.975, na.rm = TRUE),
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
catch_comp_summary_table <- data.frame(
  sp    = sp_labels,
  M_mean  = sapply(2:8, function(s) mean(post[[paste0("M[", s, "]")]])),
  M_sd    = sapply(2:8, function(s) sd(post[[paste0("M[", s, "]")]])),
  M_2.5   = sapply(2:8, function(s) quantile(post[[paste0("M[", s, "]")]], probs = 0.025)),
  M_97.5  = sapply(2:8, function(s) quantile(post[[paste0("M[", s, "]")]], probs = 0.975))
)
rownames(catch_comp_summary_table) <- NULL
catch_comp_summary_table[7,2:5] = catch_comp_summary_table[7,2:5]*scal 
print(catch_comp_summary_table)
