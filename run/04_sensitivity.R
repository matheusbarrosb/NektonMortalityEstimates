# 1. Catch curves --------------------------------------------------------------
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
require(ggpubr)

## Source functions ------------------------------------------------------------
function_dir = here::here("R")
source_files = list.files(function_dir, pattern = "\\.R$", full.names = TRUE)
for (file in source_files) source(file)

### Load data ------------------------------------------------------------------
data_dir  = here::here("data")
files     = list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
files = files[!grepl("ShrimpRegression.csv", files)]
data_list = lapply(files, read.csv) 

names(data_list) = tools::file_path_sans_ext(basename(files))

## Process data ----------------------------------------------------------------
# include dates
for (i in 1:length(data_list)) {
  data_list[[i]]$samp_date = as.Date(data_list[[i]]$samp_date, format = "%d.%m.%Y")
}

data_list$LITSET_PaP$Length = data_list$LITSET_PaP$Length * 4.781 + 5.847

# create length frequency objects
lfqs = list()
bin_sizes = c(2,1,2,2,2,2,3,1)
for (i in 1:length(data_list)) {
  lfqs[[i]] = lfqCreate(data     = data_list[[i]],
                        Lname    = "Length",
                        Dname    = "samp_date",
                        bin_size = bin_sizes[i])
}

# create catch composition dataframes
catches    = list()
sd_catches = list()
mids       = list()

for (i in 1:length(lfqs)) {
  # catches[[i]] = rowMeans(lfqs[[i]]$catch)
  catches[[i]] = rowSums(lfqs[[i]]$catch) + 1
  sd_catches[[i]] = apply(lfqs[[i]]$catch, 1, function(x) sd(x)/sqrt(length(x)))
  mids[[i]]    = lfqs[[i]]$midLengths 
}; names(catches) = names(data_list); names(mids) = names(data_list)

dfs = list()
for (i in 1:length(catches)) {
  dfs[[i]] = data.frame(catch = catches[[i]], sd_catches[[i]], mids  = mids[[i]])
  colnames(dfs[[i]])[2] = "sd_catches"
}; names(dfs) = names(data_list)

## Run - 20%--------------------------------------------------------------------
k_list    = c(0.25, 0.11/30, 0.513, 0.325/365, 0.61, 2.43/365, 1.35, 0.815) * 0.8
linf_list = c(410, 145, NA, 336.85, NA, 87.27, NA, NA) * 0.8
absolute  = ifelse(is.na(linf_list), TRUE, FALSE)
ex_points = c(30, 50, 86, 2, 108, 1, 27, 3)
bin_size  = c(2, 1, 2, 0.1, 2, 0.001, 3, 1)

cc_res = list()
N_runs = 1000
t0     = FALSE
plot   = FALSE
for (i in 1:length(data_list)) {
  
  if (absolute[i] == FALSE) {
    
    cc_res[[i]] = iLCCC(
      
      mids     = dfs[[i]]$mids,
      catch    = dfs[[i]]$catch,
      K        = k_list[i],
      Linf     = linf_list[i],
      t0       = t0,
      binsize  = bin_size[i],
      N_runs   = N_runs,
      plot     = plot,
      ex.points = ex_points[i]
    )
    
  } else {
    
    cc_res[[i]] = iLCCC_abgrowth(
      
      mids     = dfs[[i]]$mids,
      catch    = dfs[[i]]$catch,
      GR       = k_list[i],
      t0       = t0,
      binsize  = bin_size[i],
      N_runs   = N_runs,
      plot     = plot,
      ex.points = ex_points[i]
    )
    
  }
  
}; names(cc_res) = names(data_list)

## table with species and mortality column
mortality_table = data.frame(
  sp = names(cc_res),
  M = sapply(cc_res, function(x) x$Z)
)

catch_curve_summary_minus_20 = data.frame(
  sp     = titles,
  M_mean = sapply(cc_res, function(x) mean(x$Z)),
  M_sd   = sapply(cc_res, function(x) sd(x$`Z posterior distribution`, na.rm = TRUE)),
  M_2.5  = sapply(cc_res, function(x) quantile(x$`Z posterior distribution`, probs = 0.025, na.rm = TRUE)),
  M_97.5 = sapply(cc_res, function(x) quantile(x$`Z posterior distribution`, probs = 0.975, na.rm = TRUE))
);rownames(catch_curve_summary_minus_20) = NULL
print(catch_curve_summary_minus_20)

## Run - 10% --------------------------------------------------------------------
k_list    = c(0.25, 0.11/30, 0.513, 0.325/365, 0.61, 2.43/365, 1.35, 0.815) * 0.9
linf_list = c(410, 145, NA, 336.85, NA, 87.27, NA, NA) * 0.9
absolute  = ifelse(is.na(linf_list), TRUE, FALSE)
ex_points = c(30, 50, 86, 2, 108, 1, 27, 3)
bin_size  = c(2, 1, 2, 0.1, 2, 0.001, 3, 1)

cc_res = list()
N_runs = 1000
t0     = FALSE
plot   = FALSE
for (i in 1:length(data_list)) {
  
  if (absolute[i] == FALSE) {
    
    cc_res[[i]] = iLCCC(
      
      mids     = dfs[[i]]$mids,
      catch    = dfs[[i]]$catch,
      K        = k_list[i],
      Linf     = linf_list[i],
      t0       = t0,
      binsize  = bin_size[i],
      N_runs   = N_runs,
      plot     = plot,
      ex.points = ex_points[i]
    )
    
  } else {
    
    cc_res[[i]] = iLCCC_abgrowth(
      
      mids     = dfs[[i]]$mids,
      catch    = dfs[[i]]$catch,
      GR       = k_list[i],
      t0       = t0,
      binsize  = bin_size[i],
      N_runs   = N_runs,
      plot     = plot,
      ex.points = ex_points[i]
    )
    
  }
  
}; names(cc_res) = names(data_list)

## table with species and mortality column
mortality_table = data.frame(
  sp = names(cc_res),
  M = sapply(cc_res, function(x) x$Z)
)

catch_curve_summary_minus_10 = data.frame(
  sp     = titles,
  M_mean = sapply(cc_res, function(x) mean(x$Z)),
  M_sd   = sapply(cc_res, function(x) sd(x$`Z posterior distribution`, na.rm = TRUE)),
  M_2.5  = sapply(cc_res, function(x) quantile(x$`Z posterior distribution`, probs = 0.025, na.rm = TRUE)),
  M_97.5 = sapply(cc_res, function(x) quantile(x$`Z posterior distribution`, probs = 0.975, na.rm = TRUE))
);rownames(catch_curve_summary_minus_10) = NULL
print(catch_curve_summary_minus_10)

## Run + 10% --------------------------------------------------------------------
k_list    = c(0.25, 0.11/30, 0.513, 0.325/365, 0.61, 2.43/365, 1.35, 0.815) * 1.1
linf_list = c(410, 145, NA, 336.85, NA, 87.27, NA, NA) * 1.1
absolute  = ifelse(is.na(linf_list), TRUE, FALSE)
ex_points = c(30, 65, 86, 2, 108, 1, 27, 3)
bin_size  = c(2, 1, 2, 0.1, 2, 0.001, 3, 1)

cc_res = list()
N_runs = 1000
t0     = FALSE
plot   = FALSE
for (i in 1:length(data_list)) {
  
  if (absolute[i] == FALSE) {
    
    cc_res[[i]] = iLCCC(
      
      mids     = dfs[[i]]$mids,
      catch    = dfs[[i]]$catch,
      K        = k_list[i],
      Linf     = linf_list[i],
      t0       = t0,
      binsize  = bin_size[i],
      N_runs   = N_runs,
      plot     = plot,
      ex.points = ex_points[i]
    )
    
  } else {
    
    cc_res[[i]] = iLCCC_abgrowth(
      
      mids     = dfs[[i]]$mids,
      catch    = dfs[[i]]$catch,
      GR       = k_list[i],
      t0       = t0,
      binsize  = bin_size[i],
      N_runs   = N_runs,
      plot     = plot,
      ex.points = ex_points[i]
    )
    
  }
  
}; names(cc_res) = names(data_list)

## table with species and mortality column
mortality_table = data.frame(
  sp = names(cc_res),
  M = sapply(cc_res, function(x) x$Z)
)

catch_curve_summary_plus_10 = data.frame(
  sp     = titles,
  M_mean = sapply(cc_res, function(x) mean(x$Z)),
  M_sd   = sapply(cc_res, function(x) sd(x$`Z posterior distribution`, na.rm = TRUE)),
  M_2.5  = sapply(cc_res, function(x) quantile(x$`Z posterior distribution`, probs = 0.025, na.rm = TRUE)),
  M_97.5 = sapply(cc_res, function(x) quantile(x$`Z posterior distribution`, probs = 0.975, na.rm = TRUE))
);rownames(catch_curve_summary_plus_10) = NULL
print(catch_curve_summary_plus_10)

## Run + 20% --------------------------------------------------------------------
k_list    = c(0.25, 0.11/30, 0.513, 0.325/365, 0.61, 2.43/365, 1.35, 0.815) * 1.2
linf_list = c(410, 145, NA, 336.85, NA, 87.27, NA, NA) * 1.2
absolute  = ifelse(is.na(linf_list), TRUE, FALSE)
ex_points = c(30, 65, 86, 2, 108, 3, 27, 3)
bin_size  = c(2, 1, 2, 0.1, 2, 0.001, 3, 1)

cc_res = list()
N_runs = 1000
t0     = FALSE
plot   = FALSE
for (i in 1:length(data_list)) {
  
  if (absolute[i] == FALSE) {
    
    cc_res[[i]] = iLCCC(
      
      mids     = dfs[[i]]$mids,
      catch    = dfs[[i]]$catch,
      K        = k_list[i],
      Linf     = linf_list[i],
      t0       = t0,
      binsize  = bin_size[i],
      N_runs   = N_runs,
      plot     = plot,
      ex.points = ex_points[i]
    )
    
  } else {
    
    cc_res[[i]] = iLCCC_abgrowth(
      
      mids     = dfs[[i]]$mids,
      catch    = dfs[[i]]$catch,
      GR       = k_list[i],
      t0       = t0,
      binsize  = bin_size[i],
      N_runs   = N_runs,
      plot     = plot,
      ex.points = ex_points[i]
    )
    
  }
  
}; names(cc_res) = names(data_list)

## table with species and mortality column
mortality_table = data.frame(
  sp = names(cc_res),
  M = sapply(cc_res, function(x) x$Z)
)

catch_curve_summary_plus_20 = data.frame(
  sp     = titles,
  M_mean = sapply(cc_res, function(x) mean(x$Z)),
  M_sd   = sapply(cc_res, function(x) sd(x$`Z posterior distribution`, na.rm = TRUE)),
  M_2.5  = sapply(cc_res, function(x) quantile(x$`Z posterior distribution`, probs = 0.025, na.rm = TRUE)),
  M_97.5 = sapply(cc_res, function(x) quantile(x$`Z posterior distribution`, probs = 0.975, na.rm = TRUE))
);rownames(catch_curve_summary_plus_20) = NULL
print(catch_curve_summary_plus_20)

# Combine summaries into one table and plot
catch_curve_summary = rbind(
  catch_curve_summary_table    %>% mutate(sensitivity = "base"),
  catch_curve_summary_minus_20 %>% mutate(sensitivity = "- 20%"),
  catch_curve_summary_minus_10 %>% mutate(sensitivity = "- 10%"),
  catch_curve_summary_plus_10  %>% mutate(sensitivity = "+ 10%"),
  catch_curve_summary_plus_20  %>% mutate(sensitivity = "+ 20%")
)

# 2. Catch comps ---------------------------------------------------------------
## - 10 % ----------------------------------------------------------------------
data_dir  = here::here("data")
files     = list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
files = files[!grepl("ShrimpRegression.csv", files)]
data_list = lapply(files, read.csv) 
names(data_list) = tools::file_path_sans_ext(basename(files))

k_list    = c(0.25/365, 0.11/30, 0.513, 0.325/365, 0.61, 2.43/365, 1.35, 0.815) * 0.9
linf_list = c(410, 145, NA, 336.85, NA, 87.27, NA, NA) * 0.9
absolute  = ifelse(is.na(linf_list), TRUE, FALSE)

for (i in 1:length(data_list)) {
  data_list[[i]]$samp_date = as.Date(data_list[[i]]$samp_date, format = "%d.%m.%Y")
}
data_list$LITSET_PaP$Length = data_list$LITSET_PaP$Length * 4.781 + 5.847

# cutting underrepresented lengths
data_list$CALSAP_PaP = data_list$CALSAP_PaP[data_list$CALSAP_PaP$Length <= 40, ]
data_list$CYNNEB_PaP = data_list$CYNNEB_PaP[data_list$CYNNEB_PaP$Length <= 100, ]
data_list$BAICHR_PaP = data_list$BAICHR_PaP[data_list$BAICHR_PaP$Length <= 70, ]

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
data_list$M_prior_sd[6] = 0.1
data_list$a50_prior_mean = c(400, 80, 15, 120, 60, 80, 40, 60)
data_list$a50_prior_sd   = c(100, 10,   3,   5,  2, 5, 40,  2)

# fit model --------------------------------------------------------------------
stan_model_1 = cmdstanr::cmdstan_model("stan/cc.stan")
n_chains = 2
fit = stan_model_1$sample(
  data          = data_list,
  iter_sampling = 5000,
  iter_warmup   = 1000,
  chains        = n_chains,
  seed          = 2027
)

post = fit$draws(variables = c("k", "a50", "M"), format = "df")

sp_labels = c("Silver perch", "Blue crab", "White trout", "Spotted seatrout",
              "Gulf killifish", "Pinfish", "White shrimp")

catch_comp_summary_minus_10 = data.frame(
  sp    = sp_labels,
  M_mean  = sapply(2:8, function(s) mean(post[[paste0("M[", s, "]")]])),
  M_sd    = sapply(2:8, function(s) sd(post[[paste0("M[", s, "]")]])),
  M_2.5   = sapply(2:8, function(s) quantile(post[[paste0("M[", s, "]")]], probs = 0.025)),
  M_97.5  = sapply(2:8, function(s) quantile(post[[paste0("M[", s, "]")]], probs = 0.975))
)
rownames(catch_comp_summary_minus_10) <- NULL
catch_comp_summary_minus_10[7,2:5] = catch_comp_summary_minus_10[7,2:5]*scal

## - 20 % ----------------------------------------------------------------------
data_dir  = here::here("data")
files     = list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
files = files[!grepl("ShrimpRegression.csv", files)]
data_list = lapply(files, read.csv) 
names(data_list) = tools::file_path_sans_ext(basename(files))

k_list    = c(0.25/365, 0.11/30, 0.513, 0.325/365, 0.61, 2.43/365, 1.35, 0.815) * 0.8
linf_list = c(410, 145, NA, 336.85, NA, 87.27, NA, NA) * 0.8
absolute  = ifelse(is.na(linf_list), TRUE, FALSE)

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
data_list$M_prior_sd[6] = 0.1
data_list$a50_prior_mean = c(400, 80, 15, 120, 60, 80, 40, 60)
data_list$a50_prior_sd   = c(100, 10,   3,   5,  2, 5, 40,  100)

# fit model --------------------------------------------------------------------
stan_model_1 = cmdstanr::cmdstan_model("stan/cc.stan")
n_chains = 2
fit = stan_model_1$sample(
  data          = data_list,
  iter_sampling = 5000,
  iter_warmup   = 1000,
  chains        = n_chains,
  seed          = 2027
)

post = fit$draws(variables = c("k", "a50", "M"), format = "df")

sp_labels = c("Silver perch", "Blue crab", "White trout", "Spotted seatrout",
              "Gulf killifish", "Pinfish", "White shrimp")

catch_comp_summary_minus_20 = data.frame(
  sp    = sp_labels,
  M_mean  = sapply(2:8, function(s) mean(post[[paste0("M[", s, "]")]])),
  M_sd    = sapply(2:8, function(s) sd(post[[paste0("M[", s, "]")]])),
  M_2.5   = sapply(2:8, function(s) quantile(post[[paste0("M[", s, "]")]], probs = 0.025)),
  M_97.5  = sapply(2:8, function(s) quantile(post[[paste0("M[", s, "]")]], probs = 0.975))
)
rownames(catch_comp_summary_minus_20) <- NULL
catch_comp_summary_minus_20[7,2:5] = catch_comp_summary_minus_20[7,2:5]*scal

## + 10 % ----------------------------------------------------------------------
data_dir  = here::here("data")
files     = list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
files = files[!grepl("ShrimpRegression.csv", files)]
data_list = lapply(files, read.csv) 
names(data_list) = tools::file_path_sans_ext(basename(files))

k_list    = c(0.25/365, 0.11/30, 0.513, 0.325/365, 0.61, 2.43/365, 1.35, 0.815) * 1.1
linf_list = c(410, 145, NA, 336.85, NA, 87.27, NA, NA) * 1.1
absolute  = ifelse(is.na(linf_list), TRUE, FALSE)

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
data_list$M_prior_sd[6] = 0.1
data_list$a50_prior_mean = c(400, 80, 15, 120, 60, 80, 40, 60)
data_list$a50_prior_sd   = c(100, 10,   3,   5,  2, 5, 40,  100)

# fit model --------------------------------------------------------------------
stan_model_1 = cmdstanr::cmdstan_model("stan/cc.stan")
n_chains = 2
fit = stan_model_1$sample(
  data          = data_list,
  iter_sampling = 5000,
  iter_warmup   = 1000,
  chains        = n_chains,
  seed          = 2027
)

post = fit$draws(variables = c("k", "a50", "M"), format = "df")

sp_labels = c("Silver perch", "Blue crab", "White trout", "Spotted seatrout",
              "Gulf killifish", "Pinfish", "White shrimp")

catch_comp_summary_plus_10 = data.frame(
  sp    = sp_labels,
  M_mean  = sapply(2:8, function(s) mean(post[[paste0("M[", s, "]")]])),
  M_sd    = sapply(2:8, function(s) sd(post[[paste0("M[", s, "]")]])),
  M_2.5   = sapply(2:8, function(s) quantile(post[[paste0("M[", s, "]")]], probs = 0.025)),
  M_97.5  = sapply(2:8, function(s) quantile(post[[paste0("M[", s, "]")]], probs = 0.975))
)
rownames(catch_comp_summary_plus_10) <- NULL
catch_comp_summary_plus_10[7,2:5] = catch_comp_summary_plus_10[7,2:5]*scal


## + 20 % ----------------------------------------------------------------------
data_dir  = here::here("data")
files     = list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
files = files[!grepl("ShrimpRegression.csv", files)]
data_list = lapply(files, read.csv) 
names(data_list) = tools::file_path_sans_ext(basename(files))

k_list    = c(0.25/365, 0.11/30, 0.513, 0.325/365, 0.61, 2.43/365, 1.35, 0.815) * 1.2
linf_list = c(410, 145, NA, 336.85, NA, 87.27, NA, NA) * 1.2
absolute  = ifelse(is.na(linf_list), TRUE, FALSE)

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
data_list$M_prior_sd[6] = 0.1
data_list$a50_prior_mean = c(400, 80, 15, 120, 60, 80, 40, 60)
data_list$a50_prior_sd   = c(100, 10,   3,   5,  2, 5, 40,  100)

# fit model --------------------------------------------------------------------
stan_model_1 = cmdstanr::cmdstan_model("stan/cc.stan")
n_chains = 2
fit = stan_model_1$sample(
  data          = data_list,
  iter_sampling = 5000,
  iter_warmup   = 1000,
  chains        = n_chains,
  seed          = 2027
)

post = fit$draws(variables = c("k", "a50", "M"), format = "df")

sp_labels = c("Silver perch", "Blue crab", "White trout", "Spotted seatrout",
              "Gulf killifish", "Pinfish", "White shrimp")

catch_comp_summary_plus_20 = data.frame(
  sp    = sp_labels,
  M_mean  = sapply(2:8, function(s) mean(post[[paste0("M[", s, "]")]])),
  M_sd    = sapply(2:8, function(s) sd(post[[paste0("M[", s, "]")]])),
  M_2.5   = sapply(2:8, function(s) quantile(post[[paste0("M[", s, "]")]], probs = 0.025)),
  M_97.5  = sapply(2:8, function(s) quantile(post[[paste0("M[", s, "]")]], probs = 0.975))
)
rownames(catch_comp_summary_plus_20) <- NULL
catch_comp_summary_plus_20[7,2:5] = catch_comp_summary_plus_20[7,2:5]*scal

# Plotting ---------------------------------------------------------------------

# Combine summaries into one table and plot
catch_comp_summary = rbind(
  catch_comp_summary_table    %>% mutate(sensitivity = "base"),
  catch_comp_summary_minus_20 %>% mutate(sensitivity = "- 20%"),
  catch_comp_summary_minus_10 %>% mutate(sensitivity = "- 10%"),
  catch_comp_summary_plus_10  %>% mutate(sensitivity = "+ 10%"),
  catch_comp_summary_plus_20  %>% mutate(sensitivity = "+ 20%")
)
catch_comp_summary %>%
  filter(sp != "Hardhead catfish") %>%
  pivot_wider(
    names_from = sensitivity,
    values_from = c(M_mean, M_sd, M_2.5, M_97.5)
  ) %>%
  mutate(
    M_mean_diff_20 = M_mean_base - `M_mean_- 20%`,
    M_mean_diff_10 = M_mean_base - `M_mean_- 10%`,
    M_mean_diff_plus_10 = M_mean_base - `M_mean_+ 10%`,
    M_mean_diff_plus_20 = M_mean_base - `M_mean_+ 20%`
  ) %>%
  select(sp, M_mean_base, M_mean_diff_20, M_mean_diff_10, M_mean_diff_plus_10, M_mean_diff_plus_20) %>%
  pivot_longer(
    cols = starts_with("M_mean_diff"),
    names_to = "sensitivity",
    values_to = "M_mean_diff"
  ) %>%
  
  ggplot(aes(x = sp, y = 100*M_mean_diff, fill = sensitivity)) +
    geom_bar(stat = "identity", position = position_dodge(0.9), color = "black") +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    coord_flip() +
    custom_theme() +
    xlab("") +
    ylab("% difference")

# make a combined plot
sens_df = rbind(
  catch_curve_summary %>% mutate(type = "Catch curve"),
  catch_comp_summary %>% mutate(type = "IBCaL")
)

sens_df %>%
  filter(sp != "Hardhead catfish") %>%
  pivot_wider(
    names_from = sensitivity,
    values_from = c(M_mean, M_sd, M_2.5, M_97.5)
  ) %>%
  mutate(
    M_mean_diff_20 = M_mean_base - `M_mean_- 20%`,
    M_mean_diff_10 = M_mean_base - `M_mean_- 10%`,
    M_mean_diff_plus_10 = M_mean_base - `M_mean_+ 10%`,
    M_mean_diff_plus_20 = M_mean_base - `M_mean_+ 20%`
  ) %>%
  select(sp, type, M_mean_base, M_mean_diff_20, M_mean_diff_10, M_mean_diff_plus_10, M_mean_diff_plus_20) %>%
  pivot_longer(
    cols = starts_with("M_mean_diff"),
    names_to = "sensitivity",
    values_to = "M_mean_diff"
  ) %>%
  # rename 'sensitivity' column for better labels
  mutate(sensitivity = recode(sensitivity,
                              M_mean_diff_20 = "- 20%",
                              M_mean_diff_10 = "- 10%",
                              M_mean_diff_plus_10 = "+ 10%",
                              M_mean_diff_plus_20 = "+ 20%")) %>%
  
  ggplot(aes(x = sp, y = 100*M_mean_diff, fill = sensitivity)) +
    geom_bar(stat = "identity", position = position_dodge(0.9), color = "black") +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    coord_flip() +
    custom_theme() +
    xlab("") +
    ylab("% difference") +
    facet_wrap(~type, nrow = 2) +
    theme(legend.title = element_blank(),
          legend.position = "right") +
    scale_fill_manual(values = c("red", "black", "honeydew4", "steelblue1"))

# save figure
res_dir = here::here("res/figures")
ggsave(
  last_plot(),
  filename = file.path(res_dir, "sensitivity.pdf"),
  width = 8, height = 6, dpi = 300
)

