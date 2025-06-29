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
data_list = lapply(files, read.csv) 

names(data_list) = tools::file_path_sans_ext(basename(files))

## Process data ----------------------------------------------------------------
# include dates
for (i in 1:length(data_list)) {
  data_list[[i]]$samp_date = as.Date(data_list[[i]]$samp_date, format = "%d.%m.%Y")
}

# create length frequency objects
lfqs = list()
bin_sizes = c(2,2,2,2,2,2,3,1)
for (i in 1:length(data_list)) {
  lfqs[[i]] = lfqCreate(data     = data_list[[i]],
                        Lname    = "Length",
                        Dname    = "samp_date",
                        bin_size = bin_sizes[i])
}

# create catch composition dataframes
catches = list()
sd_catches = list()
mids    = list()

for (i in 1:length(lfqs)) {
 # catches[[i]] = rowMeans(lfqs[[i]]$catch)
 catches[[i]] = rowSums(lfqs[[i]]$catch)
  sd_catches[[i]] = apply(lfqs[[i]]$catch, 1, function(x) sd(x)/sqrt(length(x)))
  mids[[i]]    = lfqs[[i]]$midLengths 
}; names(catches) = names(data_list); names(mids) = names(data_list)

dfs = list()
for (i in 1:length(catches)) {
  dfs[[i]] = data.frame(catch = catches[[i]], sd_catches[[i]], mids  = mids[[i]])
  colnames(dfs[[i]])[2] = "sd_catches"
}; names(dfs) = names(data_list)

dfs$LITSET_PaP$mids = dfs$LITSET_PaP$mids * 10

## Run catch curves ------------------------------------------------------------
k_list    = c(0.25, 0.11/30, 0.513, 0.325/365, 0.61, 2.43/365, 1.35, 0.815)
linf_list = c(410, 145, NA, 336.85, NA, 87.27, NA, NA)
absolute  = ifelse(is.na(linf_list), TRUE, FALSE)
ex_points = c(50, 10, 75, 1, 105, 3, 0, 3)
bin_size  = c(0.1, 2, 2, 0.1, 2, 0.1, 2, 0.1)

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

## Plotting --------------------------------------------------------------------
titles = c("Hardhead catfish", "Silver perch", "Blue crab", "White trout", "Spotted seatrout",
           "Gulf killifish", "Pinfish", "White shrimp")

plots = list(
  plot_catch_curve(dfs$BAICHR_PaP, cc_res$BAICHR_PaP, pars = c(0.11/30, 145), absolute = FALSE, bin_size = 2) +
    xlim(50,210) + ylim(0,3) + xlab("") + ylab("") + ggtitle(titles[2]),
  plot_catch_curve(dfs$CALSAP_PaP, cc_res$CALSAP_PaP, pars = 0.513, absolute = TRUE, bin_size = 2) +
    xlim(0,80) + ylim(1,4) + xlab("") + ylab("") + ggtitle(titles[3]),
  plot_catch_curve(dfs$CYNARE_PaP, cc_res$CYNARE_PaP, pars = c(k_list[4], linf_list[4]), absolute = FALSE, bin_size = 0.1) +
    xlim(100, 270) + ylim(0,4) + xlab("") + ylab("") + ggtitle(titles[4]),
  plot_catch_curve(dfs$CYNNEB_PaP, cc_res$CYNNEB_PaP, pars = k_list[5], absolute = TRUE, bin_size = 2) +
    xlim(60,130) + ylim(1.5,3.5) + xlab("") + ylab("") + ggtitle(titles[4]),
  plot_catch_curve(dfs$FUNGRA_PaP, cc_res$FUNGRA_PaP, pars = c(k_list[6], linf_list[6]), absolute = FALSE, bin_size = 0.1) +
    xlim(60, 210) + ylim(1.5, 4) + xlab("") + ylab("") + ggtitle(titles[6]),
  plot_catch_curve(dfs$LAGRHO_PaP, cc_res$LAGRHO_PaP, pars = k_list[7], absolute = TRUE, bin_size = 2) +
    xlim(10, 130) + ylim(1, 3.5) + xlab("") + ylab("") + ggtitle(titles[7]),
  plot_catch_curve(dfs$LITSET_PaP, cc_res$LITSET_PaP, pars = k_list[8], absolute = TRUE) +
    xlim(40, 300) + xlab("") + ylab("") + ggtitle(titles[8])
)

plots = lapply(
  plots,
  function(p) p + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 9),
      plot.margin = margin(1, 0.5, 1, 0.5)  # left, top, right, bottom (in pts)
    )
)

# Arrange the plots with ggarrange
combined_plot = ggarrange(
  plotlist = plots,
  ncol = 2, nrow = 4, align = "hv",
  common.legend = TRUE, legend = "bottom")

annotate_figure(combined_plot,
         bottom = text_grob("Relative age (days)", size = 12, vjust = -1),
         left = text_grob("log Catch", rot = 90, size = 12, vjust = 1))

fig_dir = here::here("res", "figures")
ggsave(
  filename = file.path(fig_dir, "catch_curves.pdf"),
  plot = combined_plot,
  width = 4.5, height = 5.5, units = "in", dpi = 300
)

# summary of M posterior distributions (means, sds, quantiles) for each species
catch_curve_summary_table = data.frame(
  sp = titles,
  M_mean = sapply(cc_res, function(x) mean(x$Z)),
  M_sd = sapply(cc_res, function(x) sd(x$`Z posterior distribution`, na.rm = TRUE)),
  M_2.5 = sapply(cc_res, function(x) quantile(x$`Z posterior distribution`, probs = 0.025)),
  M_97.5 = sapply(cc_res, function(x) quantile(x$`Z posterior distribution`, probs = 0.975))
);rownames(summary_table) = NULL




