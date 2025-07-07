for (i in 1:length(data_list)) {
  data_list[[i]]$samp_date = as.Date(data_list[[i]]$samp_date, format = "%d.%m.%Y")
}

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

# plot size distribution histograms
sp_labels = c("Silver perch", "Blue crab", "White trout", "Spotted seatrout",
              "Gulf killifish", "Pinfish", "White shrimp")

# join dfs and add species names column
dfs_joined = do.call(rbind, lapply(names(dfs), function(x) {
  df = dfs[[x]]
  df$sp = x
  return(df)
})); dfs_joined = dfs_joined %>% filter(catch > 0) %>% filter(sp != "ARIFEL_PaP") 
# replace sp labels with sp_labels
dfs_joined$sp = factor(dfs_joined$sp, levels = names(data_list[-1]), labels = sp_labels)

# transform to long format for histograms 
dfs_joined %>%
  select(sp, mids, catch) %>%
  uncount(catch) %>%
  rename(length = mids) %>%
  
  ggplot(aes(x = length)) +
  geom_histogram(bins = 20, position = "identity", color = "black", fill = "red", alpha = 0.7) +
  facet_wrap(~ sp, scales = "free", nrow = 4) +
  custom_theme() +
  xlab("Length (mm)") +
  ylab("Count")
  
fig_dir = here::here("res", "figures")
ggsave(
  filename = file.path(fig_dir, "lfqs.pdf"),
  width = 4, height = 5, dpi = 300
)

# plot M posterior distribution across methods
M_df = rbind(
  catch_curve_summary_table %>% mutate(method = "Catch curves"),
    catch_comp_summary_table %>% mutate(method = "IBCL")
)

M_df[4,4] = 0
M_df %>%
  # filter catchfish
  filter(sp %in% sp_labels) %>%
  ggplot(aes(x = sp, y = M_mean)) +
  geom_point(aes(shape = method, color = method), position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = M_2.5, ymax = M_97.5, color = method),
                width = 0, position = position_dodge(width = 0.5)) +
  coord_flip() +
  custom_theme() +
  xlab("") +
  ylab(expression(M~"(day"^{-1}*")")) +
  scale_color_manual(values = c("IBCL" = "red", "Catch curves" = "black")) +
  theme(legend.title = element_blank(),
        legend.position = "top")

fig_dir = here::here("res", "figures")
ggsave(
  filename = file.path(fig_dir, "M_posterior_distribution.pdf"),
  width = 4, height = 4, dpi = 300
)

library(knitr)
library(kableExtra)
# Combine the two summary tables
combined_summary_table <- merge(
  catch_curve_summary_table %>% mutate(method = "Catch curves"),
  catch_comp_summary_table %>% mutate(method = "IBCL"),
  by = "sp", suffixes = c("_iLCCC", "_CatchComp")
); catch_curve_summary_table[7,4] = 0
# Select relevant columns and rename them
pretty_table <- combined_summary_table %>%
  select(sp, 
         M_mean_iLCCC, M_sd_iLCCC, M_2.5_iLCCC, M_97.5_iLCCC,
         M_mean_CatchComp, M_sd_CatchComp, M_2.5_CatchComp, M_97.5_CatchComp) %>%
  rename(
    Species = sp,
    `i-LCCC M mean` = M_mean_iLCCC,
    `i-LCCC M sd` = M_sd_iLCCC,
    `i-LCCC M 2.5%` = M_2.5_iLCCC,
    `i-LCCC M 97.5%` = M_97.5_iLCCC,
    `Catch Comp M mean` = M_mean_CatchComp,
    `Catch Comp M sd` = M_sd_CatchComp,
    `Catch Comp M 2.5%` = M_2.5_CatchComp,
    `Catch Comp M 97.5%` = M_97.5_CatchComp
  )

# Print the and export the table with kableExtra for better formatting
table = kable(pretty_table, format = "pipe", escape = TRUE, digits = 3) %>%
  kable_styling(full_width = TRUE, position = "left") %>%
  column_spec(1, bold = TRUE) %>%
  add_header_above(c(" " = 1, "i-LCCC" = 4, "Catch Composition" = 4)) %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(nrow(pretty_table), bold = FALSE)


kableExtra::save_kable(
  table,
  file = file.path(res_dir, "M_summary_table.pdf"),
  self_contained = TRUE
)






