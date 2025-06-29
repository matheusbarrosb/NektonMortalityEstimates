# plot M posterior distribution across methods
M_df = rbind(
  catch_curve_summary_table %>% mutate(method = "iLCCC"),
  catch_comp_summary_table %>% mutate(method = "Catch composition")
); M_df[nrow(M_df),2:5] = M_df[nrow(M_df),2:5]/10 

M_df %>%
  # filter catchfish
  filter(sp %in% sp_labels) %>%
  ggplot(aes(x = sp, y = M_mean)) +
  geom_point(aes(shape = method, color = method), position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = M_2.5, ymax = M_97.5, color = method),
                width = 0, position = position_dodge(width = 0.3)) +
  coord_flip() +
  custom_theme() +
  xlab("") +
  ylab(expression(M~"(day"^{-1}*")")) +
  scale_color_manual(values = c("Catch composition" = "red", "iLCCC" = "black")) +
  theme(legend.title = element_blank(),
        legend.position = "top")

fig_dir = here::here("res", "figures")
ggsave(
  filename = file.path(fig_dir, "M_posterior_distribution.pdf"),
  width = 4, height = 4, dpi = 300
)

# Unified pretty table with means, sds, and quantiles for M for each species and method. No repetition of species names.
library(knitr)
library(kableExtra)
# Combine the two summary tables
combined_summary_table <- merge(
  catch_curve_summary_table %>% mutate(method = "iLCCC"),
  catch_comp_summary_table %>% mutate(method = "Catch composition"),
  by = "sp", suffixes = c("_iLCCC", "_CatchComp")
)
# Select relevant columns and rename them
pretty_table <- combined_summary_table %>%
  select(sp, 
         M_mean_iLCCC, M_sd_iLCCC, M_2.5_iLCCC, M_97.5_iLCCC,
         M_mean_CatchComp, M_sd_CatchComp, M_2.5_CatchComp, M_97.5_CatchComp) %>%
  rename(
    Species = sp,
    `Mean ` = M_mean_iLCCC,
    `SD ` = M_sd_iLCCC,
    `2.5% ` = M_2.5_iLCCC,
    `97.5% ` = M_97.5_iLCCC,
    `Mean` = M_mean_CatchComp,
    `SD` = M_sd_CatchComp,
    `2.5%` = M_2.5_CatchComp,
    `97.5%` = M_97.5_CatchComp
  )
# Print the table with kableExtra for better formatting
kable(pretty_table, format = "pipe", escape = TRUE, digits = 3) %>%
  kable_styling(full_width = TRUE, position = "left") %>%
  column_spec(1, bold = TRUE) %>%
  add_header_above(c(" " = 1, "i-LCCC" = 4, "Catch Composition" = 4)) %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(nrow(pretty_table), bold = FALSE)

# export as .csv
res_dir = here::here("res", "data")
write.csv(
  pretty_table,
  file = file.path(res_dir, "M_summary_table.csv"),
  row.names = FALSE
)


