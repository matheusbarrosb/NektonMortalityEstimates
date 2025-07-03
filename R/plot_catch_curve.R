plot_catch_curve = function(raw_data, cc_res, pars, absolute, bin_size) {
  
  if (absolute) {
  
    raw_data$ages = raw_data$mids / pars[1] # convert to relative ages
    reg = cc_res$`Regression input`
    
    # class.min <- mids - (bin_size/2)
    # class.max <- mids + (bin_size/2)
    # dt = (class.max - class.min)/pars[1] # amount of time from one size class to another
    
    p = raw_data %>%
      ggplot(aes(x = ages, y = log(catch + 1))) +
      geom_point() +
      geom_smooth(data = reg, aes(x = xvar, y = yvar), method = "lm", alpha = 0.5, color = "red") +
      custom_theme() +
      xlim(0, max(raw_data$ages) * 1.1) +
      labs(x = "Relative age (days)", y = "log Catch")
    
  } else {
    
    min = raw_data$mids - (bin_size / 2)
    max = raw_data$mids + (bin_size / 2)
    dt = (-(1/pars[1]) * log(1 - max/pars[2])) - (-(1/pars[1]) * log(1 - min/pars[2]))
    raw_data$ages = -(1/pars[1]) * log(1 - (raw_data$mids / pars[2]))# convert to relative ages
    reg = cc_res$`Regression input`
    
  p = raw_data %>%
      ggplot(aes(x = ages, y = log((catch+1)/dt))) +
      geom_point() +
      geom_smooth(data = reg, aes(x = xvar, y = yvar), method = "lm", alpha = 0.5, color = "red") +
      custom_theme() +
      xlim(0, max(raw_data$ages) * 1.1) +
      labs(x = "Relative age (days)", y = "log Catch")
    
  }
  
  return(p)
  
}

# 
# df = dfs$LITSET_PaP
# df$age = df$mids / 0.815
# 
# reg = cc_res$LITSET_PaP$`Regression input`
# 
# df %>%
#   ggplot(aes(x = age, y = log(catch+1))) +
#   geom_point() +
#   geom_smooth(data = reg, aes(x = xvar, y = yvar), method = "lm") +
#   custom_theme() +
#   xlim()
# 
# pars = 0.815
# plot_catch_curve(dfs$LITSET_PaP, cc_res$LITSET_PaP, pars, absolute = TRUE)
# 
# 




