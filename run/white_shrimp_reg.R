data_dir  = here::here("data")
shrimp_reg = read.csv(file.path(data_dir, "ShrimpRegression.csv"))

#exclude columns > index 3
shrimp_reg = shrimp_reg[, 1:3]

names(shrimp_reg) = c("Haul", "CL", "TL")
shrimp_reg = na.exclude(shrimp_reg)


# fit regression model
model = lm(TL ~ CL, data = shrimp_reg)

# Call:
#   lm(formula = TL ~ CL, data = shrimp_reg)
# 
# Coefficients:
#   (Intercept)       (Slope)  
#    5.847            4.781 

# plot relationship
library(ggplot2)
ggplot(shrimp_reg, aes(x = CL, y = TL)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Total Length vs Carapace Length",
       x = "Carapace Length (mm)",
       y = "Total Length (mm)") +
  theme_minimal() +
  xlim(0,30) + ylim(0,110)

