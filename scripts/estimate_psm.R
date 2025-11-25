
# Load necessary libraries
library(eurostat)
library(tidyverse)

# Load data from eurostat ----

#Number of pigs (tag00018) 	
pigs_eurostat <- get_eurostat(id = "tag00018", time_format = "date") %>% 
  label_eurostat()
  
# Production of pigmeat in slaughterhouses (tag00042)
pigmeat_eurostat <- get_eurostat(id = "tag00042", time_format = "date") %>% 
  label_eurostat()
# Slaughter of pigs (apro_mt_pann)
slaughter_eurostat <- get_eurostat(id = 'apro_mt_pann', time_format = "date") %>% 
  label_eurostat()

## Estimate PSM ----
# Filter data for 2024
pigs_slaughter_2024 <- slaughter_eurostat %>% 
  filter(unit == 'Thousand heads (animals)') %>%
  filter(meat == 'Pigmeat') %>% 
  filter(TIME_PERIOD == '2024-01-01') 
# Filter pig population data for 2024
pigs_eurostat_2024 <- pigs_eurostat %>% 
  filter(unit == 'Thousand heads (animals)') %>%
  filter(TIME_PERIOD == '2024-01-01')
# Calculate monthly probability of slaughter
psm_table <- pigs_slaughter_2024 %>% 
  select(geo, values) %>%
  rename(slaughter = values) %>% 
  mutate(slaughter_month = slaughter / 12) %>%
  left_join(pigs_eurostat_2024 %>%
              select(geo, values) %>%
              rename(total_pigs = values),
            by = 'geo') %>%
  mutate(prob_slaughter = slaughter_month / total_pigs) 

# Calculate mean and standard deviation of slaughter probability
psm_estimated <- psm_table %>%
  summarise(mean_prob = mean(prob_slaughter, na.rm = TRUE),
            sd_prob = sd(prob_slaughter, na.rm = TRUE),
            min = min(prob_slaughter, na.rm = TRUE),
            max = max(prob_slaughter, na.rm = TRUE))

# checking distribution and fitting normal truncated distribution
hist(psm_table$prob_slaughter,
     breaks = 30, freq = FALSE, col = "lightgray",
     main = "Prob_slaughter",
     xlab = "prob_slaughter",
     xlim = c(0,0.25))

# plot normal curve over histogram
curve_vals <- dnorm(x_vals, mean = psm_estimated$mean_prob, sd = psm_estimated$sd_prob)
lines(x_vals, curve_vals, col = "blue", lwd = 2)

# Function to generate truncated normal distribution
truncnorm <- function(n, mean, sd, min, max) {
  q <- pnorm(c(min, max), mean = mean, sd = sd)
  u <- runif(n, q[1], q[2])
  return(qnorm(u, mean = mean, sd = sd))
}

# Generate random values from truncated normal distribution
random_values <- truncnorm(n = 1000, mean = psm_estimated$mean_prob , 
                           sd = psm_estimated$sd_prob,
                           min = psm_estimated$min , max = psm_estimated$max)
# Plot histogram of random values
hist(random_values, breaks = 30, freq = FALSE, col = rgb(0,0,1,0.5),
     add = TRUE)

# Save the psm_estimated table for future use in .csv format
write.csv(psm_estimated, "data/processed/psm_estimated.csv", row.names = FALSE)

## Estimate how many tons each slaughtered pig produces ----
# Filter data for 2024
pigmeat_2024 <- pigmeat_eurostat %>%
  filter(unit == 'Thousand tonnes') %>%
  filter(meat == 'Pigmeat') %>% 
  filter(TIME_PERIOD == '2024-01-01')
# Filter slaughter data for 2024
slaughter_2024 <- slaughter_eurostat %>%
  filter(unit == 'Thousand heads (animals)') %>%
  filter(meat == 'Pigmeat') %>% 
  filter(TIME_PERIOD == '2024-01-01')
# Join datasets and calculate tons per pig
tons_per_pig_table <- pigmeat_2024 %>%
  select(geo, values) %>%
  rename(pigmeat_tons = values) %>%
  left_join(slaughter_2024 %>%
              select(geo, values) %>%
              rename(slaughter = values),
            by = 'geo') %>%
  mutate(tons_per_pig = pigmeat_tons / slaughter)
# Calculate mean tons per pig
tons_per_pig_estimated <- tons_per_pig_table %>%
  summarise(mean_tons_per_pig = mean(tons_per_pig, na.rm = TRUE),
            sd_tons_per_pig = sd(tons_per_pig, na.rm = TRUE),
            min = min(tons_per_pig, na.rm = TRUE),
            max = max(tons_per_pig, na.rm = TRUE))
# Save the tons_per_pig_estimated table for future use in .csv format
write.csv(tons_per_pig_estimated, "data/processed/tons_per_pig_estimated.csv", row.names = FALSE)

# checking distribution and fitting normal truncated distribution
hist(tons_per_pig_table$tons_per_pig,
     breaks = 30, freq = FALSE, col = "lightgray",
     main = "Tons per pig",
     xlab = "tons per pig")
# plot normal curve over histogram
x_vals <- seq(0, 0.2, length.out = 100)
curve_vals <- dnorm(x_vals, mean = tons_per_pig_estimated$mean_tons_per_pig, sd = tons_per_pig_estimated$sd_tons_per_pig)
lines(x_vals, curve_vals, col = "blue", lwd = 2)
# Generate random values from truncated normal distribution
random_values_tons <- truncnorm(n = 1000, mean = tons_per_pig_estimated$mean_tons_per_pig , 
                               sd = tons_per_pig_estimated$sd_tons_per_pig,
                               min = tons_per_pig_estimated$min , max = tons_per_pig_estimated$max)
# Plot histogram of random values
hist(random_values_tons, breaks = 30, freq = FALSE, col = rgb(0,0,1,0.5),
     add = TRUE)




