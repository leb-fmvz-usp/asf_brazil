
#load packages
library(tidyverse)

#loading data
table_meat_production <- read.csv("data/raw/nm_meat_pig_production_exporting_country_faostat.csv")

# remove M49 == 159, which is the aggregate of China mainland and China Taiwan
table_meat_production <- table_meat_production %>% 
  rename(M49 = Area.Code..M49.) %>% 
  filter(M49 != 159)

# filtering tonnes of product
table_meat_production <- table_meat_production %>% 
  filter(Unit == "t")

# sum of meat production by country and year
production_sum <- table_meat_production %>% 
  group_by(Area, M49, Year) %>% 
  summarise(total_production = sum(Value, na.rm = TRUE))

# estimate mean and standard deviation of meat production by country
production_stats <- production_sum %>% 
  group_by(Area, M49) %>% 
  summarise(mean_production = mean(total_production, na.rm = TRUE),
            sd_production = sd(total_production, na.rm = TRUE))


# generating random values from mean and sd for all countries
production_random_values <- rnorm(n = 1000*nrow(production_stats), 
                                  mean = production_stats$mean_production,
                                  sd = production_stats$sd_production) %>%
  matrix(nrow = 1000, ncol = nrow(production_stats), byrow = TRUE) %>%
  as.data.frame() %>%
  setNames(production_stats$Area) %>% 
  pivot_longer(cols = everything(),
               names_to = "Area",
               values_to = "simulated_production")
  
#visualize production distribution by country
ggplot(production_sum, aes(x = total_production)) +
  geom_histogram(bins = 5, fill = "lightblue", color = "black") +
  facet_wrap(~ Area, scales = "free") +
  labs(title = "Distribution of Meat Production by Country",
       x = "Total Meat Production (tonnes)",
       y = "Frequency")
  
#visualizing histograms overlapped 
ggplot() +
  geom_histogram(data = production_sum, aes(x = total_production, y = ..density..),
                 bins = 30, fill = "lightgray", color = "black") +
  geom_density(data = production_random_values, aes(x = simulated_production),
               color = "blue", size = 1) +
  facet_wrap(~ Area, scales = "free") +
  labs(title = "Meat Production Distribution with Simulated Density Overlay",
       x = "Total Meat Production (tonnes)",
       y = "Density")

# Save production stats for future use
write.csv(production_stats, "data/processed/meat_production_stats.csv", row.names = FALSE)
