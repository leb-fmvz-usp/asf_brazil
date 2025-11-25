# Script to calculate mean and sd of pig population from FAOSTAT data (2013-2023)
library(tidyverse)

# load data
tabela_no <- read.csv("data/raw/No_pig_population_2013-2023_faostat.csv")

#For each year in the tabela_no, add an entry with M49 = 156, representing China mainland, 
#the Diff between China (M49 = 159) and China Taiwan (M49 = 158)
tabela_no <- tabela_no %>% 
  filter(Area.Code..M49. %in% c(158,159)) %>%
  group_by(Year) %>%
  summarise(Value = as.numeric(Value[Area.Code..M49. == 159]) - as.numeric(Value[Area.Code..M49. == 158])) %>%
  mutate(Area.Code..M49. = 156,
         Area..M49. = "China, mainland") %>%
  select(Area..M49., Area.Code..M49., Year, Value) %>%
  bind_rows(tabela_no, .)

#recalculate mean and sd including China mainland
tabela_no <- tabela_no %>% 
  rename(M49 = Area.Code..M49.) %>% 
  group_by(M49) %>%
  summarise(mean_animals = mean(as.numeric(Value),na.rm = TRUE),
            sd_animals = sd(as.numeric(Value),na.rm = TRUE))
#save table
write.csv(tabela_no, "data/processed/pig_population_country.csv", row.names = FALSE)


