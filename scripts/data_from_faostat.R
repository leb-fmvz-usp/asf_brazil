#load libraries
library(tidyverse)
library(readr)

#show files in zip file in downloads folder
unz <- "data/raw/Production_Crops_Livestock_E_All_Data.zip"
files <- unzip(unz, list = TRUE)
print(files)
#read in a specific file from the zip
faodata <- read.csv(unzip(unz, 'Production_Crops_Livestock_E_All_Data.csv'))

#filter for pig products
pig_products <- faodata %>% 
  filter(grepl("pig", Item, ignore.case = TRUE)) %>% 
  filter(!grepl("pigeon", Item, ignore.case = TRUE)) %>% 
  filter(Element == 'Production') %>% 
  rename(M49 = Area.Code..M49.) %>% 
  mutate(M49 = parse_number(M49))

#select columns M49 to Unit and years 2019 to 2023
pig_products <- pig_products %>% 
  select(M49:Unit, num_range("Y", 2019:2023))

#sum of all pork production each year
pig_products <- pig_products %>% 
  group_by(Area, M49) %>% 
  summarise(across(starts_with("Y"), ~sum(.x, na.rm = TRUE)))

#calculate mean and sd of pork production by country
pig_products_stats <- pig_products %>%
  group_by(Area, M49) %>% 
  summarise(mean_production = mean(c_across(starts_with("Y")), na.rm = TRUE),
            sd_production = sd(c_across(starts_with("Y")), na.rm = TRUE))

#save table
write.csv(pig_products_stats, "data/processed/meat_production_faostat.csv", row.names = FALSE)

