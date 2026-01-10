#Merging probabilities for live imports and product imports

#import tables
live_animals_risk <- read.csv('results/live_animals_risk.csv')
products_risk <- read.csv('results/products_risk.csv')

#calculating union of probabilities
# PRPSA=PRfL ∪ PRfP= PRfL+PRfP-(PRfL*PRfP).      
pr <- live_animals_risk %>% 
  filter(Países == 'All') %>% 
  left_join(products_risk %>% filter(Países == 'All'), by = 'Países') %>%
  mutate(pr_mean = prL_mean + prP_mean - (prL_mean * prP_mean),
         pr_025 = prL_025 + prP_025 - (prL_025 * prP_025),
         pr_975 = prL_975 + prP_975 - (prL_975 * prP_975)) %>% 
  select(Países, pr_mean, pr_025, pr_975)
  
write.csv(pr, 'results/union_risk.csv', row.names = FALSE)
