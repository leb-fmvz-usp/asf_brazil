#Scenario-tree for imports of live animals

#loading packages
library(tidyverse)
library(readxl)
library(mc2d)
library(lme4)

## loading data ----
p1_table <- read.csv("data/processed/p1_table.csv")
imports <- read_xlsx("data/processed/imports_table_2023.xlsx") %>%#translate this table
  filter(grepl("produto", tipo, ignore.case = TRUE))
tabela_hp <- read.csv("data/processed/Hp_EU_PSA_prevalence_2013-2023.csv") #translate this table
tabela_so <- read_xlsx("data/processed/So_population_establishment_exporting_countries_WOAH.xlsx")
tabela_ou <- read_xlsx("data/processed/Ou_high_risk_period_cases.xlsx")
table_psm <- read.csv('data/processed/psm_estimated.csv')
tons_per_pig_estimated <- read.csv("data/processed/tons_per_pig_estimated.csv")
meat_production_stats <- read.csv("data/processed/meat_production_faostat.csv")
pig_population_country <- read.csv("data/processed/pig_population_country.csv")

# Function to generate truncated normal distribution
truncnorm <- function(n, mean, sd, min, max) {
  q <- pnorm(c(min, max), mean = mean, sd = sd)
  u <- runif(n, min = q[1], max = q[2])
  return(qnorm(u, mean = mean, sd = sd))
}

#import P1
imports <- imports %>% 
  left_join(p1_table)

## So - Premises
imports <- imports %>% 
  left_join(tabela_so, by = join_by(ISO3 == 'Country(ISO3)') )
imports$so <- imports$`Nb animals premises`

#no - number of animals
imports <- imports %>% 
  left_join(pig_population_country)

## number of imported tons of productss
imports <- imports %>% 
  rename(ncP = `2023 - Quilograma Líquido`) %>% 
  mutate(ncP = ncP / 1000)

## meat production stats
imports <- imports %>% 
  left_join(meat_production_stats %>% 
               select(M49, mean_production, sd_production),
             by = join_by(M49)
  )

### simulation parameters ----
n_sim <- 1000
parameters <- c('P1','hp','ou','no','to','NI','P3','P4',
                'Pus','Psm','Pm','Mp','Nm','Qim','alpha1p','alpha2p',
                'P2P','pcP','prP')
countries <- imports$M49


### generating simulation array ----
simulation_array <- array(data = NA, dim = c(length(countries), length(parameters),n_sim ), 
                          dimnames = list(countries,parameters,1:n_sim))

#P1: countries with no history of ASF, probability from SEIR model
countries_no_ASF <- imports %>% 
  filter(is.na(casos_PSA))
simulation_array[ as.character(countries_no_ASF$M49), 'P1', ] <- 
  rpert(n = n_sim*nrow(countries_no_ASF), 
        min = countries_no_ASF$min, 
        mode = countries_no_ASF$mais_prov, 
        max = countries_no_ASF$max) 

#P1: countries with history of ASF, probability from Poisson model
countries_with_ASF <- imports %>% 
  filter(casos_PSA == 'P')
simulation_array[ as.character(countries_with_ASF$M49), 'P1', ] <- 
  rpert(n = n_sim*nrow(countries_with_ASF), 
        min = countries_with_ASF$min, 
        mode = countries_with_ASF$mais_prov, 
        max = countries_with_ASF$max) 

#HP - Herd Prevalence
simulation_array[ as.character(countries), 'hp', ] <- 
  rpert(n = n_sim*length(countries), 
        min = min(tabela_hp$prevalecia), 
        mode =sum(tabela_hp$prevalecia)/nrow(tabela_hp),
        max = max(tabela_hp$prevalecia))

#OU - Outbreaks
simulation_array[ as.character(countries), 'ou', ] <- 
  rpert(n = n_sim*length(countries), 
        min = min(tabela_ou$Casos), 
        mode = 1, #mode of outbreaks, check with table(tabela_ou$Casos)
        max = 231) #exclding extreme value of max outbreaks in Romenia

# number of animals
simulation_array[ as.character(countries), 'no', ] <- 
  rnorm(n = n_sim*length(countries), 
        mean = imports$mean_animals, 
        sd = imports$sd_animals)

##To - herd size: Population / Number of Herds
simulation_array[ as.character(countries), 'to', ] <- 
  simulation_array[ as.character(countries), 'no', ] / imports$so

# NI - Number of Infected: Outbreaks * herd size * prevalence
simulation_array[ as.character(countries), 'NI', ] <- 
  simulation_array[ as.character(countries), 'ou', ] *
  simulation_array[ as.character(countries), 'to', ] *
  simulation_array[ as.character(countries), 'hp', ]

## P3 #chance of an infected pig survive the infection period
simulation_array[ as.character(countries), 'P3', ] <- 
  rpert(n = n_sim*length(countries), 
        min = 0.206, 
        mode = 0.633, 
        max = 1)

## P4 #chance of survival during transit
simulation_array[ as.character(countries), 'P4', ] <- 
  1 - rpert(n = n_sim*length(countries), 
            min=0.0005, 
            mode=0.0027, 
            max=0.092)

## Pus (check needed) ##
# Probability of an infected pig pass unreported slaughtering process
simulation_array[ as.character(countries), 'Pus', ] <- 
  rbeta(n = n_sim*length(countries),
        shape1 = 1.34, 
        shape2 = 34.17)

## Psm ##
# Probability of a pig being slaughtered for meat production
simulation_array[ as.character(countries), 'Psm', ] <- 
  truncnorm(n = n_sim*length(countries),
            mean = table_psm$mean_prob,
            sd = table_psm$sd_prob,
            min = table_psm$min,
            max = table_psm$max)

## Pm ##
# Probability of an infected pig being slaughtered and used for meat production
simulation_array[ as.character(countries), 'Pm', ] <- 
  simulation_array[ as.character(countries), 'P3', ] *
  simulation_array[ as.character(countries), 'P4', ] *
  simulation_array[ as.character(countries), 'Psm', ] *
  simulation_array[ as.character(countries), 'Pus', ]

## Mp ##
# Average kg meat obtained from pig (in tons)
simulation_array[ as.character(countries), 'Mp', ] <- 
  truncnorm(n = n_sim*length(countries),
            mean = tons_per_pig_estimated$mean_tons_per_pig,
            sd = tons_per_pig_estimated$sd_tons_per_pig,
            min = tons_per_pig_estimated$min,
            max = tons_per_pig_estimated$max)

# Qim
# Estimated infected meat produced
simulation_array[ as.character(countries), 'Qim', ] <- 
  simulation_array[ as.character(countries), 'NI', ] *
  simulation_array[ as.character(countries), 'Pm', ] *
  simulation_array[ as.character(countries), 'Mp', ]

## Nm
# Total meat production in exporting country
simulation_array[ as.character(countries), 'Nm', ] <- 
  rnorm(n = n_sim*length(countries),
        mean = imports$mean_production,
        sd = imports$sd_production)

## P2P ##
# Beta distribution for P2P Beta(α1p, α2p)
#alpha1  = QIM + 1
simulation_array[ as.character(countries), 'alpha1p', ] <- 
  simulation_array[ as.character(countries), 'Qim', ] + 1
#alpha2 = NM - (QIM + 1)
simulation_array[ as.character(countries), 'alpha2p', ] <- 
  simulation_array[ as.character(countries), 'Nm', ] - 
  simulation_array[ as.character(countries), 'alpha1p', ]

#P2P
simulation_array[ as.character(countries), 'P2P', ] <- 
  rbeta(n = n_sim*length(countries), 
        shape1 = simulation_array[ as.character(countries), 'alpha1p', ], 
        shape2 = simulation_array[ as.character(countries), 'alpha2p', ])

#PCP
simulation_array[ as.character(countries), 'pcP', ] <- 
  simulation_array[ as.character(countries), 'P1', ] *
  simulation_array[ as.character(countries), 'P2P', ]

#PRp
# probability of release of at least one contaminated meat product
simulation_array[ as.character(countries), 'prP', ] <- 
  1 - (1 - simulation_array[ as.character(countries), 'pcP', ])^(imports$ncP)


#prP by country
prP_mean <- apply(FUN = mean, X = simulation_array[,'prP',], MARGIN = 1)
prP_025 <- apply(FUN = quantile, X = simulation_array[,'prP',], MARGIN = 1, probs = 0.025)
prP_975 <- apply(FUN = quantile, X = simulation_array[,'prP',], MARGIN = 1, probs = 0.975)
prP <- data.frame(M49 = as.numeric(names(prP_mean)), prP_mean, prP_025, prP_975)

#merge results with imports table
results_meat <- select(imports, M49, ISO3, `Países`) %>% 
  left_join(prP)

#visualize Prp and its interval by country
ggplot(results_meat, aes(x = reorder(`Países`, -prP_mean), y = prP_mean)) +
  geom_point(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = prP_025, ymax = prP_975), width = 0.2) +
  coord_flip() +
  labs(title = "Probability of release of at least one contaminated meat product by country",
       x = "Country",
       y = "Probability of release (prP)") +
  theme_minimal()


#probability of release of at least one contaminated meat product combining all countries
print_results_meat <- results_meat %>%
  bind_rows(
    results_meat %>% 
      summarise(across(prP_mean:prP_975, function(.x) {1 - prod( 1- .x)})) %>%
      mutate(M49 = NA, ISO3 = '', `Países` = 'All')
  )

#format results in scientific notation
print_results_meat <- print(print_results_meat %>% 
                              mutate(across(prP_mean:prP_975, 
                                            ~ formatC(.x, format = "e", digits = 2))))

#save results
write.csv(print_results_meat, "results/products_risk.csv", row.names = FALSE)

#save simulation array
saveRDS(simulation_array, "data/processed/simulation_products_array.rds")

#analyze p2p for italy
p2p_italy <- simulation_array[ as.character(380), 'P2P', ]
hist(p2p_italy, breaks = 50, main = "P2P Italy", xlab = "P2P")
#analyze pcp for italy
pcp_italy <- simulation_array[ as.character(380), 'pcP', ]
hist(pcp_italy, breaks = 50, main = "PcP Italy", xlab = "PcP")
#analyze prp for italy
prp_italy <- simulation_array[ as.character(380), 'prP', ]
hist(prp_italy, breaks = 50, main = "PRp Italy", xlab = "prP")
#analyze ncP from imports table
hist(imports$ncP)
#visualize ncp by country
imports %>%
  ggplot(aes(x = reorder(`Países`, -ncP), y = ncP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Number of imported pigs by country",
       x = "Country",
       y = "Number of imported pigs")
