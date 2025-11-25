#Scenario-tree for imports of live animals

#loading packages
library(tidyverse)
library(readxl)
library(mc2d)
library(lme4)

## loading data ----
p1_table <- read.csv("data/processed/p1_table.csv")
imports <- read_xlsx("data/processed/imports_table_2023.xlsx") %>%#translate this table
  filter(grepl("vivo", tipo, ignore.case = TRUE))
tabela_hp <- read.csv("data/processed/Hp_EU_PSA_prevalence_2013-2023.csv") #translate this table
tabela_so <- read_xlsx("data/processed/So_population_establishment_exporting_countries_WOAH.xlsx")
tabela_no <- read.csv("data/raw/No_pig_population_2013-2023_faostat.csv")
tabela_ou <- read_xlsx("data/processed/Ou_high_risk_period_cases.xlsx")

#import P1
imports <- imports %>% 
  left_join(p1_table)

## So - Premises
imports <- imports %>% 
  left_join(tabela_so, by = join_by(ISO3 == 'Country(ISO3)') )
imports$so <- imports$`Nb animals premises`

#no - number of animals
imports <- imports %>% 
  left_join( tabela_no %>% 
               mutate(M49 = Area.Code..M49.) %>% 
               group_by(M49) %>% 
               summarise(mean_animals = mean(as.numeric(Value),na.rm = TRUE),
                         sd_animals = sd(as.numeric(Value),na.rm = TRUE))
  )
## number of imported pigs
imports <- imports %>% 
  rename(ncL = `2023 - Quilograma Líquido`)

### simulation parameters ----
n_sim <- 1000
parameters <- c('P1','hp','ou','no','to','NI','alpha1','alpha2','P2L','P3','P4','pcL','prL')
countries <- imports$M49


### generating simulation array ----
simulation_array <- array(data = NA, dim = c(length(countries), length(parameters),n_sim ), 
                          dimnames = list(countries,parameters,1:n_sim))

#P1
simulation_array[ as.character(countries), 'P1', ] <- 
  rpert(n = n_sim*length(countries), min = imports$min, mode = imports$mais_prov, max = imports$max) 

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
        mode = mean(tabela_ou$Casos),
        max = max(tabela_ou$Casos))
simulation_array[ as.character(countries), 'ou', ] <- 
  rpert(n = n_sim*length(countries), 
        min = min(tabela_ou$Casos), 
        mode = 2,
        max = 10)
#ou <- rpert(n = 1, min = 1, mode = 1.28,max = 6) #Beatriz data

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

# Beta distribution for P2L Beta(α1, α2)
#alpha1  = NI + 1
simulation_array[ as.character(countries), 'alpha1', ] <- 
  simulation_array[ as.character(countries), 'NI', ] + 1
#alpha2 = no - (NI + 1)
simulation_array[ as.character(countries), 'alpha2', ] <- 
  simulation_array[ as.character(countries), 'no', ] - 
  simulation_array[ as.character(countries), 'alpha1', ]

#P2L
simulation_array[ as.character(countries), 'P2L', ] <- 
  rbeta(n = n_sim*length(countries), 
        shape1 = simulation_array[ as.character(countries), 'alpha1', ], 
        shape2 = simulation_array[ as.character(countries), 'alpha2', ])

## P3 
simulation_array[ as.character(countries), 'P3', ] <- 
  rpert(n = n_sim*length(countries), 
        min = 0.206, 
        mode = 0.633, 
        max = 1)

## P4 
simulation_array[ as.character(countries), 'P4', ] <- 
  rpert(n = n_sim*length(countries), 
        min=0.0005, 
        mode=0.0027, 
        max=0.092)


### Calculating probabilities ----

#probability of release for each live pig imported
simulation_array[ as.character(countries), 'pcL', ] <- 
  simulation_array[ as.character(countries), 'P1', ] *
  simulation_array[ as.character(countries), 'P2L', ] *
  simulation_array[ as.character(countries), 'P3', ] *
  simulation_array[ as.character(countries), 'P4', ] 

# probability of release of at least one live pig
simulation_array[ as.character(countries), 'prL', ] <- 
  1 - (1 - simulation_array[ as.character(countries), 'pcL', ])^(imports$ncL)
  
#prL by country
prL_mean <- apply(FUN = mean, X = simulation_array[,'prL',], MARGIN = 1)
prL_025 <- apply(FUN = quantile, X = simulation_array[,'prL',], MARGIN = 1, probs = 0.025)
prL_975 <- apply(FUN = quantile, X = simulation_array[,'prL',], MARGIN = 1, probs = 0.975)
prL <- data.frame(M49 = as.numeric(names(prL_mean)), prL_mean, prL_025, prL_975)

#merge results with imports table
results <- select(imports, M49, ISO3, `Países`) %>% 
  left_join(prL)

#probability of release of at least one live pig combining all countries
print_results <- results %>%
  bind_rows(
    results %>% 
      summarise(across(prL_mean:prL_975, function(.x) {1 - prod( 1- .x)})) %>%
      mutate(M49 = NA, ISO3 = '', `Países` = 'All')
  )

print_results <- print(print_results %>% 
        mutate(across(prL_mean:prL_975, 
                      ~ formatC(.x, format = "e", digits = 2))))
                      

#save simulation array
saveRDS(simulation_array, "data/processed/simulation_live_animals_array.rds")

#save results
write.csv(print_results, "results/live_animals_risk.csv", row.names = FALSE)



