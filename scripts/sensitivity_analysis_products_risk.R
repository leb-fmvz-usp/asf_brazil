#Sensitivity analysis - 

#loading packages
library(tidyverse)
library(readxl)
library(mc2d)
library(lme4)
library(sensitivity)

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
  u <- runif(n, q[1], q[2])
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

## number of imported pigs
imports <- imports %>% 
  rename(ncP = `2023 - Quilograma Líquido`)

## meat production stats
imports <- imports %>% 
  left_join(meat_production_stats %>% 
              select(M49, mean_production, sd_production),
            by = join_by(M49)
  )

### simulation parameters ----
n_sim <- 1000
parameters <- c('P1','hp','ou','no','to','NI','alpha1','alpha2','P2L','P3','P4','pcL','prL',
                'Pus','Psm','Pm','Mp','Nm','Qim','alpha1p','alpha2p','P2P','pcP','PRp')
countries <- imports$M49

### function to generate simulation array ----
generate_array <- function(){
  
  ### generating simulation array ----
  simulation_array <- array(data = NA, 
                            dim = c(length(countries), length(parameters),n_sim ), 
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
  
  
  ## Pus (check needed) ##
  
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
    simulation_array[ as.character(countries), 'alpha1', ]
  
  #P2P
  simulation_array[ as.character(countries), 'P2P', ] <- 
    rbeta(n = n_sim*length(countries), 
          shape1 = simulation_array[ as.character(countries), 'alpha1p', ], 
          shape2 = simulation_array[ as.character(countries), 'alpha2p', ])
  
  # #PCP
  # simulation_array[ as.character(countries), 'pcP', ] <- 
  #   simulation_array[ as.character(countries), 'P1', ] *
  #   simulation_array[ as.character(countries), 'P2P', ]

  return(simulation_array)
}
### preparing data for sensitivity analysis ----
#drop first dimension of the simulation_array [countries]
X1 <- apply(generate_array(), 2, c)
X2 <- apply(generate_array(), 2, c)
#convert to data frame
X1 <- as.data.frame(X1)
X2 <- as.data.frame(X2)

#prepare a Sensitivity analysis with Sobol
# Define the model function
model_function <- function(X) {
  alpha1p <- X$Qim + 1
  alpha2p <- X$Nm - alpha1p
  pcP <- X$P1 * X$P2P
  return(pcP)
}

# Perform Sobol sensitivity analysis
sobol_result <- sobolmartinez(model = model_function, X1 = X1, X2 = X2, nboot = 100)
#graphical representation
ggplot(sobol_result)

# Print the results
print(sobol_result)

S  <- as.data.frame(sobol_result$S)
ST <- as.data.frame(sobol_result$T)

S  <- S %>% 
  rownames_to_column("param") %>% 
  select(param, S = original)

ST <- ST %>% 
  rownames_to_column("param") %>% 
  select(param, ST = original)

df_plot <- S %>%
  left_join(ST, by = "param")

df_long <- df_plot %>%
  pivot_longer(cols = c(S, ST),
               names_to = "tipo",
               values_to = "valor")

#save df_plot in table in csv format, with two decimals and percentage format
df_csv <- df_plot %>%
  mutate(
    S = sprintf("%.2f%%", S * 100),
    ST = sprintf("%.2f%%", ST * 100)
  )

write.csv(df_csv, "results/sobol_indices_products.csv", row.names = FALSE)

# Tornado plot
ggplot(df_long, aes(x = reorder(param, valor), y = valor, fill = tipo)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(x = "Parameter", y = "Sobol Indices", fill = "",
       title = "Tornado plot for S and ST (live animals)") +
  theme_bw()
# save plot in results
ggsave("results/sobol_tornado_plot_products.png", width = 8, height = 6)

##  Sensitivy analysis with mixed effects model ----
#load simulation array
simulation_array <- readRDS('data/processed/simulation_products_array.rds')

#Rearrange array into a data frame
df <- as.data.frame.table(simulation_array, responseName = "value") %>%
  rename(
    country = Var1,
    parameter = Var2,
    simulation = Var3
  ) %>%
  pivot_wider(names_from = parameter, values_from = value)

#model
model_lmer <- lmer(pcP ~  P2P + Qim + Nm + P4 + P3 + ou + to + hp + no + P1 + (1|country), data = df)
summary(model_lmer)

# Extract model output
s <- summary(model_lmer)

# Fixed effect table
tab <- as.data.frame(s$coefficients) %>%
  rownames_to_column("Parameter") %>%
  mutate(
    Est  = Estimate,
    SE   = `Std. Error`,
    Lower = Est - 1.96 * SE,
    Upper = Est + 1.96 * SE
  )

# Graph
ggplot(tab, aes(x = reorder(Parameter, Est), y = Est)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15) +
  coord_flip() +
  labs(
    x = "",
    y = "Coefficient (± 1.96 * SE)",
    title = "Fixed effects"
  ) +
  theme_bw(base_size = 14)

