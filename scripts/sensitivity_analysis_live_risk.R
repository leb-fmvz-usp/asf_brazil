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
parameters <- c('P1','hp','ou','no','to','P3','P4')
countries <- imports$M49
### function to generate simulation array ----
generate_array <- function(){
  simulation_array <- array(data = NA, 
                            dim = c(length(countries), length(parameters),n_sim ), 
                            dimnames = list(countries,parameters,1:n_sim))
  
  #P1
  simulation_array[ as.character(countries), 'P1', ] <- 
    rpert(n = n_sim*length(countries), 
          min = imports$min, 
          mode = imports$mais_prov, 
          max = imports$max) 
  
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
  # simulation_array[ as.character(countries), 'ou', ] <- 
  #   rpert(n = n_sim*length(countries), 
  #         min = min(tabela_ou$Casos), 
  #         mode = 2,
  #         max = 10)
  #ou <- rpert(n = 1, min = 1, mode = 1.28,max = 6) #Beatriz data
  
  simulation_array[ as.character(countries), 'no', ] <- 
    rnorm(n = n_sim*length(countries), 
          mean = imports$mean_animals, 
          sd = 2*imports$sd_animals)
  
  ##To - herd size: Population / Number of Herds
  simulation_array[ as.character(countries), 'to', ] <- 
    simulation_array[ as.character(countries), 'no', ] / imports$so
  
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
  # NI - Number of Infected: Outbreaks * herd size * prevalence
  NI <- X$ou * X$to * X$hp
  # Beta distribution for P2L Beta(α1, α2)
  #alpha1  = NI + 1
  alpha1 <- NI + 1
  #alpha2 = no - (NI + 1)
  alpha2 <- X$no - alpha1
  #P2L
  P2L <- rbeta(n = 1, shape1 = alpha1, shape2 = alpha2)
  pcL <- X$P1 * P2L * X$P3 * X$P4
  return(pcL)
}

# Perform Sobol sensitivity analysis
sobol_result <- sobol2007(model = model_function, X1 = X1, X2 = X2, nboot = 100)
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

write.csv(df_csv, "results/sobol_indices_live_animals.csv", row.names = FALSE)

# Tornado plot
ggplot(df_long, aes(x = reorder(param, valor), y = valor, fill = tipo)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(x = "Parameter", y = "Sobol Indices", fill = "",
       title = "Tornado plot for S and ST (live animals)") +
  theme_bw()
# save plot in results
ggsave("results/sobol_tornado_plot_live_animals.png", width = 8, height = 6)

##  Sensitivy analysis with mixed effects model ----
#load simulation array
simulation_array <- readRDS('data/processed/simulation_live_animals_array.rds')

#Rearrange array into a data frame
df <- as.data.frame.table(simulation_array, responseName = "value") %>%
  rename(
    country = Var1,
    parameter = Var2,
    simulation = Var3
  ) %>%
  pivot_wider(names_from = parameter, values_from = value)

#model
model_lmer <- lmer(prL ~  P4 + P3 + ou + to + hp + no + P1 + (1|country), data = df)
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

