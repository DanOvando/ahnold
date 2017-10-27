rm(list = ls())
set.seed(123)
library(tidyverse)
library(rstan)
library(stringr)
library(purrr)

run_name <- 'Working'

run_dir <- file.path('results', run_name)

load(file = paste0(run_dir, '/abundance_indices.Rdata'))

data <- abundance_indices %>%
  slice(1:4) %>%
  select(classcode, data) %>%
  unnest()

seen_data <- data %>%
  filter(any_seen == T) %>%
  mutate(mean_vis = (mean_vis - mean(mean_vis)) / (2 * sd(mean_vis))) %>%
  mutate(mean_canopy = (mean_canopy - mean(mean_canopy, na.rm = T)) / (2 * sd(mean_canopy, na.rm = T))) %>%
  mutate(classcode = as.factor(classcode)) %>%
  select(log_density,
         classcode,
         year,
         mean_canopy,
         mean_vis,
         region) %>%
  na.omit()


wtf <-
  rstanarm::stan_glmer(log_density ~ ((factor_year |
                                        classcode) + mean_canopy + mean_vis, data = seen_data))

# lme4::glmer('log_density ~ (factor_year|classcode) + mean_canopy + mean_vis', data = seen_data)

year_species_data <- seen_data %>%
  # filter(classcode %in% c('hcar','hros')) %>%
  select(year, classcode) %>%
  arrange(year, classcode) %>%
  mutate(year_classcode = paste(classcode, year, sep = '-')) %>%
  select(year_classcode) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(year_classcode, marker, fill = 0) %>%
  select(-index)

years_per_species <- data_frame(species = colnames(year_species_data)) %>%
  mutate(position = 1:nrow(.),
         year = str_replace_all(species,'(\\D)',''),
         classcode = str_replace_all(species,'(\\d)|(-)','')) %>%
  group_by(classcode) %>%
  summarise(nyears = length(year)) %>% {
    .$nyears
  }

region_data <- seen_data %>%
  select(region) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(region, marker, fill = 0) %>%
  select(-index)

visibility_data <- seen_data %>%
  select(mean_vis)

canopy_data <- seen_data %>%
  select(mean_canopy)

x <-
  bind_cols(intercept = rep(1, nrow(year_species_data)),
            visibility_data,
            canopy_data,
            year_species_data)

# check <- rstanarm::stan_glmer('log_density ~ (1|factor_year) + (1|region) + mean_canopy + mean_vis', data = seen_data)
# check_binom <- rstanarm::stan_glmer('any_seen ~ (1|factor_year) + (1|region) + mean_canopy + mean_vis', data = seeing_data, family = 'binomial')


stan_data <- list(
  n_parameters = ncol(x),
  n_observations = nrow(x),
  n_year_species = ncol(year_species_data),
  n_species = length(years_per_species),
  years_per_species = years_per_species,
  year_species_positions = which(str_detect(colnames(x), '\\d')),
  x = x,
  log_density = seen_data$log_density
)


ahnold_stan_fit <- stan(
  file = 'scripts/fit-ahnold-allspecies.stan',
  data = stan_data,
  chains = 4,
  warmup = 2000,
  iter = 5000,
  cores = 4,
  refresh = 500
)

seen_data <- seen_data %>%
  mutate(factor_year = year %>% as.character())

check <-
  rstanarm::stan_glmer(log_density ~ (factor_year | classcode) + mean_canopy + mean_vis, data = seen_data)

# control = list(adapt_delta = 0.8)
