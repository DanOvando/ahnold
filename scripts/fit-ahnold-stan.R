rm(list = ls())
set.seed(123)
library(tidyverse)
library(rstan)
library(stringr)
library(purrr)

run_name <- 'Working'

run_dir <- file.path('results', run_name)

load(file = paste0(run_dir, '/abundance_indices.Rdata'))

data <- abundance_indices$data[[1]]

seen_data <- data %>%
  filter(any_seen == T) %>%
  mutate(mean_vis = (mean_vis - mean(mean_vis)) / (2 * sd(mean_vis))) %>%
  mutate(mean_canopy = (mean_canopy - mean(mean_canopy)) / (2 * sd(mean_canopy)))

year_data <- seen_data %>%
  select(year) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(year, marker, fill = 0) %>%
  select(-index)

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

x <- bind_cols(intercept = rep(1, nrow(year_data)),visibility_data,canopy_data,year_data, region_data)

seeing_data <- data %>%
  mutate(mean_vis = (mean_vis - mean(mean_vis)) / (2 * sd(mean_vis))) %>%
  mutate(mean_canopy = (mean_canopy - mean(mean_canopy)) / (2 * sd(mean_canopy)))

year_data <- seeing_data %>%
  select(year) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(year, marker, fill = 0) %>%
  select(-index)

region_data <- seeing_data %>%
  select(region) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(region, marker, fill = 0) %>%
  select(-index)

visibility_data <- seeing_data %>%
  select(mean_vis)

canopy_data <- seeing_data %>%
  select(mean_canopy)

x_obs <- bind_cols(intercept = rep(1, nrow(year_data)),visibility_data,canopy_data,year_data, region_data)

# check <- rstanarm::stan_glmer('log_density ~ (1|factor_year) + (1|region) + mean_canopy + mean_vis', data = seen_data)
# check_binom <- rstanarm::stan_glmer('any_seen ~ (1|factor_year) + (1|region) + mean_canopy + mean_vis', data = seeing_data, family = 'binomial')


stan_data <- list(
  n_parameters = ncol(year_data) + ncol(region_data) + 3,
  n_observations = nrow(seen_data),
  n_observations_obs = nrow(seeing_data),
  n_years = ncol(year_data),
  n_regions = ncol(region_data),
  year_positions = which(str_detect(colnames(x),'20')),
  region_positions = which(colnames(x) %in% colnames(region_data)),
  intercept_position = 1,
  vis_position = 2,
  kelp_position = 3,
  x_obs = x_obs,
  x = x,
  log_density = seen_data$log_density,
  observed = as.numeric(seeing_data$any_seen)
)


ahnold_stan_fit<- stan(
  file = 'scripts/fit-ahnold-joint.stan',
  data = stan_data,
  chains = 4,
  warmup = 2000,
  iter = 5000,
  cores = 4,
  refresh = 500
)
# control = list(adapt_delta = 0.8)
