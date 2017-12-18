set.seed(666)
library(tidyverse)
library(rstan)
library(stringr)
library(purrr)
# rstan_options(auto_write = TRUE)
demons::load_functions()
run_name <- 'Working'

run_dir <- file.path('results', run_name)

load(file = paste0(run_dir, '/abundance_indices.Rdata'))

subspecies <- abundance_indices %>%
  select(classcode, targeted) %>%
  mutate(targeted = as.numeric(targeted > 0)) %>%
  mutate(a = 1:nrow(.)) %>%
  group_by(targeted) %>%
  top_n(2,a)

data <- abundance_indices %>%
  filter(population_structure == 'one-pop',
         population_filtering == 'all',
         data_source == 'length_to_density') %>%
  select(classcode, data) %>%
  unnest() %>%
  # filter(classcode %in% subspecies$classcode) %>%
  mutate(targeted = as.numeric(targeted > 0)) %>%
  group_by(classcode) %>%
  mutate(mean_vis = (mean_vis - mean(mean_vis)) / (2 * sd(mean_vis))) %>%
  mutate(mean_canopy = (mean_canopy - mean(mean_canopy, na.rm = T)) / (2 * sd(mean_canopy, na.rm = T))) %>%
  ungroup()

numeric_species_key <- data_frame(classcode = unique(data$classcode)) %>%
  arrange(classcode) %>%
  mutate(numeric_classcode = 1:nrow(.))

raw_length_covars <-
  c('classcode',
    'factor_year',
    'region',
    'targeted',
    'zone',
    'site_side',
    'level',
    'mean_vis',
    'surge',
    'factor_month',
    'trunc_observer',
    'cumulative_n_obs',
    'cumulative_n_obs_2')

seen_data <- data %>%
  filter(any_seen == T) %>%
  left_join(numeric_species_key, by = 'classcode')

seeing_data <- data %>%
  left_join(numeric_species_key, by = 'classcode')


non_nested_variables <- c(
  'zone',
  'site_side',
  'level',
  'mean_vis',
  'surge',
  'factor_month',
  'trunc_observer',
  'cumulative_n_obs'
)


# seen data ---------------------------------------------------------------

x_seen_non_nested <- seen_data %>%
  select(non_nested_variables) %>%
  mutate(intercept = 1) %>%
  spread_factors(drop_one = T) %>%
  purrrlyr::dmap(center_scale) %>%
  mutate(cumulative_n_obs_2 = cumulative_n_obs^2)

year_species_data <- seen_data %>%
  mutate(year_classcode = paste(classcode, year, sep = '-')) %>%
  select(year_classcode) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(year_classcode, marker, fill = 0) %>%
  select(-index)

year_fixed_effects <- seen_data %>%
  select(year) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(year, marker, fill = 0) %>%
  select(-index) %>%
  select(-1)

years_per_species <-
  data_frame(species = colnames(year_species_data)) %>%
  mutate(
    position = 1:nrow(.),
    year = str_replace_all(species, '(\\D)', ''),
    classcode = str_replace_all(species, '(\\d)|(-)', '')
  ) %>%
  group_by(classcode) %>%
  summarise(nyears = length(year)) %>% {
    .$nyears
  }

region_species_data <- seen_data %>%
  select(region, classcode) %>%
  mutate(region_classcode = paste(classcode, region, sep = '-')) %>%
  select(region_classcode) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(region_classcode, marker, fill = 0) %>%
  select(-index)

region_geographic_cluster_data <- seen_data %>%
  select(region, geographic_cluster) %>%
  mutate(region_cluster = paste(geographic_cluster, region, sep = '-')) %>%
  select(region_cluster) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(region_cluster, marker, fill = 0) %>%
  select(-index)

regions_per_species <- seen_data %>%
  group_by(classcode) %>%
  summarise(nr = length(unique(region))) %>% {
    .$nr
  }

regions_per_cluster <- seen_data %>%
  group_by(geographic_cluster) %>%
  summarise(nr = length(unique(region))) %>% {
    .$nr
  }

n_intercepts <- length(years_per_species)

intercept_data <- seen_data %>%
  select(classcode) %>%
  mutate(seen = 1, id = 1:nrow(.)) %>%
  spread(classcode, seen, fill = 0) %>%
  select(-id) %>%
  set_names(glue::glue('{colnames(.)}_intercept')) %>%
  select(-1)

x_seen <-
  bind_cols(
    x_seen_non_nested,
    year_species_data,
    region_geographic_cluster_data
  )


# seeing data -------------------------------------------------------------

x_seeing_non_nested <- seeing_data %>%
  select(non_nested_variables) %>%
  mutate(intercept = 1) %>%
  spread_factors(drop_one = T) %>%
  purrrlyr::dmap(center_scale) %>%
  mutate(cumulative_n_obs_2 = cumulative_n_obs^2)

year_species_data <- seeing_data %>%
  mutate(year_classcode = paste(classcode, year, sep = '-')) %>%
  select(year_classcode) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(year_classcode, marker, fill = 0) %>%
  select(-index)

year_fixed_effects <- seeing_data %>%
  select(year) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(year, marker, fill = 0) %>%
  select(-index) %>%
  select(-1)

years_per_species <-
  data_frame(species = colnames(year_species_data)) %>%
  mutate(
    position = 1:nrow(.),
    year = str_replace_all(species, '(\\D)', ''),
    classcode = str_replace_all(species, '(\\d)|(-)', '')
  ) %>%
  group_by(classcode) %>%
  summarise(nyears = length(year)) %>% {
    .$nyears
  }

region_species_data <- seeing_data %>%
  select(region, classcode) %>%
  mutate(region_classcode = paste(classcode, region, sep = '-')) %>%
  select(region_classcode) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(region_classcode, marker, fill = 0) %>%
  select(-index)

region_geographic_cluster_data <- seeing_data %>%
  select(region, geographic_cluster) %>%
  mutate(region_cluster = paste(geographic_cluster, region, sep = '-')) %>%
  select(region_cluster) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(region_cluster, marker, fill = 0) %>%
  select(-index)

regions_per_species <- seeing_data %>%
  group_by(classcode) %>%
  summarise(nr = length(unique(region))) %>% {
    .$nr
  }

regions_per_cluster <- seeing_data %>%
  group_by(geographic_cluster) %>%
  summarise(nr = length(unique(region))) %>% {
    .$nr
  }

n_intercepts <- length(years_per_species)

intercept_data <- seeing_data %>%
  select(classcode) %>%
  mutate(seen = 1, id = 1:nrow(.)) %>%
  spread(classcode, seen, fill = 0) %>%
  select(-id) %>%
  set_names(glue::glue('{colnames(.)}_intercept')) %>%
  select(-1)

x_seeing <-
  bind_cols(
    x_seeing_non_nested,
    year_species_data,
    region_geographic_cluster_data
  )


# prepare data for stan ---------------------------------------------------

x_seen_species_index <-  seen_data$classcode %>% as.factor() %>%  as.numeric()

log_density_species_index <- seen_data %>%
  select(classcode) %>%
  mutate(index = 1:nrow(.)) %>%
  group_by(classcode) %>%
  summarise(start_position = min(index),
            end_position = max(index)) %>%
  arrange(start_position) %>%
  select(-classcode) %>%
  as.matrix()

year_species_positions <- which(str_detect(colnames(x_seen), '20') & str_detect(colnames(x_seen),'-'))

region_cluster_positions <- which(str_detect(
  colnames(x_seen), paste(unique(data$region), collapse = '|')
) & !str_detect(colnames(x_seen),'site_side'))

non_nested_positions <- which(!str_detect(colnames(x_seen), '-') & !str_detect(
  colnames(x_seen), paste(unique(data$region), collapse = '|')
))

stan_data <- list(
  x_seen = x_seen,
  n_parameters = ncol(x_seen),
  n_observations_seen = nrow(x_seen),
  n_species = length(years_per_species),
  n_clusters = length(regions_per_cluster),
  n_year_species = ncol(year_species_data),
  n_region_clusters = ncol(region_geographic_cluster_data),
  years_per_species = years_per_species,
  regions_per_cluster = regions_per_cluster,
  year_species_positions = year_species_positions,
  region_cluster_positions = region_cluster_positions,
  non_nested_positions = non_nested_positions,
  n_non_nested = length(non_nested_positions),
  log_density = seen_data$log_density,
  log_density_species_index = log_density_species_index,
  x_seen_species_index = x_seen_species_index,
  n_observations_seeing = nrow(x_seeing),
  x_seeing = x_seeing,
  observed = as.numeric(seeing_data$any_seen)
)

stan_data <- map_if(stan_data, is.data.frame,~as.matrix(.x))


set.seed(666)
a <- Sys.time()
ahnold_stan_fit <- stan(
  file = 'scripts/fit-ahnold-hurdle.stan',
  data = stan_data,
  chains = 1,
  warmup = 1000,
  iter = 2000,
  cores = 1,
  refresh = 25,
  control = list(max_treedepth = 9))

Sys.time() - a

# ,
# control = list(max_treedepth = 7))
save(file = 'cluster-seen-only.Rdata', ahnold_stan_fit)





