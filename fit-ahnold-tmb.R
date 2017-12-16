set.seed(666)
library(TMB)
library(stringr)
library(purrr)
library(tidyverse)

# rstan_options(auto_write = TRUE)

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


# prepare seen ------------------------------------------------------------



seen_data <- data %>%
  filter(any_seen == T) %>%
  select(log_density,
         classcode,
         geographic_cluster,
         year,
         month,
         mean_vis,
         region,
         trunc_observer,
         cumulative_n_obs,
         method,
         level,
         surge) %>%
  na.omit() %>%
  left_join(numeric_species_key, by = 'classcode')

x_seen_non_nested <- seen_data %>%
  select(mean_vis) %>%
  mutate(intercept = 1)

x_seen_year_species <- seen_data %>%
  mutate(year_classcode = paste(classcode, year, sep = '-')) %>%
  select(year_classcode) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(year_classcode, marker, fill = 0) %>%
  select(-index)

x_seen_region_cluster <- seen_data %>%
  select(region, geographic_cluster) %>%
  mutate(region_cluster = paste(geographic_cluster, region, sep = '-')) %>%
  select(region_cluster) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(region_cluster, marker, fill = 0) %>%
  select(-index)

log_density <-  seen_data$log_density

seen_species_index <- seen_data$numeric_classcode

seen_year_species_index <- data_frame(classcode = colnames(x_seen_year_species)) %>%
  mutate(classcode = str_replace_all(classcode,'(\\d)|-','')) %>%
  left_join(numeric_species_key, by = 'classcode') %>% {
    .$numeric_classcode
  }

seen_region_cluster_index <- data_frame(region = colnames(x_seen_region_cluster)) %>%
  mutate(region = str_replace_all(region,'(\\D)','') %>% as.factor() %>% as.numeric()) %>% {
    .$region
  }


# prepare seeing ----------------------------------------------------------

seeing_data <- data %>%
  select(any_seen,
         classcode,
         geographic_cluster,
         year,
         month,
         mean_vis,
         region,
         trunc_observer,
         cumulative_n_obs,
         method,
         level,
         surge) %>%
  na.omit() %>%
  left_join(numeric_species_key, by = 'classcode')

x_seeing_non_nested <- seeing_data %>%
  select(mean_vis) %>%
  mutate(intercept = 1)

x_seeing_year_species <- seeing_data %>%
  mutate(year_classcode = paste(classcode, year, sep = '-')) %>%
  select(year_classcode) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(year_classcode, marker, fill = 0) %>%
  select(-index)

x_seeing_region_cluster <- seeing_data %>%
  select(region, geographic_cluster) %>%
  mutate(region_cluster = paste(geographic_cluster, region, sep = '-')) %>%
  select(region_cluster) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(region_cluster, marker, fill = 0) %>%
  select(-index)

any_seen <-  seeing_data$any_seen %>% as.numeric()

seeing_species_index <- seeing_data$numeric_classcode

seeing_year_species_index <- data_frame(classcode = colnames(x_seeing_year_species)) %>%
  mutate(classcode = str_replace_all(classcode,'(\\d)|-','')) %>%
  left_join(numeric_species_key, by = 'classcode') %>% {
    .$numeric_classcode
  }

seeing_region_cluster_index <- data_frame(region = colnames(x_seeing_region_cluster)) %>%
  mutate(region = str_replace_all(region,'(\\D)','') %>% as.factor() %>% as.numeric()) %>% {
    .$region
  }

# fit TMB -----------------------------------------------------------------

n_species <- length(unique(seen_species_index))

n_clusters <- length(unique(seen_region_cluster_index))

ahnold_data <- list(
  x_seen_non_nested = x_seen_non_nested,
  x_seen_year_species = x_seen_year_species,
  x_seen_region_cluster = x_seen_region_cluster,
  log_density = log_density,
  seen_species_index = seen_species_index,
  seen_year_species_index = seen_year_species_index,
  seen_region_cluster_index = seen_region_cluster_index,
  x_seeing_non_nested = x_seeing_non_nested,
  x_seeing_year_species = x_seeing_year_species,
  x_seeing_region_cluster = x_seeing_region_cluster,
  any_seen = any_seen,
  seeing_species_index = seeing_species_index,
  seeing_year_species_index = seeing_year_species_index,
  seeing_region_cluster_index = seeing_region_cluster_index
)

ahnold_data <- map_if(ahnold_data, is.data.frame,~as.matrix(.x))

ahnold_params <- list(
  seen_non_nested_betas = rep(0,ncol(x_seen_non_nested)),
  seen_year_species_betas = rep(0, ncol(x_seen_year_species)),
  seen_region_cluster_betas = rep(0, ncol(x_seen_region_cluster)),
  seen_year_species_sigmas = rep(log(1), n_species),
  seen_region_cluster_sigmas = rep(log(1), n_clusters),
  seen_density_species_sigma = rep(log(1), n_species),
  seeing_non_nested_betas = rep(0,ncol(x_seeing_non_nested)),
  seeing_year_species_betas = rep(0, ncol(x_seeing_year_species)),
  seeing_region_cluster_betas = rep(0, ncol(x_seeing_region_cluster)),
  seeing_year_species_sigmas = rep(log(1), n_species),
  seeing_region_cluster_sigmas = rep(log(1), n_clusters),
  seeing_density_species_sigma = rep(log(1), n_species)
)

compile(here::here('scripts','fit_ahnold.cpp'),"-O0") # what is the -O0?

dyn.load(dynlib(here::here('scripts','fit_ahnold')))

ahnold_model <- MakeADFun(ahnold_data, ahnold_params, DLL="fit_ahnold")

ahnold_fit <- nlminb(ahnold_model$par, ahnold_model$fn, ahnold_model$gr,control = list(iter.max=1000, eval.max = 5000))

ahnold_fit <- nlminb(jitter(ahnold_fit$par), ahnold_model$fn, ahnold_model$gr,control = list(iter.max=1000, eval.max = 5000))

ahnold_fit <- nlminb(jitter(ahnold_fit$par), ahnold_model$fn, ahnold_model$gr,control = list(iter.max=1000, eval.max = 5000))


a <- ahnold_model$report()$log_density_hat

seen_data$log_density_hat <- a %>% as.numeric()

seen_data %>%
  ggplot(aes(log_density,log_density_hat, color = geographic_cluster %>% factor())) +
  geom_point()


seen_data %>%
  ggplot(aes(log_density_hat, log_density_hat - log_density, color = geographic_cluster %>% factor())) +
  geom_point()

a <- ahnold_model$report()$prob_seeing

seeing_data$prob_seen <- a %>% as.numeric()

seeing_data %>%
  ggplot(aes(prob_seen, any_seen, color = geographic_cluster %>% factor())) +
  geom_point()


