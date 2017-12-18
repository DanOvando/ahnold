set.seed(666)
library(TMB)
library(stringr)
library(purrr)
library(tidyverse)
demons::load_functions()

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


seen_data <- data %>%
  filter(any_seen == T) %>%
  left_join(numeric_species_key, by = 'classcode')

seeing_data <- data %>%
  left_join(numeric_species_key, by = 'classcode')


x_seen_non_nested <- seen_data %>%
  select(non_nested_variables) %>%
  mutate(intercept = 1) %>%
  spread_factors(drop_one = T) %>%
  purrrlyr::dmap(center_scale) %>%
  mutate(cumulative_n_obs_2 = cumulative_n_obs^2)

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


x_seeing_non_nested <- seeing_data %>%
  select(non_nested_variables) %>%
  mutate(intercept = 1) %>%
  spread_factors(drop_one = T) %>%
  purrrlyr::dmap(center_scale) %>%
  mutate(cumulative_n_obs_2 = cumulative_n_obs^2)

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


# create standard matrix --------------------------------------------------


standard_non_nested <- x_seeing_non_nested %>%
  mutate(classcode = seeing_data$classcode) %>%
  gather(variable, value, -classcode) %>%
  group_by(classcode, variable) %>%
  summarise(mean_value = mean(value)) %>%
  spread(variable, mean_value)

standard_non_nested <- map_df(seq_along(unique(seeing_data$year)), function(x,sdata){sdata}, sdata = standard_non_nested) %>%
  arrange(classcode) %>%
  ungroup() %>%
  select(-classcode)

standard_years <- expand.grid(classcode = unique(data$classcode), year = unique(data$year), stringsAsFactors = F) %>%
  mutate(year_classcode = paste(classcode, year, sep = '-')) %>%
  arrange(classcode) %>%
  select(year_classcode) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(year_classcode, marker, fill = 0) %>%
  select(-index)

standard_regions <- x_seeing_region_cluster %>%
  colMeans() %>%
  as.matrix() %>%
  t() %>%
  as_data_frame()

standard_regions <- map_df(1:nrow(standard_years), function(x,sdata){sdata}, sdata = standard_regions)


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
  seeing_region_cluster_index = seeing_region_cluster_index,
  standard_years = standard_years,
  standard_non_nested = standard_non_nested,
  standard_regions = standard_regions
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
  seeing_region_cluster_sigmas = rep(log(1), n_clusters)
)

compile(here::here('scripts','fit_ahnold.cpp'),"-O0") # what is the -O0?

dyn.load(dynlib(here::here('scripts','fit_ahnold')))

ahnold_model <- MakeADFun(ahnold_data, ahnold_params, DLL="fit_ahnold", random = c('seen_year_species_betas','seen_region_cluster_betas','seeing_year_species_betas','seeing_region_cluster_betas'))


ahnold_fit <- nlminb(ahnold_model$par, ahnold_model$fn, ahnold_model$gr,control = list(iter.max=1000, eval.max = 5000))

# ahnold_fit <- nlminb(jitter(ahnold_fit$par), ahnold_model$fn, ahnold_model$gr,control = list(iter.max=1000, eval.max = 5000))
#
# ahnold_fit <- nlminb(jitter(ahnold_fit$par), ahnold_model$fn, ahnold_model$gr,control = list(iter.max=1000, eval.max = 5000))
#

ahnold_model$gr(ahnold_fit$par)

rep <- ahnold_model$report()
sd_report <- sdreport(ahnold_model)
save(file = paste0(run_dir, '/ahnold-tmb-model.Rdata'), ahnold_model)
save(file = paste0(run_dir, '/ahnold-tmb-fit.Rdata'), ahnold_fit)
save(file = paste0(run_dir, '/ahnold-tmb-report.Rdata'), sd_report)

SD = summary(sd_report) # per millar example this should have the actual damn parameter names



ahnold_fit$par

ahnold_estimates <- data_frame(variable = names(ahnold_fit$par),
                          value = ahnold_fit$par)

year_terms <- ahnold_estimates %>%
  filter(str_detect())


a <- ahnold_model$report()$log_density_hat


ahnold_betas <- data_frame(beta = ahnold_fit$par, variable = names(ahnold_fit$par))

seen_abundance_trends <- ahnold_betas %>%
  filter(variable == 'seen_year_species_betas')

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
  ggplot(aes(any_seen, prob_seen)) +
  geom_boxplot()

seeing_data %>%
  group_by(any_seen) %>%
  summarise(a =mean(prob_seen))


