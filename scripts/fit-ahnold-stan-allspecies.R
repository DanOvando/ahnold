rm(list = ls())
set.seed(123)
library(tidyverse)
library(rstan)
library(stringr)
library(purrr)
rstan_options(auto_write = TRUE)

run_name <- 'Working'

run_dir <- file.path('results', run_name)

load(file = paste0(run_dir, '/abundance_indices.Rdata'))

subspecies <- abundance_indices %>%
  select(classcode, targeted) %>%
  mutate(a = 1:nrow(.)) %>%
  group_by(targeted) %>%
  top_n(3,a)

data <- abundance_indices %>%
  filter(population_structure == 'one-pop',
         population_filtering == 'all',
         data_source == 'length_to_density') %>%
  select(classcode, data) %>%
  unnest() #%>%
  # filter(classcode %in% subspecies$classcode)

wtf <- data %>%
  filter(any_seen == T)

# View(data$data[[1]] %>% filter(any_seen == T))

seen_data <- data %>%
  filter(any_seen == T) %>%
  mutate(mean_vis = (mean_vis - mean(mean_vis)) / (2 * sd(mean_vis))) %>%
  mutate(mean_canopy = (mean_canopy - mean(mean_canopy, na.rm = T)) / (2 * sd(mean_canopy, na.rm = T))) %>%
  mutate(classcode = as.factor(classcode)) %>%
  select(log_density,
         classcode,
         year,
         month,
         mean_vis,
         region,
         trunc_observer,
         cumulative_n_obs,
         method,
         level,
         surge) %>%
  na.omit()

arm_data <- seen_data %>%
  mutate(factor_year = as.factor(year),
         factor_month = as.factor(month))

year_effects <- test$stan_summary %>%
  as.data.frame() %>%
  mutate(variable = rownames(.)) %>%
  as_data_frame() %>%
  filter(str_detect(variable,'factor_year')) %>%
  mutate(year = str_replace_all(variable,'\\D','') %>% as.numeric(),
         species = str_replace_all(variable,'.*(?=classcode)','')) %>%
  mutate(species = str_replace_all(species, 'classcode',''))

year_effects %>%
  ggplot() +
  geom_line(aes(year, exp(mean + se_mean^2/2), color = species)) +
  facet_wrap(~species, scales = 'free_y') +
  theme_classic()

save(file = 'rstanarm-test.Rdata', test)


abundance_comparison <-did_models %>%
  filter(data_source == 'length_to_density',
         population_structure == 'one-pop',
         abundance_source == 'glm_abundance_index' |abundance_source == 'raw_abundance_index' ,
         population_filtering == 'all',
         did_term_names == 'years-protected') %>%
  select(abundance_source, data) %>%
  unnest()


abundance_comparison %>%
  group_by(commonname, abundance_source) %>%
  mutate(abundance_index = abundance_index / max(abundance_index)) %>%
  ungroup() %>%
  ggplot(aes(year, abundance_index, color = abundance_source)) +
  geom_line() +
  facet_wrap(~classcode) +
  theme_classic() +
  theme(strip.text = element_text(size = 8))
# test <- rstanarm::stan_glm(
#   'log_density ~
#   factor_year:classcode  + mean_vis + factor_month + trunc_observer + cumulative_n_obs + method + level + surge',
#   data = arm_data,
#   refresh = 1,
#   chains = 1
# )


# lme4::glmer('log_density ~ (factor_year|classcode) + mean_canopy + mean_vis', data = seen_data)

year_species_data <- seen_data %>%
  # filter(classcode %in% c('hcar','hros')) %>%
  mutate(year_classcode = paste(classcode, year, sep = '-')) %>%
  select(year_classcode) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(year_classcode, marker, fill = 0) %>%
  select(-index)

observer_data <- seen_data %>%
  select(trunc_observer) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(trunc_observer, marker, fill = 0) %>%
  select(-index)

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

regions_per_species <- seen_data %>%
  group_by(classcode) %>%
  summarise(nr = length(unique(region))) %>% {
    .$nr
  }

visibility_data <- seen_data %>%
  select(mean_vis)

# canopy_data <- seen_data %>%
#   select(mean_canopy)

n_intercepts <- length(years_per_species)

# intercept_data <-
#   matrix(1, nrow = nrow(year_species_data), ncol = n_intercepts) %>%
#   as_data_frame() %>%
#   set_names(paste0(unique(seen_data$classcode)[1:n_intercepts], '_intercept'))

intercept_data <- seen_data %>%
  select(classcode) %>%
  mutate(seen = 1, id = 1:nrow(.)) %>%
  spread(classcode, seen, fill = 0) %>%
  select(-id) %>%
  set_names(glue::glue('{colnames(.)}_intercept'))

x_seen <-
  bind_cols(
    intercept_data,
    visibility_data,
    year_species_data,
    region_species_data
  )

# seeing data -------
#
seeing_data <- data %>%
  mutate(mean_vis = (mean_vis - mean(mean_vis)) / (2 * sd(mean_vis))) %>%
  mutate(mean_canopy = (mean_canopy - mean(mean_canopy, na.rm = T)) / (2 * sd(mean_canopy, na.rm = T))) %>%
  mutate(classcode = as.factor(classcode)) %>%
  select(any_seen,
         classcode,
         year,
         mean_canopy,
         mean_vis,
         region) %>%
  na.omit()

# lme4::glmer('log_density ~ (factor_year|classcode) + mean_canopy + mean_vis', data = seen_data)
year_species_data <- seeing_data %>%
  # filter(classcode %in% c('hcar','hros')) %>%
  mutate(year_classcode = paste(classcode, year, sep = '-')) %>%
  select(year_classcode) %>%
  mutate(marker = 1,
         index = 1:nrow(.)) %>%
  spread(year_classcode, marker, fill = 0) %>%
  select(-index)

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

regions_per_species <- seeing_data %>%
  group_by(classcode) %>%
  summarise(nr = length(unique(region))) %>% {
    .$nr
  }

visibility_data <- seeing_data %>%
  select(mean_vis)

# canopy_data <- seeing_data %>%
#   select(mean_canopy)

n_intercepts <- length(years_per_species)
#
# intercept_data <-
#   matrix(1, nrow = nrow(year_species_data), ncol = n_intercepts) %>%
#   as_data_frame() %>%
#   set_names(paste0(unique(seeing_data$classcode)[1:n_intercepts], '_intercept'))

intercept_data <- seeing_data %>%
  select(classcode) %>%
  mutate(seen = 1, id = 1:nrow(.)) %>%
  spread(classcode, seen, fill = 0) %>%
  select(-id) %>%
  set_names(glue::glue('{colnames(.)}_intercept'))

x_seeing <-
  bind_cols(
    intercept_data,
    visibility_data,
    year_species_data,
    region_species_data
  )


# check <- rstanarm::stan_glmer('log_density ~ (1|factor_year) + (1|region) + mean_canopy + mean_vis', data = seen_data)
# check_binom <- rstanarm::stan_glmer('any_seen ~ (1|factor_year) + (1|region) + mean_canopy + mean_vis', data = seeing_data, family = 'binomial')

create_standard_mat <- function(classcode, species_data, x) {
  environment <- species_data %>%
    summarise(mean_vis = mean(mean_vis),
              mean_canopy = mean(mean_canopy))

  n_years <- length(unique(species_data$year))

  common_region <- species_data %>%
    group_by(region) %>%
    count() %>%
    arrange(desc(n)) %>% {
      .$region[1]
    }

  sframe <- matrix(0, nrow = n_years, ncol = ncol(x)) %>%
    as_data_frame() %>%
    set_names(colnames(x))

  sframe$mean_vis <- environment$mean_vis

  sframe$mean_canopy <- environment$mean_canopy

  sframe[, paste0(classcode, '_intercept')] <- 1

  sframe[, paste0(classcode, '-', common_region)] <- 1

  years <- unique(species_data$year) %>% sort()

  for (y in 1:length(years)) {
    sframe[y, paste0(classcode, '-', years[y])] <- 1

  }

  return(sframe)

}


species_data <- seeing_data %>%
  nest(-classcode)

standard_matrix <- species_data %>%
  mutate(smat = map2(classcode, data, create_standard_mat, x = x_seeing)) %>%
  select(smat) %>%
  unnest() %>%
  select(-mean_canopy)

# a <- standard_matrix %>%
#   select(contains('intercept')) %>%
#   gather(variable, value)
#

# prep did data -----------------------------------------------------------

load(file = paste0(run_dir, '/did_models.Rdata'))

dat <- did_models %>%
  filter(
    population_filtering == 'all',
    population_structure == 'one-pop',
    data_source == 'length_to_density'
  )

x_did <- dat$data[[1]] %>%
  filter(classcode %in% unique(seen_data$classcode)) %>%
  select(classcode,
         year,
         targeted,
         contains('enso'),
         contains('pdo'),
         mean_annual_temp) %>%
  mutate(fished = targeted) %>%
  mutate(post_mpa = as.numeric(year > 2003)) %>%
  spread(year, targeted, fill = 0) %>%
  mutate(intercept = 1)


x_did <- x_did %>%
  select(-classcode) #,-contains('enso'),-contains('pdo'),-contains('temp'))

#
# x_seen <- x_seen %>%
#   select(-contains('intercept')) %>%
#   mutate(intercept = 1,
#          species_intercept = 1)

# x_seeing <-  x_seeing %>%
#   select(-contains('intercept')) %>%
#   mutate(intercept = 1,
#          species_intercept = 1)

# standard_matrix <- standard_matrix %>%
#   select(-contains('intercept')) %>%
#   mutate(intercept = 1,
#          species_intercept = 1)


did_names <- colnames(x_did)

enso_locations <- which(str_detect(did_names, 'enso'))

pdo_locations <- which(str_detect(did_names, 'pdo'))

did_positions <- which(str_detect(did_names, '20'))

stan_data <- list(
  n_parameters = ncol(x_seen),
  n_did_parameters = ncol(x_did),
  n_observations_seen = nrow(x_seen),
  n_observations_seeing = nrow(x_seeing),
  n_observations_did = nrow(x_did),
  n_species = length(years_per_species),
  n_year_species = ncol(year_species_data),
  n_region_species = ncol(region_species_data),
  n_standard = nrow(standard_matrix),
  years_per_species = years_per_species,
  regions_per_species = regions_per_species,
  year_species_positions = which(str_detect(colnames(x_seen), '\\d')),
  region_species_positions = which(str_detect(
    colnames(x_seen), paste(unique(data$region), collapse = '|')
  )),
  # species_intercepts_positions = which(str_detect(colnames(x_seen), 'intercept')),
  # species_intercept_position = which(colnames(x_seen) == 'species_intercept'),
  non_nested_positions = which(!str_detect(colnames(x_seen), '\\d') & !str_detect(
    colnames(x_seen), paste(unique(data$region), collapse = '|')
  )),
  n_non_nested = length(which(!str_detect(colnames(x_seen), '\\d') & !str_detect(
    colnames(x_seen), paste(unique(data$region), collapse = '|')
  ))),
  did_positions = did_positions,
  n_did = length(did_positions),
  x_seen = x_seen,
  x_seeing = x_seeing,
  x_did = x_did,
  log_density = seen_data$log_density,
  observed = seeing_data$any_seen %>% as.numeric(),
  standard_matrix = standard_matrix,
  cauchy_2 = 2.5
)


# a <- Sys.time()
# ahnold_stan_fit <- stan(
#   file = 'scripts/fit-ahnold-allspecies-joint.stan',
#   data = stan_data,
#   chains = 1,
#   warmup = 100,
#   iter = 200,
#   cores = 1,
#   refresh = 1
# )
#
# Sys.time() - a



a <- Sys.time()
ahnold_stan_fit <- stan(
  file = 'scripts/fit-ahnold-abundance.stan',
  data = stan_data,
  chains = 1,
  warmup = 100,
  iter = 200,
  cores = 1,
  refresh = 1
)

Sys.time() - a

save(file = 'woop.Rdata', ahnold_stan_fit)

# abundance_names <- colnames(standard_matrix)[ str_detect(colnames(standard_matrix), '\\d')]
#
#
# abundance_indicies <- ahnold_stan_fit %>%
#   as.data.frame() %>%
#   as_data_frame() %>%
#   select(contains('standardized_abundance_index')) %>%
#   set_names(abundance_names) %>%
#   gather(name, abundance) %>%
#   mutate(species = str_replace_all(name,'\\d|-',''),
#          year = str_replace_all(name,'\\D|-','') %>% as.numeric())
#
# abundance_indicies %>%
#   ggplot(aes(year, abundance, color = species)) +
#   geom_point() +
#   geom_smooth() +
#   facet_wrap(~species)
#
# abundance_indicies %>%
#   group_by(species) %>%
#   mutate(abundance = abundance / max(abundance)) %>%
#   ggplot(aes(factor(year), abundance, fill = species)) +
# geom_violin() +
#   facet_wrap(~species)


# check2 <-
#   rstanarm::stan_glmer(log_density ~ (1 |classcode) + mean_canopy + mean_vis,
#                        data = seen_data,chains = 1)

# control = list(adapt_delta = 0.8)
