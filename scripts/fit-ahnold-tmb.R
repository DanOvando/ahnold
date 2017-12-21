set.seed(666)
library(TMB)
library(stringr)
library(purrr)
library(rstan)
library(tidyverse)
demons::load_functions()

rstan_options(auto_write = TRUE)

run_tmb <-  TRUE

tmb_to_stan <-  FALSE # fit the model in stan instead of TMB

run_name <- 'Working'

run_dir <- file.path('results', run_name)

load(file = paste0(run_dir, '/abundance_indices.Rdata'))

subspecies <- abundance_indices %>%
  select(classcode, targeted) %>%
  mutate(targeted = as.numeric(targeted > 0)) %>%
  mutate(a = 1:nrow(.)) %>%
  group_by(targeted) %>%
  top_n(2, a)

data <- abundance_indices %>%
  filter(
    population_structure == 'one-pop',
    population_filtering == 'all',
    data_source == 'length_to_density'
  ) %>%
  select(classcode, data) %>%
  unnest() %>%
  # filter(classcode %in% subspecies$classcode) %>%
  mutate(targeted = as.numeric(targeted > 0)) %>%
  group_by(classcode) %>%
  mutate(mean_vis = (mean_vis - mean(mean_vis)) / (2 * sd(mean_vis))) %>%
  mutate(mean_canopy = (mean_canopy - mean(mean_canopy, na.rm = T)) / (2 * sd(mean_canopy, na.rm = T))) %>%
  ungroup()

numeric_species_key <-
  data_frame(classcode = unique(data$classcode)) %>%
  arrange(classcode) %>%
  mutate(numeric_classcode = 1:nrow(.))


# prepare seen ------------------------------------------------------------

raw_length_covars <-
  c(
    'classcode',
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
    'cumulative_n_obs_2'
  )

non_nested_variables <- c(# 'zone',
  # 'site_side',
  'level',
  'mean_vis',
  'surge',
  'factor_month',
  'trunc_observer',
  'cumulative_n_obs')


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
  mutate(cumulative_n_obs_2 = cumulative_n_obs ^ 2)

x_seen_year_species <- seen_data %>%
  mutate(year_classcode = paste(classcode, year, sep = '-')) %>%
  select(year_classcode) %>%
  spread_factors(drop_one = F)

x_seen_region_cluster <- seen_data %>%
  select(region, geographic_cluster) %>%
  mutate(region_cluster = paste(geographic_cluster, region, sep = '-')) %>%
  select(region_cluster) %>%
  spread_factors(drop_one = T)


log_density <-  seen_data$log_density

seen_species_index <- seen_data$numeric_classcode

seen_year_species_index <-
  data_frame(classcode = colnames(x_seen_year_species)) %>%
  mutate(classcode = str_extract_all(classcode, '(?<=-).*(?=-)', '')) %>%
  left_join(numeric_species_key, by = 'classcode') %>% {
    .$numeric_classcode
  }

seen_region_cluster_index <-
  data_frame(region = colnames(x_seen_region_cluster)) %>%
  mutate(region = str_replace_all(region, '(\\D)', '') %>% as.factor() %>% as.numeric()) %>% {
    .$region
  }


# prepare seeing ----------------------------------------------------------


x_seeing_non_nested <- seeing_data %>%
  select(non_nested_variables) %>%
  mutate(intercept = 1) %>%
  spread_factors(drop_one = T) %>%
  purrrlyr::dmap(center_scale) %>%
  mutate(cumulative_n_obs_2 = cumulative_n_obs ^ 2)

x_seeing_year_species <- seeing_data %>%
  mutate(year_classcode = paste(classcode, year, sep = '-')) %>%
  select(year_classcode) %>%
  spread_factors(drop_one = F)

x_seeing_region_cluster <- seeing_data %>%
  select(region, geographic_cluster) %>%
  mutate(region_cluster = paste(geographic_cluster, region, sep = '-')) %>%
  select(region_cluster) %>%
  spread_factors(drop_one = T)

any_seen <-  seeing_data$any_seen %>% as.numeric()

seeing_species_index <- seeing_data$numeric_classcode

seeing_year_species_index <-
  data_frame(classcode = colnames(x_seeing_year_species)) %>%
  mutate(classcode = str_extract(classcode, '(?<=-).*(?=-)')) %>%
  left_join(numeric_species_key, by = 'classcode') %>% {
    .$numeric_classcode
  }

seeing_region_cluster_index <-
  data_frame(region = colnames(x_seeing_region_cluster)) %>%
  mutate(region = str_replace_all(region, '(\\D)', '') %>% as.factor() %>% as.numeric()) %>% {
    .$region
  }


# create standard matrix --------------------------------------------------


standard_non_nested <- x_seeing_non_nested %>%
  mutate(classcode = seeing_data$classcode) %>%
  gather(variable, value, -classcode) %>%
  group_by(classcode, variable) %>%
  summarise(mean_value = mean(value)) %>%
  spread(variable, mean_value)

standard_non_nested <-
  map_df(seq_along(unique(seeing_data$year)), function(x, sdata) {
    sdata
  }, sdata = standard_non_nested) %>%
  arrange(classcode) %>%
  ungroup() %>%
  select(-classcode) %>%
  select(colnames(x_seeing_non_nested)) # make sure things are in the same column order as they used to be

standard_years <-
  expand.grid(
    classcode = unique(data$classcode),
    year = unique(data$year),
    stringsAsFactors = F
  ) %>%
  mutate(year_classcode = paste(classcode, year, sep = '-')) %>%
  arrange(classcode) %>%
  select(year_classcode) %>%
  spread_factors(drop_one = F)


standard_regions <- x_seeing_region_cluster %>%
  colMeans() %>%
  as.matrix() %>%
  t() %>%
  as_data_frame()

standard_regions <-
  map_df(1:nrow(standard_years), function(x, sdata) {
    sdata
  }, sdata = standard_regions)


# prepare difference in difference ----------------------------------------


load(file = here::here(run_dir, 'did_models.Rdata'))

did_data <- did_models %>%
  filter(
    timing == 'years',
    complexity == 'kitchen_sink',
    population_structure == 'one-pop',
    abundance_source == 'glm_abundance_index',
    data_source == 'length_to_density',
    population_filtering == 'all'
  ) %>% {
    .$data[[1]]
  }

x_did_non_nested <- did_data %>%
  select(targeted, factor_year, loo) %>%
  spread_factors(drop_one = T)

x_did_did_term <- did_data %>%
  select(targeted, factor_year) %>%
  mutate(index = 1:nrow(.)) %>%
  spread(factor_year, targeted, fill = 0) %>%
  arrange(index) %>%
  select(-(1:2)) #drop index and one factor level

x_did_non_nested <- bind_cols(x_did_non_nested, x_did_did_term) %>%
  mutate(intercept = 1)


species_effect_variables <-
  c('mean_enso', 'mean_annual_kelp', 'temp_deviation')

x_did_species_effects <-
  as_data_frame(matrix(NA, nrow = nrow(did_data), ncol = 0))

for (i in 1:length(species_effect_variables)) {
  temp <- did_data %>%
    select(species_effect_variables[i], classcode) %>%
    mutate(index = 1:nrow(.)) %>%
    spread(classcode, species_effect_variables[i], fill = 0) %>%
    arrange(index) %>%
    select(-index) %>%
    set_names(paste(colnames(.), species_effect_variables[i], sep = '-'))

  x_did_species_effects <-  x_did_species_effects %>%
    bind_cols(temp)

}

x_did_species_effects_index <- colnames(x_did_species_effects)

x_did_species_effects_index <-
  str_extract_all(x_did_species_effects_index, '.*(?=-)', simplify = T)[, 1]

x_did_species_effects_index <-
  as.numeric(as.factor(x_did_species_effects_index))

x_did <- bind_cols(x_did_non_nested, x_did_species_effects)

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
  standard_regions = standard_regions,
  x_did_non_nested = x_did_non_nested,
  x_did_species_effects = x_did_species_effects,
  x_did_species_effects_index = x_did_species_effects_index
)

ahnold_data <- map_if(ahnold_data, is.data.frame,  ~ as.matrix(.x))

any_na <- map_lgl(ahnold_data,  ~ any(is.na(.x))) %>% any()

if (any_na) {
  stop("NAs in ahnold_data")
}

ahnold_params <- list(
  seen_non_nested_betas = rep(0, ncol(x_seen_non_nested)),
  seen_year_species_betas = rep(0, ncol(x_seen_year_species)),
  seen_region_cluster_betas = rep(0, ncol(x_seen_region_cluster)),
  seen_year_species_sigmas = rep(log(1), n_species),
  seen_density_species_sigma = rep(log(1), n_species),
  seeing_non_nested_betas = rep(0, ncol(x_seeing_non_nested)),
  seeing_year_species_betas = rep(0, ncol(x_seeing_year_species)),
  seeing_region_cluster_betas = rep(0, ncol(x_seeing_region_cluster)),
  seeing_year_species_sigmas = rep(log(1), n_species),
   did_non_nested_betas = rep(0, ncol(x_did_non_nested)),
  # did_species_betas = rep(0, ncol(x_did_species_effects)) ,
   did_sigma = log(1)
)

any_na <- map_lgl(ahnold_params,  ~ any(is.na(.x))) %>% any()

if (any_na) {
  stop("NAs in ahnold_data")
}



if (run_tmb == T) {
  script_name <-  'fit_ahnold_hurdle'
  compile(here::here('scripts', paste0(script_name, '.cpp')), "-O0") # what is the -O0?

  dyn.load(dynlib(here::here('scripts', script_name)))

  ahnold_model <-
    MakeADFun(
      ahnold_data,
      ahnold_params,
      DLL = script_name,
      random = c('seen_year_species_betas',
                 'seeing_year_species_betas')
    )

  if (tmb_to_stan == F) {
    set.seed(42)
    ahnold_fit <-
      nlminb(
        ahnold_model$par,
        ahnold_model$fn,
        ahnold_model$gr,
        control = list(iter.max = 4000, eval.max = 5000)
      )


    save(file = here::here(run_dir, 'ahnold-tmb-model.Rdata'),
         ahnold_model)

    save(file = here::here(run_dir, 'ahnold-tmb-fit.Rdata'), ahnold_fit)
    print('made it here')
    sd_report <- sdreport(ahnold_model)

    # sd_report <- sdreport(ahnold_model,getReportCovariance = TRUE, skip.delta.method = TRUE)

    save(file = here::here(run_dir, 'ahnold-tmb-report.Rdata'),
         sd_report)
  } else{
    ahnold_fit <- tmbstan::tmbstan(ahnold_model, chains = 1)

    save(file = here::here(run_dir, 'ahnold-tmbtostan-model.Rdata'),
         ahnold_model)

    save(file = here::here(run_dir, 'ahnold-tmbtostan-fit.Rdata'),
         ahnold_fit)


  }

} else{
  load(file = here::here(run_dir, 'ahnold-tmb-model.Rdata'))

  load(file = here::here(run_dir, 'ahnold-tmb-fit.Rdata'))

  load(file = here::here(run_dir, 'ahnold-tmb-report.Rdata'))

}


# pull out estimates ------------------------------------------------------


ahnold_estimates <-
  summary(sd_report) %>%
  as.data.frame() %>%
  mutate(variable = rownames(.)) %>%
  set_names(tolower) %>%
  rename(std_error = `std. error`) %>%
  mutate(lower = estimate - 1.96 * std_error,
         upper = estimate + 1.96 * std_error)

seen_non_nested_betas <- ahnold_estimates %>%
  filter(str_detect(variable, 'seen_non_nested_betas')) %>%
  mutate(group = variable) %>%
  mutate(variable = colnames(x_seen_non_nested))

seeing_non_nested_betas <- ahnold_estimates %>%
  filter(str_detect(variable, 'seeing_non_nested_betas')) %>%
  mutate(group = variable) %>%
  mutate(variable = colnames(x_seeing_non_nested))

seen_region_betas <- ahnold_estimates %>%
  filter(str_detect(variable, 'seen_region_cluster_betas')) %>%
  mutate(group = variable) %>%
  mutate(variable = colnames(x_seen_region_cluster))

seeing_region_betas <- ahnold_estimates %>%
  filter(str_detect(variable, 'seeing_region_cluster_betas')) %>%
  mutate(group = variable) %>%
  mutate(variable = colnames(x_seeing_region_cluster))

seen_year_species_betas <- ahnold_estimates %>%
  filter(str_detect(variable, 'seen_year_species_betas')) %>%
  mutate(group = variable) %>%
  mutate(variable = colnames(x_seen_year_species))

seeing_year_species_betas <- ahnold_estimates %>%
  filter(str_detect(variable, 'seeing_year_species_betas')) %>%
  mutate(group = variable) %>%
  mutate(variable = colnames(x_seeing_year_species))

betas <- bind_rows(
  seen_non_nested_betas,
  seeing_non_nested_betas,
  seen_region_betas,
  seeing_region_betas,
  seen_year_species_betas,
  seeing_year_species_betas
) %>%
  as_data_frame()



betas %>%
  ggplot() +
  geom_pointrange(aes(
    variable,
    y = estimate,
    ymin = lower,
    ymax = upper
  )) +
  geom_hline(aes(yintercept = 0)) +
  coord_flip() +
  facet_wrap(~ group, scales = 'free')


abundance_indices <- ahnold_estimates %>%
  filter(str_detect(variable, 'standardized_abundance')) %>%
  mutate(classcode_year = colnames(x_seen_year_species)) %>%
  mutate(classcode_year = str_replace(classcode_year, 'year_classcode-','')) %>%
  separate(classcode_year, c('classcode', 'year'), '-') %>%
  mutate(year = as.numeric(year))


abundance_indices %>%
  ggplot() +
  geom_line(aes(year, estimate, color = classcode), show.legend = F) +
  geom_ribbon(
    aes(
      x = year,
      ymin = lower,
      ymax = upper,
      fill = classcode
    ),
    alpha = 0.25,
    show.legend = F
  ) +
  facet_wrap( ~ classcode, scales = 'free_y')


# diagnostics -------------------------------------------------------------


a <- ahnold_model$report()$log_density_hat


ahnold_betas <-
  data_frame(beta = ahnold_fit$par, variable = names(ahnold_fit$par))

seen_abundance_trends <- ahnold_betas %>%
  filter(variable == 'seen_year_species_betas')

seen_data$log_density_hat <- a %>% as.numeric()

seen_data %>%
  ggplot(aes(log_density, log_density_hat, color = geographic_cluster %>% factor())) +
  geom_point()


lm(log_density ~ log_density_hat, data = seen_data) %>% summary()

seen_data %>%
  ggplot(
    aes(
      log_density_hat,
      log_density_hat - log_density,
      color = geographic_cluster %>% factor()
    )
  ) +
  geom_point()

a <- ahnold_model$report()$prob_seeing

seeing_data$prob_seen <- a %>% as.numeric()

seeing_data %>%
  ggplot(aes(any_seen, prob_seen)) +
  geom_boxplot()

seeing_data %>%
  group_by(any_seen) %>%
  summarise(a = mean(prob_seen))
