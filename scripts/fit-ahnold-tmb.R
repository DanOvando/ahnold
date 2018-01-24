set.seed(666)
library(TMB)
library(stringr)
library(purrr)
library(rstan)
library(tidyverse)
demons::load_functions()

rstan_options(auto_write = TRUE)

run_tmb <- TRUE

tmb_to_stan <- FALSE # fit the model in stan instead of TMB

run_tmb <-  T

max_generations <- 4

run_name <- "Working"

script_name <- "fit_ahnold"

run_dir <- file.path("results", run_name)

load(file = paste0(run_dir, "/ahnold_model_data.Rdata"))

mpa_year <-  2003

data <- abundance_models %>%
  filter(
    population_structure == "one-pop",
    population_filtering == "all",
    data_source == "length_to_density"
  ) %>%
  select(classcode, data) %>%
  unnest() %>%
  ungroup() %>%
  mutate(
         mean_temp = ifelse(is.na(mean_temp), mean(mean_temp, na.rm = T), mean_temp)) %>%
  mutate(temp_deviation = (mean_temp - temperature)^2) %>%
  mutate(generations_protected = pmin(round((year - mpa_year - 1) / tm), max_generations),
  targeted = as.numeric(targeted > 0)) %>%
  mutate(temp_deviation = center_scale(temp_deviation))


# test <-
#   rstanarm::stan_glmer(
#     log_density ~ (factor_year  - 1|  classcode) + temp_deviation - 1,
#     data = data %>% filter(any_seen == T),
#     cores = 1,
#     chains = 1
#   )

numeric_species_key <-
  data_frame(classcode = unique(data$classcode)) %>%
  arrange(classcode) %>%
  mutate(numeric_classcode = 1:nrow(.))


# prepare seen ------------------------------------------------------------

non_nested_variables <- c(
  'targeted',
  'factor_year',
  'level',
  'factor_month',
  'cumulative_n_obs',
  'temp_deviation',
  'surge',
  'mean_canopy',
  'mean_depth'
)

seen_has_important <- data %>%
  filter(any_seen == T) %>%
  select(non_nested_variables) %>%
  mutate(index = 1:nrow(.)) %>%
  na.omit()

seeing_has_important <- data %>%
  select(non_nested_variables) %>%
  mutate(index = 1:nrow(.)) %>%
  na.omit()


numeric_species_key <-
  data_frame(classcode = unique(data$classcode)) %>%
  arrange(classcode) %>%
  mutate(numeric_classcode = 1:nrow(.))

seen_data <- data %>%
  filter(any_seen == T) %>%
  left_join(numeric_species_key, by = "classcode") %>%
  slice(seen_has_important$index)

log_density <- seen_data$log_density

seen_data <- seen_data %>%
  select(-log_density)

seeing_data <- data %>%
  left_join(numeric_species_key, by = "classcode") %>%
  slice(seeing_has_important$index)

any_seen <- seeing_data$any_seen

seeing_data <- seeing_data %>%
  select(-any_seen)


x_seen_non_nested <- seen_data %>%
  select(non_nested_variables) %>%
  mutate(intercept = 1) %>%
  spread_factors(drop_one = T) %>%
  purrrlyr::dmap(center_scale) %>%
  mutate(cumulative_n_obs_2 = cumulative_n_obs ^ 2) %>%
  mutate(targeted = seen_data$targeted)


x_seen_enso <- seen_data %>%
  select(classcode, enso) %>%
  mutate(index = 1:nrow(.)) %>%
  spread(classcode, enso, fill = 0) %>%
  arrange(index) %>%
  select(-index) %>%
  set_names(paste0("enso_",colnames(.)))

x_seen_pdo <- seen_data %>%
  select(classcode, pdo) %>%
  mutate(index = 1:nrow(.)) %>%
  spread(classcode, pdo, fill = 0) %>%
  arrange(index) %>%
  select(-index) %>%
  set_names(paste0("pdo_",colnames(.)))

x_seen_non_nested <- x_seen_non_nested %>%
  bind_cols(x_seen_enso,
            x_seen_pdo)

x_seen_did <- seen_data %>%
  select(targeted, factor_year) %>%
  mutate(index = 1:nrow(.)) %>%
  spread(factor_year, targeted, fill = 0) %>%
  arrange(index) %>%
  select(-(1:2)) # drop index and one factor level

x_seen_year_species <- seen_data %>%
  mutate(year_classcode = paste(classcode, year, sep = "-")) %>%
  select(year_classcode) %>%
  spread_factors(drop_one = F)

seen_species_index <- seen_data$numeric_classcode

seen_year_species_index <-
  data_frame(classcode = colnames(x_seen_year_species)) %>%
  mutate(classcode = str_extract_all(classcode, "(?<=-).*(?=-)", "")) %>%
  left_join(numeric_species_key, by = "classcode") %>%
  {
    .$numeric_classcode
  }

x_seen_region_cluster <- seen_data %>%
  select(region, geographic_cluster) %>%
  mutate(region_cluster = paste(geographic_cluster, region, sep = "-")) %>%
  select(region_cluster) %>%
  spread_factors(drop_one = T)

seen_region_cluster_index <-
  data_frame(region = colnames(x_seen_region_cluster)) %>%
  mutate(region = str_replace_all(region, "(\\D)", "") %>% as.factor() %>% as.numeric()) %>%
  {
    .$region
  }

# x_seen <-
#   bind_cols(seen_data %>% select(log_density),
#             x_seen_non_nested, x_seen_did) #,
#             # x_seen_region_cluster)


# prepare seeing ----------------------------------------------------------

x_seeing_non_nested <- seeing_data %>%
  select(non_nested_variables) %>%
  mutate(intercept = 1) %>%
  spread_factors(drop_one = T) %>%
  purrrlyr::dmap(center_scale) %>%
  mutate(cumulative_n_obs_2 = cumulative_n_obs ^ 2)


x_seeing_enso <- seeing_data %>%
  select(classcode, enso) %>%
  mutate(index = 1:nrow(.)) %>%
  spread(classcode, enso, fill = 0) %>%
  arrange(index) %>%
  select(-index) %>%
  set_names(paste0("enso_",colnames(.)))

x_seeing_pdo <- seeing_data %>%
  select(classcode, pdo) %>%
  mutate(index = 1:nrow(.)) %>%
  spread(classcode, pdo, fill = 0) %>%
  arrange(index) %>%
  select(-index) %>%
  set_names(paste0("pdo_",colnames(.)))

x_seeing_non_nested <- x_seeing_non_nested %>%
  bind_cols(x_seeing_enso,
            x_seeing_pdo)

x_seeing_did <- seeing_data %>%
  select(targeted, factor_year) %>%
  mutate(index = 1:nrow(.)) %>%
  spread(factor_year, targeted, fill = 0) %>%
  arrange(index) %>%
  select(-(1:2)) # drop index and one factor level

x_seeing_year_species <- seeing_data %>%
  mutate(year_classcode = paste(classcode, year, sep = "-")) %>%
  select(year_classcode) %>%
  spread_factors(drop_one = F)

seeing_species_index <- seeing_data$numeric_classcode

year_species_index <-
  data_frame(classcode = colnames(x_seeing_year_species)) %>%
  mutate(classcode = str_extract_all(classcode, "(?<=-).*(?=-)", "")) %>%
  left_join(numeric_species_key, by = "classcode") %>%
  {
    .$numeric_classcode
  }


x_seeing_region_cluster <- seeing_data %>%
  select(region, geographic_cluster) %>%
  mutate(region_cluster = paste(geographic_cluster, region, sep = "-")) %>%
  select(region_cluster) %>%
  spread_factors(drop_one = T)

seeing_region_cluster_index <-
  data_frame(region = colnames(x_seeing_region_cluster)) %>%
  mutate(region = str_replace_all(region, "(\\D)", "") %>% as.factor() %>% as.numeric()) %>%
  {
    .$region
  }

# x_seeing <-
#   bind_cols(
#     seeing_data %>% select(any_seen),
#     x_seeing_non_nested,
#     x_seeing_did,
#     x_seeing_region_cluster
#   )

# create standard matrix --------------------------------------------------


standard_non_nested <- x_seeing_non_nested %>%
  gather(variable, value) %>%
  group_by(variable) %>%
  summarise(mean_value = mean(value, na.rm = T)) %>%
  spread(variable, mean_value)

standard_non_nested <-
  map_df(seq_along(unique(seeing_data$factor_year)), function(x, sdata) {
    sdata
  }, sdata = standard_non_nested) %>%
  ungroup() %>%
  select(colnames(x_seeing_non_nested))

standard_year_species <- x_seeing_year_species %>%
  gather(variable, value) %>%
  group_by(variable) %>%
  summarise(mean_value = mean(value, na.rm = T)) %>%
  spread(variable, mean_value)

standard_region_cluster <- x_seeing_region_cluster %>%
  gather(variable, value) %>%
  group_by(variable) %>%
  summarise(mean_value = mean(value, na.rm = T)) %>%
  spread(variable, mean_value)

standard_region_cluster <- standard_region_cluster[rep(1, n_distinct(seeing_data$factor_year)),]


standard_year_species <-
  map_df(seq_along(unique(seeing_data$factor_year)), function(x, sdata) {
    sdata
  }, sdata = standard_year_species) %>%
  ungroup() %>%
  select(colnames(x_seeing_year_species))


standard_did_with_mpa <- data_frame(year = unique(seeing_data$factor_year)) %>%
  spread_factors(drop_one = T)

standard_did_without_mpa <- standard_did_with_mpa

standard_did_without_mpa[standard_did_without_mpa > 0] = 0

# fit TMB -----------------------------------------------------------------


seen_species_index <- data %>%
  filter(any_seen == T) %>%
  select(non_nested_variables, classcode) %>%
  na.omit() %>%
  left_join(numeric_species_key, by = "classcode") %>%  {
    .$numeric_classcode
  }

n_species <- length(unique(seen_species_index))

ahnold_data <- list(
  x_seen_non_nested = x_seen_non_nested,
  x_seen_did = x_seen_did,
  x_seen_year_species = x_seen_year_species,
  x_seen_region_cluster = x_seen_region_cluster,
  x_seeing_non_nested = x_seeing_non_nested,
  x_seeing_did = x_seeing_did,
  x_seeing_year_species = x_seeing_year_species,
  x_seeing_region_cluster = x_seeing_region_cluster,
  region_cluster_index = seeing_region_cluster_index,
  log_density = log_density,
  any_seen = any_seen,
  standard_non_nested = standard_non_nested,
  standard_year_species = standard_year_species,
  standard_did_with_mpa = standard_did_with_mpa,
  standard_did_without_mpa = standard_did_without_mpa,
  standard_region_cluster = standard_region_cluster,
  seen_species_index = seen_species_index,
  year_species_index = year_species_index,
  seen_weights = rep(1, nrow(x_seen_non_nested)),
  seeing_weights = rep(1, nrow(x_seeing_non_nested))
)

ahnold_data <- map_if(ahnold_data, is.data.frame, ~ as.matrix(.x))

any_na <- map_lgl(ahnold_data, ~ any(is.na(.x))) %>% any()

if (any_na) {
  stop("NAs in ahnold_data")
}

ahnold_params <- list(
  seen_non_nested_betas = rep(0, ncol(x_seen_non_nested)),
  seen_did_betas = rep(0, ncol(x_seen_did)),
  seen_year_species_betas = rep(0, ncol(x_seen_year_species)),
  seen_year_species_sigmas = rep(log(1), n_species),
  seen_region_cluster_betas = rep(0, ncol(x_seen_region_cluster)),
  seen_region_cluster_sigmas = rep(log(1), n_distinct(seeing_region_cluster_index)),
  seen_density_species_sigma = rep(log(1), n_species),
  seeing_non_nested_betas = rep(0, ncol(x_seeing_non_nested)),
  seeing_did_betas = rep(0, ncol(x_seeing_did)),
  seeing_year_species_betas = rep(0, ncol(x_seen_year_species)),
  seeing_year_species_sigmas = rep(log(1), n_species),
  seeing_region_cluster_betas = rep(0, ncol(x_seeing_region_cluster)),
  seeing_region_cluster_sigmas = rep(log(1), n_distinct(seeing_region_cluster_index))
)

any_na <- map_lgl(ahnold_params, ~ any(is.na(.x))) %>% any()

if (any_na) {
  stop("NAs in ahnold_params")
}



if (run_tmb == T) {

  compile(here::here("scripts", paste0(script_name, ".cpp")), "-O0") # what is the -O0?

  dyn.load(dynlib(here::here("scripts", script_name)))

  # ahnold_model <-
  #   MakeADFun(
  #     ahnold_data,
  #     ahnold_params,
  #     DLL = script_name,
  #     random = c(
  #       "seen_region_cluster_betas",
  #       "seeing_region_cluster_betas"
  #     ))

  ahnold_model <-
    MakeADFun(
      ahnold_data,
      ahnold_params,
      DLL = script_name,
      random = c(
        "seen_year_species_betas",
        "seeing_year_species_betas",
        "seen_region_cluster_betas",
        "seeing_region_cluster_betas"
      ))

  if (tmb_to_stan == F) {
    a <- Sys.time()
    set.seed(42)
    ahnold_fit <-
      nlminb(
        ahnold_model$par,
        ahnold_model$fn,
        ahnold_model$gr,
        control = list(iter.max = 4000, eval.max = 5000)
      )
    Sys.time() - a

    save(
      file = here::here(run_dir, "ahnold-onestage-tmb-model.Rdata"),
      ahnold_model
    )

    save(file = here::here(run_dir, "ahnold-onestage-tmb-fit.Rdata"), ahnold_fit)


    ahnold_report <- ahnold_model$report()

    sd_report <- sdreport(ahnold_model)

    # sd_report <- sdreport(ahnold_model,getReportCovariance = TRUE, skip.delta.method = TRUE)

    save(
      file = here::here(run_dir, "ahnold-tmb-onestage-sdreport.Rdata"),
      sd_report
    )

    save(
      file = here::here(run_dir, "ahnold-tmb-onestage-report.Rdata"),
      ahnold_report
    )
  } else {
    ahnold_fit_stan <- tmbstan::tmbstan(ahnold_model, chains = 1,refresh = 25)

    save(
      file = here::here(run_dir, "ahnold-tmbtostan-model.Rdata"),
      ahnold_model
    )

    save(
      file = here::here(run_dir, "ahnold-tmbtostan-fit.Rdata"),
      ahnold_fit_stan
    )
  }
} else {
  load(file = here::here(run_dir, "ahnold-tmb-model.Rdata"))

  load(file = here::here(run_dir, "ahnold-tmb-fit.Rdata"))

  load(file = here::here(run_dir, "ahnold-tmb-report.Rdata"))
}


# pull out estimates ------------------------------------------------------


ahnold_estimates <-
  summary(sd_report) %>%
  as.data.frame() %>%
  mutate(variable = rownames(.)) %>%
  set_names(tolower) %>%
  rename(std_error = `std. error`) %>%
  mutate(
    lower = estimate - 1.96 * std_error,
    upper = estimate + 1.96 * std_error
  )

seen_non_nested_betas <- ahnold_estimates %>%
  filter(str_detect(variable, "seen_non_nested_betas")) %>%
  mutate(group = variable) %>%
  mutate(variable = colnames(x_seen_non_nested))

seeing_non_nested_betas <- ahnold_estimates %>%
  filter(str_detect(variable, "seeing_non_nested_betas")) %>%
  mutate(group = variable) %>%
  mutate(variable = colnames(x_seeing_non_nested))


seen_region_cluster_betas <- ahnold_estimates %>%
  filter(str_detect(variable, "seen_region_cluster_betas")) %>%
  mutate(group = variable) %>%
  mutate(variable = colnames(x_seen_region_cluster))

seeing_region_cluster_betas <- ahnold_estimates %>%
  filter(str_detect(variable, "seeing_region_cluster_betas")) %>%
  mutate(group = variable) %>%
  mutate(variable = colnames(x_seeing_region_cluster))

seen_year_species_betas <- ahnold_estimates %>%
  filter(str_detect(variable, "seen_year_species_betas")) %>%
  mutate(group = variable) %>%
  mutate(variable = colnames(x_seen_year_species))

seeing_year_species_betas <- ahnold_estimates %>%
  filter(str_detect(variable, "seeing_year_species_betas")) %>%
  mutate(group = variable) %>%
  mutate(variable = colnames(x_seeing_year_species))


seen_trend <- seen_year_species_betas %>%
  mutate(classcode = str_replace_all(variable, "(year_classcode-)|(-)|(\\d)",'')) %>%
  mutate(type = 'expected abundance')

seeing_trend <- seeing_year_species_betas %>%
  mutate(classcode = str_replace_all(variable, "(year_classcode-)|(-)|(\\d)",'')) %>%
  mutate(type = 'probability seen')

seen_trend %>%
  filter(classcode == 'hsem') %>%
  mutate(year = str_replace_all(variable,'\\D','') %>% as.numeric()) %>%
  ggplot() +
  geom_line(aes(year, estimate, color = classcode), show.legend = F) +
  facet_wrap(~classcode)


did_betas <- ahnold_estimates %>%
  filter(str_detect(variable, "net_did")) %>%
  mutate(group = variable) %>%
  mutate(variable = unique(data$year) %>% as.character()) %>%
  ungroup()


betas <- bind_rows(
  seen_non_nested_betas,
  seeing_non_nested_betas,
  did_betas
) %>%
  as_data_frame()

betas %>%
  filter(variable != 'intercept') %>%
  ggplot() +
  geom_pointrange(aes(
    variable,
    y = estimate,
    ymin = lower,
    ymax = upper
  )) +
  geom_hline(aes(yintercept = 0)) +
  coord_flip() +
  facet_wrap(~ group, scales = "free")

seen_region_cluster_betas %>%
  ggplot() +
  geom_pointrange(aes(
    variable,
    y = estimate,
    ymin = lower,
    ymax = upper
  )) +
  geom_hline(aes(yintercept = 0)) +
  coord_flip() +
  facet_wrap(~ group, scales = "free")


did_betas %>%
  ggplot() +
  geom_pointrange(aes(
    variable,
    y = estimate,
    ymin = lower,
    ymax = upper
  )) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 2003))




# diagnostics -------------------------------------------------------------

ahnold_betas <-
  data_frame(beta = ahnold_fit$par, variable = names(ahnold_fit$par))

seen_data$log_density_hat <- ahnold_report$log_density_hat %>% as.numeric()

seen_data %>%
  ggplot(aes(log_density, log_density_hat)) +
  geom_point()


lm(log_density ~ log_density_hat, data = seen_data) %>% summary()

seen_data %>%
  ggplot(
    aes(
      log_density_hat,
      log_density_hat - log_density    )
  ) +
  geom_point()

a <- ahnold_model$report()$prob_seeing

seeing_data$prob_seen <- a %>% as.numeric()

seeing_data %>%
  ggplot(aes(any_seen, prob_seen)) +
  geom_boxplot()

