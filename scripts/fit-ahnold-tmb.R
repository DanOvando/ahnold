set.seed(666)
library(TMB)
library(stringr)
library(purrr)
library(rstan)
library(tidyverse)
demons::load_functions()

rstan_options(auto_write = TRUE)

run_tmb <- FALSE

tmb_to_stan <- FALSE # fit the model in stan instead of TMB

run_tmb <-  T

run_dir <- file.path('results', run_name)

run_name <- "Working"

run_dir <- file.path("results", run_name)

load(file = paste0(run_dir, "/abundance_indices.Rdata"))


data <- abundance_indices %>%
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
  mutate(temp_deviation = (mean_temp - temperature)^2)

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



numeric_species_key <-
  data_frame(classcode = unique(data$classcode)) %>%
  arrange(classcode) %>%
  mutate(numeric_classcode = 1:nrow(.))

seen_data <- data %>%
  filter(any_seen == T) %>%
  left_join(numeric_species_key, by = "classcode") %>%
  select(non_nested_variables, log_density) %>%
  na.omit()

log_density <- seen_data$log_density

seen_data <- seen_data %>%
  select(-log_density)

seeing_data <- data %>%
  left_join(numeric_species_key, by = "classcode") %>%
  select(non_nested_variables, any_seen) %>%
  na.omit()

any_seen <- seeing_data$any_seen

seeing_data <- seeing_data %>%
  select(-any_seen)


x_seen_non_nested <- seen_data %>%
  mutate(intercept = 1) %>%
  spread_factors(drop_one = T) %>%
  purrrlyr::dmap(center_scale) %>%
  mutate(cumulative_n_obs_2 = cumulative_n_obs ^ 2) %>%
  mutate(targeted = seen_data$targeted)

x_seen_did <- seen_data %>%
  select(targeted, factor_year) %>%
  mutate(index = 1:nrow(.)) %>%
  spread(factor_year, targeted, fill = 0) %>%
  arrange(index) %>%
  select(-(1:2)) # drop index and one factor level

# x_seen <-
#   bind_cols(seen_data %>% select(log_density),
#             x_seen_non_nested, x_seen_did) #,
#             # x_seen_region_cluster)


# prepare seeing ----------------------------------------------------------

x_seeing_non_nested <- seeing_data %>%
  mutate(intercept = 1) %>%
  spread_factors(drop_one = T) %>%
  purrrlyr::dmap(center_scale) %>%
  mutate(cumulative_n_obs_2 = cumulative_n_obs ^ 2)

x_seeing_did <- seeing_data %>%
  select(targeted, factor_year) %>%
  mutate(index = 1:nrow(.)) %>%
  spread(factor_year, targeted, fill = 0) %>%
  arrange(index) %>%
  select(-(1:2)) # drop index and one factor level

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
  select(colnames(x_seeing_non_nested)) %>%
  slice(-1) # make sure things are in the same column order as they used to be

standard_did <- data_frame(year = unique(seeing_data$factor_year)) %>%
  spread_factors(drop_one = T) %>%
  slice(-1)


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
  x_seeing_non_nested = x_seeing_non_nested,
  x_seeing_did = x_seeing_did,
  log_density = log_density,
  any_seen = any_seen,
  standard_non_nested = standard_non_nested,
  standard_did = standard_did,
  seen_species_index = seen_species_index
)

ahnold_data <- map_if(ahnold_data, is.data.frame, ~ as.matrix(.x))

any_na <- map_lgl(ahnold_data, ~ any(is.na(.x))) %>% any()

if (any_na) {
  stop("NAs in ahnold_data")
}

ahnold_params <- list(
  seen_non_nested_betas = rep(0, ncol(x_seen_non_nested)),
  seen_did_betas = rep(0, ncol(x_seen_did)),
  seen_density_species_sigma = rep(log(1), n_species),
  seeing_non_nested_betas = rep(0, ncol(x_seeing_non_nested)),
  seeing_did_betas = rep(0, ncol(x_seeing_did))
)

any_na <- map_lgl(ahnold_params, ~ any(is.na(.x))) %>% any()

if (any_na) {
  stop("NAs in ahnold_params")
}



if (run_tmb == T) {
  script_name <- "fit_ahnold"
  compile(here::here("scripts", paste0(script_name, ".cpp")), "-O0") # what is the -O0?

  dyn.load(dynlib(here::here("scripts", script_name)))

  ahnold_model <-
    MakeADFun(
      ahnold_data,
      ahnold_params,
      DLL = script_name
    )

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
    print("made it here")

    ahnold_report <- ahnold_model$report()

    sd_report <- sdreport(ahnold_model)

    # sd_report <- sdreport(ahnold_model,getReportCovariance = TRUE, skip.delta.method = TRUE)

    save(
      file = here::here(run_dir, "ahnold-tmb-sdreport.Rdata"),
      sd_report
    )

    save(
      file = here::here(run_dir, "ahnold-tmb-report.Rdata"),
      ahnold_report
    )
  } else {
    ahnold_fit_stan <- tmbstan::tmbstan(ahnold_model, chains = 1)

    save(
      file = here::here(run_dir, "ahnold-tmbtostan-model.Rdata"),
      ahnold_model
    )

    save(
      file = here::here(run_dir, "ahnold-tmbtostan-fit.Rdata"),
      ahnold_fit
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

did_betas <- ahnold_estimates %>%
  filter(str_detect(variable, "net_did")) %>%
  mutate(group = variable) %>%
  mutate(variable = colnames(standard_did))


betas <- bind_rows(
  seen_non_nested_betas,
  seeing_non_nested_betas,
  did_betas
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

