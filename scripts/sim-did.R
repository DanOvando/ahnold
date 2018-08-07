library(scales)
library(scales)
library(viridis)
library(ggmap)
library(forcats)
library(stringr)
library(lubridate)
library(purrr)
library(lme4)
library(TMB)
library(FishLife)
library(patchwork)
library(spasm)
library(doParallel)
library(tidyverse)

functions <- list.files(here::here("functions"))

walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions

run_name <- 'v3.0'

run_dir <- file.path('results', run_name)

simulate_samples <- T

burn_years <- 25

sim_years <- 75

year_mpa <- 30

num_patches <-  2

mpa_size <- 0.5

n_samples <- 5

time_step <-  1


plot_theme <- hrbrthemes::theme_ipsum(base_size = 14,
                                      axis_title_size = 16)


theme_set(plot_theme)

# prepare data ------------------------------------------------------------

load(paste0(run_dir, '/abundance_data.Rdata'))

load(paste0(run_dir, '/rawish_zissou_data.Rdata'))

enso <- read_csv('data/enso.csv')

pisco <- abundance_data$data[[1]]

pisco_fish <- life_history_data %>%
  filter(classcode %in% unique(pisco$classcode)) %>%
  mutate(enviro_effect = ifelse(geographic_cluster > 1,-1, 1)) %>%
  mutate(enviro = list(enso$enso[(nrow(enso) - (sim_years + burn_years)):(nrow(enso) - 1)])) %>%
  mutate(enviro = map2(enviro, enviro_effect, ~ .x * .y))

n_groups <- 5

simple_fish <-
  data_frame(
    loo = c(rnorm(n_groups, 120, 10), rnorm(n_groups, 100, 10)),
    k = 0.4,
    lm = .75 * loo,
    m = 0.2,
    targeted = rep(c(1, 0), each = n_groups),
    classcode = fruit[1:(n_groups * 2)],
    commonname = colors()[1:(n_groups * 2)],
    enviro = NA
  )


diver <- list(q = .1,
              sel_size_50 = 2 ,
              sel_size_delta = 2)



pisco_divers <- data_frame(diver = fruit[1:3],
                           diver_stats = list(
                             list(
                               q = .01,
                               sel_size_50 = 2 ,
                               sel_size_delta = 2
                             ),
                             b = list(
                               q = .067,
                               sel_size_50 = 10 ,
                               sel_size_delta = 2
                             ),
                             d = list(
                               q = .1,
                               sel_size_50 = 4 ,
                               sel_size_delta = 2
                             )
                           ))

if (simulate_samples == T) {
  simple_fish <- create_samples(
    fishes = simple_fish,
    divers = pisco_divers %>% slice(1),
    mpa_size = mpa_size,
    burn_years = burn_years,
    sim_years = sim_years,
    year_mpa = year_mpa,
    num_patches = num_patches,
    samples = n_samples,
    rec_driver = 'stochastic',
    enviro_strength = 1,
    sigma_r = 0,
    cores = 1,
    time_step = time_step
  )

  pisco_fish <- create_samples(
    fishes = pisco_fish,
    divers = pisco_divers,
    mpa_size = mpa_size,
    burn_years = burn_years,
    sim_years = sim_years,
    year_mpa = year_mpa,
    num_patches = num_patches,
    samples = n_samples,
    rec_driver = 'environment',
    enviro_strength = 1,
    sigma_r = 0.25,
    cores = 4,
    time_step = time_step
  )

  save(file = paste0(run_dir, '/simulated-data.Rdata'),
       simple_fish,
       pisco_fish)

} else {
  load(file = paste0(run_dir, '/simulated-data.Rdata'))
}


# fit simple model --------------------------------------------------------


simple_performance <-
  test_performance(
    simple_fish,
    year_mpa = year_mpa + burn_years,
    min_year = 50,
    time_step = time_step
  )

pisco_performance <-
  test_performance(pisco_fish,
                   year_mpa + burn_years,
                   min_year = 50,
                   time_step = time_step)

save(file = 'about_time.Rdata', simple_performance, pisco_performance)


bayes_model <-
  rstanarm::stan_glmer(
    log_density ~ logical_targeted + factor_year +  logical_targeted:factor_year + (1 |
                                                                                      classcode),
    data = simple_data %>% slice(1:400) ,
    chains = 1,
    cores = 4
  )


check_stan_did <- broom::tidy(test) %>%
  filter(str_detect(term, 'logical_targetedTRUE:')) %>%
  mutate(year = str_replace_all(term, '\\D', '') %>% as.numeric()) %>%
  ggplot() +
  geom_pointrange(
    aes(
      x = year,
      y = estimate,
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    )
  ) +
  geom_point(data = true_effect, aes(year, mpa_effect), color = 'red')



script_name <- "fit_simmpa"

sfa <- safely(fit_sim_mpa)


tmb_runs <-
  data_frame(
    data = list(simple_data),
    non_nested_variables = list(c('targeted')),
    include_intercept = c(TRUE),
    fixed_did = c(TRUE)
  )


tmb_runs <- tmb_runs %>%
  mutate(tmb_fit = pmap(
    list(
      data = data,
      non_nested_variables = non_nested_variables,
      include_intercept = include_intercept,
      fixed_did = fixed_did
    ),
    sfa,
    run_dir = run_dir,
    script_name = script_name,
    fixed_regions = T
  ))

huh <- tmb_runs$tmb_fit[[1]]$result

a <- huh$ahnold_estimates %>%
  filter(variable == 'net_did')

a %>%
  mutate(year = simple_data$year %>% unique()) %>%
  ggplot(aes(year, estimate)) +
  geom_point()
