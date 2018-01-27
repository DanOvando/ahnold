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

if (("demons" %in% installed.packages()) == F) {
  devtools::install_github('danovando/demons')
}

demons::load_functions('functions')

run_name <- 'v1.0'

run_dir <- file.path('results', run_name)

simulate_samples <- T

burn_years <- 10

sim_years <- 100

year_mpa <- 75

num_patches <-  2

mpa_size <- 0.5

n_samples <- 5

time_step <-  1
# prepare data ------------------------------------------------------------

load(paste0(run_dir, '/abundance_data.Rdata'))

load(paste0(run_dir, '/rawish_ahnold_data.Rdata'))

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
    cores = 4,
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

test_performance <-
  function(fishes,
           year_mpa,
           min_year = 75,
           max_year = 100,
           time_step = 1) {
    simple_data <- fishes %>%
      select(loo, k, lm, m, targeted, classcode, commonname, pisco_samples) %>%
      unnest() %>%
      filter(year > min_year, year < max_year) %>%
      select(-pop, -sampled_lengths, -diver_stats) %>%
      mutate(year = year * time_step) %>%
      mutate(subyear = year - floor(year),
             year = floor(year)) %>%
      mutate(
        factor_year = as.factor(year),
        log_density = log(density),
        logical_targeted = targeted > 0,
        any_seen = density > 0,
        region = patches,
        post_mpa = year >= year_mpa
      )


    true_effect <- fishes %>%
      select(classcode, targeted, mpa_effect) %>%
      unnest() %>%
      filter(year > min_year, targeted == 1) %>%
      mutate(post_mpa = year > year_mpa) %>%
      mutate(year = year * time_step) %>%
      mutate(subyear = year - floor(year),
             year = floor(year))

    bare_bones_model <-
      lm(log_density ~ targeted + factor_year + targeted:factor_year + loo,
         data = simple_data)

    if (n_distinct(simple_data$diver) > 1) {
      # mixed_effect_model <- lme4::lmer(log_density ~ targeted + factor_year + targeted:factor_year + diver + (factor_year |classcode), data = simple_data)
      mixed_effect_model <-
        rstanarm::stan_glmer(
          log_density ~  targeted:factor_year + diver + (factor_year  - 1|classcode),
          data = simple_data,
          iter = 1000,
          chains = 4,
          cores = 4,
          refresh = 1,
          QR = TRUE,
          control = list(max_treedepth = 8)
        )

    } else
      {
      # mixed_effect_model <-
      #   lme4::lmer(
      #     log_density ~ targeted + factor_year + targeted:factor_year + (factor_year |
      #                                                                      classcode),
      #     data = simple_data
      #   )

      mixed_effect_model <-
        rstanarm::stan_glmer(
          log_density ~ targeted + factor_year + targeted:factor_year + (factor_year | classcode),
          data = simple_data %>% select(log_density, targeted, factor_year, classcode),
          iter = 1000,
          chains = 4,
          cores = 4,
          refresh = 1,
          QR = TRUE,
          control = list(max_treedepth = 8)
        )

    }

    pre_post_model <-
      lm(log_density ~ targeted + post_mpa + targeted:post_mpa, data = simple_data)


    bare_bones_did <- broom::tidy(bare_bones_model) %>%
      filter(str_detect(term, 'targeted:')) %>%
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
      geom_smooth(data = true_effect, aes(year, mpa_effect), color = 'red') +
      labs(title = 'bare bones')

    mixed_effect_did <- broom::tidy(mixed_effect_model) %>%
      filter(str_detect(term, 'targeted:')) %>%
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
      geom_smooth(data = true_effect, aes(year, mpa_effect), color = 'red') +
      labs(title = 'mixed effects')

    check_block_did <- broom::tidy(pre_post_model) %>%
      filter(str_detect(term, 'targeted:')) %>%
      mutate(post_mpa = TRUE) %>%
      ggplot() +
      geom_boxplot(data = true_effect,
                   aes(post_mpa, mpa_effect),
                   color = 'red') +
      geom_pointrange(
        aes(
          x = post_mpa,
          y = estimate,
          ymin = estimate - 1.96 * std.error,
          ymax = estimate + 1.96 * std.error
        )
      )

    out_plot <-
    {
      bare_bones_did + mixed_effect_did + plot_layout(ncol = 1)
    } + check_block_did + plot_layout(ncol = 2)

    out <- list(
      out_plot = out_plot,
      bare_bones_model = bare_bones_model,
      mixed_effect_model = mixed_effect_model,
      pre_post_model = pre_post_model
    )

    return(out)
  }


simple_performance <-
  test_performance(
    simple_fish,
    year_mpa = year_mpa,
    min_year = 75,
    time_step = time_step
  )

pisco_performance <-
  test_performance(pisco_fish,
                   year_mpa,
                   min_year = 75,
                   time_step = time_step)

save(file = 'about_time.Rdata', simple_performance, pisco_performance)

# tmb side ----------------------------------------------------------------


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
