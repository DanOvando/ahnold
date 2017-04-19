rm(list = ls())
library(tidyverse)
library(spasm)
library(rfishbase)
library(multidplyr)
library(trelliscopejs)
library(hrbrthemes)
library(scales)

comp_foo <-
  function(fish,
           fleet,
           year_mpa,
           mpa_size,
           sim_years,
           num_patches,
           burn_years,
           run_id,
           num_runs) {

    write(paste('run', run_id, 'out of',num_runs), file = 'simprog.txt',append = T,
          ncolumns = 1)

    no_mpa <-
      sim_fishery(
        fish = fish,
        fleet = fleet,
        manager = create_manager(year_mpa = year_mpa, mpa_size = 0),
        sim_years = sim_years,
        num_patches = num_patches,
        burn_years = burn_years
      ) %>%
      mutate(experiment = 'no-mpa')

    wi_mpa <-
      sim_fishery(
        fish = fish,
        fleet = fleet,
        manager = create_manager(year_mpa = year_mpa, mpa_size = mpa_size),
        sim_years = sim_years,
        num_patches = num_patches,
        burn_years = burn_years
      ) %>%
      mutate(experiment = 'with-mpa')

    outcomes <- no_mpa %>%
      bind_rows(wi_mpa) %>%
      group_by(year, experiment) %>%
      summarise(
        ssb = sum(ssb),
        percent_mpa = mean(mpa),
        catch = sum(biomass_caught),
        profits = sum(profits),
        effort = sum(effort)
      )



  }

sim_grid <- expand.grid(
  max_age = 25,
  steepness = seq(0.3, .9, length.out = 2),
  adult_movement = seq(1, 100, length.out = 3),
  larval_movement = seq(1, 100, length.out = 3),
  density_dependence_form = 1:5,
  m = seq(0.3, .8, by = .3),
  lhi_type = 1:3,
  target_catch = seq(1, 1000, length.out = 5),
  mpa_size = seq(0, 1, length.out = 4),
  fleet_model = c('constant-catch'),
  stringsAsFactors = F
)

sim_grid <- sim_grid %>%
  slice(1:10) %>%
  mutate(
    fish = pmap(
      list(
        max_age = max_age,
        steepness = steepness,
        adult_movement = adult_movement,
        larval_movement = larval_movement,
        m = m,
        lhi_type = lhi_type,
        density_dependence_form = density_dependence_form
      ),
      create_fish,
      scientific_name = 'fakeish fishis',
      query_fishbase = F,
      linf = 100,
      t0 = 0
    ),
    initial_effort = target_catch
  )

sim_grid <- sim_grid %>%
  mutate(fleet = pmap(
    list(
      target_catch = target_catch,
      fish = fish,
      fleet_model = fleet_model,
      initial_effort = initial_effort
    ),
    create_fleet,
    price = 100
  ))


year_mpa <- 25


cluster <- create_cluster(cores = 4)

file.remove('simprog.txt')

a <- proc.time()


sim_grid_outcomes <- sim_grid %>%
  mutate(run_id = 1:dim(.)[1],
         num_runs = dim(.)[1]) %>%
  mutate(year_mpa = year_mpa)  %>%
  partition(run_id, cluster = cluster) %>%
  cluster_library("tidyverse") %>%
  cluster_library("spasm") %>%
  cluster_assign_value("comp_foo", comp_foo) %>%
  mutate(outcomes = pmap(
    list(
      fish = fish,
      fleet = fleet,
      year_mpa = year_mpa,
      mpa_size = mpa_size,
      run_id = run_id,
      num_runs = num_runs
    ),
    comp_foo,
    sim_years = 50,
    burn_years = 25,
    num_patches = 50
  )) %>%
  collect() %>%
  as_tibble()

proc.time() - a



plot_foo <- function(outcomes, year_mpa) {
  min_year <- min(outcomes$year)

  mpa_size <-  max(outcomes$percent_mpa)

  outcomes %>%
    filter(year > (min_year + year_mpa - 3) &
             year < (min_year + year_mpa + 10)) %>%
    ggplot(aes(year, ssb)) +
    geom_line(aes(color = experiment)) +
    geom_point(size = 2, shape = 21, aes(fill = experiment)) +
    geom_vline(aes(xintercept = year_mpa + min_year - 1)) +
    theme_ipsum() +
    labs(caption = paste0(100 * mpa_size, '% MPA'))
}

sim_grid_plots <- sim_grid_outcomes %>%
  mutate(sim_plot = map2_plot(outcomes, year_mpa, plot_foo)) %>%
  mutate(eq_f = map_dbl(fleet, 'eq_f')) %>%
  select(-fish, -fleet, -outcomes)

trelliscope(sim_grid_plots, name = 'arg')