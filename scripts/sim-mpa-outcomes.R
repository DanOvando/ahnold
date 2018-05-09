
# setup -------------------------------------------------------------------

library(tidyverse)
library(spasm)
library(FishLife)
library(hrbrthemes)
library(scales)
functions <- list.files(here::here("functions"))

walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions


# options -----------------------------------------------------------------

sim_years <-  50

burn_years <- 25

num_patches <- 10


# prepare data -----------------------------------------------------

run_name <- 'v2.1'

run_dir <- file.path('results', run_name)

load(file = file.path(run_dir, "rawish_ahnold_data.Rdata"))

load(file = file.path(run_dir, "tmb_runs.Rdata"))


fitted_data <- tmb_runs$data[[1]]

cdfw_data <- read_csv(file = file.path('data','cfdw-catches.csv'))


has_timeseries <- cdfw_data %>%
  group_by(sci_name) %>%
  summarise(has_catch = min(year) <=2000 & max(year) == 2015) %>%
  filter(has_catch == T)

fill_catches <- function(min_year, max_year, catches){

  full_frame <- data_frame(year = min_year:max_year) %>%
    left_join(catches, by = 'year') %>%
    mutate(catch = zoo::na.approx(catch_lbs, rule = 2)) %>%
    select(-catch_lbs)

}

cdfw_catches <- cdfw_data %>%
  filter(sci_name %in% unique(has_timeseries$sci_name)) %>%
  group_by(sci_name,year) %>%
  summarise(catch_lbs = sum(pounds_caught)) %>%
  group_by(sci_name) %>%
  nest(-sci_name, .key = 'catches') %>%
  filter(is.na(sci_name) == F) %>%
  mutate(catches = map(catches, fill_catches, min_year = 2000, max_year = 2015))

seen_species <- life_history_data %>%
  filter(classcode %in% (fitted_data$classcode %>% unique())) %>%
  rename(sci_name = taxa,
         linf = vbgf.linf,
         common_name = commonname) %>%
  mutate(sci_name = tolower(sci_name)) %>%
  filter(sci_name %in% has_timeseries$sci_name, is.na(linf) == F)



# prepare experiments -----------------------------------------------------


sim_grid <- expand.grid(
  scientific_name = unique(fitted_data$taxa)[1],
  steepness = seq(0.6, 1, by = .2),
  adult_movement = seq(1, 20, length.out = 3),
  larval_movement = seq(1, 20, length.out = 3),
  density_movement_modifier = c(0,1),
  density_dependence_form = 1:5,
  mpa_size = c(.1,.3,.75),
  f_v_m = seq(.2,2, by = 0.5),
  fleet_model = c('constant-effort'),
  effort_allocation = c("profit-gravity","simple"),
  stringsAsFactors = F
)

# create fish objects
sim_grid <- sim_grid %>%
  sample_n(200) %>%
  mutate(fish = pmap(list(
    scientific_name = scientific_name,
    steepness = steepness,
    adult_movement = adult_movement,
    larval_movement = larval_movement,
    density_dependence_form = density_dependence_form,
    density_movement_modifier = density_movement_modifier
  ), create_fish))

# create fleet objects

sim_grid <- sim_grid %>%
  mutate(fleet = pmap(
    list(
      fish = fish,
      fleet_model = fleet_model,
      effort_allocation = effort_allocation
    ),
    create_fleet))

# tune fleet objects
sim_grid <- sim_grid %>%
  mutate(fleet =  pmap(
    list(
      f_v_m = f_v_m,
      fish = fish,
      fleet = fleet
    ),
    tune_fleet))


# run experiments ---------------------------------------------------------

a <- Sys.time()

sim_grid <- sim_grid %>%
  mutate(mpa_experiment = pmap(
    list(
      fish = fish,
      fleet = fleet,
      mpa_size = mpa_size
    ),
    comp_foo,
   year_mpa = 15,
   sim_years = 50,
   burn_years = 10,
   num_patches = 20
  ))

Sys.time() - a

calc_mpa_effect <- function(outcomes){


  mpa_effect <- outcomes %>%
    group_by(year) %>%
    mutate(mpa_size = max(percent_mpa)) %>%
    ungroup() %>%
    select(year, experiment, biomass,mpa_size) %>%
    spread(experiment, biomass) %>%
    mutate(mpa_effect = `with-mpa`/`no-mpa` - 1) # %>%
    # select(year, mpa_size, mpa_effect)

}

sim_grid <- sim_grid %>%
  mutate(mpa_effect = map(map(mpa_experiment,"outcomes"), calc_mpa_effect))

outcomes <-  sim_grid %>%
  select(-fish, -fleet,-mpa_experiment) %>%
  unnest()

outcomes %>%
  ggplot(aes(year, mpa_effect,group = year, fill = mpa_size)) +
  geom_boxplot() +
  # coord_cartesian(ylim = c(0,2)) +
  scale_y_continuous(labels = scales::percent)


# process outcomes --------------------------------------------------------


# save outcomes -----------------------------------------------------------


