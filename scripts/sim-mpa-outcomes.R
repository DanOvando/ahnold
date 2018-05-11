
# setup -------------------------------------------------------------------

library(spasm)
library(FishLife)
library(hrbrthemes)
library(scales)
library(doParallel)
library(tidyverse)

functions <- list.files(here::here("functions"))

walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions


# options -----------------------------------------------------------------

sim_years <- 50

burn_years <- 20

num_patches <- 25

run_experiments <- TRUE

n_cores <- 3

samps <- 5000

grid_search <-  FALSE

# prepare data -----------------------------------------------------

run_name <- "v2.1"

run_dir <- file.path("results", run_name)

experiment_dir <- here::here(run_dir, "experiments")

load(file = file.path(run_dir, "rawish_ahnold_data.Rdata"))

load(file = file.path(run_dir, "tmb_runs.Rdata"))

fitted_data <- tmb_runs$data[[1]]

cdfw_data <- read_csv(file = file.path("data", "cfdw-catches.csv"))

has_timeseries <- cdfw_data %>%
  group_by(sci_name) %>%
  summarise(has_catch = min(year) <= 2000 & max(year) == 2015) %>%
  filter(has_catch == T)

fill_catches <- function(min_year, max_year, catches) {
  full_frame <- data_frame(year = min_year:max_year) %>%
    left_join(catches, by = "year") %>%
    mutate(catch = zoo::na.approx(catch_lbs, rule = 2)) %>%
    select(-catch_lbs)
}

cdfw_catches <- cdfw_data %>%
  filter(sci_name %in% unique(has_timeseries$sci_name)) %>%
  group_by(sci_name, year) %>%
  summarise(catch_lbs = sum(pounds_caught)) %>%
  group_by(sci_name) %>%
  nest(-sci_name, .key = "catches") %>%
  filter(is.na(sci_name) == F) %>%
  mutate(catches = map(catches, fill_catches, min_year = 2000, max_year = 2015))

seen_species <- life_history_data %>%
  filter(classcode %in% (fitted_data$classcode %>% unique())) %>%
  rename(
    sci_name = taxa,
    linf = vbgf.linf,
    common_name = commonname
  ) %>%
  mutate(sci_name = tolower(sci_name)) %>%
  filter(sci_name %in% has_timeseries$sci_name, is.na(linf) == F)



# prepare experiments -----------------------------------------------------




# run experiments ---------------------------------------------------------

if (run_experiments == T) {


  if (grid_search == T) {

  sim_grid <- expand.grid(
    scientific_name = unique(fitted_data$taxa),
    steepness = seq(0.6, 1, by = .2),
    adult_movement = seq(1, 20, length.out = 3),
    larval_movement = seq(1, 20, length.out = 3),
    density_movement_modifier = c(0, 1),
    density_dependence_form = 1:5,
    mpa_size = c(.1, .3, .75),
    f_v_m = seq(.01, 1.25, by = 0.5),
    fleet_model = c("constant-effort"),
    effort_allocation = c("profit-gravity", "simple"),
    stringsAsFactors = F
  )
  } else{


    sim_grid <- data_frame(scientific_name = sample(unique(fitted_data$taxa), samps, replace = T),
                           steepness = runif(samps, min = 0.6, max = 1),
                           adult_movement = sample(1:num_patches,samps, replace = T),
                           larval_movement = sample(1:num_patches,samps, replace = T),
                           density_movement_modifier = sample(c(0,1),samps, replace = T),
                           density_dependence_form = sample(1:5,samps, replace = T),
                           mpa_size = runif(samps, min = 0.05, max = 1),
                           f_v_m = runif(samps, min = 0.01, max = 1.25),
                           fleet_model = c("constant-effort"),
                           effort_allocation = sample(c("profit-gravity", "simple"), samps, replace = T),
                           year_mpa = sim_years / 2
                           # year_mpa = sample(sim_years/2, samps, replace = T)
                           )

  }

  # create fish objects
  sim_grid <- sim_grid %>%
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
      create_fleet
    ))

  # tune fleet objects
  sim_grid <- sim_grid %>%
    mutate(fleet = pmap(
      list(
        f_v_m = f_v_m,
        fish = fish,
        fleet = fleet
      ),
      tune_fleet,
      num_patches = num_patches
    ))


  doParallel::registerDoParallel(cores = n_cores)

  foreach::getDoParWorkers()

  sim_grid$experiment <- 1:nrow(sim_grid)


  if (dir.exists(experiment_dir) == F) {
    dir.create(experiment_dir, recursive = T)
  }

  mpa_experiments <- foreach::foreach(i = 1:nrow(sim_grid)) %dopar% {
    results <- sim_grid %>%
      slice(i) %>%
      mutate(mpa_experiment = pmap(
        list(
          fish = fish,
          fleet = fleet,
          mpa_size = mpa_size,
          year_mpa = year_mpa
        ),
        comp_foo,
        sim_years = sim_years,
        burn_years = burn_years,
        num_patches = num_patches
      ))


    calc_mpa_effect <- function(outcomes) {
      mpa_effect <- outcomes %>%
        group_by(year) %>%
        mutate(mpa_size = max(percent_mpa)) %>%
        ungroup() %>%
        select(year, experiment, biomass, mpa_size) %>%
        spread(experiment, biomass) %>%
        mutate(mpa_effect = `with-mpa` / `no-mpa` - 1) # %>%
      # select(year, mpa_size, mpa_effect)
    }

    results <- results %>%
      mutate(mpa_effect = map(map(mpa_experiment, "outcomes"), calc_mpa_effect))


    filename <- glue::glue("experiment_{i}.rds")

    saveRDS(results, file = glue::glue("{experiment_dir}/{filename}"))

    # started at 9:41 pm
  } # close dopar
}

loadfoo <- function(experiment, experiment_dir, output = "mpa-effect") {
  ex <- readRDS(glue::glue("{experiment_dir}/experiment_{experiment}.rds"))

  if (output == "mpa-effect") {
    ex <- ex %>%
      select(-mpa_experiment)
  }
  return(ex)
}

processed_grid <- map_df(1:samps, loadfoo, experiment_dir = experiment_dir)

save(processed_grid,file = paste0(run_dir,"/processed_grid.Rdata"))

# outcomes <- processed_grid %>%
#   select(-fish, -fleet) %>%
#   unnest()
#
# outcomes %>%
#   ggplot(aes(
#     year,
#     mpa_effect,
#     group = interaction(year, factor(density_movement_modifier)),
#     fill = factor(density_movement_modifier),
#     color = factor(density_movement_modifier)
#   )) +
#   geom_boxplot() +
#   scale_y_continuous(labels = scales::percent) +
#   scale_x_continuous(limits = c(40, NA))


# process outcomes --------------------------------------------------------


# save outcomes -----------------------------------------------------------
