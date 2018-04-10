library(tidyverse)
library(spasm)
library(FishLife)
library(multidplyr)
# library(trelliscopejs)
library(hrbrthemes)
library(scales)

sim_years <-  50

burn_years <- 25

num_patches <- 10

comp_foo <- function(fish,
                     fleet,
                     year_mpa,
                     mpa_size,
                     sim_years,
                     num_patches,
                     burn_years,
                     run_id,
                     num_runs) {
  write(
    paste('run', run_id, 'out of', num_runs),
    file = 'simprog.txt',
    append = T,
    ncolumns = 1
  )

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


# load data ---------------------------------------------------------------

run_name <- 'Working'

run_dir <- file.path('results', run_name)

load(file = file.path(run_dir, "/rawish_ahnold_data.Rdata"))

load(file = file.path(run_dir, "/abundance_indices.Rdata"))


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

cdfw_data %>%
  filter(sci_name == 'semicossyphus pulcher') %>%
  group_by(year) %>%
  summarise(catch = sum(pounds_caught))

cdfw_catches <- cdfw_data %>%
  filter(sci_name %in% unique(has_timeseries$sci_name)) %>%
  group_by(sci_name,year) %>%
  summarise(catch_lbs = sum(pounds_caught)) %>%
  group_by(sci_name) %>%
  nest(-sci_name, .key = 'catches') %>%
  filter(is.na(sci_name) == F) %>%
  mutate(catches = map(catches, fill_catches, min_year = 2000, max_year = 2015))

seen_species <- life_history_data %>%
  filter(classcode %in% (abundance_indices$classcode %>% unique())) %>%
  rename(sci_name = taxa,
         linf = vbgf.linf,
         common_name = commonname) %>%
  mutate(sci_name = tolower(sci_name)) %>%
  filter(sci_name %in% has_timeseries$sci_name, is.na(linf) == F)


sim_grid <- expand.grid(
  common_name = seen_species$common_name,
  depletion = c(0.2, .4, .8),
  adult_movement = c(.1,10),
  larval_movement = c(.1,10),
  density_dependence_form = 1:5,
  steepness = c(0.4,0.9),
  stringsAsFactors = F
) %>%
  left_join(seen_species, by = 'common_name') %>%
  as_tibble() %>%
  filter(is.na(linf) == F) %>%
  arrange(linf, sci_name) %>%
  mutate(m = 0.2)

sim_grid <- sim_grid %>%
  left_join(cdfw_catches, by = 'sci_name') %>%
  mutate(target_catch = map_dbl(catches,~dplyr::first(.x$catch)))



tune_f <-
  function(params,
           target_depletion,
           linf,
           scientific_name,
           adult_movement,
           larval_movement,
           steepness,
           density_dependence_form,
           target_catch = 100,
           alpha = 0.5,
           sim_years = sim_years,
           burn_years = burn_years,
           num_patches = num_patches) {

    scientific_name <- gsub("(^)([[:alpha:]])", "\\1\\U\\2", scientific_name, perl=TRUE)

    fish <- create_fish(
      scientific_name = scientific_name,
      linf = linf,
      m = NA,
      r0 = exp(params[1]),
      adult_movement = adult_movement,
      larval_movement = larval_movement,
      density_dependence_form = density_dependence_form,
      steepness = steepness,
      query_fishlife = T
    )

    manager <- create_manager(year_mpa = 999)
    unfished <-
      sim_fishery(
        fish = fish,
        fleet = create_fleet(initial_effort = 0, fish = fish, fleet_model = "constant-effort"),
        manager = manager,
        num_patches = num_patches,
        sim_years = sim_years,
        burn_year = burn_years
      )

    ssb0 <- unfished %>%
      filter(year == max(year)) %>%
      summarise(ssb = sum(ssb)) %>% {
        .$ssb
      }

    fished <-
      sim_fishery(
        fish = fish,
        fleet = create_fleet(initial_effort = exp(params[2]), fish = fish, fleet_model = "constant-effort"),
        manager = manager,
        sim_years = sim_years,
        burn_year = burn_years,
        num_patches = num_patches
      )

    ssb <- fished %>%
      filter(year == max(year)) %>%
      summarise(ssb = sum(ssb)) %>% {
        .$ssb
      }

    fish_caught <- fished %>%
      filter(year == max(year)) %>%
      summarise(caught = sum(biomass_caught)) %>% {
        .$caught
      }


    depletion <- ssb / ssb0

    ss <-
      alpha * (log(depletion) - log(target_depletion)) ^ 2 + (1 - alpha) * (log(fish_caught) - log(target_catch)) ^
      2

    print(ss)
    return(ss)
  }


fitfoo <- function(target_depletion,scientific_name, target_catch, linf,sim_years,burn_years,
                   num_patches,
                   larval_movement,
                   adult_movement,
                   density_dependence_form,
                   steepness)
{
  nlminb(
    c(log(1000*target_catch), log(200)),
    tune_f,
    target_depletion = target_depletion,
    target_catch = target_catch,
    scientific_name = scientific_name,
    linf = NA,
    alpha = 0.5,
    sim_years = sim_years,
    burn_years = burn_years,
    num_patches = num_patches,
    adult_movement = adult_movement,
    larval_movement = larval_movement,
    steepness = steepness,
    density_dependence_form = density_dependence_form
  )$par

}

slice_grid <- sim_grid %>%
  # slice(1:5) %>%
  mutate(tuned_pars = pmap(
    list(
      target_depletion = depletion,
      target_catch = target_catch,
      linf = linf,
      scientific_name = sci_name,
      adult_movement = adult_movement,
      larval_movement = larval_movement,
      density_dependence_form = density_dependence_form,
      steepness = steepness
    ),
    fitfoo,
    sim_years = sim_years,
    burn_years = burn_years,
    num_patches = num_patches
  ))

fishfoo <- function(common_name,
                    scientific_name,
                    linf,
                    m,
                    adult_movement,
                    larval_movement,
                    density_dependence_form,
                    steepness,
                    pars) {

  scientific_name <- gsub("(^)([[:alpha:]])", "\\1\\U\\2", scientific_name, perl=TRUE)

  create_fish(
    common_name = common_name,
    scientific_name = scientific_name,
    linf = NA,
    m = NA,
    r0 = exp(pars[1]),
    adult_movement = adult_movement,
    larval_movement = larval_movement,
    density_dependence_form = density_dependence_form,
    steepness = steepness
  )
}

slice_grid <- slice_grid %>%
  mutate(
    fish = pmap(
      list(
        common_name = common_name,
        scientific_name = sci_name,
        linf = linf,
        m = m,
        pars = tuned_pars,
        adult_movement = adult_movement,
        larval_movement = larval_movement,
        density_dependence_form = density_dependence_form,
        steepness = steepness
      ),
      fishfoo
    ),
    fleet = map2(tuned_pars, fish, ~ create_fleet(
      initial_effort = exp(.x[2]), fish = .y
    ))
  )

get_catches <- function(fish, fleet) {
  fished <-
    sim_fishery(
      fish = fish,
      fleet = fleet,
      manager = create_manager(year_mpa = 999),
      num_patches = 10,
      sim_years = 50,
      burn_year = 25
    )

  catches <- fished %>%
    group_by(year) %>%
    summarise(catch = sum(biomass_caught))

  return(as.numeric(catches$catch))

}

slice_grid <- slice_grid %>%
  mutate(historic_catches = map2(fish, fleet, ~ get_catches(fish = .x, fleet = .y))) %>%
  mutate(comp_catches = map2(historic_catches, catches, ~c(.x,.y$catch)))


slice_grid <- slice_grid %>%
  mutate(
    fleet = map2(
      comp_catches,
      fish,
      ~ create_fleet(
        catches = .x,
        fish = .y,
        fleet_model = 'supplied-catch'
      )
    )
  )


slice_grid <- slice_grid %>%
  mutate(sim = map2(
    fish,
    fleet,
    ~ sim_fishery(
      fish = .x,
      fleet = .y,
      manager = create_manager(year_mpa = (sim_years - burn_years) + 4),
      num_patches = num_patches,
      sim_years = burn_years + length(.y$catches),
      burn_year = burn_years
    )
  ),
  mpa_experiment = map2(
    fish,
    fleet,
    ~ comp_foo(fish = .x,
               fleet = .y,
               year_mpa = (sim_years - burn_years) + 4,
               sim_years = sim_years,
               burn_years = burn_years,
               num_patches = num_patches,
               mpa_size = 0.25,
               run_id = 1,
               num_runs = 1
               )

  ))


year_mpa = (sim_years - burn_years) + 4

slice_grid$sim[[1]] %>%
  group_by(year) %>%
  summarise(
    ssb = sum(ssb),
    catch = sum(biomass_caught),
    mpa = mean(mpa)
  ) %>%
  ungroup() %>%
  mutate(depletion = ssb / ssb[1]) %>%
  ggplot(aes(year, depletion, color = mpa)) +
  geom_point()

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



sim_grid_plots <- slice_grid %>%
  mutate(sim_plot = map(mpa_experiment, ~plot_foo(.x, year_mpa = year_mpa))) %>%
  select(-fish,-fleet,-mpa_experiment)

# a <- slice_grid %>%
#   mutate(run = 1:nrow(.)) %>%
#   select(run,mpa_experiment) %>%
#   unnest() %>%
#   group_by(run,year) %>%
#   mutate(post_mpa = any(percent_mpa > 0))

# trelliscope(sim_grid_plots, name = 'a')

save.image(file.path(run_dir,'mpa-sims.Rdata'))
