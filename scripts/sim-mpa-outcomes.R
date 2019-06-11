



# setup -------------------------------------------------------------------

library(spasm)
library(FishLife)
library(hrbrthemes)
library(scales)
library(doParallel)
library(furrr)
library(tidyverse)

functions <- list.files(here::here("functions"))

walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions


# options -----------------------------------------------------------------

sim_years <- 50

burn_years <- 25

num_patches <- 50

run_experiments <- TRUE

save_experiment <- TRUE

create_grid <- TRUE

n_cores <- 6

samps <- 100

grid_search <-  FALSE

in_clouds <- F

if (in_clouds == T) {
  system("umount results/zissou-results")

  system("rm -r results/zissou-results")

  if (dir.exists("results/zissou-results") == F) {
    system("mkdir results/zissou-results")

  }

  system("gcsfuse zissou-results results/zissou-results")

  system("umount data/zissou-data")

  system("rm -r data/zissou-data")

  if (dir.exists("results/zissou-data") == F) {
    system("mkdir data/zissou-data")

  }

  # system("mkdir data/scrooge-data")

  system("gcsfuse zissou-data data/zissou-data")

  cloud_dir <- here::here("results", "zissou-results", run_name)

  if (dir.exists(cloud_dir) == F) {
    dir.create(cloud_dir)

  }

}


# prepare data -----------------------------------------------------

run_name <- "v4.1"

run_dir <- here::here("results", run_name)

experiment_dir <- file.path(run_dir, "experiments")

load(file = file.path(run_dir, "rawish_zissou_data.Rdata"))

load(file = file.path(run_dir, "model_runs.Rdata"))

load(file = file.path(run_dir, "abundance_data.Rdata"))


fitted_data <- model_runs$data[[1]]

seen_species <- life_history_data %>%
  filter(classcode %in% (fitted_data$classcode %>% unique())) %>%
  rename(sci_name = taxa,
         linf = vbgf.linf,
         common_name = commonname) %>%
  mutate(sci_name = tolower(sci_name)) %>%
  filter(is.na(linf) == F)

fitted_data <- fitted_data %>%
  left_join(life_history_data %>% select(classcode, taxa) %>% unique(), by = "classcode")

# prepare experiments -----------------------------------------------------


# run experiments ---------------------------------------------------------

if (run_experiments == T) {
  if (create_grid == TRUE) {
    fun <- function() {
      ANSWER <-
        readline("STOP!!! Are you sure you want to overwrite the experiment grid? y/n ")
      if (substr(ANSWER, 1, 1) == "y") {
        cat("OK, sit back, this might take a while")

        fun <- function() {
          ANSWER <-
            readline("Delete old experiments? Recommend if tuning grid again")
          if (substr(ANSWER, 1, 1) == "y") {
            unlink(here::here("results", run_name, "experiments"),
                   recursive = TRUE)

          } else{

          }
        }

      } else{
        cat("Probably a good idea, stopping")
        stop()

      }
    }
    fun()

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
        fleet_model = c("constant-catch"),
        effort_allocation = c("profit-gravity", "simple"),
        stringsAsFactors = F
      )
    } else{
      sim_grid <-
        tibble(
          scientific_name = sample(unique(fitted_data$taxa), samps, replace = T),
          steepness = runif(samps, min = 0.6, max = 0.95),
          adult_movement = sample(0:(0.5 * num_patches), samps, replace = T),
          larval_movement = sample(0:(0.5 * num_patches), samps, replace = T),
          density_movement_modifier = sample(c(0, 0.25), samps, replace = T),
          density_dependence_form = sample(1:3, samps, replace = T),
          mpa_size = runif(samps, min = 0.01, max = 1),
          f_v_m = runif(samps, min = 0.01, max = 4),
          fleet_model = sample(
            c("open-access", "constant-effort", "constant-catch"),
            samps,
            replace = T
          ),
          effort_allocation = sample(
            c("profit-gravity", "simple", "gravity"),
            samps,
            replace = T
          ),
          year_mpa = sample(5:(sim_years / 1.5), samps, replace = T),
          sprinkler = sample(c(TRUE, FALSE), samps, replace = TRUE),
          mpa_reaction   =  sample(c("stay", "leave"), samps, replace = TRUE),
          min_size = runif(samps, min = 0.01, max = 0.75),
          mpa_habfactor = sample(c(1, 4), samps, replace = TRUE),
          size_limit = runif(samps, 0.1, 1.25),
          random_mpa = sample(c(TRUE, FALSE), samps, replace = TRUE)
        )


      # sim_grid <-
      #   tibble(
      #     scientific_name = sample(unique(fitted_data$taxa), samps, replace = T),
      #     steepness = runif(samps, min = 0.6, max = 1),
      #     adult_movement = sample(0:num_patches, samps, replace = T),
      #     larval_movement = sample(0:num_patches, samps, replace = T),
      #     density_movement_modifier = sample(c(0, 0.5), samps, replace = T),
      #     density_dependence_form = sample(1:3, samps, replace = T),
      #     mpa_size = runif(samps, min = 0.05, max = 1),
      #     f_v_m = runif(samps, min = 0.01, max = 2),
      #     fleet_model = sample(c("open-access","constant-effort","constant-catch"), samps, replace = T),
      #     effort_allocation = sample(c("profit-gravity", "simple","gravity"), samps, replace = T),
      #     year_mpa = sample(5:(sim_years/2), samps, replace = T),
      #     sprinkler = sample(c(TRUE,FALSE), samps, replace = TRUE),
      #     mpa_reaction   =  sample(c("stay","leave"), samps, replace = TRUE),
      #     min_size = runif(samps, min = 0.01, max = 1),
      #     mpa_habfactor = sample(c(1,4), samps, replace = TRUE)
      # )

    }

    # create fish objects
    sim_grid <- sim_grid %>%
      mutate(fish = pmap(
        list(
          scientific_name = scientific_name,
          steepness = steepness,
          adult_movement = adult_movement,
          larval_movement = larval_movement,
          density_dependence_form = density_dependence_form,
          density_movement_modifier = density_movement_modifier
        ),
        create_fish,
        price = 10
      ))


    # create fleet objects
    sim_grid <- sim_grid %>%
      mutate(fleet = pmap(
        list(
          fish = fish,
          fleet_model = fleet_model,
          effort_allocation = effort_allocation,
          mpa_reaction = mpa_reaction,
          length_50_sel = size_limit * map_dbl(sim_grid$fish, "length_50_mature")
        ),
        create_fleet,
        q = .1
      ))

    # tune fleet objects

    future::plan(future::multiprocess, workers = n_cores)

    message("starting fishery tuning")
    sim_grid <- sim_grid %>%
      mutate(
        tuned_fishery = future_pmap(
          list(
            f_v_m = f_v_m,
            fish = fish,
            fleet = fleet,
            sprinkler = sprinkler,
            mpa_habfactor = mpa_habfactor
          ),
          safely(tune_fishery),
          num_patches = num_patches,
          sim_years = sim_years,
          burn_years = burn_years,
          .progress = T
        )
      )
    # mutate(tune_worked = map_lgl(map(tune_fishery,"error"), is.null)) %>%
    # filter(tune_worked) %>%
    # mutate(tune_fishery = map(tune_fishery,"result"))


    save(sim_grid, file = paste0(run_dir, "/sim_grid.Rdata"))
    message("finished fishery tuning")

  } else{
    load(file = paste0(run_dir, "/sim_grid.Rdata"))

  }

  tuning_worked <-
    map(sim_grid$tuned_fishery, "error") %>% map_lgl(is.null)

  sim_grid <- sim_grid %>%
    filter(tuning_worked) %>%
    mutate(fish = map(tuned_fishery, c("result", "fish")),
           fleet = map(tuned_fishery, c("result", "fleet")))

  sim_grid$tuned_fishery <- map(sim_grid$tuned_fishery, "result")

  doParallel::registerDoParallel(cores = n_cores)

  foreach::getDoParWorkers()

  sim_grid$experiment <- 1:nrow(sim_grid)

  if (dir.exists(experiment_dir) == F) {
    dir.create(experiment_dir, recursive = T)
  }

  # library(progress)

  message("starting mpa experiments")

  # sim_grid <- sample_n(sim_grid, 100)

  mpa_experiments <-
    foreach::foreach(i = 1:nrow(sim_grid)) %dopar% {
      # pb$tick()
      results <- sim_grid %>%
        slice(i) %>%
        mutate(
          mpa_experiment = pmap(
            list(
              fish = fish,
              fleet = fleet,
              mpa_size = mpa_size,
              year_mpa = year_mpa,
              sprinkler = sprinkler,
              mpa_habfactor = mpa_habfactor,
              min_size = min_size,
              random_mpa = random_mpa
            ),
            run_mpa_experiment,
            sim_years = sim_years,
            burn_years = burn_years,
            num_patches = num_patches
          )
        )

      # results$mpa_experiment[[1]]$raw_outcomes %>% filter(year == max(year)) -> a
      #
      # a %>% filter(experiment == "with-mpa") %>% group_by(patch) %>% summarise(mpa = unique(mpa)) %>% ggplot(aes(patch, mpa)) + geom_col()

      filename <- glue::glue("experiment_{i}.rds")
      #
      # out <- results

      if (save_experiment == TRUE) {
        saveRDS(results, file = glue::glue("{experiment_dir}/{filename}"))

      } else {
        out <- results

      }
    } # close dopar

  # mpa_experiments[[1]]$mpa_experiment[[1]]$raw_outcomes %>% filter(year == max(year)) -> a
  #
  # a %>% filter(experiment == "with-mpa") %>% group_by(patch) %>% summarise(mpa = unique(mpa)) %>% ggplot(aes(patch, mpa)) + geom_col()
  #
  # filename <- glue::glue("experiment_{i}.rds")

  # out <- results


} # close run experiments

message("finished mpa experiments")

# process outcomes --------------------------------------------------------


loadfoo <-
  function(experiment, experiment_dir, output = "mpa-effect") {
    ex <-
      readRDS(glue::glue("{experiment_dir}/experiment_{experiment}.rds"))

    # if (output == "mpa-effect") {
    ex$msy <- ex$tuned_fishery[[1]]$fish$msy

    ex$b_msy <- ex$tuned_fishery[[1]]$fish$b_msy$b_msy

    ex <- purrr::map_df(list(ex), study_mpa)

    ex <- ex %>%
      select(-mpa_experiment, -fish,-fleet,-tuned_fishery)

    return(ex)
  }


# if (logged == TRUE) {
future::plan(future::multiprocess, workers = n_cores)

processed_grid <-
  future_map(1:samps,
             safely(loadfoo),
             experiment_dir = experiment_dir,
             .progress = T)

grid_worked <- map(processed_grid, "error") %>% map_lgl(is_null)

processed_grid <- processed_grid %>%
  keep(grid_worked)

processed_grid <- map(processed_grid, "result") %>%
  bind_rows()


save(processed_grid, file = paste0(run_dir, "/processed_grid.Rdata"))

outcomes <- processed_grid %>%
  unnest()

outcomes %>%
  ggplot(aes(mpa_effect, fill = factor(density_movement_modifier))) +
  geom_histogram() +
  facet_wrap(~factor(density_movement_modifier))




# save outcomes -----------------------------------------------------------
