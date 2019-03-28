# zissou controlfile. Use to set global options for all runs

set.seed(42)
library(scales)
library(viridis)
library(ggmap)
library(forcats)
library(stringr)
library(lubridate)
library(lme4)
library(TMB)
library(FishLife)
library(patchwork)
library(rstan)
library(extrafont)
library(hrbrthemes)
library(caret)
library(here)
library(doParallel)
library(tidyverse)
library(recipes)
# library(furrr)

functions <- list.files(here::here("functions"))

walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions

run_name <- 'v4.0'

run_description <- "post defense improvements and author feedback. Ideally publication version"

in_clouds <- F

if (in_clouds == F){

  run_dir <- here::here("results", run_name)

} else{

  run_dir <- here::here("results","zissou-results", run_name)

}

if (dir.exists(run_dir) == F) {
  dir.create(run_dir, recursive = T)
}


if (in_clouds == T){


  system("umount results/zissou-results")

  system("rm -r results/zissou-results")

  if (dir.exists("results/zissou-results") == F){

    system("mkdir results/zissou-results")

  }

  system("gcsfuse zissou-results results/zissou-results")

  system("umount data/zissou-data")

  system("rm -r data/zissou-data")

  if (dir.exists("results/zissou-data") == F){

    system("mkdir data/zissou-data")

  }

  # system("mkdir data/scrooge-data")

  system("gcsfuse zissou-data data/zissou-data")

  cloud_dir <- here::here("results","zissou-results",run_name)

  if (dir.exists(cloud_dir) == F){

    dir.create(cloud_dir)

  }

}

if (in_clouds == F ){

  data_dir <- "data"
} else{

  data_dir <- "data/zissou-data"
}

write(run_description,
      file = paste(run_dir, 'RUN_DESCRIPTION.txt', sep = '/'))


# fit difference in difference options -----------------------------------------------------------------

rstan_options(auto_write = TRUE)

run_tmb <- FALSE

n_cores <- 5

tmb_to_stan <- FALSE # fit the model in stan instead of TMB

max_generations <- 4

run_length_to_density <-  FALSE

run_vast <- FALSE # run VAST, best to leave off for now

num_knots <-  10

channel_islands_only <- T # only include channel islands, leave T

min_year <- 1999 # included years must be greater than this

mpa_only <- FALSE

occurance_ranking_cutoff <- 0.5

small_num <-  0 # no idea

use_mpa_site_effects <- F # no idea

use_cfdw_fished <-  F

year_mpa <- 2003

rank_targeting <- F

max_generations <- 5

max_year <- 2017

plot_theme <- hrbrthemes::theme_ipsum(base_size = 14,
                                      axis_title_size = 16)


theme_set(plot_theme)
