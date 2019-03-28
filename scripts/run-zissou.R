# Run Zissou -------
# Author: Dan Ovando
# Project: zissou
# Summary: Run a revised version of zissou, redoing some of the database filtering and
# prep, and comparing raw and standardized approaches, making it easier to run a bunch of
# different model configs all at once. With end goal: accepting or rejecting null hypotheses,
# inclurind zero, but more importantly including MPA literature motivated outcomes.


# setup ---------------------------------------------------------------
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



# options -----------------------------------------------------------------

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

# load data ---------------------------------------------------------------

length_data <- read.csv(glue::glue("{data_dir}/UCSB_FISH.csv")) %>%
  magrittr::set_colnames(., tolower(colnames(.))) %>%
  mutate(classcode = tolower(classcode)) %>%
  mutate(observer = as.character(observer)) %>%
  mutate(observer = ifelse(observer == '', 'unknown', observer))

# Filter data per operations in Fish size biomass processing CIMPA.sas file

length_data <- length_data %>%
  filter(
    level != 'CAN',
    campus == 'UCSB',
    method %in%  c(
      'SBTL_FISH'
    ),!(site == 'SCI_PELICAN' & side == 'FAR WEST'),!(toupper(classcode) %in% c('NO_ORG', 'LDAL', 'CNIC'))
  )

yoy_foo <- function(classcode, fish_tl) {
  new_classcode <- classcode
  if (fish_tl <= 5 &
      is.na(fish_tl) == F &
      classcode %in% c('cpun', 'bfre', 'ocal', 'hsem')) {
    new_classcode <- glue::glue('{classcode}_yoy')

  } # close less than 5cm

  if (fish_tl <= 10 & is.na(fish_tl) == F) {
    if (toupper(classcode) %in% c('PCLA', 'SMYS', 'SPAU', 'SPIN', 'SPUL')) {
      new_classcode <- glue::glue('{classcode}_yoy')
    }

    if (toupper(classcode)  %in% c('SATR', 'SCAR', 'SCAU', 'SCHR', 'GBY')) {
      new_clascode = 'kgb'
    }

    if (toupper(classcode) %in% c('SFLA', 'SSER', 'SMEL', 'OYT')) {
      new_classcode = 'oyb'
    }

    if (toupper(classcode) %in% c('SMIN',
                                  'SNEB',
                                  'SRAS',
                                  'STRE',
                                  'SSAX',
                                  'SDAL',
                                  'SDIP',
                                  'SEBSPP')) {
      new_classcode = 'r_yoy'
    }


  } #close less than than 10cm

  return(new_classcode)

} # close yoy_foo


length_data <- length_data %>%
  mutate(classcode = map2_chr(classcode, fish_tl, yoy_foo))

# flip situations where min_tl is greater than max_tl

flipped <- length_data$max_tl < length_data$min_tl & !is.na(length_data$max_tl) & !is.na(length_data$min_tl)

length_data[flipped,c("min_tl","max_tl")] <-  length_data[flipped,c("max_tl","min_tl")]


# add in covariates -------------------------------------------------------

life_history_data <-
  read_csv(here::here(
    data_dir,
    'VRG Fish Life History in MPA_04_08_11_12 11-Mar-2014.csv'
  )) %>%
  rename(classcode = pisco_classcode) %>%
  mutate(classcode = tolower(classcode)) %>%
  magrittr::set_colnames(., tolower(colnames(.)))

fish_life <-
  life_history_data$taxa %>% str_split(' ', simplify = T) %>%
  as_data_frame() %>%
  select(1:2) %>%
  set_names(c('genus', 'species'))

get_fish_life <- function(genus, species) {
  Predict = Plot_taxa(
    Search_species(Genus = genus, Species = species)$match_taxonomy,
    mfrow = c(2, 2),
    partial_match = T,
    verbose = F
  )
  out <- Predict[[1]]$Mean_pred %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()

  out[colnames(out) != 'Temperature'] <-
    exp(out[colnames(out) != 'Temperature'])

  return(out)

}

fish_life <- fish_life %>%
  mutate(life_traits = map2(genus, species, safely(get_fish_life)))

fish_life <- fish_life %>%
  mutate(fish_life_worked = map(life_traits, 'error') %>% map_lgl(is.null)) %>%
  filter(fish_life_worked) %>%
  mutate(life_traits = map(life_traits, 'result')) %>%
  unnest() %>%
  mutate(taxa = glue::glue('{genus} {species}')) %>%
  set_names(tolower)

life_history_data <- life_history_data %>%
  left_join(fish_life, by = 'taxa')

# check age at recruitment to survey

smallest_length_seen <- length_data %>%
  group_by(classcode) %>%
  summarise(smallest_length_seen = min(fish_tl, na.rm = T))

life_history_data <- life_history_data %>%
  left_join(smallest_length_seen, by = 'classcode') %>%
  mutate(first_age_seen = log(1 - smallest_length_seen / loo) / -k)

#Convert to fished species things with CDFW catches

caselle_fished_species <-
  read_csv(file = here::here(data_dir, 'caselle-2015-fished-species-list.csv')) %>%
  mutate(commonname = tolower(commonname),
         mod_commonname = tolower(mod_commonname)) %>%
  select(mod_commonname, caselle_targeted)

# make sure targeted list matches that from Caselle 2015
life_history_data <- life_history_data %>%
  left_join(caselle_fished_species, by = c('commonname' = 'mod_commonname')) %>%
  mutate(targeted = ifelse(is.na(caselle_targeted), targeted, caselle_targeted))


cdfw_catches <-
  read_csv(file = here::here("processed_data", 'cdfw-catches.csv')) %>%
  group_by(sci_name, year) %>%
  summarise(catch = sum(pounds_caught, na.rm = T)) %>%
  left_join(
    life_history_data %>% select(taxa, classcode) %>% mutate(taxa = tolower(taxa)),
    by = c("sci_name" = "taxa")
  ) %>%
  filter(!is.na(classcode))




targeting_rank <- cdfw_catches %>%
  group_by(sci_name) %>%
  summarise(total_catch = sum(catch)) %>%
  arrange(desc(total_catch)) %>%
  ungroup() %>%
  mutate(catch_rank = percent_rank(total_catch))

fished_species <-
  data_frame(classcode = unique(cdfw_catches$classcode),
             fished = 1)

c99 <- cdfw_catches %>%
  ungroup() %>%
  filter(year == min(year)) %>%
  mutate(year = 1999)

lag_catches <-  cdfw_catches %>%
  bind_rows(c99) %>%
  group_by(classcode) %>%
  arrange(year) %>%
  mutate(year = year + 1,
         lag_catch = catch) %>%
  select(-catch)

if (use_cfdw_fished == T) {
  life_history_data <- life_history_data %>%
    left_join(fished_species, by = 'classcode') %>%
    mutate(targeted = ifelse(fished == 1 &
                               is.na(fished) == F,
                             'Targeted',
                             targeted)) %>%
    select(-fished)
}

life_history_data$targeted <-
  as.numeric(life_history_data$targeted == 'Targeted')

if (rank_targeting == T) {
  life_history_data <- life_history_data %>%
    left_join(targeting_rank, by = 'classcode') %>%
    mutate(targeted = ifelse(targeted == 1 &
                               (is.na(catch_rank) == F), catch_rank, targeted)) %>%
    mutate(targeted = ifelse(is.na(catch_rank) &
                               targeted == 1, 0.5, targeted))

}


site_data <- read_csv(here::here(data_dir,'Final_Site_Table_UCSB.csv')) %>%
  magrittr::set_colnames(., tolower(colnames(.))) %>%
  select(
    site,
    side,
    mpagroup,
    mpa_status,
    reserve,
    region,
    year_mpa,
    mpaareanm2,
    lat_wgs84,
    lon_wgs84
  ) %>%
  rename(lat_wgs84 = lon_wgs84,
         lon_wgs84 = lat_wgs84) %>%
  unique() %>%
  mutate(eventual_mpa = (year_mpa > 0))

# create geographic clustering of species

sightings <- length_data %>%
  filter(is.na(count) == F) %>%
  select(classcode, site, side) %>%
  left_join(site_data, by = c('site', 'side')) %>%
  filter(is.na(lat_wgs84) == F & is.na(lon_wgs84) == F,
         region %in% c('ANA', 'SCI', 'SMI', 'SRI')) %>%
  group_by(classcode) %>%
  summarise(
    mean_lat = mean(lat_wgs84),
    mean_long = mean(lon_wgs84),
    min_lat = min(lat_wgs84),
    max_lat = max(lat_wgs84),
    min_long = min(lon_wgs84),
    max_long = max(lon_wgs84)
  ) %>%
  ungroup()


num_clusters <- data_frame(clusters = 1:20) %>%
  mutate(within_ss = map_dbl(clusters, ~ sum(
    kmeans(
      sightings %>% select(contains('_')),
      centers = .x,
      nstart = 25,
      iter.max = 1000
    )$withinss
  )))

cluster_classcodes <-
  kmeans(
    sightings %>% select(contains('_')),
    centers = 5,
    nstart = 25,
    iter.max = 1000
  )

sightings <- sightings %>%
  mutate(geographic_cluster = cluster_classcodes$cluster)

# geographic_cluster_plot <-
#   ggmap::qmplot(mean_long,
#                 mean_lat ,
#                 color = factor(geographic_cluster),
#                 data = sightings) + theme_classic()

life_history_data <- life_history_data %>%
  left_join(sightings %>% select(classcode, geographic_cluster),
            life_history_data,
            by = 'classcode')

# add life history data into length data
length_data <- length_data %>%
  left_join(life_history_data, by = 'classcode')


# load kelp data

kelp_conn <- ncdf4::nc_open(here::here("data/LandsatKelpBiomass_2017.nc"))

kelp_vars <- c('lat','lon','year','quarter', 'biomass')

kelp <- ncdf4::ncvar_get(kelp_conn, varid = c("biomass")) %>%
  as_data_frame() %>%
  mutate(lat = ncdf4::ncvar_get(kelp_conn, varid = c("lat")),
         lon = ncdf4::ncvar_get(kelp_conn, varid = c("lon"))) %>%
  gather(tdex, kelp_biomass, -lat,-lon) %>%
  mutate(tdex = str_replace_all(tdex,"\\D","") %>% as.numeric()) %>%
  na.omit()

time <-
  data_frame(
    year = ncdf4::ncvar_get(
      kelp_conn,
      varid = c("year")),
      quarter = ncdf4::ncvar_get(kelp_conn, varid = c("quarter"))) %>%
  mutate(tdex = 1:nrow(.))

kelp <- kelp %>%
  left_join(time, by = "tdex") %>%
  select(-tdex) %>%
  mutate(rounded_lat = plyr::round_any(lat, 0.1),
         rounded_lon = plyr::round_any(lon, 0.1)) %>%
  group_by(rounded_lat, rounded_lon, year, quarter) %>%
  summarise(mean_kelp = mean(kelp_biomass, na.rm = T))


kelp_recipe <- recipes::recipe(mean_kelp ~ ., data = kelp)


## process kfm dat

kfm_data <-
  read_csv(here::here(data_dir,'kfm_data/SBCMBON_integrated_fish_20170520.csv'))

kfm_locations <-
  read_csv(here::here(data_dir,'kfm_data/SBCMBON_site_geolocation_20170520.csv'))

conditions_data <- length_data %>%
  group_by(site, side, classcode, year) %>%
  summarise(
    mean_temp = mean(temp, na.rm = T),
    mean_kelp = mean(pctcnpy, na.rm = T),
    mean_vis = mean(vis, na.rm = T)
  ) %>%
  group_by(year) %>%
  mutate(
    mean_temp = ifelse(is.na(mean_temp), mean(mean_temp, na.rm = T), mean_temp),
    mean_kelp = ifelse(is.na(mean_kelp), mean(mean_kelp, na.rm = T), mean_kelp)
  ) %>%
  ungroup() %>%
  mutate(
    mean_temp = ifelse(is.na(mean_temp), mean(mean_temp, na.rm = T), mean_temp),
    mean_kelp = ifelse(is.na(mean_kelp), mean(mean_kelp, na.rm = T), mean_kelp)
  )

observer_experience <- length_data %>%
  group_by(year, month, observer) %>%
  summarise(n_obs = length(side)) %>%
  arrange(observer, year) %>%
  group_by(observer) %>%
  mutate(cumulative_n_obs = cumsum(n_obs),
         lifetime_obs = sum(n_obs)) %>%
  ungroup() %>%
  mutate(observer_ranking = percent_rank(lifetime_obs)) %>%
  mutate(
    trunc_observer = ifelse(
      observer_ranking > 0.5 | observer == 'unknown',
      observer,
      'infrequent'
    )
  ) %>%
  arrange(desc(observer_ranking)) %>%
  mutate(
    cumulative_n_obs = ifelse(
      trunc_observer == 'infrequent' |
        trunc_observer == 'unknown',
      0,
      cumulative_n_obs
    ),
    n_obs = ifelse(
      trunc_observer == 'infrequent' |
        trunc_observer == 'unknown',
      0,
      n_obs
    ),
    lifetime_obs = ifelse(
      trunc_observer == 'infrequent' |
        trunc_observer == 'unknown',
      0,
      lifetime_obs
    ),
    observer_ranking = ifelse(
      trunc_observer == 'infrequent' |
        trunc_observer == 'unknown',
      0,
      observer_ranking
    )
  ) %>%
  mutate(cumulative_n_obs_2 = cumulative_n_obs ^ 2)


length_data <- length_data %>%
  left_join(observer_experience, by = c('year', 'month', 'observer'))


density_data <- read_csv(here::here(data_dir,"ci_reserve_data_final3 txt.csv")) %>%
  magrittr::set_colnames(., tolower(colnames(.))) %>%
  gather('concat.name', 'value', grep('_', colnames(.)), convert = T) %>%
  mutate(
    data.type = gsub('\\_.*', '', concat.name),
    classcode = gsub('.*\\_', '', concat.name)
  ) %>%
  mutate(value = as.numeric(value),
         index = 1:nrow(.)) %>%
  spread(data.type, value) %>%
  rename(site_side = site.side) %>%
  arrange(index) %>%
  select(-index)

site_coords <- density_data %>%
  group_by(site, side) %>%
  summarise(latitude = mean(lon.wgs84, na.rm = T),
            longitude = mean(lat.wgs84, na.rm = T))

species_distributions <- length_data %>%
  left_join(site_coords, by = 'site') %>%
  ungroup() %>%
  filter(is.na(latitude) == F & is.na(longitude) == F) %>%
  group_by(classcode) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  summarise(
    mean_latitude = sum(count * latitude, na.rm = T) / sum(count, na.rm = T),
    mean_longitude = sum(count * longitude, na.rm = T) / sum(count, na.rm = T),
    max_latitude = max(latitude, na.rm = T),
    min_latitude = min(latitude, na.rm = T),
    max_longitude = max(longitude, na.rm = T),
    min_longitude = min(longitude, na.rm = T)
  )

# ggmap::qmplot(min_longitude,min_latitude, data = species_distributions)


# ggmap::qmplot(longitude,latitude, color = side, data = site_coords)


if (file.exists(here::here(data_dir,'enso.csv'))) {
  enso <- read_csv(here(data_dir,'enso.csv'))
} else {

  scrape_enso(outdir = paste0(data_dir,'/'))

  enso <- read_csv(here(data_dir,'enso.csv'))
}

if (file.exists(here(data_dir,'pdo.csv'))) {
  pdo <- read_csv(here(data_dir,'pdo.csv'))

} else {
  scrape_pdo(outdir = paste0(data_dir,'/'))

  pdo <- read_csv(here(data_dir,'pdo.csv'))
}



# convert transect data to density estimates ------------------------------


if (file.exists('processed_data/pisco-data.Rdata') == F |
    run_length_to_density == T) {


  length_example <-   length_data %>%
    filter(is.na(commonname) == F) %>%
    mutate(biomass_g = pmap_dbl(
      list(
        mean_length = fish_tl,
        min_length = min_tl,
        max_length = max_tl,
        count = count,
        weight_a = wl_a,
        weight_b = wl_b,
        length_type_for_weight = wl_input_length,
        length_for_weight_units = wl_l_units,
        weight_units = wl_w_units,
        tl_sl_a = lc.a._for_wl,
        tl_sl_b = lc.b._for_wl,
        tl_sl_type = lc_type_for_wl,
        tl_sl_formula = ll_equation_for_wl
      ),
      length_to_weight
    ))

  pisco_data <- length_example %>%
    mutate(
      observer = ifelse(is.na(observer), 'unknown', observer),
      surge = ifelse(is.na(surge), 'unknown', surge)
    ) %>%
    rename(
      total_biomass_g = biomass_g,
      mean_temp = temp,
      mean_vis = vis,
      mean_depth = depth,
      mean_canopy = pctcnpy
    )

  species_sightings <- pisco_data %>%
    left_join(site_data, by = 'site') %>%
    group_by(site) %>%
    summarise(species_seen = list(unique(classcode)))

  pisco_data <- pisco_data %>%
    ungroup() %>%
    left_join(site_data %>% select(site, region), by = 'site') %>%
    select(region, site, side, year, month, day, zone, level, transect) %>%
    unique() %>%  {
      pmap(
        list(
          this_region = .$region,
          this_site = .$site,
          this_side = .$side,
          this_year = .$year,
          this_month = .$month,
          this_day = .$day,
          this_transect = .$transect,
          this_zone = .$zone,
          this_level = .$level
        ),
        add_missing_fish,
        observations = pisco_data,
        species_sightings = species_sightings,
        life_history_vars = colnames(life_history_data)
      )
    } %>%
    bind_rows()


  classcodes <- unique(pisco_data$classcode)

  life_history_vars <-
    which(colnames(pisco_data) %in% colnames(life_history_data))

  for (i in 1:length(classcodes)) {
    where_class <-
      pisco_data$classcode == classcodes[i] &
      pisco_data$count == 0

    temp_life_history <- life_history_data %>%
      filter(classcode == classcodes[i])

    temp_life_history <-
      temp_life_history[, colnames(pisco_data)[life_history_vars]]

    if (nrow(temp_life_history) > 1) {
      stop('multiple classcodes')
    }

    pisco_data[where_class, life_history_vars] <-
      temp_life_history

    if (any((
      colnames(temp_life_history) == colnames(pisco_data[1, life_history_vars])
    ) == F)) {
      stop()
    }

    if (any(is.na(pisco_data$classcode))) {
      stop()
    }


  }

  check_life_history_fill <- pisco_data %>%
    group_by(classcode) %>%
    summarise(a = n_distinct(commonname)) %>%
    arrange(desc(a))

  if (any(check_life_history_fill$a > 1)) {
    stop('multiple species per classcode')
  }

  save(file = paste0('processed_data/pisco-data.Rdata'),
       pisco_data)

} else {
  load('processed_data/pisco-data.Rdata')

}


# sum biomass across all observed sizes ---------------------------------


transect_covariates <-
  c(
    'mean_depth',
    'mean_vis',
    'mean_temp',
    'surge',
    'mean_canopy',
    'observer',
    'n_obs',
    'cumulative_n_obs'
  )

mean_foo <- function(var) {
  if (class(var) != 'character') {
    out <- mean(var, na.rm = T)
  } else{
    out <- paste(unique(var), collapse = '-')
  }

}


channel_sides <-
  data_frame(region = c("ANA", "SCI", "SRI", "SMI"),
             midline = c(34, 34, 34, 34.03))

pisco_data <- pisco_data %>%
  mutate(year_month = year + (month / 12 - .1))

rolling_mean_temperatures <- pisco_data %>%
  mutate(year_month = year + (month / 12 - .1)) %>%
  left_join(site_data, by = c('site','side')) %>%
  filter(region %in% c("ANA", "SCI","SRI","SMI")) %>%
  group_by(year) %>%
  summarise(mean_temperature = mean(mean_temp)) %>%
  mutate(rolling_mean = RcppRoll::roll_mean(mean_temperature,4, align = "right",
                                            fill = NA)) %>%
  mutate(rolling_mean3 = RcppRoll::roll_mean(mean_temperature,3, align = "right",
                                             fill = NA))  %>%
  mutate(rolling_mean2 = RcppRoll::roll_mean(mean_temperature,2, align = "right",
                                             fill = NA)) %>%
  mutate(rolling_mean1 = RcppRoll::roll_mean(mean_temperature,1, align = "right",
                                             fill = NA)) %>%
  mutate(rolling_mean = ifelse(is.na(rolling_mean), rolling_mean3, rolling_mean)) %>%
  mutate(rolling_mean = ifelse(is.na(rolling_mean), rolling_mean2, rolling_mean)) %>%
  mutate(rolling_mean = ifelse(is.na(rolling_mean), rolling_mean1, rolling_mean)) %>%
  select(year, rolling_mean) %>%
  mutate(rolling_mean = zoo::na.approx(rolling_mean)) %>%
  ungroup()


# sum biomass across all lengths
pisco_data <- pisco_data %>%
  mutate(sampling_event = paste(campus, method, year, month, day, site, side, zone, transect, level) %>%
           as.factor() %>% as.numeric()) %>%
  group_by(sampling_event,
           classcode) %>%
  mutate(
    total_biomass_g = sum(total_biomass_g),
    density_g_m2 = total_biomass_g / 60,
    total_count = sum(count),
    mean_length = mean(fish_tl, na.rm = T),
    counter = 1:length(total_biomass_g)
  ) %>%
  ungroup() %>%
  filter(counter == 1) %>%
  select(-counter)

pisco_data <- pisco_data %>%
  mutate(
    any_seen = total_biomass_g > 0,
    factor_year = factor(year),
    log_density = log(density_g_m2),
    factor_month = factor(month),
    site_side = glue::glue('{site}-{side}')
  ) %>%
  group_by(site, side, month, year) %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_canopy = ifelse(is.na(mean_canopy), mean(mean_canopy, na.rm = T), mean_canopy)
  ) %>%
  group_by(site, side, year) %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_canopy = ifelse(is.na(mean_canopy), mean(mean_canopy, na.rm = T), mean_canopy)
  ) %>%
  group_by(site, year) %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_canopy = ifelse(is.na(mean_canopy), mean(mean_canopy, na.rm = T), mean_canopy)
  ) %>%
  group_by(year) %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_canopy = ifelse(is.na(mean_canopy), mean(mean_canopy, na.rm = T), mean_canopy)
  ) %>%
  ungroup() %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_canopy = ifelse(is.na(mean_canopy), mean(mean_canopy, na.rm = T), mean_canopy)
  ) %>% mutate(
    mean_temp = ifelse(is.na(mean_temp), mean(mean_temp, na.rm = T), mean_temp)) %>%
  mutate(temp_deviation = abs(mean_temp - temperature)) %>%
  mutate(generations_protected = pmin(round((year - year_mpa - 1) / tm), max_generations))


# Process kfm count data

kfm_sites <- kfm_locations %>%
  select(latitude, longitude) %>%
  as.matrix() %>%
  na.omit()

pisco_sites <- site_coords %>%
  ungroup() %>%
  select(latitude, longitude) %>%
  as.matrix()

nearest_pisco_site <- RANN::nn2(pisco_sites, kfm_sites)

nearest_pisco_site <- site_coords[nearest_pisco_site$nn.idx[, 1], ]

kfm_pisco_locations <- kfm_locations %>%
  na.omit() %>%
  bind_cols(nearest_pisco_site)


kfm_data <- kfm_data %>%
  filter(sample_method != "visualfish" &
           sample_method != "crypticfish",
         data_source == 'kfm') %>%
  left_join(life_history_data, by = c('taxon_name' = 'taxa')) %>%
  left_join(kfm_pisco_locations %>% select(-data_source), by = 'site_id') %>%
  filter(str_detect(geolocation, '_island'),
         data_source != 'lter',
         is.na(geolocation) == F) %>%
  mutate(region = str_extract(geolocation, '.*(?=_island)')) %>%
  filter(region %in% c('anacapa', 'santa_cruz', 'santa_rosa', 'san_miguel'))

kfm_data <- kfm_data %>%
  mutate(
    count = as.numeric(count),
    density = count / area,
    year = lubridate::year(date),
    month = lubridate::month(date)
  )

region_table <- tribble(
  ~ region,
  ~ short_region,
  'anacapa',
  'ANA',
  'santa_cruz',
  'SCI',
  'santa_rosa',
  'SRI',
  'san_miguel',
  'SMI'
)

kfm_data <- kfm_data %>%
  left_join(region_table, by = 'region') %>%
  select(-region) %>%
  rename(region = short_region)

kfm_data <- kfm_data %>%
  mutate(
    log_density = log(density),
    any_seen = density > 0,
    factor_year = as.factor(year),
    factor_month = as.factor(month),
    total_biomass_g = density
  ) %>%
  select(-region) #for compatibility with things later on


# save raw-isa data -------------------------------------------------------

save(
  file = glue::glue("{run_dir}/rawish_zissou_data.Rdata"),
  life_history_data,
  pisco_data,
  kfm_data
)

# filter data -------------------------------------------------------------

consistent_sites <- pisco_data %>%
  group_by(site) %>%
  summarise(
    num_years = length(unique(year)),
    min_year = min(year),
    max_year = max(year)
  ) %>%
  arrange(desc(num_years)) %>%
  filter(min_year <= 2000,
         num_years > 10)

consistent_regions <- pisco_data %>%
  left_join(site_data, by = 'site') %>%
  group_by(region) %>%
  summarise(
    num_years = length(unique(year)),
    min_year = min(year),
    max_year = max(year)
  ) %>%
  arrange(desc(num_years)) %>%
  filter(min_year <= 2000,
         num_years > 10)

nobs_quantiles <- pisco_data %>%
  filter(any_seen == T) %>%
  group_by(year, classcode) %>%
  summarise(nobs = sum(any_seen)) %>%
  ungroup() %>%
  {
    quantile(.$nobs)
  }

well_observed_species <- pisco_data %>%
  filter(year > 1999) %>%
  group_by(year, classcode, commonname, targeted) %>%
  summarise(nseen = sum(any_seen, na.rm = T)) %>%
  group_by(commonname, classcode, targeted) %>%
  summarise(min_seen = min(nseen, na.rm = T)) %>%
  arrange(desc(min_seen)) %>%
  ungroup() %>%
  filter(min_seen > 2) %>%
  # filter(min_seen > nobs_quantiles[3]) %>%
  mutate(classcode = (classcode))


filterfoo <-
  function(x,
           y,
           min_seen_years = 14,
           mpa_start_year = 2003,
           min_year,
           filter_level) {
    # only years above 1999 at the main channel islands
    x <- x %>%
      mutate(classcode = y) %>%
      filter(is.na(log_density) == F) %>%
      filter(year > min_year &
               region %in% c('ANA', 'SCI', 'SMI', 'SRI')) %>%
      group_by(!!filter_level) %>%
      mutate(
        num_years_observed = length(unique(year[any_seen == T])),
        min_year = min(year[any_seen == T], na.rm = T),
        max_year = max(year[any_seen == T], na.rm = T)
      ) %>%
      ungroup() %>%
      filter(
        num_years_observed >= min_seen_years,
        min_year <= mpa_start_year,
        max_year >= mpa_start_year
      ) %>%
      select(-classcode)# actually observed for at least 10 years
  }


# create dropped pairs


site_data$rounded_lat <- site_data$lat_wgs84

site_data$rounded_lon <- site_data$lon_wgs84

abundance_data <- pisco_data %>%
  left_join(site_data, by = c("site","side")) %>%
  left_join(lag_catches %>% select(year, classcode, lag_catch), by = c('year', 'classcode')) %>%
  mutate(data_source = 'pisco') %>%
  nest(-data_source,-classcode) %>%
  bind_rows(kfm_data %>%  mutate(data_source = 'kfm') %>%
              left_join(site_data, by = c("site","side")) %>%
              left_join(lag_catches %>% select(year, classcode, lag_catch), by = c('year', 'classcode')) %>%
              mutate(year_month = year + (month / 12 - .1),
                     site_side = site_id) %>%
              nest(-data_source,-classcode)
              ) %>%
  filter(classcode %in% well_observed_species$classcode) %>%
  # mutate(data = map(data, ~ left_join(.x, site_data, by = c('site', 'side')))) %>%
  # mutate(data = map(data, ~ left_join(.x, channel_sides, by = c("region")))) %>%
  # mutate(data = map(data, ~ mutate(.x, channel_zone = ifelse(lat_wgs84 > midline, "inner", "outer")))) %>%
  mutate(data = map(data, ~ left_join(.x, rolling_mean_temperatures, by = c("year")))) %>%
  mutate(data = map(data, ~ mutate(.x, temp_deviation = abs(rolling_mean - temperature)))) %>%
  mutate(data = map(data, ~ mutate(.x, month = as.numeric(as.character(factor_month))))) %>%
  mutate(data = map(data, ~ left_join(.x, enso, by = c('year', 'month')))) %>%
  mutate(data = map(data, ~ left_join(.x, pdo, by = c('year', 'month')))) %>%
  mutate(data = map(data, ~ mutate(.x, lag_catch = ifelse(is.na(lag_catch), 0, lag_catch)))) %>%
  mutate(
    data = map2(
      data,
      classcode,
      filterfoo,
      min_year = min_year,
      min_seen_years = 14,
      filter_level = quo(classcode)
    )
  ) %>% # filter out things
  mutate(dim_data = map_dbl(data, nrow)) %>%
  filter(dim_data > 0) %>%
  left_join(life_history_data %>% select(classcode, commonname, targeted),
            by = 'classcode') %>%
  filter(
    str_detect(commonname, 'YOY') == F,
    is.na(targeted) == F,
    str_detect(classcode, '_yoy') == F
  )

abundance_data <- abundance_data %>%
  select(classcode, data_source, data) %>%
  unnest() %>%
  nest(-data_source)


# add kelp data

kelp_model <- caret::knnreg(mean_kelp ~ ., data = kelp, k = 5)

abundance_data <- abundance_data %>%
  mutate(data = map2(data,list(kelp_model),  ~ mutate(.x, interp_kelp = predict(
    .y,
    newdata = .x %>% mutate(quarter = (month/3) %>% ceiling()) %>%  select(year, quarter, rounded_lat, rounded_lon)
  ))))

save(
  file = glue::glue("{run_dir}/abundance_data.Rdata"),
  abundance_data
)


# run vast ----------------------------------------------------------------

# if (run_vast == T) {
#   vast_data <- abundance_data %>%
#     filter(
#       data_source == 'pisco'
#     ) %>%
#     unnest() %>%
#     nest(-data_source, -classcode) %>%
#     mutate(survey_region = 'california_current',
#            n_x = num_knots) %>%
#     mutate(
#       vast_data = map2(
#         data,
#         classcode,
#         vast_prep,
#         site_coords = site_coords,
#         conditions_data = conditions_data
#       )
#     )
#
#
#   vast_abundance <- vast_data %>%
#     mutate(
#       vast_results = purrr::pmap(
#         list(
#           region = survey_region,
#           raw_data = vast_data,
#           n_x = n_x
#         ),
#         safely(vasterize_pisco_data),
#         run_dir = run_dir,
#         nstart = 100,
#         obs_model = c(2, 0),
#         catchability_variables_names = c('mean_vis', 'surge', 'factor_month', 'zone', 'level'),
#         vessel = 0,
#         vessel_year = 1
#       )
#     )
#
#   vast_abundance <- vast_abundance %>%
#     mutate(vast_error = map(vast_results, 'error'))
#
#   vast_abundance <- vast_abundance %>%
#     mutate(vast_index = map(vast_results, 'result')) %>%
#     mutate(no_vast_error = map_lgl(vast_error, is.null)) %>%
#     filter(no_vast_error)
#
#   save(file = paste0(run_dir, '/vast_abundance.Rdata'), vast_abundance)
# } else {
#   load(file = paste0(run_dir, '/vast_abundance.Rdata'))
# }

# fit model ---------------------------------------------------------------


abundance_data <- abundance_data %>%
  mutate(data = map(data, ~select(.,log_density,
                                  targeted,
                                  site_side,
                                  region,
                                  geographic_cluster,
                                  level,
                                  factor_month,
                                  cumulative_n_obs,
                                  surge,
                                  mean_canopy,
                                  mean_vis,
                                  mean_depth,
                                  lag_catch,
                                  factor_year,
                                  classcode,
                                  interp_kelp,
                                  eventual_mpa,
                                  temp_deviation,
                                  any_seen,
                                  year
  )))


script_name <- "fit_zissou"



var_options <-  data_frame(var_names = c("pisco_a","kfm"),
                           non_nested_variables =
                             list(c(
                               "site_side",
                               'level',
                               'factor_month',
                               'cumulative_n_obs',
                               'surge',
                               'mean_depth',
                               'mean_vis'
                             ),
                             c(
                               "site_side",
                               'factor_month'
                             )
                             ))


stupid_filter <- function(x){

  missing = x %>%
    group_by(classcode) %>%
    summarise(nyd = n_distinct(year)) %>%
    ungroup() %>%
    filter(nyd < 17)

  out <- x %>%
    filter(!(classcode %in% missing$classcode))

}

abundance_data <- abundance_data %>%
  mutate(data = map(data, stupid_filter))

model_runs <- cross_df(
  list(
    data_source = abundance_data$data_source,
    var_names = var_options$var_names,
    mpa_only = c(TRUE,FALSE),
    center_scale = c(TRUE)
   )
  ) %>%
  left_join(var_options, by = "var_names") %>%
  left_join(abundance_data, by = "data_source") %>%
  filter(!(data_source == "kfm" & str_detect(var_names,"pisco"))) %>%
  filter(!(data_source == "kfm" & mpa_only == TRUE))

# model_runs <- model_runs %>%
#   filter((data_source == "pisco" & mpa_only == FALSE & center_scale == TRUE) |
#            data_source == "kfm" & center_scale == TRUE)

# model_runs <- model_runs %>%
#   filter(data_source == "pisco" & mpa_only == FALSE & center_scale == TRUE)

if (run_tmb == T){

  if (file.exists("fit-progress.txt")){
  file.remove("fit-progress.txt")
  }

  # future::plan(future::multiprocess, workers = 4)

  doParallel::registerDoParallel(cores = n_cores)
  #
  fits <- foreach::foreach(i = 1:nrow(model_runs)) %dopar% {

  # fits <- list()

  # for (i in 1:nrow(model_runs)){

    sfz <- safely(fit_zissou)

    fits <- sfz(data = model_runs$data[[i]],
                  non_nested_variables = model_runs$non_nested_variables[[i]],
                  mpa_only = model_runs$mpa_only[[i]],
                  center_scale = model_runs$center_scale[[i]],
                  run_dir = run_dir,
                  script_name = script_name,
                  fixed_regions = FALSE,
                  include_intercept = TRUE,
                  fixed_did = FALSE,
                non_nested_did_variables = c(
                  "temp", "kelp", "lag_catch"
                )
    )

    write(glue::glue("{round(100*i/nrow(model_runs),2)}% done with model fits"), file = "fit-progress.txt",
          append = T)

    out <- fits

    } # close dopar


save(file = paste0(run_dir, '/model_fits.Rdata'),
       fits)

model_runs$tmb_fit <- fits

save(file = paste0(run_dir, '/model_runs.Rdata'),
     model_runs)


} else {

  load(file = file.path(run_dir, 'model_runs.Rdata'))
}


models_worked <- model_runs$tmb_fit %>% map("error") %>% map_lgl(is_null)

model_runs$tmb_fit %>% map_dbl(length)

model_runs <- model_runs %>%
  filter(models_worked) %>%
  mutate(tmb_fit = map(tmb_fit,"result")) %>%
  mutate(processed_fits = map(tmb_fit, process_fits)) %>%
  mutate(did_plot = map(processed_fits, "did_plot"))




