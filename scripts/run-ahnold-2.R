# Run Ahnold v2 -------
# Author: Dan Ovando
# Project: Ahnold
# Summary: Run a revised version of ahnold, redoing some of the database filtering and
# prep, and comparing raw and standardized approaches, making it easier to run a bunch of
# different model configs all at once. With end goal: accepting or rejecting null hypotheses,
# inclurind zero, but more importantly including MPA literature motivated outcomes.


# setup ---------------------------------------------------------------
rm(list = ls())
library(scales)
# library(rstanarm)
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
library(tidyverse)
demons::load_functions('functions')

run_name <- 'Working'

run_dir <- file.path('results', run_name)

run_description <-
  'Model selection process, testing STAN selection'

if (dir.exists(run_dir) == F) {
  dir.create(run_dir)
}

write(run_description,
      file = paste(run_dir, 'RUN_DESCRIPTION.txt', sep = '/'))


# options -----------------------------------------------------------------

run_length_to_density <-  F

run_vast <- F # run VAST, best to leave off for now

num_knots <-  10

aggregate_transects <- T # should transects be aggregated up (mean across observer) or left as raw

channel_islands_only <- T # only include channel islands, leave T

min_year <- 1999 # included years must be greater than this

occurance_ranking_cutoff <- 0.5

small_num <-  0 # no idea

use_mpa_site_effects <- F # no idea

# load data ---------------------------------------------------------------

length_data <- read_csv('data/UCSB_FISH raw thru 2013.csv') %>%
  magrittr::set_colnames(., tolower(colnames(.))) %>%
  mutate(classcode = tolower(classcode)) %>%
  mutate(observer = ifelse(is.na(observer), 'unknown', observer))

# Filter data per operations in Fish size biomass processing CIMPA.sas file

length_data <- length_data %>%
  filter(level != 'CAN',
         campus == 'UCSB',
         method %in%  c('SBTL_FISH', 'SBTL_FISH_NPS', 'SBTL_FISH_CRANE', 'SBTL_FISH_VRG'),
         !(site == 'SCI_PELICAN' & side == 'FAR WEST'),
         !(toupper(classcode) %in% c('NO_ORG', 'LDAL', 'CNIC'))
  )

yoy_foo <- function(classcode, fish_tl){

  new_classcode <- classcode
  if (fish_tl <= 5 & is.na(fish_tl) == F & classcode %in% c('cpun','bfre','ocal','hsem')){

    new_classcode <- glue::glue('{classcode}_yoy')

  } # close less than 5cm

  if (fish_tl <= 10 & is.na(fish_tl) == F){

    if (toupper(classcode) %in% c('PCLA','SMYS','SPAU','SPIN','SPUL')){
      new_classcode <- glue::glue('{classcode}_yoy')
    }

    if (toupper(classcode)  %in% c('SATR', 'SCAR', 'SCAU', 'SCHR', 'GBY')){new_clascode = 'kgb'}

    if (toupper(classcode)%in% c('SFLA', 'SSER', 'SMEL', 'OYT')){ new_classcode ='oyb'}

    if (toupper(classcode) %in% c('SMIN', 'SNEB', 'SRAS', 'STRE', 'SSAX', 'SDAL', 'SDIP', 'SEBSPP')){
      new_classcode = 'r_yoy'}


  } #close less than than 10cm

  return(new_classcode)

} # close yoy_foo


length_data <- length_data %>%
  mutate(classcode = map2_chr(classcode, fish_tl, yoy_foo))

# add in covariates -------------------------------------------------------

life_history_data <-
  read_csv('data/VRG Fish Life History in MPA_04_08_11_12 11-Mar-2014.csv') %>%
  rename(classcode = pisco_classcode) %>%
  mutate(classcode = tolower(classcode)) %>%
  magrittr::set_colnames(., tolower(colnames(.)))

fish_life <- life_history_data$taxa %>% str_split(' ', simplify = T) %>%
  as_data_frame() %>%
  select(1:2) %>%
  set_names(c('genus','species'))

get_fish_life <- function(genus, species){
  Predict = Plot_taxa( Search_species(Genus=genus,Species=species)$match_taxonomy, mfrow=c(2,2), partial_match = T, verbose = F )

  out <- Predict[[1]]$Mean_pred %>%
    as.matrix() %>%
    as.data.frame() %>%
    mutate(variable = row.names(.)) %>%
    spread(variable, V1)

  out[colnames(out) != 'Temperature'] <- exp(out[colnames(out) != 'Temperature'])

  return(out)

}

fish_life <- fish_life %>%
  mutate(life_traits = map2(genus, species, safely(get_fish_life)))

fish_life <- fish_life %>%
  mutate(fish_life_worked = map(life_traits,'error') %>% map_lgl(is.null)) %>%
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
  mutate(first_age_seen = log(1 - smallest_length_seen/loo)/-k)

#Convert to fished species things with CDFW catches

ci_catches <-
  read_csv(file = file.path('data', 'cfdw-channel-islands-catches.csv')) %>% group_by(classcode, year) %>%
  summarise(catch = sum(pounds_caught, na.rm = T))

fished_species <-
  data_frame(classcode = unique(ci_catches$classcode),
             fished = 1)

life_history_data <- life_history_data %>%
  left_join(fished_species, by = 'classcode') %>%
  mutate(targeted = ifelse(
    fished == 1 &
      is.na(fished) == F,
    'Targeted',
    targeted
  )) %>%
  select(-fished)

# add life history data into length data
length_data <- length_data %>%
  left_join(life_history_data, by = 'classcode')

site_data <- read_csv('data/Final_Site_Table_UCSB.csv') %>%
  magrittr::set_colnames(., tolower(colnames(.))) %>%
  select(site,
         side,
         mpagroup,
         mpa_status,
         reserve,
         region,
         year_mpa,
         mpaareanm2) %>%
  unique() %>%
  mutate(eventual_mpa = (year_mpa > 0))


kfm_data <- read_csv('data/kfm_data/SBCMBON_integrated_fish_20170520.csv')

kfm_locations <- read_csv('data/kfm_data/SBCMBON_site_geolocation_20170520.csv')

conditions_data <- length_data %>%
  group_by(site, side, year) %>%
  summarise(
    mean_temp = mean(temp, na.rm = T),
    mean_kelp = mean(pctcnpy, na.rm = T),
    mean_vis = mean(vis, na.rm = T)
  )


observer_experience <- length_data %>%
  group_by(year,month,observer) %>%
  summarise(n_obs = length(side)) %>%
  arrange(observer, year) %>%
  group_by(observer) %>%
  mutate(cumulative_n_obs = cumsum(n_obs),
         lifetime_obs = sum(n_obs)) %>%
  ungroup() %>%
  mutate(observer_ranking = percent_rank(lifetime_obs)) %>%
  mutate(trunc_observer = ifelse(observer_ranking > 0.5 | observer == 'unknown', observer,
                                 'infrequent')) %>%
  arrange(desc(observer_ranking)) %>%
  mutate(cumulative_n_obs = ifelse(trunc_observer == 'infrequent' | trunc_observer == 'unknown',0,cumulative_n_obs),
         n_obs = ifelse(trunc_observer == 'infrequent' | trunc_observer == 'unknown',0,n_obs),
         lifetime_obs = ifelse(trunc_observer == 'infrequent' | trunc_observer == 'unknown',0,lifetime_obs),
         observer_ranking = ifelse(trunc_observer == 'infrequent' | trunc_observer == 'unknown',0,observer_ranking)) %>%
  mutate(cumulative_n_obs_2 = cumulative_n_obs^2)


length_data <- length_data %>%
  left_join(observer_experience, by = c('year','month','observer'))

# ggmap::qmplot(x = mean_longitude, y = mean_latitude,data = species_distributions)

density_data <- read_csv('data/ci_reserve_data_final3 txt.csv') %>%
  magrittr::set_colnames(., tolower(colnames(.))) %>%
  gather('concat.name', 'value', grep('_', colnames(.)), convert = T) %>%
  mutate(
    data.type = gsub('\\_.*', '', concat.name),
    classcode = gsub('.*\\_', '', concat.name)
  ) %>%
  mutate(value = as.numeric(value)) %>%
  spread(data.type, value) %>%
  rename(site_side = site.side)



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
  summarise(mean_latitude = sum(count * latitude, na.rm = T) / sum(count, na.rm = T),
            mean_longitude = sum(count * longitude, na.rm = T) / sum(count, na.rm = T),
            max_latitude = max(latitude, na.rm = T),
            min_latitude = min(latitude, na.rm = T),
            max_longitude = max(longitude, na.rm = T),
            min_longitude = min(longitude, na.rm = T))

# ggmap::qmplot(min_longitude,min_latitude, data = species_distributions)


# ggmap::qmplot(longitude,latitude, color = side, data = site_coords)


if (file.exists('data/enso.csv')) {
  enso <- read_csv('data/enso.csv') %>%
    group_by(year) %>%
    summarise(mean_enso = mean(enso, na.rm = T)) %>%
    mutate(
      lag1_enso = dplyr::lag(mean_enso, 1),
      lag2_enso = dplyr::lag(mean_enso, 2),
      lag3_enso = dplyr::lag(mean_enso, 3),
      lag4_enso = dplyr::lag(mean_enso, 4)
    )

} else {
  scrape_enso(outdir = 'data/')

}

if (file.exists('data/pdo.csv')) {
  pdo <- read_csv('data/pdo.csv') %>%
    group_by(year) %>%
    summarise(mean_pdo = mean(pdo, na.rm = T)) %>%
    mutate(
      lag1_pdo = dplyr::lag(mean_pdo, 1),
      lag2_pdo = dplyr::lag(mean_pdo, 2),
      lag3_pdo = dplyr::lag(mean_pdo, 3),
      lag4_pdo = dplyr::lag(mean_pdo, 4)
    )

} else {
  scrape_pdo(outdir = 'data/')

  pdo <- read_csv('data/pdo.csv') %>%
    group_by(year) %>%
    summarise(mean_pdo = mean(pdo, na.rm = T)) %>%
    mutate(
      lag1_pdo = dplyr::lag(mean_pdo, 1),
      lag2_pdo = dplyr::lag(mean_pdo, 2),
      lag3_pdo = dplyr::lag(mean_pdo, 3),
      lag4_pdo = dplyr::lag(mean_pdo, 4)
    )

}


# deal with processed densities -------------------------------------------

reg_data <- density_data %>%
  select(biomass, site, side, site_side, year, classcode) %>%
  ungroup() %>%
  left_join(conditions_data, by = c('site', 'side', 'year')) %>%
  left_join(life_history_data %>% select(classcode, targeted, trophicgroup,
                                         commonname),
            by = 'classcode') %>%
  left_join(enso, by = 'year') %>%
  left_join(pdo, by = 'year') %>%
  mutate(
    any_seen = biomass > 0,
    log_density = log(biomass),
    targeted = as.numeric(targeted == 'Targeted'),
    post_mlpa = as.numeric(year >= 2003)
  )


# convert transect data to density estimates ------------------------------


if (file.exists('data/length-to-density-data.Rdata') == F |
    run_length_to_density == T) {

  density_example <- density_data %>%
    filter(is.na(biomass) == F & biomass > 0) %>%
    sample_n(1)

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


  length_to_density_data <- length_example %>%
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

  # species_sightings <- length_data %>%
  #   left_join(site_data, by = 'site') %>%
  #   group_by(region) %>%
  #   summarise(species_seen = list(unique(classcode)))
  #
  species_sightings <- length_to_density_data %>%
    left_join(site_data, by = 'site') %>%
    group_by(site) %>%
    summarise(species_seen = list(unique(classcode)))

  length_to_density_data <- length_to_density_data %>%
    ungroup() %>%
    left_join(site_data %>% select(site, region), by = 'site') %>%
    select(region, site, side, year,month,day,zone, level,transect) %>%
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
        observations = length_to_density_data,
        species_sightings = species_sightings,
        life_history_vars = colnames(life_history_data)
      )
    } %>%
    bind_rows()


  classcodes <- unique(length_to_density_data$classcode)

  life_history_vars <- which(colnames(length_to_density_data) %in% colnames(life_history_data))

  for (i in 1:length(classcodes)) {

    where_class <-
      length_to_density_data$classcode == classcodes[i] &
      length_to_density_data$count == 0

    temp_life_history <- life_history_data %>%
      filter(classcode == classcodes[i])

    temp_life_history <- temp_life_history[,colnames(length_to_density_data)[life_history_vars]]

    if (nrow(temp_life_history) > 1){stop('multiple classcodes')}

    length_to_density_data[where_class, life_history_vars] <- temp_life_history

    if (any(    (colnames(temp_life_history) == colnames(length_to_density_data[1,life_history_vars])
) == F)){stop()}

    if (any(is.na(length_to_density_data$classcode))){stop()}


  }

  check_life_history_fill <- length_to_density_data %>%
    group_by(classcode) %>%
    summarise(a = n_distinct(commonname)) %>%
    arrange(desc(a))

  if (any(check_life_history_fill$a > 1)){stop('multiple species per classcode')}

  save(file = paste0('data/length-to-density-data.Rdata'),
       length_to_density_data)

} else {
  load('data/length-to-density-data.Rdata')

}


# sum biomass across all observed sizes ---------------------------------


transect_covariates <- c('mean_depth','mean_vis','mean_temp','surge','mean_canopy','observer','n_obs','cumulative_n_obs')

mean_foo <- function(var){

  if (class(var) != 'character') {
    out <- mean(var, na.rm = T)
  } else{
    out <- paste(unique(var), collapse = '-')
  }

}

length_to_density_data <- length_to_density_data %>%
  group_by(campus, method, year, month, day, site, side, zone,transect,level, classcode) %>%
  mutate(total_biomass_g = sum(total_biomass_g),
         counter = 1:length(total_biomass_g)) %>%
  # purrrlyr::dmap_at(which(colnames(.) %in% transect_covariates), mean_foo) %>%
  ungroup() %>%
  filter(counter == 1) %>%
  select(-counter)

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

nearest_pisco_site <- site_coords[nearest_pisco_site$nn.idx[,1],]

kfm_pisco_locations <- kfm_locations %>%
  na.omit() %>%
  bind_cols(nearest_pisco_site)


kfm_data <- kfm_data %>%
  filter(sample_method !="visualfish" & sample_method!="crypticfish", data_source == 'kfm')%>%
  left_join(life_history_data, by = c('taxon_name' = 'taxa')) %>%
  left_join(kfm_pisco_locations %>% select(-data_source), by = 'site_id') %>%
  filter(str_detect(geolocation, '_island'),
         data_source != 'lter',
         is.na(geolocation) == F
         ) %>%
  mutate(region = str_extract(geolocation, '.*(?=_island)')) %>%
  filter(region %in% c('anacapa','santa_cruz','santa_rosa','san_miguel'))

kfm_data <- kfm_data %>%
  mutate(count = as.numeric(count),
         density = count / area,
         year = lubridate::year(date),
         month = lubridate::month(date))

region_table <- tribble(
  ~region,~short_region,
  'anacapa','ANA',
  'santa_cruz','SCI',
  'santa_rosa','SRI',
  'san_miguel','SMI'
)

kfm_data <- kfm_data %>%
  left_join(region_table, by = 'region') %>%
  select(-region) %>%
  rename(region = short_region)

kfm_data <- kfm_data %>%
  mutate(log_density = log(density),
         any_seen = density >0,
         factor_year = as.factor(year),
         factor_month = as.factor(month),
         total_biomass_g = density) %>%
  select(-region) #for compatibility with things later on




# deal with missing covariates --------------------------------------------
# pretty hacky for now, need to go back and deal with this better

length_to_density_data <- length_to_density_data %>%
  mutate(
    any_seen = total_biomass_g > 0,
    factor_year = factor(year),
    log_density = log(total_biomass_g),
    factor_month = factor(month),
    site_side = glue::glue('{site}-{side}')
  ) %>%
  group_by(site, side, month,year) %>%
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
  )


density_data <- reg_data %>%
  mutate(factor_year = factor(year),
         log_density = log(biomass),
         factor_month = factor(9)) %>%
  group_by(site, side, factor_month,year) %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_canopy = ifelse(is.na(mean_kelp), mean(mean_kelp, na.rm = T), mean_kelp)
  ) %>%
  group_by(site, side, year) %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_kelp = ifelse(is.na(mean_kelp), mean(mean_kelp, na.rm = T), mean_kelp)
  ) %>%
  group_by(site, year) %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_kelp = ifelse(is.na(mean_kelp), mean(mean_kelp, na.rm = T), mean_kelp)
  ) %>%
  group_by(year) %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_kelp = ifelse(is.na(mean_kelp), mean(mean_kelp, na.rm = T), mean_kelp)
  ) %>%
  ungroup() %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_kelp = ifelse(is.na(mean_kelp), mean(mean_kelp, na.rm = T), mean_kelp)
  )


# save raw-isa data -------------------------------------------------------

save(file = glue::glue("{run_dir}/rawish_ahnold_data.Rdata"), length_to_density_data,
     density_data, kfm_data)

# prepare data for model fitting ------------------------------------------

# length_to_density_data <- length_to_density_data %>%
#   left_join(observer_experience %>% select(observer, year,trunc_observer, cumulative_n_obs),
#             by = c('observer','year'))

consistent_sites <- length_to_density_data %>%
  group_by(site) %>%
  summarise(num_years = length(unique(year)),
            min_year = min(year),
            max_year = max(year)) %>%
  arrange(desc(num_years)) %>%
  filter(min_year <= 2000,
         num_years > 10)

consistent_regions <- length_to_density_data %>%
  left_join(site_data, by = 'site') %>%
  group_by(region) %>%
  summarise(num_years = length(unique(year)),
            min_year = min(year),
            max_year = max(year)) %>%
  arrange(desc(num_years)) %>%
  filter(min_year <= 2000,
         num_years > 10)

raw_length_covars <-
  paste(c(
    'region',
    'zone',
    'site_side',
    'level',
    'mean_vis',
    'surge',
    'factor_month',
    'trunc_observer',
    'cumulative_n_obs',
    'cumulative_n_obs_2'
  ),
  collapse = '+')

prob_raw_length_covars <-
  paste(c(
    'region',
    'zone',
    'site_side',
    'level',
    'mean_vis',
    'surge',
    'factor_month',
    'trunc_observer',
    'cumulative_n_obs',
    'cumulative_n_obs_2'
  ),
  collapse = '+')

kfm_length_covars <-
  paste(c(
    'region',
    'factor_month'
  ),
  collapse = '+')

prob_kfm_length_covars <-
  paste(c(
    'region',
    'factor_month'
  ),
  collapse = '+')

# prob_raw_length_covars <-
#   paste(c('region','mean_vis', 'mean_canopy','factor_month'),
#         collapse = '+')

supplied_density_covars <-
  paste(c('region','mean_kelp', 'mean_vis','site_side'),
        collapse = '+')

prob_supplied_density_covars <-
  paste(c('region','mean_kelp', 'mean_vis','site_side'),
        collapse = '+')


length_to_density_models <- length_to_density_data %>%
  nest(-classcode) %>%
  mutate(ind_covars = raw_length_covars,
         prob_ind_covars = prob_raw_length_covars) %>%
  mutate(data_source = 'length_to_density')


supplied_density_models <- density_data %>%
  nest(-classcode) %>%
  mutate(ind_covars = supplied_density_covars,
         prob_ind_covars = prob_supplied_density_covars) %>%
  mutate(data_source = 'supplied_density')


kfm_to_density_models <- kfm_data %>%
  nest(-classcode) %>%
  mutate(ind_covars = kfm_length_covars,
         prob_ind_covars = prob_kfm_length_covars) %>%
  mutate(data_source = 'kfm_density')

# filter data -------------------------------------------------------------

nobs_quantiles <- length_to_density_data %>%
  filter(any_seen == T) %>%
  group_by(year, classcode) %>%
  summarise(nobs = sum(any_seen)) %>%
  ungroup() %>%
  {quantile(.$nobs)}

well_observed_species <- length_to_density_data %>%
  filter(year > 1999) %>%
  group_by(year, classcode,commonname, targeted) %>%
  summarise(nseen = sum(any_seen, na.rm = T)) %>%
  group_by(commonname,classcode, targeted) %>%
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


abundance_models <- length_to_density_models %>%
  bind_rows(supplied_density_models) %>%
  bind_rows(kfm_to_density_models) %>%
  filter(classcode %in% well_observed_species$classcode) %>%
  mutate(data = map(data, ~ left_join(.x, site_data, by = c('site', 'side')))) %>%
  mutate(data = map2(data,classcode, filterfoo, min_year = min_year,
                     min_seen_years = 14,
                    filter_level = quo(classcode))) %>% # filter out things
  mutate(dim_data = map_dbl(data, nrow)) %>%
  filter(dim_data > 0)

# prepare data for abundance estimates aggregate, augment, center-scale data ----------------------------------------------

abundance_models <- abundance_models %>%
  mutate(classcode = tolower(classcode)) %>%
  left_join(life_history_data %>% select(classcode, commonname, targeted),
            by = 'classcode') %>%
  # left_join(fished_species, by = 'classcode') %>%
  # mutate(targeted = ifelse(
  #   fished == 1 &
  #     is.na(fished) == F &
  #     targeted == "Non-targeted",
  #   'Targeted',
  #   targeted
  # )) %>%
  filter(str_detect(commonname, 'YOY') == F,
         is.na(targeted) == F,
         str_detect(classcode,'_yoy') == F) #%>%
  # mutate(data = map(data, ~ purrrlyr::dmap_if(.x, is.numeric, center_scale))) # center and scale continuos data


population_structure <- c('one-pop', 'regional-pops', 'mpa-pops')

population_filtering <- c('all', 'mpa-only', 'consistent-sites')

abundance_models <- cross_df(list(data = list(abundance_models),
         population_structure = population_structure,
         population_filtering = population_filtering)) %>%
  unnest()

# run vast ----------------------------------------------------------------

if (run_vast == T) {


  vast_abundance <- abundance_models %>%
    # select(classcode, data_source, data, ) %>%
    mutate(survey_region = 'california_current',
           n_x = num_knots) %>%
    mutate(
      vast_data = map2(
        data,
        classcode,
        vast_prep,
        site_coords = site_coords,
        conditions_data = conditions_data
      )
    )


  arg <- safely(vasterize_pisco_data)

  vast_abundance <- vast_abundance %>%
    mutate(vast_results = purrr::pmap(
      list(
        region = survey_region,
        raw_data = vast_data,
        n_x = n_x
      ),
      arg,
      run_dir = run_dir,
      nstart = 100,
      obs_model = c(2, 0)
    ))

  vast_abundance <- vast_abundance %>%
    mutate(vast_error = map(vast_results, 'error'))

  save(file = paste0(run_dir, '/vast_abundance.Rdata'), vast_abundance)
} else {
  load(file = paste0(run_dir, '/vast_abundance.Rdata'))
}


vast_abundance <- vast_abundance %>%
  mutate(vast_index = map(vast_results, 'result'))

species_comp_by_dbase_plot <- abundance_models %>%
  ggplot(aes(commonname, dim_data, fill = data_source)) +
  geom_col(position = 'dodge') +
  coord_flip()

species_comp_and_targeted_by_dbase_plot <- abundance_models %>%
  ggplot(aes(commonname, dim_data, fill = targeted)) +
  geom_col(position = 'dodge') +
  coord_flip() +
  facet_wrap(~ data_source)

# estimate abundance through delta-glm ------------------------------------

safely_fit_fish <- safely(fit_fish)

abundance_models <- abundance_models %>%
  # filter(data_source == 'kfm_density') %>%
    mutate(
    seen_model = pmap(
      list(data = data,
           ind_covars = ind_covars,
           pop_structure = population_structure,
           pop_filter = population_filtering),
      safely_fit_fish,
      dep_var = 'log_density',
      fit_seen = T,
      family = 'gaussian',
      model_type = 'glm',
      consistent_sites = consistent_sites
    ),
    seeing_model = pmap(
      list(data = data,
           ind_covars = prob_ind_covars,
           pop_structure = population_structure,
           pop_filter = population_filtering),
      safely_fit_fish,
      dep_var = 'any_seen',
      fit_seen = F,
      family = 'binomial',
      model_type = 'glm',
      consistent_sites = consistent_sites
    )
  )

abundance_models <- abundance_models %>%
  mutate(
    seeing_error = map(seeing_model, 'error'),
    seen_error =  map(seen_model, 'error'),
    seeing_model = map(seeing_model, 'result'),
    seen_model =  map(seen_model, 'result')
  )


#filter out models that didn't converge for some reason
abundance_models <- abundance_models %>%
  mutate(no_error = map2_lgl(seeing_error, seen_error, ~ is.null(.x) &
                               is.null(.y))) %>%
  filter(no_error == T)

abundance_models <- abundance_models %>%
  mutate(
    seen_coefs = map(seen_model, broom::tidy),
    seen_aug = map(seen_model, broom::augment),
    seeing_coefs = map(seeing_model, broom::tidy),
    seeing_aug = map(seeing_model, broom::augment)
  )

add_covars_foo <- function(data){

  if(is.null(data$region)){data$region <- 'somewhere'}

  if(is.null(data$eventual_mpa)){data$eventual_mpa <- TRUE}

return(data)
}

abundance_models <- abundance_models %>%
  mutate(seen_aug = map(seen_aug, add_covars_foo),
         seeing_aug = map(seeing_aug, add_covars_foo))

abundance_models <- abundance_models %>%
  mutate(seen_rank_deficient = map_lgl(seen_model, ~ length(.$coefficients) > .$rank),
         seeing_rank_deficient = map_lgl(seeing_model, ~ length(.$coefficients) > .$rank))


abundance_models <- abundance_models %>%
  mutate(abundance_index = pmap(
    list(
      seen_model = seen_model,
      seeing_model = seeing_model,
      seeing_aug = seeing_aug,
      seen_aug = seen_aug,
      pop_structure = population_structure
    ),
    create_abundance_index))

abundance_plot_foo <- function(data,pop_structure){

  if (pop_structure == 'one-pop'){
  outplot <- ggplot(
    data,
    aes(
      factor_year %>% as.character() %>% as.numeric(),
      abundance_index
    )
  ) + geom_line() +
    geom_point()
  }
  if (pop_structure == 'regional-pops'){

    outplot <- ggplot(
      data,
      aes(
        factor_year %>% as.character() %>% as.numeric(),
        abundance_index,
        color = region
      )
    ) +
      geom_line() +
      geom_point()
  }
  if (pop_structure == 'mpa-pops'){

    outplot <- ggplot(
      data,
      aes(
        factor_year %>% as.character() %>% as.numeric(),
        abundance_index,
        color = eventual_mpa
      )
    ) +
      geom_line() +
      geom_point()
  }
  return(outplot)

}


abundance_models <- abundance_models %>%
  mutate(abundance_plot = map2(abundance_index, population_structure, abundance_plot_foo))

abundance_models <- abundance_models %>%
  mutate(num_years = map_dbl(abundance_index, ~ nrow(.x)))

abundance_models <- abundance_models %>%
  mutate(classcode = tolower(classcode))

save_foo <-
  function(species,
           pop_structure,
           pop_filtering,
           data_source,
           abundance_plot,
           run_dir) {
    ggsave(
      file = paste0(
        run_dir,
        '/',
        species,
        '-',
        pop_structure,
        '-',
        pop_filtering,
        '-',
        data_source,
        '.pdf'
      ),
      abundance_plot,
      width = 8,
      height = 8
    )
  }

pwalk(
  list(species = abundance_models$commonname,
       pop_structure = abundance_models$population_structure,
       pop_filtering = abundance_models$population_filtering,
       data_source = abundance_models$data_source,
       abundance_plot = abundance_models$abundance_plot),
       save_foo,
  run_dir = run_dir
)


# prepare abundance indicies assess standardization --------------------------------------------------



calc_raw_abundance <- function(data, population_structure) {
  if (any(colnames(data) == 'biomass')) {
    data$total_biomass_g <- data$biomass

  }

  if (population_structure == 'one-pop'){

  raw_trend <- data %>%
    group_by(factor_year) %>%
    summarise(abundance_index = mean(total_biomass_g, na.rm = T)) %>%
    ungroup() %>%
    # mutate(abundance_index = center_scale(abundance_index + 1e-3)) %>%
    mutate(year = factor_year %>% as.character() %>% as.numeric()) %>%
    ungroup()
  }
  if (population_structure == 'regional-pops'){

    raw_trend <- data %>%
      group_by(region,factor_year) %>%
      summarise(abundance_index = mean(total_biomass_g, na.rm = T)) %>%
      ungroup() %>%
      # mutate(abundance_index = center_scale(abundance_index + 1e-3)) %>%
      mutate(year = factor_year %>% as.character() %>% as.numeric()) %>%
      ungroup()
  }
  if (population_structure == 'mpa-pops'){

    raw_trend <- data %>%
      group_by(eventual_mpa,factor_year) %>%
      summarise(abundance_index = mean(total_biomass_g, na.rm = T)) %>%
      ungroup() %>%
      # mutate(abundance_index = center_scale(abundance_index + 1e-3)) %>%
      mutate(year = factor_year %>% as.character() %>% as.numeric()) %>%
      ungroup()
  }

  return(raw_trend)

}

process_vast_foo <- function(data){


  out <- data$time_index %>%
    select(Year, abundance) %>%
    rename(year = Year, abundance_index = abundance) %>%
    mutate(population_level = 'one-pop')
}

vast_abundance <- vast_abundance %>%
  mutate(vast_abundance_index = map(vast_index, process_vast_foo))

simplify_trend <- function(dat, population_structure){

  nobs <- nrow(dat)

  out_frame <- data_frame(year = rep(NA,nobs), abundance_index = rep(NA,nobs), population_level = rep(NA,nobs))

  out_frame$year <- dat$factor_year %>% as.character() %>% as.numeric()

  out_frame$abundance_index <- dat$abundance_index

  if (population_structure == 'one-pop'){

    out_frame$population_level <- 'channel-islands'

  }
  if (population_structure == 'regional-pops'){
    out_frame$population_level <- dat$region
  }
  if (population_structure == 'mpa-pops'){

    out_frame$population_level <- dat$eventual_mpa %>% as.character()

  }


  return(out_frame)

}


abundance_indices <- abundance_models %>%
  left_join(vast_abundance %>%   select(classcode, population_filtering, population_structure, data_source, vast_abundance_index)
, by = c('classcode','population_filtering', 'population_structure','data_source')) %>%
  mutate(raw_abundance_trend = map2(data,population_structure, calc_raw_abundance)) %>%
  mutate(raw_abundance_index = map2(raw_abundance_trend, population_structure, simplify_trend),
         glm_abundance_index = map2(abundance_index, population_structure, simplify_trend))


compare_trends <- abundance_indices %>%
  nest(-classcode,-commonname)

walk2(compare_trends$commonname,
      compare_trends$data,
      safely(compare_trend_foo),
      run_dir = run_dir)

rm(compare_trends)

# a <- object.size(abundance_indices)

abundance_indices <- abundance_indices %>%
  mutate(seen_model = map(seen_model, ~strip::strip(.x, keep = c('predict','print'))),
         seeing_model = map(seeing_model, ~strip::strip(.x, keep = c('predict','print'))))

walk(abundance_indices, ~print(object.size(.x), units = 'Gb'))

# b <- object.size(ai2)

save(file = paste0(run_dir, '/abundance_indices.Rdata'),
     abundance_indices)

# prepare candidate did data and runs ----------------------------------------------


annual_regional_conditions <- conditions_data %>%
  left_join(site_data %>% select(site, region), by = 'site') %>%
  group_by(region, year) %>%
  summarise(
    mean_annual_temp = mean(mean_temp, na.rm = T),
    mean_annual_kelp = mean(mean_kelp, na.rm  = T)
  ) %>%
  gather(variable, value,-year,-region) %>%
  group_by(variable) %>%
  mutate(value = zoo::na.approx(value),
         value = center_scale(value)) %>%
  spread(variable, value)

annual_conditions <- conditions_data %>%
  left_join(site_data %>% select(site, region), by = 'site') %>%
  group_by(year) %>%
  summarise(
    mean_annual_temp = mean(mean_temp, na.rm = T),
    mean_annual_kelp = mean(mean_kelp, na.rm  = T)
  ) %>%
  gather(variable, value,-year) %>%
  group_by(variable) %>%
  mutate(value = zoo::na.approx(value),
         value = center_scale(value)) %>%
  spread(variable, value)


did_data <- abundance_indices %>%
  select(classcode, population_filtering, population_structure, data_source, contains('_index')) %>%
  select(-abundance_index) %>%
  gather('abundance_source','abundance_index',contains('_index')) %>%
  filter(!(abundance_source == 'vast_abundance_index' & data_source == 'kfm_density')) %>%
  unnest() %>%
  left_join(life_history_data, by = 'classcode') %>%
  left_join(enso, by = 'year') %>%
  left_join(pdo, by = 'year') %>%
  left_join(ci_catches, by = c('classcode', 'year')) %>%
  left_join(species_distributions, by = 'classcode') %>%
  mutate(catch = ifelse(is.na(catch), 0, catch)) %>%
  mutate(targeted = as.numeric(targeted == 'Targeted'),
         post_mpa = as.numeric(year >= 2003)) %>%
  mutate(abundance_index = abundance_index + 1e-6)


annual_conditions_foo <- function(population_structure, data, abundance_source, annual_conditions,annual_regional_conditions){

  if (population_structure == 'regional-pops' &  abundance_source != 'vast_abundance_index'){
    data <- data %>%
      left_join(annual_regional_conditions, by = c('year','population_level' = 'region'))
  } else{

    data <- data %>%
      left_join(annual_conditions, by = 'year')
  }

return(data)
}

did_data <- did_data %>%
  nest(-population_structure, -abundance_source) %>%
  mutate(data = pmap(list(population_structure = population_structure,
                          data = data,
                          abundance_source = abundance_source), annual_conditions_foo, annual_conditions,
                     annual_regional_conditions)) %>%
           unnest() %>%
  filter(year <=2013)


compare_annual_abundance <- function(data){


  out <-  data %>%
    group_by(abundance_source) %>%
    mutate(abundance_index = (abundance_index - mean(abundance_index)) / (2 * sd(abundance_index))) %>%
    ggplot(aes(year, abundance_index, color = abundance_source)) +
    geom_line() +
    geom_point()


}



annual_abundance_trends_foo <-
  function(classcode,
           data_source,
           annual_abundance_trends,
           run_dir) {
    ggsave(file = paste0(run_dir, '/', classcode, '-', data_source, '-annual-abundance-trends.pdf'),
           annual_abundance_trends, height = 8, width = 8)
  }

trend_data <- did_data %>%
  filter(population_structure == 'one-pop',
         population_filtering ==  'all') %>%
  nest(-classcode,-data_source) %>%
  mutate(annual_abundance_trends = map(data, compare_annual_abundance))

pwalk(list(classcode = trend_data$classcode, data_source = trend_data$data_source,
           annual_abundance_trends = trend_data$annual_abundance_trends), annual_abundance_trends_foo,
      run_dir = run_dir)


did_terms <- did_data %>%
  select(year, targeted) %>%
  mutate(index = 1:nrow(.)) %>%
  spread(year, targeted, fill = 0) %>%
  select(-index) %>%
  set_names(., paste0('did_', colnames(.))) %>%
  select(-did_2000)

year_mpa <- 2003

generation_did_terms <- did_data %>%
  select(year, targeted, tm) %>%
  mutate(index = 1:nrow(.),
         years_protected = pmax(0,year - year_mpa),
         generations_protected = round(years_protected / round(tm))) %>%
  select(index,generations_protected, targeted) %>%
  spread(generations_protected, targeted, fill = 0) %>%
  select(-index) %>%
  set_names(., paste0('did_', colnames(.))) %>%
  select(-did_0)

generation_did_terms <- generation_did_terms[, colSums(generation_did_terms) >0]


recruitment_did_terms <- did_data %>%
  select(year, targeted, first_age_seen) %>%
  mutate(index = 1:nrow(.),
         years_protected = pmax(0,year - year_mpa),
         recruits_protected = round(years_protected / round(first_age_seen))) %>%
  select(index,recruits_protected, targeted) %>%
  spread(recruits_protected, targeted, fill = 0) %>%
  select(-index) %>%
  set_names(., paste0('did_', colnames(.))) %>%
  select(-did_0)

# did_terms <- generation_did_terms


did_models <-
  data_frame(
    did_data = list(did_data),
    did_terms = list(did_terms, generation_did_terms, recruitment_did_terms),
    did_term_names = c('years-protected','generations-protected','recruits-protected')
  )


fit_did <- function(did_data, did_terms){

did_data <- did_data %>%
  bind_cols(did_terms)

did_data <- did_data %>%
  nest(-population_structure, -data_source, -population_filtering,-abundance_source)

# visually inspect abundance trends ---------------------------------------


did_data <- did_data %>%
  mutate(correlation_tests = map(data, test_parallel_trends)) %>%
  mutate(pre_correlation = map_dbl(correlation_tests, ~.x$overall_correlation_test$estimate),
         pre_correlation_signif = map_dbl(correlation_tests, ~.x$overall_correlation_test$p.value))


# fit DiD estimator on abundance indicies ---------------------------------


did_reg <-
  paste0('log(abundance_index) ~', paste(
    c(
      'targeted',
      'mean_annual_kelp',
      '(1|factor_year)',
      'mean_pdo',
      # 'lag1_pdo',
      # 'lag2_pdo',
      '((1 + mean_annual_temp)|classcode)',
      colnames(did_terms)),
    collapse = '+'
  ))

did_models <- did_data %>%
  mutate(did_reg = did_reg) %>%
  mutate(did_model = map2(data, did_reg, ~ lme4::glmer(.y, data = .x %>% mutate(factor_year = factor(year)))))
} # close fit_did


did_models <- did_models %>%
  mutate(did_model = map2(did_data, did_terms, fit_did)) %>%
  select(did_model,did_term_names) %>%
  unnest()

# check_model <- did_models %>%
#   filter(data_source == 'length_to_density',
#          population_structure == 'one-pop',
#          abundance_source == 'glm_abundance_index',
#          population_filtering == 'all',
#          did_term_names == 'years-protected')
#
# check_model$did_model[[1]] %>% broom::glance()

did_plot_foo <- function(x) {
  x %>%
    broom::tidy() %>%
    filter(str_detect(term, 'did')) %>%
    mutate(year = str_replace(term, 'did_', '') %>% as.numeric()) %>%
    ggplot() +
    geom_pointrange(aes(
      year,
      estimate,
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    )) +
    # geom_vline(aes(xintercept = 2003),
    #            color = 'red',
    #            linetype = 2) +
    geom_hline(aes(yintercept = 0))
}

did_models <- did_models %>%
  mutate(did_plot = map(did_model, did_plot_foo) ,
         did_diagnostics = map(did_model, diagnostic_plots))

did_plot_foo <-
  function(data_source,
           abundance_source,
           population_filtering,
           population_structure,
           did_term_name,
           did_plot,
           run_dir) {
    ggsave(file = paste0(run_dir, '/', data_source, '-', abundance_source,'-',population_filtering,'-',population_structure,
                         '-',did_term_name, '-did.pdf'),
           did_plot, height = 8, width = 8)
  }

pwalk(
  list(
    data_source = did_models$data_source,
    abundance_source =  did_models$abundance_source,
    population_filtering = did_models$population_filtering,
    population_structure = did_models$population_structure,
    did_term_name = did_models$did_term_names,
    did_plot = did_models$did_plot
  ),
  did_plot_foo,
  run_dir = run_dir
)

check_ahnold(length_to_density_data = length_to_density_data,
             abundance_indices = abundance_indices,
             did_models = did_models)

save(file = glue::glue('{run_dir}/did_models.Rdata'), did_models)

