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

aggregate_transects  <-
  F # should transects be aggregated up (mean across observer) or left as raw

channel_islands_only <- T # only include channel islands, leave T

min_year <- 1999 # included years must be greater than this

occurance_ranking_cutoff <- 0.5

small_num <-  0 # no idea

use_mpa_site_effects <- F # no idea

# load data ---------------------------------------------------------------

length_data <- read_csv('data/UCSB_FISH raw thru 2013.csv') %>%
  magrittr::set_colnames(., tolower(colnames(.)))

life_history_data <-
  read_csv('data/VRG Fish Life History in MPA_04_08_11_12 11-Mar-2014.csv') %>%
  rename(classcode = pisco_classcode) %>%
  mutate(classcode = tolower(classcode)) %>%
  # rename(description_2 = Description) %>%
  magrittr::set_colnames(., tolower(colnames(.)))

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

length_data <- length_data %>%
  left_join(life_history_data %>% mutate(classcode = toupper(classcode)), by = 'classcode') #%>%
# left_join(site_data, by = c('site', 'side'))

ci_catches <-
  read_csv(file = file.path('data', 'cfdw-channel-islands-catches.csv')) %>% group_by(classcode, year) %>%
  summarise(catch = sum(pounds_caught, na.rm = T))

fished_species <-
  data_frame(classcode = unique(ci_catches$classcode),
             fished = 1)


conditions_data <- length_data %>%
  group_by(site, side, year) %>%
  summarise(
    mean_temp = mean(temp, na.rm = T),
    mean_kelp = mean(pctcnpy, na.rm = T),
    mean_vis = mean(vis, na.rm = T)
  )


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
# # select_(.dots = as.list(reg_vars)) %>%
# map2_df(
#   colnames(.),
#   center_scale,
#   omit_names = c('log_density','biomass', 'year', 'mean_enso', 'mean_pdo',
#                  'targeted','year_mpa',paste0('lag',1:4,'_enso'),paste0('lag',1:4,'_pdo'))
# ) %>%


# convert transect data to density estimates ------------------------------

if (file.exists('data/length-to-density-data.Rdata') == F |
    run_length_to_density == T) {
  density_example <- density_data %>%
    filter(is.na(biomass) == F & biomass > 0) %>%
    sample_n(1)

  length_example <- length_data %>%
    filter(
      classcode == toupper(density_example$classcode)
      ,
      site == density_example$site,
      side == density_example$side,
      year == density_example$year
    )

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
      surge = ifelse(is.na(observer), 'unknown', surge)
    ) %>%
    group_by(classcode, site, side, year, transect, observer) %>%
    summarise(
      total_biomass_g = sum(biomass_g),
      mean_temp = mean(temp, na.rm = T),
      mean_vis = mean(vis, na.rm = T),
      mean_depth = mean(depth, na.rm = T),
      mean_canopy = mean(pctcnpy, na.rm = T)
    )

  length_to_density_data %>%
    summarise(nobs = length(classcode))

  species_sightings <- length_data %>%
    left_join(site_data, by = 'site') %>%
    group_by(region) %>%
    summarise(species_seen = list(unique(classcode)))



  length_to_density_data <- length_to_density_data %>%
    ungroup() %>%
    left_join(site_data %>% select(site, region), by = 'site') %>%
    select(region, site, side, year, transect) %>%
    unique() %>%  {
      pmap(
        list(
          this_region = .$region,
          this_site = .$site,
          this_side = .$side,
          this_year = .$year,
          this_transect = .$transect
        ),
        add_missing_fish,
        observations = length_to_density_data,
        species_sightings = species_sightings
      )
    } %>%
    bind_rows()

  save(file = paste0('data/length_to_density_data.Rdata'),
       length_to_density_data)

} else {
  load('data/length-to-density-data.Rdata')

}


if (aggregate_transects == T) {
  length_to_density_data <- length_to_density_data %>%
    group_by(classcode, observer, site, side, year) %>%
    summarise(
      mean_biomass_g = mean(total_biomass_g, na.rm = T),
      mean_temp = mean(mean_temp, na.rm = T),
      mean_vis = mean(mean_vis, na.rm = T),
      mean_depth = mean(mean_depth, na.rm = T),
      mean_canopy = mean(mean_canopy, na.rm = T)
    ) %>%
    mutate(
      biomass_g_per_m2 = mean_biomass_g / (30 * 4),
      biomass_g_per_hectare = biomass_g_per_m2 * 10000,
      biomass_ton_per_hectare = biomass_g_per_hectare * 1e-6
    ) %>%
    ungroup()
} else {
  length_to_density_data <- length_to_density_data %>%
    rename(mean_biomass_g = total_biomass_g) %>%
    mutate(
      biomass_g_per_m2 = mean_biomass_g / (30 * 4),
      biomass_g_per_hectare = biomass_g_per_m2 * 10000,
      biomass_ton_per_hectare = biomass_g_per_hectare * 1e-6
    ) %>%
    ungroup()

}


# deal with missing covariates --------------------------------------------
# pretty hacky for now, need to go back and deal with this better


length_to_density_data <- length_to_density_data %>%
  mutate(
    any_seen = mean_biomass_g > 0,
    factor_year = factor(year),
    log_density = log(mean_biomass_g)
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
         log_density = log(biomass)) %>%
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


# prepare data for model fitting ------------------------------------------

raw_length_vars <-
  paste(c(
    'factor_year:region',
    'observer',
    'site',
    'mean_vis',
    'mean_canopy'
  ),
  collapse = '+')

supplied_density_vars <-
  paste(c('factor_year:region', 'site', 'mean_kelp', 'mean_vis'),
        collapse = '+')

prob_raw_length_vars <-
  paste(c('factor_year', 'observer', 'site', 'mean_vis', 'mean_canopy'),
        collapse = '+')

prob_supplied_density_vars <-
  paste(c('factor_year', 'region', 'site', 'mean_kelp', 'mean_vis'),
        collapse = '+')

length_to_density_models <- length_to_density_data %>%
  nest(-classcode) %>%
  mutate(ind_vars = raw_length_vars,
         prob_ind_vars = prob_raw_length_vars)

supplied_density_models <- density_data %>%
  nest(-classcode) %>%
  mutate(ind_vars = supplied_density_vars,
         prob_ind_vars = prob_supplied_density_vars)

length_to_density_models <- length_to_density_models %>%
  mutate(data_source = 'length_to_density')

supplied_density_models <- supplied_density_models %>%
  mutate(data_source = 'supplied_density')

# filter data -------------------------------------------------------------

filterfoo <-
  function(x,
           min_seen_years = 8,
           mpa_start_year = 2003) {
    # only years above 1999 at the main channel islands

    x <- x %>%
      filter(is.na(log_density) == F) %>%
      filter(year > 1999 & region %in% c('ANA', 'SCI', 'SMI', 'SRI')) %>%
      group_by(region) %>%
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
      ) # actually observed for at least 10 years
  }


abundance_models <- length_to_density_models %>%
  bind_rows(supplied_density_models) %>%
  mutate(data = map(data, ~ left_join(.x, site_data, by = c('site', 'side')))) %>%
  mutate(data = map(data, (filterfoo))) %>% # filter out things
  mutate(dim_data = map_dbl(data, nrow)) %>%
  filter(dim_data > 0)


# aggregate, augment, center-scale data ----------------------------------------------

abundance_models <- abundance_models %>%
  mutate(classcode = tolower(classcode)) %>%
  left_join(life_history_data %>% select(classcode, commonname, targeted),
            by = 'classcode') %>%
  left_join(fished_species, by = 'classcode') %>%
  mutate(targeted = ifelse(
    fished == 1 &
      is.na(fished) == F &
      targeted == "Non-targeted",
    'Targeted',
    targeted
  )) %>%
  filter(str_detect(commonname, 'YOY') == F, is.na(targeted) == F) %>%
  mutate(data = map(data, ~ purrrlyr::dmap_if(.x, is.numeric, center_scale))) # center and scale continuos data


vast_abundance <- abundance_models %>%
  select(classcode, data_source, data) %>%
  mutate(survey_region = 'california_current',
         n_x = 200) %>%
  mutate(vast_data = map2(data,classcode, vast_prep, site_coords = site_coords, conditions_data = conditions_data))

if (run_vast == T){
  arg <- safely(vasterize_pisco_data)

  vast_abundance <- vast_abundance %>%
    mutate(
      vast_results = purrr::pmap(
        list(
          region = survey_region,
          raw_data = vast_data,
          n_x = n_x
        ),
        arg,
        run_dir = run_dir,
        nstart = 100,
        obs_model = c(2, 0)
      )
    )

  vast_abundance <- vast_abundance %>%
    mutate(vast_error = map(vast_results, 'error'))

  save(file = paste0(run_dir,'/vast_abundance.Rdata'), vast_abundance)
} else {

  load(file = paste0(run_dir,'/vast_abundance.Rdata'))


}


vast_abundance <- vast_abundance %>%
  mutate(vast_index = map(vast_results,'result'))


species_comp_by_dbase_plot <- abundance_models %>%
  ggplot(aes(commonname, dim_data, fill = data_source)) +
  geom_col(position = 'dodge') +
  coord_flip()

species_comp_and_targeted_by_dbase_plot <- abundance_models %>%
  ggplot(aes(commonname, dim_data, fill = targeted)) +
  geom_col(position = 'dodge') +
  coord_flip() +
  facet_wrap( ~ data_source)


# estimate abundance through delta-glm ------------------------------------

safely_fit_fish <- safely(fit_fish)

abundance_models <- abundance_models %>%
  mutate(seen_model = map2(data,ind_vars, safely_fit_fish, dep_var = 'log_density', fit_seen = T, family = 'gaussian'),
         seeing_model =map2(data,prob_ind_vars, safely_fit_fish, dep_var = 'any_seen', fit_seen = F, family = 'binomial') )

abundance_models <- abundance_models %>%
  mutate(seeing_error = map(seeing_model, 'error'),
         seen_error =  map(seen_model, 'error'),
         seeing_model = map(seeing_model, 'result'),
         seen_model =  map(seen_model, 'result'))

abundance_models <- abundance_models %>%
  mutate(no_error = map2_lgl(seeing_error, seen_error, ~is.null(.x) & is.null(.y))) %>%
  filter(no_error == T) #filter out models that didn't converge for some reason


abundance_models <- abundance_models %>%
  mutate(seen_coefs = map(seen_model, broom::tidy),
         seen_aug = map(seen_model, broom::augment),
         seeing_coefs = map(seeing_model, broom::tidy),
         seeing_aug = map(seeing_model, broom::augment))

abundance_models <- abundance_models %>%
  mutate(abundance_index = pmap(
    list(
      seen_model = seen_model,
      seeing_model = seeing_model,
      seeing_aug = seeing_aug,
      seen_aug = seen_aug
    ),
    create_abundance_index,
    model_resolution = 'regional'
  ))

abundance_models <- abundance_models %>%
  mutate(abundance_plot = map(abundance_index, ~ ggplot(.x,aes(factor_year %>% as.character() %>% as.numeric(), abundance_index, color = region)) + geom_line()))

abundance_models <- abundance_models %>%
  mutate(num_years = map_dbl(abundance_index, ~ nrow(.x)))

abundance_models <- abundance_models %>%
  mutate(classcode = tolower(classcode))

walk2(abundance_models$commonname, abundance_models$abundance_plot,
      ~ ggsave(file = paste0(run_dir,'/',.x,'.pdf'), .y), run_dir = run_dir)


# prepare abundance indicies assess standardization --------------------------------------------------

calc_raw_abundance <- function(data) {

  if (any(colnames(data) == 'biomass')){

    data$mean_biomass_g <- data$biomass

  }

  raw_trend <- data %>%
    group_by(region,factor_year) %>%
    summarise(abundance_index = mean(mean_biomass_g, na.rm = T)) %>%
    ungroup() %>%
    mutate(abundance_index = center_scale(abundance_index + 1e-3)) %>%
    mutate(year = factor_year %>% as.character() %>% as.numeric()) %>%
    ungroup()

}

abundance_indices <- abundance_models %>%
  mutate(standardized_abundance_trend = map(abundance_index, ~select(.,factor_year,region,abundance_index) %>% mutate(year = as.character(factor_year) %>% as.numeric()))) %>%
  mutate(raw_abundance_trend = map(data,calc_raw_abundance )) %>%
  select(classcode,data_source,standardized_abundance_trend,raw_abundance_trend)


raw_model_indices <- abundance_indices %>%
  select(classcode, data_source, raw_abundance_trend) %>%
  unnest() %>%
  mutate(abundance_source = 'raw')

standardized_model_indices <- abundance_indices %>%
  select(classcode, data_source, standardized_abundance_trend) %>%
  unnest() %>%
  mutate(abundance_source = 'standardized')

vast_model_indices <- vast_abundance %>%
  mutate(time_index = map(vast_index, 'time_index')) %>%
  select(classcode, data_source,time_index) %>%
  unnest() %>%
  mutate(factor_year = as.factor(Year),
         region = 'SCI') %>%
  rename(year = Year) %>%
  select(classcode, data_source,factor_year, region, abundance,year) %>%
  rename(abundance_index = abundance) %>%
  set_names(tolower) %>%
  mutate(abundance_source = 'vast') %>%
  group_by(classcode, data_source) %>%
  mutate(abundance_index = center_scale(abundance_index + 1e-3)) %>%
  ungroup()

abundance_indices <- raw_model_indices %>%
  bind_rows(standardized_model_indices) %>%
  bind_rows(vast_model_indices) %>%
  nest(-classcode)

comp_plot_foo <- function(data){
  data %>%
    ggplot(aes(year, abundance_index, color = abundance_source)) +
    geom_line() +
    geom_point() +
    facet_grid(data_source ~ region)


}

abundance_indices <- abundance_indices %>%
  mutate(comp_plot = map(data, comp_plot_foo)) %>%
  mutate(classcode = tolower(classcode)) %>%
  left_join(life_history_data %>% select(classcode, commonname))


save(file = paste0(run_dir,'/abundance_indices.Rdata'), abundance_indices)


walk2(abundance_indices$commonname, abundance_indices$comp_plot,
      ~ ggsave(file = paste0(run_dir,'/',.x,'-modelcomp.pdf'), .y), run_dir = run_dir)


# fit DiD estimator on abundance indicies ---------------------------------

annual_conditions <- conditions_data %>%
  left_join(site_data %>% select(site, region), by = 'site') %>%
  group_by(region,year) %>%
  summarise(mean_annual_temp = mean(mean_temp, na.rm = T),
            mean_annual_kelp = mean(mean_kelp, na.rm  = T)) %>%
  gather(variable, value,-year,-region) %>%
  group_by(variable) %>%
  mutate(value = zoo::na.approx(value),
         value = center_scale(value)) %>%
  spread(variable, value)

did_data <- abundance_indices %>%
  select(classcode, data) %>%
  unnest() %>%
  left_join(life_history_data, by = 'classcode') %>%
  left_join(enso, by = 'year') %>%
  left_join(pdo, by = 'year') %>%
  left_join(annual_conditions, by = c('year','region')) %>%
  left_join(ci_catches, by = c('classcode','year')) %>%
  mutate(catch = ifelse(is.na(catch), 0, catch))

did_terms <- did_data %>%
  select(year, targeted) %>%
  mutate(index = 1:nrow(.),
         targeted = as.numeric(targeted == 'Targeted')) %>%
  spread(year, targeted, fill = 0) %>%
  select(-index) %>%
  set_names(., paste0('did_',colnames(.))) %>%
  select(-did_2002)

did_data <- did_data %>%
  mutate(targeted = as.numeric(targeted == 'Targeted'),
         post_mpa = as.numeric(year >= 2003)) %>%
  bind_cols(did_terms)


did_reg <- paste0('abundance_index ~',paste(c('targeted','post_mpa','mean_enso','mean_pdo','(mean_annual_temp|classcode)','mean_annual_kelp',colnames(did_terms)), collapse = '+'))


did_models <- did_data %>%
  nest(-data_source,-abundance_source) %>%
  mutate(did_reg = did_reg) %>%
  mutate(did_model = map2(data, did_reg, ~lme4::glmer(.y, data = .x)))


did_plot_foo <- function(x) {
  x %>%
    broom::tidy() %>%
    filter(str_detect(term, 'did')) %>%
    mutate(year = str_replace(term,'did_','') %>% as.numeric()) %>%
    ggplot() +
    geom_pointrange(aes(year, estimate, ymin = estimate - 1.96 *std.error, ymax = estimate + 1.96 * std.error)) +
    geom_vline(aes(xintercept = 2003), color = 'red', linetype = 2) +
    geom_hline(aes(yintercept = 0))
}

did_models <- did_models %>%
  mutate(did_plot = map(did_model, did_plot_foo))

did_plot_foo <- function(data_source, abundance_source,did_plot, run_dir){

  ggsave(file = paste0(run_dir, '/',data_source,'-',abundance_source,'-did.pdf'), did_plot)
}

pwalk(list(data_source = did_models$data_source,
           abundance_source =  did_models$abundance_source,
           did_plot = did_models$did_plot),did_plot_foo, run_dir = run_dir)

