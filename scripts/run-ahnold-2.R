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

run_vast <- T # run VAST, best to leave off for now

num_knots <-  10

aggregate_transects  <-
  F # should transects be aggregated up (mean across observer) or left as raw

channel_islands_only <- T # only include channel islands, leave T

min_year <- 1999 # included years must be greater than this

occurance_ranking_cutoff <- 0.5

small_num <-  0 # no idea

use_mpa_site_effects <- F # no idea

# load data ---------------------------------------------------------------

length_data <- read_csv('data/UCSB_FISH raw thru 2013.csv') %>%
  magrittr::set_colnames(., tolower(colnames(.))) %>%
  mutate(classcode = tolower(classcode))

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
  left_join(life_history_data, by = 'classcode') #%>%
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
      classcode == (density_example$classcode)
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

  save(file = paste0('data/length-to-density-data.Rdata'),
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


observer_experience <- length_to_density_data %>%
  group_by(year,observer) %>%
  summarise(n_obs = length(log_density)) %>%
  arrange(observer, year) %>%
  group_by(observer) %>%
  mutate(cumulative_n_obs = cumsum(n_obs))

# a <- length_to_density_data %>%
#   left_join(observer_experience, by = c('observer','year')) %>%
#   ungroup() %>%
#   left_join(life_history_data %>% select(classcode,vbgf.linf), by = 'classcode')
#
# b <- a %>%
#   group_by(classcode,year) %>%
#   summarise(nobs = sum(any_seen, na.rm = T))
#
# b <- glm(any_seen ~ vbgf.linf, data = a, family = 'binomial')

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
    'mean_vis',
    'mean_canopy'
  ),
  collapse = '+')


prob_raw_length_covars <-
  paste(c('region','mean_vis', 'mean_canopy'),
        collapse = '+')

supplied_density_covars <-
  paste(c('region','mean_kelp', 'mean_vis'),
        collapse = '+')

prob_supplied_density_covars <-
  paste(c('region','mean_kelp', 'mean_vis'),
        collapse = '+')

length_to_density_models <- length_to_density_data %>%
  nest(-classcode) %>%
  mutate(ind_covars = raw_length_covars,
         prob_ind_covars = prob_raw_length_covars)

supplied_density_models <- density_data %>%
  nest(-classcode) %>%
  mutate(ind_covars = supplied_density_covars,
         prob_ind_covars = prob_supplied_density_covars)

length_to_density_models <- length_to_density_models %>%
  mutate(data_source = 'length_to_density')

supplied_density_models <- supplied_density_models %>%
  mutate(data_source = 'supplied_density')

# filter data -------------------------------------------------------------

nobs_quantiles <- length_to_density_data %>%
  filter(any_seen == T) %>%
  group_by(year, classcode) %>%
  summarise(nobs = sum(any_seen)) %>%
  ungroup() %>%
  {quantile(.$nobs)}

well_observed_species <- length_to_density_data %>%
  filter(year > 1999) %>%
  left_join(life_history_data %>% select(classcode,commonname, targeted), by = c('classcode')) %>%
  group_by(year, classcode,commonname, targeted) %>%
  summarise(nseen = sum(any_seen, na.rm = T)) %>%
  group_by(commonname,classcode, targeted) %>%
  summarise(min_seen = min(nseen, na.rm = T)) %>%
  arrange(desc(min_seen)) %>%
  ungroup() %>%
  filter(min_seen > nobs_quantiles[2]) %>%
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
  filter(classcode %in% well_observed_species$classcode) %>%
  mutate(data = map(data, ~ left_join(.x, site_data, by = c('site', 'side')))) %>%
  mutate(data = map2(data,classcode, filterfoo, min_year = min_year,
                    filter_level = quo(classcode))) %>% # filter out things
  mutate(dim_data = map_dbl(data, nrow)) %>%
  filter(dim_data > 0)


# prepare data for abundance estimates aggregate, augment, center-scale data ----------------------------------------------

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
  filter(str_detect(commonname, 'YOY') == F, is.na(targeted) == F) #%>%
  # mutate(data = map(data, ~ purrrlyr::dmap_if(.x, is.numeric, center_scale))) # center and scale continuos data


population_structure <- c('one-pop', 'regional-pops', 'mpa-pops')

population_filtering <- c('all', 'mpa-only', 'consistent-sites')

abundance_models <- cross_df(list(data = list(abundance_models),
         population_structure = population_structure,
         population_filtering = population_filtering)) %>%
  unnest()



# run vast ----------------------------------------------------------------



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

if (run_vast == T) {
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


# test aggreate model -----------------------------------------------------

# joint_data <- abundance_models %>%
#   filter(population_structure == 'one-pop',
#          population_filtering == 'all',
#          data_source == 'length_to_density') %>%
#   unnest()
#
# afla <- joint_data %>%
#   filter(classcode == 'afla') %>%
#   group_by(factor_year) %>%
#   summarise(nobs = sum(any_seen),
#             not_seen = sum(!any_seen))
#
# observed <- glm('log_density ~ factor_year*classcode + observer + mean_vis',
#                 data = joint_data %>% filter(any_seen))

# estimate abundance through delta-glm ------------------------------------

safely_fit_fish <- safely(fit_fish)

abundance_models <- abundance_models %>%
  # filter(commonname == 'blue rockfish',
  #        population_structure == 'one-pop',
  #        population_filtering == 'all',
  #        data_source  == 'length_to_density') %>%
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


abundance_models <- abundance_models %>%
  mutate(no_error = map2_lgl(seeing_error, seen_error, ~ is.null(.x) &
                               is.null(.y))) %>%
  filter(no_error == T) #filter out models that didn't converge for some reason


# abundance_models %>%
#   group_by(targeted) %>%
#   summarise(nobs = length(no_error))

abundance_models <- abundance_models %>%
  mutate(
    seen_coefs = map(seen_model, broom::tidy),
    seen_aug = map(seen_model, broom::augment),
    seeing_coefs = map(seeing_model, broom::tidy),
    seeing_aug = map(seeing_model, broom::augment)
  )

add_covars_foo <- function(data){

  if(is.null(data$region)){data$region <- 'somewhere'}

  if(is.null(data$eventual_mpa)){data$eventual_mpa <- 'maybe'}

return(data)
}

abundance_models <- abundance_models %>%
  mutate(seen_aug = map(seen_aug, add_covars_foo),
         seeing_aug = map(seeing_aug, add_covars_foo))

abundance_models <- abundance_models %>%
  mutate(seen_rank_deficient = map_lgl(seen_model, ~ length(.$coefficients) > .$rank),
         seeing_rank_deficient = map_lgl(seeing_model, ~ length(.$coefficients) > .$rank))

wtf <- abundance_models %>%
  filter(seen_rank_deficient == T)


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
      abundance_plot
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
    data$mean_biomass_g <- data$biomass

  }

  if (population_structure == 'one-pop'){

  raw_trend <- data %>%
    group_by(factor_year) %>%
    summarise(abundance_index = mean(mean_biomass_g, na.rm = T)) %>%
    ungroup() %>%
    mutate(abundance_index = center_scale(abundance_index + 1e-3)) %>%
    mutate(year = factor_year %>% as.character() %>% as.numeric()) %>%
    ungroup()
  }
  if (population_structure == 'regional-pops'){

    raw_trend <- data %>%
      group_by(region,factor_year) %>%
      summarise(abundance_index = mean(mean_biomass_g, na.rm = T)) %>%
      ungroup() %>%
      mutate(abundance_index = center_scale(abundance_index + 1e-3)) %>%
      mutate(year = factor_year %>% as.character() %>% as.numeric()) %>%
      ungroup()
  }
  if (population_structure == 'mpa-pops'){

    raw_trend <- data %>%
      group_by(eventual_mpa,factor_year) %>%
      summarise(abundance_index = mean(mean_biomass_g, na.rm = T)) %>%
      ungroup() %>%
      mutate(abundance_index = center_scale(abundance_index + 1e-3)) %>%
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

save(file = paste0(run_dir, '/abundance_indices.Rdata'),
     abundance_indices)



# prepare candidate did runs ----------------------------------------------


# annual_conditions <- conditions_data %>%
#   left_join(site_data %>% select(site, region), by = 'site') %>%
#   group_by(region, year) %>%
#   summarise(
#     mean_annual_temp = mean(mean_temp, na.rm = T),
#     mean_annual_kelp = mean(mean_kelp, na.rm  = T)
#   ) %>%
#   gather(variable, value,-year,-region) %>%
#   group_by(variable) %>%
#   mutate(value = zoo::na.approx(value),
#          value = center_scale(value)) %>%
#   spread(variable, value)

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
  unnest() %>%
  left_join(life_history_data, by = 'classcode') %>%
  left_join(enso, by = 'year') %>%
  left_join(pdo, by = 'year') %>%
  left_join(annual_conditions, by = c('year')) %>%
  left_join(ci_catches, by = c('classcode', 'year')) %>%
  mutate(catch = ifelse(is.na(catch), 0, catch)) %>%
  mutate(targeted = as.numeric(targeted == 'Targeted'),
         post_mpa = as.numeric(year >= 2003))


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
  mutate(index = 1:nrow(.),
         targeted = as.numeric(targeted == 'Targeted')) %>%
  spread(year, targeted, fill = 0) %>%
  select(-index) %>%
  set_names(., paste0('did_', colnames(.))) %>%
  select(-did_2002)

# fit DiD estimator on abundance indicies ---------------------------------
# so this is where you need to do a bunch of work in setting up the model runs







did_reg <-
  paste0('abundance_index ~', paste(
    c(
      'targeted',
      'post_mpa',
      'mean_enso',
      'mean_pdo',
      '(mean_annual_temp|classcode)',
      'mean_annual_kelp',
      colnames(did_terms)
    ),
    collapse = '+'
  ))

did_models <- did_data %>%
  nest(-data_source, -abundance_source) %>%
  mutate(did_reg = did_reg) %>%
  mutate(did_model = map2(data, did_reg, ~ lme4::glmer(.y, data = .x)))


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
    geom_vline(aes(xintercept = 2003),
               color = 'red',
               linetype = 2) +
    geom_hline(aes(yintercept = 0))
}

did_models <- did_models %>%
  mutate(did_plot = map(did_model, did_plot_foo))

did_plot_foo <-
  function(data_source,
           abundance_source,
           did_plot,
           run_dir) {
    ggsave(file = paste0(run_dir, '/', data_source, '-', abundance_source, '-did.pdf'),
           did_plot)
  }

pwalk(
  list(
    data_source = did_models$data_source,
    abundance_source =  did_models$abundance_source,
    did_plot = did_models$did_plot
  ),
  did_plot_foo,
  run_dir = run_dir
)
