# run_ahnold -------
# Author: Dan Ovando
# Project: ahnold
# Summary: Master script to run the effects of the MLPA
# project


# load libraries ------------------------------
# Summary: Load things you need
rm(list = ls())
set.seed(123)
library(tidyverse)
library(forcats)
library(modelr)
library(stringr)
library(rstan)
library(rstanarm)
library(car)
library(AER)
library(broom)
library(viridis)
library(scales)

demons::load_functions()

# set options ------------------------------
# Summary: set options for model run

run_name <- 0.51

run_dir <- paste('results', run_name, sep = '/')

run_description <-
  'Model development and testing. Focusing on hierarchichal
positive observations. Tackling logit part next. Trying to get the new structure up and running.
Will move to 1.0 when fully functional delta-glm. Adding in species-region interactions'

if (dir.exists(run_dir) == F) {
  dir.create(run_dir)
}

write(run_description,
      file = paste(run_dir, 'RUN_DESCRIPTION.txt', sep = '/'))

# set run parameters ------------------------------
# Summary: set things like filters

channel_islands_only <- T

min_year <- 1999

base_theme <- hrbrthemes::theme_ipsum(base_size = 18, axis_title_size = 16)

theme_set(base_theme)

# load Data ------------------------------
# Summary: Bring in required data for project

length_data <- read_csv('data/UCSB_FISH raw thru 2013.csv') %>%
  magrittr::set_colnames(., tolower(colnames(.)))

life_history_data <-
  read_csv('data/VRG Fish Life History in MPA_04_08_11_12 11-Mar-2014.csv') %>%
  rename(classcode = pisco_classcode) %>%
  mutate(classcode = tolower(classcode)) %>%
  rename(description_2 = Description) %>%
  magrittr::set_colnames(., tolower(colnames(.)))

site_data <- read_csv('data/Final_Site_Table_UCSB.csv') %>%
  magrittr::set_colnames(., tolower(colnames(.))) %>%
  select(site,side,mpagroup, mpa_status, reserve, region, year_mpa,mpaareanm2) %>%
  unique()

length_data <- length_data %>%
  left_join(life_history_data, by = 'classcode') %>%
  left_join(site_data, by = c('site', 'side'))

conditions_data <- length_data %>%
  group_by(site,side, year) %>%
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

if (file.exists('data/enso.csv')) {
  enso <- read_csv('data/enso.csv') %>%
    group_by(year) %>%
    summarise(mean_enso = mean(enso, na.rm = T)) %>%
    mutate(lag1_enso = dplyr::lag(mean_enso,1),
           lag2_enso = dplyr::lag(mean_enso,2),
           lag3_enso = dplyr::lag(mean_enso,3),
           lag4_enso = dplyr::lag(mean_enso,4))

} else {
  scrape_enso(outdir = 'data/')

}

if (file.exists('data/pdo.csv')) {
  pdo <- read_csv('data/pdo.csv') %>%
    group_by(year) %>%
    summarise(mean_pdo = mean(pdo, na.rm = T))

} else {
  scrape_pdo(outdir = 'data/')

}

# prepare data----------------------------
# Summary: apply transformations, calculations etc.

has_all <- function(x,reg_vars)
  any(is.na(x[reg_vars])) == F


reg_vars <- c('log_density', 'year','targeted', 'region' ,
              'mean_enso','mean_pdo', 'mean_temp','classcode','site','side','post_mlpa','region',
              'trophicgroup')

reg_data <- density_data %>%
  select(biomass, site,side,site_side, year, classcode) %>%
  # group_by(site,side, year, classcode) %>%
  # summarise(mean_density = mean(biomass, na.rm = T)) %>%
  ungroup() %>%
  left_join(conditions_data, by = c('site','side', 'year')) %>%
  left_join(life_history_data %>% select(classcode, targeted, trophicgroup,
                                         commonname),
            by = 'classcode') %>%
  left_join(enso, by = 'year') %>%
  left_join(pdo, by = 'year') %>%
  left_join(site_data %>% select(site,side,region,year_mpa), by = c('site','side')) %>%
  mutate(
    any_seen = biomass > 0,
    log_density = log(biomass),
    targeted = as.numeric(targeted == 'Targeted'),
    post_mlpa = as.numeric(year >= 2003)
  ) %>%
  filter(is.na(targeted) == F,
         is.na(year) == F,
         is.na(biomass) == F) %>%
  # select_(.dots = as.list(reg_vars)) %>%
  map2_df(
    colnames(.),
    center_scale,
    omit_names = c('log_density','mean_density', 'biomass', 'year', 'mean_enso', 'mean_pdo',
                   'targeted','year_mpa')
  ) %>%
  by_row(function(x,y) any(is.na(x[,y])), y = reg_vars) %>%
  filter(.out == F)

reg_data <- reg_data %>%
  mutate(did_dummy = targeted,
         did_year = paste('did', year, sep = '_')) %>%
  spread(did_year, did_dummy, fill = 0) %>%
  mutate(temp2 = mean_temp ^ 2,
         pdo2 = mean_pdo ^ 2,
         enso2 = mean_enso ^ 2,
         site_side = paste(site,side,sep = '_'))


# filter data ------------------------------
# Summary: Apply all filters to data here for clarity


reg_data <- reg_data %>%
  filter(is.na(biomass) == F,
         is.na(commonname) == F,
         is.na(targeted) == F,
         year > min_year,
         !str_detect(commonname %>% tolower(), 'yoy'))

if (channel_islands_only == T){
  reg_data <- reg_data %>%
    filter( !region %in% c('SCSR','CCSR'))
}

seen_reg_data <- reg_data %>%
  filter(any_seen == T)

# prep regression ------------------------------
# Summary: select variables, forms, etc.
#

did_years <-
  paste('did', min(seen_reg_data$year): max(seen_reg_data$year), sep = '_')

did_year <- did_years[did_years != 'did_2002']


reg <-
  as.formula(
    paste0(
      'log_density ~',
      paste(did_year, collapse = '+'),
      ' + (1|year) + (1 + mean_temp + temp2 + region|classcode)+ mean_enso + lag1_enso  + lag2_enso +
      lag3_enso + lag4_enso + mean_pdo  + (1 | site) + targeted'
    )
    )

reg <-
  as.formula(
    paste0(
      'log_density ~',
      paste(did_year, collapse = '+'),
      ' + (1|year) + (1 + mean_temp + temp2 |classcode)+ mean_enso + lag1_enso  + lag2_enso +
      lag3_enso + lag4_enso + mean_pdo  + (site | region) + (1 | region:classcode) + targeted'
    )
    )


posreg <-
  as.formula(
    paste0(
      'biomass ~',
      paste(did_year, collapse = '+'),
      ' + (1|year) + (1 + mean_temp + temp2 |classcode)+ mean_enso + lag1_enso  + lag2_enso +
      lag3_enso + lag4_enso + mean_pdo  + (1 | site) + region  + targeted'
    )
    )

# ' + (1|year) + (1 + mean_temp + temp2 |classcode) + mean_enso +  mean_pdo + (1 | site_side) + (1 | site) + (1 | region) + targeted + post_mlpa'


# run model ------------------------------
# Summary: fit ahnold model


seen_model <- lme4::lmer(reg, data = seen_reg_data)


# seen_model <- lme4::glmer(posreg, data = seen_reg_data, family = Gamma(link = ''))

stan_seen_model <- stan_glmer(reg, data = seen_reg_data)

save(file = 'stan_check.Rdata', stan_seen_model, seen_reg_data)
