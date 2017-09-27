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
library(purrr)
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
library(trelliscopejs)


demons::load_functions()

# set options ------------------------------
# Summary: set options for model run

run_name <- 'Working'

run_dir <- file.path('results', run_name)

run_description <-
  'Model selection process, testing STAN selection'

if (dir.exists(run_dir) == F) {
  dir.create(run_dir)
}

write(run_description,
      file = paste(run_dir, 'RUN_DESCRIPTION.txt', sep = '/'))

# set run parameters ------------------------------
# Summary: set things like filters

channel_islands_only <- T

min_year <- 1999

occurance_ranking_cutoff <- 0.5

small_num <-  0

use_mpa_site_effects <- F

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
  # rename(description_2 = Description) %>%
  magrittr::set_colnames(., tolower(colnames(.)))

site_data <- read_csv('data/Final_Site_Table_UCSB.csv') %>%
  magrittr::set_colnames(., tolower(colnames(.))) %>%
  select(site,side,mpagroup, mpa_status, reserve, region, year_mpa,mpaareanm2) %>%
  unique()

length_data <- length_data %>%
  left_join(life_history_data %>% mutate(classcode = toupper(classcode)), by = 'classcode') %>%
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
    summarise(mean_pdo = mean(pdo, na.rm = T)) %>%
    mutate(lag1_pdo = dplyr::lag(mean_pdo,1),
           lag2_pdo = dplyr::lag(mean_pdo,2),
           lag3_pdo = dplyr::lag(mean_pdo,3),
           lag4_pdo = dplyr::lag(mean_pdo,4))

} else {
  scrape_pdo(outdir = 'data/')

  pdo <- read_csv('data/pdo.csv') %>%
    group_by(year) %>%
    summarise(mean_pdo = mean(pdo, na.rm = T)) %>%
    mutate(lag1_pdo = dplyr::lag(mean_pdo,1),
           lag2_pdo = dplyr::lag(mean_pdo,2),
           lag3_pdo = dplyr::lag(mean_pdo,3),
           lag4_pdo = dplyr::lag(mean_pdo,4))

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
    omit_names = c('log_density','biomass', 'year', 'mean_enso', 'mean_pdo',
                   'targeted','year_mpa',paste0('lag',1:4,'_enso'),paste0('lag',1:4,'_pdo'))
  ) %>%
  mutate(log_density = log(biomass + small_num)) %>%
  purrrlyr::by_row(function(x,y) any(is.na(x[,y])), y = reg_vars) %>%
  filter(.out == F)


# reg_data <- reg_data %>%
#   group_by(classcode) %>%
#   mutate(log_density = center_scale(log_density,'log_density'))
#
# a %>%
#   group_by(year)

# cs_reg_data <- reg_data


if (use_mpa_site_effects == F)
{

reg_data <- reg_data %>%
  mutate(did_dummy = targeted,
         did_year = paste('did', year, sep = '_')) %>%
  spread(did_year, did_dummy, fill = 0) %>%
  mutate(temp2 = mean_temp ^ 2,
         pdo2 = mean_pdo ^ 2,
         enso2 = mean_enso ^ 2,
         site_side = paste(site,side,sep = '_'),
         factor_year = as.factor(year))

} else {
reg_data <- reg_data %>%
  mutate(did_dummy = targeted * (year_mpa >0),
         did_year_inside = paste('did', year, 'inside',sep = '_')) %>%
  spread(did_year_inside, did_dummy, fill = 0) %>%
  mutate(did_dummy = targeted * (year_mpa == 0),
         did_year_outside = paste('did', year,'outside', sep = '_')) %>%
  spread(did_year_outside, did_dummy, fill = 0) %>%
  mutate(temp2 = mean_temp ^ 2,
         pdo2 = mean_pdo ^ 2,
         enso2 = mean_enso ^ 2,
         site_side = paste(site,side,sep = '_'),
         factor_year = as.factor(year))

}


# filter data ------------------------------
# Summary: Apply all filters to data here for clarity

# Evaluate species occurances

species_sample_counts = reg_data %>%
  group_by(classcode) %>%
  summarise(num_samples = sum(is.na(biomass) == F & biomass > 0)) %>%
  ungroup() %>%
  mutate(prank = percent_rank(num_samples)) %>%
  filter(prank > occurance_ranking_cutoff) %>%
  arrange(prank %>% desc()) %>%
  left_join(life_history_data %>% select(classcode, targeted), by = 'classcode')

fished_balance_plot <- species_sample_counts %>%
  ggplot(aes(targeted)) +
  geom_bar()

reg_data <- reg_data %>%
  dplyr::filter(is.na(biomass) == F,
         is.na(commonname) == F,
         is.na(targeted) == F,
         year > min_year,
         !str_detect(commonname %>% tolower(), 'yoy'),
         region != 'SBI',
         classcode %in% species_sample_counts$classcode)

if (channel_islands_only == T){
  reg_data <- reg_data %>%
    filter( !region %in% c('SCSR','CCSR'))
}

reg_data$targeted[reg_data$commonname == 'black surfperch'] <-  1

seen_reg_data <- reg_data %>%
  filter(any_seen == T)

if (use_mpa_site_effects == F){

did_years <-
  paste('did', min(seen_reg_data$year): max(seen_reg_data$year), sep = '_')

# did_year <- did_years[did_years != 'did_2002']

} else {

did_years <-
  paste('did', min(seen_reg_data$year): max(seen_reg_data$year), sep = '_')

did_years_inside <- paste(did_years, 'inside', sep = '_')

did_years_outside <- paste(did_years, 'outside', sep = '_')

did_years <- c(did_years_inside, did_years_outside)


}
did_year <- did_years[str_detect(did_years,'2003') == F]
# did_year <- did_years[did_years != 'did_2002']


# model selection ------------------------------
# Summary: Run and compare a variety of model structures

reg_fmla_1 <-
  as.formula(
    paste0(
      "log_density ~",
      paste(did_year, collapse = "+"),
      " + (1|year) + (1 + mean_temp + temp2 + region|classcode)+ mean_enso + lag1_enso  + lag2_enso +
      lag3_enso + lag4_enso + mean_pdo + lag1_pdo + lag2_pdo + lag3_pdo + lag4_pdo + (1 |region:site) + targeted"
    )
  )

reg_fmla_2 <-
  as.formula(
    paste0(
      "log_density ~",
      paste(did_year, collapse = "+"),
      " + factor_year + (1 + mean_temp + temp2 + region|classcode)+ mean_enso + lag1_enso  + lag2_enso +
      lag3_enso + lag4_enso + mean_pdo + lag1_pdo + lag2_pdo + lag3_pdo + lag4_pdo + (1 |region:site) + targeted"
    )
    )




reg_fmla_4 <-
  as.formula(
    paste0(
      "log_density ~",
      paste(did_year, collapse = "+"),
      "+ (1|year) + (1 + mean_temp + temp2 + region|classcode)+ mean_enso + mean_pdo + (1 |region:site)"
    )
    )

reg_fmla_5 <-
  as.formula(
    paste0(
      "log_density ~",
      paste(did_year, collapse = "+"),
      "+ (1|year) + (1 + mean_temp + temp2 + region|classcode)+ mean_enso + (1 |region:site)"
    )
  )

# reg_fmla_12 <-
#   as.formula(
#     paste0(
#       "log_density ~",
#       paste(did_year, collapse = "+"),
#       "+ (1|year) + (1 + mean_temp + temp2 + region|classcode)"
#     )
#   )


# wtf <- lme4::lmer(reg_fmla_12, seen_reg_data)


reg_fmla_6 <-
  as.formula(
    paste0(
      "log_density ~",
      paste(did_year, collapse = "+"),
      "+ (1|year) + (1 + mean_temp + temp2 + region|classcode)+ mean_pdo + (1 |region:site) + targeted"
    )
  )

reg_fmla_7 <-
  as.formula(
    paste0(
      "log_density ~",
      paste(did_year, collapse = "+"),
      "+ (1|year) + (1 + region|classcode)+ mean_pdo + mean_temp + (1 |region:site) + targeted"
    )
  )

reg_fmla_8 <-
  as.formula(
    paste0(
      "log_density ~",
      paste(did_year, collapse = "+"),
      "+ (1|year) + (1 + region|classcode)+ mean_enso + lag1_enso  + lag2_enso +
      lag3_enso + lag4_enso + mean_pdo + lag1_pdo + lag2_pdo + lag3_pdo + lag4_pdo + (1 |region:site) + targeted"
    )
    )

reg_fmla_9 <-
  as.formula(
    paste0(
      "log_density ~",
      paste(did_year, collapse = "+"),
      "+ (1|year) + (1 + region|classcode)+ mean_enso + lag1_enso  + lag2_enso +
      lag3_enso + lag4_enso + mean_pdo + lag1_pdo + lag2_pdo + lag3_pdo + lag4_pdo + (1 |region:site) + targeted"
    )
  )

reg_fmla_10 <-
  as.formula(
    paste0(
      "log_density ~",
      paste(did_year, collapse = "+"),
      "+ (1|year) + (1 + mean_temp + temp2 + region|classcode)+ mean_pdo + (1 |site) + targeted"
    )
  )

reg_fmla_11 <-
  as.formula(
    paste0(
      "log_density ~",
      paste(did_year, collapse = "+"),
      "+ (1|year) + (1 + mean_temp + temp2 + region|classcode)+ mean_pdo + (1 |region) + targeted"
    )
  )



regs <- ls()[str_detect(ls(),"reg_fmla")]

reg_list <- map(regs, get) %>%
  set_names(regs)


validate <- seen_reg_data %>%
  crossv_mc(1, test = 0.05)

test_data <- list(test_training = list(validate), reg_fmlas = reg_list)


test_data <- cross_d(test_data) %>%
  mutate(reg_name = regs) %>%
  unnest(test_training, .drop = F)

canditate_models <- test_data %>%
  mutate(fitted_model = map2(reg_fmlas, train,~lme4::lmer(.x, data = .y)))


canditate_models <- canditate_models %>%
  mutate(root_mean_se = map2_dbl(fitted_model, train , ~modelr::rmse(.x,.y)),
         aic = map_dbl(fitted_model, AIC),
         bic = map_dbl(fitted_model, BIC),
         r2 = map2_dbl(fitted_model, train, ~rsquare(.x,.y))) %>%
  arrange(aic)

# prep regression ------------------------------
# Summary: select variables, forms, etc.
#

# reg <-
#   as.formula(
#     paste0(
#       "log_density ~",
#       paste(did_year, collapse = "+"),
#       "+ factor_year + mean_temp + temp2 + mean_pdo + targeted"
#     )
#   )
#
#
reg <-
  as.formula(
    paste0(
      "log_density ~",
      paste(did_year, collapse = "+"),
      "+ (1|year) + ( mean_temp + temp2 + region|classcode)+ mean_pdo + (1 |site) + targeted"
    )
  )


reg <- canditate_models$reg_fmlas[[1]]
# reg <-
#   as.formula(
#     paste0(
#       'log_density ~',
#       paste(did_year, collapse = '+'),
#       ' + (1|year) + (1 + mean_temp + temp2 + region|classcode)+ mean_enso + lag1_enso  + lag2_enso +
#       lag3_enso + lag4_enso + mean_pdo + lag1_pdo + lag2_pdo + lag3_pdo + lag4_pdo + (1 |region:site) + targeted'
#     )
#   )

# gamma_reg <-
#   as.formula(
#     paste0(
#       'biomass ~',
#       paste(did_year, collapse = '+'),
#       ' + (1|year) + (1 + mean_temp + temp2 + region|classcode)+ mean_enso + lag1_enso  + lag2_enso +
#       lag3_enso + lag4_enso + mean_pdo + lag1_pdo + lag2_pdo + lag3_pdo + lag4_pdo + (1 |region:site) + targeted'
#     )
#     )

# reg <-
#   as.formula(
#     paste0(
#       'log_density ~',
#       paste(did_year, collapse = '+'),
#       ' + (1|year) + (1 + mean_temp + temp2 + region |classcode)+ mean_enso + lag1_enso  + lag2_enso +
#       lag3_enso + lag4_enso + mean_pdo + (site | region) + targeted'
#     )
#     )


# posreg <-
#   as.formula(
#     paste0(
#       "biomass ~",
#       paste(did_year, collapse = '+'),
#       ' + (1|year) + (1 + mean_temp + temp2 |classcode)+ mean_enso + lag1_enso  + lag2_enso +
#       lag3_enso + lag4_enso + mean_pdo  + (1 | site) + region  + targeted'
#     )
#     )

# ' + (1|year) + (1 + mean_temp + temp2 |classcode) + mean_enso +  mean_pdo + (1 | site_side) + (1 | site) + (1 | region) + targeted + post_mlpa'


# run model ------------------------------
# Summary: fit ahnold model

seen_model <- lm(reg, data = seen_reg_data)


seen_model <- lme4::lmer(reg, data = seen_reg_data)


# seen_model <- lme4::glmer(posreg, data = seen_reg_data, family = Gamma(link = ''))

# stan_seen_model <- stan_glmer(reg, data = seen_reg_data, cores = 4, chains = 4, iter = 4000)


consistent_sites <-  seen_reg_data %>%
  group_by(site_side) %>%
  summarise(spans_transition = all(2001:2005 %in% unique(year))) %>%
  filter(spans_transition == T)




seen_ana_sci_model <- lme4::lmer(reg, data = seen_reg_data %>%
                                   filter(site_side %in% consistent_sites$site_side))
                                # filter(region %in% c('ANA','SCI')))

seen_mpa_model <- lme4::lmer(reg, data = seen_reg_data %>%
                                   filter(year_mpa >0))

seen_nompa_model <- lme4::lmer(reg, data = seen_reg_data %>%
                               filter(year_mpa <=0))

unseen_model <- lme4::glmer(any_seen ~ mean_temp + mean_kelp + mean_vis + (1 | site_side) +
                              (1 | year) + (1 | classcode), family = binomial(link = 'logit'), data = reg_data)


# diagnose model ------------------------------
# Summary: run model diagnostics


density_cor <- reg_data %>%
  select(year, site_side,classcode,biomass) %>%
  spread(classcode, biomass) %>%
  select(-year,-site_side) %>%
  cor() %>%
  as_data_frame() %>%
  mutate(classcode = colnames(.)) %>%
  gather(classcode_2,correlation, -classcode) %>%
  left_join(life_history_data %>% select(classcode, targeted), by = 'classcode')


species_correlation_plot <- density_cor %>%
  ggplot() +
  geom_tile(aes(classcode, classcode_2, fill = correlation)) +
  scale_fill_viridis()

# reg_data %>%
#   group_by(year, targeted) %>%
#   summarise(mb = mean(biomass)) %>%
#   spread(targeted,mb) %>%
#   ggplot(aes(`0`, `1`)) +
#   geom_point()


aug_seen_model <- seen_model %>%
  augment()

resid_v_fit_plot <- aug_seen_model %>%
  ggplot(aes(.fitted, .resid, color = site)) +
  geom_ref_line(h = 0, colour = 'red') +
  geom_point(alpha = 0.5, show.legend = F)

resid_hist_plot <- aug_seen_model %>%
  ggplot(aes(.resid)) +
  geom_histogram() +
  geom_ref_line(v = 0, colour = 'red')

qq_plot <- aug_seen_model %>%
  ggplot(aes(sample = .resid)) +
  geom_abline(aes(slope = 1, intercept = 0), color = 'red', linetype = 2) +
  stat_qq() +
  labs(y = 'Deviance Residuals', title = 'Normal QQ plot')

aug_unseen_model <- unseen_model %>%
  augment() %>%
  mutate(prob_seen = 1 / (1 + exp(-.fitted)),
         wtf = predict(unseen_model, type = 'response'))

unseen_plot <- aug_unseen_model %>%
  ggplot(aes(prob_seen, any_seen)) +
  geom_point()

# unseen_plot



# Run additional analyses -------------------------------------------------
# Summary:
#


length_hists <- length_data %>%
  group_by(year, classcode, fish_tl) %>%
  summarise(count = sum(count),
            comm_name = unique(commonname)) %>%
  ungroup() %>%
  nest(-comm_name) %>%
  mutate(samples = map_dbl(data, ~sum(.x$count))) %>%
  mutate(length_hist_plot = map_plot(
    data,
    ~ ggplot(.x, aes(fish_tl, count, color = year >= 2003)) + geom_line(show.legend = F) + facet_grid(year ~ ., as.table = F, scales = 'free_y') + theme_fivethirtyeight() +
      theme(axis.text.y = element_blank(), axis.title = element_text(),
            strip.text = element_text(size = 6), panel.grid.minor = element_blank())
  )) %>%
  select(comm_name,samples, length_hist_plot) %>%
  filter(is.na(comm_name) == F) %>%
  arrange(desc(samples))

trelliscopejs::trelliscope(length_hists,name = 'huh',
                           panel_col ='length_hist_plot' )

fished_length_hists <- length_data %>%
  group_by(year, targeted, fish_tl) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  filter(is.na(targeted) == F) %>%
  ggplot(aes(fish_tl, count, color = year >= 2003)) +
  geom_line(show.legend = F) +
  facet_grid(year ~ targeted, as.table = F, scales = 'free_y') +
  theme_fivethirtyeight() +
      theme(axis.text.y = element_blank(), axis.title = element_text(),
            strip.text = element_text(size = 6), panel.grid.minor = element_blank())



trophic_effects <- seen_reg_data %>%
  nest(-trophicgroup) %>%
  mutate(ahnold_model = map(data, ~lme4::lmer(reg, data = .x)))

plot_did <- function(model,trophicgroup){

  mpa_effect_plot <-  model %>%
    tidy() %>%
    mutate(lower = estimate - 1.96 * std.error,
           upper = estimate + 1.96 * std.error) %>%
    filter(str_detect(term,'did')) %>%
    mutate(year = str_replace(term, 'did_','') %>% as.numeric()) %>%
    ggplot() +
    geom_hline(aes(yintercept = 0)) +
    geom_vline(aes(xintercept = 2003), color = 'red', linetype = 2, size = 2) +
    geom_pointrange(aes(year,estimate, ymax = upper, ymin = lower),color = 'skyblue4', size = 2) +
    geom_pointrange(data = data_frame(year = 2002, estimate = 0), aes(year, estimate,ymin = estimate, ymax = estimate),color = 'skyblue4', size = 2) +
    ylab('Estimated MLPA Effect') +
    xlab('Year') +
    coord_cartesian(ylim = c(-1.25,1.25))
}

trophic_effects <- trophic_effects %>%
  mutate(did_plot = map2_plot(ahnold_model, trophicgroup, ~plot_did(.x,.y)))
#
trophic_effects %>%
  select(trophicgroup, did_plot) %>%
  trelliscopejs::trelliscope(name = 'DiD Effects by Trophic Group')

# huh <- life_history_data %>%
#   filter(trophicgroup == 'benthic micro-invertivore') %>%
#   select(commonname, targeted)




# figures ------------------------------
# Summary: make summary figures

with_mlpa <- seen_model %>% augment()

without_mlpa <- seen_model %>% augment()

without_mlpa[, str_detect(colnames(without_mlpa),'did_')] <-  0

with_mlpa <- with_mlpa %>%
  add_predictions(model = seen_model) %>%
  mutate(world = 'With MLPA')

without_mlpa <- without_mlpa %>%
  add_predictions(model = seen_model) %>%
  mutate(world = 'Without MLPA')

mlpa_experiment <- with_mlpa %>%
  bind_rows(without_mlpa) %>%
  ggplot(aes((year),pred %>% exp(), color = world, fill = world)) +
  geom_smooth()

with_mlpa_mpas <- seen_mpa_model %>% augment()

without_mlpa_mpas <- seen_mpa_model %>% augment()

without_mlpa_mpas[, str_detect(colnames(without_mlpa_mpas),'did_')] <-  0

with_mlpa_mpas <- with_mlpa_mpas %>%
  add_predictions(model = seen_mpa_model) %>%
  mutate(world = 'With MLPA')

without_mlpa_mpas <- without_mlpa_mpas %>%
  add_predictions(model = seen_mpa_model) %>%
  mutate(world = 'Without MLPA')

mlpa_mpa_experiment <- with_mlpa_mpas %>%
  bind_rows(without_mlpa_mpas) %>%
  ggplot(aes((year),pred %>% exp(), color = world, fill = world)) +
  geom_smooth()

with_mlpa_no_mpas <- seen_nompa_model %>% augment()


without_mlpa_no_mpas <- seen_nompa_model %>% augment()

without_mlpa_no_mpas[, str_detect(colnames(without_mlpa_no_mpas),'did_')] <-  0

with_mlpa_no_mpas <- with_mlpa_no_mpas %>%
  add_predictions(model = seen_nompa_model) %>%
  mutate(world = 'With MLPA')

without_mlpa_no_mpas <- without_mlpa_no_mpas %>%
  add_predictions(model = seen_nompa_model) %>%
  mutate(world = 'Without MLPA')

mlpa_no_mpa_experiment <- with_mlpa_no_mpas %>%
  bind_rows(without_mlpa_no_mpas) %>%
  ggplot(aes((year),pred %>% exp(), color = world, fill = world)) +
  geom_smooth()


# mpa_site_effect_plot <-  seen_model %>%
#   tidy() %>%
#   mutate(lower = estimate - 1.96 * std.error,
#          upper = estimate + 1.96 * std.error) %>%
#   filter(str_detect(term,'did')) %>%
#   mutate(site_type = str_split(term,"_", simplify = T)[,3]) %>%
#   mutate(year = str_split(term,"_", simplify = T)[,2] %>% as.numeric()) %>%
#   ggplot() +
#   geom_hline(aes(yintercept = 0)) +
#   geom_vline(aes(xintercept = 2003), color = 'red', linetype = 2, size = 2) +
#   geom_pointrange(aes(year,estimate, ymax = upper, ymin = lower),color = 'skyblue4', size = 2) +
#   geom_pointrange(data = data_frame(year = 2002, estimate = 0), aes(year, estimate,ymin = estimate, ymax = estimate),color = 'skyblue4', size = 2) +
#   ylab('Estimated MLPA Effect') +
#   ggrepel::geom_text_repel(data = data_frame(x = 2003, y = 1), aes(x,y, label = 'MLPA Enacted'),nudge_x = 2) +
#   xlab('Year') +
#   coord_cartesian(ylim = c(-1.25,1.25)) +
#   labs(title = 'System-Wide Effect') +
#   facet_wrap(~site_type)


mpa_effect_plot <-  seen_model %>%
  tidy() %>%
  mutate(lower = estimate - 1.96 * std.error,
         upper = estimate + 1.96 * std.error) %>%
  filter(str_detect(term,'did')) %>%
  mutate(year = str_replace(term, 'did_','') %>% as.numeric()) %>%
  ggplot() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 2003), color = 'red', linetype = 2, size = 2) +
  geom_pointrange(aes(year,estimate, ymax = upper, ymin = lower),color = 'skyblue4', size = 2) +
  geom_pointrange(data = data_frame(year = 2003, estimate = 0), aes(year, estimate,ymin = estimate, ymax = estimate),color = 'skyblue4', size = 2) +
  ylab('Estimated MLPA Effect') +
  ggrepel::geom_text_repel(data = data_frame(x = 2003, y = 1), aes(x,y, label = 'MLPA Enacted'),nudge_x = 2) +
  xlab('Year') +
  coord_cartesian(ylim = c(-1.25,1.25)) +
  labs(title = 'System-Wide Effect')

# mpa_effect_plot


in_mpa_effect_plot <-  seen_mpa_model %>%
  tidy() %>%
  mutate(lower = estimate - 1.96 * std.error,
         upper = estimate + 1.96 * std.error) %>%
  filter(str_detect(term,'did')) %>%
  mutate(year = str_replace(term, 'did_','') %>% as.numeric()) %>%
  ggplot() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 2003), color = 'red', linetype = 2, size = 2) +
  geom_pointrange(aes(year,estimate, ymax = upper, ymin = lower), color = 'skyblue4', size = 2) +
  geom_pointrange(data = data_frame(year = 2002, estimate = 0), aes(year, estimate,ymin = estimate, ymax = estimate), color = 'skyblue4', size = 2) +
  ylab('Estimated MLPA Effect') +
  ggrepel::geom_text_repel(data = data_frame(x = 2003, y = 1), aes(x,y, label = 'MLPA Enacted'),nudge_x = 2) +
  xlab('Year') +
  coord_cartesian(ylim = c(-1.25,1.25)) +
  labs(title = 'Inside MPA Effect')

# in_mpa_effect_plot


mpa_anasci_effect_plot <-  seen_ana_sci_model %>%
  tidy() %>%
  mutate(lower = estimate - 1.96 * std.error,
         upper = estimate + 1.96 * std.error) %>%
  filter(str_detect(term,'did')) %>%
  mutate(year = str_replace(term, 'did_','') %>% as.numeric()) %>%
  ggplot() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 2003), color = 'red', linetype = 2, size = 2) +
  geom_pointrange(aes(year,estimate, ymax = upper, ymin = lower),color = 'skyblue4', size = 2) +
  geom_pointrange(data = data_frame(year = 2002, estimate = 0), aes(year, estimate,ymin = estimate, ymax = estimate),color = 'skyblue4', size = 2) +
  ylab('Estimated MLPA Effect') +
  ggrepel::geom_text_repel(data = data_frame(x = 2003, y = 1), aes(x,y, label = 'MLPA Enacted'),nudge_x = 2) +
  xlab('Year') +
  coord_cartesian(ylim = c(-1.25,1.25)) +
  labs(title = 'Anacapa and Santa Cruz Only')


# mpa_anasci_effect_plot

# save final results ------------------------------
# Summary: save final things
#
demons::save_plots(plot_dir = run_dir)

plots <- ls()[str_detect(ls(), "_plot")]

plot_list <- purrr::map(plots, get) %>%
  set_names(plots)

save(file = paste0(run_dir,"/data.Rdata"), reg_data, seen_reg_data)

save(file = paste0(run_dir,"/regressions.Rdata"), seen_model, seen_ana_sci_model, seen_mpa_model)

save(file = paste0(run_dir,"/plots.Rdata"), plot_list)

