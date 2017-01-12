
# Set up ------------------------------------------------------------------
# This script runs the development version of the MLPA effects model

rm(list = ls())
set.seed(666)
library(tidyverse)
library(demons)
library(knitr)
library(gridExtra)
library(broom)
library(coda)
library(ggmcmc)
library(LaplacesDemon)
library(scales)
library(AER)
library(mvtnorm)
library(stringr)
library(grid)
library(fishMod)
library(runjags)
library(rjags)
library(rstanarm)
library(forcats)

demons::load_functions('functions')


# Run Options -------------------------------------------------------------


runfolder <- '5.0'

description <- 'versions 5.0 signify start of conversion to STAN'


agg_level <- 'site'

MPA_Only <- TRUE

scale_numerics <- T

its <- 1e6

run_mcmc <- T

model_form <- 'Summon Demon'

runpath <- paste('results/',runfolder,'/', sep = '')

if (dir.exists(runpath) == F)
{
  dir.create(runpath, recursive = T)
}

write(description, file = paste(runpath,'RUN_DESCRIPTION.txt'))


# Load Data ---------------------------------------------------------------


rawdat <- read.csv('data/UCSB_FISH raw thru 2013.csv', stringsAsFactors = F)


rawdat %>%
  group_by(observer) %>%
  summarise(num_samps = length(fish_tl)) %>%
  mutate(observer = fct_reorder(observer, num_samps)) %>%
  ggplot(aes(observer, num_samps)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  theme(text = element_text(size = 8))

rawdat %>%
  ggplot(aes(pctcnpy,count)) +
  geom_point()

life.history <- read.csv('data/VRG Fish Life History in MPA_04_08_11_12 11-Mar-2014.csv', stringsAsFactors = F) %>%
  rename(classcode = pisco_classcode)

species.dat <- read.csv('data/master_spp_table.csv', stringsAsFactors = F)

site.dat <- read.csv('data/Final_Site_Table_UCSB.csv', stringsAsFactors = F) %>%
  rename(site = SITE)

mlpa_dat <- left_join(rawdat,life.history, by = 'classcode') %>%
  left_join(site.dat, by = 'site') %>%
  left_join(species.dat, by = 'classcode' )


conditions_dat <- mlpa_dat %>%
  group_by(year,SITE_SIDE) %>%
  summarise(mean_temp = mean(temp, na.rm = T),
            num_transects = length(site),
            mean_depth = mean(depth, na.rm = T),
            #             mean_surge = mean(surge, na.rm = T), # need to convert this to numeric
            mean_vis = mean(vis, na.rm = T),
            mpa_group = unique(MPAGROUP),
            pct_kelp_canopy = mean(pctcnpy, na.rm = T)
  ) %>%
  rename(site_side = SITE_SIDE)

# Processed PISCO data ----
# These are the processed biomass density estimates metric tons / hectare for is "biomass"

# Tidy data
# note that biomass is in units of mt/hectare and fish is un # per hectare


processed_site_dat <- read.csv('data/ci_reserve_data_final3 txt.csv', stringsAsFactors = F) %>%
  gather('concat.name','value', grep('_',colnames(.)),convert = T) %>%
  mutate(data.type = gsub('\\_.*', '', concat.name),
         classcode = gsub('.*\\_','',concat.name)) %>%
  mutate(value = as.numeric(value)) %>%
  spread(data.type,value) %>%
  rename(site_side = site.side) %>%
  left_join(conditions_dat,by = c('year','site_side')) %>%
  left_join(life.history, by = 'classcode') %>%
  rename(description2 = Description) %>%
  subset(is.na(biomass) == F)

colnames(processed_site_dat) <- tolower(colnames(processed_site_dat)) #make all variables lower case

capfun <- function(x,cap1=0,cap2 = 3000)
{
  x[x == cap1] <- cap2

  return(x)
}

processed_dat <- processed_site_dat %>%
  rename(common_name = commonname) %>%
  group_by(site_side, common_name) %>%
  mutate(ever.seen = mean(biomass, na.rm = T)>0) %>%
  subset(year != 1999 & ever.seen == T & grepl('YOY',classcode) == F  & grepl('YOY',common_name) == F) %>% #remove 1999, species never seen at a site, and young of the year. Why remove YOY again?
  ungroup() %>%
  mutate(year_mlpa_mpa = year.mpa * as.numeric(year.mpa >= 2003)) %>%
  group_by(site_side) %>%
  mutate(num.years.surveyed = length(unique(year)),
         will_be_mlpa = year_mlpa_mpa>0) %>%
  subset(num.years.surveyed > 2) %>% #Only keep sites that have been surveyed two or more years
  ungroup() %>%
  mutate(year_mpa2 = capfun(year_mlpa_mpa)) %>%
  group_by(region,year) %>%
  mutate(region_has_mpas = any(year_mpa2 <= year)) %>%
  ungroup() %>%
  group_by(region) %>%
  mutate(year_region_has_mlpa_mpas = min(year[year >= year_mpa2])) %>% #redefining 'MPA' as an MLPA created MPA
  ungroup() %>%
  group_by(common_name) %>%
  mutate(total_species_samples = sum(biomass>0, na.rm = T)) %>%
  ungroup() %>%
  subset(total_species_samples > quantile(total_species_samples,.1))

save(mlpa_dat,processed_dat, file = paste(runpath,'Raw MPA Database.Rdata',sep = ''))


regional_temptrend_plot <- processed_dat %>% group_by(region,year) %>% summarize(mt = mean(mean_temp)) %>%
  ggplot(aes(year,mt,fill = region)) +
  geom_point(shape = 21) +
  geom_smooth(aes(color = region),method = 'lm', se = F) +
  facet_grid(~region) +
  scale_fill_discrete(guide = F)+
  scale_color_discrete(guide = F) +
  ylab(expression(paste('Temperature (',degree ~ C,')', sep = ''))) +
  xlab('Year') #+
#   theme_fivethirtyeight()


temptrend_plot <- processed_dat %>% group_by(site,year) %>% summarize(mt = mean(mean_temp), region = unique(region)) %>%
  ggplot(aes(year,mt)) +
  geom_point(shape = 21, fill = 'lightseagreen') +
  stat_smooth() +
  #   geom_smooth(method = 'lm', se = F) +
  scale_color_discrete(guide = F) +
  #   theme_fivethirtyeight() +
  ylab(expression(paste('Temperature (',degree ~ C,')', sep = ''))) +
  xlab('Year')


# Format data for regression ----------------------------------------------

if (agg_level == 'siteside')
{

  # Generate temperature lags ----
  #

  temp_lags <- processed_dat %>%
    group_by(site_side,year) %>%
    summarise(mean_temp = unique(mean_temp)) %>%
    ungroup() %>%
    group_by(site_side) %>%
    mutate(mean_site_temp = mean(mean_temp, na.rm = T))

  temp_lags$mean_temp[is.na(temp_lags$mean_temp)] <- temp_lags$mean_site_temp[is.na(temp_lags$mean_temp)]

  temp_lags <- temp_lags %>%
    mutate(
           mean_temp_lag0 = (lag(mean_temp,0)),
           mean_temp_lag1 = (lag(mean_temp,1)),
           mean_temp_lag2 = lag(mean_temp,2),
           mean_temp_lag3 = lag(mean_temp,3),
           mean_temp_lag4 = lag(mean_temp,4))

  for (i in 1:dim(temp_lags)[1]){

    missing_lags <- which(is.na(temp_lags[i,]))

    temp_lags[i,missing_lags] <- temp_lags$mean_site_temp[i]

  }

  ggplot(temp_lags,aes(year,mean_temp, color = site_side)) +
    geom_line() +
    scale_color_discrete(guide = F)


  reg_data <- processed_dat %>% #group data to the species, site_side, year level
    subset(is.na(targeted) == F & targeted != '' & mpaareanm2 >0) %>%
    group_by(year,site_side,classcode) %>%
    summarise(mpa_period = unique(region_has_mpas),
              sites_checked = length(biomass),
              sites_seenat = sum(biomass >0),
              mean_density = mean(biomass, na.rm = T),
              site.type = unique(mpa.status),
              site = unique(site),
              site_resolution = unique(site_side),
              site_will_be_mlpa = unique(will_be_mlpa),
              year_mlpa = unique(year_mlpa_mpa),
              year_region_has_mlpas = unique(year_region_has_mlpa_mpas),
              years_mpa = max(0,(year - (year_mlpa_mpa-1)) * as.numeric(year_mlpa_mpa>0)), #years a site_side has been an MLPA MPA
              years_mlpa_mpas = max(0,(year - (year_region_has_mlpa_mpas-1)) * as.numeric(year_region_has_mlpa_mpas>0)), #years a region has had MLPA mpas
              region = unique(region),
              targeted = unique(targeted),
              trophic.group = unique(trophicgroup),
              mpa_area = mean(mpaareanm2, na.rm = T),
              linf = mean(c(vbgf.linf,vbgf.linf.f,vbgf.linf.m), na.rm = T),
              vbk = mean(c(vbgf.k,vbgf.k.f,vbgf.k.m), na.rm = T),
              size_mature = mean(size_mature_cm, na.rm = T),
              broadtrophic = unique(broadtrophic),
              mean_vis = mean(mean_vis, na.rm = T),
              mean_kelp = mean(pct_kelp_canopy, na.rm = T),
              species = unique(common_name)) %>%
    ungroup() %>%
    group_by(site_side, species) %>%
    mutate(anybio = sum(mean_density, na.rm = T)) %>%
    subset(anybio >0) %>%
    ungroup() %>%
    mutate(trunc_mean_density = pmax(quantile(unique(mean_density[mean_density >0]),.001, na.rm = T),mean_density)
           ,log_density = log(trunc_mean_density),
           fished = as.numeric(targeted == 'Targeted'),
           mpa_applied = as.numeric(mpa_period == 'TRUE'),
           fished_x_mpa = fished * mpa_applied,
           fished_x_yearsmlpa = fished * years_mlpa_mpas,
           fished_x_yearsmpa = fished * years_mpa,
           fishedeffect = as.factor(fished*year*mpa_applied),
           factor_year = as.factor(year)) %>%
    group_by(site_side) %>%
    mutate(eventual_mpa = as.numeric(max(years_mpa) >0)) %>%
    subset(is.na(mean_density) == F & is.na(log_density) == F) %>%
    ungroup() %>%
    mutate(fished_x_eventualmpa = fished * site_will_be_mlpa,
           fished_x_yearsmlpa_x_eventualmpa = fished * site_will_be_mlpa * years_mlpa_mpas) %>%
    left_join(temp_lags, by = c('site_side','year'))

}
if (agg_level == 'site'){

  temp_lags <- processed_dat %>%
    group_by(site,year) %>%
    summarise(mean_temp = unique(mean_temp)) %>%
    ungroup() %>%
    group_by(site) %>%
    mutate(mean_site_temp = mean(mean_temp, na.rm = T))

  temp_lags$mean_temp[is.na(temp_lags$mean_temp)] <- temp_lags$mean_site_temp[is.na(temp_lags$mean_temp)]

  temp_lags <- temp_lags %>%
    mutate(
      mean_temp_lag0 = lag(mean_temp,0),
      mean_temp_lag1 = lag(mean_temp,1),
           mean_temp_lag2 = lag(mean_temp,2),
           mean_temp_lag3 = lag(mean_temp,3),
           mean_temp_lag4 = lag(mean_temp,4))

  for (i in 1:dim(temp_lags)[1]){

    missing_lags <- which(is.na(temp_lags[i,]))

    temp_lags[i,missing_lags] <- temp_lags$mean_site_temp[i]

  }

  reg_data <- processed_dat %>% #group data to the species, site_side, year level
    subset(is.na(targeted) == F & targeted != '' & mpaareanm2 >0) %>%
    group_by(year,site,classcode) %>%
    summarise(mpa_period = unique(region_has_mpas),
              sites_checked = length(biomass),
              sites_seenat = sum(biomass >0),
              mean_density = mean(biomass, na.rm = T),
              site.type = unique(mpa.status),
              site_resolution = unique(site),
              #               site = unique(site),
              site_will_be_mlpa = unique(will_be_mlpa),
              year_mlpa = unique(year_mlpa_mpa),
              year_region_has_mlpas = unique(year_region_has_mlpa_mpas),
              years_mpa = max(0,(year - (year_mlpa_mpa-1)) * as.numeric(year_mlpa_mpa>0)), #years a site_side has been an MLPA MPA
              years_mlpa_mpas = max(0,(year - (year_region_has_mlpa_mpas-1)) * as.numeric(year_region_has_mlpa_mpas>0)), #years a region has had MLPA mpas
              region = unique(region),
              targeted = unique(targeted),
              trophic.group = unique(trophicgroup),
              mpa_area = mean(mpaareanm2, na.rm = T),
              linf = mean(c(vbgf.linf,vbgf.linf.f,vbgf.linf.m), na.rm = T),
              vbk = mean(c(vbgf.k,vbgf.k.f,vbgf.k.m), na.rm = T),
              size_mature = mean(size_mature_cm, na.rm = T),
              broadtrophic = unique(broadtrophic),
              mean_vis = mean(mean_vis, na.rm = T),
              mean_kelp = mean(pct_kelp_canopy, na.rm = T),
              species = unique(common_name)) %>%
    ungroup() %>%
    group_by(site, species) %>%
    mutate(anybio = sum(mean_density, na.rm = T)) %>%
    subset(anybio >0) %>%
    ungroup() %>%
    mutate(trunc_mean_density = pmax(quantile(unique(mean_density[mean_density >0]),.001, na.rm = T),mean_density)
           ,log_density = log(trunc_mean_density),
           fished = as.numeric(targeted == 'Targeted'),
           mpa_applied = as.numeric(mpa_period == 'TRUE'),
           fished_x_mpa = fished * mpa_applied,
           fished_x_yearsmlpa = fished * years_mlpa_mpas,
           fished_x_yearsmpa = fished * years_mpa,
           fishedeffect = as.factor(fished*year*mpa_applied),
           factor_year = as.factor(year)) %>%
    group_by(site) %>%
    mutate(eventual_mpa = as.numeric(max(years_mpa) >0)) %>%
    subset(is.na(mean_density) == F & is.na(log_density) == F) %>%
    ungroup() %>%
    mutate(fished_x_eventualmpa = fished * site_will_be_mlpa,
           fished_x_yearsmlpa_x_eventualmpa = fished * site_will_be_mlpa * years_mlpa_mpas) %>%
    left_join(temp_lags, by = c('site','year'))

}
varnames <- colnames(reg_data)

reg_data$MLPA <- ordered((reg_data$mpa_period), levels = c('BEFORE','AFTER'))

reg_data$na_temp <- (reg_data$mean_temp)

reg_data$na_temp[is.na(reg_data$mean_temp)] <- mean(reg_data$mean_temp, na.rm = T)

reg_data$na_vis <- (reg_data$mean_vis)

reg_data$na_vis[is.na(reg_data$mean_vis)] <- mean(reg_data$mean_vis, na.rm = T)

if (MPA_Only == TRUE){
  reg_data <- filter(reg_data, eventual_mpa == T)
}

# Plot some stuff ---------------------------------------------------------



density.by.time.in.mpa_plot <-
  ggplot(subset(reg_data), aes(factor(years_mpa),mean_density)) +
  geom_boxplot(aes(fill = factor(targeted))) +
  facet_grid(eventual_mpa~.) +
  scale_y_log10()



# Density vs Year by Targeted or Untargeted Inside and Outside of Eventual MPA ----

density.by.year_plot <-
  ggplot(subset(reg_data,mean_density >0), aes(factor(year),mean_density)) +
  geom_boxplot(aes(fill = factor(targeted))) +
  facet_grid(eventual_mpa~.) +
  scale_y_log10()


# Covatiate Plots ----

temp_gradient_plot <- reg_data %>%
  group_by(site_resolution) %>%
  mutate(mean_mean_temp = mean(mean_temp), na.rm = T) %>%
  subset(is.na(mean_mean_temp) == F) %>%
  ungroup() %>%
  mutate(factor_site_side = factor(site_resolution)) %>%
  ggplot(aes(reorder(factor_site_side,mean_mean_temp),mean_temp)) +
  geom_boxplot() +
  xlab('Site Side') +
  ylab("Mean Temperature ('C')") +
  coord_flip()

vis_gradient_plot <- reg_data %>%
  group_by(site_resolution) %>%
  mutate(mean_mean_vis = mean(mean_vis), na.rm = T) %>%
  subset(is.na(mean_mean_vis) == F) %>%
  ungroup() %>%
  mutate(factor_site_side = factor(site_resolution)) %>%
  ggplot(aes(reorder(factor_site_side,mean_mean_vis),mean_vis)) +
  geom_boxplot() +
  xlab('Site Side') +
  ylab("Mean Visibility") +
  coord_flip()

mean_density_plot <- reg_data %>%
  group_by(site_resolution) %>%
  mutate(mean_mean_density = mean(mean_density, na.rm = T)) %>%
  subset(is.na(mean_mean_density) == F) %>%
  ungroup() %>%
  mutate(factor_site_side = factor(site_resolution)) %>%
  ggplot(aes(reorder(factor_site_side,mean_mean_density),mean_density)) +
  geom_boxplot() +
  xlab('Site Side') +
  ylab("Mean Biomass Density") +
  scale_y_log10() +
  coord_flip()

perc_targeted_plot <- reg_data %>%
  group_by(site_resolution) %>%
  summarize(perc_targeted = mean(targeted == 'Targeted'), MPA = unique(eventual_mpa)) %>%
  ungroup() %>%
  mutate(factor_site_side = as.factor(site_resolution)) %>%
  ggplot(aes(reorder(factor_site_side,perc_targeted), perc_targeted,fill= factor(MPA) )) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(labels = percent) +
  xlab('Site Side') +
  ylab("% of Species Targeted") +
  coord_flip() +
  theme(axis.text.y = element_text(size = 4), axis.title.y = element_blank())

community_structure_plot <- reg_data %>%
  ggplot(aes(site_resolution,fill= trophic.group)) +
  geom_bar() +
  xlab('Site Side') +
  ylab("Community Structure") +
  coord_flip()


# Examine Correlations ----------------------------------------------------


correlations <- reg_data %>%
  ungroup() %>%
  dplyr::select(mean_density,site_resolution,year,mean_temp,mean_vis,mpa_area) %>%
  gather('variable','value',year:mpa_area) %>%
  ggplot(aes(value,mean_density, fill = variable)) +
  geom_point(shape = 21, alpha = 0.75) +
  scale_y_log10() +
  facet_wrap(~variable, scales = 'free')

write.csv(file = paste(runpath,'processed reg_data.csv', sep = ''),
          subset(reg_data, is.na(mean_density) == F))


# Run Regression ----------------------------------------------------------

# tobit_reg <- tobit(log_density ~
#                      factor(year) +
#                      factor(region) +
#                      years_mpa*fished +
#                      factor(trophic.group) +
#                      years_mpa +
#                      linf +
#                      vbk +
#                      mpa_applied +
#                      fished +
#                      mpa_applied*fished +
#                      na_temp +
#                      na_vis +
#                      mean_temp_lag1 +
#                      mean_temp_lag2 +
#                      mean_temp_lag3 +
#                      mean_temp_lag4,
#                    data = subset(reg_data, is.na(log_density) == F),
#                    left = min(reg_data$log_density, na.rm = T))
#
# vcov_plot <- as.data.frame(vcov(tobit_reg)) %>%
#   mutate(var2 = rownames(.)) %>%
#   gather('var1','vcov',1:(dim(.)[2]-1)) %>%
#   ggplot(aes(var1,var2,fill = vcov)) +
#   geom_raster() +
#   scale_fill_gradient2(low = 'red', high = 'green', mid = 'white',midpoint = 0)
#
# summary(tobit_reg)
#
# reg_results <- augment(tobit_reg)
#
# reg_results$observed <- reg_results$.fitted + reg_results$.resid
#
#
# write.csv(file = paste(runpath,'tobit data.csv', sep = ''),
#           reg_results)
#
# pdf(file = paste(runpath,'tobit diagnostics.pdf', sep = ''))
# par(mfrow = c(2,3))
# plot((reg_results$factor.year.),reg_results$.resid)
# boxplot(reg_results$.resid ~ reg_results$factor.region.)
# plot(reg_results$observed,reg_results$.fitted)
# abline(0,1, col = 'blue')
# abline(lm(reg_results$.fitted ~ reg_results$observed), col="red") # regression line (y~x)
# hist(reg_results$.resid)
# abline(v = mean(reg_results$.resid))
# qqnorm(reg_results$.resid)
# qqline(reg_results$.resid)
# dev.off()
#
# summary(tobit_reg)
#
#
# eval_resid <- reg_results %>%
#   ggplot(aes(observed,.fitted, fill = factor.region.)) +
#   geom_point(shape = 21, alpha = .4) +
#   geom_abline(intercept = 0, slope = 1, color = 'red')
#


# Prep Bayesian Regression ------------------------------------------------------------




# pos_vars <- c('fished','years_mpa','fished_x_yearsmpa','factor_year',
#               'region','trophic.group',
#               'linf','vbk','na_temp','na_vis', 'mean_temp_lag1', 'mean_temp_lag2','mean_temp_lag3',
#               'mean_temp_lag4')
#
# delta_vars <- c('fished','years_mpa','fished_x_yearsmpa','factor_year',
#               'trophic.group',
#               'linf','na_temp','na_vis')

any_fish <- reg_data$log_density > min(reg_data$log_density, na.rm = T)

logit_model = glm(any_fish ~  fished + mpa_applied + fished_x_mpa + trophic.group +
                    + linf + na_temp + na_vis ,family = 'binomial', data = reg_data)

logit_reg <- augment(logit_model)

logit_reg$prob <- exp(logit_reg$.fitted)/(1+exp(logit_reg$.fitted))

logit_reg_plot <- logit_reg %>%
  mutate(bin = cut(prob,100)) %>%
  group_by(bin) %>%
  summarize(seen = mean(any_fish)) %>%
  ggplot(aes(bin,seen)) +
  geom_bar(stat = 'identity', position = 'dodge')

plots <- ls()

plots <- plots[grepl('_plot',plots, fixed = T)]

save(file = paste(runpath,'MLPA Plots.Rdata', sep = ''), list = plots)


# Run Bayesian Regression -------------------------------------------------

if (run_mcmc == T){


  dep_var <- 'log_density'


  pos_vars <- c('fished','mpa_applied','fishedeffect',
                'region','linf','vbk','broadtrophic','mean_temp_lag0',
                'mean_temp_lag1','mean_temp_lag2','mean_temp_lag3','mean_temp_lag4') #, 'mean_kelp')

  delta_vars <- c('fished','mpa_applied','fishedeffect',
                  'broadtrophic',
                  'linf','mean_temp_lag0','na_vis') #,'mean_kelp')

  # Try frequentist deltaGLM ----

  density_fmla <- formula(paste('mean_density ~',paste(pos_vars, collapse = '+'), sep = ''))

  logit_fmla <- formula(paste('~ 1+',paste(pos_vars, collapse = '+'), sep = ''))

  has_all_needed <- (rowSums(is.na(reg_data[,pos_vars])) & rowSums(is.na(reg_data[,delta_vars]))) == 0

  reg_data <- reg_data[has_all_needed,]

  save(file = paste(runpath,'reg_data.Rdata', sep = ''), reg_data)

  delta_glm <- deltaLN(ln.form = density_fmla,binary.form = logit_fmla, data = reg_data)
#   reg_data$binom_predict <- exp(predict(delta_glm$binMod))/(1+exp(predict(delta_glm$binMod)))
#
#   reg_data$binom <- as.numeric(reg_data$mean_density >0)

#   quartz()
#   plot(reg_data$binom_predict, reg_data$binom)

#   reg_data$density_hat <- delta_glm$fitted
#
#   arg <- reg_data %>%
#     filter(mean_density>0) %>%
#     group_by(fished,year) %>%
#     summarize(mean_observed = mean(log_density), mean_predicted = mean(log(density_hat))) %>%
#     gather('data_type','density', mean_observed:mean_predicted)
#
#     ggplot() +
#     geom_point(data = filter(arg,data_type == 'mean_observed'),aes(year,density,fill = factor(fished)), shape = 21) +
#     geom_line(data = filter(arg,data_type == 'mean_predicted'), aes(year,density,color = factor(fished)))




#   quartz()
#
#   plot(log(reg_data$mean_density[complete_reg_data$mean_density>0]),log(delta_glm$fitted[reg_data$mean_density>0]))
#

  # Run MCMC ----------------------------------------------------------------


#   reg_fmla <- formula(paste('log_density ~',paste(pos_vars, collapse = '+'), sep = ''))
#
#   summary(lm(reg_fmla, data = check))

  browser()

 log_density_fmla <- formula(paste('log_density ~',paste(pos_vars, collapse = '+'), sep = ''))


  stan_reg <-  rstanarm::stan_glm(formula = log_density_fmla, family = gaussian, data = reg_data %>% filter(mean_density > 0))


  a <- rstanarm::posterior_linpred(stan_reg)

  rstanarm::posterior_vs_prior(stan_reg)

  rstanarm::posterior_interval(stan_reg)


  stan_reg %>%
    tidy() %>%
    filter(str_detect(term,'fishedeffect')) %>%
    mutate(year = str_replace(term,'fishedeffect','')) %>%
    ggplot() +
    geom_pointrange(aes(x = year, y = estimate, ymin = estimate - 1.96*std.error, ymax = estimate + 1.96 * std.error))

   launch_shinystan(stan_reg)

  bayes_reg <- run_delta_demon(dat = reg_data, method = model_form,dep_var = dep_var,
                               pos_vars = pos_vars, delta_vars = delta_vars,runpath = runpath,scale_numerics = scale_numerics,
                               iterations = its,status = .01, acceptance_rate = 0.43, thin = its/1e4)


  reg_results <- list(acceptance_rate = bayes_reg$demon_fit$Acceptance.Rate, post =
                        bayes_reg$demon_fit$Posterior1, Data = bayes_reg$Data)

  save(file = paste(runpath,'MCMC results.Rdata', sep = ""), reg_results)
}

processed_demon <- process_demon(runfolder = runfolder, post_sample_size = 1000, burn = 0.75)

grid.draw(processed_demon$plot_list$density_summary_plot)

ggsave(file = paste(runpath,'Posterior summary plot.pdf', sep = ""), processed_demon$plot_list$density_summary_plot)

# outlist <- list(thinned_post = processed_demon$thinned_post,data_and_predictions = processed_demon$predictions, diagnostic_plots = processed_demon$plot_list)

pres_figures <- mpa_presentation_figs(run = runfolder,
                                      out = 'PhD Symposium Figs',
                                      fig_width = 7, fig_height = 4, font_size = 14)

save(file = paste(runpath,'pres_quality_figures.Rdata', sep = ""),  pres_figures)


save(file = paste(runpath,'Processed MCMC results.Rdata', sep = ""),  processed_demon)



