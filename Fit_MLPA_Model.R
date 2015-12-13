
# Set up ------------------------------------------------------------------
# This script runs the development version of the MLPA effects model

rm(list = ls())
set.seed(666)
library(knitr)
library(gridExtra)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(broom)
library(coda)
library(ggmcmc)
library(LaplacesDemon)
library(foreach)
library(scales)
library(stargazer)
library(DT)
library(ggmap)
library(texreg)
library(AER)
library(msm)
library(mvtnorm)
devtools::load_all('MLPAFuns')


# Run Options -------------------------------------------------------------


runfolder <- '7.0'

scale_numerics <- T

its <- 1e4

run_mcmc <- T

runpath <- paste('MLPA Effects Results/',runfolder,'/', sep = '')

if (dir.exists(runpath) == F)
{
  dir.create(runpath, recursive = T)
}


# Load Data ---------------------------------------------------------------


rawdat <- read.csv('MLPA Data/UCSB_FISH raw thru 2013.csv', stringsAsFactors = F)

life.history <- read.csv('MLPA Data/VRG Fish Life History in MPA_04_08_11_12 11-Mar-2014.csv', stringsAsFactors = F) %>%
  rename(classcode = pisco_classcode)

species.dat <- read.csv('MLPA Data/master_spp_table.csv', stringsAsFactors = F)

site.dat <- read.csv('MLPA Data/Final_Site_Table_UCSB.csv', stringsAsFactors = F) %>%
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


processed_site_dat <- read.csv('MLPA Data/ci_reserve_data_final3 txt.csv', stringsAsFactors = F) %>%
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

save(mlpa_dat,processed_dat, file = paste(runpath,'Combined MLPA Database.Rdata',sep = ''))


temptrend_plot <- processed_dat %>% group_by(region,year) %>% summarize(mt = mean(mean_temp)) %>%
  ggplot(aes(year,mt,color = region)) +
  geom_line()

# Generate temperature lags ----

temp_lags <- processed_dat %>%
  group_by(site_side,year) %>%
  summarise(mean_temp = unique(mean_temp)) %>%
  ungroup() %>%
  group_by(site_side) %>%
  mutate(mean_site_temp = mean(mean_temp, na.rm = T))

temp_lags$mean_temp[is.na(temp_lags$mean_temp)] <- temp_lags$mean_site_temp[is.na(temp_lags$mean_temp)]

temp_lags <- temp_lags %>%
  mutate(mean_temp_lag1 = (lag(mean_temp,1)),
         mean_temp_lag2 = lag(mean_temp,2),
         mean_temp_lag3 = lag(mean_temp,3),
         mean_temp_lag4 = lag(mean_temp,4))

for (i in 1:dim(temp_lags)[1]){

  missing_lags <- which(is.na(temp_lags[i,]))

  temp_lags[i,missing_lags] <- temp_lags$mean_site_temp[i]

}


# Format data for regression ----------------------------------------------


species_siteside_year <- processed_dat %>% #group data to the species, site_side, year level
  subset(is.na(targeted) == F & targeted != '' & mpaareanm2 >0) %>%
  group_by(year,site_side,classcode) %>%
  summarise(mpa_period = unique(region_has_mpas),
            mean_density = mean(biomass, na.rm = T),
            site.type = unique(mpa.status),
            site = unique(site),
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
         factor_year = as.factor(year)) %>%
  group_by(site_side) %>%
  mutate(eventual_mpa = max(years_mpa) >0) %>%
  subset(is.na(mean_density) == F & is.na(log_density) == F) %>%
  ungroup() %>%
  mutate(fished_x_eventualmpa = fished * site_will_be_mlpa,
         fished_x_yearsmlpa_x_eventualmpa = fished * site_will_be_mlpa * years_mlpa_mpas) %>%
  left_join(temp_lags, by = c('site_side','year'))


varnames <- colnames(species_siteside_year)

species_siteside_year$MLPA <- ordered((species_siteside_year$mpa_period), levels = c('BEFORE','AFTER'))

species_siteside_year$na_temp <- (species_siteside_year$mean_temp)

species_siteside_year$na_temp[is.na(species_siteside_year$mean_temp)] <- mean(species_siteside_year$mean_temp, na.rm = T)

species_siteside_year$na_vis <- (species_siteside_year$mean_vis)

species_siteside_year$na_vis[is.na(species_siteside_year$mean_vis)] <- mean(species_siteside_year$mean_vis, na.rm = T)

save(file = paste(runpath,'species_siteside_year.Rdata', sep = ''), species_siteside_year)

# Plot some stuff ---------------------------------------------------------



density.by.time.in.mpa_plot <-
  ggplot(subset(species_siteside_year), aes(factor(years_mpa),mean_density)) +
  geom_boxplot(aes(fill = factor(targeted))) +
  facet_grid(eventual_mpa~.) +
  scale_y_log10()



# Density vs Year by Targeted or Untargeted Inside and Outside of Eventual MPA ----

density.by.year_plot <-
  ggplot(subset(species_siteside_year,mean_density >0), aes(factor(year),mean_density)) +
  geom_boxplot(aes(fill = factor(targeted))) +
  facet_grid(eventual_mpa~.) +
  scale_y_log10()


# Covatiate Plots ----

temp_gradient_plot <- species_siteside_year %>%
  group_by(site_side) %>%
  mutate(mean_mean_temp = mean(mean_temp), na.rm = T) %>%
  subset(is.na(mean_mean_temp) == F) %>%
  ungroup() %>%
  mutate(factor_site_side = factor(site_side)) %>%
  ggplot(aes(reorder(factor_site_side,mean_mean_temp),mean_temp)) +
  geom_boxplot() +
  xlab('Site Side') +
  ylab("Mean Temperature ('C')") +
  coord_flip()

vis_gradient_plot <- species_siteside_year %>%
  group_by(site_side) %>%
  mutate(mean_mean_vis = mean(mean_vis), na.rm = T) %>%
  subset(is.na(mean_mean_vis) == F) %>%
  ungroup() %>%
  mutate(factor_site_side = factor(site_side)) %>%
  ggplot(aes(reorder(factor_site_side,mean_mean_vis),mean_vis)) +
  geom_boxplot() +
  xlab('Site Side') +
  ylab("Mean Visibility") +
  coord_flip()

mean_density_plot <- species_siteside_year %>%
  group_by(site_side) %>%
  mutate(mean_mean_density = mean(mean_density, na.rm = T)) %>%
  subset(is.na(mean_mean_density) == F) %>%
  ungroup() %>%
  mutate(factor_site_side = factor(site_side)) %>%
  ggplot(aes(reorder(factor_site_side,mean_mean_density),mean_density)) +
  geom_boxplot() +
  xlab('Site Side') +
  ylab("Mean Biomass Density") +
  scale_y_log10() +
  coord_flip()

perc_targeted_plot <- species_siteside_year %>%
  group_by(site_side) %>%
  summarize(perc_targeted = mean(targeted == 'Targeted')) %>%
  ungroup() %>%
  mutate(factor_site_side = as.factor(site_side)) %>%
  ggplot(aes(reorder(factor_site_side,perc_targeted), perc_targeted,fill= perc_targeted)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(labels = percent) +
  xlab('Site Side') +
  ylab("% of Species Targeted") +
  coord_flip()

community_structure_plot <- species_siteside_year %>%
  ggplot(aes(site_side,fill= trophic.group)) +
  geom_bar() +
  xlab('Site Side') +
  ylab("Community Structure") +
  coord_flip()


# Examine Correlations ----------------------------------------------------


correlations <- species_siteside_year %>%
  ungroup() %>%
  dplyr::select(mean_density,site_side,year,mean_temp,mean_vis,mpa_area) %>%
  gather('variable','value',year:mpa_area) %>%
  ggplot(aes(value,mean_density, fill = variable)) +
  geom_point(shape = 21, alpha = 0.75) +
  scale_y_log10() +
  facet_wrap(~variable, scales = 'free')

write.csv(file = paste(runpath,'species_siteside_year.csv', sep = ''),
          subset(species_siteside_year, is.na(mean_density) == F))


# Run Regression ----------------------------------------------------------

tobit_reg <- tobit(log_density ~
                     factor(year) +
                     factor(region) +
                     years_mpa*fished +
                     factor(trophic.group) +
                     years_mpa +
                     linf +
                     vbk +
                     mpa_applied +
                     fished +
                     mpa_applied*fished +
                     na_temp +
                     na_vis +
                     mean_temp_lag1 +
                     mean_temp_lag2 +
                     mean_temp_lag3 +
                     mean_temp_lag4,
                   data = subset(species_siteside_year, is.na(log_density) == F),
                   left = min(species_siteside_year$log_density, na.rm = T))

vcov_plot <- as.data.frame(vcov(tobit_reg)) %>%
  mutate(var2 = rownames(.)) %>%
  gather('var1','vcov',1:(dim(.)[2]-1)) %>%
  ggplot(aes(var1,var2,fill = vcov)) +
  geom_raster() +
  scale_fill_gradient2(low = 'red', high = 'green', mid = 'white',midpoint = 0)

summary(tobit_reg)

reg_results <- augment(tobit_reg)

reg_results$observed <- reg_results$.fitted + reg_results$.resid


write.csv(file = paste(runpath,'tobit data.csv', sep = ''),
          reg_results)

pdf(file = paste(runpath,'tobit diagnostics.pdf', sep = ''))
par(mfrow = c(2,3))
plot((reg_results$factor.year.),reg_results$.resid)
boxplot(reg_results$.resid ~ reg_results$factor.region.)
plot(reg_results$observed,reg_results$.fitted)
abline(0,1, col = 'blue')
abline(lm(reg_results$.fitted ~ reg_results$observed), col="red") # regression line (y~x)
hist(reg_results$.resid)
abline(v = mean(reg_results$.resid))
qqnorm(reg_results$.resid)
qqline(reg_results$.resid)
dev.off()

summary(tobit_reg)


eval_resid <- reg_results %>%
  ggplot(aes(observed,.fitted, fill = factor.region.)) +
  geom_point(shape = 21, aes = .4) +
  geom_abline(intercept = 0, slope = 1, color = 'red')



# Prep Bayesian Regression ------------------------------------------------------------




# pos_vars <- c('fished','years_mpa','fished_x_yearsmpa','factor_year',
#               'region','trophic.group',
#               'linf','vbk','na_temp','na_vis', 'mean_temp_lag1', 'mean_temp_lag2','mean_temp_lag3',
#               'mean_temp_lag4')
#
# delta_vars <- c('fished','years_mpa','fished_x_yearsmpa','factor_year',
#               'trophic.group',
#               'linf','na_temp','na_vis')

any_fish <- species_siteside_year$log_density > min(species_siteside_year$log_density, na.rm = T)

logit_model = glm(any_fish ~  fished + mpa_applied + fished_x_mpa + trophic.group +
                    + linf + na_temp + na_vis ,family = 'binomial', data = species_siteside_year)

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

  devtools::load_all('MLPAFuns')

  check <- subset(species_siteside_year, mean_density >0)

  dep_var <- 'log_density'

#   pos_vars <- c('fished','years_mpa','fished_x_yearsmpa','factor_year',
#                 'region','linf','vbk','trophic.group','na_temp','na_vis', 'mean_temp_lag1', 'mean_temp_lag2','mean_temp_lag3',
#                 'mean_temp_lag4')

  pos_vars <- c('fished','years_mlpa_mpas','fished_x_yearsmlpa','factor_year',
                'region','linf','vbk','trophic.group','na_temp','na_vis', 'mean_temp_lag1', 'mean_temp_lag2','mean_temp_lag3',
                'mean_temp_lag4')

#   species_siteside_year %>%
#     group_by(region,years_mlpa_mpas) %>%
#     summarise() %>%
#     ggplot(aes)


  delta_vars <- c('fished','years_mlpa_mpas','fished_x_yearsmlpa','factor_year',
                  'trophic.group',
                  'linf','na_temp','na_vis')

  reg_fmla <- formula(paste('log_density ~',paste(pos_vars, collapse = '+'), sep = ''))

  summary(lm(reg_fmla, data = check))

  bayes_reg <- run_delta_demon(dat = species_siteside_year, method = 'Summon Demon',dep_var = dep_var,
                               pos_vars = pos_vars, delta_vars = delta_vars,runpath = runpath,scale_numerics = scale_numerics,
                               iterations = its,status = .05, acceptance_rate = 0.3, thin = its/1e4)


  reg_results <- list(acceptance_rate = bayes_reg$demon_fit$Acceptance.Rate, post =
                        bayes_reg$demon_fit$Posterior1, Data = bayes_reg$Data)

  save(file = paste(runpath,'MCMC results.Rdata', sep = ""), list = 'reg_results')
}

processed_demon <- process_demon(runfolder = runfolder)

outlist <- list(thinned_post = processed_demon$thinned_post,data_and_predictions = processed_demon$predictions, diagnostic_plots = processed_demon$plot_list)

save(file = paste(runpath,'Processed MCMC results.Rdata', sep = ""), list = 'outlist')



