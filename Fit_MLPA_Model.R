
# Set up ------------------------------------------------------------------
# This script runs the development version of the MLPA effects model

rm(list = ls())
library(knitr)
knitr::opts_chunk$set(fig.path='Figs/', echo=FALSE, warning=FALSE, message=FALSE)
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
# library(VGAM)
library(AER)
library(msm)
library(mvtnorm)
devtools::load_all('MLPAFuns')


devtools::load_all('~/Custom R Packages/ZurrCode')
devtools::load_all('~/Custom R Packages/RobustRegression')



# Run Options -------------------------------------------------------------


runfolder <- 'development'

runpath <- paste('MLPA Effects Results/',runfolder,'/', sep = '')

if (dir.exists(runpath) == F)
{
  dir.create(runpath, recursive = T)
}


# Load Data ---------------------------------------------------------------


rawdat <- read.csv('UCSB_FISH raw thru 2013.csv', stringsAsFactors = F)

life.history <- read.csv('VRG Fish Life History in MPA_04_08_11_12 11-Mar-2014.csv', stringsAsFactors = F) %>%
  rename(classcode = pisco_classcode)


# lh <- life.history %>%
#   group_by(classcode) %>%
#   summarise(L_inf = mean(linf),VbK = mean(vbk),size_mat = mean(size_mature)) %>%
#   ungroup() %>%
#   mutate(hasall = is.na(L_inf) == F & is.na(VbK) == F & is.na(size_mat) == F)

species.dat <- read.csv('master_spp_table.csv', stringsAsFactors = F)

site.dat <- read.csv('Final_Site_Table_UCSB.csv', stringsAsFactors = F) %>%
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


processed_site_dat <- read.csv('ci_reserve_data_final3 txt.csv', stringsAsFactors = F) %>%
  gather('concat.name','value', grep('_',colnames(.)),convert = T) %>%
  mutate(data.type = gsub('\\_.*', '', concat.name),
         classcode = gsub('.*\\_','',concat.name)) %>%
  mutate(value = as.numeric(value)) %>%
  spread(data.type,value) %>%
  rename(site_side = site.side) %>%
  left_join(conditions_dat,by = c('year','site_side')) %>%
  left_join(life.history, by = 'classcode') %>%
  rename(description2 = Description)

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
  subset(year != 1999 & ever.seen == T & grepl('YOY',classcode) == F) %>% #remove 1999, species never seen at a site, and young of the year. Why remove YOY again?
  ungroup() %>%
  group_by(site_side) %>%
  mutate(num.years.surveyed = length(unique(year)),
         will_be_mpa = year.mpa>0) %>%
  subset(num.years.surveyed > 2) %>% #Only keep sites that have been surveyed two or more years
  ungroup() %>%
  mutate(year_mpa2 = capfun(year.mpa)) %>%
  group_by(region,year) %>%
  mutate(region_has_mpas = any(year_mpa2 <= year)) %>%
  ungroup() %>%
  group_by(common_name) %>%
  mutate(num_samples = sum(biomass>0, na.rm = T)) %>%
  ungroup() %>%
  subset(num_samples > quantile(num_samples,.1))

save(mlpa_dat,processed_dat, file = paste(runpath,'Combined MLPA Database.Rdata',sep = ''))


# Format data for regression ----------------------------------------------

species_siteside_year <- processed_dat %>% #group data to the species, site_side, year level
  subset(is.na(targeted) == F & targeted != '' & mpaareanm2 >0) %>%
  group_by(year,site_side,classcode) %>%
  summarise(mpa.period = unique(region_has_mpas),
            mean_density = mean(biomass, na.rm = T),
            site.type = unique(mpa.status),
            years_mpa = max(0,(year - (year.mpa-1)) * as.numeric(year.mpa>0)),
            region = unique(region),
            targeted = unique(targeted),
            trophic.group = unique(trophicgroup),
            mpa_area = mean(mpaareanm2, na.rm = T),
            linf = mean(c(vbgf.linf,vbgf.linf.f,vbgf.linf.m), na.rm = T),
            vbk = mean(c(vbgf.k,vbgf.k.f,vbgf.k.m), na.rm = T),
            size_mature = mean(size_mature_cm, na.rm = T),
            broadtrophic = unique(broadtrophic),
            mean_temp = mean(mean_temp, na.rm = T),
            mean_vis = mean(mean_vis, na.rm = T),
            species = unique(common_name)) %>%
  #   ungroup() %>% min(mean_density[mean_density>0], na.rm = T)
  ungroup() %>%
  group_by(site_side, species) %>%
  mutate(anybio = sum(mean_density, na.rm = T)) %>%
  subset(anybio >0) %>%
  ungroup() %>%
  mutate(trunc_mean_density = pmax(quantile(unique(mean_density[mean_density >0]),.001, na.rm = T),mean_density)
         ,log_density = log(trunc_mean_density),
         fished = as.numeric(targeted == 'Targeted'),
         mpa_applied = as.numeric(mpa.period == 'TRUE'),
         fished_x_mpa = fished * mpa_applied,
         fished_x_yearsmpa = fished * years_mpa,
         factor_year = as.factor(year)) %>%
  group_by(site_side) %>%
  mutate(eventual_mpa = max(years_mpa) >0) %>%
  subset(is.na(mean_density) == F & is.na(log_density) == F) %>%
  ungroup()


varnames <- colnames(species_siteside_year)


# Plot some stuff ---------------------------------------------------------


species_siteside_year$MLPA <- ordered((species_siteside_year$mpa.period), levels = c('BEFORE','AFTER'))

density.by.time.in.mpa <-
  ggplot(subset(species_siteside_year), aes(factor(years_mpa),mean_density)) +
  geom_boxplot(aes(fill = factor(targeted))) +
  facet_grid(eventual_mpa~.) +
  scale_y_log10()

density.by.time.in.mpa


# Density vs Year by Targeted or Untargeted Inside and Outside of Eventual MPA ----

density.by.year <-
  ggplot(subset(species_siteside_year,mean_density >0), aes(factor(year),mean_density)) +
  geom_boxplot(aes(fill = factor(targeted))) +
  facet_grid(eventual_mpa~.) +
  scale_y_log10()

density.by.year

# Covatiate Plots ----

temp_gradient <- species_siteside_year %>%
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

vis_gradient <- species_siteside_year %>%
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

mean_density <- species_siteside_year %>%
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

perc_targeted <- species_siteside_year %>%
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

community_structure <- species_siteside_year %>%
  ggplot(aes(site_side,fill= trophic.group)) +
  geom_bar() +
  xlab('Site Side') +
  ylab("Community Structure") +
  coord_flip()

community_structure <- species_siteside_year %>%
  ggplot(aes(site_side,fill= trophic.group)) +
  geom_bar() +
  xlab('Site Side') +
  ylab("Community Structure") +
  coord_flip()


# Examine Correlations ----------------------------------------------------


correlations <- species_siteside_year %>%
  ungroup() %>%
  select(mean_density,site_side,year,mean_temp,mean_vis,mpa_area) %>%
  gather('variable','value',year:mpa_area) %>%
  ggplot(aes(value,mean_density, fill = variable)) +
  geom_point(shape = 21, alpha = 0.75) +
  scale_y_log10() +
  facet_wrap(~variable, scales = 'free')

write.csv(file = paste(runpath,'species_siteside_year.csv', sep = ''),
          subset(species_siteside_year, is.na(mean_density) == F))


# Run Regression ----------------------------------------------------------

# reg <- lm(log_density ~  factor(year) + mpa_applied + fished + mpa_applied*fished  , data = (species_siteside_year))

# tobit_reg <- vglm(log_density ~  factor(year) + factor(region) + mpa_applied + fished + mpa_applied*fished + fished*years_mpa + years_mpa +
#                     factor(broadtrophic) + size_mature,
#                   data = subset(species_siteside_year, is.na(log_density) == F),
#                   tobit(Lower = min(species_siteside_year$log_density, na.rm = T)),control = vglm.control(maxit = 100),trace = T)
#
# tobit_reg <- tobit(log_density ~  factor(year) + factor(region) + mpa_applied + fished + mpa_applied*fished + fished*years_mpa + years_mpa +
#                     factor(broadtrophic) + size_mature,
#                   data = subset(species_siteside_year, is.na(log_density) == F),
#                   left = min(species_siteside_year$log_density, na.rm = T))
#

species_siteside_year$na_temp <- (species_siteside_year$mean_temp)

species_siteside_year$na_temp[is.na(species_siteside_year$mean_temp)] <- mean(species_siteside_year$mean_temp, na.rm = T)

species_siteside_year$na_vis <- (species_siteside_year$mean_vis)

species_siteside_year$na_vis[is.na(species_siteside_year$mean_vis)] <- mean(species_siteside_year$mean_vis, na.rm = T)

pos_vars <- c('fished','mpa_applied','fished_x_mpa','factor_year',
              'region','broadtrophic','years_mpa','fished_x_yearsmpa',
              'linf','vbk')

species_siteside_year %>%
  group_by(species) %>%
  summarize(l = mean(linf),k = mean(vbk))

tobit_reg <- tobit(log_density ~
                     factor(year) +
                     factor(region) +
                     years_mpa*fished +
                     factor(broadtrophic) +
                     years_mpa +
                     linf +
                     vbk +
                     mpa_applied +
                     fished +
                     mpa_applied*fished,
                   data = subset(species_siteside_year, is.na(log_density) == F),
                   left = min(species_siteside_year$log_density, na.rm = T))
#
# tobit_reg <- tobit(log_density ~  (year) +
#                      years_mpa,data = subset(species_siteside_year, is.na(log_density) == F),
#                    left = min(species_siteside_year$log_density, na.rm = T))
#

# quartz()
# image(a)

vcov_plot <- as.data.frame(vcov(tobit_reg)) %>%
  mutate(var2 = rownames(.)) %>%
  gather('var1','vcov',1:(dim(.)[2]-1)) %>%
  ggplot(aes(var1,var2,fill = vcov)) +
  geom_raster() +
  scale_fill_gradient2(low = 'red', high = 'green', mid = 'white',midpoint = 0)

summary(tobit_reg)


reg_results <- augment(tobit_reg)



# predicted_log_density <- as.data.frame(predict(tobit_reg))
#
# observed_log_density <- as.data.frame(depvar(tobit_reg))
#
# observed_log_density$row <- rownames(observed_log_density)
#
# used <- species_siteside_year[as.numeric(observed_log_density$row),]
#
# used$predicted_log_density <- fitted(tobit_reg)

# used$predicted_log_density <- predicted_log_density$mu

# used$residuals <- used$predicted_log_density - used$log_density

reg_results$observed <- reg_results$.fitted + reg_results$.resid


write.csv(file = paste(runpath,'tobit data.csv', sep = ''),
          reg_results)

quartz()
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

summary(tobit_reg)


eval_resid <- reg_results %>%
  ggplot(aes(observed,.fitted, fill = factor.region.)) +
  geom_point(shape = 21, aes = .4) +
  geom_abline(intercept = 0, slope = 1, color = 'red')



# Develop Bayesian Regression ------------------------------------------------------------

dep_var <- 'log_density'


pos_vars <- c('fished','mpa_applied','fished_x_mpa','factor_year',
              'region','broadtrophic','years_mpa','fished_x_yearsmpa',
              'linf','vbk')

devtools::load_all('MLPAFuns')

a = which(species_siteside_year$log_density == min(species_siteside_year$log_density))

b <- sample(a,.99*length(a), replace = F)

# Develop delta method ----------------------------------------------------

any_fish <- species_siteside_year$log_density > min(species_siteside_year$log_density, na.rm = T)

a = glm(any_fish ~ linf + fished + fished_x_mpa + broadtrophic
        + mpa_applied + mean_vis + region + factor(year) ,family = 'binomial', data = species_siteside_year)

b <- augment(a)

b$prob <- exp(b$.fitted)/(1+exp(b$.fitted))

b %>%
  mutate(bin = cut(prob,100)) %>%
  group_by(bin) %>%
  summarize(seen = mean(any_fish)) %>%
  ggplot(aes(bin,seen)) +
  geom_bar(stat = 'identity', position = 'dodge')


bayesian_tobit <- run_mlpa_demon(dat = species_siteside_year, dep_var = dep_var,
                                 pos_vars = pos_vars,runpath = runpath,scale_numerics = T,
                                 iterations = 1e5,status = .025, acceptance_rate = 0.22, method = 'Summon Demon')

bayesian_tobit$Acceptance.Rate

# post <- thin_mcmc(bayesian_tobit$posteriors, thin_every = 10)
#
# ggmcmc(ggs(mcmc(post)), file = paste(runpath,'mcmc diagnostics.pdf', sep = ''))

post <- bayesian_tobit$Posterior1

post <- thin_mcmc(post[2.5e5:dim(post)[1],], thin_every = 4*bayesian_tobit$Rec.Thinning)

ggmcmc(ggs(mcmc(post)), file = paste(runpath,'mcmc diagnostics.pdf', sep = ''))

