## ---- include=FALSE------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, dev = "CairoPDF", dpi = 600,
                      cache = FALSE)
# knitr::opts_knit$set(root.dir = "/Users/danovando/PhD/zissou") # not entirely reprodcible hack to deal with knit child in dissertation.. not idealk



## ---- include=FALSE------------------------------------------------------
library(sf)
library(ggmap)
library(viridis)
library(hrbrthemes)
library(scales)
library(patchwork)
library(rpart)
library(extrafont)
library(caret)
library(ggsci)
library(rEDM)
library(tidyverse)
extrafont::loadfonts()
functions <- list.files(here::here("functions"))

walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions
save_plots <- TRUE

sim_years <- 50

burn_years <- 20

num_patches <- 25

  samps <- 20000
 
  run_name <- "v4.0"

run_dir <- here::here("results", run_name)

experiment_dir <- here::here("results",run_name, "experiments")

load(file = here::here("results",run_name, "rawish_zissou_data.Rdata"))

load(file = here::here("results",run_name, "model_runs.Rdata"))

load(file = here::here("results",run_name,"processed_grid.Rdata"))

load(file = here::here("results",run_name,"about_time.Rdata"))

channel_islands <- readRDS(here::here("data","channel_islands_map.RDS"))

ca_mpas <- sf::st_read(here::here("data","MPA_CA_Existing_160301")) %>% 
  rmapshaper::ms_simplify() %>% 
  sf::st_transform(crs = 4326)



zissou_theme <-
  theme_ipsum(
    base_size = 12,
    axis_title_size = 12,
    strip_text_size = 12
  )

theme_set(zissou_theme)


models_worked <- model_runs$tmb_fit %>% map("error") %>% map_lgl(is_null)

model_runs <- model_runs %>%
  filter(models_worked) %>%
  mutate(tmb_fit = map(tmb_fit,"result")) %>%
  mutate(processed_fits = map(tmb_fit, process_fits)) %>%
  mutate(did_plot = map(processed_fits, "did_plot"))

base_run <- model_runs %>% 
filter(var_names == "pisco_a", mpa_only == FALSE, center_scale == TRUE)

mpa_run <- model_runs %>% 
  filter(data_source == "pisco",
         var_names == "pisco_a",
         mpa_only == TRUE)

kfm_run <- model_runs %>% 
  filter(data_source == "kfm")
  
fitted_data <- base_run$data[[1]]

zissou_fit <- base_run$tmb_fit[[1]]

report <- base_run$tmb_fit[[1]]

site_data <- read_csv(here::here("data",'Final_Site_Table_UCSB.csv')) %>%
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


seen_non_nested_betas <- zissou_fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "seen_non_nested_betas")) %>%
  rename(group = variable) %>%
  mutate(variable  = zissou_fit$zissou_data$x_seen_non_nested %>% colnames())

seeing_non_nested_betas <- zissou_fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "seeing_non_nested_betas")) %>%
  rename(group = variable) %>%
  mutate(variable  = zissou_fit$zissou_data$x_seen_non_nested %>% colnames())

seen_year_species_betas <- zissou_fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "seen_year_species_betas")) %>%
  rename(group = variable) %>%
  mutate(variable = zissou_fit$zissou_data$x_seen_year_species %>% colnames())

seeing_year_species_betas <- zissou_fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "seeing_year_species_betas")) %>%
  # filter(variable == "seeing_year_species_betas") %>%
  rename(group = variable) %>%
  mutate(variable = zissou_fit$zissou_data$x_seen_year_species %>% colnames())

seen_region_cluster_betas <- zissou_fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "seen_region_cluster_betas")) %>%
  # filter(variable == "seen_region_cluster_betas") %>%
  rename(group = variable) %>%
  mutate(variable = zissou_fit$zissou_data$x_seen_region_cluster %>% colnames())

seeing_region_cluster_betas <- zissou_fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "seeing_region_cluster_betas")) %>%
  # filter(variable == "seeing_region_cluster_betas") %>%
  rename(group = variable) %>%
  mutate(variable = zissou_fit$zissou_data$x_seen_region_cluster %>% colnames())

did_betas <- zissou_fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "mpa_effect")) %>%
  # filter(variable == "mpa_effect") %>%
  mutate(group = variable) %>%
  mutate(year = zissou_fit$did_data$year %>% unique())


 betas <- bind_rows(
   seen_non_nested_betas,
   seeing_non_nested_betas,
   seen_year_species_betas,
   seeing_year_species_betas,
   seen_region_cluster_betas,
   seeing_region_cluster_betas,
   did_betas %>% select(-year)
 ) %>%
   as_data_frame()

 non_nested_beta_plot <- betas %>%
   filter(str_detect(group, "non_nested")) %>%
   ggplot() +
   geom_hline(aes(yintercept = 0), linetype = 2, color = "red") +
   geom_pointrange(aes(x = variable,
                       y = estimate,
                       ymin = lower,
                       ymax = upper)) +
   facet_wrap(~group) +
   coord_flip() +
   theme(axis.text.y = element_text(size = 10))


 year_species_effects_plots <- betas %>%
   filter(str_detect(group, "year_species_betas")) %>%
   mutate(year = str_replace_all(variable,"\\D","") %>% as.numeric()) %>%
   mutate(classcode = str_split(variable,'-', simplify = T)[,2]) %>%
   ggplot() +
   geom_hline(aes(yintercept = 0), linetype = 2, color = "red") +
   geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = classcode), alpha = 0.25) +
   geom_line(aes(x = year, y = estimate, color = classcode)) +
   facet_wrap(~group)

  region_cluster_plots <- betas %>%
   filter(str_detect(group, "region_cluster_betas")) %>%
   mutate(cluster = str_replace_all(variable,"\\D","")) %>%
   mutate(region = str_split(variable,'-', simplify = T)[,3]) %>%
   ggplot() +
    geom_hline(aes(yintercept = 0), linetype = 2, color = "red") +
    geom_pointrange(aes(x = region,
                        y = estimate,
                        ymin = lower,
                        ymax = upper, color = cluster)) +
    facet_wrap(~group)



 did_plot <- did_betas %>%
   ggplot() +
   geom_vline(aes(xintercept = 2003), color = 'red', linetype = 2, size = 2) +
 geom_hline(aes(yintercept = 0)) +
   geom_pointrange(aes(
     year,
     y = estimate,
     ymin = lower,
     ymax = upper
   ),
   size = 1.5) +
   labs(x = "Year", y = "Targeted Trend Divergence")

top_species <- base_run$data[[1]]$classcode %>% unique()

cip_data <- pisco_data %>% 
  left_join(site_data, by = c("site","side")) %>% 
  filter(is.na(eventual_mpa) == F) %>% 
  filter(region %in% c("ANA", "SCI","SRI",'SMI'),
         classcode %in% top_species)

# process outcomes

outcomes <- processed_grid %>%
  select(-fishery_effect, -density_ratio) %>% 
  slice(1:samps) %>% 
  unnest() %>% 
  group_by(experiment) %>% 
  mutate(year = 1:length(year)) %>% 
  ungroup() %>% 
  mutate(years_protected = year - year_mpa + 1) %>% 
  mutate(mpa_effect = pmin(mpa_effect,2)) %>% 
  group_by(experiment) %>% 
  mutate(b0 = `no-mpa`[year == min(year)]) %>% 
  mutate(depletion = 1 - `no-mpa`/b0) %>% 
  ungroup() %>% 
  mutate(pop_effect = pmin(1,(`with-mpa` - `no-mpa`) / b0))

eqo <- outcomes %>% 
  filter(year == max(year))

fishery_outcomes <- processed_grid %>%
  select(-mpa_effect, -density_ratio) %>% 
  slice(1:samps) %>% 
  unnest() %>% 
  group_by(experiment) %>% 
  mutate(year = 1:length(year)) %>% 
  ungroup() %>% 
  mutate(years_protected = year - year_mpa + 1) %>% 
  mutate(fishery_effect = pmin(mpa_effect,1)) %>% 
  left_join(outcomes %>% select(experiment, year, depletion), by = c("experiment", "year")) %>% 
  filter(fleet_model != "constant-catch") %>% 
  mutate(msy_effect = pmin(1,(`with-mpa` - `no-mpa`) / msy))

density_ratios <- processed_grid %>%
  select(-fishery_effect) %>% 
  slice(1:samps) %>% 
  unnest() %>% 
  group_by(experiment) %>% 
  mutate(year = 1:length(year)) %>% 
  ungroup() %>% 
  mutate(years_protected = year - year_mpa + 1) %>% 
  left_join(outcomes %>% select(experiment, year, depletion), by = c("experiment", "year"))


fishery_eqo <- fishery_outcomes %>% 
  filter(year == max(year))


facet_labels <- c(
  mpa_size = "Range in MPA",
  depletion = "Pre-MPA Depletion"
)



## ------------------------------------------------------------------------
short_term <- outcomes %>% 
  filter(years_protected == 15, years_protected >= 0)


fishery_short_term <- fishery_outcomes %>% 
  filter(years_protected <= 15, years_protected >= 0)

median_effect <- median(short_term$mpa_effect[short_term$year == max(short_term$year)])

median_pop_effect <- median(short_term$pop_effect[short_term$year == max(short_term$year)])

median_fish_effect <- median(fishery_short_term$mpa_effect[short_term$year == max(short_term$year)], na.rm = TRUE)

median_fish_effect <- median(fishery_short_term$msy_effect[short_term$year == max(short_term$year)], na.rm = TRUE)


pos_runs <- round(mean(short_term$pop_effect[short_term$year == max(short_term$year)] > 0),2)

pos_fish_runs <- round(mean(fishery_short_term$msy_effect[fishery_short_term$year == max(fishery_short_term$year)] > 0, na.rm = TRUE),2)

pos_fish_effect <- round(mean(short_term$f[short_term$year == max(short_term$year)] > 0),2)

small_effect <- round( 100* median(short_term$mpa_effect[short_term$mpa_size <= 0.25 & short_term$depletion <= 0.5]))

min_small_effect <- round( 100* min(short_term$mpa_effect[short_term$mpa_size <= 0.25 & short_term$depletion <= 0.5]))

max_small_effect <- round( 100* max(short_term$mpa_effect[short_term$mpa_size <= 0.25 & short_term$depletion <= 0.5]))

medium_effect <- round(100*median(short_term$mpa_effect[short_term$depletion > 0.5 & short_term$depletion <= 0.75]))

big_effect <- round(100*median(short_term$mpa_effect[short_term$depletion > 0.5 & short_term$depletion > 0.75]))


## ----conservation-effect, fig.cap = "Median (A) and range (B) regional MPA conservation effect (expressed as percent changes in biomass with MPAs relative to what would have occured without MPAs) after 15 years of protection. For (A), X-axes indicate the pre-MPA depletion of the fishery, where depletion is the percentage of unfished biomass that has been removed from the population, and Y-axes is the percent of the population's range encompasssed inside an MPA. For B), y-axes show the regional conservation effect.", include = FALSE----


pop_depletion_plot <- outcomes %>%
  filter(years_protected == 15) %>%
  ggplot() +
  geom_bin2d(aes(depletion, mpa_effect), binwidth = c(.05, .05),
             show.legend = FALSE) +
  scale_fill_viridis(
  option = "A",
  trans = "log10",
  guide = guide_colorbar(frame.colour = "black",
  frame.linewidth = 1,),
  name = "Median Effect"
  ) +
  scale_x_continuous(labels = percent, name = "Pre-MPA Depletion") + 
  scale_y_continuous(labels = percent, name = "")


pop_size_plot <- outcomes %>%
  filter(years_protected == 15) %>%
  ggplot() +
  geom_bin2d(aes(mpa_size, mpa_effect), binwidth = c(.05, .05), show.legend = TRUE) +
  scale_fill_viridis(
  option = "A",
  trans = "log10",
  guide = guide_colorbar(frame.colour = "black",
  frame.linewidth = 1,),
  name = "# of Sims"
  ) +
  scale_x_continuous(labels = percent, name = "Range in MPA") + 
    scale_y_continuous(labels = percent, name = "Conservation Effect") + 
  theme(legend.position = "right")


 pop_facet_effect_plot <- outcomes %>%
  filter(years_protected == 15) %>%
   select(mpa_effect, mpa_size, depletion) %>% 
   gather(measure, value, -mpa_effect) %>% 
   ggplot() +
   geom_bin2d(aes(value, mpa_effect), binwidth = c(.05, .05), show.legend = TRUE) +
  scale_fill_viridis(
  option = "A",
  trans = "log10",
  guide = guide_colorbar(frame.colour = "black",
  frame.linewidth = 1,barheight = unique(10, "lines")),
  name = "# of Sims"
  ) + 
   facet_wrap(~measure, labeller = labeller(measure = facet_labels), strip.position = "bottom") + 
   scale_x_percent(name = "") + 
   scale_y_percent(name = "Conservation Effect")

 
  pop_depletion_and_size_plot <- outcomes %>% 
    filter(years_protected >= 0) %>%
    group_by(experiment) %>% 
    mutate(depletion = plyr::round_any(depletion[years_protected == 0], .05),
           mpa_size = plyr::round_any(mpa_size, .05)) %>% 
    filter(years_protected == 15) %>% 
    group_by(depletion, mpa_size) %>% 
    summarise(median_mpa_effect = median(mpa_effect)) %>% 
    ggplot(aes(depletion, mpa_size, fill = median_mpa_effect)) + 
    geom_tile() + 
    # geom_contour(aes(z = median_mpa_effect)) +
    scale_fill_viridis(labels = percent,
                       guide = guide_colorbar(frame.colour = "black",
                                              frame.linewidth = 1,
                                              barheight = unit(10,"lines")), 
                       name = "Median Effect") +
  labs(x = "Pre-MPA Depletion", y = "Range in MPA") + 
    scale_y_continuous(labels = percent) + 
    scale_x_continuous(labels = percent)
  
  pop_combo_plot <-
    (pop_depletion_and_size_plot + labs(title = "A")) + ((pop_size_plot + labs(title = "B")) / pop_depletion_plot)  + plot_layout(widths = c(1.5, 1)) & theme(
    plot.margin = unit(c(0, 0, 0, 0), units = "lines"),
    axis.text.x = element_text(size = 8),
    legend.box.margin = unit(c(0, 0, 0, 0), units = "lines"),
    axis.text.y = element_text(size = 10)
    )
    
  
  pdf(file =  "figs/conservation-effect.pdf",width = 8, height = 5, useDingbats = TRUE)
print(pop_combo_plot)
dev.off()
  


## ----fishery-effects,fig.cap = "Median (A) and range (B) MPA fishery effects, expressed as the difference in catch with and without MPAs  as a proportion of MSY, after 15 years of protection. For (A), X-axes indicate the pre-MPA depletion of the fishery, where depletion is the percentage of unfished biomass that has been removed from the population, and Y-axes is the percent of the population's range encompasssed inside an MPA. For B), y-axes show the regional conservation effect. Constant-catch scenarios are not included in this plot since by definition catches are equal with or without MPAs", include = FALSE----

msy_depletion_plot <- fishery_outcomes %>%
  filter(years_protected == 15) %>%
  ggplot() +
  geom_bin2d(aes(depletion, msy_effect), binwidth = c(.05, .05),
             show.legend = FALSE) +
  scale_fill_viridis(
  option = "A",
  guide = guide_colorbar(frame.colour = "black",
  frame.linewidth = 1,),
  name = "Median Effect"
  ) +
  scale_x_continuous(labels = percent, name = "Pre-MPA Depletion") + 
  scale_y_continuous(labels = percent, name = "Fishery Effect")


msy_size_plot <- fishery_outcomes %>%
  filter(years_protected == 15) %>%
  ggplot() +
  geom_bin2d(aes(mpa_size, msy_effect), binwidth = c(.05, .05), show.legend = TRUE) +
  scale_fill_viridis(
  option = "A",
  guide = guide_colorbar(frame.colour = "black",
  frame.linewidth = 1,),
  name = "# of Sims"
  ) +
  scale_x_continuous(labels = percent, name = "Range in MPA") + 
    scale_y_continuous(labels = percent, name = "Fishery Effect") + 
  theme(legend.position = "right")

msy_facet_effect_plot <- fishery_outcomes %>%
  filter(years_protected == 15) %>%
   select(msy_effect, mpa_size, depletion) %>% 
   gather(measure, value, -msy_effect) %>% 
   ggplot() +
   geom_bin2d(aes(value, msy_effect), binwidth = c(.05, .05), show.legend = TRUE) +
  scale_fill_viridis(
  option = "A",
  guide = guide_colorbar(frame.colour = "black",
  frame.linewidth = 1,barheight = unique(10, "lines")),
  name = "# of Sims"
  ) + 
   facet_wrap(~measure, labeller = labeller(measure = facet_labels), strip.position = "bottom") + 
   scale_x_percent(name = "") + 
   scale_y_percent(name = "Fishery Effect")


  msy_depletion_and_size_plot <- fishery_outcomes %>% 
    filter(years_protected >= 0,!is.na(msy)) %>%
    group_by(experiment) %>% 
    mutate(depletion = plyr::round_any(depletion[years_protected == 0], .05),
           mpa_size = plyr::round_any(mpa_size, .05)) %>% 
    filter(years_protected == 15) %>% 
    group_by(depletion, mpa_size) %>% 
    summarise(median_mpa_effect = median(msy_effect, na.rm = T)) %>% 
    ggplot(aes(depletion, mpa_size, fill = median_mpa_effect)) + 
    geom_raster(interpolate = TRUE) + 
    # geom_contour(aes(z = median_mpa_effect)) +
    scale_fill_gradient2(midpoint = 0,
                         low = "tomato",
                         high = "steelblue",
                         mid = "white",
                         labels = percent,
                       guide = guide_colorbar(frame.colour = "black",
                                              frame.linewidth = 1,
                                              barheight = unit(10,"lines")), 
                       name = "Median Effect") +
  labs(x = "Pre-MPA Depletion", y = "Range in MPA") + 
    scale_y_continuous(labels = percent) + 
    scale_x_continuous(labels = percent, limits = c(0,NA))
  
  msy_fishery_combo_plot <-
    (msy_depletion_and_size_plot + labs(title = "A")) + ((msy_size_plot + labs(title = "B")) / msy_depletion_plot)  + plot_layout(widths = c(1.5, 1)) & theme(
    plot.margin = unit(c(0, 0, 0, 0), units = "lines"),
    axis.text.x = element_text(size = 8),
    legend.box.margin = unit(c(0, 0, 0, 0), units = "lines"),
    axis.text.y = element_text(size = 10)
    )
    
pdf(file =  "figs/fishery-effect.pdf",width = 8, height = 5, useDingbats = TRUE)
print(msy_fishery_combo_plot)
dev.off()
  


## ----density-ratio-plot--------------------------------------------------


unbiased_dr <- density_ratios %>% 
  group_by(experiment) %>% 
  filter(years_protected == max(years_protected)) %>% 
  ggplot(aes(1 - depletion, pmin(10,true_density_ratio - 1))) + 
  geom_point(size = 2, alpha = 0.75)  + 
  labs(x = "Percent of Unfished Biomass", y = "Percent Difference in Inside/Outside MPA Densities", caption = "Density ratios calculated as density inside MPAs relative to density of identical scenario without MPAs",
       title = "'Unbiased' Density Ratio") + 
  scale_x_percent() + 
  scale_y_percent()

biased_dr <- density_ratios %>% 
  group_by(experiment) %>% 
  filter(years_protected == max(years_protected)) %>% 
  ggplot(aes(1 - depletion, pmin(10,biased_density_ratio - 1))) + 
  geom_point() + 
    labs(x = "Percent of Unfished Biomass", y = "Percent Difference in Inside/Outside MPA Densities", caption = "Density ratios calculated as density inside MPAs relative to density outside MPAs, weighted by distance from MPA border",
         title = "'Biased' Density Ratio") + 
  scale_x_percent() + 
  scale_y_percent()




## ----dr-vs-mpa-effect-plot-----------------------------------------------

unbiased_dr_plot <- density_ratios %>% 
  group_by(experiment) %>% 
  filter(years_protected == max(years_protected)) %>% 
  ggplot(aes(pmin(4,mpa_effect), pmin(4,true_density_ratio - 1))) + 
  geom_point(aes(color = adult_movement),size = 2, alpha = 0.75)  + 
  labs(x = "Percent difference in Total Biomass With MPAs vs Without", y = "Percent Difference in Inside/Outside MPA Densities", caption = "Density ratios calculated as density inside MPAs relative to density of identical scenario without MPAs",
       title = "'SUTVA Satisfied' Density Ratio") + 
  scale_x_percent() + 
  scale_y_percent() + 
  geom_smooth(method = "lm") + 
  geom_abline(aes(slope = 1, intercept = 0)) +
  scale_color_viridis()

biased_dr_plot <- density_ratios %>% 
  group_by(experiment) %>% 
  filter(years_protected == max(years_protected)) %>% 
  ggplot(aes( pmin(4,biased_density_ratio - 1),pmin(4, mpa_effect))) + 
  geom_point(aes(color = adult_movement),size = 2, alpha = 0.75)  + 
    labs(y = "Percent difference in Total Biomass With MPAs vs Without", x = "Percent Difference in Inside/Outside MPA Densities", caption = "Density ratios calculated as density inside MPAs relative to density outside MPAs, weighted by distance from MPA border",
         title = "'SUTVA Violated' Density Ratio") + 
  scale_x_percent() + 
  scale_y_percent() + 
  geom_smooth(method = "lm") + 
  geom_abline(aes(slope = 1, intercept = 0)) +
  scale_color_viridis()




## ------------------------------------------------------------------------
fish_v_fishing <- processed_grid %>%
  select(-density_ratio) %>%
  slice(1:samps) %>%
  unnest() %>%
  rename(fishery_effect = mpa_effect1) %>%
  group_by(experiment) %>%
  mutate(year = 1:length(year)) %>%
  ungroup() %>%
  mutate(years_protected = year - year_mpa + 1) %>%
  left_join(outcomes %>% select(experiment, year, depletion),
            by = c("experiment", "year"))

fish_v_fishing %>%
  filter(years_protected == max(years_protected)) %>%
  ggplot(aes(pmin(2, mpa_effect), pmin(2, fishery_effect))) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_smooth(method = "lm") +
  geom_hex()




## ----ci-map, fig.cap = "Map of study region and sampling locations. Shaded polygons indicate location of MPAs. Points represent sampling locations, and color indicates the number of observations recorded at a given point", include = FALSE----

ci_map <- pisco_data %>% 
  left_join(site_data, by = c("site","side")) %>% 
  filter(is.na(eventual_mpa) == F) %>% 
  filter(region %in% c("ANA", "SCI","SRI",'SMI'),
         classcode %in% top_species) %>% 
  group_by(lat_wgs84, lon_wgs84, site, region, eventual_mpa) %>% 
  count()



ci_map <-  ci_map %>%
  dplyr::mutate(geometry = purrr::map2(lon_wgs84, lat_wgs84, ~ sf::st_point(x = c(.x, .y), dim = 'XY'))) %>%
  ungroup() %>%
  mutate(geometry = sf::st_sfc(geometry, crs = 4326)) %>%
  sf::st_sf()


bbox <- sf::st_bbox(ci_map)

pisco_ci_map_plot <-  ggmap(channel_islands) +
  geom_sf(data = ca_mpas, inherit.aes = FALSE, alpha = 0.5) +
  geom_sf(data = ci_map,inherit.aes = FALSE,
  aes(color = n),
  alpha = 0.75,
  size = 1) +
  coord_sf(xlim = c(bbox['xmin'] - .3, bbox['xmax'] + .06),
  ylim = c(bbox['ymin'] - .1, bbox['ymax'] + .4)) +
  scale_color_gradient(low = "white", high = "orangered",
  guide = guide_colorbar(frame.colour = "black",frame.linewidth = 1,
  barheight = unit(13, "lines")),
  name = "Samples")  + 
  labs(x = "Longitude", y = "Latitude")

pdf(file =  "figs/channel-islands.pdf",width = 8, height = 5, useDingbats = TRUE)
print(pisco_ci_map_plot)
dev.off()




## ------------------------------------------------------------------------
ci_effects <- short_term %>% 
  filter(years_protected == max(years_protected),
         mpa_size  <= 0.25, mpa_size >= 0.15) 

med_ci_effect <- median(ci_effects$mpa_effect)

# rethinking::HPDI(ci_effects$mpa_effect, prob = 0.5)

quart_ci_effect <- quantile(ci_effects$mpa_effect)



## ----raw-trend, fig.cap="Centered and scaled mean annual density of included species (faded lines) and smoothed means of targeted and non-targeted groups, and mean (darker lines) and 95% confidence interval of the mean (ribbon) over time", include = FALSE----
raw_did_plot <- base_run$data[[1]] %>% 
  select(year, targeted,classcode, log_density, any_seen) %>% 
  mutate(density = exp(log_density) * any_seen) %>% 
  group_by(classcode,targeted, year) %>% 
  summarise(md = mean(density)) %>% 
  group_by(classcode,targeted) %>% 
  mutate(smd = (md - mean(md)) / sd(md)) %>% 
  ggplot(aes(year, smd, color = targeted == 1, fill = targeted == 1)) + 
  geom_vline(aes(xintercept = 2003), linetype = 2, color = "red") +
  geom_line(aes(group = interaction(targeted, classcode)), alpha = 0.5) +
  geom_smooth(size = 2) + 
  labs(y = "Centered and Scaled Mean Density", x = "Year") +
    scale_color_manual(values = c("steelblue", "tomato") ,name = "Targeted?") +
      scale_fill_manual(values = c("steelblue", "tomato") ,name = "Targeted?") +
    labs(y = "Centered and Scaled Mean Density", x = "Year")

  
raw_did_plot <- raw_did_plot + 
      scale_color_manual(values = c("steelblue", "tomato"), labels = c("Non-Targeted", "Targeted"), name = "") +
        scale_fill_manual(values = c("steelblue", "tomato"), labels = c("Non-Targeted", "Targeted"), name = "") 

pdf(file =  "figs/pop-trends.pdf",width = 8, height = 5, useDingbats = TRUE)
print(raw_did_plot)
dev.off()





## ----did-plot, fig.cap = "Estimated divergence in biomass densities of targeted and non-targeted fishes throughout the Channel Islands (i.e. integrated across inside and outside of MPAs). MPAs are implemented in 2003 (red dashed line). Estimates are from a regression on log(abundance index), so estimated effects roughly correspond to percentage changes", include = FALSE----

did_plot <- did_plot + 
  scale_y_continuous(name = "~Divergence from Non-Targeted", labels = percent)

pdf(file =  "figs/did_plot.pdf",width = 8, height = 5, useDingbats = TRUE)
print(did_plot)
dev.off()


