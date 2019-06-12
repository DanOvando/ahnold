

## ---- include=FALSE------------------------------------------------------
library(sf)
library(ggmap)
library(viridis)
library(hrbrthemes)
library(scales)
library(patchwork)
# library(rpart)
library(extrafont)
# library(caret)
library(ggsci)
library(rEDM)
library(rstanarm)
library(bayesplot)
library(tidyverse)
library(caret)
library(recipes)
library(scales)
library(spasm)
extrafont::loadfonts()
functions <- list.files(here::here("functions"))

walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions

sim_years <- 50

burn_years <- 1

num_patches <- 50

samps <- 20000

short_frame <- 15

run_name <- "v4.1"

fig_name <- "presentations"

fig_width <- 12

fig_height <- fig_width / 1.333

device <- "pdf"

run_dir <- here::here("results", run_name)

fig_dir <- file.path(run_dir,"figures",fig_name)

if (dir.exists(fig_dir) == FALSE){
  dir.create(fig_dir, recursive = TRUE)
}

write(glue::glue("Figures generated on {Sys.Date()} using results {run_name}"),
      file = file.path(fig_dir,"README.txt"))

experiment_dir <- here::here("results",run_name, "experiments")

load(file = here::here("results",run_name, "rawish_zissou_data.Rdata"))

load(file = here::here("results",run_name, "model_runs.Rdata"))

load(file = here::here("results",run_name,"processed_grid.Rdata"))

load(file = here::here("results",run_name,"simulated_did.Rdata"))

load(file = file.path(run_dir, "abundance_data.Rdata"))


channel_islands <- readRDS(here::here("data","channel_islands_map.RDS"))

ca_mpas <- sf::st_read(here::here("data","MPA_CA_Existing_160301")) %>%
  rmapshaper::ms_simplify() %>%
  sf::st_transform(crs = 4326)



zissou_theme <-
  theme_ipsum(
    base_size = 22,
    axis_title_size = 26,
    strip_text_size = 26
  )

theme_set(zissou_theme)

gc <- guide_colorbar(frame.colour = "black",
                     ticks.colour = "black",
                     barheight = 35)


plot_trans <- "identity"

# filter results ----------------------------------------------------------


colfoo <- function(sim, collapsed = 0.1){
  # sim <- processed_grid$mpa_effect[[1]]
  prepop <- sim$`no-mpa`[sim$mpa_size == 0]

  prepop <- prepop / max(prepop)

  collapse <- any(prepop <= collapsed)

}

processed_grid <- processed_grid %>%
  mutate(collapsed = map_lgl(mpa_effect, colfoo)) %>%
  filter(collapsed == FALSE)

# make examples -----------------------------------------------------------

fish <- create_fish(r0 = 100)

fleet <- create_fleet(fleet_model = "constant-effort", fish = fish,
                      initial_effort = 1, q  = 1)


sizes <- tibble(mpa_size = seq(0,0.75, by = 0.25))

year_mpa <- 20
sizes <- sizes %>%
  mutate(foo = map(mpa_size, ~ run_mpa_experiment(fish = fish,
                                        fleet = fleet,
                                        year_mpa = year_mpa,
                                        mpa_size = .x,
                                        sim_years = 50,
                                        num_patches = 25,
                                        burn_years = 1,
                                        sprinkler = FALSE,
                                        mpa_habfactor = 100,
                                        enviro = NA,
                                        enviro_strength = NA,
                                        rec_driver = 'stochastic',
                                        simseed = 42)))

compare <- sizes %>%
  mutate(delta = map(foo,"outcomes")) %>%
  select(-foo) %>%
  unnest() %>%
  mutate(years_protected = year - year_mpa) %>%
  filter(years_protected <=short_frame)

ribbon <- compare %>%
  select(years_protected, mpa_size, experiment, biomass) %>%
  spread(experiment, biomass) %>%
  group_by(mpa_size) %>%
  mutate(delta = mean(`with-mpa` - `no-mpa`))

labelfoo <- function(x){
  paste0("MPA:", percent(as.numeric(x)))

}


effect_example_plot <- compare %>%
  ggplot() +
  geom_ribbon(
    data = ribbon,
    aes(
      years_protected,
      ymin = `no-mpa`,
      ymax = `with-mpa`,
      fill = delta
    ),
    alpha = 0.5,
    show.legend = FALSE
  ) +
  geom_line(aes(years_protected , biomass, color = experiment), size = 1.5) +
  geom_vline(aes(xintercept = 0), linetype = 2, color = "red") +
  facet_wrap(~ mpa_size, labeller = labeller(mpa_size = labelfoo)) +
  labs(x = "Years with MPA", y = "Biomass") +
  scale_color_aaas(labels = c("Without MPAs", "With MPAs"), name  = '') +
  scale_fill_gradient(low = "white", high = "steelblue",
                      guide = gc)



# fig_width <- 12
#
# fig_height <- fig_width / 1.333
#
#
# savefoo <- function(fig,
#                     device = "pdf",
#                     fig_width = 6,
#                     fig_height = 5) {
#   ggsave(
#     filename =  file.path(fig_dir, paste(fig, device, sep = '.')),
#     plot  = get(fig),
#     width = fig_width,
#     height = fig_height
#   )
#
# }
#
# walk(plots, savefoo, device = device, fig_height = fig_height,
#      fig_width = fig_width)


# make other plots --------------------------------------------------------



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

 doh_did_plot <- did_betas %>%
   filter(year < 2015) %>%
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
   labs(x = "Year", y = "Targeted Trend Divergence") +
  xlim(2000,2017)


top_species <- base_run$data[[1]]$classcode %>% unique()

cip_data <- pisco_data %>%
  left_join(site_data, by = c("site","side")) %>%
  filter(is.na(eventual_mpa) == F) %>%
  filter(region %in% c("ANA", "SCI","SRI",'SMI'),
         classcode %in% top_species)


processed_grid$adult_movement <- (2 * processed_grid$adult_movement) / num_patches

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

# density ratio is being saved as a list of a
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





# run super simple counterfactual -----------------------------------------

year_mpa <- 2002

data <- abundance_data$data[abundance_data$data_source == "pisco"][[1]] %>%
  mutate(density = exp(log_density) * any_seen) %>%
  mutate(targ = case_when(targeted == 1 ~ "Targeted",  TRUE ~ "Non-Targeted"))

sortasimple <- data %>%
  group_by(year, classcode) %>%
  summarise(mean_density = mean(density),
            targ = unique(targ)) %>%
  mutate(mpa = year > year_mpa) %>%
  ungroup() %>%
  mutate(fyear = factor(year)) %>%
  mutate(fyear = relevel(fyear, ref = "2003"))

sortasimple %>%
  filter(targ == "Targeted") %>%
  ggplot(aes(year, mean_density, color = classcode)) +
  geom_point() +
  geom_label(aes(year,mean_density, label = classcode) ) +
  facet_wrap(~targ)


# sortasimple_reg <- stan_glm(log(mean_density + 1e-3) ~ targ*fyear, data = sortasimple, cores = 4)
#
# bayesplot::mcmc_areas(as.matrix(sortasimple_reg),
#                       regex_pars = ":fyear")



sortasimple2_reg <-
  stan_glmer(
    log(mean_density + 1e-3) ~ targ * fyear + (1 |
                                                 classcode),
    data = sortasimple,
    cores = 4
  )

supersimple_plot <- bayesplot::mcmc_areas(as.matrix(sortasimple2_reg),
                      regex_pars = ":fyear") +
  geom_vline(aes(xintercept = 0), linetype = 2, color = "red") +
  scale_y_discrete(labels = c(2000:2001, 2003:2017)) +
  coord_flip() +
  scale_x_percent()


## ------------------------------------------------------------------------
short_term <- outcomes %>%
  filter(years_protected == short_frame, years_protected >= 0)


fishery_short_term <- fishery_outcomes %>%
  filter(years_protected <= short_frame, years_protected >= 0)

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
  filter(years_protected == short_frame) %>%
  ggplot() +
  geom_bin2d(aes(depletion, mpa_effect), binwidth = c(.05, .05),
             show.legend = FALSE) +
  scale_fill_viridis(
  option = "A",
  trans = plot_trans,
  guide = gc,
  name = "Median Effect"
  ) +
  scale_x_continuous(labels = percent, name = "Pre-MPA Depletion") +
  scale_y_continuous(labels = percent, name = "")

fleet_importance_plot <- outcomes %>%
  group_by(experiment) %>%
  filter(years_protected == max(years_protected)) %>%
  ggplot(aes(mpa_effect, fill = fleet_model)) +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  geom_histogram(show.legend = FALSE, color = "black") +
  facet_wrap(~fleet_model) +
  scale_fill_npg() +
  labs(y = "# of Sims") +
  scale_x_percent(name = "MPA Conservation Effect")


# fleet_importance_plot <- outcomes %>%
#   group_by(experiment) %>%
#   filter(years_protected == max(years_protected)) %>%
#   ggplot(aes(mpa_effect, fill = fleet_model)) +
#   geom_vline(aes(xintercept = 0), linetype = 2) +
#   geom_histogram(show.legend = FALSE, color = "black", aes(alpha = mpa_effect)) +
#   facet_wrap(~fleet_model) +
#   scale_fill_npg() +
#   labs(y = "# of Sims") +
#   scale_x_percent(name = "MPA Conservation Effect")
#


cc_outcomes <- outcomes %>%
  group_by(experiment) %>%
  filter(years_protected == max(years_protected), fleet_model == "constant-catch") %>%
  mutate(conservation_loss = as.factor(mpa_effect < 0)) %>%
  select(-mpa_effect,-pop_effect)


# negative_tree = train(
#   conservation_loss ~  adult_movement + larval_movement + size_limit + f_v_m + mpa_size + factor(density_dependence_form) + density_movement_modifier,
#   method = "rpart",
#   data =  cc_outcomes
# )


# class(rpart.plot::rpart.plot(negative_tree$finalModel))

cc_importance_plot <- outcomes %>%
  group_by(experiment) %>%
  filter(years_protected == max(years_protected), fleet_model == "constant-catch", mpa_effect < 0) %>%
  ggplot(aes(mpa_size, adult_movement)) +
  geom_point(alpha = 0.5, size = 4) +
  geom_density2d(size = 2, color = "tomato",alpha = 0.75) +
  scale_y_percent(name = "Adult Movement") +
  scale_x_percent(name = "MPA Size")


pop_size_plot <- outcomes %>%
  filter(years_protected == short_frame) %>%
  ggplot() +
  geom_bin2d(aes(mpa_size, mpa_effect), binwidth = c(.05, .05), show.legend = TRUE) +
  scale_fill_viridis(
  option = "A",
  trans = plot_trans,
  guide = gc,
  name = "# of Sims"
  ) +
  scale_x_continuous(labels = percent, name = "Range in MPA") +
    scale_y_continuous(labels = percent, name = "Conservation Effect") +
  theme(legend.position = "right")


 pop_facet_effect_plot <- outcomes %>%
  filter(years_protected == short_frame) %>%
   select(mpa_effect, mpa_size, depletion) %>%
   gather(measure, value, -mpa_effect) %>%
   ggplot() +
   geom_bin2d(aes(value, mpa_effect), binwidth = c(.05, .05), show.legend = TRUE) +
  scale_fill_viridis(
  option = "A",
  trans = plot_trans,
  guide = gc,
  name = "# of Sims"
  ) +
   facet_wrap(~measure, labeller = labeller(measure = facet_labels), strip.position = "bottom") +
   scale_x_percent(name = "") +
   scale_y_percent(name = "Conservation Effect")

 pop_log10_facet_effect_plot <- outcomes %>%
   group_by(experiment) %>%
   filter(years_protected == max(years_protected)) %>%
   ungroup() %>%
   select(mpa_effect, mpa_size, depletion) %>%
   gather(measure, value, -mpa_effect) %>%
   ggplot() +
   geom_bin2d(aes(value, mpa_effect), binwidth = c(.05, .05), show.legend = TRUE) +
   scale_fill_viridis(
     option = "A",
     trans = "log10",
     guide = gc,
     name = "# of Sims"
   ) +
   facet_wrap(~measure, labeller = labeller(measure = facet_labels), strip.position = "bottom") +
   scale_x_percent(name = "") +
   scale_y_percent(name = "Conservation Effect")


 pop_f_facet_effect_plot <- outcomes %>%
   filter(years_protected == short_frame) %>%
   select(mpa_effect, mpa_size, f_v_m) %>%
   gather(measure, value,-mpa_effect) %>%
   ggplot() +
   geom_bin2d(aes(value, mpa_effect),
              binwidth = c(.05, .05),
              show.legend = TRUE) +
   scale_fill_viridis(
     option = "A",
     trans = plot_trans,
     guide = gc,
     name = "# of Sims"
   ) +
   facet_wrap( ~ measure,
               strip.position = "bottom") +
   scale_x_percent(name = "") +
   scale_y_percent(name = "Conservation Effect")



 pop_depletion_and_size_plot <- outcomes %>%
   filter(years_protected >= 0) %>%
   group_by(experiment) %>%
   mutate(depletion = plyr::round_any(depletion[years_protected == 0], .05),
          mpa_size = plyr::round_any(mpa_size, .05)) %>%
   filter(years_protected == short_frame) %>%
   group_by(depletion, mpa_size) %>%
   summarise(median_mpa_effect = median(mpa_effect)) %>%
   ggplot(aes(depletion, mpa_size, fill = median_mpa_effect)) +
   geom_tile() +
   # geom_contour(aes(z = median_mpa_effect)) +
   scale_fill_viridis(labels = percent,
                      guide = gc,
                      name = "Median Effect") +
   labs(x = "Pre-MPA Depletion", y = "Range in MPA") +
   scale_y_continuous(labels = percent) +
   scale_x_continuous(labels = percent)

 eq_pop_depletion_and_size_plot <- outcomes %>%
   filter(years_protected >= 0) %>%
   group_by(experiment) %>%
   mutate(depletion = plyr::round_any(depletion[years_protected == max(years_protected)], .05),
          mpa_size = plyr::round_any(mpa_size, .05)) %>%
   filter(years_protected == max(years_protected)) %>%
   group_by(depletion, mpa_size) %>%
   summarise(median_mpa_effect = median(mpa_effect)) %>%
   ggplot(aes(depletion, mpa_size, fill = median_mpa_effect)) +
   geom_tile() +
   # geom_contour(aes(z = median_mpa_effect)) +
   scale_fill_viridis(labels = percent,
                      guide = gc,
                      name = "Median Pop. Effect") +
   labs(x = "Equilibrium Depletion", y = "Range in MPA") +
   scale_y_continuous(labels = percent) +
   scale_x_continuous(labels = percent)


 eq_pop_f_and_size_plot <- outcomes %>%
   filter(years_protected >= 0) %>%
   group_by(experiment) %>%
   mutate(f = plyr::round_any(f_v_m[years_protected == 0], .1),
          mpa_size = plyr::round_any(mpa_size, .05)) %>%
   filter(years_protected == max(years_protected)) %>%
   group_by(f, mpa_size) %>%
   summarise(median_mpa_effect = median(mpa_effect)) %>%
   ggplot(aes(f, mpa_size, fill = median_mpa_effect)) +
   geom_tile() +
   # geom_contour(aes(z = median_mpa_effect)) +
   scale_fill_viridis(labels = percent,
                      guide = gc,
                      name = "Median Effect") +
   labs(x = "EQ F/M", y = "Range in MPA") +
   scale_y_continuous(labels = percent) +
   scale_x_continuous()



  pop_depletion_and_size_plot <- outcomes %>%
    filter(years_protected >= 0) %>%
    group_by(experiment) %>%
    mutate(depletion = plyr::round_any(depletion[years_protected == 0], .05),
           mpa_size = plyr::round_any(mpa_size, .05)) %>%
    filter(years_protected == short_frame) %>%
    group_by(depletion, mpa_size) %>%
    summarise(median_mpa_effect = median(mpa_effect)) %>%
    ggplot(aes(depletion, mpa_size, fill = median_mpa_effect)) +
    geom_tile() +
    # geom_contour(aes(z = median_mpa_effect)) +
    scale_fill_viridis(labels = percent,
                       guide = gc,
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


## ----fishery-effects,fig.cap = "Median (A) and range (B) MPA fishery effects, expressed as the difference in catch with and without MPAs  as a proportion of MSY, after 15 years of protection. For (A), X-axes indicate the pre-MPA depletion of the fishery, where depletion is the percentage of unfished biomass that has been removed from the population, and Y-axes is the percent of the population's range encompasssed inside an MPA. For B), y-axes show the regional conservation effect. Constant-catch scenarios are not included in this plot since by definition catches are equal with or without MPAs", include = FALSE----

msy_depletion_plot <- fishery_outcomes %>%
  filter(years_protected == short_frame) %>%
  ggplot() +
  geom_bin2d(aes(depletion, msy_effect), binwidth = c(.05, .05),
             show.legend = FALSE) +
  scale_fill_viridis(
  option = "A",
  guide = gc,
  name = "Median Effect"
  ) +
  scale_x_continuous(labels = percent, name = "Pre-MPA Depletion") +
  scale_y_continuous(labels = percent, name = "Fishery Effect")


msy_size_plot <- fishery_outcomes %>%
  filter(years_protected == short_frame) %>%
  ggplot() +
  geom_bin2d(aes(mpa_size, msy_effect), binwidth = c(.05, .05), show.legend = TRUE) +
  scale_fill_viridis(
  option = "A",
  guide = gc,
  name = "# of Sims"
  ) +
  scale_x_continuous(labels = percent, name = "Range in MPA") +
    scale_y_continuous(labels = percent, name = "Fishery Effect") +
  theme(legend.position = "right")

msy_facet_effect_plot <- fishery_outcomes %>%
  filter(years_protected == short_frame) %>%
   select(msy_effect, mpa_size, depletion) %>%
   gather(measure, value, -msy_effect) %>%
   ggplot() +
   geom_bin2d(aes(value, msy_effect), binwidth = c(.05, .05), show.legend = TRUE) +
  scale_fill_viridis(
  option = "A",
  guide = gc,
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
    filter(years_protected == short_frame) %>%
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
                       guide = gc,
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



  catch_depletion_plot <- fishery_outcomes %>%
    filter(years_protected == short_frame) %>%
    ggplot() +
    geom_bin2d(aes(depletion, fishery_effect), binwidth = c(.05, .05),
               show.legend = FALSE) +
    scale_fill_viridis(
      option = "A",
      guide = gc,
      name = "Median Effect"
    ) +
    scale_x_continuous(labels = percent, name = "Pre-MPA Depletion") +
    scale_y_continuous(labels = percent, name = "Fishery Effect")


  catch_size_plot <- fishery_outcomes %>%
    filter(years_protected == short_frame) %>%
    ggplot() +
    geom_bin2d(aes(mpa_size, fishery_effect), binwidth = c(.05, .05), show.legend = TRUE) +
    scale_fill_viridis(
      option = "A",
      guide = gc,
      name = "# of Sims"
    ) +
    scale_x_continuous(labels = percent, name = "Range in MPA") +
    scale_y_continuous(labels = percent, name = "Fishery Effect") +
    theme(legend.position = "right")

  catch_facet_effect_plot <- fishery_outcomes %>%
    filter(years_protected == short_frame) %>%
    select(fishery_effect, mpa_size, depletion) %>%
    gather(measure, value, -fishery_effect) %>%
    ggplot() +
    geom_bin2d(aes(value, fishery_effect), binwidth = c(.05, .05), show.legend = TRUE) +
    scale_fill_viridis(
      option = "A",
      guide = gc,
      name = "# of Sims"
    ) +
    facet_wrap(~measure, labeller = labeller(measure = facet_labels), strip.position = "bottom") +
    scale_x_percent(name = "") +
    scale_y_percent(name = "Fishery Effect")


  catch_depletion_and_size_plot <- fishery_outcomes %>%
    filter(years_protected >= 0,!is.na(msy)) %>%
    group_by(experiment) %>%
    mutate(depletion = plyr::round_any(depletion[years_protected == 0], .05),
           mpa_size = plyr::round_any(mpa_size, .05)) %>%
    filter(years_protected == short_frame) %>%
    group_by(depletion, mpa_size) %>%
    summarise(median_mpa_effect = median(fishery_effect, na.rm = T)) %>%
    ggplot(aes(depletion, mpa_size, fill = median_mpa_effect)) +
    geom_raster(interpolate = TRUE) +
    # geom_contour(aes(z = median_mpa_effect)) +
    scale_fill_gradient2(midpoint = 0,
                         low = "tomato",
                         high = "steelblue",
                         mid = "white",
                         labels = percent,
                         guide = gc,
                         name = "Median Effect") +
    labs(x = "Pre-MPA Depletion", y = "Range in MPA") +
    scale_y_continuous(labels = percent) +
    scale_x_continuous(labels = percent, limits = c(0,NA))

  catch_fishery_combo_plot <-
    (catch_depletion_and_size_plot + labs(title = "A")) + ((catch_size_plot + labs(title = "B")) / catch_depletion_plot)  + plot_layout(widths = c(1.5, 1)) & theme(
      plot.margin = unit(c(0, 0, 0, 0), units = "lines"),
      axis.text.x = element_text(size = 8),
      legend.box.margin = unit(c(0, 0, 0, 0), units = "lines"),
      axis.text.y = element_text(size = 10)
    )




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
  ggplot(aes(y = pmin(4,mpa_effect), x = pmin(4,true_density_ratio - 1))) +
  geom_point(aes(color = adult_movement),size = 2, alpha = 0.75)  +
  labs(y = "True % Biomass Difference", x = "% Difference in Inside/Outside Biomass") +
  scale_x_percent() +
  scale_y_percent() +
  geom_smooth(method = "lm") +
  geom_abline(aes(slope = 1, intercept = 0)) +
  scale_color_viridis(guide = gc, name = "Adult Movement", labels = percent)


baci_plot <- density_ratios %>%
  group_by(experiment) %>%
  filter(years_protected == max(years_protected)) %>%
  ggplot(aes(y = mpa_density - nompa_density, x = baci)) +
  geom_point(aes(color = adult_movement),size = 2, alpha = 0.75)  +
  labs(y = "True Change in Density", x = "BACI Estimate of Density Change") +
  geom_smooth(method = "lm") +
  geom_abline(aes(slope = 1, intercept = 0)) +
  scale_color_viridis(guide = gc, name = "Adult Movement", labels = percent)

biased_dr_plot <- density_ratios %>%
  group_by(experiment) %>%
  filter(years_protected == max(years_protected)) %>%
  ggplot(aes(x =  pmin(4,biased_density_ratio - 1),y = pmin(4, mpa_effect))) +
  geom_point(aes(color = adult_movement),size = 2, alpha = 0.75)  +
  labs(y = "True % Biomass Difference", x = "% Difference in Inside/Outside Biomass") +
  scale_x_percent() +
  scale_y_percent() +
  geom_smooth(method = "lm") +
  geom_abline(aes(slope = 1, intercept = 0)) +
  scale_color_viridis(guide = gc, name = "Adult Movement", labels = percent)




## ------------------------------------------------------------------------
fish_v_fishing <- processed_grid %>%
  select(-density_ratio) %>%
  slice(1:samps) %>%
  unnest() %>%
  group_by(experiment) %>%
  mutate(pop_effect = pmin(1,(`with-mpa` - `no-mpa`) / `no-mpa`[year == min(year)])) %>%
  rename(fishery_effect = mpa_effect1) %>%
  group_by(experiment) %>%
  mutate(year = 1:length(year)) %>%
  ungroup() %>%
  mutate(years_protected = year - year_mpa + 1) %>%
  left_join(outcomes %>% select(experiment, year, depletion),
            by = c("experiment", "year")) #%>%
  # left_join(fishery_outcomes %>% select(experiment, year, msy_effect))

fish_v_fishing %>%
  group_by(experiment) %>%
  filter(years_protected == max(years_protected)) %>%
  ungroup() %>%
  ggplot(aes(pmin(2,fishery_effect), pmin(2,mpa_effect))) +
  geom_bin2d() +
  scale_fill_viridis(guide = gc)


labelfoo <- function(x){

  paste0("%MPA:",x)

}

fish_v_fishing_plot <- fish_v_fishing %>%
  group_by(experiment) %>%
  filter(years_protected == max(years_protected),
         fleet_model != "constant-catch") %>%
  ungroup() %>%
  # mutate(mpa_bins = cut_width(100*mpa_size, width = 25)) %>%
  mutate(mpa_bins = cut(100*mpa_size, breaks = seq(0,100, by = 25))) %>%
  ggplot(aes(x =pmin(1, pop_effect), y = pmin(1, fishery_effect))) +
  geom_hline(aes(yintercept = 0), size = 1) +
  geom_vline(aes(xintercept = 0), size = 1) +
  # geom_smooth(method = "lm") +
  geom_point(alpha = 0.25, aes(fill = mpa_bins), shape = 21, show.legend = FALSE, size = 2) +
  # geom_hex(bins = 10) +
  # geom_abline(aes(slope = 1, intercept = 0), color = "red", size = 1) +
  geom_abline(aes(slope = -1, intercept = 0), color = "steelblue", size = 1) +
  geom_smooth(aes(color = mpa_bins),method = "lm", linetype = 1, se = TRUE) +

  scale_x_percent(name = "% of Unfished Biomass Recovered") +
  scale_y_percent(name = "% Change in Catch") +
  # scale_(labels = percent, guide = gc) +
  scale_color_brewer(palette = "Reds",labels = percent, guide = gc) +
  scale_fill_brewer(palette = "Reds",labels = percent, guide = gc) +
  facet_wrap(~mpa_bins,
             labeller = labeller(mpa_bins = labelfoo))




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

ci_mpa_plot <-  ggmap(channel_islands) +
  geom_sf(data = ca_mpas, inherit.aes = FALSE, alpha = 0.5) +
  coord_sf(xlim = c(bbox['xmin'] - .3, bbox['xmax'] + .06),
           ylim = c(bbox['ymin'] - .1, bbox['ymax'] + .4)) +
  labs(x = "Longitude", y = "Latitude")


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

a = outcomes %>%
  group_by(experiment) %>%
  filter(years_protected == max(years_protected)) %>%
  ungroup() %>%
  select(
    -experiment,
    -msy,
    -b_msy,
    -year,
    -mpa_size1,-`no-mpa`,
    -`with-mpa`,
    -years_protected,
    -depletion,
    -scientific_name,
    -pop_effect,
    -b0,
    -collapsed
    ) %>%
  na.omit()


d <- recipes::recipe(mpa_effect ~ ., data = a) %>%
  step_dummy(all_nominal())

d <- prep(d, data = a) %>%
  juice() %>%
  mutate(rando = rnorm(nrow(.)))


b <- caret::train(mpa_effect ~ .,
           data = d,
           importance = "impurity_corrected",
           method = "ranger")


imp <- b$finalModel$variable.importance

i <- imp %>%
  broom::tidy() %>%
  mutate(names = fct_reorder(names, x))


importance_plot <- i %>%
  filter(x >= x[names == "rando"]) %>%
  top_n(15,x) %>%
  ggplot(aes(names,x)) +
  geom_col(fill = "steelblue", color = "black") +
  coord_flip() +
  labs(x = "", y = "Importance to MPA Conservation Effect")

# ggsave(
#   filename =  file.path(fig_dir, "channel-islands.pdf"),
#   plot  = pisco_ci_map_plot,
#   width = fig_width,
#   height = fig_height,
#   useDingbats = TRUE
# )

# pdf(
#   file =  file.path(fig_dir, "channel-islands.pdf")
#   ,
#   width = 8,
#   height = 5,
#   useDingbats = TRUE
# )
# print(pisco_ci_map_plot)
# dev.off()



## ------------------------------------------------------------------------
ci_effects <- short_term %>%
  filter(years_protected == max(years_protected),
         mpa_size  <= 0.25, mpa_size >= 0.15)

med_ci_effect <- median(ci_effects$mpa_effect)

# rethinking::HPDI(ci_effects$mpa_effect, prob = 0.5)

quart_ci_effect <- quantile(ci_effects$mpa_effect)



# plot the actual data ---------------------------------------------------------------

used_data <- base_run$data[[1]] %>%
  left_join(life_history_data %>% select(classcode, commonname), by = "classcode")


species_plot <- used_data %>%
  group_by(commonname) %>%
  summarise(obs = sum(any_seen),
            targeted = unique(targeted)) %>%
  mutate(commonname = fct_reorder(commonname, obs)) %>%
  ggplot(aes(commonname,obs, fill = targeted == 1)) +
  geom_col(color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("steelblue", "tomato"), labels = c("Non-Targeted", "Targeted"), name = "") +
  labs(x= "", y = "Observations")



## ----raw-trend, fig.cap="Centered and scaled mean annual density of included species (faded lines) and smoothed means of targeted and non-targeted groups, and mean (darker lines) and 95% confidence interval of the mean (ribbon) over time", include = FALSE----
raw_did_plot <- base_run$data[[1]] %>%
  select(year, targeted,classcode, log_density, any_seen, region) %>%
  mutate(density = exp(log_density) * any_seen) %>%
  # group_by(region, classcode) %>%
  # mutate(density = scale(density)) %>%
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


# ggsave(
#   filename =  file.path(fig_dir, "pop-trends.pdf"),
#   plot  = raw_did_plot,
#   width = fig_width,
#   height = fig_height,
#   useDingbats = TRUE
# )

# pdf(
#   file =  file.path(fig_dir, "pop-trends.pdf"),
#   width = 8,
#   height = 5,
#   useDingbats = TRUE
# )
# print(raw_did_plot)
# dev.off()




## ----did-plot, fig.cap = "Estimated divergence in biomass densities of targeted and non-targeted fishes throughout the Channel Islands (i.e. integrated across inside and outside of MPAs). MPAs are implemented in 2003 (red dashed line). Estimates are from a regression on log(abundance index), so estimated effects roughly correspond to percentage changes", include = FALSE----

did_plot <- did_plot +
  scale_y_continuous(name = "~Divergence from Non-Targeted", labels = percent, limits = c(-.75,1))

doh_did_plot <- doh_did_plot +
  scale_y_continuous(name = "~Divergence from Non-Targeted", labels = percent, limits = c(-.75,1))

# EDM ---------------------------------------------------------------------



cip_data <- pisco_data %>%
  left_join(site_data, by = c("site","side")) %>%
  filter(is.na(eventual_mpa) == F) %>%
  filter(region %in% c("ANA", "SCI","SRI",'SMI'),
         classcode %in% top_species)

dat <- cip_data %>%
  group_by(broadtrophic, region, site, side, transect, year, eventual_mpa) %>%
  summarise(density = sum(density_g_m2)) %>%
  group_by(broadtrophic, region, site, year) %>%
  summarise(mean_density = mean(density, na.rm = T)) %>%
  group_by(broadtrophic,region, year) %>%
  summarise(regional_density = mean(mean_density, na.rm = T)) %>%
  ungroup() %>%
  spread(broadtrophic, regional_density) %>%
  group_by(region) %>%
  arrange(year) %>%
  mutate(lag4_carnivore = lag(carnivore,4),
         lag4_piscovore = lag(piscivore,4)) %>%
  na.omit() %>%
  ungroup() %>%
  filter(year > 1999)

dat <- dat %>%
  select(-lag4_carnivore,-lag4_piscovore)%>%
  group_by(region)%>%
  mutate_at(vars(carnivore:planktivore),funs(norm=(.-mean(.))/sd(.)))->dat

group_trends_plot <- dat %>%
  select(region,year,contains("norm"))%>%
  gather("group","density",-region,-year)%>%
  ggplot(aes(year,density,col=group))+
  geom_line()+
  scale_color_locuszoom(name="Trophic Troup",labels=c("carnivore","herbivore","piscivore","planktivore"))+
  labs(x="year",y="Centered and Scaled Density")+
  facet_wrap(~region)

## split data by region to analye timeseries from each island separately
ana.dat <- dat %>% filter(region=="ANA")
sci.dat <- dat %>% filter(region=="SCI")
smi.dat <- dat %>% filter(region=="SMI")
sri.dat <- dat %>% filter(region=="SRI")

datnest <- dat %>% group_by(region) %>% nest()

datnorm <- dat %>% select(region,year,contains('norm')) %>% ungroup()

# have to record the segments corresponding to each "replicate" so simplex algorithm does not try to make predictions crossing time barriers
segs <- datnorm %>% mutate(ind=row_number()) %>% group_by(region) %>% summarise(first=first(ind),last=last(ind)) %>%
  select(-region)

var_names <- c("carnivore_norm","herbivore_norm","piscivore_norm","planktivore_norm")


regions.combined.simp.list <- map(var_names,function(x){
  temp <- datnorm %>% ungroup() %>% select(matches(x)) %>% as.data.frame()
  out <- simplex(as.numeric(temp[,1]),E=1:10,lib=as.matrix(segs),silent=T) %>%
    mutate(trophic=x)
  out
})


embed_plot <- bind_rows(regions.combined.simp.list) %>%
  ggplot(aes(E,rho,color=trophic))+
  geom_line(size=2)+
  facet_wrap(~trophic,nrow=2,scales="free_y")+
  geom_hline(yintercept = 0,color="black")+
  labs(x="Embedding Dimension (E)",y=expression(paste("Skill, ",rho)))+
  scale_x_continuous(breaks=seq(0,12,by=2))+
  scale_color_locuszoom()+
  guides(color=F)

tempE <- 8

temp <- ccm(datnorm,lib=as.matrix(segs),pred=as.matrix(segs),E=tempE,lib_column= 'carnivore_norm',target_column = 'herbivore_norm',lib_sizes = c(10,25,50,75),num_samples=100,replace=T,silent=T,RNGseed = 41389)

c_xmap_h <- temp %>%
  group_by(lib_size)%>%
  summarise(rhomean=mean(rho,na.rm=T),upper=quantile(rho, 0.975),lower=quantile(rho, 0.025))%>%
  ungroup()%>%
  ggplot(aes(lib_size,rhomean))+
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3,fill="red")+
  geom_line(color="darkorchid3")+
  geom_hline(aes(yintercept = 0), linetype = 2) +
  labs(x="",y=expression(paste(rho, " (predictive skill)")),title="Carnivores on Herbivores")

# cross map herbivore to carnivore
# inspect the output of simplex from the previous step and use the best embedding dimension (highest rho) for the carnivore time series

tempE <- 6

temp <- ccm(datnorm,lib=as.matrix(segs),pred=as.matrix(segs),E=tempE,lib_column= 'herbivore_norm',target_column = 'carnivore_norm',lib_sizes = c(10,25,50,75),num_samples=100,replace=T,silent=T,RNGseed = 41389)

h_xmap_c <- temp %>%
  group_by(lib_size)%>%
  summarise(rhomean=mean(rho,na.rm=T),upper=quantile(rho, 0.975),lower=quantile(rho, 0.025))%>%
  ungroup()%>%
  ggplot(aes(lib_size,rhomean))+
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3,fill="red")+
  geom_line(color="darkorchid3")+
  geom_hline(aes(yintercept = 0), linetype = 2) +
  labs(x="Library Size",y=expression(paste(rho, " (predictive skill)")),title="Herbivores on Carnivores")

cross_map_plot <- h_xmap_c + c_xmap_h



# environment effects and MPA only-------------------------------------------------------------------------

facet_labels <- c(
  `TRUE` = "MPA",
  `FALSE` = "Non-MPA"
)


species_mpa_trend_plot <- cip_data %>%
  group_by(year, classcode, eventual_mpa, targeted) %>%
  summarise(mean_density = mean(density_g_m2, na.rm = T)) %>%
  ungroup() %>%
  mutate(targeted = targeted == 1) %>%
  group_by(classcode, eventual_mpa) %>%
  mutate(cs_density = (mean_density - mean(mean_density)) / sd(mean_density)) %>%
  ungroup() %>%
  rename(`Eventual MPA?` = eventual_mpa) %>%
  ggplot(aes(year, cs_density, color = targeted, fill = targeted)) +
  geom_vline(aes(xintercept = 2003), linetype = 2, color = "red") +
  geom_line(aes(group = interaction(classcode, targeted)),alpha = 0.5) +
  geom_smooth() +
  scale_color_manual(values = c("steelblue", "tomato") ,name = "", labels = c("Non-Targeted","Targeted")) +
  scale_fill_manual(values = c("steelblue", "tomato") ,name = "",labels = c("Non-Targeted","Targeted")) +

  facet_wrap(~`Eventual MPA?`, labeller = labeller(`Eventual MPA?` = facet_labels)) +
  labs(y = "Centered and Scaled Mean Density", x = "Year")

species_mpa_trend_plot



mpa_did_plot <- mpa_run$did_plot[[1]] +
  labs(x = "Year", y = "Estimate of Regional Effect",
       caption = "")

temp_trends_plot <- cip_data %>%
  mutate(year_month = year + (month / 12 - .1)) %>%
  filter(region %in% c("ANA", "SCI","SRI","SMI")) %>%
  group_by(year,region) %>%
  summarise(mean_temperature = mean(mean_temp)) %>%
  ggplot(aes(year, mean_temperature, color = region)) +
  geom_line(size = 1) +
  scale_color_locuszoom(name = "") +
  scale_y_continuous(name =bquote("Mean Temperature"~(degree~ "C"))) +
  labs(x = "Year")



# kfm example -------------------------------------------------------------

kfm_did_plot <- kfm_run$did_plot[[1]] +
  scale_y_continuous(name = "~Divergence from Non-Targeted", labels = percent) +
  labs(caption = "")


# regional catches --------------------------------------------------------

cdfw_catches <-
  read_csv(file = here::here("processed_data", 'cdfw-catches.csv')) %>%
  group_by(sci_name, year) %>%
  summarise(catch = sum(pounds_caught, na.rm = T)) %>%
  left_join(
    life_history_data %>% select(taxa, classcode) %>% mutate(taxa = tolower(taxa)),
    by = c("sci_name" = "taxa")
  ) %>%
  filter(!is.na(classcode)) %>%
  left_join(life_history_data %>% select(classcode, commonname), by = "classcode") %>%
  filter(classcode %in% abundance_data$data[[1]]$classcode)

catch_trend_plot <- cdfw_catches %>%
  group_by(classcode) %>%
  mutate(catch = catch * 0.000453592) %>%
  mutate(scatch = scale(catch),
         ncatch = length(catch),
         mcatch = mean(catch)) %>%
  ungroup() %>%
  filter(ncatch > 10) %>%
  ggplot(aes(year, scatch)) +
  geom_line(show.legend = TRUE, aes(size = mcatch)) +
  facet_wrap(~commonname) +
  theme_minimal() +
  labs(y = "Centered and Scaled Catch", x = "Year") +
  scale_size_continuous(name = "Mean Catch (tons)")




# simulation testing ------------------------------------------------------

pisco_sim_plot <- pisco_performance$mixed_effect_did +
  labs(title = "", x= "Year") +
  scale_y_percent(name = "MPA Effect")


simple_sim_plot <- simple_performance$mixed_effect_did +
  labs(title = "", x= "Year") +
  scale_y_percent(name = "MPA Effect")



# case studies ------------------------------------------------------------

fish <- create_fish(r0 = 100, adult_movement = 2, larval_movement = 10)

fleet <- create_fleet(fish = fish,
                      initial_effort = 10,
                      q = 1,
                      target_catch = 25,
                      fleet_model = "constant-effort")

experiment <- spasm::mpa_counterfactual(
  fish = fish,
  fleet = fleet,
  year_mpa = 25,
  mpa_size = .5,
  sim_years = 50,
  burn_years = 1,
  num_patches = 50,
  random_mpas = FALSE,
  min_size = 0.1
)

raw <- experiment$raw_outcomes %>%
  select(year,
         patch,
         biomass,
         biomass_caught,
         profits,
         effort,
         experiment,
         mpa) %>%
  rename(catch = biomass_caught) %>%
  gather(metric, value, -year,-patch,-experiment,-mpa) %>%
  group_by(year, patch, metric, experiment) %>%
  summarise(value = sum(value),
            mpa = unique(mpa)) %>%
  group_by(year, patch) %>%
  mutate(mpa = unique(mpa[experiment == "with-mpa"])) %>%
  ungroup() %>%
  group_by(metric) %>%
  mutate(value = value / max(value)) %>%
  ungroup() %>%
  spread(experiment, value) %>%
  ungroup() %>%
  mutate(delta = `with-mpa` / `no-mpa`,
         ref = 1)

doh <- raw %>%
  filter(metric %in% c("biomass","effort")) %>%
  ggplot() +
  geom_col(
    aes(patch, `with-mpa`, fill = mpa),
    width = 1
  ) +
  geom_line(aes(patch, `no-mpa`), size = 2,
            color = "black") +
  facet_wrap(~ metric) +
  coord_polar() +
  gganimate::transition_time(year) +
  gganimate::ease_aes('linear') +
  labs(title = 'Year: {frame_time}', x = "'",  y = "") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top") +
  scale_fill_npg(labels = c("Witout MPAs","With MPAs"), name = "")


file.path(fig_dir,"simple_cs.gif")

gganimate::anim_save(animation = doh,filename = file.path(fig_dir,"simple_cs.gif")
)

raw_sum <- raw %>%
  group_by(year, metric) %>%
  summarise(no_mpa = sum(`no-mpa`),
            with_mpa = sum(`with-mpa`)) %>%
  ungroup() %>%
  gather(experiment, value, no_mpa:with_mpa) %>%
  filter(year > 20)

ribbon <- raw_sum %>%
  group_by(year, metric) %>%
  summarise(ymin = min(value),
            ymax = max(value),
            delta = value[experiment == "with_mpa"] - value[experiment == "no_mpa"]) %>%
  ungroup() %>%
  group_by(metric) %>%
  mutate(max_delta = last(delta)) %>%
  ungroup()

simple_cs_plot <- raw_sum %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(year, value, color = experiment), size = 1.5) +
  geom_ribbon(data = ribbon, aes(
    year,
    ymin = ymin,
    ymax = ymax,
    fill = max_delta
  ),
  alpha = 0.5,
  show.legend = FALSE) +
  facet_wrap( ~ metric, scales = "free_y",
              labeller = labeller(metric = tools::toTitleCase)) +
  scale_fill_gradient(low = "tomato", high = "steelblue") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.spacing = unit(1,"lines")) +
  labs(x = "Year") +
  scale_color_npg(labels = c("Without MPAs","With MPAs"), name = '')


# trickier case study

  fish <-
    create_fish(r0 = 100,
                adult_movement = 2,
                larval_movement = 10,
                density_dependence_form = 1)

  fleet <- create_fleet(fish = fish,
                        initial_effort = 10,
                        q = 1,
                        target_catch = 25,
                        fleet_model = "open-access",
                        max_cr_ratio = 0.5,
                        effort_allocation = "profit-gravity"
                        )

  experiment <- spasm::mpa_counterfactual(
    fish = fish,
    fleet = fleet,
    year_mpa = 25,
    mpa_size = .5,
    sim_years = 50,
    burn_years = 1,
    num_patches = 50,
    random_mpas = TRUE,
    min_size = 0.1,
    sprinkler = FALSE,
    mpa_habfactor = 1
  )

  raw <- experiment$raw_outcomes %>%
    select(year,
           patch,
           biomass,
           biomass_caught,
           profits,
           effort,
           experiment,
           mpa) %>%
    rename(catch = biomass_caught) %>%
    gather(metric, value, -year,-patch,-experiment,-mpa) %>%
    group_by(year, patch, metric, experiment) %>%
    summarise(value = sum(value),
              mpa = unique(mpa)) %>%
    group_by(year, patch) %>%
    mutate(mpa = unique(mpa[experiment == "with-mpa"])) %>%
    ungroup() %>%
    group_by(metric) %>%
    mutate(value = value / max(value)) %>%
    ungroup() %>%
    spread(experiment, value) %>%
    ungroup() %>%
    mutate(delta = `with-mpa` / `no-mpa`,
           ref = 1)

  doh <- raw %>%
    filter(metric %in% c("biomass","effort")) %>%
    ggplot() +
    geom_col(
      aes(patch, `with-mpa`, fill = mpa),
      width = 1
    ) +
    geom_line(aes(patch, `no-mpa`), size = 2,
              color = "black") +
    facet_wrap(~ metric) +
    coord_polar() +
    gganimate::transition_time(year) +
    gganimate::ease_aes('linear') +
    labs(title = 'Year: {frame_time}', x = "'",  y = "") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "top") +
    scale_fill_npg(labels = c("Fished","MPA"), name = "")


  gganimate::anim_save(animation = doh,filename = file.path(fig_dir,"complex_cs.gif")
  )

  raw_sum <- raw %>%
    group_by(year, metric) %>%
    summarise(no_mpa = sum(`no-mpa`),
              with_mpa = sum(`with-mpa`)) %>%
    ungroup() %>%
    gather(experiment, value, no_mpa:with_mpa) %>%
    filter(year > 20)

  ribbon <- raw_sum %>%
    group_by(year, metric) %>%
    summarise(ymin = min(value),
              ymax = max(value),
              delta = value[experiment == "with_mpa"] - value[experiment == "no_mpa"]) %>%
    ungroup() %>%
    group_by(metric) %>%
    mutate(max_delta = last(delta)) %>%
    ungroup()

  complex_cs_plot <- raw_sum %>%
    ungroup() %>%
    ggplot() +
    geom_line(aes(year, value, color = experiment), size = 1.5) +
    geom_ribbon(data = ribbon, aes(
      year,
      ymin = ymin,
      ymax = ymax,
      fill = max_delta
    ),
    alpha = 0.5,
    show.legend = FALSE) +
    facet_wrap( ~ metric, scales = "free_y",
                labeller = labeller(metric = tools::toTitleCase)) +
    scale_fill_gradient(low = "tomato", high = "steelblue") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.spacing = unit(1,"lines")) +
    labs(x = "Year") +
    scale_color_npg(labels = c("Without MPAs","With MPAs"), name = '')



# save plots --------------------------------------------------------------


rm(diagnostic_plots)

plots <- ls()[str_detect(ls(),"_plot")]



savefoo <- function(fig,
                    device = "pdf",
                    fig_width = 6,
                    fig_height = 5) {
  ggsave(
    filename =  file.path(fig_dir, paste(fig, device, sep = '.')),
    plot  = get(fig),
    width = fig_width,
    height = fig_height
  )

}

walk(plots, safely(savefoo), device = device, fig_height = fig_height,
     fig_width = fig_width)


# pdf(file =  file.path(fig_dir,"did_plot.pdf")
# ,width = 8, height = 5, useDingbats = TRUE)
# print(did_plot)
# dev.off()


