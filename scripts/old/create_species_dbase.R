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
library(spasm)
library(rfishbase)

demons::load_functions()


fish <- create_fish(adult_movement = 15,
                    density_dependence_form = 2)

fleet <-
  create_fleet(
    eq_f = .2,
    length_50_sel = 1,
    length_95_sel = 2,
    fish = fish,
    mpa_reaction = 'concentrate'
  )


with_mpa <-
  spasm::sim_fishery(
    fish = fish,
    fleet = fleet,
    manager  = create_manager(mpa_size = 0.5, year_mpa = 3),
    num_patches = 20,
    burn_year = 25,
    sim_years = 60
  ) %>%
  mutate(run = 'with_mpa')

with_mpa %>%
  group_by(year) %>%
  summarise(prop_mpa = mean(mpa))

life_history_data <-
  read_csv('data/VRG Fish Life History in MPA_04_08_11_12 11-Mar-2014.csv') %>%
  rename(classcode = pisco_classcode) %>%
  mutate(classcode = tolower(classcode)) %>%
  rename(description_2 = Description) %>%
  magrittr::set_colnames(., tolower(colnames(.)))

# Goals
# get data you can get
#

arg <- life_history_data %>%
  mutate(fook = str_detect(spelling.old.name, '( sp)') == F |
           str_detect(taxa, '( sp.)')) %>%
  select(spelling.old.name, commonname, taxa, fook)


fish_data <- life_history_data %>%
  mutate(has_both = map_dbl(taxa, ~ (
    str_split(.x, pattern = ' ') %>% unlist() %>% length()
  ))) %>%
  filter(
    str_detect(commonname, 'YOY') == F,
    has_both == 2,
    str_detect(commonname, ' sp.') == F,
    str_detect(taxa, ' sp.') == F,
    str_detect(taxa, '(unidentified)') == F,
    is.na(family) == F,
    is.na(classcode) == F
  ) %>%
  select(
    classcode,
    commonname,
    taxa,
    wl_a,
    wl_b,
    wl_w_units,
    wl_l_units,
    vbgf.linf,
    vbgf.k,
    vbgf.t0,
    vbgf_input_length,
    size_mature_cm
  ) %>%
  rename(scientific_name = taxa)

for (i in 1:dim(fish_data)[1]){

  if (is.na(fish_data$vbgf.linf[i]) == F &
       is.na(fish_data$wl_a[i]) == F &
       is.na(fish_data$wl_b[i]) == F &
      fish_data$wl_l_units[i] == 'mm') {


        new_wl <- refit_lw(linf = fish_data$vbgf.linf[i],
                       a = fish_data$wl_a[i],
                       b = fish_data$wl_b[i], linf_factor = 1/10)

    fish_data$wl_a[i] <- new_wl$a

    fish_data$wl_b[i] <- new_wl$b

  }

}

fish_data <- fish_data %>%
   mutate(wl_a = wl_a / 1000,
          infweight = wl_a * vbgf.linf ^ wl_b)

fish_data %>%
  ggplot(aes(vbgf_input_length)) +
  geom_bar()

fish_data %>%
  ggplot(aes(wl_l_units)) +
  geom_bar()


simufish <- fish_data %>%
  mutate(fish = pmap(
    list(
      scientific_name = scientific_name,
      linf = vbgf.linf,
      vbk = vbgf.k,
      t0 = vbgf.t0,
      length_mature = size_mature_cm,
      weight_a = wl_a,
      weight_b = wl_b
    ),
    create_fish
  )) %>%
  select(fish) %>%
  mutate(has_enough = map(fish, ~ is.na(.$linf) == F &
                            any(is.na(.$maturity_at_age) == F))) %>%
  filter(has_enough == T)

foo <- function(fish,fleet){

  length(fish) + length(fleet)
}

simulations <- simufish %>%
  mutate(fleet = pmap(list(fish = fish, eq_f = rlnorm(1, log(.1),.5)),
                      create_fleet, length_50_sel = 1,
                      length_95_sel = 2,
                      mpa_reaction = 'concentrate')) %>%
mutate(simulation = pmap(list(fish = fish,
                         fleet =fleet),
                  sim_fishery,
                  manager = create_manager(
                    mpa_size = .25,
                    year_mpa = 35),
                  num_patches = 20,
                  burn_year = 25,
                 sim_years = 100))



plot_simulation <- function(sim){

sim %>%
    group_by(year) %>%
    summarise(ssb = sum(ssb), prop_mpa = mean(mpa)) %>%
    ggplot(aes(year,ssb, fill = prop_mpa)) +
    geom_point(shape = 21) +
    scale_fill_viridis(labels = percent, name = '% MPA') +
    xlab('Year') +
    ylab('SSB')

}

simulations <- simulations %>%
  mutate(sim_plot = map_plot(simulation,plot_simulation)) %>%
  mutate(species = map_chr(fish,'scientific_name'))

trelliscope(simulations %>% select(species, sim_plot), panel_col = 'sim_plot', name = 'testing')

hmm <- simulations %>%
  filter(species == 'Aulorhynchus flavidus')

create_fish(scientific_name = 'Aulorhynchus flavidus',t0 = -2.31)

