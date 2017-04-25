library(tidyverse)
library(spasm)
library(stringr)
library(rfishbase)
library(modelr)
library(Rcpp)
library(RcppEigen)

# sourceCpp("scripts/eigen_mat_mult.cpp")

fish <- create_fish(adult_movement = 1e-3,
                    larval_movement = 1e-3,
                    density_dependence_form = 1,
                    price = 10,
                    query_fishbase = F,
                    linf = 100,
                    age_mature = 2,
                    vbk = .2,
                    t0 = 0,
                    max_age = 20,
                    weight_a = 1e-3,
                    weight_b = 3)

fleet <-
  create_fleet(
    eq_f = .2,
    length_50_sel = 1,
    length_95_sel = 2,
    fish = fish,
    mpa_reaction = 'concentrate',
    price = 10,
    cost = 0.01,
    beta = 1.3,
    theta = 1e-1,
    q = 1e-3,
    fleet_model = 'constant-effort',
    effort_allocation = 'gravity',
    initial_effort = 1000,
    target_catch = 1
  )

# fish$price <- 10
#
# fleet$cost <- .01
#
# fleet$beta <- 1.3
#
# fleet$theta <- 1e-1
#
# fleet$q <- 1e-3
#
# fleet$fleet_model <- 'open-access'
#
# fleet$effort_allocation <- 'gravity'
#
# fleet$initial_effort <-  1000


without_mpa <-
  spasm::sim_fishery(
    fish = fish,
    fleet = fleet,
    manager  = create_manager(mpa_size = 0, year_mpa = 10),
    num_patches = 100,
    burn_year = 25,
    sim_years = 200
  ) %>%
  mutate(run = 'without_mpa')

without_mpa %>%
  group_by(year) %>%
  summarise(effort = sum(effort),
            profits = sum(profits),
            f = max(f)) %>%
  ggplot(aes(year, effort)) +
  geom_point()

with_mpa <-
  spasm::sim_fishery(
    fish = fish,
    fleet = fleet,
    manager  = create_manager(mpa_size = 0.5, year_mpa = 10),
    num_patches = 10,
    burn_year = 25,
    sim_years = 100
  ) %>%
  mutate(run = 'with_mpa')


with_mpa %>%
  filter(year == max(year)) %>%
  ggplot(aes(age,numbers)) +
  geom_col()

with_mpa %>%
  filter(year == max(year)) %>%
  group_by(patch) %>%
summarise(ssb = sum(ssb), mpa = unique(mpa)) %>%
  ggplot(aes(patch,ssb, fill = mpa)) +
  geom_col()


with_mpa %>%
  filter(year == max(year)) %>%
  group_by(mpa,age) %>%
  summarise(numbers = sum(numbers)) %>%
  ggplot(aes(x = age,y =numbers, fill = mpa)) +
  geom_col(position = 'dodge')





with_mpa %>%
  group_by(year) %>%
  summarise(catch = sum(biomass_caught)) %>%
  ggplot(aes(year,catch)) +
  geom_point()

with_mpa %>%
  group_by(year, eventual_mpa) %>%
  summarise(effort = max(effort),
            f = max(f),
            total_f = sum(f),
            catch = sum(biomass_caught),
            numbers = sum(numbers)) %>%
  ggplot(aes(year, catch, color = eventual_mpa)) +
  geom_line()

with_mpa %>%
  group_by(year,patch) %>%
  summarise(effort = sum(effort),
            pop = sum(ssb),
            is_mpa = unique(mpa)) %>%
  ungroup() %>%
  ggplot(aes(pop, effort, color = is_mpa )) +
  geom_point()

with_mpa %>%
  group_by(year) %>%
  summarise(effort = mean(effort),
            profits = sum(profits),
            f = max(f),
            total_b = sum(ssb),
            catch = sum(biomass_caught)) %>%
  ggplot(aes(year, profits)) +
  geom_point()


with_mpa %>%
  group_by(year) %>%
  summarise(effort = sum(effort),
            profits = sum(profits),
            f = max(f)) %>%
  ggplot(aes(year, effort)) +
  geom_point()


with_mpa %>%
  group_by(year) %>%
  summarise(prop_mpa = mean(mpa)) %>%
  ggplot(aes(year, prop_mpa)) +
  geom_point()


without_mpa %>%
  group_by(year, mpa) %>%
  summarise(numbers = sum(numbers), ssb = sum(ssb)) %>%
  ggplot(aes(year, ssb, color = mpa)) +
  geom_line()

with_mpa %>%
  bind_rows(without_mpa) %>%
  group_by(year, run) %>%
  summarise(ssb = sum(ssb),
            catch = sum(biomass_caught)) %>%
  ggplot(aes(year,ssb, color = run)) +
  geom_line()

