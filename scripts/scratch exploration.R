
life_history <- pisco %>%
  select(classcode, targeted, loo) %>%
  unique() %>%
  mutate(numeric_classcode = as.numeric(as.factor(classcode)))


numeric_species_key <-
  data_frame(classcode = unique(pisco$classcode)) %>%
  arrange(classcode) %>%
  mutate(numeric_classcode = 1:nrow(.))

raw <- pisco %>%
  group_by(year, classcode, targeted) %>%
  summarise(md = mean(density_g_m2)) %>%
  group_by(classcode) %>%
  mutate(smd = (md - mean(md))/(sd(md))) %>%
  ungroup() %>%
  left_join(numeric_species_key, by = "classcode")


ahnold_estimates <- tmb_runs$tmb_fit[[1]]$ahnold_estimates

ahnold_estimates %>%
  filter(variable == "mpa_effect") %>%
  mutate(year = 2000:2013) %>%
  ggplot() +
  geom_pointrange(aes(x = year, y = estimate, ymin = lower, ymax = upper)) +
  geom_hline(aes(yintercept = 0),linetype = 2) +
  geom_vline(aes(xintercept = 2003), linetype = 2, color= "red")

ahnold_report <- tmb_runs$tmb_fit[[1]]$ahnold_report


abundance_trends <- data_frame(abundance_hat = ahnold_report$abundance_hat,
                               classcode = rep(1:n_distinct(pisco$classcode), each = length(2000:2013))) %>%
  mutate(log_abundance_hat = log(abundance_hat)) %>%
  group_by(classcode) %>%
  mutate(year = 1999 + 1:length(abundance_hat)) %>%
  mutate(scaled_abundance_hat = (abundance_hat - mean(abundance_hat)) / sd(abundance_hat)) %>%
  ungroup() %>%
  rename(numeric_classcode = classcode) %>%
  left_join(life_history, by = "numeric_classcode" )

abundance_trends %>%
  ggplot(aes(year, scaled_abundance_hat, color = factor(classcode))) +
  geom_line(show.legend = F) +
  geom_vline(aes(xintercept = 2003), linetype = 2, color = "red") +
  geom_point(data = raw,aes(year, smd, color = factor(classcode)), show.legend = F)+
  facet_wrap(~classcode) +
  theme_minimal()


abundance_trends %>%
  ggplot() +
  geom_line(aes(
    year,
    abundance_hat,
    color = factor(targeted),
    group = interaction(targeted, classcode)
  ),
  show.legend = F,
  alpha = 0.5) +
  geom_smooth(aes(year, abundance_hat, color = factor(targeted))) +
  geom_vline(aes(xintercept = 2003), linetype = 2, color = "red")



did_data <- tmb_runs$tmb_fit[[1]]$did_data %>%
  mutate(log_abundance_hat = log(ahnold_report$abundance_hat)) %>%
  mutate(abundance_hat = exp(log_abundance_hat)) %>%
  group_by(classcode) %>%
  mutate(scaled_abundance_hat = (abundance_hat - mean(abundance_hat)) / sd(abundance_hat)) %>%
  ungroup() %>%
  mutate(factor_year = as.factor(year))

a <- lm(abundance_hat ~ targeted + factor_year  + targeted:factor_year, data = did_data %>% filter(year > 1999)) %>%
  broom::tidy() %>%
  filter(str_detect(term,"targeted:"))

test <- stan_glmer(log_abundance_hat ~ (1|classcode) + (factor_year - 1|targeted) - 1, data = did_data, chains = 1)

test <- stan_glmer(scaled_abundance_hat ~ (1|classcode) + (factor_year - 1|targeted) - 1, data = did_data, chains = 1)

of_interest <- test$coefficients[ str_detect(names(test$coefficients),"targeted:")]

mpa_effect <- data_frame(coef = of_interest, variable = names(of_interest)) %>%
  mutate(targeted = map_dbl(variable, ~str_extract(.x, pattern = '(?<=:).*(?=])') %>% as.numeric())) %>%
  mutate(year =  map_dbl(variable, ~str_extract(.x, pattern = "(?<=year).*(?=\\s)") %>% as.numeric()))


mpa_effect %>%
  select(coef, targeted, year) %>%
  spread(targeted, coef) %>%
  mutate(delta = `1` - `0`) %>%
  ggplot(aes(year, delta)) +
  geom_vline(aes(xintercept = 2003), linetype = 2, color = "red") +
  geom_line() +
  geom_point(data = a %>% mutate(year = 2001:2013),aes(year, estimate))



mpa_effect <- data_frame(mpa_effect = ahnold_report$mpa_effect,
                         year = unique(data$year))

mpa_effect %>%
  ggplot(aes(year, mpa_effect)) +
  geom_line() +
  geom_hline(aes(yintercept = 0),linetype = 2) +
  geom_vline(aes(xintercept = 2003), linetype = 2, color= "red")
