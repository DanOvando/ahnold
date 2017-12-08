plot_enviro_effects <- function(model)



enso_terms <- model$stan_summary %>%
  as.data.frame() %>%
  mutate(term = rownames(.)) %>%
  filter(str_detect(term, 'enso') & (str_detect(term,'Sigma') == F))

enso_terms %>%
  ggplot() +
  geom_hline(aes(yintercept = 0), color = 'red', linetype = 2) +
  geom_pointrange(aes(
    x = term,
    y = mean,
    ymin = `2.5%`,
    ymax = `97.5%`
  )) +
  coord_flip()

pdo_terms <- a$stan_summary %>%
  as.data.frame() %>%
  mutate(term = rownames(.)) %>%
  filter(str_detect(term, 'pdo') & (str_detect(term,'Sigma') == F)) %>%
  mutate(year = map_dbl(term, extract_year))

pdo_terms %>%
  ggplot() +
  geom_hline(aes(yintercept = 0), color = 'red', linetype = 2) +
  geom_pointrange(aes(
    x = term,
    y = mean,
    ymin = `2.5%`,
    ymax = `97.5%`
  )) +
  coord_flip()