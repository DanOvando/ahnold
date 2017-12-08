max_generations <- 5

max_year <- 2013

did_data <- abundance_indices %>%
  select(
    classcode,
    population_filtering,
    population_structure,
    data_source,
    contains('_index')
  ) %>%
  select(-abundance_index) %>%
  gather('abundance_source', 'abundance_index', contains('_index')) %>%
  mutate(abundance_index_passed = !map_lgl(abundance_index, is.null)) %>%
  filter(abundance_index_passed == T) %>%
  unnest() %>%
  filter(year <= max_year) %>%
  left_join(life_history_data, by = 'classcode') %>%
  left_join(enso, by = 'year') %>%
  left_join(pdo, by = 'year') %>%
  left_join(ci_catches, by = c('classcode', 'year')) %>%
  left_join(species_distributions, by = 'classcode') %>%
  mutate(catch = ifelse(is.na(catch), 0, catch)) %>%
  mutate(
    post_mpa = as.numeric(year >= 2003),
    factor_year = as.factor(year),
    generations_protected = pmin(round(pmax(
      0, (year - year_mpa - 1) / tm
    )), max_generations)
  ) %>%
  mutate(abundance_index = abundance_index + 1e-6) %>%
  mutate(
    log_abundance_index = log(abundance_index),
    generations_protected = as.factor(generations_protected),
    targeted = targeted / max(targeted)
  ) %>%
  filter(population_filtering == 'all',
         data_source == 'length_to_density',
         population_structure == 'one-pop',
         abundance_source == 'glm_abundance_index')


annual_conditions_foo <- function(population_structure, data, abundance_source, annual_conditions,annual_regional_conditions){

  if (population_structure == 'regional-pops' &  abundance_source != 'vast_abundance_index'){
    data <- data %>%
      left_join(annual_regional_conditions, by = c('year','population_level' = 'region', 'classcode'))
  } else{

    data <- data %>%
      left_join(annual_conditions, by = c('year','classcode'))
  }

  return(data)
}

did_data <- did_data %>%
  nest(-population_structure, -abundance_source) %>%
  mutate(data = pmap(list(population_structure = population_structure,
                          data = data,
                          abundance_source = abundance_source), annual_conditions_foo, annual_conditions,
                     annual_regional_conditions)) %>%
  unnest() %>%
  mutate(temp_deviation = abs(mean_annual_temp - temperature))


compare_annual_abundance <- function(data) {
  out <-  data %>%
    group_by(abundance_source) %>%
    mutate(abundance_index = (abundance_index - mean(abundance_index)) / (2 * sd(abundance_index))) %>%
    ggplot(aes(year, abundance_index, color = abundance_source)) +
    geom_line() +
    geom_point()
}



annual_abundance_trends_foo <-
  function(classcode,
           data_source,
           annual_abundance_trends,
           run_dir) {
    ggsave(
      file = paste0(
        run_dir,
        '/',
        classcode,
        '-',
        data_source,
        '-annual-abundance-trends.pdf'
      ),
      annual_abundance_trends,
      height = 8,
      width = 8
    )
  }

trend_data <- did_data %>%
  filter(population_structure == 'one-pop',
         population_filtering ==  'all') %>%
  nest(-classcode, -data_source) %>%
  mutate(annual_abundance_trends = map(data, compare_annual_abundance))

# pwalk(
#   list(
#     classcode = trend_data$classcode,
#     data_source = trend_data$data_source,
#     annual_abundance_trends = trend_data$annual_abundance_trends
#   ),
#   annual_abundance_trends_foo,
#   run_dir = run_dir
# )

did_models <-
  cross_df(
    list(
      did_data = list(did_data),
      timing = c('years', 'generations'),
      complexity = c('bare_bones', 'kitchen_sink'),
      dirty_dishes = 'loo + (mean_enso + mean_annual_kelp + temp_deviation +mean_pdo + lag1_pdo + lag2_pdo + lag3_pdo + lag4_pdo |classcode)'
    )
  )


did_models <- did_models %>%
  slice(3) %>%
  mutate(did_model = pmap(list(did_data = did_data,
                               timing = timing,
                               complexity = complexity,
                               dirty_dishes = dirty_dishes), fit_did,
                          cores = 1,
                          chains = 1)) %>%
  select(-did_data) %>%
  unnest()


did_models <- did_models %>%
  mutate(did_plot = map(did_model, plot_did),
         did_diagnostics = map(did_model, diagnostic_plots))


did_plot_foo <-
  function(data_source,
           abundance_source,
           population_filtering,
           population_structure,
           timing,
           complexity,
           did_plot,
           run_dir) {
    ggsave(
      file = paste0(
        run_dir,
        '/',
        data_source,
        '-',
        abundance_source,
        '-',
        population_filtering,
        '-',
        population_structure,
        '-',
        timing,
        '-',
        complexity,
        '-did.pdf'
      ),
      did_plot,
      height = 8,
      width = 8
    )
  }

pwalk(
  list(
    data_source = did_models$data_source,
    abundance_source =  did_models$abundance_source,
    population_filtering = did_models$population_filtering,
    population_structure = did_models$population_structure,
    timing = did_models$timing,
    complexity = did_models$complexity,
    did_plot = did_models$did_plot
  ),
  did_plot_foo,
  run_dir = run_dir
)

check_ahnold(
  length_to_density_data = length_to_density_data,
  abundance_indices = abundance_indices,
  did_models = did_models
)
