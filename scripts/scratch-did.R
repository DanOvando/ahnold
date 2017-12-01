did_data <- abundance_indices %>%
  select(classcode, population_filtering, population_structure, data_source, contains('_index')) %>%
  select(-abundance_index) %>%
  gather('abundance_source','abundance_index',contains('_index')) %>%
  mutate(abundance_index_passed = !map_lgl(abundance_index, is.null)) %>%
  filter(abundance_index_passed == T) %>%
  unnest() %>%
  left_join(life_history_data, by = 'classcode') %>%
  left_join(enso, by = 'year') %>%
  left_join(pdo, by = 'year') %>%
  left_join(ci_catches, by = c('classcode', 'year')) %>%
  left_join(species_distributions, by = 'classcode') %>%
  mutate(catch = ifelse(is.na(catch), 0, catch)) %>%
  mutate(targeted = as.numeric(targeted == 'Targeted'),
         post_mpa = as.numeric(year >= 2003)) %>%
  mutate(abundance_index = abundance_index + 1e-6)


annual_conditions_foo <- function(population_structure, data, abundance_source, annual_conditions,annual_regional_conditions){

  if (population_structure == 'regional-pops' &  abundance_source != 'vast_abundance_index'){
    data <- data %>%
      left_join(annual_regional_conditions, by = c('year','population_level' = 'region'))
  } else{

    data <- data %>%
      left_join(annual_conditions, by = 'year')
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



compare_annual_abundance <- function(data){


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
    ggsave(file = paste0(run_dir, '/', classcode, '-', data_source, '-annual-abundance-trends.pdf'),
           annual_abundance_trends, height = 8, width = 8)
  }

trend_data <- did_data %>%
  filter(population_structure == 'one-pop',
         population_filtering ==  'all') %>%
  nest(-classcode,-data_source) %>%
  mutate(annual_abundance_trends = map(data, compare_annual_abundance))

pwalk(list(classcode = trend_data$classcode, data_source = trend_data$data_source,
           annual_abundance_trends = trend_data$annual_abundance_trends), annual_abundance_trends_foo,
      run_dir = run_dir)


did_terms <- did_data %>%
  select(year, targeted) %>%
  mutate(index = 1:nrow(.)) %>%
  spread(year, targeted, fill = 0) %>%
  select(-index) %>%
  set_names(., paste0('did_', colnames(.))) %>%
  select(-did_2000)

year_mpa <- 2003

generation_did_terms <- did_data %>%
  select(year, targeted, tm) %>%
  mutate(index = 1:nrow(.),
         years_protected = pmax(0,year - year_mpa),
         generations_protected = round(years_protected / round(tm))) %>%
  select(index,generations_protected, targeted) %>%
  spread(generations_protected, targeted, fill = 0) %>%
  select(-index) %>%
  set_names(., paste0('did_', colnames(.))) %>%
  select(-did_0)

generation_did_terms <- generation_did_terms[, colSums(generation_did_terms) >0]


recruitment_did_terms <- did_data %>%
  select(year, targeted, first_age_seen) %>%
  mutate(index = 1:nrow(.),
         years_protected = pmax(0,year - year_mpa),
         recruits_protected = round(years_protected / round(first_age_seen))) %>%
  select(index,recruits_protected, targeted) %>%
  spread(recruits_protected, targeted, fill = 0) %>%
  select(-index) %>%
  set_names(., paste0('did_', colnames(.))) %>%
  select(-did_0)

# did_terms <- generation_did_terms


did_models <-
  data_frame(
    did_data = list(did_data),
    did_terms = list(did_terms, generation_did_terms, recruitment_did_terms),
    did_term_names = c('years-protected','generations-protected','recruits-protected')
  )


fit_did <- function(did_data, did_terms){

  did_data <- did_data %>%
    bind_cols(did_terms)

  did_data <- did_data %>%
    nest(-population_structure, -data_source, -population_filtering,-abundance_source)

  # visually inspect abundance trends ---------------------------------------


  did_data <- did_data %>%
    mutate(correlation_tests = map(data, test_parallel_trends)) %>%
    mutate(pre_correlation = map_dbl(correlation_tests, ~.x$overall_correlation_test$estimate),
           pre_correlation_signif = map_dbl(correlation_tests, ~.x$overall_correlation_test$p.value))


  # fit DiD estimator on abundance indicies ---------------------------------



  did_reg <-
    paste0('log(abundance_index) ~', paste(
      c('targeted','mean_annual_kelp','temp_deviation','mean_pdo','(1|factor_year:targeted)'), collapse = '+'
    ))

  # did_reg <-
  #   paste0('log(abundance_index) ~', paste(
  #     c(
  #       'targeted',
  #       'mean_annual_kelp',
  #       '(1|factor_year)',
  #       'mean_pdo',
  #       # 'lag1_pdo',
  #       # 'lag2_pdo',
  #       '(1 + temp_deviation|classcode)',
  #       colnames(did_terms)),
  #     collapse = '+'
  #   ))

  did_models <- did_data %>%
    mutate(did_reg = did_reg) %>%
    mutate(did_model = map2(data, did_reg, ~ rstanarm::stan_glmer(.y, data = .x %>% mutate(factor_year = factor(
      year)), chains = 1)))


} # close fit_did


did_models <- did_models %>%
  slice(1) %>%
  mutate(did_model = map2(did_data, did_terms, fit_did)) %>%
  select(did_model,did_term_names) %>%
  unnest()



a <- did_models$did_model[[1]]


did_models$data[[1]] %>%
  group_by(year, targeted) %>%
  summarise(mean_log_abundance = mean(log(abundance_index/1000))) %>%
  ggplot(aes(year, mean_log_abundance, color = targeted == 1)) +
  geom_point() +
  geom_line()

extract_year <- function(term){

  term <- str_split(term,':')

  year <- term[[1]][3] %>% as.numeric()
}

did_terms <- a$stan_summary %>%
  as.data.frame() %>%
  mutate(term = rownames(.)) %>%
  filter(str_detect(term,':1')) %>%
  mutate(year = map_dbl(term, extract_year))

did_terms %>%
  ggplot() +
  geom_hline(aes(yintercept = 0), color = 'red') +
  geom_pointrange(aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`))

unfished_terms <- a$stan_summary %>%
  as.data.frame() %>%
  mutate(term = rownames(.)) %>%
  filter(str_detect(term,':0')) %>%
  mutate(year = map_dbl(term, extract_year))

unfished_terms %>%
  ggplot() +
  geom_hline(aes(yintercept = 0), color = 'red') +
  geom_pointrange(aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`))



# check_model <- did_models %>%
#   filter(data_source == 'length_to_density',
#          population_structure == 'one-pop',
#          abundance_source == 'glm_abundance_index',
#          population_filtering == 'all',
#          did_term_names == 'years-protected')
#
# check_model$did_model[[1]] %>% broom::glance()

did_plot_foo <- function(x) {
  x %>%
    broom::tidy() %>%
    filter(str_detect(term, 'did')) %>%
    mutate(year = str_replace(term, 'did_', '') %>% as.numeric()) %>%
    ggplot() +
    geom_pointrange(aes(
      year,
      estimate,
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    )) +
    # geom_vline(aes(xintercept = 2003),
    #            color = 'red',
    #            linetype = 2) +
    geom_hline(aes(yintercept = 0))
}

did_models <- did_models %>%
  mutate(did_plot = map(did_model, did_plot_foo) ,
         did_diagnostics = map(did_model, diagnostic_plots))




did_plot_foo <-
  function(data_source,
           abundance_source,
           population_filtering,
           population_structure,
           did_term_name,
           did_plot,
           run_dir) {
    ggsave(file = paste0(run_dir, '/', data_source, '-', abundance_source,'-',population_filtering,'-',population_structure,
                         '-',did_term_name, '-did.pdf'),
           did_plot, height = 8, width = 8)
  }

pwalk(
  list(
    data_source = did_models$data_source,
    abundance_source =  did_models$abundance_source,
    population_filtering = did_models$population_filtering,
    population_structure = did_models$population_structure,
    did_term_name = did_models$did_term_names,
    did_plot = did_models$did_plot
  ),
  did_plot_foo,
  run_dir = run_dir
)

check_ahnold(length_to_density_data = length_to_density_data,
             abundance_indices = abundance_indices,
             did_models = did_models)

