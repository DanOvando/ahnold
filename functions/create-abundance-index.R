create_abundance_index <- function(seen_model, seeing_model, seeing_aug, seen_aug){
  # find seen omitted terms
  omitted_levels <- data_frame(factor_terms = names(seen_model$xlevels),
                               omitted_term = NA)

  seen_coefs <- broom::tidy(seen_model)
  for (i in seq_along(1:nrow(omitted_levels))){

    possible_levels <- seen_model$xlevels[[i]]

    estimated_levels <- seen_coefs %>%
      filter(str_detect(term, omitted_levels$factor_terms[i])) %>% {
        .$term
      }

    estimated_levels <- str_replace(estimated_levels,omitted_levels$factor_terms[i],'' )

    omitted_level <- possible_levels[!possible_levels %in% estimated_levels]

    omitted_levels$omitted_term[i] <- omitted_level

  }

  omitted_levels <- omitted_levels %>%
    mutate(omitted_name = map2_chr(factor_terms, omitted_term, ~paste(.x,.y, sep = '_')))

  intercept_term <- paste(omitted_levels$omitted_name, collapse = '-') # a whole ton of work to double check what the intercept is

# extract and convert intercept and year terms the old fashioned way

  abundance_index <- seen_coefs %>%
    filter(term == "(Intercept)" | str_detect(term, 'factor_year')) %>%
    mutate(term = str_replace(term, 'factor_year',''),
           abundance_index = exp(estimate + std.error ^ 2 / 2)) %>%
    select(term, abundance_index) %>%
    filter(term != '(Intercept)') %>%
    ungroup() %>%
    mutate(abundance_index = abundance_index / max(abundance_index))

  # in order to not drop a year use model to predict abundance over time for the most frequent factor levels held constant

  # seen_aug <-   seen_model %>%
  #   broom::augment()

  seen_factors <- colnames(seen_aug)[map_lgl(seen_aug, ~class(.x) == 'factor' |class(.x) == 'character')]

top_factors <- seen_aug %>%
    group_by_at(vars(one_of(seen_factors[seen_factors != 'factor_year']))) %>%
    count() %>%
    arrange(desc(n)) %>%
  ungroup() %>%
  slice(1) %>%
  select(-n)

seen_series <- seen_aug %>%
  purrrlyr::dmap_if(is.character, as.factor)

factor_names <- colnames(top_factors)

for (i in factor_names){ # filter down to the most common thing

  where_seen <- (seen_series[,i] == (top_factors[,i][[1]]))

  seen_series <- seen_series[where_seen,]

}

seen_series <- seen_series %>%
  slice(1)

num_original_years <- seen_model$xlevels$factor_year

seen_series <- seen_series[rep(1, length(num_original_years)),] #replicate the number of years

seen_series$factor_year <- seen_model$xlevels$factor_year %>% as.factor()

smearing_term <- mean(exp(seen_aug $.resid))

seen_series <- seen_series %>%
  modelr::add_predictions(seen_model) %>%
  mutate(linear_predictions = exp(pred) * smearing_term) %>%
  ungroup() %>%
  mutate(linear_predictions = linear_predictions / max(linear_predictions))

# ggplot() +
#   geom_point(data = abundance_index, aes(term, abundance_index)) +
#  geom_point(data = seen_series, aes(factor_year, linear_predictions), color = 'red')


  # create reference frame for probabilities and calculate abundance index

prob_seen <- predict(seeing_model,newdata = seen_series, type = 'response') # calculate probabilty of seeing

  seen_series <- seen_series %>%
    mutate(prob_seen = prob_seen) %>%
    mutate(abundance_index = linear_predictions * prob_seen,
           abundance_index = abundance_index / max(abundance_index))

  # seen_series %>%
  #   ggplot(aes(factor_year, abundance_index)) +
  #   geom_point()
    return(seen_series)

  # predict observation probabilities for reference frame

  # create composite index

}