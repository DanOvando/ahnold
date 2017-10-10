fit_fish <-
  function(data,
           pop_structure,
           pop_filter,
           ind_covars,
           dep_var,
           family,
           fit_seen = T,
           consistent_sites) {
    if (pop_structure == 'one-pop') {
      pop_term <- 'factor_year'

    }
    if (pop_structure == 'regional-pops') {
      pop_term <- 'factor_year*region'

    }
    if (pop_structure == 'mpa-pops') {
      pop_term <- 'factor_year*eventual_mpa'

    }

    ind_vars <- paste(c(pop_term, ind_covars), collapse = '+')

    if (pop_filter == 'mpa-only') {
      data <- data %>%
        filter(eventual_mpa == T)

    }
    if (pop_filter == 'consistent-sites'){
      data <- data %>%
        filter(site %in% consistent_sites$site)
    }

      if (fit_seen == T) {
        data <- data %>%
          filter(any_seen == fit_seen & is.na(any_seen) == F)
      } else{
        data <- data %>%
          filter(is.na(any_seen) == F)
      }

    reg_fmla <- paste(dep_var, ind_vars, sep = '~') %>% as.formula()
    model <- glm(reg_fmla, data  = data, family = family)
  }