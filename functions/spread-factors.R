spread_factors <- function(data){


  factors <- colnames(data)[map_lgl(data, ~class(.x) %in% c('character','factor'))]

  for (i in seq_along(factors)) {

  var <- factors[i]

  data <- data %>%
    mutate(dummy = 1, index = 1:nrow(.)) %>%
    spread_(var, 'dummy', fill = 0, sep = 'dummy') %>%
    select(-index,-contains('dummy')[1]) %>%
    set_names(str_replace(colnames(.),'dummy',''))

  }

  return(data)
}