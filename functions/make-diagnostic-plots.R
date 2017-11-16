diagnostic_plots <- function(model){


  augmod <- broom::augment(model)
  normal_qq_plot <- augmod %>%
    ggplot(aes(sample = .resid)) +
    stat_qq() +
    stat_qq_line(color = 'red')

  hist_resid_plot <- augmod %>%
    ggplot(aes(.resid)) +
    geom_histogram() +
    geom_vline(aes(xintercept = 0), color = 'blue', linetype = 2) +
    geom_vline(aes(xintercept = mean(.resid)), color = 'red', linetype = 3)


  resid_v_fitted_plot <- augmod %>%
    ggplot(aes(.fitted, .resid)) +
    geom_point() +
    geom_hline(aes(yintercept = 0), color = 'red')

  # qq_plot <- augmod %>%
  #   ggplot(aes(sample = .resid)) +
  #   stat_qq(distribution = qpois,
  #           dparams = list(lambda = mean(model$y))) +
  #   stat_qq_line(color = 'red',
  #                distribution = qpois,
  #                dparams = list(lambda = mean(model$y)))


  if (class(model)[[1]] == 'lmerMod'){

    ind_vars <- model@frame %>% colnames()

  }else{

  ind_vars <- (model$terms %>% as.character())[3]

  ind_vars <- str_split(ind_vars, pattern = '\\*|\\+', simplify = T)

  ind_vars <-  map_chr(ind_vars, ~str_replace_all(.x, ' ',''))
}
  numeric_vars <- colnames(augmod)[map_lgl(augmod, is.numeric)]

  factor_vars <-  colnames(augmod)[map_lgl(augmod, ~!is.numeric(.))]


  numeric_coef_v_resid <- augmod %>%
    select(ind_vars[ind_vars %in% numeric_vars], .resid) %>%
    gather(variable, value, -.resid)

  factor_coef_v_resid <- augmod %>%
    select(ind_vars[ind_vars %in% factor_vars], .resid) %>%
    gather(variable, value, -.resid)


  numeric_coef_v_resid_plot <- numeric_coef_v_resid %>%
    ggplot(aes(value, .resid)) +
    geom_point() +
    facet_wrap(~variable, scales = 'free_x')


  factor_coef_v_resid_plot <- factor_coef_v_resid %>%
    ggplot(aes(value, .resid)) +
    geom_boxplot() +
    facet_wrap(~variable, scales = 'free_y') +
    coord_flip()



out <- list(factor_coef_v_resid_plot = factor_coef_v_resid_plot,
            numeric_coef_v_resid_plot = numeric_coef_v_resid_plot,
            resid_v_fitted_plot = resid_v_fitted_plot,
            hist_resid_plot = hist_resid_plot,
            normal_qq_plot = normal_qq_plot)

}