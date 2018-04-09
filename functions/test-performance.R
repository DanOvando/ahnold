test_performance <-
  function(fishes,
           year_mpa,
           min_year = 75,
           max_year = 100,
           time_step = 1) {

    simple_data <- fishes %>%
      select(loo, k, lm, m, targeted, classcode, commonname, pisco_samples) %>%
      unnest() %>%
      filter(year > min_year, year < max_year) %>%
      select(-pop, -sampled_lengths, -diver_stats) %>%
      mutate(year = year * time_step) %>%
      mutate(subyear = year - floor(year),
             year = floor(year)) %>%
      mutate(
        factor_year = as.factor(year),
        log_density = log(density),
        logical_targeted = targeted > 0,
        any_seen = density > 0,
        region = patches,
        post_mpa = year >= year_mpa
      )

    years_protected <- unique(simple_data$year) - year_mpa

    bins <- c(seq(min(years_protected), -1, by = 5), seq(0, max(years_protected), by = 5), max(years_protected))

    year_block <- data_frame(year = unique(simple_data$year)) %>%
      mutate(years_protected = year - year_mpa) %>%
      mutate(protected_block = cut(
        years_protected,
        breaks = bins,
        include.lowest = T
      )) %>%
      select(-years_protected)

    simple_data <- simple_data %>%
      left_join(year_block, by = "year")

    true_effect <- fishes %>%
      select(classcode, targeted, mpa_effect) %>%
      unnest() %>%
      filter(year > min_year, targeted == 1) %>%
      mutate(post_mpa = year > year_mpa) %>%
      mutate(year = year * time_step) %>%
      mutate(subyear = year - floor(year),
             year = floor(year))

    bare_bones_model <-
      lm(log_density ~ targeted + protected_block + targeted:protected_block,
         data = simple_data)

    if (n_distinct(simple_data$diver) > 1) {
      # mixed_effect_model <- lme4::lmer(log_density ~ targeted + factor_year + targeted:factor_year + diver + (factor_year |classcode), data = simple_data)
      # mixed_effect_model <-
      #   rstanarm::stan_glm(
      #     log_density ~  targeted + protected_block + targeted:protected_block + loo,
      #     data = simple_data,
      #     iter = 4000,
      #     warmup = 2000,
      #     cores = 1,
      #     refresh = 25,
      #     QR = TRUE,
      #     control = list(max_treedepth = 10)
      #   )

      mixed_effect_model <-
        rstanarm::stan_glmer(
          log_density ~  (1|classcode) + (protected_block - 1 | targeted:protected_block),
          data = simple_data,
          iter = 4000,
          warmup = 2000,
          cores = 1,
          refresh = 25,
          control = list(max_treedepth = 10)
        )
browser()
    } else
    {
      # mixed_effect_model <-
      #   lme4::lmer(
      #     log_density ~ targeted + factor_year + targeted:factor_year + (factor_year |
      #                                                                      classcode),
      #     data = simple_data
      #   )
      mixed_effect_model <-
        rstanarm::stan_glm(
          log_density ~ targeted + protected_block + targeted:protected_block + loo,
          data = simple_data %>% select(log_density, targeted, protected_block,factor_year, classcode,loo),
          iter = 4000,
          warmup = 2000,
          chains = 4,
          cores = 4,
          refresh = 1,
          QR = TRUE,
          control = list(max_treedepth = 8)
        )

    }

    pre_post_model <-
      lm(log_density ~ targeted + post_mpa + targeted:post_mpa, data = simple_data)

    get_range <- function(bin){

      bin_range <- str_extract(bin, pattern = '(?<=\\().*(?=])')

      mean(str_split(bin_range, ',', simplify = T) %>% as.numeric())

    }



    bare_bones_did <- broom::tidy(bare_bones_model) %>%
      filter(str_detect(term, 'targeted:')) %>%
      mutate(year = map_dbl(term, get_range)) %>%
      ggplot() +
      geom_pointrange(
        aes(
          x = year,
          y = estimate,
          ymin = estimate - 1.96 * std.error,
          ymax = estimate + 1.96 * std.error
        )
      ) +
      geom_line(data = true_effect, aes(year - year_mpa, mpa_effect, color = classcode), show.legend = F) +
      labs(title = 'bare bones')

    mixed_effect_did <- broom::tidy(mixed_effect_model) %>%
      filter(str_detect(term, 'targeted:')) %>%
      mutate(year = map_dbl(term, get_range)) %>%
      ggplot() +
      geom_pointrange(
        aes(
          x = year,
          y = estimate,
          ymin = estimate - 1.96 * std.error,
          ymax = estimate + 1.96 * std.error
        )
      ) +
      geom_line(data = true_effect, aes(year - year_mpa, mpa_effect, color = classcode), show.legend = F) +
      labs(title = 'mixed effects')

    check_block_did <- broom::tidy(pre_post_model) %>%
      filter(str_detect(term, 'targeted:')) %>%
      mutate(post_mpa = TRUE) %>%
      ggplot() +
      geom_boxplot(data = true_effect,
                   aes(post_mpa, mpa_effect),
                   color = 'red') +
      geom_pointrange(
        aes(
          x = post_mpa,
          y = estimate,
          ymin = estimate - 1.96 * std.error,
          ymax = estimate + 1.96 * std.error
        )
      )

    out_plot <-
    {bare_bones_did + mixed_effect_did + plot_layout(ncol = 1)} + check_block_did + plot_layout(ncol = 2)

    out <- list(
      bare_bones_did = bare_bones_did,
      mixed_effect_did = mixed_effect_did,
      out_plot = out_plot,
      bare_bones_model = bare_bones_model,
      mixed_effect_model = mixed_effect_model,
      pre_post_model = pre_post_model
    )

    return(out)
  }
