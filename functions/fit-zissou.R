fit_zissou <- function(data,
                       non_nested_variables =  c(
                         'targeted',
                         'factor_year',
                         'level',
                         'factor_month',
                         'cumulative_n_obs',
                         'temp_deviation',
                         'surge',
                         'mean_canopy',
                         'mean_depth'
                       ),
                       script_name,
                       seed = 42,
                       run_dir,
                       use_tmb = T,
                       fixed_did = T,
                       include_intercept = T,
                       center_and_scale = T,
                       fixed_regions = F) {


  numeric_species_key <-
    data_frame(classcode = unique(data$classcode)) %>%
    arrange(classcode) %>%
    mutate(numeric_classcode = 1:nrow(.))


  seen_has_important <- data %>%
    filter(any_seen == T) %>%
    select(non_nested_variables) %>%
    mutate(index = 1:nrow(.)) %>%
    na.omit()

  seeing_has_important <- data %>%
    select(non_nested_variables) %>%
    mutate(index = 1:nrow(.)) %>%
    na.omit()

  seen_data <- data %>%
    filter(any_seen == T) %>%
    left_join(numeric_species_key, by = "classcode") %>%
    slice(seen_has_important$index)

  log_density <- seen_data$log_density

  seen_data <- seen_data %>%
    select(-log_density)

  seeing_data <- data %>%
    left_join(numeric_species_key, by = "classcode") %>%
    slice(seeing_has_important$index)

  any_seen <- seeing_data$any_seen

  seeing_data <- seeing_data %>%
    select(-any_seen)


  # prepare data for c++ ----------------------------------------------------


  seen_cdata <-
    make_c_worthy(seen_data,
                  non_nested_vars = non_nested_variables,
                  fixed_regions = fixed_regions,
                  include_intercept = include_intercept,
                  fixed_did = fixed_did,
                  numeric_species_key = numeric_species_key)

  seeing_cdata <-
    make_c_worthy(seeing_data,
                  non_nested_vars = non_nested_variables,
                  fixed_regions = fixed_regions,
                  include_intercept = include_intercept,
                  fixed_did = fixed_did,
                  numeric_species_key = numeric_species_key)

  # prepare standardized matrices -------------------------------------------

  standard_year_species <- expand.grid(year = unique(seeing_data$factor_year), classcode = unique(seeing_data$classcode), stringsAsFactors = F) %>%
    as_data_frame() %>%
    mutate(marker = 1,
           classcode_year = paste(classcode, year, sep = "-")) %>%
    spread(classcode_year, marker, fill = 0) %>%
    arrange(classcode, year) %>%
    left_join(numeric_species_key, by = "classcode")

  standard_species <- standard_year_species$numeric_classcode

  standard_year_species <-  standard_year_species %>%
    select(-year,-classcode,-numeric_classcode)

  standard_non_nested <- seeing_cdata$x_non_nested %>%
    gather(variable, value) %>%
    group_by(variable) %>%
    summarise(mean_value = mean(value, na.rm = T)) %>%
    spread(variable, mean_value)

  standard_non_nested <-
    standard_non_nested[rep(1, nrow(standard_year_species)), ]

  standard_non_nested <-
    standard_non_nested[colnames(seeing_cdata$x_non_nested)]

  standard_region_cluster <- seeing_cdata$x_region_cluster %>%
    gather(variable, value) %>%
    group_by(variable) %>%
    summarise(mean_value = mean(value, na.rm = T)) %>%
    spread(variable, mean_value)

  standard_region_cluster <-
    standard_region_cluster[rep(1, nrow(standard_year_species)), ]

  standard_region_cluster <-
    standard_region_cluster[colnames(seeing_cdata$x_region_cluster)]


  # fit model ---------------------------------------------------------------

  seen_species_index <- seen_data$numeric_classcode

  n_species <- length(unique(seen_species_index))

  zissou_data <- list(
    x_seen_non_nested = seen_cdata$x_non_nested,
    x_seen_year_species = seen_cdata$x_year_species,
    x_seen_region_cluster = seen_cdata$x_region_cluster,
    x_seeing_non_nested = seeing_cdata$x_non_nested,
    x_seeing_year_species = seeing_cdata$x_year_species,
    x_seeing_region_cluster = seeing_cdata$x_region_cluster,
    region_cluster_index = seeing_cdata$region_cluster_index,
    log_density = log_density,
    any_seen = any_seen,
    standard_non_nested = standard_non_nested,
    standard_year_species = standard_year_species,
    standard_region_cluster = standard_region_cluster,
    seen_species_index = seen_species_index,
    year_species_index = seen_cdata$year_species_index,
    seen_weights = rep(1, nrow(seen_cdata$x_non_nested)),
    seeing_weights = rep(1, nrow(seeing_cdata$x_non_nested))
  )

  zissou_data <- map_if(zissou_data, is.data.frame, ~ as.matrix(.x))

  any_na <- map_lgl(zissou_data, ~ any(is.na(.x))) %>% any()

  if (any_na) {
    stop("NAs in zissou_data")
  }

  zissou_params <- list(
    seen_non_nested_betas = rep(0, ncol(zissou_data$x_seen_non_nested)),
    seen_year_species_betas = rep(0, ncol(zissou_data$x_seen_year_species)),
    seen_year_species_sigmas = rep(log(1), n_species),
    seen_region_cluster_betas = rep(0, ncol(zissou_data$x_seen_region_cluster)),
    seen_region_cluster_sigmas = rep(log(1), n_distinct(zissou_data$region_cluster_index)),
    seen_density_species_sigma = rep(log(1), n_species),
    seeing_non_nested_betas = rep(0, ncol(zissou_data$x_seeing_non_nested)),
    seeing_year_species_betas = rep(0, ncol(zissou_data$x_seen_year_species)),
    seeing_year_species_sigmas = rep(log(1), n_species),
    seeing_region_cluster_betas = rep(0, ncol(zissou_data$x_seeing_region_cluster)),
    seeing_region_cluster_sigmas = rep(log(1), n_distinct(zissou_data$region_cluster_index))
  )

  any_na <- map_lgl(zissou_params, ~ any(is.na(.x))) %>% any()

  if (any_na) {
    stop("NAs in ahnold_params")
  }
browser()
  if (use_tmb == T) {

    compile(here::here("src", paste0(script_name, ".cpp")), "-O0") # what is the -O0?

    dyn.load(dynlib(here::here("src", script_name)))

    if (fixed_regions == T){

      randos <- c(
        "seen_year_species_betas",
        "seeing_year_species_betas"
      )
    } else {
      randos = c(
        "seen_year_species_betas",
        "seeing_year_species_betas",
        "seen_region_cluster_betas",
        "seeing_region_cluster_betas"
      )    }

    ahnold_model <-
      MakeADFun(
        zissou_data,
        zissou_params,
        DLL = script_name,
        random = randos
      )

    # save(file = here::here(run_dir, "ahnold-onestage-tmb-model.Rdata"),
    #      ahnold_model)

browser()
    a <- Sys.time()
    set.seed(seed)
    ahnold_fit <-
      nlminb(
        ahnold_model$par,
        ahnold_model$fn,
        ahnold_model$gr,
        control = list(iter.max = 4000, eval.max = 5000)
      )
    Sys.time() - a

    save(file = here::here(run_dir, "ahnold-onestage-tmb-fit.Rdata"),
         ahnold_fit)

    ahnold_report <- ahnold_model$report()

    ahnold_sd_report <- sdreport(ahnold_model)


    save(file = here::here(run_dir, "ahnold-tmb-onestage-sdreport.Rdata"),
         ahnold_sd_report)

    save(file = here::here(run_dir, "ahnold-tmb-onestage-report.Rdata"),
         ahnold_report)


  }


  ahnold_estimates <-
    summary(ahnold_sd_report) %>%
    as.data.frame() %>%
    mutate(variable = rownames(.)) %>%
    set_names(tolower) %>%
    rename(std_error = `std. error`) %>%
    mutate(lower = estimate - 1.96 * std_error,
           upper = estimate + 1.96 * std_error)

  ahnold_fits <- summary(ahnold_report) %>%
    as.data.frame() %>%
    mutate(variable = rownames(.)) %>%
    set_names(tolower)

  abundance_trends <- data_frame(log_abundance_hat = ahnold_report$log_abundance_hat,
                                 classcode = standard_species) %>%
    mutate(abundance_hat = exp(log_abundance_hat)) %>%
    group_by(classcode) %>%
    mutate(year = 1999 + 1:length(abundance_hat)) %>%
    mutate(scaled_abundance_hat = (abundance_hat - mean(abundance_hat)) / sd(abundance_hat)) %>%
    ungroup()

  raw <- pisco %>%
    group_by(year, classcode) %>%
    summarise(md = mean(density_g_m2)) %>%
    group_by(classcode) %>%
    mutate(smd = (md - mean(md))/(sd(md))) %>%
    ungroup() %>%
    left_join(numeric_species_key, by = "classcode") %>%
    select(-classcode) %>%
    rename(classcode = numeric_classcode)


  abundance_trends %>%
    ggplot(aes(year, scaled_abundance_hat, color = factor(classcode))) +
    geom_line(show.legend = F) +
    geom_vline(aes(xintercept = 2003), linetype = 2, color = "red") +
    geom_point(data = raw,aes(year, smd, color = factor(classcode)), show.legend = F)+
    facet_wrap(~classcode) +
    theme_minimal()





  return(
    list(
      seen_cdata = seen_cdata,
      seeing_cdata = seeing_cdata,
      ahnold_fit = ahnold_fit,
      ahmold_model = ahnold_model,
      ahnold_report = ahnold_report,
      ahnold_sd_report = ahnold_sd_report,
      ahnold_estimates = ahnold_estimates
    )
  )
}