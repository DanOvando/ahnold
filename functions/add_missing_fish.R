#' add missing fish to data
#'
#' @param this_site
#' @param this_side
#' @param this_year
#' @param this_transect
#' @param observations
#' @param species_sightings
#'
#' @return
#' @export
#'
add_missing_fish <-
  function(this_site,
           this_side,
           this_year,
           this_transect,
           observations,
           species_sightings) {

    sampling_event <- observations %>%
      ungroup() %>%
      filter(site == this_site,
             side == this_side,
             year == this_year,
             transect == this_transect)

    species_seen <- unique(sampling_event$classcode)

    species_possible <- species_sightings %>%
      filter(site == this_site) %>%
      select(species_seen) %>%
      unlist() %>%
      as.character()

    species_missing <-
      species_possible[!species_possible %in% species_seen]

    if (length(species_missing) > 0) {
      blank <- sampling_event[rep(1, length(species_missing)),] %>%
        mutate(classcode = species_missing) %>%
        gather('metric', 'value', dplyr::contains('biomass_g')) %>%
        mutate(value = 0) %>%
        spread(metric, value)

      missing_fish <- blank[, colnames(sampling_event)]


      sampling_event <- sampling_event %>%
        bind_rows(missing_fish)
    }
    return(sampling_event)

  }