library(tidyverse)
library(lubridate)
library(stringr)
library(tabulizer)
library(tabulizerjars)
library(hrbrthemes)
library(rfishbase)
library(scales)
library(forcats)
library(trelliscopejs)

theme_set(theme_ipsum())

file_names <- list.files('data/cdfw-data/')


process_cdfw <- function(file_name) {
  dat <-
    tabulizer::extract_tables(paste0("data/cdfw-data/", file_name))

  year <- str_extract(file_name, pattern = '(?<=landings).*(?=_)')

  year <- as.numeric(paste0('20', year))

  collapse_pages <- function(x) {
    missing_all <-
      map_lgl(x %>% as_tibble, ~ all(str_count(.x) == 0)) == F

    x <- x[, missing_all]

    # if (dim(x)[2] == 15){ #super hacky fix for random blank column in some years
    #
    #   x <- x[,-14]
    #
    # }

    x <- as_tibble(x)

    x <- slice(x,-(1:2))

    species <- x[, 1]
    species <-
      str_replace_all(species$V1, "(\\.)|^[ \t]|[ \t]$", '')

    species_mat <-
      str_split(species, pattern = ',|:', simplify = T)[, c(2, 1)]

    species <-
      unite(species_mat %>% as_data_frame(),
            col = species,
            V1,
            V2,
            sep = ' ') %>%
      mutate(species = str_replace(species, "^[ \t]|[ \t]$", '') %>% tolower())

    x[, 1] <- species$species
    numfoo <- function(z) {
      z <- str_replace_all(z, '\\.|\\,', '')

      if (any(is.na(as.numeric(z)) == F)) {
        z = as.numeric(z)
      } else{
        z = z
      }

    }

    x <- map_df(x, numfoo)


    return(x)
  }

  flat_dat <- map_df(dat, collapse_pages)


  flat_dat <- flat_dat %>%
    set_names(
      c(
        'species',
        'january',
        'february',
        'march',
        'april',
        'may',
        'june',
        'july',
        'august',
        'september',
        'october',
        'november',
        'december',
        'total'
      )
    ) %>%
    select(-total) %>%
    gather('month', 'pounds_caught',-species) %>%
    mutate(year = year) %>%
    filter(str_detect(species, '(total)|(Total)|(ishes)|(aters)') == F)
}

spc <- safely(process_cdfw)

flat_cdfw <- map_df(file_names, process_cdfw)

write_csv(flat_cdfw, path = "processed_data/flat_cdfw.csv")

flat_cdfw <- read_csv( "processed_data/clean_flat_cdfw.csv")

# out <- bind_rows(flat_cdfw)

sci_names <- map_chr((flat_cdfw$common_name %>% unique()), ~common_to_sci(.x)[1])

comm_to_sci <- data_frame(common_name = flat_cdfw$common_name %>% unique(), sci_name = sci_names)


flat_cdfw <- flat_cdfw %>%
  select(-sci_name) %>%
  left_join(comm_to_sci, by = 'common_name') %>%
  mutate(sci_name = tolower(sci_name))

write_csv(flat_cdfw, path = file.path('processed_data','cdfw-catches.csv'))


life_history_data <-
  read_csv('data/VRG Fish Life History in MPA_04_08_11_12 11-Mar-2014.csv') %>%
  rename(classcode = pisco_classcode) %>%
  mutate(classcode = tolower(classcode)) %>%
  # rename(description_2 = Description) %>%
  magrittr::set_colnames(., tolower(colnames(.)))

# seen_species <- seen_reg_data %>%
#   select(commonname, classcode) %>%
#   unique() %>%
#   left_join(life_history_data %>% select(classcode, taxa), by = 'classcode') %>%
#   rename(sci_name = taxa) %>%
#   mutate(sci_name = tolower(sci_name))
#
#
# ci_catches <- flat_cdfw %>%
#   left_join(seen_species, by = 'sci_name') %>%
#   filter(!is.na(classcode))
#
# plotfoo <- function(x){
#
#   x %>%
#     ggplot(aes(year,catch)) +
#     geom_vline(aes(xintercept = 2003), color = 'red',
#                linetype = 2) +
#     geom_line() +
#     geom_point() +
#     scale_y_continuous(labels = comma, name = 'catch (lbs)')
#
# }
#
#
# ci_catches %>%
#   group_by(year) %>%
#   summarise(catch = sum(pounds_caught)) %>%
#   ggplot(aes(year,catch)) +
#   geom_vline(aes(xintercept = 2003), color = 'red',
#              linetype = 2) +
#   geom_line() +
#   geom_point() +
#   scale_y_continuous(labels = comma, name = 'catch (lbs)')
#
#
# ci_catches %>%
#   filter(common_name != 'pacific sardine') %>%
#   group_by(common_name, year) %>%
#   summarise(catch = sum(pounds_caught))  %>%
#   group_by(common_name) %>%
#   mutate(total_catch = sum(catch)) %>%
#   ungroup() %>%
#   mutate(common_name = fct_reorder(common_name,total_catch )) %>%
#   ggplot(aes(year,catch, fill = common_name)) +
#   geom_vline(aes(xintercept = 2003), color = 'red',
#              linetype = 2) +
#   geom_area( show.legend = T,position = 'stack', alpha = .85)
#
# ci_catch_plots <- ci_catches %>%
#   group_by(common_name, year) %>%
#   summarise(catch = sum(pounds_caught)) %>%
#   nest(-common_name) %>%
#   mutate(total_catch = map_dbl(data, ~sum(.x$catch))) %>%
#   arrange(desc(total_catch)) %>%
#   mutate(catch_plot = map_plot(data, plotfoo))
#
# trelliscope(ci_catch_plots, name = 'Channel Islands Catches')
#
#
# a <- seen_reg_data %>%
#   group_by(commonname) %>%
#   summarise(samples = length(biomass),
#             fished = unique(targeted)) %>%
#   arrange(desc(samples))
#
# write_csv(ci_catches, path = file.path('data','cdfw-channel-islands-catches.csv'))



