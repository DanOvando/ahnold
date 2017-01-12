#' \code{scrape_enso} extracts El Niño data from
# [here]('https://www.esrl.noaa.gov/psd/enso/mei/table.html)
#'
#' @return a saved csv of ENSO idices
#' @export
#' @import readr
#'
#' @examples
#' scrape_enso(outdir = '../data/')
scrape_enso <- function(outdir = '../data/'){

  library(tidyverse)

  enso <- readr::read_lines('https://www.esrl.noaa.gov/psd/enso/mei/table.html')

  enso <- enso[stringr::str_detect(enso,'\t|(YEAR)')] %>%
    write('enso.txt')

  enso <-  read.csv('enso.txt', sep = '\t', header = F)

  table_names <- enso$V1[1] %>%
    as.character() %>%
    stringr::str_split(stringr::boundary('word'), simplify = T) %>%
    as.character() %>%
    tolower()

  enso <- enso %>%
    slice(-1) %>%
    as_data_frame()

  colnames(enso) <-  table_names

  enso <- enso %>%
    gather('bimonth','enso',-year)

  # enso <- read_table("http://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/Data/nino34.long.anom.data",
  #                    na = c("-99.99", "99.99",'-99'), skip = 1, n_max = lubridate::year(Sys.time()) - 1870 + 1,
  #                    col_names = c("year", 1:12)) %>%
  #   gather(month, enso, -year) %>%
  #   mutate(month = as.double(month))

readr::write_csv(enso, path = paste0(outdir,'enso.csv'))


}