library(dplyr)
library(ggplot2)
library(tidyr)

rawdat <- read.csv('UCSB_FISH raw thru 2013.csv', stringsAsFactors = F)

life.history <- read.csv('VRG Fish Life History in MPA_04_08_11_12 11-Mar-2014.csv', stringsAsFactors = F) %>%
  rename(classcode = pisco_classcode)

species.dat <- read.csv('master_spp_table.csv', stringsAsFactors = F)

site.dat <- read.csv('Final_Site_Table_UCSB.csv', stringsAsFactors = F) %>%
  rename(site = SITE)

mlpa.dat <- left_join(rawdat,life.history, by = 'classcode') %>%
  left_join(site.dat, by = 'site') %>%
  left_join(species.dat, by = 'classcode' )


head(mlpa.dat)

mlpa.summary <- mlpa.dat %>%
  group_by(year,MPA_STATUS,Targeted) %>%
  summarise(total.count = sum(count), mean.total.length = mean(fish_tl, na.rm =T))


ggplot(subset(mlpa.summary,MPA_STATUS == 'SMCA') , aes(year,total.count, shape = Targeted, color = Targeted)) + geom_point()

ggplot(mlpa.dat, aes(vis,log(count))) + geom_point() + geom_smooth()

