
library(tidyverse)

sites <- read.csv("./MLPA Data/Final_Site_Table_UCSB.csv", stringsAsFactors = F) %>%
  select(SITE_SIDE, RESERVE)

master_spp <- read.csv("./MLPA Data/master_spp_table.csv", stringsAsFactors = F) %>%
  group_by(classcode, genus, species) %>%
  summarize(N = n()) %>%
  select(classcode, genus, species, N) %>%
  ungroup()

duplicated_sp <- filter(master_spp, N>1)

master_spp <- select(master_spp, -N)

UCSB_fish <- read.csv("./MLPA Data/UCSB_FISH raw thru 2013.csv", stringsAsFactors = F)

data1 <- UCSB_fish %>%
  left_join(master_spp, by = "classcode") %>%
  mutate(SITE_SIDE = paste(site, side, sep = "_")) %>%
  left_join(sites, by = "SITE_SIDE") %>%
  mutate(sciname = paste(genus, species)) %>%
  filter(!sciname == " ") %>%
  group_by(day, month, year, site, side, RESERVE, zone, level, transect) %>%
  summarize(N = sum(count)) %>%
  mutate(key = paste(day, month, year, site, side, RESERVE, zone, level, transect)) %>%
  ungroup() %>%
  select(key, N)

data2 <- UCSB_fish %>%
  left_join(master_spp, by = "classcode") %>%
  mutate(SITE_SIDE = paste(site, side, sep = "_")) %>%
  left_join(sites, by = "SITE_SIDE") %>%
  mutate(sample = paste(day, month, year, site, side, RESERVE, zone, level, transect, sep = "-"),
         key = paste(day, month, year, site, side, RESERVE, zone, level, transect),
         sciname = paste(genus, species)) %>%
  filter(!sciname == " ") %>%
  group_by(day, month, year, site, side, RESERVE, zone, level, transect, sample, key, sciname) %>%
  summarize(count = sum(count)) %>%
  group_by(month, year, site, side, RESERVE, sample, key, sciname) %>%
  summarize(count = mean(count)) %>%
  left_join(data1, by = "key") %>%
  mutate(Ni = count/N,
         date = as.Date(paste(1, month, year, sep = "/"), format = "%d/%m/%Y")) %>%
  ungroup()

ggplotly(ggplot(data2, aes(x = date, y = Ni, color = sciname)) +
           geom_line() +
           geom_point() +
           theme(legend.position="none"))

data3 <- ungroup(data2) %>%
  select(sciname, key, N) %>%
  spread(key, N)

rows <- data3$sciname
data3<-select(data3, -sciname)
data3 <- as.data.frame(data3)
rownames(data3) <- rows
data3 <- as.data.frame(data3)
data3[is.na(data3)] = 0
JC <- bvi()

