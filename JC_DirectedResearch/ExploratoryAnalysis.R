
# Load packages
library(tidyverse)
library(plotly)
library(Rbvi)

# Load the list of sites and their categories
sites <- read.csv("./MLPA Data/Final_Site_Table_UCSB.csv", stringsAsFactors = F) %>%
  select(SITE_SIDE, RESERVE)

# Load the master spreadsheet for species
master_spp <- read.csv("./MLPA Data/master_spp_table.csv", stringsAsFactors = F) %>%
  group_by(classcode, genus, species) %>%
  summarize(N = n()) %>%
  select(classcode, genus, species, N) %>%
  ungroup()

 #Identify duplicated spp
duplicated_sp <- filter(master_spp, N>1)

# Get rid of the N column
master_spp <- select(master_spp, -N)

# Load the transect database, and filter it to keep only Bottom surveys
UCSB_fish <- read.csv("./MLPA Data/UCSB_FISH raw thru 2013.csv", stringsAsFactors = F) %>%
  filter(level == "BOT")

# Calculate averages per survey (average across transects)
data2 <- UCSB_fish %>%
  left_join(master_spp, by = "classcode") %>%
  mutate(SITE_SIDE = paste(site, side, sep = "_")) %>%
  left_join(sites, by = "SITE_SIDE") %>%
  mutate(sample = paste(day, month, year, site, side, RESERVE, zone, level, transect, sep = "-"),
         sciname = paste(genus, species)) %>%
  filter(!sciname == " ") %>%
  group_by(day, month, year, site, side, RESERVE, zone, level, transect, sample, sciname) %>%
  summarize(count = sum(count)) %>%
  group_by(day, month, year, site, side, RESERVE, zone, sciname) %>%
  summarize(N = mean(count)) %>%
  ungroup() %>%
  mutate(key = paste(day, month, year, site, side, RESERVE, zone))

# Keep the identifier, abundance and scientific name
data3 <- ungroup(data2) %>%
  select(sciname, key, N) %>%
  spread(key, N)

# Format this data to match the bvi input format
rows <- data3$sciname
data3<-select(data3, -sciname)
data3 <- as.data.frame(data3)
rownames(data3) <- rows
data3 <- as.data.frame(data3)
data3[is.na(data3)] = 0

# Calculate bvi
BVI <- bvi(data3)

# Keep only BVI score and %BVI scores

BVI2 <- BVI[,c(1,2837:2838)]
BVI2
  

