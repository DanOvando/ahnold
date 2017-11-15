rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(reshape2)

data<-read.csv("fish_integrated.csv",na.strings = "",stringsAsFactors = FALSE)

data$date <- as.Date(data$date)

data1 <- data %>% 
  mutate(year = format(date,'%Y')) %>%
  filter( year %in% (2005:2014))

#fish diversity calculation
# data2<-data1 %>%
#   filter(!count==0) %>%
#   group_by(site_id,year,auth_taxon_id) %>%
#   summarise(freq=n()) %>%
#   ungroup()%>%
#   dcast(site_id+year~auth_taxon_id,length,value.var="freq")
# 
# diversity<-rowSums(data2[,3:ncol(data2)])
# 
# data3<-cbind(data2[,1:2],diversity)

data4<-data1 %>%
  mutate(count=as.numeric(count)) %>%
  mutate(density=count/area) %>% 
  #these two methods were excluded because there are larger scale fish survey in the same program
  filter(!sample_method=="visualfish"|sample_method=="crypticfish")%>% 
  group_by(site_id,subsite_id,transect_id,replicate_id,date,year) %>%
  summarise(density=sum(density)) %>% #sum up fish density in the same plot
  ungroup() %>%
  group_by(site_id,year) %>%
  summarise(density=mean(density)) %>% #average fish density over different plots
  ungroup() 

site<-read.csv("site_geolocation.csv") 

data5<-data4 %>%
  left_join(select(site,site_id,geolocation),by="site_id") %>%
  select(site_id,geolocation,year,density)

write.csv(data5,"fish_density_mapping.csv",row.names = F)