huh <- processed_dat %>%
  group_by(site,year,common_name) %>%
  summarize(total_bio = sum(biomass, na.rm = T)) %>%
  ungroup() %>%
  group_by(site,common_name) %>%
  summarize(num_obs = sum(total_bio>0), possible = length(total_bio)) %>%
  mutate(perc_obs = num_obs/possible) %>%
  ungroup() %>%
  filter(perc_obs > 0.25) %>%
  arrange(num_obs)