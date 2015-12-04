species_siteside_year <- processed_dat %>% #group data to the species, site_side, year level
  subset(is.na(targeted) == F & targeted != '' & mpaareanm2 >0) %>%
  group_by(year,site_side,classcode) %>%
  summarise(mpa.period = unique(region_has_mpas),
            mean_density = mean(biomass, na.rm = T),
            site.type = unique(mpa.status),
            years_mpa = max(0,(year - (year.mpa-1)) * as.numeric(year.mpa>0)),
            region = unique(region),
            targeted = unique(targeted),
            trophic.group = unique(trophicgroup),
            mpa_area = mean(mpaareanm2, na.rm = T),
            linf = mean(c(vbgf.linf,vbgf.linf.f,vbgf.linf.m), na.rm = T),
            vbk = mean(c(vbgf.k,vbgf.k.f,vbgf.k.m), na.rm = T),
            size_mature = mean(size_mature_cm, na.rm = T),
            broadtrophic = unique(broadtrophic),
            mean_temp = mean(mean_temp, na.rm = T),
            mean_vis = mean(mean_vis, na.rm = T),
            species = unique(common_name)) %>%
  #   ungroup() %>% min(mean_density[mean_density>0], na.rm = T)
  ungroup() %>%
  group_by(site_side, species)  %>%
  #   mutate(anybio = sum(mean_density, na.rm = T)) %>%
  #   subset(anybio >0) %>%
  #   mutate(mean_density = mean_density + 0.01
  mutate(mean_density = pmax(quantile(unique(mean_density[mean_density >0]),.01, na.rm = T),mean_density)
         ,log_density = log(mean_density),
         fished = as.numeric(targeted == 'Targeted'),
         mpa_applied = as.numeric(mpa.period == 'TRUE'),
         fished_x_mpa = fished * mpa_applied,
         fished_x_yearsmpa = fished * years_mpa,
         factor_year = as.factor(year)) %>%
  group_by(site_side) %>%
  mutate(eventual_mpa = max(years_mpa) >0) %>%
  subset(is.na(mean_density) == F & is.na(log_density) == F) %>%
  ungroup()