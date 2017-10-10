fit_fish <- function(data, ind_vars, dep_var, family,fit_seen = T){

  if (fit_seen == T){
    data <- data %>%
      filter(any_seen == fit_seen & is.na(any_seen) == F)
  } else{
    data <- data %>%
      filter(is.na(any_seen) == F)
  }
  reg_fmla <- paste(dep_var, ind_vars, sep = '~') %>% as.formula()
  model <- glm(reg_fmla, data  = data, family = family)
}