check <- lm(log_density ~ ., data = x_seen)

check

if (sum(is.na(check$coefficients)) > 0){warning('some model terms are perfectly collinear')}

check$coefficients[is.na(check$coefficients)]

test <- broom::tidy(check) %>%
  mutate(variable = colnames(x_seen)) %>%
  mutate(lower = estimate - std.error * 1.96,
         upper = estimate + std.error * 1.96) %>%
  filter(variable %in% as.character(2001:2013)) %>%
  mutate(year = as.numeric(variable))

test %>%
  ggplot() +
  geom_pointrange(aes(x = year, y = estimate, ymin = lower, ymax = upper))+
  geom_vline(aes(xintercept = 2003)) +
  geom_hline(aes(yintercept = 0), color = 'red')


check <- glm(any_seen ~ ., data = x_seeing, family = binomial)

check

sum(is.na(check$coefficients))

check$coefficients[is.na(check$coefficients)]

test <- broom::tidy(check) %>%
  slice(31:43) %>%
  mutate(lower = estimate - std.error * 1.96,
         upper = estimate + std.error * 1.96,
         year = 2001:2013)

test %>%
  ggplot() +
  geom_pointrange(aes(x = year, y = estimate, ymin = lower, ymax = upper))+
  geom_vline(aes(xintercept = 2003))

