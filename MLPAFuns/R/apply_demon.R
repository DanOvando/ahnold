#' Function to estimate predicted values from
#'
#' \code{apply_demon} applies the posterior from
#' a Laplaces Demon call to a supplied and appropriate
#' dataframe to obtain predictions and posterior predictions
#' @param demonpost the burned and thinned posterior
#' @param demon an object of class demonoid produced by
#' LaplacesDemon
#' @param dat the dataframe associated with the
#' demonoid object

apply_demon <- function(demonpost,demon,dat)
{

  #   demonpost = trimmed_post
  #
  #   demon <- bayes_reg$demon_fit
  #
  #   dat <- bayes_reg$Data

  # Obtain predicted values ----

  betas <- demonpost[,dat$pos_den_beta]

  ind_vars <- dat$reg_dat[,dat$beta_to_use_binom == F]

  predicted_density <- as.data.frame(t(betas %*% t(ind_vars)))

  predicted_density$obs_log_density <- dat$dep_var

  sigma_density <- as.data.frame(demonpost[,'sigma_density'])

  colnames(sigma_density) <- 'sigma_density'

  sigma_density$chain <- 1:dim(sigma_density)[1]

  obs <- dim(predicted_density)[1]

  predicted_density$observation <- 1:obs

  posterior <- predicted_density %>%
    gather('chain','pred_log_density',which(grepl('V', colnames(.), fixed = T))) %>%
    mutate(chain = as.numeric(gsub('V', '',chain))) %>%
    left_join(sigma_density, by = 'chain') %>%
    mutate(post_predict = rnorm(length(pred_log_density),pred_log_density, sigma_density),
           resid = pred_log_density - obs_log_density)

  # Obtain predicted zero probs ---------------------------------------------

  betas <- demonpost[,dat$beta_to_use_binom]

  ind_vars <- dat$reg_dat[,dat$beta_to_use_binom == T]

  prob_zero <- as.data.frame(t(betas %*% t(ind_vars)))

  prob_zero <- exp(prob_zero) / (1+exp(prob_zero))

  prob_zero$log_density <- dat$dep_var

  prob_zero$is_zero <- prob_zero$log_density > min(prob_zero$log_density)

  obs <- dim(prob_zero)[1]

  prob_zero$observation <- 1:obs

  post_prob_zero <- prob_zero %>%
    gather('chain','prob_zero',which(grepl('V', colnames(.), fixed = T))) %>%
    mutate(chain = as.numeric(gsub('V', '',chain)))





  # Scratch Plots -----------------------------------------------------------


#   resid_plot <- posterior %>%
#     group_by(observation) %>%
#     summarise(obs = mean(obs_log_density), pred = mean(pred_log_density),
#               post_pred = mean(post_predict), mean_resid = mean(resid)) %>%
#     ggplot(aes(pred,mean_resid)) +
#     geom_point(shape = 21, alpha = 0.5, fill = 'blue') +
#     geom_abline(intercept = 0, slope = 0, color = 'red')
#
#
#   obs_pred_plot <- posterior %>%
#     group_by(observation) %>%
#     summarise(obs = mean(obs_log_density), pred = mean(pred_log_density),
#               post_pred = mean(post_predict)) %>%
#     subset(obs>min(obs)) %>%
#     ggplot(aes(obs,pred)) +
#     geom_point(shape = 21, alpha = 0.5, fill = 'blue') +
#     geom_abline(intercept = 0, slope = 1, color = 'red') +
#     geom_smooth(method = 'lm')
#
#   please <- posterior %>%
#     group_by(observation) %>%
#     summarise(obs = mean(obs_log_density), pred = mean(pred_log_density),
#               post_pred = mean(post_predict)) %>%
#     ggplot(aes(obs,post_pred)) +
#     geom_point(shape = 21, alpha = 0.5, fill = 'blue') +
#     geom_abline(intercept = 0, slope = 1, color = 'red')
#
#   plots <- list(zero_plot = zero_plot, resid_plot = resid_plot, obs_pred_plot = obs_pred_plot)

  return(list(posterior = posterior, post_prob_zero = post_prob_zero))
}