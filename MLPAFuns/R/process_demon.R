#' analyze demon object
#'
#' \code{process_demon} loads saved data
#' from run_demon and produces diagnostics
#' fits etc
#' @param runfolder the folder where results are stored
#' @param fontsize numeric fontsize to use
#' @param post_sample_size the desired final number of samples
#' from the MCMC after burning and thinning
#' @param burn proportion 0 to 1 of iterations to burn

process_demon <- function(run_name,fontsize = 14,post_sample_size = 1000, burn = 0.6)
{

  runpath <- paste('MLPA Effects Results/',run_name,'/', sep = '')

  load(paste(runpath,'species_siteside_year.Rdata', sep = ''))

  load(paste(runpath,'MCMC results.Rdata', sep = ''))

  plot_theme <- theme_classic() + theme(text = element_text(size = fontsize, family = 'Helvetica'))

  post <- bayes_reg$demon_fit$Posterior1

  its <- dim(post)[1]

  burn <- burn

  thinned_post <- thin_mcmc(post[(burn * its):its, ], thin_every = (its*(1-burn))/post_sample_size)

  predictions <- apply_demon(demonpost = thinned_post, demon = bayes_reg$demon_fit, dat = bayes_reg$Data)

  resids <- predictions$posterior %>%
    subset(obs_log_density > min(obs_log_density)) %>%
    group_by(observation) %>%
    summarise(mean_resid = mean(resid), mean_obs = mean(obs_log_density),mean_pred = mean(pred_log_density), mean_post_pred = mean(post_predict))

  resid_plot <- resids %>%
    ggplot(aes(mean_resid)) +
    geom_histogram(fill = 'steelblue4',color = 'black') +
    plot_theme +
    xlab('Residual') +
    ylab('Count')

  qq_plot <- resids %>%
    ggplot(aes(sample = mean_resid)) +
    stat_qq(shape = 21, size = 4, alpha = 0.6, fill = 'steelblue4') +
    plot_theme +
    xlab('Theoretical') +
    ylab('Sample')

  misspec_plot <- resids %>%
    ggplot(aes(mean_pred,mean_resid)) +
    geom_point(shape = 21, fill = 'steelblue4', size = 3, alpha = 0.5) +   geom_hline(aes(yintercept = 0), linetype = 'longdash') +
    geom_smooth(method = 'lm', color = 'red', se = F) +
    plot_theme +
    xlab('Predicted Density') +
    ylab('Residual')


  lm_eqn = function(m) {

    l <- list(a = format(coef(m)[1], digits = 2),
              b = format(abs(coef(m)[2]), digits = 2),
              r2 = format(summary(m)$r.squared, digits = 3));

    if (coef(m)[2] >= 0)  {
      eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
    } else {
      eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
    }

    as.character(as.expression(eq));
  }

  fit_plot <- resids %>%
    ggplot(aes(mean_obs,mean_pred)) +
    geom_point(shape = 21, fill = 'steelblue4', size = 3, alpha = 0.5) +
    geom_abline(aes(slope = 1, intercept = 0)) +
    geom_smooth(method = 'lm', color = 'red',linetype = 'longdash') +
    plot_theme +
    xlab('Observed') +
    ylab('Predicted')
  #+ geom_text(x = -8,y= -3, label = lm_eqn(lm(mean_pred ~ mean_obs,resids )), parse = T)

  post_pred_plot <- resids %>%
    ggplot(aes(mean_obs,mean_post_pred)) +
    geom_point(shape = 21, fill = 'steelblue4', size = 3, alpha = 0.5) +
    geom_abline(aes(slope = 1, intercept = 0)) +
    geom_smooth(method = 'lm', color = 'red',linetype = 'longdash') +
    plot_theme +
    xlab('Observed') +
    ylab('Posterior Predicted')

  r2 <- round(summary(lm(mean_post_pred ~ mean_obs,resids ))$r.squared, digits = 2)

  gweke_plot <- ggs_geweke(ggs(mcmc(thinned_post))) +
    plot_theme

  eff_sample_plot <- data.frame(effective.ss = effectiveSize(mcmc(thinned_post))) %>%
    mutate(parameter = rownames(.)) %>%
    subset(!parameter %in% c('ll','deviance')) %>%
    ggplot(aes(factor(parameter),effective.ss)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    coord_flip() +
    xlab('Parameter') +
    ylab('Effective Sample Size') +
    plot_theme


  lag_1 <- rep(NA,dim(thinned_post)[2])
  for (j in 1:dim(thinned_post)[2])
  {
    lag_1[j] <- acf(thinned_post[,j], plot = F)$acf[2]
  }


  bino_plot <- predictions$post_prob_zero %>%
    group_by(observation) %>%
    summarize(is_zero = mean(is_zero), prob_zero = mean(prob_zero)) %>%
    mutate(bin = cut(prob_zero,5, dig.lab = 2)) %>%
    group_by(bin) %>%
    summarize(seen = mean(is_zero)) %>%
    ggplot(aes(bin,seen)) +
    geom_bar(stat = 'identity', position = 'dodge',
             fill = 'steelblue4', color = 'black') +
    plot_theme +
    xlab('Proportion Observed > 0') +
    ylab('Mean P(D)')


  acf_hist_plot <- ggplot(data.frame(acf = lag_1), aes(acf)) +
    geom_histogram(binwidth = .05, fill = 'steelblue4', color = 'black') +
    xlab('Lag 1 ACF') +
    ylab('Count') +
    plot_theme

  crosscor_plot <- ggs_crosscorrelation(ggs(mcmc(thinned_post))) +
    plot_theme


  outs_of_interest_plot <- as.data.frame(thinned_post) %>%
    dplyr::select(fished,mpa_applied,fished_x_mpa,fished_x_yearsmpa) %>%
    rename(Fished = fished,'MPA Applied' = mpa_applied, 'Fished and Applied' = fished_x_mpa,
           'Fished by Years Protected' = fished_x_yearsmpa) %>%
    gather('Variable','Coefficient') %>%
    group_by(Variable) %>%
    mutate(lower95 = quantile(Coefficient, 0.025), upper95 = quantile(Coefficient, 0.975)) %>%
    ggplot(aes(Coefficient, fill = Variable)) +
    scale_fill_brewer(guide = F, palette = 'Spectral') +
    geom_density() +
    geom_vline(aes(xintercept = c(lower95,upper95)),alpha = 0.75) +
    geom_vline(aes(xintercept = 0), linetype = 'longdash', color = 'red')+
    facet_wrap(~Variable) +
    theme_light() +
    theme(text = element_text(size = 18)) +
    ylab('Density')

  outs_of_interest_binomial_plot <- as.data.frame(thinned_post) %>%
    dplyr::select(bi.fished,bi.mpa_applied,bi.fished_x_mpa) %>%
    rename(Fished = bi.fished,'MPA Applied' = bi.mpa_applied, 'Fished and Applied' = bi.fished_x_mpa) %>%
    gather('Variable','Coefficient') %>%
    group_by(Variable) %>%
    mutate(lower95 = quantile(Coefficient, 0.025), upper95 = quantile(Coefficient, 0.975)) %>%
    ggplot(aes(Coefficient, fill = Variable)) +
    scale_fill_brewer(guide = F, palette = 'Spectral') +
    geom_density() +
    geom_vline(aes(xintercept = c(lower95,upper95)),alpha = 0.75) +
    geom_vline(aes(xintercept = 0), linetype = 'longdash', color = 'red')+
    facet_grid(Variable~.) +
    theme_light() +
    theme(text = element_text(size = 18)) +
    ylab('Density')

  year_trend_plot <- as.data.frame(thinned_post) %>%
    select(contains('factor_year')) %>%
    select(which(!grepl('bi.',colnames(.), fixed = T))) %>%
    gather('factor_year','Coef') %>%
    mutate(Year = as.numeric(factor_year) + 2001 -1) %>%
    ggplot(aes(factor(Year),Coef)) +
    geom_boxplot(fill = 'steelblue4', color = 'black') +
    geom_hline(aes(yintercept = 0))+
    plot_theme +
    theme(axis.text.x = element_text(size = 12)) +
    xlab('Year') +
    ylab('Effect Relative to 2000')

  local_files <- ls()

  plot_files <- local_files[grep('_plot', local_files, fixed = T)]

  plot_list <- list()

  for (i in 1:length(plot_files))
  {
    eval(parse( text = paste('plot_list$',plot_files[i],' <- ',plot_files[i], sep = '')))
  }


return(list(plot_list = plot_list, resids = resids, thinned_post = thinned_post,
            predictions = predictions))
}