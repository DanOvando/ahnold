
parma <- Initial.Values


cppFunction('NumericVector floookys(NumericVector parm, NumericVector pos_den_beta,
            NumericVector pos_sigma_year,NumericVector pos_sigma_bi_year, NumericVector pos_sigma_region,
            NumericVector pos_sigma_density, NumericVector pos_den_time_terms,
            NumericVector pos_bi_time_terms,NumericVector pos_den_region_terms,
            NumericVector pos_beta_to_use_binom, NumericMatrix reg_dat, NumericVector binom_dep_var,
            NumericMatrix bi_reg_mat, NumericMatrix den_reg_mat, NumericVector dep_var)
            {

            // init things you need

            NumericVector beta(pos_den_beta.size());

            NumericVector sigma_year(pos_sigma_year.size());

            NumericVector sigma_bi_year(pos_sigma_bi_year.size());

            NumericVector sigma_region(pos_sigma_region.size());

            NumericVector sigma_density(pos_sigma_density.size());

            beta = parm[pos_den_beta];

            sigma_year = parm[pos_sigma_year];

            sigma_bi_year = parm[pos_sigma_bi_year];

            sigma_region = parm[pos_sigma_region];

            sigma_density = parm[pos_sigma_density];

            // Log priors

            NumericVector den_time_terms = parm[pos_den_time_terms];

            NumericVector bi_time_temrs = parm[pos_bi_time_terms];

            NumericVector den_region_terms = parm[pos_den_region_terms];

            double year_priors = sum(dnorm(den_time_terms, 0.0, sigma_year[0], true)) ;

            double bi_year_priors = sum(dnorm(bi_time_temrs,0,sigma_bi_year[0], true));

            double sigma_year_prior = sum(dnorm(sigma_year, 0.1,0.2, true));

            double sigma_bi_year_prior = sum(dnorm(sigma_bi_year, .1,.2, true));

            double region_priors = sum(dnorm(den_region_terms,0,sigma_region[0], true));

            double sigma_region_prior = sum(dnorm(sigma_region, 0.1,0.2, true));

            double sigma_density_prior = sum(dnorm(sigma_density, 0.1, 0.2, true));

            //   ### Hurdle Log-Likelihood ----

            NumericVector bi_beta = parm[pos_beta_to_use_binom];

            int nrows = reg_dat.nrow();

            int ncols = pos_beta_to_use_binom.size();

            //NumericMatrix bi_reg_mat(nrows, ncols);

            NumericVector bi_hat(nrows);

            NumericVector temp_hat(ncols);

            for(int i = 0; i < nrows; ++i) {

            for(int j = 0; j <ncols; ++j)
            {
            //bi_reg_mat(i,j) = reg_dat(i,j);

            temp_hat[j] = bi_reg_mat(i,j) * bi_beta[j];
            }
            bi_hat[i] = sum(temp_hat);

            }

            NumericVector prob_hat = exp(bi_hat)/(1 + exp(bi_hat));

            double bi_loglike =  sum(dbinom(binom_dep_var,1,prob_hat[0], true));

            // Density likelihood

            NumericVector den_beta = parm[pos_den_beta];

            int den_ncols = den_beta.size();

            NumericVector mu(nrows);

            NumericVector temp_mu(ncols);

            for(int k = 0; k < nrows; ++k) {

            for(int l = 0; l <den_ncols; ++l)
            {
            temp_mu[l] = den_reg_mat(k,l) * den_beta[l];
            }
            mu[k] = sum(temp_mu);

            }

            double density_loglike = sum(pow(dnorm(dep_var, mu[0],sigma_density[0], true),binom_dep_var[0]));

            // Posterior

            double LP = density_loglike + bi_loglike  + year_priors + bi_year_priors +
            sigma_year_prior + sigma_bi_year_prior + region_priors + sigma_region_prior +
            sigma_density_prior;

            NumericVector out(2);

            out[1] = LP;

            out[2] = -2*(density_loglike + bi_loglike);

            return(out);
            }')

# a <- proc.time()
damnit = floookys(parm = parma, pos_den_beta = (pos_den_beta - 1),
                  pos_sigma_year = pos_sigma_year - 1,pos_sigma_bi_year = pos_sigma_bi_year - 1,
                  pos_sigma_region = pos_sigma_region - 1,pos_sigma_density = pos_sigma_density - 1,
                  pos_den_time_terms = pos_den_time_terms - 1, pos_bi_time_terms = pos_bi_time_terms - 1,
                  pos_den_region_terms = pos_den_region_terms - 1,
                  pos_beta_to_use_binom = pos_beta_to_use_binom-1,
                  reg_dat = Data$reg_dat, binom_dep_var = Data$binom_dep_var,
                  bi_reg_mat = bi_reg_mat, den_reg_mat = den_reg_mat, dep_var = Data$dep_var)
# proc.time() - a

