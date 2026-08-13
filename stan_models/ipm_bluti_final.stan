// Integrated population model for blue tits
functions {
  real real_poisson_lpdf(real n, real lambda) {
    real lp;
    if (lambda < 0) {
      reject("lambda must be non-negative; found lambda=", lambda);
    } else if (n < 0) {
      reject("n must not be negative; found n=", n);
    } else {
      lp = n * log(lambda) - lambda - lgamma(n + 1);
    }
    return lp;
  }

  real real_binomial_lpdf(real n, real N, real theta) {
    real lp;
    if (N < 0) {
      reject("N must be non-negative; found N=", N);
    } else if (theta < 0 || theta > 1) {
      reject("theta must be in [0,1]; found theta=", theta);
    } else if (n < 0) {
      reject("n must not be negative; found n=", n);
    } else {
      real n_capped = fmin(n, N);
      lp = lchoose(N, n_capped) + n_capped * log(theta) + (N - n_capped) * log(1 - theta);
    }
    return lp;
  }

  // marray_adults 
  array[] vector marray_adults(int nyears, vector phia, vector p_s) {
    int ny_minus_1 = nyears - 1;
    array[ny_minus_1] vector[nyears] pr_a;
    vector[nyears - 1] q = 1 - p_s;
    real prod_phi;
    real prod_q;

    for (t in 1:(nyears - 1)) {
      pr_a[t, t] = phia[t] * p_s[t];

      prod_phi = phia[t];
      prod_q   = 1;
      for (j in (t + 1):(nyears - 1)) {
        prod_phi   = prod_phi * phia[j];
        prod_q     = prod_q   * q[j - 1];
        pr_a[t, j] = prod_phi * prod_q * p_s[j];
      }

      for (j in 1:(t - 1))
        pr_a[t, j] = 0;

      pr_a[t, nyears] = 1 - sum(pr_a[t, 1:(nyears - 1)]);
    }
    return pr_a;
  }
}

data {
  int nyears;
  int nsites;
  array[nsites, nyears]             int y;        // pop counts or max number of nestboxes occupied [site, year]
  array[nsites, nyears - 1]         int J;        // number of fledglings [site, year]
  array[nsites, nyears - 1]         int R;        // broods surveyed [site, year]
  array[nsites, nyears - 1, nyears] int marray_a; // adult m-array 
}

transformed data {
  int ny_minus_1 = nyears - 1;
}

parameters {
  // Latent state variables (N1, NadSurv, Nadimm) varying by site x year
  array[nsites] vector<lower=0>[nyears] N1;
  array[nsites] vector<lower=0>[nyears] NadSurv;
  array[nsites] vector<lower=0>[nyears] Nadimm;

  // Grand means (logit or log scale)
  real l_mphij;
  real l_mphia;
  real l_mprod;
  real l_mim;
  real l_p;             

  // Year random effects (non-centered)
  vector[ny_minus_1] epsilon_phia_raw;
  vector[ny_minus_1] epsilon_prod_raw;
  vector[ny_minus_1] epsilon_im_raw;
  vector[ny_minus_1] epsilon_p_raw;    

  // Site random effects (non-centered)
  vector[nsites] zeta_phia_raw;
  vector[nsites] zeta_prod_raw;
  vector[nsites] zeta_im_raw;
  vector[nsites] zeta_p_raw;           

  // Site×year interaction random effects (non-centered)
  array[nsites] vector[ny_minus_1] eta_phia_raw;
  array[nsites] vector[ny_minus_1] eta_prod_raw;
  array[nsites] vector[ny_minus_1] eta_im_raw;
  array[nsites] vector[ny_minus_1] eta_p_raw;   

  // Standard deviations — year level
  real<lower=0> sig_year_phia;
  real<lower=0> sig_year_prod;
  real<lower=0> sig_year_im;
  real<lower=0> sig_year_p;                 

  // Standard deviations — site level
  real<lower=0> sig_site_phia;
  real<lower=0> sig_site_prod;
  real<lower=0> sig_site_im;
  real<lower=0> sig_site_p;           

  // Standard deviations — site×year level
  real<lower=0> sig_sy_phia;
  real<lower=0> sig_sy_prod;
  real<lower=0> sig_sy_im;
  real<lower=0> sig_sy_p;             
}

transformed parameters {
  // Scaled year deviations
  vector[ny_minus_1] epsilon_phia = sig_year_phia * epsilon_phia_raw;
  vector[ny_minus_1] epsilon_prod  = sig_year_prod  * epsilon_prod_raw;
  vector[ny_minus_1] epsilon_im   = sig_year_im   * epsilon_im_raw;
  vector[ny_minus_1] epsilon_p    = sig_year_p    * epsilon_p_raw;   

  // Scaled site deviations
  vector[nsites] zeta_phia = sig_site_phia * zeta_phia_raw;
  vector[nsites] zeta_prod  = sig_site_prod  * zeta_prod_raw;
  vector[nsites] zeta_im   = sig_site_im   * zeta_im_raw;
  vector[nsites] zeta_p    = sig_site_p    * zeta_p_raw;       

  // Site×year interaction deviations
  array[nsites] vector[ny_minus_1] eta_phia;
  array[nsites] vector[ny_minus_1] eta_prod;
  array[nsites] vector[ny_minus_1] eta_im;
  array[nsites] vector[ny_minus_1] eta_p;                   

  // Demographic rates [site, year]
  array[nsites] vector<lower=0, upper=1>[ny_minus_1] phij;
  array[nsites] vector<lower=0, upper=1>[ny_minus_1] phia;
  array[nsites] vector<lower=0>[ny_minus_1]          f;
  array[nsites] vector<lower=0>[ny_minus_1]          omega;

  // Recapture probability
  array[nsites] vector<lower=0, upper=1>[ny_minus_1] p;        

  // Derived measures
  array[nsites] vector<lower=0>[nyears]     Ntot;
  array[nsites] vector<lower=0>[ny_minus_1] rho;
  array[nsites, ny_minus_1] simplex[nyears] pr_a;

  for (s in 1:nsites) {
    // Scale interaction terms
    eta_phia[s] = sig_sy_phia * eta_phia_raw[s];
    eta_prod[s]  = sig_sy_prod  * eta_prod_raw[s];
    eta_im[s]   = sig_sy_im   * eta_im_raw[s];
    eta_p[s]    = sig_sy_p    * eta_p_raw[s];                  

    for (t in 1:ny_minus_1) {
      // phij: grand mean + assuming adult year/site/site×year effects 
      phij[s, t] = inv_logit(l_mphij
                              + epsilon_phia[t]
                              + zeta_phia[s]
                              + eta_phia[s][t]);
                              
 // grand mean + year + site + site×year random effects
      phia[s, t] = inv_logit(l_mphia
                              + epsilon_phia[t]
                              + zeta_phia[s]
                              + eta_phia[s][t]);
 // grand mean + year + site + site×year random effects
      f[s, t]    = exp(l_mprod
                       + epsilon_prod[t]
                       + zeta_prod[s]
                       + eta_prod[s][t]);
                       
 // grand mean + year + site + site×year random effects
      omega[s, t] = exp(l_mim
                        + epsilon_im[t]
                        + zeta_im[s]
                        + eta_im[s][t]);

      // grand mean + year + site + site×year random effects
      p[s, t]    = inv_logit(l_p                               
                              + epsilon_p[t]
                              + zeta_p[s]
                              + eta_p[s][t]);
    }

    // Total population
    Ntot[s] = NadSurv[s] + Nadimm[s] + N1[s];

    // Adult m-array cell probabilities 
    pr_a[s] = marray_adults(nyears, phia[s], p[s]);             

    // Productivity or expected juveniles: R[s,t] broods x f
    for (t in 1:ny_minus_1)
      rho[s, t] = R[s, t] * f[s, t];
  }
}

model {
  // Priors: initial population sizes
  for (s in 1:nsites) {
    N1[s, 1]      ~ normal(50, 100);
    NadSurv[s, 1] ~ normal(50, 100);
    Nadimm[s, 1]  ~ normal(50, 100);
  }

  // Prior on juvenile survival (moderately informative) - on logit scale 
  //(use plogis to understand what the numbers mean on probability scale, easier to interpret)
  l_mphij ~ normal(-1.1, 0.7);

  // Priors: grand means 
  //Back-transform to help with interpretation of what below numbers mean
  l_mphia ~ normal(0, 1.5); //logit scale
  l_mprod  ~ normal(1.5, 0.7); //log scale
  l_mim   ~ normal(-3, 3); //log scale
  l_p     ~ normal(0, 1.5); //logit scale
 

  // Priors: year random effects (non-centered)
  epsilon_phia_raw ~ normal(0, 1);
  epsilon_prod_raw  ~ normal(0, 1);
  epsilon_im_raw   ~ normal(0, 1);
  epsilon_p_raw    ~ normal(0, 1);   

  // Priors: site random effects (non-centered)
  zeta_phia_raw ~ normal(0, 1);
  zeta_prod_raw  ~ normal(0, 1);
  zeta_im_raw   ~ normal(0, 1);
  zeta_p_raw    ~ normal(0, 1);      

  // Priors: site×year interaction random effects (non-centered)
  for (s in 1:nsites) {
    eta_phia_raw[s] ~ normal(0, 1);
    eta_prod_raw[s]  ~ normal(0, 1);
    eta_im_raw[s]   ~ normal(0, 1);
    eta_p_raw[s]    ~ normal(0, 1);  
  }

  // Weakly informative priors for std dev
  sig_year_phia      ~ normal(0, 1);
  sig_year_prod       ~ normal(0, 1);
  sig_year_im        ~ normal(0, 1);
  sig_year_p         ~ normal(0, 1);      
  sig_site_phia ~ normal(0, 1);
  sig_site_prod  ~ normal(0, 1);
  sig_site_im   ~ normal(0, 1);
  sig_site_p    ~ normal(0, 1);      
  sig_sy_phia   ~ normal(0, 1);
  sig_sy_prod    ~ normal(0, 1);
  sig_sy_im     ~ normal(0, 1);
  sig_sy_p      ~ normal(0, 1);    

  // Likelihood: state process 
  for (s in 1:nsites) {
    for (t in 2:nyears) {
      real mean1 = 0.5 * f[s, t-1] * phij[s, t-1] * Ntot[s, t-1];
      real mpo   = Ntot[s, t-1] * omega[s, t-1];

      N1[s, t]      ~ real_poisson(mean1);
      NadSurv[s, t] ~ real_binomial(Ntot[s, t-1], phia[s, t-1]);
      Nadimm[s, t]  ~ real_poisson(mpo);
    }
  }

  // Likelihood: population counts
  for (s in 1:nsites)
    y[s] ~ poisson(Ntot[s]);

  // Likelihood: adult CJS m-array
  for (s in 1:nsites) {
    for (t in 1:ny_minus_1) {
      marray_a[s, t] ~ multinomial(pr_a[s, t]);
    }
  }

  // Likelihood: productivity
  for (s in 1:nsites)
    J[s] ~ poisson(rho[s]);
}

generated quantities {
  // Grand mean demographic rates (back-transformed)
  real<lower=0, upper=1> mphij = inv_logit(l_mphij);  // grand mean juv survival
  real<lower=0, upper=1> mphia = inv_logit(l_mphia); // grand mean adult survival
  real<lower=0, upper=1> mp    = inv_logit(l_p); // grand mean recapture
  real<lower=0>          mprod  = exp(l_mprod); // grand mean productivity
  real<lower=0>          mim   = exp(l_mim); // grand immigration rate

  // Population growth rates per site
  array[nsites] vector<lower=0>[ny_minus_1] lambda; //annual pop growth rate
  array[nsites] vector[ny_minus_1]          logla; //instantaneous pop growth rate (i.e. r=log(lambda))
  array[nsites] real<lower=0>               mlam; //long-term pop growth rate (geometric mean of lambda)

  for (s in 1:nsites) {
    lambda[s] = Ntot[s, 2:nyears] ./ Ntot[s, 1:ny_minus_1]; //annual pop growth rate per site per year
    logla[s]  = log(lambda[s]); //instantaneous pop growth rate (i.e. r=log(lambda))
    mlam[s]   = exp(1.0 / ny_minus_1 * sum(logla[s])); //long-term pop growth rate (geometric mean of lambda)
  }

  // Variance components and ICC
  
  //phia = adult survival
  real var_phia_year     = square(sig_year_phia);
  real var_phia_site     = square(sig_site_phia);
  real var_phia_siteyear = square(sig_sy_phia);
  real icc_phia          = var_phia_year / (var_phia_year + var_phia_siteyear);
  
  // prod = productivity
  real var_prod_year      = square(sig_year_prod);
  real var_prod_site     = square(sig_site_prod);
  real var_prod_siteyear  = square(sig_sy_prod);
  real icc_prod           = var_prod_year / (var_prod_year + var_prod_siteyear);
 
 // im = immigration
  real var_im_year       = square(sig_year_im);
  real var_im_site     = square(sig_site_im);
  real var_im_siteyear   = square(sig_sy_im);
  real icc_im            = var_im_year / (var_im_year + var_im_siteyear);

  //p = recapture probability (unecessary but whatevs)...
  real var_p_year        = square(sig_year_p);
  real var_p_site     = square(sig_site_p);
  real var_p_siteyear    = square(sig_sy_p);
  real icc_p             = var_p_year / (var_p_year + var_p_siteyear);
}

