// Integrated population model
// Site, year, and site×year random effects

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
    } else if (n < 0 || n > N) {
      reject("n must be in [0,N]; found n=", n);
    } else {
      lp = lchoose(N, n) + n * log(theta) + (N - n) * log(1 - theta);
    }
    return lp;
  }


  array[] vector marray_adults(int nyears, vector phia, vector p) {
    int ny_minus_1 = nyears - 1;
    array[ny_minus_1] vector[nyears] pr_a;
    vector[nyears - 1] q = 1 - p;
    real prod_phi;
    real prod_q;

    for (t in 1:(nyears - 1)) {
      pr_a[t, t] = phia[t] * p[t];

      prod_phi = phia[t];
      prod_q   = 1;
      for (j in (t + 1):(nyears - 1)) {
        prod_phi   = prod_phi * phia[j];
        prod_q     = prod_q   * q[j - 1];
        pr_a[t, j] = prod_phi * prod_q * p[j];
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
  array[nsites, nyears] int y;         // counts [site, year]
  array[nsites, nyears] int J;         // juveniles caught [site, year]
  array[nsites, nyears] int R;         // adults monitored [site, year]
  array[nsites, nyears - 1, nyears] int marray_a;
}

transformed data {
  int ny_minus_1 = nyears - 1;
}

parameters {
  // Latent state variables — site × year
  array[nsites] vector<lower=0>[nyears] N1;
  array[nsites] vector<lower=0>[nyears] NadSurv;
  array[nsites] vector<lower=0>[nyears] Nadimm;

  // Grand means (logit / log scale)
  real l_mphia;
  real l_mfec;
  real l_mim;
  real l_p;

  // Year random effects (non-centred)
  vector[ny_minus_1] epsilon_phia_raw;
  vector[ny_minus_1] epsilon_fec_raw;
  vector[ny_minus_1] epsilon_im_raw;

  // Site random effects (non-centred)
  vector[nsites] zeta_phia_raw;
  vector[nsites] zeta_fec_raw;
  vector[nsites] zeta_im_raw;

  // Site×year interaction random effects (non-centred)
  array[nsites] vector[ny_minus_1] eta_phia_raw;
  array[nsites] vector[ny_minus_1] eta_fec_raw;
  array[nsites] vector[ny_minus_1] eta_im_raw;

  // Standard deviations — year level
  real<lower=0> sig_phia;
  real<lower=0> sig_fec;
  real<lower=0> sig_im;

  // Standard deviations — site level
  real<lower=0> sig_site_phia;
  real<lower=0> sig_site_fec;
  real<lower=0> sig_site_im;

  // Standard deviations — site×year level
  real<lower=0> sig_sy_phia;
  real<lower=0> sig_sy_fec;
  real<lower=0> sig_sy_im;
}

transformed parameters {
  // Scale year deviations 
  vector[ny_minus_1] epsilon_phia = sig_phia * epsilon_phia_raw;
  vector[ny_minus_1] epsilon_fec  = sig_fec  * epsilon_fec_raw;
  vector[ny_minus_1] epsilon_im   = sig_im   * epsilon_im_raw;

  // Scale site deviations
  vector[nsites] zeta_phia = sig_site_phia * zeta_phia_raw;
  vector[nsites] zeta_fec  = sig_site_fec  * zeta_fec_raw;
  vector[nsites] zeta_im   = sig_site_im   * zeta_im_raw;

  //Site×year interaction deviations
  array[nsites] vector[ny_minus_1] eta_phia;
  array[nsites] vector[ny_minus_1] eta_fec;
  array[nsites] vector[ny_minus_1] eta_im;

  //Demographic rates [site, year]
  array[nsites] vector<lower=0, upper=1>[ny_minus_1] phia;
  array[nsites] vector<lower=0>[ny_minus_1] f;
  array[nsites] vector<lower=0>[ny_minus_1] omega;
  vector<lower=0, upper=1>[ny_minus_1] p;  // recapture: year-only

  //Derived measures
  array[nsites] vector<lower=0>[nyears]      Ntot;
  array[nsites] vector<lower=0>[ny_minus_1]  rho;
  array[nsites, ny_minus_1] simplex[nyears]  pr_a;

  // Recapture probability — no random effects
  for (t in 1:ny_minus_1)
    p[t] = inv_logit(l_p);

  for (s in 1:nsites) {
    // Scale interaction terms
    eta_phia[s] = sig_sy_phia * eta_phia_raw[s];
    eta_fec[s]  = sig_sy_fec  * eta_fec_raw[s];
    eta_im[s]   = sig_sy_im   * eta_im_raw[s];

    for (t in 1:ny_minus_1) {
      // Linear predictor: grand mean + year[t] + site[s] + site×year[s,t]
      phia[s, t]  = inv_logit(l_mphia
                               + epsilon_phia[t]
                               + zeta_phia[s]
                               + eta_phia[s][t]);
      f[s, t]     = exp(l_mfec
                        + epsilon_fec[t]
                        + zeta_fec[s]
                        + eta_fec[s][t]);
      omega[s, t] = exp(l_mim
                        + epsilon_im[t]
                        + zeta_im[s]
                        + eta_im[s][t]);
    }

    // Total population
    Ntot[s] = NadSurv[s] + Nadimm[s] + N1[s];

    // m-array cell probabilities
    pr_a[s] = marray_adults(nyears, phia[s], p);

    // Expected juveniles caught: R[s,t] adults monitored × fecundity rate
    for (t in 1:ny_minus_1)
      rho[s, t] = R[s, t] * f[s, t];
  }
}

model {
  //Priors: initial population sizes
  for (s in 1:nsites) {
    N1[s, 1]      ~ normal(100, 100);
    NadSurv[s, 1] ~ normal(100, 100);
    Nadimm[s, 1]  ~ normal(100, 100);
  }

  //Priors: grand means
  l_mphia ~ normal(0, 100);
  l_mfec  ~ normal(0, 100);
  l_mim   ~ normal(0, 100);
  l_p     ~ normal(0, 100);

  //Priors: year random effects (non-centred)
  epsilon_phia_raw ~ normal(0, 1);
  epsilon_fec_raw  ~ normal(0, 1);
  epsilon_im_raw   ~ normal(0, 1);

  //Priors: site random effects (non-centred)
  zeta_phia_raw ~ normal(0, 1);
  zeta_fec_raw  ~ normal(0, 1);
  zeta_im_raw   ~ normal(0, 1);

  //Priors: site×year interaction random effects (non-centred) 
  for (s in 1:nsites) {
    eta_phia_raw[s] ~ normal(0, 1);
    eta_fec_raw[s]  ~ normal(0, 1);
    eta_im_raw[s]   ~ normal(0, 1);
  }

// Weakly informative half-normal priors on sigma parameters
sig_phia      ~ normal(0, 1);
sig_fec       ~ normal(0, 1);
sig_im        ~ normal(0, 1);
sig_site_phia ~ normal(0, 1);
sig_site_fec  ~ normal(0, 1);
sig_site_im   ~ normal(0, 1);
sig_sy_phia   ~ normal(0, 0.5);
sig_sy_fec    ~ normal(0, 0.5);
sig_sy_im     ~ normal(0, 0.5);

  //Likelihood
  for (s in 1:nsites) {
    for (t in 2:nyears) {
      real mean1 = 0.5 * f[s, t-1] * Ntot[s, t-1];
      real mpo   = Ntot[s, t-1] * omega[s, t-1];

      N1[s, t]      ~ real_poisson(mean1);
      NadSurv[s, t] ~ real_binomial(Ntot[s, t-1], phia[s, t-1]);
      Nadimm[s, t]  ~ real_poisson(mpo);
    }
  }

  //Likelihood (counts)
  for (s in 1:nsites)
    y[s] ~ poisson(Ntot[s]);

  //Likelihood (CJS)
  for (s in 1:nsites) {
    for (t in 1:ny_minus_1) {
      marray_a[s, t] ~ multinomial(pr_a[s, t]);
    }
  }

  //Likelihood (productivity)
  for (s in 1:nsites)
    J[s] ~ poisson(rho[s]);
}

generated quantities {
  //Grand mean demographic rates (back-transformed)
  real<lower=0, upper=1> mphia = inv_logit(l_mphia);
  real<lower=0>          mfec  = exp(l_mfec);
  real<lower=0>          mim   = exp(l_mim);

  //Population growth rates per site
  array[nsites] vector<lower=0>[ny_minus_1] lambda;
  array[nsites] vector[ny_minus_1]          logla;
  array[nsites] real<lower=0>               mlam;

  for (s in 1:nsites) {
    lambda[s] = Ntot[s, 2:nyears] ./ Ntot[s, 1:ny_minus_1];
    logla[s]  = log(lambda[s]);
    mlam[s]   = exp(1.0 / ny_minus_1 * sum(logla[s]));
  }
  
  real var_phia_year       = square(sig_phia);
  real var_phia_siteyear   = square(sig_sy_phia);
  real icc_phia = var_phia_year / (var_phia_year + var_phia_siteyear);
  
  real var_fec_year       = square(sig_fec);
  real var_fec_siteyear   = square(sig_sy_fec);
  real icc_fec = var_fec_year / (var_fec_year + var_fec_siteyear);
  
  real var_im_year       = square(sig_im);
  real var_im_siteyear   = square(sig_sy_im);
  real icc_im = var_im_year / (var_im_year + var_im_siteyear);
}
