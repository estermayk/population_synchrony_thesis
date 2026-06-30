# Dataset: Hoopoe dataset in the SW Swiss Alps (Valais) studied from 2002 to 2009 
# (More details of population in Chapter 11 of book, "Bayesian Population Analysis using WinBUGS - A Hierarchical Perspective" (2012) by Marc Kéry and Michael Schaub and also published in Arlettaz et al., 2010, Schaub et al., 2011).

# Description of Integrated population model from Kery and Schaub 2012 book 
# - Age structured model with 2 age classes:
# 1-year old and adults (at least 2-years old)
# - Age at first breeding = 1 year
# - Prebreeding census, female-based
# - All vital rates are assumed to be time-dependent (random) (NB: I have now added site and site:year interaction as random effects for all demographic parameters, also estimating spatial synchrony ICCs)
# - Explicit estimation of immigration
# - Data provided to IPM are:
#   A) Recapture histories (specified as m-arrays. More details on m-arrays and code to convert individual capture histories to m-arrays can be found in Chapter 7.10 in book, "Bayesian Population Analysis using WinBUGS --- A Hierarchical Perspective" (2012) by Marc Kéry and Michael Schaub.)
# B) Population count data (the maximal number of simultaneously occupied nest boxes in each year for each site)
# C) Number of offspring (per year per site)
# D) Number of surveyed broods (per year per site)
# NB: More detailed description of IPM code available in Chapter 11, Kery & Schaub 2012 book

# Further details on immigration:
# How is immigration estimated? From Chapter 11, Kery & Schaub 2012 book: "By modeling apparent survival in a CJS model as a part of the integrated population model,
# we automatically account for emigration, even if we cannot estimate it. On the other hand, information
# on immigration is very elusive. In our study, the population counts contain information about
# immigration because the annual change in population size is a result of all four demographic processes
# operating in a population: survival, productivity, immigration, and emigration. We have independent
# data for the other demographic processes (capture–recapture for apparent survival, reproductive output
# for productivity); hence, the combined analysis should enable us to obtain an estimate of immigration.
# This neat idea was brought up by Abadi et al. (2010b)."

# BUGS code for IPM from book "Bayesian Population Analysis using WinBUGS - A Hierarchical Perspective" (2012) by Marc Kéry and Michael Schaub 
# was translated to Stan by Hiroki Itô (https://researchmap.jp/read0208767/) and is available at https://github.com/stan-dev/example-models/tree/master/BPA) .
# Hiroki Itô's Stan model code only fit year as a random effect for all demographic parameters (1-year old and 2+ years old survival, productivity, immigration) estimated in IPM.
# I simulated 2 more datasets based on original hoopoe dataset (site 1) and also updated Hiroki Itô's Stan model to include site and site:year as random effects for juvenile and adult survival, productivity and immigration, also estimating spatial synchrony ICCs

# Load libraries
library(rstan)
library(tidyverse)
library(ggplot2)
library(patchwork)

# Simulate datasets for two more sites
nsites <- 3
nyears <- 9

# Note: CJS model is specified as a multinomial likelihood (i.e. individual capture histories are collapsed into m-arrays) and were not specified as individual-level state-space formulation I had done previously for blue tits
# More details on m-arrays and code to convert individual capture histories to m-arrays can be found in Chapter 7.10 in book, "Bayesian Population Analysis using WinBUGS --- A Hierarchical Perspective" (2012) by Marc Kéry and Michael Schaub.

# All site-1 associated data is the hoopoe dataset from Chapter 11, Kery & Schaub 2012 book
# Capture-recapture data: m-array of juveniles and adults (these are males and females together)
marray_j_site1 <- structure(
  c(15, 0, 0, 0, 0, 0, 0, 0, 3, 34, 0, 0, 0, 0, 0, 0,
    0, 9, 56, 0, 0, 0, 0, 0, 0, 1, 8, 48, 0, 0, 0, 0, 0, 0, 1, 3,
    45, 0, 0, 0, 0, 0, 0, 1, 13, 27, 0, 0, 0, 0, 0, 0, 2, 7, 37,
    0, 0, 0, 0, 0, 0, 0, 3, 39, 198, 287, 455, 518, 463, 493, 434, 405),
  .Dim = 8:9)

marray_a_site1 <- structure(
  c(14, 0, 0, 0, 0, 0, 0, 0, 2, 22, 0, 0, 0, 0, 0, 0,
    0, 4, 34, 0, 0, 0, 0, 0, 0, 0, 2, 51, 0, 0, 0, 0, 0, 0, 0, 3,
    45, 0, 0, 0, 0, 0, 0, 0, 3, 44, 0, 0, 0, 0, 0, 0, 0, 3, 48, 0,
    0, 0, 0, 0, 0, 0, 2, 51, 43, 44, 79, 94, 118, 113, 99, 90),
  .Dim = 8:9)

y_site1 <- c(32, 42, 64, 85, 82, 78, 73, 69, 79) # Population count data  (the maximal number of simultaneously occupied nest boxes in each year)
J_site1 <- c(189, 274, 398, 538, 520, 476, 463, 438) # Number of offspring
R_site1 <- c(28, 36, 57, 77, 81, 83, 77, 72) # Number of surveyed broods

# Site 2 
set.seed(42)
y_site2 <- round(y_site1 * 0.60 * exp(rnorm(nyears, 0, 0.10)))
J_site2 <- round(J_site1 * 0.60 * exp(rnorm(nyears - 1, 0, 0.10)))
R_site2 <- round(R_site1 * 0.60 * exp(rnorm(nyears - 1, 0, 0.10)))

view(y_site2)

marray_j_site2 <- structure(
  c(12, 0, 0, 0, 0, 0, 0, 0, 2, 28, 0, 0, 0, 0, 0, 0,
    0, 7, 45, 0, 0, 0, 0, 0, 0, 2, 6, 38, 0, 0, 0, 0, 0, 0, 2, 4,
    36, 0, 0, 0, 0, 0, 0, 1, 10, 22, 0, 0, 0, 0, 0, 0, 1, 5, 29,
    0, 0, 0, 0, 0, 0, 0, 2, 31, 210, 265, 430, 490, 440, 470, 410, 385),
  .Dim = 8:9)

marray_a_site2 <- structure(
  c(11, 0, 0, 0, 0, 0, 0, 0, 3, 19, 0, 0, 0, 0, 0, 0,
    0, 3, 29, 0, 0, 0, 0, 0, 0, 0, 3, 44, 0, 0, 0, 0, 0, 0, 0, 2,
    38, 0, 0, 0, 0, 0, 0, 0, 2, 38, 0, 0, 0, 0, 0, 0, 0, 2, 41, 0,
    0, 0, 0, 0, 0, 0, 1, 44, 40, 38, 68, 82, 105, 100, 88, 78),
  .Dim = 8:9)

# Site 3 
set.seed(43)
y_site3 <- round(y_site1 * 1.20 * exp(rnorm(nyears, 0, 0.12)))
J_site3 <- round(J_site1 * 1.20 * exp(rnorm(nyears - 1, 0, 0.12)))
R_site3 <- round(R_site1 * 1.20 * exp(rnorm(nyears - 1, 0, 0.12)))

marray_j_site3 <- structure(
  c(18, 0, 0, 0, 0, 0, 0, 0, 4, 40, 0, 0, 0, 0, 0, 0,
    0, 11, 63, 0, 0, 0, 0, 0, 0, 1, 10, 55, 0, 0, 0, 0, 0, 0, 1, 4,
    52, 0, 0, 0, 0, 0, 0, 2, 15, 31, 0, 0, 0, 0, 0, 0, 2, 9, 43,
    0, 0, 0, 0, 0, 0, 0, 4, 45, 185, 305, 478, 540, 480, 510, 452, 420),
  .Dim = 8:9)

marray_a_site3 <- structure(
  c(16, 0, 0, 0, 0, 0, 0, 0, 2, 25, 0, 0, 0, 0, 0, 0,
    0, 5, 39, 0, 0, 0, 0, 0, 0, 0, 2, 58, 0, 0, 0, 0, 0, 0, 0, 4,
    51, 0, 0, 0, 0, 0, 0, 0, 4, 50, 0, 0, 0, 0, 0, 0, 0, 4, 54, 0,
    0, 0, 0, 0, 0, 0, 3, 57, 46, 50, 88, 102, 128, 120, 108, 97),
  .Dim = 8:9)

# Combine data from all sites
# Combining pop counts
(y_mat <- rbind(y_site1, y_site2, y_site3))

# Combining offpsring numbers (J) and surveyed broods (R)
# J, R: [nsites, nyears-1] format
(J_mat <- rbind(J_site1, J_site2, J_site3))
(R_mat <- rbind(R_site1, R_site2, R_site3))

# Combining recapture histories (m-array format)
# marray: [nsites, nyears-1, nyears] format
marray_j_3d <- array(NA_integer_, dim = c(nsites, nyears - 1, nyears))
marray_a_3d <- array(NA_integer_, dim = c(nsites, nyears - 1, nyears))

marray_j_3d[1,,] <- marray_j_site1
marray_j_3d[2,,] <- marray_j_site2
marray_j_3d[3,,] <- marray_j_site3

marray_a_3d[1,,] <- marray_a_site1
marray_a_3d[2,,] <- marray_a_site2
marray_a_3d[3,,] <- marray_a_site3

view(marray_a_3d)

# Put together all data to feed into stan model
stan_data <- list(
  nyears    = nyears,
  nsites    = nsites,
  y         = y_mat,      
  J         = J_mat,      
  R         = R_mat,        
  marray_j  = marray_j_3d, 
  marray_a  = marray_a_3d  
)

# Give Stan parameters to estimate
params <- c(
  "phij", "phia", "f", "omega", "p",
  "mphij", "mphia", "mfec", "mim",
  "lambda", "mlam",
  "sig_phij",      "sig_phia",      "sig_fec",      "sig_im",
  "sig_site_phij", "sig_site_phia", "sig_site_fec", "sig_site_im",
  "sig_sy_phij",   "sig_sy_phia",   "sig_sy_fec",   "sig_sy_im",
  "N1", "NadSurv", "Nadimm", "Ntot",
  "var_phij_year", "var_phij_siteyear", "icc_phij",
  "var_phia_year", "var_phia_siteyear", "icc_phia",
  "var_fec_year", "var_fec_siteyear", "icc_fec",
  "var_im_year", "var_im_siteyear", "icc_im"
)

# Setting stan model options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

# MCMC settings
ni <- 20000   
nb <- 10000   
nt <- 5       
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i) {
  list(l_mphij = rnorm(1, 0.2, 0.5),
       l_mphia = rnorm(1, 0.2, 0.5),
       l_mfec = rnorm(1, 0.2, 0.5),
       l_mim = rnorm(1, 0.2, 0.5),
       l_p = rnorm(1, 0.2, 1),
       sig_phij = runif(1, 0.1, 10),
       sig_phia = runif(1, 0.1, 10),
       sig_fec = runif(1, 0.1, 10),
       sig_im = runif(1, 0.1, 10),
       sig_site_phij = runif(1, 0.1, 10),
       sig_site_phia = runif(1, 0.1, 10),
       sig_site_fec = runif(1, 0.1, 10),
       sig_site_im = runif(1, 0.1, 10),
       sig_sy_phij = runif(1, 0.1, 10),
       sig_sy_phia = runif(1, 0.1, 10),
       sig_sy_fec = runif(1, 0.1, 10),
       sig_sy_im = runif(1, 0.1, 10),
       N1      = lapply(1:nsites, function(s) round(runif(nyears,  1, 50))),
       NadSurv = lapply(1:nsites, function(s) round(runif(nyears,  5, 50))),
       Nadimm  = lapply(1:nsites, function(s) round(runif(nyears,  1, 50)))
)})

# Run stan model
ipm_hoopoe <- stan(
  file    = "ipm_hoopoe.stan",
  data    = stan_data,
  init    = inits,
  pars    = params,
  chains  = nc,
  iter    = ni,
  warmup  = nb,
  thin    = nt,
  seed    = 2,
  control = list(
    adapt_delta   = 0.99,
    max_treedepth = 15    
  )
)

rds_file_path <- "ipm_hoopoe.rds"

saveRDS(ipm_hoopoe, rds_file_path)

ipm_hoopoe <- readRDS("ipm_hoopoe.rds")

# Print stan model output
print(ipm_hoopoe,
      pars = c("mphij", "mphia", "mfec", "mim",
               "sig_phij",      "sig_phia",      "sig_fec",      "sig_im",
               "sig_site_phij", "sig_site_phia", "sig_site_fec", "sig_site_im",
               "sig_sy_phij",   "sig_sy_phia",   "sig_sy_fec",   "sig_sy_im",
               "mlam",   "var_phij_year", "var_phij_siteyear", "icc_phij",
               "var_phia_year", "var_phia_siteyear", "icc_phia",
               "var_fec_year", "var_fec_siteyear", "icc_fec",
               "var_im_year", "var_im_siteyear", "icc_im"),
      digits = 3)

# Get posterior distributions
posterior_ipm <- as.data.frame(ipm_hoopoe)

# extract summary (mean + 95% CrI) for parameters
extract_summary <- function(posterior, pattern) {
  cols <- grep(pattern, names(posterior), value = TRUE)
  data.frame(
    param = cols,
    mean  = colMeans(posterior[, cols, drop = FALSE]),
    lower = apply(posterior[, cols, drop = FALSE], 2, quantile, 0.025),
    upper = apply(posterior[, cols, drop = FALSE], 2, quantile, 0.975)
  )
}

# Get nyears and nsites
nyears  <- stan_data$nyears
nsites  <- stan_data$nsites
years   <- 1:nyears
years_m1 <- 1:(nyears - 1)  # for parameters indexed over nyears-1

site_labels <- paste0("Site ", 1:nsites)

# Estimated Ntot per site per year
ntot_list <- lapply(1:nsites, function(s) {
  pattern <- paste0("Ntot\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- years
  df$site <- site_labels[s]
  df
})
ntot_df <- bind_rows(ntot_list)

# Observed counts per site per year
obs_counts <- expand.grid(site = site_labels, year = years) %>%
  arrange(site, year) %>%
  mutate(observed = as.vector(t(stan_data$y)))  # y is [nsites, nyears]

# Merge both dfs
ntot_df <- left_join(ntot_df, obs_counts, by = c("site", "year"))

# Plot
p_ntot <- ggplot(ntot_df, aes(x = year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = mean, colour = "Estimated Ntot"), linewidth = 0.9) +
  geom_point(aes(y = mean, colour = "Estimated Ntot"), size = 1.5) +
  geom_line(aes(y = observed, colour = "Observed count"),
            linewidth = 0.9, linetype = "dashed") +
  geom_point(aes(y = observed, colour = "Observed count"), size = 1.5) +
  scale_colour_manual(values = c("Estimated Ntot" = "steelblue",
                                 "Observed count"  = "firebrick")) +
  facet_wrap(~ site, scales = "free_y") +
  labs(title = "Population size", x = "Year", y = "N",
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(n.breaks=8)

p_ntot

# Exercises for Ester -

print(dim(y_mat))
print(dim(J_mat))
print(dim(R_mat))
print(dim(marray_j_3d))
print(dim(marray_a_3d))

typeof(marray_a_3d)
typeof(marray_list)

# Exercise1: Plot and investigate temporal trends in other demographic parameters (survival, productivity, immigration, population growth rate etc) across the different sites (doing this will be handy for blue tit data analysis)

#productivity/fecundity 
# Estimated f per site per year
f_list <- lapply(1:nsites, function(s) {
  pattern <- paste0("f\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- years_m1
  df$site <- site_labels[s]
  df
})
f_est <- bind_rows(f_list)

# Observed f per site per year
obs_f <- expand.grid(site = site_labels, year = years_m1) %>%
  arrange(site, year) %>%
  mutate(observed = as.vector(t(stan_data$J / stan_data$R)))

# Merge both dfs
f_df <- left_join(f_est, obs_f, by = c("site", "year"))

# Plot
p_f <- ggplot(f_df, aes(x = year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = mean, colour = "Estimated f"), linewidth = 0.9) +
  geom_point(aes(y = mean, colour = "Estimated f"), size = 1.5) +
  geom_line(aes(y = observed, colour = "Observed f"),
            linewidth = 0.9, linetype = "dashed") +
  geom_point(aes(y = observed, colour = "Observed f"), size = 1.5) +
  scale_colour_manual(values = c("Estimated f" = "steelblue",
                                 "Observed f"  = "firebrick")) +
  facet_wrap(~ site, scales = "free_y") +
  labs(title = "Fecundity per Nest", x = "Year", y = "N Fledglings/female",
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(n.breaks=8)

p_f


#survival
#juveniles
# Estimated phij per site per year
phij_list <- lapply(1:nsites, function(s) {
  pattern <- paste0("phij\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- years_m1
  df$site <- site_labels[s]
  df
})
phij_est <- bind_rows(phij_list)

#adults
phia_list <- lapply(1:nsites, function(s) {
  pattern <- paste0("phia\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- years_m1
  df$site <- site_labels[s]
  df
})
phia_est <- bind_rows(phia_list)


# Observed f per site per year
#obs_phij <- expand.grid(site = site_labels, year = years_m1) %>%
 # arrange(site, year) %>%
  #mutate(observed = as.vector(t()))

# Merge both dfs
phi_df <- left_join(phij_est, phia_est, by = c("site", "year"))

# Plot
p_phi <- ggplot(phi_df, aes(x = year)) +
  geom_ribbon(aes(ymin = lower.x, ymax = upper.x), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = mean.x, colour = "Estimated juvenile survival"), linewidth = 0.9) +
  geom_point(aes(y = mean.x, colour = "Estimated juvenile survival"), size = 1.5) +
  geom_ribbon(aes(ymin = lower.y, ymax = upper.y), fill = "firebrick", alpha = 0.25) +
  geom_line(aes(y = mean.y, colour = "Estimated adult survival"), linewidth = 0.9) +
  geom_point(aes(y = mean.y, colour = "Estimated adult survival"), size = 1.5) +
  scale_colour_manual(values = c("Estimated juvenile survival" = "steelblue", "Estimated adult survival" = "firebrick")) +
  facet_wrap(~ site, scales = "free_y") +
  labs(title = "Estimated Survival", x = "Year", y = "Phi",
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(n.breaks=8)

p_phi

#population growth rate

lambda_list <- lapply(1:nsites, function(s) {
  pattern <- paste0("lambda\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- years_m1
  df$site <- site_labels[s]
  df
})
lambda_df <- bind_rows(lambda_list)

# Observed counts per site per year
#obs_counts <- expand.grid(site = site_labels, year = years) %>%
 # arrange(site, year) %>%
  #mutate(observed = as.vector(t(stan_data$y)))  # y is [nsites, nyears]

# Merge both dfs
#lambda_df <- left_join(lambda_df, obs_counts, by = c("site", "year"))

# Plot
p_lambda <- ggplot(lambda_df, aes(x = year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = mean, colour = "Estimated lambda"), linewidth = 0.9) +
  geom_point(aes(y = mean, colour = "Estimated lambda"), size = 1.5) +
  #geom_line(aes(y = observed, colour = "Observed count"),
   #         linewidth = 0.9, linetype = "dashed") +
  #geom_point(aes(y = observed, colour = "Observed count"), size = 1.5) +
  scale_colour_manual(values = c("Estimated lambda" = "steelblue")) +
  facet_wrap(~ site, scales = "free_y") +
  labs(title = "Population growth rate", x = "Year", y = "Lambda",
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(n.breaks=8)

p_lambda

#immigration
Nadimm_list <- lapply(1:nsites, function(s) {
  pattern <- paste0("Nadimm\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- years
  df$site <- site_labels[s]
  df
})
Nadimm_df <- bind_rows(Nadimm_list)

# Observed counts per site per year
obs_counts <- expand.grid(site = site_labels, year = years) %>%
  arrange(site, year) %>%
  mutate(observed = as.vector(t(stan_data$y)))  # y is [nsites, nyears]

# Merge both dfs
Nadimm_df <- left_join(Nadimm_df, obs_counts, by = c("site", "year"))

# Plot
p_Nadimm <- ggplot(Nadimm_df, aes(x = year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = mean, colour = "Estimated Nadimm"), linewidth = 0.9) +
  geom_point(aes(y = mean, colour = "Estimated Nadimm"), size = 1.5) +
#  geom_line(aes(y = observed, colour = "Observed count"),
#            linewidth = 0.9, linetype = "dashed") +
#  geom_point(aes(y = observed, colour = "Observed count"), size = 1.5) +
  scale_colour_manual(values = c("Estimated Nadimm" = "steelblue")) +
  facet_wrap(~ site, scales = "free_y") +
  labs(title = "N Immigrants", x = "Year", y = "N",
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(n.breaks=8)

p_Nadimm

rates_plots <- (p_ntot | p_lambda | p_f) / (p_phi | p_Nadimm)

rates_plots

ggsave("figs/rates_plot.png", plot = rates_plots)

# Exercise2: Plot and investigate posterior distributions of random effect variances and spatial synchrony (ICCs) for all demographic parameters (doing this will be handy for blue tit data analysis)

#lambda
lambda_pds <- posterior_ipm %>%
  select(contains("lambda")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

pds_lambda <- ggplot(lambda_pds, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Pooled Posterior Distribution of Lambda",
       x = "Lambda Value",
       y = "Density") +
  theme_minimal()

#phia icc

phia_icc <- posterior_ipm %>%
  select(contains("icc_phia")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

icc_phia_p <- ggplot(phia_icc, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Adult Survival ICC",
       x = "phia ICC",
       y = "Density") +
  theme_minimal()

icc_phia_p

#var_phia_siteyear

var_phia_siteyear <- posterior_ipm %>%
  select(contains("var_phia_siteyear")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_phia_siteyear_p <- ggplot(var_phia_siteyear, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Adult Survival Random Effect Variances",
       x = "phia var",
       y = "Density") +
  theme_minimal()

var_phia_siteyear_p

#phij icc

phij_icc <- posterior_ipm %>%
  select(contains("icc_phij")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

icc_phij_p <- ggplot(phij_icc, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Juvenile Survival ICC",
       x = "phij ICC",
       y = "Density") +
  theme_minimal()

icc_phij_p

#var_phij_siteyear

var_phij_siteyear <- posterior_ipm %>%
  select(contains("var_phij_siteyear")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_phij_siteyear_p <- ggplot(var_phij_siteyear, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Juvenile Survival Random Effect Variances",
       x = "phij var",
       y = "Density") +
  theme_minimal()

var_phij_siteyear_p

#fec icc

fec_icc <- posterior_ipm %>%
  select(contains("icc_fec")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

icc_fec_p <- ggplot(fec_icc, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Fecundity ICC",
       x = "fec ICC",
       y = "Density") +
  theme_minimal()

icc_fec_p

#var_fec_siteyear

var_fec_siteyear <- posterior_ipm %>%
  select(contains("var_fec_siteyear")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_fec_siteyear_p <- ggplot(var_fec_siteyear, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Fecundity Random Effect Variances",
       x = "fec var",
       y = "Density") +
  theme_minimal()

var_fec_siteyear_p

#im icc

im_icc <- posterior_ipm %>%
  select(contains("icc_im")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

icc_im_p <- ggplot(im_icc, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Immigration ICC",
       x = "im ICC",
       y = "Density") +
  theme_minimal()

icc_im_p

#var_im_siteyear

var_im_siteyear <- posterior_ipm %>%
  select(contains("var_im_siteyear")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_im_siteyear_p <- ggplot(var_im_siteyear, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Immigration Random Effect Variances",
       x = "im var",
       y = "Density") +
  theme_minimal()

var_im_siteyear_p

icc_var_plots <- (var_phia_siteyear_p | icc_phia_p) / (var_phij_siteyear_p | icc_phij_p) / (var_fec_siteyear_p | icc_fec_p) / (var_im_siteyear_p | icc_im_p)

icc_var_plots

ggsave("figs/icc_var_plots.png", plot = icc_var_plots)
