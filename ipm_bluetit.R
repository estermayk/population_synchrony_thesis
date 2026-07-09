library(rstan)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(viridis)
library(ClustGeo)
library(sf)

# Exercise3: Transform blue tit data into similar format used here for hoopoe dataset, namely:


# a) Collapse individual capture histories to m-array format for each site (code to convert individual capture histories to m-arrays can be found in Chapter 7.10 in book, "Bayesian Population Analysis using WinBUGS --- A Hierarchical Perspective" (2012) by Marc Kéry and Michael Schaub)

#from kery and schaub 2012:
# Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

#lets test with ALN first

ALN_adults <- adults[(adults$site == "ALN"),]

ALN_adults_CH <- data.frame(individuals = c(unique(ALN_adults$ring)))

ALN_years <- sort(unique(ALN_adults$year))

for (year in ALN_years) {
  ALN_adults_CH[[as.character(year)]] <- NA
}

for (i in 1:nrow(ALN_adults_CH)) {
  individual <- ALN_adults_CH$individuals[i]
  capture_years <- unique(ALN_adults$year[ALN_adults$ring == individual])
  ALN_adults_CH[i, as.character(capture_years)] <- 1
}

ALN_adults_CH[is.na(ALN_adults_CH)] <- 0
ALN_adults_CH_matrix <- as.matrix(ALN_adults_CH[,-1])
ALN_adults_CH_matrix <- (apply(ALN_adults_CH_matrix, 2, as.numeric))

ALN_adults_CH_matrix

ALN_adults_marray <- marray(ALN_adults_CH_matrix)

ALN_adults_marray


adults <- adults %>%
  mutate(zone = case_when(site == "EDI" | site == "RSY" | site == "FOF" ~ "A",
                          site == "BAD" | site == "LVN" | site == "DOW" | site == "GLF" ~ "B",
                          site == "SER" | site == "MCH" | site == "PTH" | site == "STY" ~ "C",
                          site == "BIR" | site == "DUN" | site == "BLG"  ~ "D",
                          site == "PIT" | site == "KCZ" | site == "KCK" | site == "BLA" | site == "CAL" ~ "E",
                          site == "DNM" | site == "DNC" | site == "DNS" | site == "DLW"  ~ "F",
                          site == "CRU" | site == "NEW" | site == "HWP" | site == "INS" | site == "FSH" ~ "G",
                          site == "RTH" | site == "AVI" | site == "AVN" | site == "CAR" ~ "H",
                          site == "SLS" | site == "TOM" | site == "DAV" ~ "I",
                          site == "ART" | site == "MUN"  ~ "J",
                          site == "FOU" | site == "ALN" | site == "DEL" ~ "K",
                          site == "TAI" | site == "SPD" | site == "OSP" | site == "DOR" ~ "L"))


sitedat <- sitedat %>%
  mutate(zone = case_when(site == "EDI" | site == "RSY" | site == "FOF" ~ "A",
                          site == "BAD" | site == "LVN" | site == "DOW" | site == "GLF" ~ "B",
                          site == "SER" | site == "MCH" | site == "PTH" | site == "STY" ~ "C",
                          site == "BIR" | site == "DUN" | site == "BLG"  ~ "D",
                          site == "PIT" | site == "KCZ" | site == "KCK" | site == "BLA" | site == "CAL" ~ "E",
                          site == "DNM" | site == "DNC" | site == "DNS" | site == "DLW"  ~ "F",
                          site == "CRU" | site == "NEW" | site == "HWP" | site == "INS" | site == "FSH" ~ "G",
                          site == "RTH" | site == "AVI" | site == "AVN" | site == "CAR" ~ "H",
                          site == "SLS" | site == "TOM" | site == "DAV" ~ "I",
                          site == "ART" | site == "MUN"  ~ "J",
                          site == "FOU" | site == "ALN" | site == "DEL" ~ "K",
                          site == "TAI" | site == "SPD" | site == "OSP" | site == "DOR" ~ "L"))


phendat <- phendat %>%
  mutate(zone = case_when(site == "EDI" | site == "RSY" | site == "FOF" ~ "A",
                          site == "BAD" | site == "LVN" | site == "DOW" | site == "GLF" ~ "B",
                          site == "SER" | site == "MCH" | site == "PTH" | site == "STY" ~ "C",
                          site == "BIR" | site == "DUN" | site == "BLG"  ~ "D",
                          site == "PIT" | site == "KCZ" | site == "KCK" | site == "BLA" | site == "CAL" ~ "E",
                          site == "DNM" | site == "DNC" | site == "DNS" | site == "DLW"  ~ "F",
                          site == "CRU" | site == "NEW" | site == "HWP" | site == "INS" | site == "FSH" ~ "G",
                          site == "RTH" | site == "AVI" | site == "AVN" | site == "CAR" ~ "H",
                          site == "SLS" | site == "TOM" | site == "DAV" ~ "I",
                          site == "ART" | site == "MUN"  ~ "J",
                          site == "FOU" | site == "ALN" | site == "DEL" ~ "K",
                          site == "TAI" | site == "SPD" | site == "OSP" | site == "DOR" ~ "L"))


site_codes


#making a function to iterate over the sites
create_marray <- function(site_code, data) {
  site_data <- data[data$zone == site_code,]
  
  capture_history <- data.frame(individuals = unique(site_data$ring))
  site_years <- 2014:2025  
  #change the above line to sort(unique(site_data$year))
  for (year in site_years) {
    capture_history[[as.character(year)]] <- NA
  }
  
  for (i in 1:nrow(capture_history)) {
    individual <- capture_history$individuals[i]
    capture_years <- unique(site_data$year[site_data$ring == individual])
    capture_history[i, as.character(capture_years)] <- 1
  }
  
  capture_history[is.na(capture_history)] <- 0
  
  capture_history_matrix <- as.matrix(capture_history[,-1])
  capture_history_matrix <- apply(capture_history_matrix, 2, as.numeric)
  
  marray_result <- marray(capture_history_matrix)
  return(marray_result)
}


site_codes <- unique(adults$zone)

marray_list <- lapply(site_codes, function(zone) create_marray(zone, adults))

names(marray_list) <- site_codes
dim(marray_list)

nsites_bt <- 12
nyears_bt <- 12

marray_list <- array(unlist(marray_list), dim = c(nsites_bt, nyears_bt - 1, nyears_bt))

dim(marray_list)

# b) Also for each site, obtain observed population count (maximal number of simultaneously occupied nest boxes in each year), number of offspring and number of surveyed broods and investigate temporal trends.

all_combinations <- expand.grid(zone = site_codes, year = 2014:2025)

#max occupancy
phendat$ID <- as.factor(paste(phendat$year, phendat$site, phendat$box, sep = "_")) 
phendat$suc[phendat$suc == "-999"] <- 0
phendat <- phendat[!duplicated(phendat$ID, fromLast = TRUE), ]

ord_dates_phendat <- select(phendat, c('year', 'site', 'zone', 'ID', 'n1', 'nl', 'latestfed', 'latestcc', 'fed', 'cc', 'cs', 'fki', 'hatching_first_recorded', 'v1date', 'v2date'))


transform_dates <- function(df) {
  date_columns <- 70:180
  occupancy_df <- matrix(0, ncol = length(date_columns), nrow = nrow(df))
  colnames(occupancy_df) <- as.character(date_columns)
  
  for (i in 1:nrow(df)) {
    dates <- suppressWarnings(as.numeric(df[i, c('n1', 'nl', 'latestfed', 'latestcc', 'fed', 'cc', 
                                                 'cs', 'fki', 'hatching_first_recorded', 'v1date', 
                                                 'v2date')]))    
    date_range <- dates[!is.na(dates)]
    if (length(date_range) > 0) {
      min_date <- max(min(date_range), 70)
      max_date <- min(max(date_range), 180)
      valid_dates <- min_date:max_date
      valid_dates <- valid_dates[valid_dates %in% date_columns]
      occupancy_df[i, as.character(valid_dates)] <- 1
    }
  }
  
  occupancy_df <- as.data.frame(occupancy_df)
  occupancy_df <- cbind(df[, c('year', 'zone', 'ID')], occupancy_df)
  
  return(occupancy_df)
}

occupancy_data <- transform_dates(ord_dates_phendat)

max_occupancy <- function(df) {
  numeric_columns <- as.character(70:180)
  df[numeric_columns] <- lapply(df[numeric_columns], as.numeric)
  results <- df %>%
    group_by(year, zone) %>%
    summarise(max_occupancy = max(colSums(across(all_of(numeric_columns)), na.rm = TRUE), na.rm = TRUE)) %>%
    ungroup()
  
  return(results)
}

max_occupancy_df <- max_occupancy(occupancy_data)

#max_occupancy_df_lat <- max_occupancy_df %>%
  #left_join(sitedat %>% select(site, Mean.Lat), by = "site")

occupancy_p <- ggplot(max_occupancy_df, aes(x = factor(year), y = max_occupancy, colour = zone, fill = zone)) +
  #scale_colour_gradient(high = 'purple3', low = 'orange') +
  #scale_fill_gradient(high = 'purple3', low = 'orange') +
  geom_violin(fill = "lightgray", alpha = 0.3, color = NA) +
  geom_path(
    aes(group = zone), 
    size = 0.4, alpha = 0.5, 
    position = position_jitter(width = 0.1, seed = 3922)
  ) + 
  #geom_point(size = 1, alpha = 0.7, position = position_jitter(width = 0.1, seed = 3922)) +  
  theme_minimal() +
  scale_y_continuous(n.breaks = 9) +
  labs(title = "Annual Occupancy by Site",
       x = "Year",
       y = "Occupancy",
       fill = "Mean Latitude",
       colour = 'Mean Latitude') +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_line(color = "gray", size = 0.5)
  )

max_occupancy_df <- all_combinations %>%
  left_join(max_occupancy_df, by = c("zone", "year")) %>%
  replace_na(list(max_occupancy = 0))

#number of offspring

count_fledglings <- function(df) {
  results <- df %>%
    group_by(year, zone) %>%
    summarise(fledgling_count = sum(suc, na.rm = TRUE)) %>%
    ungroup()
  
  return(results)
}

n_offspring_df <- count_fledglings(phendat)

#n_offspring_df_lat <- n_offspring_df %>%
 # left_join(sitedat %>% select(site, Mean.Lat), by = "site")

n_fledge_p <- ggplot(n_offspring_df, aes(x = factor(year), y = fledgling_count, fill = zone, colour = zone)) +
  #scale_colour_gradient(high = 'purple3', low = 'orange') +
  #scale_fill_gradient(high = 'purple3', low = 'orange') +
  geom_violin(fill = "lightgray", alpha = 0.3, color = NA) +
  geom_path(
    aes(group = zone), 
    size = 0.4, alpha = 0.5, 
    position = position_jitter(width = 0.1, seed = 3922)
  ) + 
  geom_point(size = 1, alpha = 0.7, position = position_jitter(width = 0.1, seed = 3922)) +  
  #geom_line(aes(group = site, colour = Mean.Lat)) +
  theme_minimal() +
  labs(title = "Annual fledgling count by Site",
       x = "Year",
       y = "N fledglings",
       fill = "Mean Latitude",
       colour = "Mean Latitude") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_line(color = "gray", size = 0.5)
  )

n_offspring_df <- all_combinations %>%
  left_join(n_offspring_df, by = c("zone", "year")) %>%
  replace_na(list(fledgling_count = 0))

count_broods_surveyed <- function(df) {
  numeric_columns <- as.character(70:180)
  df[numeric_columns] <- lapply(df[numeric_columns], as.numeric)
  results <- df %>%
    group_by(year, zone) %>%
    summarise(brood_count = sum(rowSums(across(all_of(numeric_columns))) > 0, na.rm = TRUE)) %>%    ungroup()
  
  return(results)
}

broods_surveyed_df <- count_broods_surveyed(occupancy_data)

view(occupancy_data)

#broods_surveyed_df_lat <- broods_surveyed_df %>%
 # left_join(sitedat %>% select(site, Mean.Lat), by = "site")

broods_surveyed_df <- broods_surveyed_df %>%
  group_by(year, zone, brood_count) %>%
  summarise(count = n()) %>%
  ungroup()


broods_p <- ggplot(broods_surveyed_df, aes(x = factor(year), y = brood_count, colour = zone, fill = zone)) +
  #scale_colour_gradient(high = 'purple3', low = 'orange') +
  #scale_fill_gradient(high = 'purple3', low = 'orange') +
  geom_violin(fill = "darkgray", alpha = 0.3, color = NA) +
  geom_path(
    aes(group = zone), 
    size = 0.5, 
    position = position_jitter(width = 0.1, height = 0, seed = 3922)
  ) +
#  geom_path(
#    aes(group = zone, linewidth = count),   # Use 'count' for line thickness
#    colour = "black"
#  ) +
  geom_point(size = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0, seed = 3922)) +  
  theme_minimal() +
  labs(title = "Annual broods surveyed by Site",
       x = "Year",
       y = "N broods") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_line(color = "gray", size = 0.5)
  )

broods_surveyed_df <- all_combinations %>%
  left_join(broods_surveyed_df, by = c("zone", "year")) %>%
  replace_na(list(brood_count = 0))

for (i in site_codes) {
  y_list[[i]] <- max_occupancy_df %>% filter(zone == i) %>% arrange(year) %>% pull(max_occupancy)
  J_list[[i]] <- n_offspring_df   %>% filter(zone == i) %>% arrange(year) %>% pull(fledgling_count)
  R_list[[i]] <- broods_surveyed_df %>% filter(zone == i) %>% arrange(year) %>% pull(brood_count)
}

(occupancy_p / broods_p / n_fledge_p)

# Exercise 4: Once blue tit data are in format similar to here, run the Stan model on blue tit dataset and evaluate.

list_to_mat <- function(lst, site_codes, nyears) {
  mat <- matrix(0, nrow = length(site_codes), ncol = nyears,
                dimnames = list(site_codes, NULL))
  for (i in seq_along(site_codes)) {
    vals <- unlist(lst[[site_codes[i]]])
    mat[i, seq_along(vals)] <- vals
  }
  mat
}

y_mat_bt <- list_to_mat(y_list, site_codes, nyears_bt)
J_mat_bt <- list_to_mat(J_list, site_codes, nyears_bt)
R_mat_bt <- list_to_mat(R_list, site_codes, nyears_bt)

#Y - pop counts (max occ)
#J - nestling counts
#R - surveyed broods

print(dim(y_mat_bt))
print(dim(J_mat_bt))
print(dim(R_mat_bt))
print(dim(marray_list))

view(marray_list)

y_mat_bt    <- apply(y_mat_bt,    2, as.integer)
J_mat_bt <- apply(J_mat_bt[,1:11], 2, as.integer)
R_mat_bt <- apply(R_mat_bt[,1:11], 2, as.integer)

marray_list <- vector("list", length(site_codes))
names(marray_list) <- site_codes

for (i in seq_along(site_codes)) {
  ma <- create_marray(site_codes[i], adults)
  cat("Zone", site_codes[i], "- dim:", dim(ma), "\n")  # verify each one
  marray_list[[i]] <- ma
}

marray_array <- array(NA_integer_, dim = c(nsites_bt, nyears_bt - 1, nyears_bt))
for (i in 1:nsites_bt) {
  marray_array[i,,] <- marray_list[[i]]
}

marray_list <- apply(marray_array, c(1,2,3), as.integer)

for(s in 1:nsites_bt) {
  row_sums <- apply(marray_list[s,,], 1, sum)
  if(any(row_sums == 0)) cat("Zone", s, "has zero m-array rows:", which(row_sums == 0), "\n")
}

print(dim(y_mat_bt))
print(dim(J_mat_bt))
print(dim(R_mat_bt))
print(dim(marray_list))
nyears_bt
nsites_bt

stan_data_bt <- list(
  nyears    = nyears_bt,
  nsites    = nsites_bt,
  y         = y_mat_bt,      
  J         = J_mat_bt,      
  R         = R_mat_bt,        
  marray_a  = marray_list
)

params_bt <- c(
  "phia", "omega", "p",
  "mphia", "mim",
  "lambda", "mlam",
  "sig_site_phia", "sig_site_im",
  "sig_sy_phia",   "sig_sy_im",
  "N1", "NadSurv", "Nadimm", "Ntot",
  "var_phia_year", "var_phia_siteyear", "icc_phia",
  "var_im_year", "var_im_siteyear", "icc_im"
)

 # Give Stan parameters to estimate
 params_bt_s <- c(
   # Demographic rates and recapture prob [site, year]
   "phij", "phia", "f", "omega", "p",
   
   # Grand mean demographic rates and recapture prob (back-transformed)
   "mphij", "mphia", "mprod", "mim", "mp",
   
   # Population growth rates
   "lambda", "mlam",
   
   # Standard deviations: year level
   "sig_year_phia", "sig_year_prod", "sig_year_im", "sig_year_p",
   
   # Standard deviations: site level
   "sig_site_phia", "sig_site_prod", "sig_site_im", "sig_site_p",
   
   # Standard deviations: site×year level
   "sig_sy_phia", "sig_sy_prod", "sig_sy_im", "sig_sy_p",
   
   # Latent  state variables
   "N1", "NadSurv", "Nadimm", "Ntot",
   
   # Variance components and ICC
   "var_phia_year", "var_phia_site", "var_phia_siteyear", "icc_phia",
   "var_prod_year",  "var_prod_site", "var_prod_siteyear",  "icc_prod",
   "var_im_year",   "var_im_site", "var_im_siteyear",   "icc_im",
   "var_p_year", "var_p_site", "var_p_siteyear", "icc_p"
 )


# MCMC settings
ni <- 20000   
nb <- 10000   
nt <- 5       
nc <- 4

inits_bt <- lapply(1:nc, function(i) {
     list(
       # Grand means (on logit scale; 
       # to convert to probability scale to understand what the numbers mean, run plogis(-1.1), plogis(0) etc)
       l_mphij = -1.1,   
       l_mphia =  0,
       l_mprod  =  0,
       l_mim   = 0,
       l_p     =  0,
       
       # sigmas (std dev)
       sig_phia      = 1,
       sig_prod       = 1,
       sig_im        = 1,
       sig_p         = 1,
       sig_site_phia = 1,
       sig_site_prod  = 1,
       sig_site_im   = 1,
       sig_site_p    = 1,
       sig_sy_phia   = 1,
       sig_sy_prod    = 1,
       sig_sy_im     = 1,
       sig_sy_p      = 1,
       
       # Raw random effects 
       epsilon_phia_raw = rep(0, nyears_bt - 1),
       epsilon_prod_raw  = rep(0, nyears_bt - 1),
       epsilon_im_raw   = rep(0, nyears_bt - 1),
       epsilon_p_raw = rep(0, nyears_bt - 1),
       zeta_phia_raw    = rep(0, nsites_bt),
       zeta_prod_raw     = rep(0, nsites_bt),
       zeta_im_raw      = rep(0, nsites_bt),
       zeta_p_raw    = rep(0, nsites_bt),
       eta_phia_raw     = lapply(1:nsites_bt, function(s) rep(0, nyears_bt - 1)),
       eta_prod_raw      = lapply(1:nsites_bt, function(s) rep(0, nyears_bt - 1)),
       eta_im_raw       = lapply(1:nsites_bt, function(s) rep(0, nyears_bt - 1)),
       eta_p_raw     = lapply(1:nsites_bt, function(s) rep(0, nyears_bt - 1)),
       
       # Inits for pop sizes in the components of Ntot
       N1      = lapply(1:nsites_bt, function(s) rep(40, nyears_bt)),
       NadSurv = lapply(1:nsites_bt, function(s) rep(40, nyears_bt)),
       Nadimm  = lapply(1:nsites_bt, function(s) rep(40, nyears_bt))
     )
   })
  
#run N1 etc separately 
N2 <- lapply(1:nsites_bt, function(s) round(runif(nyears,  0, 3)))
N2
#tweak to edit and get a more reasonable distribution 


closeAllConnections()

# Setting stan model options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)


ipm_bt_debug <- stan(
  file   = "ipm_bluti_final.stan",
  data   = stan_data_bt,
  #init   = inits_bt[1],   
  chains = 1,
  iter   = 10
)

?stan

ipm_bt <- stan(
  file    = "ipm_bluti_final.stan",
  data    = stan_data_bt,
  init    = inits_bt,
  pars    = params_bt_s,
  chains  = nc,
  iter    = ni,
  warmup  = nb,
  thin    = nt,
  seed    = 2,
  control = list(
    adapt_delta   = 0.9,
    max_treedepth = 10    
  )
)

ipm_bt_test <- stan(
  file    = "ipm_bluti_final.stan",
  data    = stan_data_bt,
  init    = inits_bt,
  pars    = params_bt,
  chains  = 2,
  iter    = 500,
  warmup  = 250,
  thin    = 1,
  seed    = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 10)
)

print(ipm_bt,
      pars = c("mphia", "mim", "mprod",
               "sig_site_im", "sig_site_prod",
               "sig_sy_phia",   "sig_sy_im", "sig_sy_prod",
               "mlam",
               "var_phia_year", "var_phia_siteyear", "icc_phia",
               "var_im_year", "var_im_siteyear", "icc_im",
               "var_prod_year", "var_prod_siteyear", "icc_prod"),
      digits = 3)

# Get posterior distributions
posterior_ipm <- as.data.frame(ipm_bt)

posterior_ipm_test <- as.data.frame(ipm_bt_test)

rds_file_path <- "ipm_bt.rds"

saveRDS(ipm_bt, "~/population_synchrony_thesis/ipm_bt_v2.rds")

ipm_bt <- readRDS("ipm_bt.rds")

view(posterior_ipm)

