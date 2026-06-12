library(rstan)
library(tidyverse)
library(ggplot2)
library(patchwork)

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

#making a function to iterate over the sites
create_marray <- function(site_code, data) {
  site_data <- data[data$site == site_code,]
  
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

site_codes <- unique(adults$site)

marray_list <- lapply(site_codes, function(site) create_marray(site, adults))

names(marray_list) <- site_codes

# b) Also for each site, obtain observed population count (maximal number of simultaneously occupied nest boxes in each year), number of offspring and number of surveyed broods and investigate temporal trends.

#max occupancy
phendat$ID <- as.factor(paste(phendat$year, phendat$site, phendat$box, sep = "_")) 
phendat$suc[phendat$suc == "-999"] <- 0
ord_dates_phendat <- select(phendat, c('year', 'site', 'ID', 'n1', 'nl', 'latestfed', 'latestcc', 'fed', 'cc', 'cs', 'fki', 'hatching_first_recorded', 'v1date', 'v2date'))


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
  occupancy_df <- cbind(df[, c('year', 'site', 'ID')], occupancy_df)
  
  return(occupancy_df)
}

occupancy_data <- transform_dates(ord_dates_phendat)

max_occupancy <- function(df) {
  numeric_columns <- as.character(70:180)
  df[numeric_columns] <- lapply(df[numeric_columns], as.numeric)
  results <- df %>%
    group_by(year, site) %>%
    summarise(max_occupancy = max(colSums(across(all_of(numeric_columns)), na.rm = TRUE), na.rm = TRUE)) %>%
    ungroup()
  
  return(results)
}

max_occupancy_df <- max_occupancy(occupancy_data)

max_occupancy_df_elev <- max_occupancy_df %>%
  left_join(sitedat %>% select(site, Mean.Elev), by = "site")

ggplot(max_occupancy_df_elev, aes(x = factor(year), y = max_occupancy, fill = Mean.Elev)) +
  scale_fill_gradient(low = "white", high = "black") +
  geom_violin(fill = "lightgray", alpha = 0.3, color = NA) +
  geom_jitter(shape = 21, size = 2, width = 0.3, height = 0.3) +
  #geom_line(aes(group = site, colour = site)) +
  theme_minimal() +
  labs(title = "Annual Occupancy by Site",
       x = "Year",
       y = "Occupancy",
       fill = "Elevation") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_line(color = "gray", size = 0.5)
  )

y_list <- list()
for (i in site_codes) {
  y_list[[i]] <- max_occupancy_df %>%
    filter(site == i) %>%
    select(max_occupancy) %>%
    as.list()
}

#number of offspring

count_fledglings <- function(df) {
  results <- df %>%
    group_by(year, site) %>%
    summarise(fledgling_count = sum(suc, na.rm = TRUE)) %>%
    ungroup()
  
  return(results)
}

n_offspring_df <- count_fledglings(phendat)

n_offspring_df_elev <- n_offspring_df %>%
  left_join(sitedat %>% select(site, Mean.Elev), by = "site")

ggplot(n_offspring_df_elev, aes(x = factor(year), y = fledgling_count, fill = Mean.Elev)) +
  scale_fill_gradient(low = "white", high = "black") +
  geom_violin(fill = "lightgray", alpha = 0.3, color = NA) +
  geom_jitter(shape = 21, size = 2, width = 0.3, height = 0.3) +
  #geom_line(aes(group = site, colour = site)) +
  theme_minimal() +
  labs(title = "Annual fledgling count by Site",
       x = "Year",
       y = "N fledglings",
       fill = "Elevation") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_line(color = "gray", size = 0.5)
  )

J_list <- list()
for (i in site_codes) {
  J_list[[i]] <- n_offspring_df %>%
  filter(site == i) %>%
  select(fledgling_count) %>%
  as.list()
}

count_broods_surveyed <- function(df) {
  numeric_columns <- as.character(70:180)
  df[numeric_columns] <- lapply(df[numeric_columns], as.numeric)
  results <- df %>%
    group_by(year, site) %>%
    summarise(brood_count = sum(rowSums(across(all_of(numeric_columns))) > 0, na.rm = TRUE)) %>%    ungroup()
  
  return(results)
}

broods_surveyed_df <- count_broods_surveyed(occupancy_data)

broods_surveyed_df_elev <- broods_surveyed_df %>%
  left_join(sitedat %>% select(site, Mean.Elev), by = "site")

ggplot(broods_surveyed_df_elev, aes(x = factor(year), y = brood_count, fill = Mean.Elev)) +
  scale_fill_gradient(low = "white", high = "black") +
  geom_violin(fill = "lightgray", alpha = 0.3, color = NA) +
  geom_jitter(shape = 21, size = 2, width = 0.3, height = 0.3) +
  #geom_line(aes(group = site, colour = site)) +
  theme_minimal() +
  labs(title = "Annual broods surveyed by Site",
       x = "Year",
       y = "N broods",
       fill = "Elevation") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_line(color = "gray", size = 0.5)
  )

R_list <- list()
for (i in site_codes) {
  R_list[[i]] <- broods_surveyed_df %>%
    filter(site == i) %>%
    select(brood_count) %>%
    as.list()
}


# Exercise 4: Once blue tit data are in format similar to here, run the Stan model on blue tit dataset and evaluate.

y_mat_bt <- rbind(y_list)

y_mat_bt <- as.data.frame(y_mat_bt, stringsAsFactors = FALSE)

y_mat_bt <- y_mat_bt %>%
  pivot_longer(cols = everything(), names_to = "site", values_to = "yearly_data") %>%
  rowwise() %>%
  mutate(yearly_data = list(as.numeric(unlist(yearly_data)))) %>%
  mutate(yearly_data = list(c(yearly_data, rep(NA, 12 - length(yearly_data))))) %>%
  unnest_wider(yearly_data, names_sep = "_") %>%
  rename_at(vars(starts_with("yearly_data")), ~ as.character(1:12)) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

J_mat_bt <- rbind(J_list)

J_mat_bt <- as.data.frame(J_mat_bt, stringsAsFactors = FALSE)

J_mat_bt <- J_mat_bt %>%
  pivot_longer(cols = everything(), names_to = "site", values_to = "yearly_data") %>%
  rowwise() %>%
  mutate(yearly_data = list(as.numeric(unlist(yearly_data)))) %>%
  mutate(yearly_data = list(c(yearly_data, rep(NA, 12 - length(yearly_data))))) %>%
  unnest_wider(yearly_data, names_sep = "_") %>%
  rename_at(vars(starts_with("yearly_data")), ~ as.character(1:12)) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

R_mat_bt <- rbind(R_list)

R_mat_bt <- as.data.frame(R_mat_bt, stringsAsFactors = FALSE)

R_mat_bt <- R_mat_bt %>%
  pivot_longer(cols = everything(), names_to = "site", values_to = "yearly_data") %>%
  rowwise() %>%
  mutate(yearly_data = list(as.numeric(unlist(yearly_data)))) %>%
  mutate(yearly_data = list(c(yearly_data, rep(NA, 12 - length(yearly_data))))) %>%
  unnest_wider(yearly_data, names_sep = "_") %>%
  rename_at(vars(starts_with("yearly_data")), ~ as.character(1:12)) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

y_mat_bt <- column_to_rownames(y_mat_bt, var = "site")
J_mat_bt <- column_to_rownames(J_mat_bt, var = "site")
R_mat_bt <- column_to_rownames(R_mat_bt, var = "site")

print(dim(y_mat_bt))
print(dim(J_mat_bt))
print(dim(R_mat_bt))
print(dim(marray_list))

stan_data_bt <- list(
  nyears    = 12,
  nsites    = 44,
  y         = y_mat_bt,      
  J         = J_mat_bt,      
  R         = R_mat_bt,        
  marray_a  = marray_list
)

params_bt <- c(
  "phia", "f", "omega", "p",
  "mphia", "mfec", "mim",
  "lambda", "mlam",
  "sig_phia",      "sig_fec",      "sig_im",
  "sig_site_phia", "sig_site_fec", "sig_site_im",
  "sig_sy_phia",   "sig_sy_fec",   "sig_sy_im",
  "N1", "NadSurv", "Nadimm", "Ntot",
  "var_phia_year", "var_phia_siteyear", "icc_phia",
  "var_fec_year", "var_fec_siteyear", "icc_fec",
  "var_im_year", "var_im_siteyear", "icc_im"
)

inits_bt <- lapply(1:nc, function(i) {
  list(l_mphia = rnorm(1, 0.2, 0.5),
       l_mfec = rnorm(1, 0.2, 0.5),
       l_mim = rnorm(1, 0.2, 0.5),
       l_p = rnorm(1, 0.2, 1),
       sig_phia = runif(1, 0.1, 10),
       sig_fec = runif(1, 0.1, 10),
       sig_im = runif(1, 0.1, 10),
       sig_site_phia = runif(1, 0.1, 10),
       sig_site_fec = runif(1, 0.1, 10),
       sig_site_im = runif(1, 0.1, 10),
       sig_sy_phia = runif(1, 0.1, 10),
       sig_sy_fec = runif(1, 0.1, 10),
       sig_sy_im = runif(1, 0.1, 10),
       N1      = lapply(1:nsites, function(s) round(runif(nyears,  1, 50))),
       NadSurv = lapply(1:nsites, function(s) round(runif(nyears,  5, 50))),
       Nadimm  = lapply(1:nsites, function(s) round(runif(nyears,  1, 50)))
  )})

closeAllConnections()

ipm_bt <- stan(
  file    = "ipm_bt.stan",
  data    = stan_data_bt,
  init    = inits_bt,
  pars    = params_bt,
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

print(ipm_bt,
      pars = c("mphia", "mfec", "mim",
               "sig_phia",      "sig_fec",      "sig_im",
               "sig_site_phia", "sig_site_fec", "sig_site_im",
               "sig_sy_phia",   "sig_sy_fec",   "sig_sy_im",
               "mlam",
               "var_phia_year", "var_phia_siteyear", "icc_phia",
               "var_fec_year", "var_fec_siteyear", "icc_fec",
               "var_im_year", "var_im_siteyear", "icc_im"),
      digits = 3)

# Get posterior distributions
posterior_ipm <- as.data.frame(ipm_bt)
