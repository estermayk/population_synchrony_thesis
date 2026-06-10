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

ALN_adults_marray <- marray(ALN_adults_CH_matrix)

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
ord_dates_phendat <- select(phendat, c('year', 'site', 'ID', 'n1', 'nl', 'latestfed', 'latestcc', 'fed', 'cc', 'cs', 'fki', 'hatching_first_recorded', 'v1date', 'v2date'))

transform_dates <- function(df) {
  date_columns <- 100:180
  occupancy_df <- matrix(0, ncol = length(date_columns), nrow = nrow(df))
  colnames(occupancy_df) <- as.character(date_columns)
  
  for (i in 1:nrow(df)) {
    dates <- suppressWarnings(as.numeric(df[i, c('n1', 'nl', 'latestfed', 'latestcc', 'fed', 'cc', 
                                                 'cs', 'fki', 'hatching_first_recorded', 'v1date', 
                                                 'v2date')]))    
    date_range <- dates[!is.na(dates)]
    if (length(date_range) > 0) {
      min_date <- max(min(date_range), 100)
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
  numeric_columns <- as.character(100:180)
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

J_list <- list()
for (i in site_codes) {
  J_list[[i]] <- n_offspring_df %>%
  filter(site == i) %>%
  select(fledgling_count) %>%
  as.list()
}

count_broods_surveyed <- function(df) {
  numeric_columns <- as.character(100:180)
  df[numeric_columns] <- lapply(df[numeric_columns], as.numeric)
  results <- df %>%
    group_by(year, site) %>%
    summarise(brood_count = sum(rowSums(across(all_of(numeric_columns))) > 0, na.rm = TRUE)) %>%    ungroup()
  
  return(results)
}

broods_surveyed_df <- count_broods_surveyed(occupancy_data)

R_list <- list()
for (i in site_codes) {
  R_list[[i]] <- broods_surveyed_df %>%
    filter(site == i) %>%
    select(brood_count) %>%
    as.list()
}


# Exercise 4: Once blue tit data are in format similar to here, run the Stan model on blue tit dataset and evaluate.


