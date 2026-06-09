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



# Exercise 4: Once blue tit data are in format similar to here, run the Stan model on blue tit dataset and evaluate.

