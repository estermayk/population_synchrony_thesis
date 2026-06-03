#---- packages ----

library(lme4)
library(insight)
library(dplyr)
library(ggplot2)
library(tidyverse)

#---- data import and duplicate removal ----

#importing bird phenology datasheet
phendat <- read.csv("data/Bird_Phenology.csv")

#importing site details
sitedat <- read.csv("data/Site_Details.csv")

#removing non focal species
bludat <- phendat[phendat$species == "bluti",]

#changing -999s (0s in fledgling success caused by predation) to 0s
bludat$suc[bludat$suc == "-999"] <- 0

#checking for duplicates
bludat$ID <- as.factor(paste(bludat$year, bludat$site, bludat$box, sep = "_")) 
length(bludat$ID)
length(unique(bludat$ID))
bludat$ID[duplicated(bludat$ID)]
#2230 entries, 5 of which not unique

#---- occupancy plots ----

#creating a loop to create a df containing occupancy by site 
site_codes <- unique(phendat$site)
occupancy <- numeric(length(site_codes))

for (i in seq_along(site_codes)) {
  site <- site_codes[i]
  num_bluti <- nrow(bludat[bludat$site == site, ])
  num_total <- nrow(phendat[phendat$site == site, ])
  occupancy[i] <- num_bluti / num_total
}

occbysite <- data.frame(site = site_codes, occupancy = occupancy)

ggplot(occbysite, aes(x = site, y = occupancy)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Occupancy of Blue Tits by Site",
       x = "Site",
       y = "Occupancy") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#and now by year
years <- unique(phendat$year)
occyear <- numeric(length(years))

for (i in seq_along(years)) {
  year <- years[i]
  num_blutiyear <- nrow(bludat[bludat$year == year, ])
  num_totalyear <- nrow(phendat[phendat$year == year, ])
  occyear[i] <- num_blutiyear / num_totalyear
}

occbyyear <- data.frame(year = years, occupancy = occyear)

ggplot(occbyyear, aes(x = year, y = occupancy)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Occupancy of Blue Tits by Year",
       x = "Year",
       y = "Occupancy") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#year and site
occupancy_by_site_year <- bludat %>%
  group_by(site, year) %>%
  summarise(occupancy = n() / sum(phendat$site == first(site) & phendat$year == first(year))) %>%
  ungroup()

#plotting with box and whiskers
ggplot(occupancy_by_site_year, aes(x = factor(year), y = occupancy)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 21) +
  geom_text(data = occupancy_by_site_year %>%
              filter(occupancy > (quantile(occupancy, 0.75) + 1.5 * IQR(occupancy)) |
                       occupancy < (quantile(occupancy, 0.25) - 1.5 * IQR(occupancy))),
            aes(label = site), vjust = -0.5, hjust = 0.5, size = 3) +
  theme_minimal() +
  expand_limits(y = 0,) +
  labs(title = "Annual Occupancy",
       x = "Year",
       y = "Occupancy")

?scale_fill_gradient

occupancy_by_site_year <- occupancy_by_site_year %>%
  left_join(sitedat %>% select(site, Mean.Elev), by = "site")

ggplot(occupancy_by_site_year, aes(x = factor(year), y = occupancy, fill = Mean.Elev)) +
  scale_fill_gradient(low = "white", high = "black") +
  geom_violin(fill = "lightgray", alpha = 0.3, color = NA) +
  geom_jitter(shape = 21, size = 2, width = 0.4, height = 0) +
  theme_minimal() +
  labs(title = "Annual Occupancy by Site",
       x = "Year",
       y = "Occupancy",
       fill = "Elevation") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_line(color = "gray", size = 0.5)
  )
  
?geom_jitter

#---- initial suc mod ----

#very early model attempt on fledgling success
bludat$site_year <- as.factor(paste(bludat$site, bludat$year, sep = "_")) 
sucmod <- lmer(suc ~ (1|year) + (1|site) + (1|site_year), data = bludat) 
summary(sucmod)

get_variance_intercept(sucmod)

#extracting the variance components 
var_components <- VarCorr(sucmod)
var_intercept_site_year <- attr(var_components$site_year, "stddev")^2
var_intercept_site <- attr(var_components$site, "stddev")^2
var_intercept_year <- attr(var_components$year, "stddev")^2

#a crude look at ICC synchrony from the model 
sucsync <- var_intercept_year/(var_intercept_year + var_intercept_site_year)

sucsync

#---- clutch swap variance ----

#lets look at the variance caused by clutch swap
bludatswap <- bludat %>%
  mutate(clutch.swap.treatment = ifelse(clutch.swap.treatment %in% 1:9, "swap", clutch.swap.treatment))

bludatswap <- bludatswap[bludatswap$year %in% c("2017", "2018", "2019"), ]

bludatswap$clutch.swap.treatment <- substr(bludatswap$clutch.swap.treatment, 1, 3)

clutchswapmod <- lm(suc ~ clutch.swap.treatment, data = bludatswap)
summary(clutchswapmod)

swapmodmixed <- lmer(suc ~ clutch.swap.treatment + (1|site) + (1|year), data = bludatswap)
summary(swapmodmixed)

#relevel to make unm the reference
str(bludatswap)
bludatswap$clutch.swap.treatment <- as.factor(bludatswap$clutch.swap.treatment)
bludatswap$clutch.swap.treatment <- relevel(bludatswap$clutch.swap.treatment, ref = "unm")

#making year a factor to reduce noise
bludatswap$year <- as.factor(bludatswap$year)

#removing nas in swap
bludatswap <- bludatswap %>%
  mutate(clutch.swap.treatment = na_if(clutch.swap.treatment, "")) %>%
  drop_na(clutch.swap.treatment)

#popping odd and not possible together 
bludatswap <- bludat %>%
  mutate(clutch.swap.treatment = as.character(clutch.swap.treatment)) %>%
  mutate(clutch.swap.treatment = ifelse(clutch.swap.treatment == "not", "odd", clutch.swap.treatment))

swapmodmixed <- lmer(suc ~ (1|clutch.swap.treatment) + year + (1|site), data = bludatswap)
summary(swapmodmixed)

#think about nesting - year should not be swapped 

tapply(bludatswap$clutch.swap.treatment, bludatswap$clutch.swap.treatment, length)
tapply(bludatswap$suc, bludatswap$clutch.swap.treatment, mean, na.rm=TRUE)

