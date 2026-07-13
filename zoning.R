library(rstan)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(viridis)
library(ClustGeo)
library(sf)

site_codes_2 <- unique(phendat$site)

site_means <- phendat %>%
  group_by(site) %>%                            
  summarise(mean.fed = mean(fed, na.rm = TRUE))

sitedat <- left_join(sitedat, site_means, by = "site")

site_vars <- scale(sitedat[, c("Mean.Long", "Mean.Lat", "Mean.Elev")]) 

D0 <- dist(site_vars)

coords <- st_as_sf(sitedat, coords = c("Mean.Long", "Mean.Lat"), crs = 4326)

D1 <- st_distance(coords) |> as.dist()

range_alpha <- choicealpha(D0, D1, range.alpha = seq(0, 1, 0.1), 
                           K = 12, graph = TRUE)

tree <- hclustgeo(D0, D1, alpha = 0.1)

sitedat$zone <- cutree(tree, k = 16)
