library(rstan)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(viridis)
library(ClustGeo)
library(sf)
library(zoom)
library(ggdendro)

site_codes_2 <- unique(phendat$site)

site_means <- phendat %>%
  group_by(site) %>%                            
  summarise(mean.fed = mean(fed, na.rm = TRUE))

sitedat <- left_join(sitedat, site_means, by = "site")

sitedat$management..rough. <- as.factor(sitedat$management..rough.)

site_vars <- scale(sitedat[, c("Mean.Long", "Mean.Lat", "Mean.Elev")]) 

D0 <- dist(site_vars)

coords <- st_as_sf(sitedat, coords = c("Mean.Long", "Mean.Lat"), crs = 4326)

D1 <- st_distance(coords) |> as.dist()

range_alpha <- choicealpha(D0, D1, range.alpha = seq(0, 1, 0.1), 
                           K = 12, graph = TRUE)

testtree <- hclustgeo(D0, D1, alpha = 0.8)


colours_12 <- viridis::viridis(12)
cut_h <- sort(testtree$height, decreasing = TRUE)[12]
cap_h <- cut_h * 4

testtree_viz         <- testtree
testtree_viz$height  <- pmin(testtree$height, cap_h)
testtree_viz$labels  <- sitedat$site

clusters  <- cutree(testtree, k = 12)   
leaf_ord  <- order.dendrogram(as.dendrogram(testtree_viz))
clust_lr  <- clusters[testtree_viz$labels[leaf_ord]]  

xmids  <- tapply(seq_along(leaf_ord), clust_lr, mean)  
xlefts <- tapply(seq_along(leaf_ord), clust_lr, min)   

lr_rank     <- order(xlefts)
plot_x      <- as.numeric(xmids)[lr_rank]           
plot_labels <- LETTERS[1:12]                          
plot_cols   <- colours_12                            

plot(testtree_viz,
     hang = -1, cex = 0.65,
     xlab = "", sub  = "",
     main = "ClustGeo Dendrogram (k = 12)",
     ylab = "Height")

rect.hclust(testtree_viz, k = 12, border = colours_17)  

plot(testtree)

sitedat$zone <- cutree(testtree, k = 12)

sitedat$zone <- LETTERS[as.integer(sitedat$zone)]

head(sitedat)

adults <- adults[, -which(names(adults) == "zone.y")]
phendat <- phendat[, -which(names(phendat) == "zone")]
head(phendat)

adults  <- adults  %>% left_join(sitedat[, c("site", "zone")], by = "site")
phendat <- phendat %>% left_join(sitedat[, c("site", "zone")], by = "site")
head(phendat)

#no matter how i split it, i get zeros in my m-array from J and K - which are closest to each other so even dropping below 10 sites would still cause this
#giving up on this and using the simple pairwise version which ensures sufficient sampling at each zone 