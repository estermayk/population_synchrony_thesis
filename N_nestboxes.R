boxes <- read.csv("data/Boxes.csv")

names(boxes) <- c("Zone", "Site", "2014", "2015", "2016", "2017")

boxes_long <- boxes %>% 
  pivot_longer(cols = c("2014", "2015", "2016", "2017"), 
               names_to = "year", values_to = "nestboxes") %>% 
  mutate(year = recode(year, "2017-" = "2017+")) %>% 
  group_by(year, Zone) %>% 
  summarise(total_nestboxes = sum(nestboxes, na.rm = TRUE), 
            .groups = "drop")

ggplot(boxes_long, aes(x = year, y = total_nestboxes, fill = Zone)) + 
  geom_bar(stat = "identity", colour = "white", width = 0.6) + 
  geom_text(aes(label = after_stat(y), group = year),
                                                                                                                                              stat = "summary", fun = sum, vjust = -0.5, size = 3.5, colour = "black") + scale_fill_viridis_d(option = "turbo") + labs(title = "Total nestboxes across transect by year", x = "Year", y = "N nestboxes", fill = "Zone") + theme_bw(base_size = 12) + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2))
ggplot(boxes_long, aes(x = Zone, y = total_nestboxes, fill = year)) +
  geom_bar(stat     = "identity",
           position = "dodge",
           colour   = "white",
           width    = 0.7) +
  scale_fill_viridis_d(option = "turbo") +
  labs(title = "Nestboxes per zone by year",
       x     = "Zone",
       y     = "N nestboxes",
       fill  = "Year") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

nnestboxes <- ggplot(boxes_long, aes(x = year, y = total_nestboxes, fill = Zone)) +
  geom_bar(stat     = "identity",
           position = "dodge",
           colour   = "white",
           width    = 0.7) +
  scale_fill_viridis_d(option = "turbo") +
  labs(title = "Nestboxes per zone by year",
       x     = "Year",
       y     = "N nestboxes",
       fill  = "Zone") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))

ggsave("figs/N_Nestboxes.png", nnestboxes, width = 9, height = 6)
