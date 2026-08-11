#---- posterior site year means ----

# Get nyears and nsites
nyears_bt  <- stan_data_bt$nyears
nsites_bt  <- stan_data_bt$nsites
years_bt   <- 1:nyears_bt
years_m1_bt <- 1:(nyears_bt - 1)  # for parameters indexed over nyears-1

site_labels <- paste0(unique(adults$zone))

site_labels

survey_years <- c(2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025)
survey_years_m1 <- c(2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024)

extract_summary <- function(posterior, pattern) {
  cols <- grep(pattern, names(posterior), value = TRUE)
  data.frame(
    param = cols,
    mean  = colMeans(posterior[, cols, drop = FALSE]),
    lower = apply(posterior[, cols, drop = FALSE], 2, quantile, 0.025),
    upper = apply(posterior[, cols, drop = FALSE], 2, quantile, 0.975)
  )
}

# Estimated Ntot per site per year
ntot_list <- lapply(1:nsites_bt, function(s) {
  pattern <- paste0("Ntot\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- survey_years
  df$site <- site_labels[s]
  df
})
ntot_df <- bind_rows(ntot_list)

obs_counts$site[1:12]
survey_years
nyears_bt

# Observed counts per site per year
obs_counts <- expand.grid(site = site_labels, year = survey_years) %>%
  arrange(site, year) %>%
  mutate(observed = as.vector(t(stan_data_bt$y)))  # y is [nsites, nyears]

obs_counts <- data.frame(
  site     = rep(site_labels, each = nyears_bt),   # K×12, J×12, H×12...
  year     = rep(survey_years, times = nsites_bt),  # 2014:2025 for each site
  observed = as.vector(t(stan_data_bt$y))           # matches site_codes order
)

# Merge both dfs
ntot_df <- left_join(ntot_df, obs_counts, by = c("site", "year"))

mean(ntot_df$mean)

# Plot
p_ntot <- ggplot(ntot_df, aes(x = year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = mean, colour = "Estimated Population Size"), linewidth = 0.9) +
  geom_point(aes(y = mean, colour = "Estimated Population Size"), size = 1.5) +
  geom_line(aes(y = observed, colour = "Observed count"),
            linewidth = 0.9, linetype = "dashed") +
  geom_point(aes(y = observed, colour = "Observed count"), size = 1.5) +
  scale_colour_manual(values = c("Estimated Population Size" = "steelblue",
                                 "Observed count"  = "firebrick")) +
  facet_wrap(~ site, scales = "free_y") +
  labs(title = "Population size", x = "Year", y = "N",
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(n.breaks=8)

p_ntot

ggsave("figs/pntot.png", p_ntot, width = 9, height = 4)

p_ntot2 <- ggplot(ntot_df, aes(x = year, group = site, colour = site, fill = site)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_line(aes(y = mean,     linetype = "Estimated Population Size"), linewidth = 0.9) +
  geom_point(aes(y = mean,    shape    = "Estimated Population Size"), size = 1.5) +
 # geom_line(aes(y = observed, linetype = "Observed count"), linewidth = 0.9) +
#  geom_point(aes(y = observed, shape   = "Observed count"), size = 1.5) +
#  scale_linetype_manual(values = c("Estimated Population Size" = "solid",
 #                                  "Observed count"  = "dashed")) +
#  scale_shape_manual(values   = c("Estimated Population Size" = 16,
 #                                 "Observed count"  = 1)) +
  scale_colour_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option   = "turbo") +
  labs(title    = "Population size", 
       x        = "Year", 
       y        = "N",
       linetype = NULL, 
       shape    = NULL,
       colour   = "Zone", 
       fill     = "Zone") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(n.breaks = 8) +
  guides(colour   = guide_legend(nrow = 2),
         fill     = guide_legend(nrow = 2),
        # linetype = guide_legend(nrow = 2),
         shape    = guide_legend(nrow = 2))

ggsave("figs/pntot2.png", p_ntot2, width = 8, height = 8)

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
f_list <- lapply(1:nsites_bt, function(s) {
  pattern <- paste0("f\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- survey_years_m1
  df$site <- site_labels[s]
  df
})
f_est <- bind_rows(f_list)

# Observed f per site per year
obs_f <- expand.grid(site = site_labels, year = survey_years_m1) %>%
  arrange(site, year) %>%
  mutate(observed = as.vector(t(stan_data_bt$J / stan_data_bt$R)))

# Merge both dfs
f_df <- left_join(f_est, obs_f, by = c("site", "year"))

mean(f_df$mean)

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
  labs(title = "Productivity per Nest", x = "Year", y = "N Fledglings/female",
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(n.breaks=8)

p_f

p_f2 <- ggplot(f_df, aes(x = year, group = site, colour = site, fill = site)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_line(aes(y = mean,     linetype = "Estimated Productivity"), linewidth = 0.9) +
  geom_point(aes(y = mean,    shape    = "Estimated Productivity"), size = 1.5) +
  # geom_line(aes(y = observed, linetype = "Observed count"), linewidth = 0.9) +
  #  geom_point(aes(y = observed, shape   = "Observed count"), size = 1.5) +
  #  scale_linetype_manual(values = c("Estimated Population Size" = "solid",
  #                                  "Observed count"  = "dashed")) +
  #  scale_shape_manual(values   = c("Estimated Population Size" = 16,
  #                                 "Observed count"  = 1)) +
  scale_colour_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option   = "turbo") +
  labs(title    = "Productivity", 
       x        = "Year", 
       y        = "N",
       linetype = NULL, 
       shape    = NULL,
       colour   = "Zone", 
       fill     = "Zone") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(n.breaks = 8) +
  guides(colour   = guide_legend(nrow = 2),
         fill     = guide_legend(nrow = 2),
         # linetype = guide_legend(nrow = 2),
         shape    = guide_legend(nrow = 2))
#plot on one for main text
#include original to show estimated vs observed tracking 

ggsave("figs/pprod2.png", p_f2, width = 8, height = 8)


#survival
#juveniles
# Estimated phij per site per year
phij_list <- lapply(1:nsites_bt, function(s) {
  pattern <- paste0("phij\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- survey_years_m1
  df$site <- site_labels[s]
  df
})
phij_est <- bind_rows(phij_list)

#adults
phia_list <- lapply(1:nsites_bt, function(s) {
  pattern <- paste0("phia\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- survey_years_m1
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

mean(phi_df$mean.x)
mean(phi_df$mean.y)

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

p_phi2 <- ggplot(phi_df, aes(x = year, group = site, colour = site, fill = site)) +
  geom_ribbon(aes(ymin = lower.x, ymax = upper.x), alpha = 0.12, colour = NA) +
  geom_ribbon(aes(ymin = lower.y, ymax = upper.y), alpha = 0.12, colour = NA) +
  geom_line(aes(y = mean.y,     linetype = "Estimated Adult Survival"), linewidth = 0.9) +
  geom_point(aes(y = mean.y,    shape    = "Estimated Adult Survival"), size = 1.5) +
   geom_line(aes(y = mean.x, linetype = "Estimated Juvenile Survival"), linewidth = 0.9) +
    geom_point(aes(y = mean.x, shape   = "Estimated Juvenile Survival"), size = 1.5) +
    scale_linetype_manual(values = c("Estimated Adult Survival" = "solid",
                                    "Estimated Juvenile Survival"  = "dashed")) +
    scale_shape_manual(values   = c("Estimated Adult Survival" = 16,
                                   "Estimated Juvenile Survival"  = 1)) +
  scale_colour_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option   = "turbo") +
  labs(title    = "Survival", 
       x        = "Year", 
       y        = "Phi",
       linetype = NULL, 
       shape    = NULL,
       colour   = "Zone", 
       fill     = "Zone") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(n.breaks = 8) +
  guides(colour   = guide_legend(nrow = 2),
         fill     = guide_legend(nrow = 2),
         # linetype = guide_legend(nrow = 2),
         shape    = guide_legend(nrow = 2))

ggsave("figs/pphi2.png", p_phi2, width = 8, height = 8)


#population growth rate

lambda_list <- lapply(1:nsites_bt, function(s) {
  pattern <- paste0("lambda\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- survey_years_m1
  df$site <- site_labels[s]
  df
})
lambda_df <- bind_rows(lambda_list)

mean(lambda_df$mean)

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

lambda_df <- lambda_df %>%
  filter(!year %in% c("2014", "2015", "2016"))

p_lambda2 <- ggplot(lambda_df, aes(x = year, group = site, colour = site, fill = site)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_line(aes(y = mean,     linetype = "Estimated lambda"), linewidth = 0.9) +
  geom_point(aes(y = mean,    shape    = "Estimated lambda"), size = 1.5) +
  # geom_line(aes(y = observed, linetype = "Observed count"), linewidth = 0.9) +
  #  geom_point(aes(y = observed, shape   = "Observed count"), size = 1.5) +
  #  scale_linetype_manual(values = c("Estimated Population Size" = "solid",
  #                                  "Observed count"  = "dashed")) +
  #  scale_shape_manual(values   = c("Estimated Population Size" = 16,
  #                                 "Observed count"  = 1)) +
  scale_colour_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option   = "turbo") +
  labs(title    = "Population growth rate", 
       x        = "Year", 
       y        = "Lambda",
       linetype = NULL, 
       shape    = NULL,
       colour   = "Zone", 
       fill     = "Zone") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(n.breaks = 8) +
  guides(colour   = guide_legend(nrow = 2),
         fill     = guide_legend(nrow = 2),
         # linetype = guide_legend(nrow = 2),
         shape    = guide_legend(nrow = 2)) +
  coord_cartesian(ylim = c(0.5, 2))

p_lambda2

ggsave("figs/plambda2_ylim.png", p_lambda2, width = 8, height = 8)


#immigration
Nadimm_list <- lapply(1:nsites_bt, function(s) {
  pattern <- paste0("Nadimm\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- survey_years
  df$site <- site_labels[s]
  df
})
Nadimm_df <- bind_rows(Nadimm_list)

mean(Nadimm_df$mean)

# Observed counts per site per year
obs_counts <- expand.grid(site = site_labels, year = survey_years) %>%
  arrange(site, year) %>%
  mutate(observed = as.vector(t(stan_data_bt$y)))  # y is [nsites, nyears]

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

p_Nadimm2 <- ggplot(Nadimm_df, aes(x = year, group = site, colour = site, fill = site)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_line(aes(y = mean,     linetype = "Estimated Nadimm"), linewidth = 0.9) +
  geom_point(aes(y = mean,    shape    = "Estimated Nadimm"), size = 1.5) +
  # geom_line(aes(y = observed, linetype = "Observed count"), linewidth = 0.9) +
  #  geom_point(aes(y = observed, shape   = "Observed count"), size = 1.5) +
  #  scale_linetype_manual(values = c("Estimated Population Size" = "solid",
  #                                  "Observed count"  = "dashed")) +
  #  scale_shape_manual(values   = c("Estimated Population Size" = 16,
  #                                 "Observed count"  = 1)) +
  scale_colour_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option   = "turbo") +
  labs(title    = "Immigration", 
       x        = "Year", 
       y        = "N",
       linetype = NULL, 
       shape    = NULL,
       colour   = "Zone", 
       fill     = "Zone") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(n.breaks = 8) +
  guides(colour   = guide_legend(nrow = 2),
         fill     = guide_legend(nrow = 2),
         # linetype = guide_legend(nrow = 2),
         shape    = guide_legend(nrow = 2))

p_Nadimm2

ggsave("figs/p_Nadimm2.png", p_Nadimm2, width = 8, height = 8)


rates_plots <- (p_ntot2 | p_lambda2 | p_f2) / (p_phi2 | p_Nadimm2)

rates_plots

ggsave("figs/rates_plot.png", plot = rates_plots)


p_list <- lapply(1:nsites_bt, function(s) {
  pattern <- paste0("p\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- survey_years_m1
  df$site <- site_labels[s]
  df
})
p_df <- bind_rows(p_list)

mean(p_df$mean)

# Observed counts per site per year
#obs_counts <- expand.grid(site = site_labels, year = years) %>%
# arrange(site, year) %>%
#mutate(observed = as.vector(t(stan_data$y)))  # y is [nsites, nyears]

# Merge both dfs
#p_df <- left_join(p_df, obs_counts, by = c("site", "year"))

# Plot
p_p <- ggplot(p_df, aes(x = year, group = site, colour = site, fill = site)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_line(aes(y = mean,     linetype = "Estimated p"), linewidth = 0.9) +
  geom_point(aes(y = mean,    shape    = "Estimated p"), size = 1.5) +
  # geom_line(aes(y = observed, linetype = "Observed count"), linewidth = 0.9) +
  #  geom_point(aes(y = observed, shape   = "Observed count"), size = 1.5) +
  #  scale_linetype_manual(values = c("Estimated Population Size" = "solid",
  #                                  "Observed count"  = "dashed")) +
  #  scale_shape_manual(values   = c("Estimated Population Size" = 16,
  #                                 "Observed count"  = 1)) +
  scale_colour_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option   = "turbo") +
  labs(title    = "Detection Probability", 
       x        = "Year", 
       y        = "p",
       linetype = NULL, 
       shape    = NULL,
       colour   = "Zone", 
       fill     = "Zone") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(n.breaks = 8) +
  guides(colour   = guide_legend(nrow = 2),
         fill     = guide_legend(nrow = 2),
         # linetype = guide_legend(nrow = 2),
         shape    = guide_legend(nrow = 2))


p_p

ggsave("figs/p_p.png", p_p, width = 8, height = 4)


#---- ICC and var distributions ----

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
       x = "phia ICC",,
       tag = "c)",
       y = "Density") +
  theme_minimal()

icc_phia_p

#var_phia_year

var_phia_year <- posterior_ipm %>%
  select(contains("var_phia_year")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_phia_year_p <- ggplot(var_phia_year, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Adult Survival Year Variances",,
       tag = "a)",
       x = "phia year var",
       y = "Density") +
  theme_minimal()

var_phia_year_p

var_phia_siteyear <- posterior_ipm %>%
  select(contains("var_phia_siteyear")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_phia_siteyear_p <- ggplot(var_phia_siteyear, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Adult Survival Zone Year Variances",,
       tag = "b)",
       x = "phia zone year var",
       y = "Density") +
  theme_minimal()

var_phia_siteyear_p

#im icc

im_icc <- posterior_ipm %>%
  select(contains("icc_im")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

icc_im_p <- ggplot(im_icc, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Immigration ICC",,
       tag = "i)",
       x = "im ICC",
       y = "Density") +
  theme_minimal()

icc_im_p

#var_im_year

var_im_year <- posterior_ipm %>%
  select(contains("var_im_year")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_im_year_p <- ggplot(var_im_year, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Immigration Year Variances",,
       tag = "g)",
       x = "im year var",
       y = "Density") +
  theme_minimal()

var_im_year_p

var_im_siteyear <- posterior_ipm %>%
  select(contains("var_im_siteyear")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_im_siteyear_p <- ggplot(var_im_siteyear, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Immigration Zone Year Variances",,
       tag = "h)",
       x = "im zone year var",
       y = "Density") +
  theme_minimal()

var_im_siteyear_p

#prod
prod_icc <- posterior_ipm %>%
  select(contains("icc_prod")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

icc_prod_p <- ggplot(prod_icc, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Productivity ICC",,
       tag = "f)",
       x = "prod ICC",
       y = "Density") +
  theme_minimal()

icc_prod_p

#var_prod_year

var_prod_year <- posterior_ipm %>%
  select(contains("var_prod_year")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_prod_year_p <- ggplot(var_prod_year, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Productivity Year Variances",,
       tag = "d)",
       x = "prod year var",
       y = "Density") +
  theme_minimal()

var_prod_year_p

var_prod_siteyear <- posterior_ipm %>%
  select(contains("var_prod_siteyear")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_prod_siteyear_p <- ggplot(var_prod_siteyear, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Productivity Zone Year Variances",
       tag = "e)",
       x = "prod zone year var",
       y = "Density") +
  theme_minimal()

var_prod_siteyear_p

icc_var_plots <- (var_phia_year_p | var_phia_siteyear_p | icc_phia_p) / (var_prod_year_p | var_prod_siteyear_p | icc_prod_p) / (var_im_year_p | var_im_siteyear_p |icc_im_p)

icc_var_plots

icc_var_phia <- (var_phia_year_p | icc_phia_p)

icc_var_prod <- (var_prod_year_p | icc_prod_p)

icc_var_imm <- (var_im_year_p | icc_im_p)

ggsave("figs/icc_var_plots_bt.png", plot = icc_var_plots, width = 12, height = 6)
ggsave("figs/icc_var_phia.png", plot = icc_var_phia, width = 8, height = 4)
ggsave("figs/icc_var_prod.png", plot = icc_var_prod, width = 8, height = 4)
ggsave("figs/icc_var_imm.png", plot = icc_var_imm, width = 8, height = 4)


#p icc

p_icc <- posterior_ipm %>%
  select(contains("icc_p")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

icc_p_p <- ggplot(p_icc, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Detection Probability ICC",
       x = "p ICC",,
       tag = "c)",
       y = "Density") +
  theme_minimal()

icc_p_p

#var_p_year

var_p_year <- posterior_ipm %>%
  select(contains("var_p_year")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_p_year_p <- ggplot(var_p_year, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Detection Probability Year Variances",,
       tag = "a)",
       x = "p year var",
       y = "Density") +
  theme_minimal()

var_p_year_p

var_p_siteyear <- posterior_ipm %>%
  select(contains("var_p_siteyear")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_p_siteyear_p <- ggplot(var_p_siteyear, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Detection Probability Zone Year Variances",,
       tag = "b)",
       x = "p zone year var",
       y = "Density") +
  theme_minimal()

var_p_siteyear_p

#---- lambda correlations ----

lambda_df$site_year <- as.factor(paste(lambda_df$year, lambda_df$site, sep = "_")) 

lambda_siteyear <- lambda_df %>%
  group_by(site_year) %>%
  summarise(mean_lambda = mean,
            lower_lambda = lower,
            upper_lambda = upper,
            year = year)

phia_est$site_year <- as.factor(paste(phia_est$year, phia_est$site, sep = "_")) 

phia_siteyear <- phia_est %>%
  group_by(site_year) %>%
  summarise(mean_phia = mean,
            lower_phia = lower,
            upper_phia = upper,
            year = year)

phia_lambda <- left_join(lambda_siteyear, phia_siteyear, by = "site_year")

phia_lambda <- phia_lambda %>%
  filter(!year.x %in% c("2014", "2015", "2016"))

phia_lambda_p <- ggplot(phia_lambda, aes(x = mean_phia, y = mean_lambda)) +
  #geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
                #colour = "grey70", width = 0) +
  #geom_errorbarh(aes(xmin = lower_phia, xmax = upper_phia), 
                 #colour = "grey70", height = 0) +
  geom_point(size = 1.5, colour = "steelblue") +
  #geom_text(nudge_y = 0.02, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
  labs(tag = "a)",
       x = "Mean adult survival (φa)",
       y = "Mean λ",
       title = "Population growth rate vs adult survival by zone year") +
  theme_bw(base_size = 12) +
  coord_cartesian(ylim = c(0.7, 1.3)) 


phij_est$site_year <- as.factor(paste(phij_est$year, phij_est$site, sep = "_")) 

phij_siteyear <- phij_est %>%
  group_by(site_year) %>%
  summarise(mean_phij = mean,
            lower_phij = lower,
            upper_phij = upper,
            year = year)

phij_lambda <- left_join(lambda_siteyear, phij_siteyear, by = "site_year")

phij_lambda <- phij_lambda %>%
  filter(!year.x %in% c("2014", "2015", "2016"))

phij_lambda_p <- ggplot(phij_lambda, aes(x = mean_phij, y = mean_lambda)) +
  #geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
  #colour = "grey70", width = 0) +
  #geom_errorbarh(aes(xmin = lower_phij, xmax = upper_phij), 
  #colour = "grey70", height = 0) +
  geom_point(size = 1.5, colour = "steelblue") +
  #geom_text(nudge_y = 0.02, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
  labs(tag = "b)",
       x = "Mean juvenile survival (φj)",
       y = "Mean λ",
       title = "Population growth rate vs juvenile survival by zone year") +
  theme_bw(base_size = 12)+
  coord_cartesian(ylim = c(0.7, 1.3))



f_est$site_year <- as.factor(paste(f_est$year, f_est$site, sep = "_")) 

prod_siteyear <- f_est %>%
  group_by(site_year) %>%
  summarise(mean_prod = mean,
            lower_prod = lower,
            upper_prod = upper,
            year = year)

prod_lambda <- left_join(lambda_siteyear, prod_siteyear, by = "site_year")

prod_lambda <- prod_lambda %>%
  filter(!year.x %in% c("2014", "2015", "2016"))

prod_lambda_p <- ggplot(prod_lambda, aes(x = mean_prod, y = mean_lambda)) +
  #geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
  #colour = "grey70", width = 0) +
  #geom_errorbarh(aes(xmin = lower_prod, xmax = upper_prod), 
  #colour = "grey70", height = 0) +
  geom_point(size = 1.5, colour = "steelblue") +
  #geom_text(nudge_y = 0.02, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
  labs(tag = "c)",
       x = "Mean productivity",
       y = "Mean λ",
       title = "Population growth rate vs productivity by zone year") +
  theme_bw(base_size = 12)+
  coord_cartesian(ylim = c(0.7, 1.3))


#ntot_df$site_year <- as.factor(paste(ntot_df$year, ntot_df$site, sep = "_")) 

#ntot_siteyear <- ntot_df %>%
#  group_by(site_year) %>%
#  summarise(mean_ntot = mean,
#            lower_ntot = lower,
#            upper_ntot = upper,
#            year = year)

#ntot_lambda <- left_join(lambda_siteyear, ntot_siteyear, by = "site_year")

#ntot_lambda <- ntot_lambda %>%
#  filter(!year.x %in% c("2014", "2015", "2016"))

#ntot_lambda_p <- ggplot(ntot_lambda, aes(x = mean_ntot, y = mean_lambda)) +
#  #geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
#  #colour = "grey70", width = 0) +
#  #geom_errorbarh(aes(xmin = lower_ntot, xmax = upper_ntot), 
#  #colour = "grey70", height = 0) +
#  geom_point(size = 1.5, colour = "steelblue") +
#  #geom_text(nudge_y = 0.02, size = 3) +
#  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
#  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
#  labs(tag = "b)",
#       x = "Mean Population Size",
#       y = "Mean λ",
#       title = "Population growth rate vs population size by zone year") +
#  theme_bw(base_size = 12)

Nadimm_df$site_year <- as.factor(paste(Nadimm_df$year, Nadimm_df$site, sep = "_")) 

im_siteyear <- Nadimm_df %>%
  group_by(site_year) %>%
  summarise(mean_im = mean,
            lower_im = lower,
            upper_im = upper,
            year = year)

im_lambda <- left_join(lambda_siteyear, im_siteyear, by = "site_year")

im_lambda <- im_lambda %>%
  filter(!year.x %in% c("2014", "2015", "2016"))

im_lambda_p <- ggplot(im_lambda, aes(x = mean_im, y = mean_lambda)) +
  #geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
  #colour = "grey70", width = 0) +
  #geom_errorbarh(aes(xmin = lower_im, xmax = upper_im), 
  #colour = "grey70", height = 0) +
  geom_point(size = 1.5, colour = "steelblue") +
  #geom_text(nudge_y = 0.02, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
  labs(tag = "d)",
       x = "Mean immigrants",
       y = "Mean λ",
       title = "Population growth rate vs immigration by zone year") +
  theme_bw(base_size = 12)+
  coord_cartesian(ylim = c(0.7, 1.3))

cor_ps <- (phia_lambda_p / prod_lambda_p | phij_lambda_p / im_lambda_p)

(phia_lambda_p / phij_lambda_p)

cor.test(phia_lambda$mean_lambda, phia_lambda$mean_phia)

cor.test(phij_lambda$mean_lambda, phij_lambda$mean_phij)

cor.test(prod_lambda$mean_lambda, prod_lambda$mean_prod)

#cor.test(ntot_lambda$mean_lambda, ntot_lambda$mean_ntot)

cor.test(im_lambda$mean_lambda, im_lambda$mean_im)

ggsave("figs/cor_ps.png", cor_ps, width = 12, height = 10)

#---- lambda cor with error bars ----

phia_lambda_p_CrI <- ggplot(phia_lambda, aes(x = mean_phia, y = mean_lambda)) +
  geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
  colour = "grey70", width = 0) +
  geom_errorbarh(aes(xmin = lower_phia, xmax = upper_phia), 
  colour = "grey70", height = 0) +
  geom_point(size = 1.5, colour = "steelblue") +
  #geom_text(nudge_y = 0.02, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
  labs(tag = "a)",
       x = "Mean adult survival (φa)",
       y = "Mean λ",
       title = "Population growth rate vs adult survival by zone year") +
  theme_bw(base_size = 12)

phij_lambda_p_CrI <- ggplot(phij_lambda, aes(x = mean_phij, y = mean_lambda)) +
  geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
                colour = "grey70", width = 0) +
  geom_errorbarh(aes(xmin = lower_phij, xmax = upper_phij), 
                 colour = "grey70", height = 0) +
  geom_point(size = 1.5, colour = "steelblue") +
  #geom_text(nudge_y = 0.02, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
  labs(tag = "b)",
       x = "Mean juvenile survival (φj)",
       y = "Mean λ",
       title = "Population growth rate vs juvenile survival by zone year") +
  theme_bw(base_size = 12)

prod_lambda_p_CrI <- ggplot(prod_lambda, aes(x = mean_prod, y = mean_lambda)) +
  geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
  colour = "grey70", width = 0) +
  geom_errorbarh(aes(xmin = lower_prod, xmax = upper_prod), 
  colour = "grey70", height = 0) +
  geom_point(size = 1.5, colour = "steelblue") +
  #geom_text(nudge_y = 0.02, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
  labs(tag = "c)",
       x = "Mean productivity",
       y = "Mean λ",
       title = "Population growth rate vs productivity by zone year") +
  theme_bw(base_size = 12)

#ntot_lambda_p_CrI <- ggplot(ntot_lambda, aes(x = mean_ntot, y = mean_lambda)) +
#  geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
#  colour = "grey70", width = 0) +
#  geom_errorbarh(aes(xmin = lower_ntot, xmax = upper_ntot), 
#  colour = "grey70", height = 0) +
#  geom_point(size = 1.5, colour = "steelblue") +
#  #geom_text(nudge_y = 0.02, size = 3) +
#  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
#  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
#  labs(tag = "b)",
#       x = "Mean Population Size",
#       y = "Mean λ",
#       title = "Population growth rate vs population size by zone year") +
#  theme_bw(base_size = 12)

im_lambda_p_CrI <- ggplot(im_lambda, aes(x = mean_im, y = mean_lambda)) +
  geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
  colour = "grey70", width = 0) +
  geom_errorbarh(aes(xmin = lower_im, xmax = upper_im), 
  colour = "grey70", height = 0) +
  geom_point(size = 1.5, colour = "steelblue") +
  #geom_text(nudge_y = 0.02, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
  labs(tag = "d)",
       x = "Mean immigrants",
       y = "Mean λ",
       title = "Population growth rate vs immigration by zone year") +
  theme_bw(base_size = 12)

cor_ps_CrI <- (phia_lambda_p_CrI / prod_lambda_p_CrI | phij_lambda_p_CrI / im_lambda_p_CrI)
cor_ps_CrI
ggsave("figs/cor_ps_CrI.png", cor_ps_CrI, width = 12, height = 10)

#---- ICC bar chart ----

icc_df <- posterior_ipm %>%
  select(icc_phia, icc_prod, icc_im) %>%
  pivot_longer(cols      = everything(),
               names_to  = "rate",
               values_to = "icc") %>%
  group_by(rate) %>%
  summarise(mean  = mean(icc),
            lower = quantile(icc, 0.025),
            upper = quantile(icc, 0.975)) %>%
  mutate(rate = recode(rate,
                       "icc_phia" = "Adult survival",
                       "icc_prod" = "Productivity",
                       "icc_im"   = "Immigration"),
         source = "IPM posterior")

write.csv(icc_df, "icc_table.csv", row.names = FALSE)

lmer_icc <- data.frame(
  rate   = "Population size",
  mean   = 0.373,
  lower  = 0.206,    #no posterior uncertainty from lmer
  upper  = 0.558,
  source = "lmer"
)

icc_df <- bind_rows(icc_df, lmer_icc) %>%
  mutate(rate = factor(rate, 
                       levels = c("Adult survival", 
                                  "Productivity",
                                  "Immigration",
                                  "Population size")))

ICC_barplot <- ggplot(icc_df, aes(x = rate, y = mean, fill = source)) +
  geom_bar(stat   = "identity",
           width  = 0.6,
           colour = "white") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width    = 0.15,
                linewidth = 0.8,
                colour   = "grey30") +
  scale_fill_manual(values = c("IPM posterior" = "steelblue",
                               "lmer"          = "firebrick"),
                    labels = c("IPM posterior (with CrI)",
                               "lmer (with CI)")) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2)) +
  labs(title = "ICC estimates",
       x     = NULL,
       y     = "ICC",
       fill  = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("figs/ICC_barplot.png", ICC_barplot, height = 4, width = 8)

#---- nestbox/sites sample size ----
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


#---- stats for discussion ----
mean(phendat$cs)
csdf <- phendat %>% filter(!is.na(cs))
mean(csdf$cs)
csdf %>% 
  group_by(year) %>% 
  summarise(mean_cs = mean(cs, na.rm = TRUE))

ntot_df %>% 
  group_by(site) %>% 
  reframe(range_ntot = range(mean, na.rm = TRUE))
print(n=24, ntot_df %>% 
        group_by(site) %>% 
        reframe(range_ntot = range(mean, na.rm = TRUE)))

#---- model validation stuff ----
library(bayesplot)

posterior_array <- as.array(ipm_bt)

mcmc_labels <- as_labeller(c('mphia' = "Mean Adult Survival",
                             'mim' = "Mean Immigration",
                             'mprod' = "Mean Productivity",
                             'mlam[1]' = "Mean Lambda (Zone K)"))

mcmc_trace <- mcmc_trace(posterior_array, 
           pars = c("mphia", "mim", "mprod", "mlam[1]"),
           facet_args = list(labeller = mcmc_labels)) +
  theme_bw() +
  labs(x = "Iteration", 
       y = "Parameter Value") 

ggsave("figs/mcmc_trace.png", mcmc_trace, height = 4, width = 8)
