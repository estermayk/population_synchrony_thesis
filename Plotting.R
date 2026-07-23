# Get nyears and nsites
nyears_bt  <- stan_data_bt$nyears
nsites_bt  <- stan_data_bt$nsites
years_bt   <- 1:nyears_bt
years_m1_bt <- 1:(nyears_bt - 1)  # for parameters indexed over nyears-1

site_labels <- paste0("Zone ", unique(adults$zone))

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

# Observed counts per site per year
obs_counts <- expand.grid(site = site_labels, year = survey_years) %>%
  arrange(site, year) %>%
  mutate(observed = as.vector(t(stan_data_bt$y)))  # y is [nsites, nyears]

# Merge both dfs
ntot_df <- left_join(ntot_df, obs_counts, by = c("site", "year"))

mean(ntot_df$mean)

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

p_ntot2 <- ggplot(ntot_df, aes(x = year, group = site, colour = site, fill = site)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_line(aes(y = mean,     linetype = "Estimated Ntot"), linewidth = 0.9) +
  geom_point(aes(y = mean,    shape    = "Estimated Ntot"), size = 1.5) +
 # geom_line(aes(y = observed, linetype = "Observed count"), linewidth = 0.9) +
#  geom_point(aes(y = observed, shape   = "Observed count"), size = 1.5) +
#  scale_linetype_manual(values = c("Estimated Ntot" = "solid",
 #                                  "Observed count"  = "dashed")) +
#  scale_shape_manual(values   = c("Estimated Ntot" = 16,
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

ggsave("figs/pntot2.png", p_ntot2, width = 8, height = 4)

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
  #  scale_linetype_manual(values = c("Estimated Ntot" = "solid",
  #                                  "Observed count"  = "dashed")) +
  #  scale_shape_manual(values   = c("Estimated Ntot" = 16,
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

ggsave("figs/pprod2.png", p_f2, width = 8, height = 4)


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
  geom_line(aes(y = mean.y,     linetype = "Estimated phia"), linewidth = 0.9) +
  geom_point(aes(y = mean.y,    shape    = "Estimated phia"), size = 1.5) +
   geom_line(aes(y = mean.x, linetype = "Estimated phij"), linewidth = 0.9) +
    geom_point(aes(y = mean.x, shape   = "Estimated phij"), size = 1.5) +
    scale_linetype_manual(values = c("Estimated phia" = "solid",
                                    "Estimated phij"  = "dashed")) +
    scale_shape_manual(values   = c("Estimated phia" = 16,
                                   "Estimated phij"  = 1)) +
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

ggsave("figs/pphi2.png", p_phi2, width = 8, height = 4)


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
  #  scale_linetype_manual(values = c("Estimated Ntot" = "solid",
  #                                  "Observed count"  = "dashed")) +
  #  scale_shape_manual(values   = c("Estimated Ntot" = 16,
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
         shape    = guide_legend(nrow = 2))

ggsave("figs/plambda2.png", p_lambda2, width = 8, height = 4)


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
  #  scale_linetype_manual(values = c("Estimated Ntot" = "solid",
  #                                  "Observed count"  = "dashed")) +
  #  scale_shape_manual(values   = c("Estimated Ntot" = 16,
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

ggsave("figs/p_Nadimm2.png", p_Nadimm2, width = 8, height = 4)


rates_plots <- (p_ntot2 | p_lambda2 | p_f2) / (p_phi2 | p_Nadimm2)

rates_plots

ggsave("figs/rates_plot.png", plot = rates_plots)


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

#var_phia_year

var_phia_year <- posterior_ipm %>%
  select(contains("var_phia_year")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_phia_year_p <- ggplot(var_phia_year, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Adult Survival Random Effect Variances",
       x = "phia var",
       y = "Density") +
  theme_minimal()

var_phia_year_p

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

#var_im_year

var_im_year <- posterior_ipm %>%
  select(contains("var_im_year")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_im_year_p <- ggplot(var_im_year, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Immigration Random Effect Variances",
       x = "im var",
       y = "Density") +
  theme_minimal()

var_im_year_p

#prod
prod_icc <- posterior_ipm %>%
  select(contains("icc_prod")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

icc_prod_p <- ggplot(prod_icc, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Productivity ICC",
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
  labs(title = "Productivity Random Effect Variances",
       x = "prod var",
       y = "Density") +
  theme_minimal()

var_prod_year_p

icc_var_plots <- (var_phia_year_p | icc_phia_p) / (var_prod_year_p | icc_prod_p) / (var_im_year_p | icc_im_p)

icc_var_plots

icc_var_phia <- (var_phia_year_p | icc_phia_p)

icc_var_prod <- (var_prod_year_p | icc_prod_p)

icc_var_imm <- (var_im_year_p | icc_im_p)

ggsave("figs/icc_var_plots_bt.png", plot = icc_var_plots, width = 8, height = 10)
ggsave("figs/icc_var_phia.png", plot = icc_var_phia, width = 8, height = 4)
ggsave("figs/icc_var_prod.png", plot = icc_var_prod, width = 8, height = 4)
ggsave("figs/icc_var_imm.png", plot = icc_var_imm, width = 8, height = 4)


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
  labs(x = "Mean adult survival (φa)",
       y = "Mean λ",
       title = "Population growth rate vs adult survival by zone_year") +
  theme_bw(base_size = 12)


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
  labs(x = "Mean productivity",
       y = "Mean λ",
       title = "Population growth rate vs productivity by zone_year") +
  theme_bw(base_size = 12)


ntot_df$site_year <- as.factor(paste(ntot_df$year, ntot_df$site, sep = "_")) 

ntot_siteyear <- ntot_df %>%
  group_by(site_year) %>%
  summarise(mean_ntot = mean,
            lower_ntot = lower,
            upper_ntot = upper,
            year = year)

ntot_lambda <- left_join(lambda_siteyear, ntot_siteyear, by = "site_year")

ntot_lambda <- ntot_lambda %>%
  filter(!year.x %in% c("2014", "2015", "2016"))

ntot_lambda_p <- ggplot(ntot_lambda, aes(x = mean_ntot, y = mean_lambda)) +
  #geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
  #colour = "grey70", width = 0) +
  #geom_errorbarh(aes(xmin = lower_ntot, xmax = upper_ntot), 
  #colour = "grey70", height = 0) +
  geom_point(size = 1.5, colour = "steelblue") +
  #geom_text(nudge_y = 0.02, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
  labs(x = "Mean ntot",
       y = "Mean λ",
       title = "Population growth rate vs population size by zone_year") +
  theme_bw(base_size = 12)

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
  labs(x = "Mean immigrants",
       y = "Mean λ",
       title = "Population growth rate vs immigration by zone_year") +
  theme_bw(base_size = 12)

cor_ps <- (phia_lambda_p / prod_lambda_p | ntot_lambda_p / im_lambda_p)

cor.test(phia_lambda$mean_lambda, phia_lambda$mean_phia)

cor.test(prod_lambda$mean_lambda, prod_lambda$mean_prod)

cor.test(ntot_lambda$mean_lambda, ntot_lambda$mean_ntot)

cor.test(im_lambda$mean_lambda, im_lambda$mean_im)

ggsave("figs/cor_ps.png", cor_ps, width = 12, height = 6)
