# Get nyears and nsites
nyears_bt  <- stan_data_bt$nyears
nsites_bt  <- stan_data_bt$nsites
years_bt   <- 1:nyears_bt
years_m1_bt <- 1:(nyears_bt - 1)  # for parameters indexed over nyears-1

site_labels <- paste0("Zone ", unique(adults$zone))

site_labels

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
  df$year <- years_bt
  df$site <- site_labels[s]
  df
})
ntot_df <- bind_rows(ntot_list)

# Observed counts per site per year
obs_counts <- expand.grid(site = site_labels, year = years_bt) %>%
  arrange(site, year) %>%
  mutate(observed = as.vector(t(stan_data_bt$y)))  # y is [nsites, nyears]

# Merge both dfs
ntot_df <- left_join(ntot_df, obs_counts, by = c("site", "year"))

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

ggplot(ntot_df, aes(x = year, group = site, colour = site, fill = site)) +
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
  df$year <- years_m1_bt
  df$site <- site_labels[s]
  df
})
f_est <- bind_rows(f_list)

# Observed f per site per year
obs_f <- expand.grid(site = site_labels, year = years_m1_bt) %>%
  arrange(site, year) %>%
  mutate(observed = as.vector(t(stan_data_bt$J / stan_data_bt$R)))

# Merge both dfs
f_df <- left_join(f_est, obs_f, by = c("site", "year"))

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


#survival
#juveniles
# Estimated phij per site per year
phij_list <- lapply(1:nsites_bt, function(s) {
  pattern <- paste0("phij\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- years_m1_bt
  df$site <- site_labels[s]
  df
})
phij_est <- bind_rows(phij_list)

#adults
phia_list <- lapply(1:nsites_bt, function(s) {
  pattern <- paste0("phia\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- years_m1_bt
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

#population growth rate

lambda_list <- lapply(1:nsites_bt, function(s) {
  pattern <- paste0("lambda\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- years_m1_bt
  df$site <- site_labels[s]
  df
})
lambda_df <- bind_rows(lambda_list)

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

#immigration
Nadimm_list <- lapply(1:nsites_bt, function(s) {
  pattern <- paste0("Nadimm\\[", s, ",")
  df <- extract_summary(posterior_ipm, pattern)
  df$year <- years_bt
  df$site <- site_labels[s]
  df
})
Nadimm_df <- bind_rows(Nadimm_list)

# Observed counts per site per year
obs_counts <- expand.grid(site = site_labels, year = years_bt) %>%
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

p_Nadimm

rates_plots <- (p_ntot | p_lambda | p_f) / (p_phi | p_Nadimm)

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

#var_phia_siteyear

var_phia_siteyear <- posterior_ipm %>%
  select(contains("var_phia_siteyear")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_phia_siteyear_p <- ggplot(var_phia_siteyear, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Adult Survival Random Effect Variances",
       x = "phia var",
       y = "Density") +
  theme_minimal()

var_phia_siteyear_p

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

#var_im_siteyear

var_im_siteyear <- posterior_ipm %>%
  select(contains("var_im_siteyear")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_im_siteyear_p <- ggplot(var_im_siteyear, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Immigration Random Effect Variances",
       x = "im var",
       y = "Density") +
  theme_minimal()

var_im_siteyear_p

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

#var_prod_siteyear

var_prod_siteyear <- posterior_ipm %>%
  select(contains("var_prod_siteyear")) %>%
  pivot_longer(cols = everything(), values_to = "Value")

var_prod_siteyear_p <- ggplot(var_prod_siteyear, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Productivity Random Effect Variances",
       x = "prod var",
       y = "Density") +
  theme_minimal()

var_prod_siteyear_p

icc_var_plots <- (var_phia_siteyear_p | icc_phia_p) / (var_prod_siteyear_p | icc_prod_p) / (var_im_siteyear_p | icc_im_p)

icc_var_plots

ggsave("figs/icc_var_plots_bt.png", plot = icc_var_plots)

lambda_site <- lambda_df %>%
  group_by(site) %>%
  summarise(mean_lambda = mean(mean),
            lower_lambda = mean(lower),
            upper_lambda = mean(upper))

phia_site <- phia_est %>%
  group_by(site) %>%
  summarise(mean_phia = mean(mean),
            lower_phia = mean(lower),
            upper_phia = mean(upper))

phia_lambda <- left_join(lambda_site, phia_site, by = "site")

ggplot(phia_lambda, aes(x = mean_phia, y = mean_lambda, label = site)) +
  #geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
                #colour = "grey70", width = 0) +
  #geom_errorbarh(aes(xmin = lower_phia, xmax = upper_phia), 
                 #colour = "grey70", height = 0) +
  geom_point(size = 3, colour = "steelblue") +
  geom_text(nudge_y = 0.02, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
  labs(x = "Mean adult survival (φa)",
       y = "Mean λ",
       title = "Population growth rate vs adult survival by zone") +
  theme_bw(base_size = 12)


prod_site <- f_est %>%
  group_by(site) %>%
  summarise(mean_prod = mean(mean),
            lower_prod = mean(lower),
            upper_prod = mean(upper))

prod_lambda <- left_join(lambda_site, prod_site, by = "site")

ggplot(prod_lambda, aes(x = mean_prod, y = mean_lambda, label = site)) +
  #geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
  #colour = "grey70", width = 0) +
  #geom_errorbarh(aes(xmin = lower_prod, xmax = upper_prod), 
  #colour = "grey70", height = 0) +
  geom_point(size = 3, colour = "steelblue") +
  geom_text(nudge_y = 0.02, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
  labs(x = "Mean productivity",
       y = "Mean λ",
       title = "Population growth rate vs productivity by zone") +
  theme_bw(base_size = 12)


ntot_site <- ntot_df %>%
  group_by(site) %>%
  summarise(mean_ntot = mean(mean),
            lower_ntot = mean(lower),
            upper_ntot = mean(upper))

ntot_lambda <- left_join(lambda_site, ntot_site, by = "site")

ggplot(ntot_lambda, aes(x = mean_ntot, y = mean_lambda, label = site)) +
  #geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
  #colour = "grey70", width = 0) +
  #geom_errorbarh(aes(xmin = lower_ntot, xmax = upper_ntot), 
  #colour = "grey70", height = 0) +
  geom_point(size = 3, colour = "steelblue") +
  geom_text(nudge_y = 0.02, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
  labs(x = "Mean ntot",
       y = "Mean λ",
       title = "Population growth rate vs population size by zone") +
  theme_bw(base_size = 12)

im_site <- Nadimm_df %>%
  group_by(site) %>%
  summarise(mean_im = mean(mean),
            lower_im = mean(lower),
            upper_im = mean(upper))

im_lambda <- left_join(lambda_site, im_site, by = "site")

ggplot(im_lambda, aes(x = mean_im, y = mean_lambda, label = site)) +
  #geom_errorbar(aes(ymin = lower_lambda, ymax = upper_lambda), 
  #colour = "grey70", width = 0) +
  #geom_errorbarh(aes(xmin = lower_im, xmax = upper_im), 
  #colour = "grey70", height = 0) +
  geom_point(size = 3, colour = "steelblue") +
  geom_text(nudge_y = 0.02, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", alpha = 0.2) +
  labs(x = "Mean immigrants",
       y = "Mean λ",
       title = "Population growth rate vs immigration by zone") +
  theme_bw(base_size = 12)


cor.test(phia_lambda$mean_lambda, phia_lambda$mean_phia)

cor.test(prod_lambda$mean_lambda, prod_lambda$mean_prod)

cor.test(ntot_lambda$mean_lambda, ntot_lambda$mean_ntot)

cor.test(im_lambda$mean_lambda, im_lambda$mean_im)
