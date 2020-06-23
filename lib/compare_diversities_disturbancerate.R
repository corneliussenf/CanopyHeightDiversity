
library(tidyverse)
library(brms)
library(ggthemes)
library(patchwork)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

colors_landscape <- c("#DDAA33", "#004488")
colors_metrics <- c("#BB5566", "#228833", "#66CCEE")

#### Load and merge data ####

load(file = "data/diversity_boehmerwald.RData")
load(file = "data/diversity_berchtesgaden.RData")

dat <- list(diversity_berchtesgaden, 
            diversity_boehmerwald) %>%
  set_names(c("berchtesgaden", "boehmerwald")) %>%
  bind_rows(.id = "landscape") %>%
  gather(., key = metric, value = diversity, -landscape, -(disturbance_rate:sublandscape)) %>%
  mutate(q = parse_number(metric)) %>%
  mutate(metric = substring(metric, 1, (nchar(as.character(metric)) - 1)))

dat$disturbance_rate_annual <- NA
dat[dat$landscape == "boehmerwald", "disturbance_rate_annual"] <- dat[dat$landscape == "boehmerwald", "disturbance_rate"] / length(1985:2012)
dat[dat$landscape == "berchtesgaden", "disturbance_rate_annual"] <- dat[dat$landscape == "berchtesgaden", "disturbance_rate"] / length(1985:2009)

dat$disturbance_rate_annual_scaled <- dat$disturbance_rate_annual * 10

dat$diversity_log <- log(dat$diversity)
dat <- dat %>% 
  filter(forest_rate > 0.4) # Forest threshold

dat <- dat %>%
  group_by(landscape) %>%
  mutate(elevation_mean_scaled = (elevation_mean - mean(elevation_mean)) / sd(elevation_mean)) %>%
  ungroup() %>% 
  mutate(., landscape = factor(landscape))

#### Calculate averages ###

avg_forest_rate <- dat %>% 
  filter(metric == "alpha" & q == 1) %>% 
  split(.$landscape) %>% 
  map(., ~sample_n(., 100)) %>% 
  bind_rows() %>% 
  .$forest_rate %>% 
  mean(.)

avg_elevation_varcof <- dat %>% 
  filter(metric == "alpha") %>% 
  split(.$landscape) %>% 
  map(., ~sample_n(., 100)) %>% 
  bind_rows() %>% 
  .$elevation_varcof %>% 
  mean(.)

avg_disturbance_rate <- dat %>% 
  filter(metric == "alpha") %>% 
  split(.$landscape) %>% 
  map(., ~sample_n(., 100)) %>% 
  bind_rows() %>% 
  .$disturbance_rate %>% 
  mean(.)

avg_disturbance_rate_annual <- dat %>% 
  filter(metric == "alpha") %>% 
  split(.$landscape) %>% 
  map(., ~sample_n(., 100)) %>% 
  bind_rows() %>% 
  .$disturbance_rate_annual %>% 
  mean(.)

#### Exploration ####

rawplot <- pmap(.l = list(a = dat %>%
                 mutate(landscape = ifelse(landscape == "berchtesgaden", "Berchtesgaden", "Bavarian Forest")) %>%
                 mutate(metric = ifelse(metric == "alpha", "Within-patch", metric)) %>%
                 mutate(metric = ifelse(metric == "beta", "Between-patch", metric)) %>%
                 mutate(metric = ifelse(metric == "gamma", "Overall", metric)) %>%
                 mutate(metric = forcats::fct_relevel(metric, c("Within-patch", "Between-patch", "Overall"))) %>%
                 split(.$q),
               q = list(expression(paste({}^{0}, "D")),
                        expression(paste({}^{1}, "D")),
                        expression(paste({}^{2}, "D"))),
               legendpos = list("none", "right", "none"),
               xtxt = list(NULL, NULL, expression(paste("Disturbance rate (% ", yr^{-1}, ")"))),
               striptxt = list(element_text(size = 12, hjust = 0),
                               element_blank(),
                               element_blank())),
       .f = function(a, q, legendpos, xtxt, striptxt) {
         ggplot(a, aes(x = disturbance_rate_annual, y = diversity, col = landscape, shape = factor(windowsize))) +
           geom_point(alpha = 0.1) +
           facet_wrap(~metric, scales = "free") +
           scale_y_log10() +
           theme_few() +
           scale_fill_manual(values = c(colors_metrics[c(1, 3)], "#BBBBBB"), breaks = c("Bavarian Forest", "Berchtesgaden")) +
           scale_color_manual(values = c(colors_metrics[c(1, 3)], "black"), breaks = c("Bavarian Forest", "Berchtesgaden")) +
           theme(axis.title = element_text(size = 12),
                 axis.text = element_text(size = 9, colour = "darkgrey"),
                 strip.text = striptxt,
                 legend.position = legendpos,
                 legend.title = element_text(size = 9),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 8),
                 panel.spacing.x = unit(1, "lines")) +
           labs(x = xtxt,
                y = q,
                col = "", 
                shape = "Landscape extent (n)") +
           guides(col = guide_legend(ncol = 1,
                                     keywidth = 0.1,
                                     keyheight = 0.1,
                                     default.unit = "inch")) +
           guides(fill = guide_legend(ncol = 1,
                                      keywidth = 0.1,
                                      keyheight = 0.1,
                                      default.unit = "inch")) +
           guides(linetype = guide_legend(ncol = 1,
                                          keywidth = 0.1,
                                          keyheight = 0.1,
                                          default.unit = "inch"))
       }
  ) %>%
  wrap_plots(ncol = 1)

ggsave("figures/rawplot.pdf", rawplot, width = 7.5, height = 6.5)

#### Fit models ####

iteration_names <- cross2(c(0, 1, 2), c("Within", "Between", "Overall")) %>% 
  map(~ paste(unlist(.), collapse = "."))

fits_combined <- vector("list", 3 * 3)

k <- 0

for (m in c("alpha", "beta", "gamma")) {
  
  for (i in c(0, 1, 2)) {
    
    k <- k + 1
    
    print(k)
    
    fit_full_interaction <- brm(bf(diversity_log ~ landscape +
                                     s(disturbance_rate_annual_scaled, by = landscape) +
                                     s(elevation_varcof, by = landscape) + 
                                     s(forest_rate, by = landscape) +
                                     (1 | windowsize)),
                                data = filter(dat, metric == m & q == i), 
                                family = gaussian(), 
                                cores = 4, 
                                control = list(adapt_delta = 0.99, max_treedepth = 11))
    
    loo_full_interaction <- loo(fit_full_interaction)
    
    fit_forest_topo_interaction <- brm(bf(diversity_log ~ landscape +
                                            s(disturbance_rate_annual_scaled) +
                                            s(elevation_varcof, by = landscape) + 
                                            s(forest_rate, by = landscape) +
                                            (1 | windowsize)),
                                       data = filter(dat, metric == m & q == i), 
                                       family = gaussian(), 
                                       cores = 4, 
                                       control = list(adapt_delta = 0.99, max_treedepth = 11))
    
    loo_forest_topo_interaction <- loo(fit_forest_topo_interaction)
    
    fits_combined[[k]] <- list(fit_full_interaction, 
                               loo_full_interaction, 
                               fit_forest_topo_interaction, 
                               loo_forest_topo_interaction) 
    
  }
  
}

save(fits_combined, file = "fits_combined_170620.RData")
load(file = "fits_combined_170620.RData")

#### Choose final model ####

compare_models <- function(inp, m) {
  
  comp <- loo::compare(inp[[2]], inp[[4]])
  
  if (comp[1] + comp[2] * m < 0) {
    return(inp[[1]])
  } else {
    return(inp[[3]])
  }
}

fits_combined_models <- fits_combined %>% 
  map(compare_models, m = 2.5)

formulas <- fits_combined_models %>% map(~ .$formula)

fits_combined %>% 
  map(~ loo_compare(.[[2]], .[[4]]))

#### Model checks ####

# Posterior predictive checks

n_draws <- 30

post_pred <- fits_combined_models %>%
  map(~ posterior_predict(., nsamples = n_draws) %>%
        as.data.frame(.) %>%
        mutate(iter = 1:n_draws) %>%
        gather(key = data_point, value = value, -iter)) %>%
  map2(.y = fits_combined_models, 
       ~ mutate(., windowsize = rep(.y$data$windowsize, each = n_draws))) %>%
  set_names(iteration_names) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("q", "metric"), "\\.") %>%
  mutate(metric = ifelse(metric == "Within", "Within-patch", metric)) %>%
  mutate(metric = ifelse(metric == "Between", "Between-patch", metric)) %>%
  mutate(metric = forcats::fct_relevel(metric, c("Within-patch", "Between-patch", "Overall")))

mod_dat <- fits_combined_models %>%
  map(~ .$data) %>%
  set_names(iteration_names) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("q", "metric"), "\\.") %>%
  mutate(metric = ifelse(metric == "Within", "Within-patch", metric)) %>%
  mutate(metric = ifelse(metric == "Between", "Between-patch", metric)) %>%
  mutate(metric = forcats::fct_relevel(metric, c("Within-patch", "Between-patch", "Overall")))

p <- ggplot() +
  geom_density(data = post_pred, aes(x = exp(value),
                                     group = interaction(iter, q),
                                     fill = q),
               alpha = 0.01, col = NA) +
  geom_density(data = mod_dat, aes(x = exp(diversity_log),
                                   col = q)) +
  facet_wrap(~metric, scales = "free") +
  theme_few() +
  labs(x = "Structural diversity", y = "Density", fill = NULL, color = NULL) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, name = "Dark2"),
                     labels = c(expression(paste({}^{0}, "D")),
                                expression(paste({}^{1}, "D")),
                                expression(paste({}^{2}, "D")))) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, name = "Dark2"),
                    labels = c(expression(paste({}^{0}, "D")),
                               expression(paste({}^{1}, "D")),
                               expression(paste({}^{2}, "D")))) +
  guides(colour = FALSE)

ggsave("figures/posterior_checks.pdf", p, width = 6, height = 2)

# Predicted versus observed

nd <- dat %>% 
  na.omit(.) %>%
  split(list(.$q, .$metric))

post_pred <- fits_combined_models %>%
  map2(.y = nd, ~ posterior_predict(., newdata = .y, nsamples = n_draws) %>%
        as.data.frame() %>%
        mutate(iter = 1:n_draws) %>%
        gather(key = data_point, value = value, -iter) %>%
        mutate(., windowsize = rep(nd$windowsize, each = n_draws)) %>%
        mutate(., landscape = rep(.y$landscape, each = n_draws)) %>%
        mutate(., observed = rep(.y$diversity, each = n_draws))) %>%
  set_names(iteration_names) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("q", "metric"), "\\.") %>%
  mutate(metric = ifelse(metric == "Within", "Within-patch", metric)) %>%
  mutate(metric = ifelse(metric == "Between", "Between-patch", metric)) %>%
  group_by(q, metric, landscape, data_point) %>%
  summarize(predicted = median(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975),
            observed = unique(observed)) %>%
  ungroup(.) %>%
  mutate(landscape = ifelse(landscape == "berchtesgaden", "Berchtesgaden", "Bavarian Forest")) %>%
  mutate(metric = forcats::fct_relevel(metric, c("Within-patch", "Between-patch", "Overall")))

p <- pmap(.l = list(a = post_pred %>% split(.$q),
               b = list(expression(paste("Observed ", {}^{0}, "D")),
                        expression(paste("Observed ", {}^{1}, "D")),
                        expression(paste("Observed ", {}^{2}, "D"))),
               c = list(expression(paste("Predicted ", {}^{0}, "D")),
                        expression(paste("Predicted ", {}^{1}, "D")),
                        expression(paste("Predicted ", {}^{2}, "D"))),
               d = list("right", "none", "none")),
       .f = function(a, b, c, d) {
         ggplot(a, aes(x = observed, y = exp(predicted), col = landscape, shape = landscape)) +
           geom_point(alpha = 0.25) +
           facet_wrap(~metric, scales = "free") +
           geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
           geom_smooth(method = "lm", se = FALSE) +
           theme_few() +
           scale_color_manual(values = c(colors_metrics[c(1, 3)]), 
                              breaks = c("Bavarian Forest", "Berchtesgaden")) +
           scale_shape_manual(values = 1:2, 
                              breaks = c("Bavarian Forest", "Berchtesgaden")) +
           theme(legend.position = d,
                 axis.title = element_text(size = 12),
                 axis.text = element_text(size = 9, colour = "darkgrey"),
                 legend.text = element_text(size = 8),
                 strip.text = element_text(size = 8),
                 panel.spacing.x = unit(1, "lines")) +
           labs(x = b,
                y = c,
                col = NULL,
                shape = NULL)
       }
    ) %>%
  wrap_plots(ncol = 1)
  
ggsave("figures/predicted_versus_observed.pdf", p, width = 7, height = 6.5)

# R-squares

r_squares <- fits_combined_models %>%
  map(bayes_R2) %>%
  map(~ data.frame(mean = median(.), 
                   lower = quantile(., c(0.025)), 
                   upper = quantile(., c(0.975)))) %>%
  set_names(iteration_names) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("q", "metric"), "\\.")

r_squares %>%
  mutate(label = paste0(round(mean, 2), " (", round(lower, 2), "-", round(upper, 2), ")")) %>%
  dplyr::select(q, metric, label) %>%
  spread(key = q, value = label) %>%
  write_csv("results/model_rsquares.csv")

#### Plot disturbance-structural diversity relationships ####

newdata <- vector("list", length = 9)

for (i in 1:length(newdata)) {
  
  if (grepl("landscape", as.character(formulas[[i]])[[1]])) {
    newdat <- expand.grid(landscape = unique(dat$landscape),
                          disturbance_rate_annual_scaled = seq(0, 0.4, length.out = 100),
                          forest_rate = avg_forest_rate,
                          elevation_varcof = avg_elevation_varcof)
    newdat$key <- 1:nrow(newdat)
  } else {
    newdat <- expand.grid(disturbance_rate_annual_scaled = seq(0, 0.4, length.out = 100),
                          forest_rate = avg_forest_rate,
                          elevation_varcof = avg_elevation_varcof)
    newdat$key <- 1:nrow(newdat)
  }
  
  newdata[[i]] <- newdat
    
}

ndraws <- 4000

plotdat <- fits_combined_models %>%
  map2(.y = newdata, 
       ~ posterior_epred(., newdata = .y, nsamples = ndraws, re_formula = NA) %>%
         as.data.frame(.) %>%
         mutate(draw = 1:ndraws) %>%
         gather(., key = key, value = value, -draw) %>%
         mutate(key = as.integer(gsub("V", "", key))) %>%
         group_by(key) %>%
         summarize(pred = median(value),
                   lwr = quantile(value, 0.025),
                   upr = quantile(value, 0.975)) %>%
         mutate(landscape = .y$landscape,
                disturbance_rate_annual_scaled = .y$disturbance_rate_annual_scaled)) %>%
  set_names(iteration_names) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("q", "metric"), "\\.") %>%
  mutate(landscape = ifelse(landscape == "berchtesgaden", "Berchtesgaden", "Bavarian Forest")) %>%
  mutate(metric = ifelse(metric == "Between", "Between-patch", metric)) %>%
  mutate(metric = ifelse(metric == "Within", "Within-patch", metric)) %>%
  mutate(landscape = ifelse(is.na(landscape), "Both", landscape)) %>%
  mutate(metric = forcats::fct_relevel(metric, c("Within-patch", "Between-patch", "Overall")))

dat_plotdat <- dat %>%
  mutate(landscape = ifelse(landscape == "berchtesgaden", "Berchtesgaden", "Bavarian Forest")) %>%
  mutate(metric = ifelse(metric == "alpha", "Within-patch", metric)) %>%
  mutate(metric = ifelse(metric == "beta", "Between-patch", metric)) %>%
  mutate(metric = ifelse(metric == "gamma", "Overall", metric)) %>%
  mutate(metric = forcats::fct_relevel(metric, c("Within-patch", "Between-patch", "Overall")))
  
p_response_curve_disturbances <- pmap(
  .l = list(d = plotdat %>% 
              split(.$q),
            q = list(expression(paste({}^{0}, "D")),
                     expression(paste({}^{1}, "D")),
                     expression(paste({}^{2}, "D"))),
            legendpos = list("right", "none", "none"),
            xtxt = list(NULL, NULL, expression(paste("Disturbance rate (% ", yr^{-1}, ")"))),
            striptxt = list(element_text(size = 12, hjust = 0),
                            element_blank(),
                            element_blank())),
  .f = function(d, maxpoint, q, xtxt, legendpos, striptxt) {
    ggplot(d) +
      geom_ribbon(data = d,
                  aes(x = disturbance_rate_annual_scaled * 10, 
                      ymin = exp(lwr), ymax = exp(upr),
                      fill = landscape),
                  alpha = 0.35) +
      geom_line(data = d,
                aes(x = disturbance_rate_annual_scaled * 10, y = exp(pred), 
                    col = landscape)) +
      facet_wrap(~metric, ncol = 3, scales = "free") +
      theme_few() +
      scale_fill_manual(values = c(colors_metrics[c(1, 3)], "#BBBBBB"), breaks = c("Bavarian Forest", "Berchtesgaden")) +
      scale_color_manual(values = c(colors_metrics[c(1, 3)], "black"), breaks = c("Bavarian Forest", "Berchtesgaden")) +
      theme(axis.title = element_text(size = 12),
            axis.text = element_text(size = 9, colour = "darkgrey"),
            legend.position = legendpos,
            #legend.justification = c(0, 0),
            legend.background = element_blank(),
            legend.text = element_text(size = 8),
            strip.text = striptxt,
            panel.spacing.x = unit(1, "lines")) +
      labs(x = xtxt,
           y = q,
           col = NULL, fill = NULL, linetype = NULL) +
      guides(col = guide_legend(ncol = 1,
                                keywidth = 0.1,
                                keyheight = 0.1,
                                default.unit = "inch")) +
      guides(fill = guide_legend(ncol = 1,
                                keywidth = 0.1,
                                keyheight = 0.1,
                                default.unit = "inch"))
    }
  ) %>%
  wrap_plots(ncol = 1)

ggsave("response_curves_disturbance.pdf", p_response_curve_disturbances, 
       path = "figures/", width = 7, height = 5.5)


#### Scaling among windowsizes ####

newdata <- vector("list", length = 9)

for (i in 1:length(newdata)) {
  
  if (grepl("landscape", as.character(formulas[[i]])[[1]])) {
    newdat <- expand.grid(landscape = unique(dat$landscape),
                          disturbance_rate_annual_scaled = seq(0, 0.4, length.out = 100),
                          forest_rate = avg_forest_rate,
                          elevation_varcof = avg_elevation_varcof,
                          windowsize = unique(dat$windowsize))
    newdat$key <- 1:nrow(newdat)
  } else {
    newdat <- expand.grid(disturbance_rate_annual_scaled = seq(0, 0.4, length.out = 100),
                          forest_rate = avg_forest_rate,
                          elevation_varcof = avg_elevation_varcof,
                          windowsize = unique(dat$windowsize))
    newdat$key <- 1:nrow(newdat)
  }
  
  newdata[[i]] <- newdat
  
}

ndraws <- 4000

plotdat_scaling <- fits_combined_models %>%
  map2(.y = newdata, 
       ~ posterior_epred(., newdata = .y, nsamples = ndraws) %>%
         as.data.frame(.) %>%
         mutate(draw = 1:ndraws) %>%
         gather(., key = key, value = value, -draw) %>%
         mutate(key = as.integer(gsub("V", "", key))) %>%
         group_by(key) %>%
         summarize(pred = median(value),
                   lwr = quantile(value, 0.025),
                   upr = quantile(value, 0.975)) %>%
         mutate(landscape = .y$landscape,
                windowsize = .y$windowsize,
                disturbance_rate_annual_scaled = .y$disturbance_rate_annual_scaled)) %>%
  set_names(iteration_names) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("q", "metric"), "\\.") %>%
  mutate(landscape = ifelse(landscape == "berchtesgaden", "Berchtesgaden", "Bavarian Forest")) %>%
  mutate(metric = ifelse(metric == "Between", "Between-patch", metric)) %>%
  mutate(metric = ifelse(metric == "Within", "Within-patch", metric)) %>%
  mutate(landscape = ifelse(is.na(landscape), "Both", landscape)) %>%
  mutate(metric = forcats::fct_relevel(metric, c("Within-patch", "Between-patch", "Overall")))

p_response_curve_disturbances_scaling <- pmap(
  .l = list(d = plotdat_scaling %>% 
              split(.$q),
            q = list(expression(paste({}^{0}, "D")),
                     expression(paste({}^{1}, "D")),
                     expression(paste({}^{2}, "D"))),
            legendpos = list("none", "right", "none"),
            xtxt = list(NULL, NULL, expression(paste("Disturbance rate (% ", yr^{-1}, ")"))),
            striptxt = list(element_text(size = 12, hjust = 0),
                            element_blank(),
                            element_blank())),
  .f = function(d, maxpoint, q, xtxt, legendpos, striptxt) {
    ggplot(d) +
      geom_ribbon(data = d,
                  aes(x = disturbance_rate_annual_scaled * 10, 
                      ymin = exp(lwr), ymax = exp(upr),
                      fill = landscape, group = interaction(landscape, windowsize)),
                  alpha = 0.1) +
      geom_line(data = d,
                aes(x = disturbance_rate_annual_scaled * 10, y = exp(pred), 
                    col = landscape, linetype = factor(windowsize))) +
      facet_wrap(~metric, ncol = 3, scales = "free") +
      theme_few() +
      scale_fill_manual(values = c(colors_metrics[c(1, 3)], "#BBBBBB"), breaks = c("Bavarian Forest", "Berchtesgaden")) +
      scale_color_manual(values = c(colors_metrics[c(1, 3)], "black"), breaks = c("Bavarian Forest", "Berchtesgaden")) +
      theme(axis.title = element_text(size = 12),
            axis.text = element_text(size = 9, colour = "darkgrey"),
            legend.position = legendpos,
            legend.background = element_blank(),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 9),
            strip.text = striptxt,
            panel.spacing.x = unit(1, "lines")) +
      labs(x = xtxt,
           y = q,
           col = NULL, fill = NULL, linetype = "Landscape extent (n)") +
      guides(col = guide_legend(ncol = 1,
                                keywidth = 0.1,
                                keyheight = 0.1,
                                default.unit = "inch")) +
      guides(fill = guide_legend(ncol = 1,
                                 keywidth = 0.1,
                                 keyheight = 0.1,
                                 default.unit = "inch")) +
      guides(linetype = guide_legend(ncol = 1,
                                 keywidth = 0.1,
                                 keyheight = 0.1,
                                 default.unit = "inch"))
  }
) %>%
  wrap_plots(ncol = 1)

ggsave("response_curves_disturbance_scaling.pdf", p_response_curve_disturbances_scaling, 
       path = "figures/", width = 7, height = 5.5)
