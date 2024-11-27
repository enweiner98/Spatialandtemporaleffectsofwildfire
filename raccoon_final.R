#### RACCOON FINAL ANALYSIS #####
# single-species, multi-factor models
# approach: temporal & spatial global models

#### LIBRARIES ####
library(tidyverse)
library(caret)
library(stats)
library(autoOcc)
library(beepr)

#### OCCUPANCY DATA ####
# temporal
all_occu <- read.csv("allOccu_weekly_092723.csv") # complete pre/post dataset

all_occu <- dplyr::distinct(all_occu)  # remove duplicate sites/seasons

all_occu <- split(all_occu,
                  all_occu$Species)  # create a list split by species 

# spatial 
post_occu <- read.csv("postOccu_weekly_100223.csv")  # post only dataset v3

post_occu <- dplyr::distinct(post_occu)

post_occu <- split(post_occu,
                   post_occu$Species)

#### RACCOON DATAFRAMES ####
r_all <- all_occu$raccoon  # all raccoon detections
r_post <- post_occu$raccoon  # post-fire raccoon detections

## summary statistics
# temporal: total detections
r_all %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 64 total detections

# temporal: detections at burned v. unburned sites 
r_bs <- r_all
r_bs$burn_status <- var_cov$burned

bs_count <- r_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 24 detections 
# unburned: 13 detections 
# prefire: 27 detections

# temporal: total sites 
tot_sites <- r_all %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 22 sites 

# spatial: total detections
r_post %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 4 total detections

# detections at burned v. unburned sites 
r_bs <- r_post
r_bs$burn_status <- post_varCov$burned

bs_count <- r_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 0 detections 
# unburned: 4 detections 

# temporal: total sites 
tot_sites <- r_post %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 4 sites 

## creating dataframes for autoOcc
# post-fire only
raccoon_S <- format_y(
  x = r_post,
  site_column = "Site",
  time_column = "Season",
  history_column = "Week"
)

# includes pre-fire and post-fire data
raccoon_T <- format_y(        
  x = r_all,
  site_column = "Site",
  time_column = "Season",
  history_columns = "Week" 
)

#### SPATIAL ANALYSIS ####
#### null model ####
racnS_null <- auto_occ(
  formula = ~1~1,
  y = raccoon_S
)

# overall expected occupancy: 
(
  intercept_preds_psi <- predict(
    racnS_null,
    type = "psi"
  )
)

# overall expected detection: 
(
  intercept_preds_rho <- predict(
    racnS_null, 
    type = "rho"
  )
)

#### p1: fire heterogeneity ####
racnS_het <- auto_occ(
  formula = ~onroad + burn_status * rdnbr_het
  ~burn_status * rdnbr_het,
  y = raccoon_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(opoS_het)

#### p2: fire severity ####
racnS_rdnbr <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad
  ~burn_status + rdnbr + rdnbrQuad,
  y = raccoon_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(racnS_rdnbr)

#### p3: distance to fire perimeter ####
racnS_dp <- auto_occ(
  formula = ~onroad + burn_status * dist_perim
  ~burn_status * dist_perim,
  y = raccoon_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(opoS_dp)

#### TEMPORAL ANALYSIS ####

#### null model ####
racnT_null <- auto_occ(
  formula = ~1~1,
  y = raccoon_T
)

# overall expected occupancy: 0.236
(
  intercept_preds_psi <- predict(
    racnT_null,
    type = "psi"
  )
)

# overall expected detection: 0.067
(
  intercept_preds_rho <- predict(
    racnT_null, 
    type = "rho"
  )
)

#### global model ####
racnT_global <- auto_occ(
  formula = ~onroad + burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi
  ~burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi,
  y = raccoon_T,
  det_covs = covFrame_T,
  occ_covs = covFrame_T
)
summary(racnT_global)

#### TEMPORAL MODEL SELECTION ####
compare_models(list(racnT_null, racnT_global), digits = 2)

#### TEMPORAL MODEL PREDICTIONS ####

#### burn status ####
burned_seq <- c("prefire", "burned_post", "unburned_post")
season_levels <- levels(cov_frame$season$V1)

racn_burn_real <- data.frame(
  burn_status = factor(burned_seq, 
                       levels = burned_seq),
  tsf = -3,
  ndvi = 0, 
  onroad = 0
)

# the model prediction
burn_racn <- predict(
  object = racnT_global,
  type = "rho",
  newdata = racn_burn_real
)

# add on the covariate data
burn_racn <- data.frame(
  burn_racn,
  racn_burn_real
)

# plot
r2 <- ggplot(burn_racn, aes(x = factor(burn_status,
                                levels = c("prefire", "burned_post", "unburned_post")), 
                     y = estimate,
                     color = burn_status)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, show.legend = FALSE) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_vline(aes(xintercept = 1.5), color = "gray40", linetype = "longdash", size = 1) +
  labs(x = "Burn Status", y = "Detection") +
  scale_x_discrete(labels = c("Pre-Fire", "Burned", "Unburned")) +
  scale_color_manual(values = c("darkgreen", "darkorange", "darkgreen")) +
  ylim(0, 1) +
  theme_classic(18)

r1 # occupancy
r2 # detection

#### burn_status*tsf ####
racn_int <- data.frame(
  burn_status = c(
    rep("prefire", 23),
    rep("burned_post", 23),
    rep("unburned_post", 23)
  ),
  tsf = c(
    -22:0,
    0:22, 
    0:22
  ),
  ndvi = 0,
  onroad = 0
)
racn_int$burn_status <- factor(racn_int$burn_status,
                              levels = c("prefire", "burned_post", "unburned_post"))

racn_int_pred <- predict(racnT_global,
                        type = 'rho',
                        newdata = racn_int,
                        interval = "confidence")

racn_int_pred$tsf <- c(-22:0, 0:22, 0:22)
racn_int_pred$burn_status <- racn_int$burn_status

racn_int_pred <- split(racn_int_pred, racn_int_pred$burn_status)

#plot
plot(1~1, type = "n", xlab = "Seasons since fire", ylab = "Detection",
     xlim = c(-22, 23), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(racn_int_pred$prefire$tsf, rev(racn_int_pred$prefire$tsf)),
  y = c(racn_int_pred$prefire$lower, rev(racn_int_pred$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = racn_int_pred$prefire$tsf, y = racn_int_pred$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(racn_int_pred$burned_post$tsf, rev(racn_int_pred$burned_post$tsf)),
  y = c(racn_int_pred$burned_post$lower, rev(racn_int_pred$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = racn_int_pred$burned_post$tsf, y = racn_int_pred$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(racn_int_pred$unburned_post$tsf, rev(racn_int_pred$unburned_post$tsf)),
  y = c(racn_int_pred$unburned_post$lower, rev(racn_int_pred$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = racn_int_pred$unburned_post$tsf, y = racn_int_pred$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topleft", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

r3 <- recordPlot() # occupancy
r4 <- recordPlot() # detection

#### ndvi ####
ndvi_seq <- seq(0, 0.30, 0.01)

racn_ndvi_real <- data.frame(
  burn_status = factor("burned_post",
                       levels = c("prefire", "burned_post", "unburned_post")
  ),
  tsf = -3,
  ndvi = ndvi_seq,
  onroad = 0
)

# scaling
racn_ndvi_scaled <- racn_ndvi_real
racn_ndvi_scaled$ndvi <- (
  racn_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_racn <- predict(
  object = racnT_global,
  type = "psi",
  newdata = racn_ndvi_scaled
)

# add on covariate data
ndvi_racn <- data.frame(
  ndvi_racn,
  racn_ndvi_real
)

# plot
ggplot(ndvi_racn, aes(x = ndvi, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Vegetation Greenness", y = "Detection", title = "Raccoon") + 
  ylim(0,1) +
  theme_classic(18)

#### burn_status*ndvi ####
ndvi_seq <- seq(0, 0.35, 0.01)

racn_ndvi_real <- data.frame(
  burn_status = c(
    rep("prefire", 36),
    rep("burned_post", 36),
    rep("unburned_post", 36)
  ),
  tsf = -3,
  ndvi = c(
    ndvi_seq,
    ndvi_seq,
    ndvi_seq
  ),
  onroad = 0
)

# scaling
racn_ndvi_scaled <- racn_ndvi_real
racn_ndvi_scaled$ndvi <- (
  racn_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_racn <- predict(
  object = racnT_global,
  type = "rho",
  newdata = racn_ndvi_scaled
)

# add on covariate data
ndvi_racn <- data.frame(
  ndvi_racn,
  racn_ndvi_real
)

ndvi_racn <- split(ndvi_racn, ndvi_racn$burn_status)

#plot
plot(1~1, type = "n", xlab = "NDVI", ylab = "Detection",
     xlim = c(0, 0.3), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(ndvi_racn$prefire$ndvi, rev(ndvi_racn$prefire$ndvi)),
  y = c(ndvi_racn$prefire$lower, rev(ndvi_racn$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = ndvi_racn$prefire$ndvi, y = ndvi_racn$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(ndvi_racn$burned_post$ndvi, rev(ndvi_racn$burned_post$ndvi)),
  y = c(ndvi_racn$burned_post$lower, rev(ndvi_racn$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = ndvi_racn$burned_post$ndvi, y = ndvi_racn$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(ndvi_racn$unburned_post$ndvi, rev(ndvi_racn$unburned_post$ndvi)),
  y = c(ndvi_racn$unburned_post$lower, rev(ndvi_racn$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = ndvi_racn$unburned_post$ndvi, y = ndvi_racn$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topleft", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

r5 <- recordPlot() # occupancy
r6 <- recordPlot() # detection

#### SUMMARY TABLES ####

#### temporal ####
parameter <- c(
  paste(greeks("psi"), "- Intercept"), paste(greeks("psi"), "- Burn Status (Unburned)"), 
  paste(greeks("psi"), "- Burn Status (Pre-Fire)"), paste(greeks("psi"), "- Time Since Fire"),
  paste(greeks("psi"), "- NDVI"), paste(greeks("psi"), "- Unburned x Time Since Fire"), 
  paste(greeks("psi"), "- Pre-Fire x Time Since Fire"), paste(greeks("psi"), "- Unburned x NDVI"),
  paste(greeks("psi"), "- Pre-Fire x NDVI"), paste(greeks("psi"), "-", greeks("theta")), 
  paste(greeks("rho"), "- Intercept"), paste(greeks("rho"), "- Microsite Attractant"), paste(greeks("rho"), "- Burn Status (Unburned)"), 
  paste(greeks("rho"), "- Burn Status (Pre-Fire)"), paste(greeks("rho"), "- Time Since Fire"),
  paste(greeks("rho"), "- NDVI"), paste(greeks("rho"), "- Unburned x Time Since Fire"), 
  paste(greeks("rho"), "- Pre-Fire x Time Since Fire"), paste(greeks("rho"), "- Unburned x NDVI"),
  paste(greeks("rho"), "- Pre-Fire x NDVI")
)

estimate <- racnT_global@estimates$Est
estimate <- round(estimate, 3)

SE <- racnT_global@estimates$SE
SE <- round(SE, 3)

p <- racnT_global@estimates$p
p <- round(p, 3)

racnGlobal_tab <- bind_cols(parameter, estimate, SE, p)
racnGlobal_tab <- racnGlobal_tab %>% 
  rename("Parameter" = ...1, 
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
racnGlobal_tab <- racnGlobal_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

racn_g <- racnGlobal_tab %>% 
  gt() %>% 
  cols_align(align = "left", 
             columns = "Parameter") %>%
  tab_style(style = cell_text(weight = "normal"),
            locations = cells_column_labels(c("Parameter",
                                              "Estimate",
                                              "SE",
                                              "p"))) %>% 
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(
              columns = p,
              rows = p < 0.05
            )) %>% 
  tab_spanner(label = md("Global Model"),
              columns = everything())

racn_g %>% 
  gtsave(filename = "racn_g.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")




