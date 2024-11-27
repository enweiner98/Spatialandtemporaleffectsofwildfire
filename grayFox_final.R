#### GRAY FOX FINAL ANALYSIS #####
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

#### GRAY FOX DATAFRAMES ####
gf_all <- all_occu$gray_fox  # all gray fox detections
gf_post <- post_occu$gray_fox  # post-fire gray fox detections

#### summary statistics ####
# temporal: total detections
total_count <- gf_all %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE)))
rowSums(total_count[,]) # 1078 total detections

# temporal: detections at burned v. unburned sites 
gf_bs <- gf_all
gf_bs$burn_status <- var_cov$burned

bs_count <- gf_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm = TRUE)))
rowSums(bs_count[,2:6])
# burned: 249 detections
# unburned: 206 detections
# prefire: 623 detections

# temporal: total sites 
tot_sites <- gf_all %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 63 sites 

# spatial: total detections
gf_post %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 370 total detections

# spatial: detections by burn status
gf_bs <- gf_post
gf_bs$burn_status <- post_varCov$burned

bs_count <- gf_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 212 detections 
# unburned: 158 detections

# spatial: total sites 
tot_sites <- gf_post %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 49 sites 

## creating dataframes for autoOcc
# post-fire only
grayFox_S <- format_y(
  x = gf_post,
  site_column = "Site",
  time_column = "Season",
  history_column = "Week"
)

# includes pre-fire and post-fire data
grayFox_T <- format_y(        
  x = gf_all,
  site_column = "Site",
  time_column = "Season",
  history_columns = "Week" 
)

#### SPATIAL ANALYSIS ####
#### null model ####
gfS_null <- auto_occ(
  formula = ~1~1,
  y = grayFox_S
)

# overall expected occupancy: 0.413
(
  intercept_preds_psi <- predict(
    gfS_null,
    type = "psi"
  )
)

# overall expected detection: 0.510
(
  intercept_preds_rho <- predict(
    gfS_null, 
    type = "rho"
  )
)

#### p1: fire heterogeneity ####
gfS_het <- auto_occ(
  formula = ~onroad + burn_status + rdnbr_het + burn_status:rdnbr_het
  ~burn_status + rdnbr_het + burn_status:rdnbr_het,
  y = grayFox_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(gfS_het)

#### p2: fire severity ####
gfS_rdnbr <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad
  ~burn_status + rdnbr + rdnbrQuad,
  y = grayFox_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(gfS_rdnbr)

#### p3: distance to fire perimeter ####
gfS_dp <- auto_occ(
  formula = ~onroad + burn_status + dist_perim + burn_status:dist_perim
  ~burn_status + dist_perim + burn_status:dist_perim,
  y = grayFox_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(gfS_dp)

#### 2 factor: severity ####
gfS_sev_dp <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim
  ~burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim,
  y = grayFox_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(gfS_sev_dp)

#### 2 factor: heterogeneity ####
gfS_het_dp <- auto_occ(
  formula = ~onroad + burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim
  ~burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim,
  y = grayFox_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(gfS_het_dp)

#### SPATIAL MODEL SELECTION ####
compare_models(list(gfS_null, gfS_het, gfS_rdnbr, gfS_dp, gfS_sev_dp, gfS_het_dp),
               digits = 2)  # best fit model: m2 - fire heterogeneity

#### SPATIAL MODEL PREDICTIONS ####

#### fire heterogeneity ####
gf_rdnbrHet_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  
  rdnbr_het = rdnbrHet,
  dist_perim = 0,
  onroad = 0
)

# scaling
gf_rdnbrHet_scaled <- gf_rdnbrHet_real
gf_rdnbrHet_scaled$rdnbr_het <- (
  gf_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_gf <- predict(
  object = gfS_het,
  type = "rho",
  newdata = gf_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_gf <- data.frame(
  rdnbrHet_gf,
  gf_rdnbrHet_real
)

# plot 
ggplot(rdnbrHet_gf, aes(rdnbr_het, estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Fire Heterogeneity", y = "Detection") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*fire heterogeneity ####
gf_rdnbrHet_real <- data.frame(
  burn_status = c(
    rep("burned", 174),
    rep("unburned", 174)
  ),
  rdnbr_het = c(
    rdnbrHet,
    rdnbrHet
  ),
  dist_perim = 0,
  onroad = 0
)

# scaling
gf_rdnbrHet_scaled <- gf_rdnbrHet_real
gf_rdnbrHet_scaled$rdnbr_het <- (
  gf_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_gf <- predict(
  object = gfS_het,
  type = "rho",
  newdata = gf_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_gf <- data.frame(
  rdnbrHet_gf,
  gf_rdnbrHet_real
)

rdnbrHet_gf <- split(rdnbrHet_gf, rdnbrHet_gf$burn_status)

#plot
plot(1~1, type = "n", xlab = "Fire Heterogeneity", ylab = "Detection",
     xlim = c(20, 1750), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(rdnbrHet_gf$burned$rdnbr_het, rev(rdnbrHet_gf$burned$rdnbr_het)),
  y = c(rdnbrHet_gf$burned$lower, rev(rdnbrHet_gf$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = rdnbrHet_gf$burned$rdnbr_het, y = rdnbrHet_gf$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(rdnbrHet_gf$unburned$rdnbr_het, rev(rdnbrHet_gf$unburned$rdnbr_het)),
  y = c(rdnbrHet_gf$unburned$lower, rev(rdnbrHet_gf$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = rdnbrHet_gf$unburned$rdnbr_het, y = rdnbrHet_gf$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topleft", c("Unburned", "Burned"), lwd = 2, lty = c(2, 1), col = c("darkgreen", "darkorange"), bty = "n", cex=1.25)

gf1 <- recordPlot() # occupancy
gf2 <- recordPlot() # detection

#### fire severity ####
gf_rdnbr_real <- data.frame(
  burn_status = factor("burned",
                      levels = c("unburned", "burned")
  ),
  rdnbr = rdnbr_seq,
  rdnbrQuad = 0,
  dist_perim = 0,
  onroad = 0
)

# scaling
gf_rdnbr_scaled <- gf_rdnbr_real
gf_rdnbr_scaled$rdnbr <- (
  gf_rdnbr_scaled$rdnbr - mean(post_varCov$rdnbr_abs)
) / sd(post_varCov$rdnbr_abs)

# the model prediction
rdnbr_gf <- predict(
  object = gfS_sev_dp,
  type = "rho",
  newdata = gf_rdnbr_scaled
)

# add on covariate data
rdnbr_gf <- data.frame(
  rdnbr_gf,
  gf_rdnbr_real
)

ggplot(rdnbr_gf, aes(x = rdnbr, y = estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Fire Severity (RdNBR)", y = "Detection") + 
  ylim(0,1) +
  theme_classic(18)

gf3 # occupancy
gf4 # detection

#### distance to fire perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)
dpQuad <- seq(-0.60, 3.4, 0.1)

gf_dp_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  dist_perim = dp_seq,
  onroad = 0
)

# scaling
gf_dp_scaled <- gf_dp_real
gf_dp_scaled$dist_perim <- (
  gf_dp_scaled$dist_perim - mean(post_varCov$dist_perim)
) / sd(post_varCov$dist_perim)

# the model prediction
dp_gf <- predict(
  object = gfS_dp,
  type = "rho",
  newdata = gf_dp_real
)

# add on covariate data
dp_gf <- data.frame(
  dp_gf,
  gf_dp_real
)

ggplot(dp_gf, aes(x = dist_perim, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Distance to Fire Perimeter (km)", y = "Occupancy") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*distance to perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

gf_dp_real <- data.frame(
  burn_status = c(
    rep("burned", 7),
    rep("unburned", 10)
  ),
  dist_perim = c(
    -6:0,
    0:9
  ),
  rdnbr = 0,
  rdnbrQuad = 0,
  onroad = 0
)

# the model prediction
dp_gf <- predict(
  object = gfS_sev_dp,
  type = "rho",
  newdata = gf_dp_real
)

dp_gf$dist_perim <- c(-6:0, 0:9)
dp_gf$burn_status <- gf_dp_real$burn_status

dp_gf <- split(dp_gf, gf_dp_real$burn_status)

#plot
plot(1~1, type = "n", xlab = "Distance to Fire Perimeter (km)", ylab = "Detection",
     xlim = c(-5, 5), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(dp_gf$burned$dist_perim, rev(dp_gf$burned$dist_perim)),
  y = c(dp_gf$burned$lower, rev(dp_gf$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = dp_gf$burned$dist_perim, y = dp_gf$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(dp_gf$unburned$dist_perim, rev(dp_gf$unburned$dist_perim)),
  y = c(dp_gf$unburned$lower, rev(dp_gf$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = dp_gf$unburned$dist_perim, y = dp_gf$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topright", c("Burned", "Unburned", "Fire Perimeter"), lwd = 2, lty = c(1, 2, 1), col = c("darkorange", "darkgreen", "black"), bty = "n", cex=1.25)

gf5 <- recordPlot() # occupancy
gf6 <- recordPlot() # detection

#### TEMPORAL ANALYSIS ####

#### null model ####
gfT_null <- auto_occ(
  formula = ~1~1,
  y = grayFox_T
)

# overall expected occupancy: 0.391
(
  intercept_preds_psi <- predict(
    gfT_null,
    type = "psi"
  )
)

# overall expected detection: 0.526
(
  intercept_preds_rho <- predict(
    gfT_null, 
    type = "rho"
  )
)

#### global model ####
gfT_global <- auto_occ(
  formula = ~onroad + burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi
  ~burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi,
  y = grayFox_T,
  det_covs = covFrame_T,
  occ_covs = covFrame_T
)
summary(gfT_global)

#### TEMPORAL MODEL SELECTION ####
compare_models(list(gfT_null, gfT_global), digits = 2)

#### TEMPORAL MODEL PREDICTIONS ####

#### burn status ####
burned_seq <- c("prefire", "burned_post", "unburned_post")
season_levels <- levels(cov_frame$season$V1)

gf_burn_real <- data.frame(
  burn_status = factor(burned_seq, 
                       levels = burned_seq),
  tsf = -3,
  ndvi = 0, 
  ndvi_het = 0,
  onroad = 0
)

# the model prediction
burn_gf <- predict(
  object = gfT_global,
  type = "rho",
  newdata = gf_burn_real
)

# add on the covariate data
burn_gf <- data.frame(
  burn_gf,
  gf_burn_real
)

# plot
gf8 <- ggplot(burn_gf, aes(x = factor(burn_status,
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

gf7 # occupancy
gf8 # detection

#### burn_status*tsf ####
gf_int <- data.frame(
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
  ndvi_het = 0,
  onroad = 0
)
gf_int$burn_status <- factor(gf_int$burn_status,
                               levels = c("prefire", "burned_post", "unburned_post"))

gf_int_pred <- predict(gfT_global,
                         type = 'rho',
                         newdata = gf_int,
                         interval = "confidence")

gf_int_pred$tsf <- c(-22:0, 0:22, 0:22)
gf_int_pred$burn_status <- gf_int$burn_status

gf_int_pred <- split(gf_int_pred, gf_int_pred$burn_status)

#plot
plot(1~1, type = "n", xlab = "Seasons since fire", ylab = "Detection",
     xlim = c(-22, 23), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(gf_int_pred$prefire$tsf, rev(gf_int_pred$prefire$tsf)),
  y = c(gf_int_pred$prefire$lower, rev(gf_int_pred$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = gf_int_pred$prefire$tsf, y = gf_int_pred$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(gf_int_pred$burned_post$tsf, rev(gf_int_pred$burned_post$tsf)),
  y = c(gf_int_pred$burned_post$lower, rev(gf_int_pred$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = gf_int_pred$burned_post$tsf, y = gf_int_pred$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(gf_int_pred$unburned_post$tsf, rev(gf_int_pred$unburned_post$tsf)),
  y = c(gf_int_pred$unburned_post$lower, rev(gf_int_pred$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = gf_int_pred$unburned_post$tsf, y = gf_int_pred$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topright", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

gf9 <- recordPlot() # occupancy
gf10 <- recordPlot() # detection

#### ndvi ####
ndvi_seq <- seq(0, 0.30, 0.01)

gf_ndvi_real <- data.frame(
  burn_status = factor("burned_post",
                       levels = c("prefire", "burned_post", "unburned_post")
  ),
  tsf = -3,
  ndvi = ndvi_seq,
  ndvi_het = 0,
  onroad = 0
)

# scaling
gf_ndvi_scaled <- gf_ndvi_real
gf_ndvi_scaled$ndvi <- (
  gf_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_gf <- predict(
  object = gfT_global,
  type = "psi",
  newdata = gf_ndvi_scaled
)

# add on covariate data
ndvi_gf <- data.frame(
  ndvi_gf,
  gf_ndvi_real
)

# plot
ggplot(ndvi_gf, aes(x = ndvi, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "NDVI", y = "Occupancy") + 
  ylim(0,1) +
  theme_classic(18)

#### burn_status*ndvi ####
ndvi_seq <- seq(0, 0.35, 0.01)

gf_ndvi_real <- data.frame(
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
gf_ndvi_scaled <- gf_ndvi_real
gf_ndvi_scaled$ndvi <- (
  gf_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)
view(gf_ndvi_scaled)

# the model prediction
ndvi_gf <- predict(
  object = gfT_global,
  type = "psi",
  newdata = gf_ndvi_scaled
)

# add on covariate data
ndvi_gf <- data.frame(
  ndvi_gf,
  gf_ndvi_real
)

ndvi_gf <- split(ndvi_gf, ndvi_gf$burn_status)

#plot
plot(1~1, type = "n", xlab = "Vegetation Biomass", ylab = "Occupancy",
     xlim = c(0, 0.3), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(ndvi_gf$prefire$ndvi, rev(ndvi_gf$prefire$ndvi)),
  y = c(ndvi_gf$prefire$lower, rev(ndvi_gf$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = ndvi_gf$prefire$ndvi, y = ndvi_gf$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(ndvi_gf$burned_post$ndvi, rev(ndvi_gf$burned_post$ndvi)),
  y = c(ndvi_gf$burned_post$lower, rev(ndvi_gf$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = ndvi_gf$burned_post$ndvi, y = ndvi_gf$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(ndvi_gf$unburned_post$ndvi, rev(ndvi_gf$unburned_post$ndvi)),
  y = c(ndvi_gf$unburned_post$lower, rev(ndvi_gf$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = ndvi_gf$unburned_post$ndvi, y = ndvi_gf$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topright", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

gf11 <- recordPlot() # occupancy
gf12 <- recordPlot() # detection

#### not using - heterogeneity ####
ndviHet_seq <- seq(0, 0.15, 0.01)

gf_ndviHet_real <- data.frame(
  burn_status = factor("burned_post",
                       levels = c("prefire", "burned_post", "unburned_post")
  ),
  tsf = -3,
  ndvi = 0,
  ndvi_het = ndviHet_seq,
  onroad = 0
)

# scaling
gf_ndviHet_scaled <- gf_ndviHet_real
gf_ndviHet_scaled$ndvi_het <- (
  gf_ndviHet_scaled$ndvi_het - mean(var_cov$ndvi_het)
) / sd(var_cov$ndvi_het)

# the model prediction
ndviHet_gf <- predict(
  object = gfT_global,
  type = "rho",
  newdata = gf_ndviHet_scaled
)

# add on covariate data
ndviHet_gf <- data.frame(
  ndviHet_gf,
  gf_ndviHet_real
)

# plot
ggplot(ndviHet_gf, aes(x = ndvi_het, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Vegetation Heterogeneity", y = "Detection Probability", title = "Gray Fox") + 
  ylim(0,1) +
  theme_classic(18)

#### not using - burn_status*heterogeneity ####
ndviHet_seq <- seq(0, 0.15, 0.01)

gf_ndviHet_real <- data.frame(
  burn_status = c(
    rep("prefire", 16),
    rep("burned_post", 16),
    rep("unburned_post", 16)
  ),
  tsf = -3,
  ndvi = 0,
  ndvi_het = c(
    ndviHet_seq,
    ndviHet_seq,
    ndviHet_seq
  ),
  onroad = 0
)

# scaling
gf_ndviHet_scaled <- gf_ndviHet_real
gf_ndviHet_scaled$ndvi_het <- (
  gf_ndviHet_scaled$ndvi_het - mean(var_cov$ndvi_het)
) / sd(var_cov$ndvi_het)
view(gf_ndviHet_scaled)

# the model prediction
ndviHet_gf <- predict(
  object = gfT_global,
  type = "rho",
  newdata = gf_ndviHet_scaled
)

# add on covariate data
ndviHet_gf <- data.frame(
  ndviHet_gf,
  gf_ndviHet_real
)

ndviHet_gf <- split(ndviHet_gf, ndviHet_gf$burn_status)

plot(1~1, type = "n", xlab = "Vegetation Heterogeneity", ylab = "Detection",
     xlim = c(0, 0.15), ylim = c(0, 1), bty = "l", las = 1)

#prefire
polygon(
  x = c(ndviHet_gf$prefire$ndvi_het, rev(ndviHet_gf$prefire$ndvi_het)),
  y = c(ndviHet_gf$prefire$lower, rev(ndviHet_gf$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = ndviHet_gf$prefire$ndvi_het, y = ndviHet_gf$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)


#burned
polygon(
  x = c(ndviHet_gf$burned_post$ndvi_het, rev(ndviHet_gf$burned_post$ndvi_het)),
  y = c(ndviHet_gf$burned_post$lower, rev(ndviHet_gf$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = ndviHet_gf$burned_post$ndvi_het, y = ndviHet_gf$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(ndviHet_gf$unburned_post$ndvi_het, rev(ndviHet_gf$unburned_post$ndvi_het)),
  y = c(ndviHet_gf$unburned_post$lower, rev(ndviHet_gf$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = ndviHet_gf$unburned_post$ndvi_het, y = ndviHet_gf$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("bottomleft", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n")

#### SUMMARY TABLES ####

#### spatial ####
# setup
param_h <- c(
  paste(greeks("psi"), "- Intercept"), paste(greeks("psi"), "- Burn Status (Burned)"), paste(greeks("psi"), "- Fire Heterogeneity"), 
  paste(greeks("psi"), "- Burned x Heterogeneity"), paste(greeks("psi"), "-", greeks("theta")), paste(greeks("rho"), "- Intercept"), paste(greeks("rho"), "- Microsite Attractant"),
  paste(greeks("rho"), "- Burn Status (Burned)"),  paste(greeks("rho"), "- Fire Heterogeneity"),  paste(greeks("rho"), "- Burned x Heterogeneity")
)
param_r <- c(
  paste(greeks("psi"), "- Intercept"), paste(greeks("psi"), "- Burn Status (Burned)"), paste(greeks("psi"), "- Fire Severity (RdNBR)"), paste(greeks("psi"), "- RdNBR^2"),
  paste(greeks("psi"), "-", greeks("theta")), paste(greeks("rho"), "- Intercept"), paste(greeks("rho"), "- Microsite Attractant"), paste(greeks("rho"), "- Burn Status (Burned)"),
  paste(greeks("rho"), "- Fire Severity (RdNBR)"), paste(greeks("rho"), "- RdNBR^2")
)
param_d <- c(
  paste(greeks("psi"), "- Intercept"), paste(greeks("psi"), "- Burn Status (Burned)"), paste(greeks("psi"), "- Distance to Fire Perimeter"), paste(greeks("psi"), "- Burned x Distance to Perimeter"),
  paste(greeks("psi"), "-", greeks("theta")), paste(greeks("rho"), "- Intercept"), paste(greeks("rho"), "- Microsite Attractant"), 
  paste(greeks("rho"), "- Burn Status (Burned)"), paste(greeks("rho"), "- Distance to Fire Perimeter"), paste(greeks("rho"), "- Burned x Distance to Perimeter")
)

estimate_h <- round(gfS_het@estimates$Est, 3)
estimate_r <- round(gfS_rdnbr@estimates$Est, 3)
estimate_d <- round(gfS_dp@estimates$Est, 3)

SE_h <- round(gfS_het@estimates$SE, 3)
SE_r <- round(gfS_rdnbr@estimates$SE, 3)
SE_d <- round(gfS_dp@estimates$SE, 3)

p_h <- round(gfS_het@estimates$p, 3)
p_r <- round(gfS_rdnbr@estimates$p, 3)
p_d <- round(gfS_dp@estimates$p, 3)

# fire heterogeneity 
gfHet_tab <- bind_cols(param_h, estimate_h, SE_h, p_h)
gfHet_tab <- gfHet_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
gfHet_tab <- gfHet_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

gf_h <- gfHet_tab %>% 
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
  tab_spanner(label = md("Fire Heterogeneity Model (AIC = 1409.6)"),
              columns = everything())

gf_h %>% 
  gtsave(filename = "gf_h.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# fire severity
gfR_tab <- bind_cols(param_r, estimate_r, SE_r, p_r)
gfR_tab <- gfR_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
gfR_tab <- gfR_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

gf_r <- gfR_tab %>% 
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
  tab_spanner(label = md("Fire Severity Model (AIC = 1403.2)"),
              columns = everything())

gf_r %>%  
  gtsave(filename = "gf_r.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# distance to fire perimeter
gfD_tab <- bind_cols(param_d, estimate_d, SE_d, p_d)
gfD_tab <- gfD_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
gfD_tab <- gfD_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

gf_d <- gfD_tab %>% 
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
  tab_spanner(label = md("Distance to Fire Perimeter Model (AIC = 1411.5)"),
              columns = everything()) 

gf_d %>% 
  gtsave(filename = "gf_d.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

#### temporal ####
parameter <- coyT_global@estimates$parameter
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

estimate <- gfT_global@estimates$Est
estimate <- round(estimate, 3)

SE <- gfT_global@estimates$SE
SE <- round(SE, 3)

p <- gfT_global@estimates$p
p <- round(p, 3)

gfGlobal_tab <- bind_cols(parameter, estimate, SE, p)
gfGlobal_tab <- gfGlobal_tab %>% 
  rename("Parameter" = ...1, 
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
gfGlobal_tab <- gfGlobal_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

gf_g <- gfGlobal_tab %>% 
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

gf_g %>% 
  gtsave(filename = "gf_g.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")



