#### MOUNTAIN LION FINAL ANALYSIS #####
# single-species, multi-factor models
# approach: temporal & spatial global models

#### LIBRARIES ####
library(tidyverse)
library(caret)
library(stats)
library(autoOcc)
library(beepr)
library(greekLetters)
library(gt)

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

#### MOUNTAIN LION DATAFRAMES ####
ml_all <- all_occu$mountain_lion  # all mt lion detections
ml_post <- post_occu$mountain_lion  # post-fire mt lion detections

## summary statistics
# temporal: total detections
ml_all %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 740 total detections

# temporal: detections at burned v. unburned sites 
ml_bs <- ml_all
ml_bs$burn_status <- var_cov$burned

bs_count <- ml_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 226 detections 
# unburned: 297 detections 
# prefire: 217 detections

# temporal: total sites 
tot_sites <- ml_all %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 63 sites 

# spatial: total detections
ml_post %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 195 total detections

# detections at burned v. unburned sites 
ml_bs <- ml_post
ml_bs$burn_status <- post_varCov$burned

bs_count <- ml_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 119 detections 
# unburned: 76 detections 

# spatial: total sites 
tot_sites <- ml_post %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 40 sites 

## creating dataframes for autoOcc
# post-fire only
mtLion_S <- format_y(
  x = ml_post,
  site_column = "Site",
  time_column = "Season",
  history_column = "Week"
)

# includes pre-fire and post-fire data
mtLion_T <- format_y(        
  x = ml_all,
  site_column = "Site",
  time_column = "Season",
  history_columns = "Week" 
)

#### SPATIAL ANALYSIS ####
#### null model ####
mlS_null <- auto_occ(
  formula = ~1~1,
  y = mtLion_S
)

# overall expected occupancy: 0.492
(
  intercept_preds_psi <- predict(
    mlS_null,
    type = "psi"
  )
)

# overall expected detection: 0.234
(
  intercept_preds_rho <- predict(
    mlS_null, 
    type = "rho"
  )
)

#### p1: fire heterogeneity ####
mlS_het <- auto_occ(
  formula = ~onroad + imperv + burn_status * rdnbr_het
  ~imperv + burn_status * rdnbr_het,
  y = mtLion_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(mlS_het)

#### p2: fire severity ####
mlS_rdnbr <- auto_occ(
  formula = ~onroad + imperv + burn_status + rdnbr + rdnbrQuad
  ~imperv + burn_status + rdnbr + rdnbrQuad,
  y = mtLion_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(mlS_rdnbr)

#### p3: distance to fire perimeter ####
mlS_dp <- auto_occ(
  formula = ~onroad + imperv + burn_status * dist_perim
  ~imperv + burn_status * dist_perim,
  y = mtLion_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(mlS_dp)

#### global model ####
mlS_global <- auto_occ(
  formula = ~onroad + burn_status + rdnbr_het + rdnbr + rdnbrQuad + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim
  ~burn_status + rdnbr_het + rdnbr + rdnbrQuad + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim,
  y = mtLion_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(mlS_global)

#### 2 factor: severity ####
mlS_sev_dp <- auto_occ(
  formula = ~onroad + imperv + burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim
  ~imperv + burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim,
  y = mtLion_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(mlS_sev_dp)

#### 2 factor: heterogeneity ####
mlS_het_dp <- auto_occ(
  formula = ~onroad + imperv + burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim
  ~imperv + burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim,
  y = mtLion_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(mlS_het_dp)

#### SPATIAL MODEL SELECTION ####
compare_models(list(mlS_null, mlS_het, mlS_rdnbr, mlS_dp, mlS_sev_dp, mlS_het_dp),
               digits = 2) # best fit model: fire severity/null

#### SPATIAL MODEL PREDICTIONS ####

#### fire heterogeneity ####
ml_rdnbrHet_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  
  rdnbr_het = rdnbrHet,
  onroad = 0,
  imperv = 0
)

# scaling
ml_rdnbrHet_scaled <- ml_rdnbrHet_real
ml_rdnbrHet_scaled$rdnbr_het <- (
  ml_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_ml <- predict(
  object = mlS_het,
  type = "rho",
  newdata = ml_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_ml <- data.frame(
  rdnbrHet_ml,
  ml_rdnbrHet_real
)

# plot 
ggplot(rdnbrHet_ml, aes(rdnbr_het, estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Fire Heterogeneity", y = "Detection", title = "Mountain Lion") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*fire heterogeneity ####
ml_rdnbrHet_real <- data.frame(
  burn_status = c(
    rep("burned", 174),
    rep("unburned", 174)
  ),
  rdnbr_het = c(
    rdnbrHet,
    rdnbrHet
  ),
  onroad = 0,
  imperv = 0
)

# scaling
ml_rdnbrHet_scaled <- ml_rdnbrHet_real
ml_rdnbrHet_scaled$rdnbr_het <- (
  ml_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_ml <- predict(
  object = mlS_het,
  type = "rho",
  newdata = ml_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_ml <- data.frame(
  rdnbrHet_ml,
  ml_rdnbrHet_real
)

rdnbrHet_ml <- split(rdnbrHet_ml, rdnbrHet_ml$burn_status)

#plot
plot(1~1, type = "n", xlab = "Fire Heterogeneity", ylab = "Detection",
     xlim = c(20, 1750), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(rdnbrHet_ml$burned$rdnbr_het, rev(rdnbrHet_ml$burned$rdnbr_het)),
  y = c(rdnbrHet_ml$burned$lower, rev(rdnbrHet_ml$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = rdnbrHet_ml$burned$rdnbr_het, y = rdnbrHet_ml$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(rdnbrHet_ml$unburned$rdnbr_het, rev(rdnbrHet_ml$unburned$rdnbr_het)),
  y = c(rdnbrHet_ml$unburned$lower, rev(rdnbrHet_ml$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = rdnbrHet_ml$unburned$rdnbr_het, y = rdnbrHet_ml$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topleft", c("Unburned", "Burned"), lwd = 2, lty = c(2, 1), col = c("darkgreen", "darkorange"), bty = "n", cex=1.25)

ml1 <- recordPlot()
ml2 <- recordPlot()

ml1 # occupancy
ml2 # detection

#### fire severity ####
ml_rdnbr_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("unburned", "burned")
  ),
  rdnbr = rdnbr_seq,
  rdnbrQuad = 0,
  onroad = 0,
  imperv = 0
)

# scaling
ml_rdnbr_scaled <- ml_rdnbr_real
ml_rdnbr_scaled$rdnbr <- (
  ml_rdnbr_scaled$rdnbr - mean(post_varCov$rdnbr_abs)
) / sd(post_varCov$rdnbr_abs)

# the model prediction
rdnbr_ml <- predict(
  object = mlS_rdnbr,
  type = "rho",
  newdata = ml_rdnbr_scaled
)

# add on covariate data
rdnbr_ml <- data.frame(
  rdnbr_ml,
  ml_rdnbr_real
)

ml4 <- ggplot(rdnbr_ml, aes(x = rdnbr, y = estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Fire Severity (RdNBR)", y = "Detection") + 
  ylim(0,1) +
  theme_classic(18)

ml3 # occupancy
ml4 # detection

#### distance to fire perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

ml_dp_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  dist_perim = dp_seq,
  onroad = 0
)

# scaling
ml_dp_scaled <- ml_dp_real
ml_dp_scaled$dist_perim <- (
  ml_dp_scaled$dist_perim - mean(post_varCov$dist_perim)
) / sd(post_varCov$dist_perim)

# the model prediction
dp_ml <- predict(
  object = mlS_dp,
  type = "psi",
  newdata = ml_dp_real
)

# add on covariate data
dp_ml <- data.frame(
  dp_ml,
  ml_dp_real
)

ggplot(dp_ml, aes(x = dist_perim, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Distance to Fire Perimeter (km)", y = "Occupancy", title = "Mountain Lion") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*distance to perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

ml_dp_real <- data.frame(
  burn_status = c(
    rep("burned", 7),
    rep("unburned", 10)
  ),
  dist_perim = c(
    -6:0,
    0:9
  ),
  onroad = 0,
  imperv = 0
)

# the model prediction
dp_ml <- predict(
  object = mlS_dp,
  type = "rho",
  newdata = ml_dp_real
)

dp_ml$dist_perim <- c(-6:0, 0:9)
dp_ml$burn_status <- ml_dp_real$burn_status

dp_ml <- split(dp_ml, ml_dp_real$burn_status)

#plot
plot(1~1, type = "n", xlab = "Distance to Fire Perimeter (km)", ylab = "Detection",
     xlim = c(-5, 5), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(dp_ml$burned$dist_perim, rev(dp_ml$burned$dist_perim)),
  y = c(dp_ml$burned$lower, rev(dp_ml$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = dp_ml$burned$dist_perim, y = dp_ml$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(dp_ml$unburned$dist_perim, rev(dp_ml$unburned$dist_perim)),
  y = c(dp_ml$unburned$lower, rev(dp_ml$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = dp_ml$unburned$dist_perim, y = dp_ml$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topright", c("Burned", "Unburned", "Fire Perimeter"), lwd = 2, lty = c(1, 2, 1), col = c("darkorange", "darkgreen", "black"), bty = "n", cex=1.25)

ml5 <- recordPlot()
ml6 <- recordPlot()

ml5 # occupancy
ml6 # detection

#### TEMPORAL ANALYSIS ####

#### null model ####
mlT_null <- auto_occ(
  formula = ~1~1,
  y = mtLion_T
)

# overall expected occupancy: 0.510
(
  intercept_preds_psi <- predict(
    mlT_null,
    type = "psi"
  )
)

# overall expected detection: 0.236
(
  intercept_preds_rho <- predict(
    mlT_null, 
    type = "rho"
  )
)

#### global model ####
mlT_global <- auto_occ(
  formula = ~onroad + imperv + burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi
  ~imperv + burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi,
  y = mtLion_T,
  det_covs = covFrame_T,
  occ_covs = covFrame_T
)
summary(mlT_global)

#### TEMPORAL MODEL SELECTION ####
compare_models(list(mlT_null, mlT_global),
               digits = 2) # best fit model: global

#### TEMPORAL MODEL PREDICTIONS ####

#### burn status ####
burned_seq <- c("prefire", "burned_post", "unburned_post")
season_levels <- levels(cov_frame$season$V1)

ml_burn_real <- data.frame(
  burn_status = factor(burned_seq, 
                       levels = burned_seq),
  tsf = -3,
  ndvi = 0, 
  ndvi_het = 0,
  onroad = 0,
  imperv = 0
)

# the model prediction
burn_ml <- predict(
  object = mlT_global,
  type = "rho",
  newdata = ml_burn_real
)

# add on the covariate data
burn_ml <- data.frame(
  burn_ml,
  ml_burn_real
)

# plot
ml8 <- ggplot(burn_ml, aes(x = factor(burn_status,
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

ml7 # occupancy
ml8 # detection

#### burn_status*tsf ####
ml_int <- data.frame(
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
  onroad = 0,
  imperv = 0
)
ml_int$burn_status <- factor(ml_int$burn_status,
                              levels = c("prefire", "burned_post", "unburned_post"))

ml_int_pred <- predict(mlT_global,
                        type = 'rho',
                        newdata = ml_int,
                        interval = "confidence")

ml_int_pred$tsf <- c(-22:0, 0:22, 0:22)
ml_int_pred$burn_status <- ml_int$burn_status

ml_int_pred <- split(ml_int_pred, ml_int_pred$burn_status)

#plot
plot(1~1, type = "n", xlab = "Seasons since fire", ylab = "Detection",
     xlim = c(-22, 23), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(ml_int_pred$prefire$tsf, rev(ml_int_pred$prefire$tsf)),
  y = c(ml_int_pred$prefire$lower, rev(ml_int_pred$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = ml_int_pred$prefire$tsf, y = ml_int_pred$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(ml_int_pred$burned_post$tsf, rev(ml_int_pred$burned_post$tsf)),
  y = c(ml_int_pred$burned_post$lower, rev(ml_int_pred$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = ml_int_pred$burned_post$tsf, y = ml_int_pred$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(ml_int_pred$unburned_post$tsf, rev(ml_int_pred$unburned_post$tsf)),
  y = c(ml_int_pred$unburned_post$lower, rev(ml_int_pred$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = ml_int_pred$unburned_post$tsf, y = ml_int_pred$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topright", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

ml9 <- recordPlot()
ml10 <- recordPlot()

ml9 # occupancy
ml10 # detection

#### ndvi ####
ndvi_seq <- seq(0, 0.30, 0.01)

ml_ndvi_real <- data.frame(
  burn_status = factor("burned_post",
                       levels = c("prefire", "burned_post", "unburned_post")
  ),
  tsf = -3,
  ndvi = ndvi_seq,
  ndvi_het = 0,
  onroad = 0,
  imperv = 0
)

# scaling
ml_ndvi_scaled <- ml_ndvi_real
ml_ndvi_scaled$ndvi <- (
  ml_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_ml <- predict(
  object = mlT_global,
  type = "rho",
  newdata = ml_ndvi_scaled
)

# add on covariate data
ndvi_ml <- data.frame(
  ndvi_ml,
  ml_ndvi_real
)

# plot
ggplot(ndvi_ml, aes(x = ndvi, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Vegetation Greenness", y = "Occupancy") + 
  ylim(0,1) +
  theme_classic(18)

#### burn_status*ndvi ####
ndvi_seq <- seq(0, 0.35, 0.01)

ml_ndvi_real <- data.frame(
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
  onroad = 0,
  imperv = 0
)

# scaling
ml_ndvi_scaled <- ml_ndvi_real
ml_ndvi_scaled$ndvi <- (
  ml_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_ml <- predict(
  object = mlT_global,
  type = "rho",
  newdata = ml_ndvi_scaled
)

# add on covariate data
ndvi_ml <- data.frame(
  ndvi_ml,
  ml_ndvi_real
)

ndvi_ml <- split(ndvi_ml, ndvi_ml$burn_status)

#plot
plot(1~1, type = "n", xlab = "Vegetation Greenness (NDVI)", ylab = "Detection",
     xlim = c(0, 0.3), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(ndvi_ml$prefire$ndvi, rev(ndvi_ml$prefire$ndvi)),
  y = c(ndvi_ml$prefire$lower, rev(ndvi_ml$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = ndvi_ml$prefire$ndvi, y = ndvi_ml$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(ndvi_ml$burned_post$ndvi, rev(ndvi_ml$burned_post$ndvi)),
  y = c(ndvi_ml$burned_post$lower, rev(ndvi_ml$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = ndvi_ml$burned_post$ndvi, y = ndvi_ml$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(ndvi_ml$unburned_post$ndvi, rev(ndvi_ml$unburned_post$ndvi)),
  y = c(ndvi_ml$unburned_post$lower, rev(ndvi_ml$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = ndvi_ml$unburned_post$ndvi, y = ndvi_ml$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topright", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

ml11 <- recordPlot()
ml12 <- recordPlot()

ml11 # occupancy
ml12 # detection

#### SUMMARY TABLES ####

#### spatial ####
# setup
param_h <- c(
  paste(greeks("psi"), "- Intercept"), paste(greeks("psi"), "- % Impervious Surface"), paste(greeks("psi"), "- Burn Status (Burned)"), paste(greeks("psi"), "- Fire Heterogeneity"), 
  paste(greeks("psi"), "- Burned x Heterogeneity"), paste(greeks("psi"), "-", greeks("theta")), paste(greeks("rho"), "- Intercept"), paste(greeks("rho"), "- Microsite Attractant"),
  paste(greeks("rho"), "- % Impervious Surface"), paste(greeks("rho"), "- Burn Status (Burned)"),  paste(greeks("rho"), "- Fire Heterogeneity"),  paste(greeks("rho"), "- Burned x Heterogeneity")
)
param_r <- c(
  paste(greeks("psi"), "- Intercept"), paste(greeks("psi"), "- % Impervious Surface"), paste(greeks("psi"), "- Burn Status (Burned)"), paste(greeks("psi"), "- Fire Severity (RdNBR)"), paste(greeks("psi"), "- RdNBR^2"),
  paste(greeks("psi"), "-", greeks("theta")), paste(greeks("rho"), "- Intercept"), paste(greeks("rho"), "- Microsite Attractant"), paste(greeks("rho"), "- % Impervious Surface"), paste(greeks("rho"), "- Burn Status (Burned)"),
  paste(greeks("rho"), "- Fire Severity (RdNBR)"), paste(greeks("rho"), "- RdNBR^2")
)
param_d <- c(
  paste(greeks("psi"), "- Intercept"), paste(greeks("psi"), "- % Impervious Surface"), paste(greeks("psi"), "- Burn Status (Burned)"), paste(greeks("psi"), "- Distance to Fire Perimeter"), paste(greeks("psi"), "- Burned x Distance to Perimeter"),
  paste(greeks("psi"), "-", greeks("theta")), paste(greeks("rho"), "- Intercept"), paste(greeks("rho"), "- Microsite Attractant"), paste(greeks("rho"), "- % Impervious Surface"),
  paste(greeks("rho"), "- Burn Status (Burned)"), paste(greeks("rho"), "- Distance to Fire Perimeter"), paste(greeks("rho"), "- Burned x Distance to Perimeter")
)

estimate_h <- round(mlS_het@estimates$Est, 3)
estimate_r <- round(mlS_rdnbr@estimates$Est, 3)
estimate_d <- round(mlS_dp@estimates$Est, 3)

SE_h <- round(mlS_het@estimates$SE, 3)
SE_r <- round(mlS_rdnbr@estimates$SE, 3)
SE_d <- round(mlS_dp@estimates$SE, 3)

p_h <- round(mlS_het@estimates$p, 3)
p_r <- round(mlS_rdnbr@estimates$p, 3)
p_d <- round(mlS_dp@estimates$p, 3)

# fire heterogeneity 
mlHet_tab <- bind_cols(param_h, estimate_h, SE_h, p_h)
mlHet_tab <- mlHet_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
mlHet_tab <- mlHet_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

ml_h <- mlHet_tab %>% 
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
  tab_spanner(label = md("Fire Heterogeneity Model (AIC = 1150.4)"),
              columns = everything())

ml_h %>% 
  gtsave(filename = "ml_h.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# fire severity
mlR_tab <- bind_cols(param_r, estimate_r, SE_r, p_r)
mlR_tab <- mlR_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
mlR_tab <- mlR_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

ml_r <- mlR_tab %>% 
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
  tab_spanner(label = md("Fire Severity Model (AIC = 1148.4)"),
              columns = everything())

ml_r %>%  
  gtsave(filename = "ml_r.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# distance to fire perimeter
mlD_tab <- bind_cols(param_d, estimate_d, SE_d, p_d)
mlD_tab <- mlD_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
mlD_tab <- mlD_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

ml_d <- mlD_tab %>% 
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
  tab_spanner(label = md("Distance to Fire Perimeter Model (AIC = 1153.2)"),
              columns = everything())

ml_d %>% 
  gtsave(filename = "ml_d.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

#### temporal ####
parameter <- bobT_global@estimates$parameter
parameter <- c(
  paste(greeks("psi"), "- Intercept"), paste(greeks("psi"), "- % Impervious Surface"), paste(greeks("psi"), "- Burn Status (Unburned)"), 
  paste(greeks("psi"), "- Burn Status (Pre-Fire)"), paste(greeks("psi"), "- Time Since Fire"),
  paste(greeks("psi"), "- NDVI"), paste(greeks("psi"), "- Unburned x Time Since Fire"), 
  paste(greeks("psi"), "- Pre-Fire x Time Since Fire"), paste(greeks("psi"), "- Unburned x NDVI"),
  paste(greeks("psi"), "- Pre-Fire x NDVI"), paste(greeks("psi"), "-", greeks("theta")), 
  paste(greeks("rho"), "- Intercept"), paste(greeks("rho"), "- % Impervious Surface"), paste(greeks("rho"), "- Microsite Attractant"), paste(greeks("rho"), "- Burn Status (Unburned)"), 
  paste(greeks("rho"), "- Burn Status (Pre-Fire)"), paste(greeks("rho"), "- Time Since Fire"),
  paste(greeks("rho"), "- NDVI"), paste(greeks("rho"), "- Unburned x Time Since Fire"), 
  paste(greeks("rho"), "- Pre-Fire x Time Since Fire"), paste(greeks("rho"), "- Unburned x NDVI"),
  paste(greeks("rho"), "- Pre-Fire x NDVI")
)

estimate <- mlT_global@estimates$Est
estimate <- round(estimate, 3)

SE <- mlT_global@estimates$SE
SE <- round(SE, 3)

p <- mlT_global@estimates$p
p <- round(p, 3)

mlGlobal_tab <- bind_cols(parameter, estimate, SE, p)
mlGlobal_tab <- mlGlobal_tab %>% 
  rename("Parameter" = ...1, 
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
mlGlobal_tab <- mlGlobal_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

ml_g <- mlGlobal_tab %>% 
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

ml_g %>% 
  gtsave(filename = "ml_g.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

