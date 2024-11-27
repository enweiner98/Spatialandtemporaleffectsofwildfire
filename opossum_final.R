#### OPOSSUM FINAL ANALYSIS #####
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

#### OPOSSUM DATAFRAMES ####
o_all <- all_occu$virginia_opossum  # all mule deer detections
o_post <- post_occu$virginia_opossum  # post-fire mule deer detections

#### summary statistics ####
# temporal: total detections
o_all %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 190 total detections

# temporal: detections at burned v. unburned sites 
o_bs <- o_all
o_bs$burn_status <- var_cov$burned

bs_count <- o_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 36 detections 
# unburned: 95 detections 
# prefire: 59 detections

# temporal: total sites 
tot_sites <- o_all %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 30 sites 

# spatial: total detections
o_post %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 28 total detections

# spatial: detections at burned v. unburned sites 
o_bs <- o_post
o_bs$burn_status <- post_varCov$burned

bs_count <- o_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 2 detections 
# unburned: 26 detections 

# spatial: total sites 
tot_sites <- o_post %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 7 sites 

## creating dataframes for autoOcc
# post-fire only
opossum_S <- format_y(
  x = o_post,
  site_column = "Site",
  time_column = "Season",
  history_column = "Week"
)

# includes pre-fire and post-fire data
opossum_T <- format_y(        
  x = o_all,
  site_column = "Site",
  time_column = "Season",
  history_columns = "Week" 
)

#### SPATIAL ANALYSIS ####
#### null model ####
opoS_null <- auto_occ(
  formula = ~1~1,
  y = opossum_S
)

# overall expected occupancy: 0.128
(
  intercept_preds_psi <- predict(
    opoS_null,
    type = "psi"
  )
)

# overall expected detection: 0.219
(
  intercept_preds_rho <- predict(
    opoS_null, 
    type = "rho"
  )
)

#### p1: fire heterogeneity ####
opoS_het <- auto_occ(
  formula = ~onroad + burn_status + rdnbr_het
  ~burn_status + rdnbr_het,
  y = opossum_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(opoS_het)

#### p2: fire severity ####
opoS_rdnbr <- auto_occ(
  formula = ~onroad + burn_status + rdnbr
  ~burn_status + rdnbr,
  y = opossum_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(opoS_rdnbr)

#### p3: distance to fire perimeter ####
opoS_dp <- auto_occ(
  formula = ~onroad + burn_status + dist_perim
  ~burn_status + dist_perim,
  y = opossum_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(opoS_dp)

#### SPATIAL MODEL SELECTION ####
compare_models(list(opoS_null, opoS_het, opoS_rdnbr, opoS_dp),
               digits = 2) # best fit model: rdnbr

#### SPATIAL MODEL PREDICTIONS ####

#### fire heterogeneity ####
opo_rdnbrHet_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  
  rdnbr_het = rdnbrHet,
  onroad = 0
)

# scaling
opo_rdnbrHet_scaled <- opo_rdnbrHet_real
opo_rdnbrHet_scaled$rdnbr_het <- (
  opo_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_opo <- predict(
  object = opoS_het,
  type = "psi",
  newdata = opo_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_opo <- data.frame(
  rdnbrHet_opo,
  opo_rdnbrHet_real
)

# plot 
ggplot(rdnbrHet_opo, aes(rdnbr_het, estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Fire Heterogeneity", y = "Occupancy") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*fire heterogeneity ####
opo_rdnbrHet_real <- data.frame(
  burn_status = c(
    rep("burned", 174),
    rep("unburned", 174)
  ),
  rdnbr_het = c(
    rdnbrHet,
    rdnbrHet
  ),
  onroad = 0
)

# scaling
opo_rdnbrHet_scaled <- opo_rdnbrHet_real
opo_rdnbrHet_scaled$rdnbr_het <- (
  opo_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_opo <- predict(
  object = opoS_het,
  type = "rho",
  newdata = opo_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_opo <- data.frame(
  rdnbrHet_opo,
  opo_rdnbrHet_real
)

rdnbrHet_opo <- split(rdnbrHet_opo, rdnbrHet_opo$burn_status)

#plot
plot(1~1, type = "n", xlab = "Fire Heterogeneity", ylab = "Detection",
     xlim = c(20, 1750), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(rdnbrHet_opo$burned$rdnbr_het, rev(rdnbrHet_opo$burned$rdnbr_het)),
  y = c(rdnbrHet_opo$burned$lower, rev(rdnbrHet_opo$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = rdnbrHet_opo$burned$rdnbr_het, y = rdnbrHet_opo$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(rdnbrHet_opo$unburned$rdnbr_het, rev(rdnbrHet_opo$unburned$rdnbr_het)),
  y = c(rdnbrHet_opo$unburned$lower, rev(rdnbrHet_opo$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = rdnbrHet_opo$unburned$rdnbr_het, y = rdnbrHet_opo$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topleft", c("Unburned", "Burned"), lwd = 2, lty = c(2, 1), col = c("darkgreen", "darkorange"), bty = "n", cex=1.25)

o1 <- recordPlot() # occupancy
o2 <- recordPlot() # detection

#### fire severity ####
opo_rdnbr_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("unburned", "burned")
  ),
  rdnbr = rdnbr_seq,
  onroad = 0
)

# scaling
opo_rdnbr_scaled <- opo_rdnbr_real
opo_rdnbr_scaled$rdnbr <- (
  opo_rdnbr_scaled$rdnbr - mean(post_varCov$rdnbr_abs)
) / sd(post_varCov$rdnbr_abs)

# the model prediction
rdnbr_opo <- predict(
  object = opoS_rdnbr,
  type = "rho",
  newdata = opo_rdnbr_scaled
)

# add on covariate data
rdnbr_opo <- data.frame(
  rdnbr_opo,
  opo_rdnbr_real
)

o4 <- ggplot(rdnbr_opo, aes(x = rdnbr, y = estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Fire Severity (RdNBR)", y = "Detection") + 
  ylim(0,1) +
  theme_classic(18)

o3 # occupancy
o4 # detection

#### distance to fire perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

opo_dp_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  dist_perim = dp_seq,
  onroad = 0
)

# scaling
opo_dp_scaled <- opo_dp_real
opo_dp_scaled$dist_perim <- (
  opo_dp_scaled$dist_perim - mean(post_varCov$dist_perim)
) / sd(post_varCov$dist_perim)

# the model prediction
dp_opo <- predict(
  object = opoS_dp,
  type = "psi",
  newdata = opo_dp_real
)

# add on covariate data
dp_opo <- data.frame(
  dp_opo,
  opo_dp_real
)

ggplot(dp_opo, aes(x = dist_perim, y = estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Distance to Fire Perimeter (km)", y = "Occupancy") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*distance to perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

opo_dp_real <- data.frame(
  burn_status = c(
    rep("burned", 7),
    rep("unburned", 10)
  ),
  dist_perim = c(
    -6:0,
    0:9
  ),
  dist_perimQuad = 0,
  onroad = 0
)

# the model prediction
dp_opo <- predict(
  object = opoS_dp,
  type = "rho",
  newdata = opo_dp_real
)

dp_opo$dist_perim <- c(-6:0, 0:9)
dp_opo$burn_status <- opo_dp_real$burn_status

dp_opo <- split(dp_opo, opo_dp_real$burn_status)

#plot
plot(1~1, type = "n", xlab = "Distance to Fire Perimeter (km)", ylab = "Detection",
     xlim = c(-5, 5), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(dp_opo$burned$dist_perim, rev(dp_opo$burned$dist_perim)),
  y = c(dp_opo$burned$lower, rev(dp_opo$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = dp_opo$burned$dist_perim, y = dp_opo$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(dp_opo$unburned$dist_perim, rev(dp_opo$unburned$dist_perim)),
  y = c(dp_opo$unburned$lower, rev(dp_opo$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = dp_opo$unburned$dist_perim, y = dp_opo$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topleft", c("Burned", "Unburned", "Fire Perimeter"), lwd = 2, lty = c(1, 2, 1), col = c("darkorange", "darkgreen", "black"), bty = "n", cex = 1.25)

o5 <- recordPlot() # occupancy
o6 <- recordPlot() # detection

#### TEMPORAL ANALYSIS ####

#### null model ####
opoT_null <- auto_occ(
  formula = ~1~1,
  y = opossum_T
)

# overall expected occupancy: 0.289
(
  intercept_preds_psi <- predict(
    opoT_null,
    type = "psi"
  )
)

# overall expected detection: 0.178
(
  intercept_preds_rho <- predict(
    opoT_null, 
    type = "rho"
  )
)

#### global model ####
opoT_global <- auto_occ(
  formula = ~onroad + burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi
  ~burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi,
  y = opossum_T,
  det_covs = covFrame_T,
  occ_covs = covFrame_T
)
summary(opoT_global)

#### TEMPORAL MODEL SELECTION ####
compare_models(list(opoT_null, opoT_global), digits = 2)

#### TEMPORAL MODEL PREDICTIONS ####

#### burn status ####
burned_seq <- c("prefire", "burned_post", "unburned_post")
season_levels <- levels(cov_frame$season$V1)

opo_burn_real <- data.frame(
  burn_status = factor(burned_seq, 
                       levels = burned_seq),
  tsf = -3,
  ndvi = 0, 
  onroad = 0
)

# the model prediction
burn_opo <- predict(
  object = opoT_global,
  type = "rho",
  newdata = opo_burn_real
)

# add on the covariate data
burn_opo <- data.frame(
  burn_opo,
  opo_burn_real
)

# plot
o8 <- ggplot(burn_opo, aes(x = factor(burn_status,
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

o7 # occupancy
o8 # detection

#### burn_status*tsf ####
opo_int <- data.frame(
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
opo_int$burn_status <- factor(opo_int$burn_status,
                             levels = c("prefire", "burned_post", "unburned_post"))

opo_int_pred <- predict(opoT_global,
                       type = 'psi',
                       newdata = opo_int,
                       interval = "confidence")

opo_int_pred$tsf <- c(-22:0, 0:22, 0:22)
opo_int_pred$burn_status <- opo_int$burn_status

opo_int_pred <- split(opo_int_pred, opo_int_pred$burn_status)

#plot
plot(1~1, type = "n", xlab = "Seasons Since Fire", ylab = "Occupancy",
     xlim = c(-22, 23), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(opo_int_pred$prefire$tsf, rev(opo_int_pred$prefire$tsf)),
  y = c(opo_int_pred$prefire$lower, rev(opo_int_pred$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = opo_int_pred$prefire$tsf, y = opo_int_pred$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(opo_int_pred$burned_post$tsf, rev(opo_int_pred$burned_post$tsf)),
  y = c(opo_int_pred$burned_post$lower, rev(opo_int_pred$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = opo_int_pred$burned_post$tsf, y = opo_int_pred$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(opo_int_pred$unburned_post$tsf, rev(opo_int_pred$unburned_post$tsf)),
  y = c(opo_int_pred$unburned_post$lower, rev(opo_int_pred$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = opo_int_pred$unburned_post$tsf, y = opo_int_pred$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topleft", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

o9 <- recordPlot() # occupancy
o10 <- recordPlot() # detection

#### ndvi ####
ndvi_seq <- seq(0, 0.30, 0.01)

opo_ndvi_real <- data.frame(
  burn_status = factor("burned_post",
                       levels = c("prefire", "burned_post", "unburned_post")
  ),
  tsf = -3,
  ndvi = ndvi_seq,
  onroad = 0
)

# scaling
opo_ndvi_scaled <- opo_ndvi_real
opo_ndvi_scaled$ndvi <- (
  opo_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_opo <- predict(
  object = opoT_global,
  type = "rho",
  newdata = opo_ndvi_scaled
)

# add on covariate data
ndvi_opo <- data.frame(
  ndvi_opo,
  opo_ndvi_real
)

# plot
ggplot(ndvi_opo, aes(x = ndvi, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Vegetation Greenness", y = "Detection", title = "Opossum") + 
  ylim(0,1) +
  theme_classic(18)

#### burn_status*ndvi ####
ndvi_seq <- seq(0, 0.35, 0.01)

opo_ndvi_real <- data.frame(
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
opo_ndvi_scaled <- opo_ndvi_real
opo_ndvi_scaled$ndvi <- (
  opo_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_opo <- predict(
  object = opoT_global,
  type = "psi",
  newdata = opo_ndvi_scaled
)

# add on covariate data
ndvi_opo <- data.frame(
  ndvi_opo,
  opo_ndvi_real
)

ndvi_opo <- split(ndvi_opo, ndvi_opo$burn_status)

#plot
plot(1~1, type = "n", xlab = "NDVI", ylab = "Occupancy",
     xlim = c(0, 0.3), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(ndvi_opo$prefire$ndvi, rev(ndvi_opo$prefire$ndvi)),
  y = c(ndvi_opo$prefire$lower, rev(ndvi_opo$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = ndvi_opo$prefire$ndvi, y = ndvi_opo$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(ndvi_opo$burned_post$ndvi, rev(ndvi_opo$burned_post$ndvi)),
  y = c(ndvi_opo$burned_post$lower, rev(ndvi_opo$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = ndvi_opo$burned_post$ndvi, y = ndvi_opo$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(ndvi_opo$unburned_post$ndvi, rev(ndvi_opo$unburned_post$ndvi)),
  y = c(ndvi_opo$unburned_post$lower, rev(ndvi_opo$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = ndvi_opo$unburned_post$ndvi, y = ndvi_opo$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topright", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

o11 <- recordPlot() # occupancy
o12 <- recordPlot() # detection

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
  paste(greeks("psi"), "- Intercept"), paste(greeks("psi"), "- Burn Status (Burned)"), paste(greeks("psi"), "- Distance to Fire Perimeter"),
  paste(greeks("psi"), "-", greeks("theta")), paste(greeks("rho"), "- Intercept"), paste(greeks("rho"), "- Microsite Attractant"),
  paste(greeks("rho"), "- Burn Status (Burned)"), paste(greeks("rho"), "- Distance to Fire Perimeter")
)

estimate_h <- round(opoS_het@estimates$Est, 3)
estimate_r <- round(opoS_rdnbr@estimates$Est, 3)
estimate_d <- round(opoS_dp@estimates$Est, 3)

SE_h <- round(opoS_het@estimates$SE, 3)
SE_r <- round(opoS_rdnbr@estimates$SE, 3)
SE_d <- round(opoS_dp@estimates$SE, 3)

p_h <- round(opoS_het@estimates$p, 3)
p_r <- round(opoS_rdnbr@estimates$p, 3)
p_d <- round(opoS_dp@estimates$p, 3)

# fire heterogeneity 
opoHet_tab <- bind_cols(param_h, estimate_h, SE_h, p_h)
opoHet_tab <- opoHet_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
opoHet_tab <- opoHet_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

opo_h <- opoHet_tab %>% 
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
  tab_spanner(label = md("Fire Heterogeneity Model (AIC = 191.2)"),
              columns = everything())

opo_h %>% 
  gtsave(filename = "opo_h.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# fire severity
opoR_tab <- bind_cols(param_r, estimate_r, SE_r, p_r)
opoR_tab <- opoR_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
opoR_tab <- opoR_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

opo_r <- opoR_tab %>% 
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
  tab_spanner(label = md("Fire Severity Model (AIC = 190.2)"),
              columns = everything())

opo_r %>%  
  gtsave(filename = "opo_r.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# distance to fire perimeter
opoD_tab <- bind_cols(param_d, estimate_d, SE_d, p_d)
opoD_tab <- opoD_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
opoD_tab <- opoD_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

opo_d <- opoD_tab %>% 
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
  tab_spanner(label = md("Distance to Fire Perimeter Model (AIC = 196.6)"),
              columns = everything())

opo_d %>% 
  gtsave(filename = "opo_d.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

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

estimate <- opoT_global@estimates$Est
estimate <- round(estimate, 3)

SE <- opoT_global@estimates$SE
SE <- round(SE, 3)

p <- opoT_global@estimates$p
p <- round(p, 3)

opoGlobal_tab <- bind_cols(parameter, estimate, SE, p)
opoGlobal_tab <- opoGlobal_tab %>% 
  rename("Parameter" = ...1, 
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
opoGlobal_tab <- opoGlobal_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

opo_g <- opoGlobal_tab %>% 
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

opo_g %>% 
  gtsave(filename = "opo_g.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")




