#### STRIPED SKUNK FINAL ANALYSIS #####
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

#### STRIPED SKUNK DATAFRAMES ####
ss_all <- all_occu$striped_skunk  # all mule deer detections
ss_post <- post_occu$striped_skunk  # post-fire mule deer detections

## summary statistics
# temporal: total detections
ss_all %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 545 total detections

# temporal: detections at burned v. unburned sites 
ss_bs <- ss_all
ss_bs$burn_status <- var_cov$burned

bs_count <- ss_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 77 detections 
# unburned: 348 detections 
# prefire: 120 detections

# temporal: total sites 
tot_sites <- ss_all %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 60 sites 

# spatial: total detections
ss_post %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 205 total detections

# spatial: detections at burned v. unburned sites 
ss_bs <- ss_post
ss_bs$burn_status <- post_varCov$burned

bs_count <- ss_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 106 detections 
# unburned: 99 detections 

# spatial: total sites 
tot_sites <- ss_post %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 41 sites 

## creating dataframes for autoOcc
# post-fire only
stripedSkunk_S <- format_y(
  x = ss_post,
  site_column = "Site",
  time_column = "Season",
  history_column = "Week"
)

# includes pre-fire and post-fire data
stripedSkunk_T <- format_y(        
  x = ss_all,
  site_column = "Site",
  time_column = "Season",
  history_columns = "Week" 
)

#### SPATIAL ANALYSIS ####
#### null model ####
ssS_null <- auto_occ(
  formula = ~1~1,
  y = stripedSkunk_S
)

# overall expected occupancy: 0.322
(
  intercept_preds_psi <- predict(
    ssS_null,
    type = "psi"
  )
)

# overall expected detection: 0.359
(
  intercept_preds_rho <- predict(
    ssS_null, 
    type = "rho"
  )
)

#### p1: fire heterogeneity ####
ssS_het <- auto_occ(
  formula = ~onroad + burn_status * rdnbr_het
  ~burn_status * rdnbr_het,
  y = stripedSkunk_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(ssS_het)

#### p2: fire severity ####
ssS_rdnbr <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad
  ~burn_status + rdnbr + rdnbrQuad,
  y = stripedSkunk_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(ssS_rdnbr)

#### p3: distance to fire perimeter ####
ssS_dp <- auto_occ(
  formula = ~onroad + burn_status * dist_perim
  ~burn_status * dist_perim,
  y = stripedSkunk_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(ssS_dp)

#### 2 factor: severity ####
ssS_sev_dp <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim
  ~burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim,
  y = stripedSkunk_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(ssS_sev_dp)

#### 2 factor: heterogeneity ####
ssS_het_dp <- auto_occ(
  formula = ~onroad + burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim
  ~burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim,
  y = stripedSkunk_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(ssS_het_dp)

#### SPATIAL MODEL SELECTION ####
autoOcc::compare_models(list(ssS_null, ssS_het, ssS_rdnbr, ssS_dp, ssS_sev_dp, ssS_het_dp), digits=2) # best fit model: rdnbr

#### SPATIAL MODEL PREDICTIONS ####

#### fire heterogeneity ####
ss_rdnbrHet_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  
  rdnbr_het = rdnbrHet,
  dist_perim = 0,
  onroad = 0
)

# scaling
ss_rdnbrHet_scaled <- ss_rdnbrHet_real
ss_rdnbrHet_scaled$rdnbr_het <- (
  ss_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_ss <- predict(
  object = ssS_het_dp,
  type = "rho",
  newdata = ss_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_ss <- data.frame(
  rdnbrHet_ss,
  ss_rdnbrHet_real
)

# plot 
ggplot(rdnbrHet_ss, aes(rdnbr_het, estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Fire Heterogeneity", y = "Detection", title = "Striped Skunk") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*fire heterogeneity ####
ss_rdnbrHet_real <- data.frame(
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
ss_rdnbrHet_scaled <- ss_rdnbrHet_real
ss_rdnbrHet_scaled$rdnbr_het <- (
  ss_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_ss <- predict(
  object = ssS_het_dp,
  type = "rho",
  newdata = ss_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_ss <- data.frame(
  rdnbrHet_ss,
  ss_rdnbrHet_real
)

rdnbrHet_ss <- split(rdnbrHet_ss, rdnbrHet_ss$burn_status)

#plot
plot(1~1, type = "n", xlab = "Fire Heterogeneity", ylab = "Detection",
     xlim = c(20, 1750), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(rdnbrHet_ss$burned$rdnbr_het, rev(rdnbrHet_ss$burned$rdnbr_het)),
  y = c(rdnbrHet_ss$burned$lower, rev(rdnbrHet_ss$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = rdnbrHet_ss$burned$rdnbr_het, y = rdnbrHet_ss$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(rdnbrHet_ss$unburned$rdnbr_het, rev(rdnbrHet_ss$unburned$rdnbr_het)),
  y = c(rdnbrHet_ss$unburned$lower, rev(rdnbrHet_ss$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = rdnbrHet_ss$unburned$rdnbr_het, y = rdnbrHet_ss$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topleft", c("Unburned", "Burned"), lwd = 2, lty = c(2, 1), col = c("darkgreen", "darkorange"), bty = "n", cex=1.25)

s1 <- recordPlot()
s2 <- recordPlot()

s1 # occupancy
s2 # detection

#### fire severity ####
ss_rdnbr_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("unburned", "burned")
  ),
  rdnbr = rdnbr_seq,
  rdnbrQuad = 0,
  onroad = 0
)

# scaling
ss_rdnbr_scaled <- ss_rdnbr_real
ss_rdnbr_scaled$rdnbr <- (
  ss_rdnbr_scaled$rdnbr - mean(post_varCov$rdnbr_abs)
) / sd(post_varCov$rdnbr_abs)

# the model prediction
rdnbr_ss <- predict(
  object = ssS_rdnbr,
  type = "rho",
  newdata = ss_rdnbr_scaled
)

# add on covariate data
rdnbr_ss <- data.frame(
  rdnbr_ss,
  ss_rdnbr_real
)

s4 <- ggplot(rdnbr_ss, aes(x = rdnbr, y = estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Fire Severity (RdNBR)", y = "Detection") +
  ylim(0,1) +
  theme_classic(18)

s3 # occupancy
s4 # detection

#### distance to fire perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

ss_dp_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  dist_perim = dp_seq,
  rdnbr_het = 0,
  onroad = 0
)

# scaling
ss_dp_scaled <- ss_dp_real
ss_dp_scaled$dist_perim <- (
  ss_dp_scaled$dist_perim - mean(post_varCov$dist_perim)
) / sd(post_varCov$dist_perim)

# the model prediction
dp_ss <- predict(
  object = ssS_het_dp,
  type = "psi",
  newdata = ss_dp_real
)

# add on covariate data
dp_ss <- data.frame(
  dp_ss,
  ss_dp_real
)

ggplot(dp_ss, aes(x = dist_perim, y = estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Distance to Fire Perimeter (km)", y = "Detection") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*distance to perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

ss_dp_real <- data.frame(
  burn_status = c(
    rep("burned", 7),
    rep("unburned", 10)
  ),
  dist_perim = c(
    -6:0,
    0:9
  ),
  rdnbr_het = 0,
  onroad = 0
)

# the model prediction
dp_ss <- predict(
  object = ssS_het_dp,
  type = "psi",
  newdata = ss_dp_real
)

dp_ss$dist_perim <- c(-6:0, 0:9)
dp_ss$burn_status <- ss_dp_real$burn_status

dp_ss <- split(dp_ss, ss_dp_real$burn_status)

#plot
plot(1~1, type = "n", xlab = "Distance to Fire Perimeter (km)", ylab = "Occupancy",
     xlim = c(-5, 5), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(dp_ss$burned$dist_perim, rev(dp_ss$burned$dist_perim)),
  y = c(dp_ss$burned$lower, rev(dp_ss$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = dp_ss$burned$dist_perim, y = dp_ss$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(dp_ss$unburned$dist_perim, rev(dp_ss$unburned$dist_perim)),
  y = c(dp_ss$unburned$lower, rev(dp_ss$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = dp_ss$unburned$dist_perim, y = dp_ss$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topleft", c("Burned", "Unburned", "Fire Perimeter"), lwd = 2, lty = c(1, 2, 1), col = c("darkorange", "darkgreen", "black"), bty = "n", cex=1.25)

s5 <- recordPlot()
s6 <- recordPlot()

s5 # occupancy
s6 # detection

#### TEMPORAL ANALYSIS ####

#### null model ####
ssT_null <- auto_occ(
  formula = ~1~1,
  y = stripedSkunk_T
)

# overall expected occupancy: 0.301
(
  intercept_preds_psi <- predict(
    ssT_null,
    type = "psi"
  )
)

# overall expected detection: 0.349
(
  intercept_preds_rho <- predict(
    ssT_null, 
    type = "rho"
  )
)

#### global model ####
ssT_global <- auto_occ(
  formula = ~onroad + burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi
  ~burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi,
  y = stripedSkunk_T,
  det_covs = covFrame_T,
  occ_covs = covFrame_T
)
summary(ssT_global)

#### TEMPORAL MODEL SELECTION ####
compare_models(list(ssT_null, ssT_global), digits = 2)

#### TEMPORAL MODEL PREDICTIONS ####

#### burn status ####
burned_seq <- c("prefire", "burned_post", "unburned_post")
season_levels <- levels(cov_frame$season$V1)

ss_burn_real <- data.frame(
  burn_status = factor(burned_seq, 
                       levels = burned_seq),
  tsf = -3,
  ndvi = 0, 
  onroad = 0
)

# the model prediction
burn_ss <- predict(
  object = ssT_global,
  type = "psi",
  newdata = ss_burn_real
)

# add on the covariate data
burn_ss <- data.frame(
  burn_ss,
  ss_burn_real
)

# plot
s7 <- ggplot(burn_ss, aes(x = factor(burn_status,
                               levels = c("prefire", "burned_post", "unburned_post")), 
                    y = estimate,
                    color = burn_status)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, show.legend = FALSE) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_vline(aes(xintercept = 1.5), color = "gray40", linetype = "longdash", size = 1) +
  labs(x = "Burn Status", y = "Occupancy") +
  scale_x_discrete(labels = c("Pre-Fire", "Burned", "Unburned")) +
  scale_color_manual(values = c("darkgreen", "darkorange", "darkgreen")) +
  ylim(0, 1) +
  theme_classic(18)

s7 # occupancy
s8 # detection

#### burn_status*tsf ####
ss_int <- data.frame(
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
ss_int$burn_status <- factor(ss_int$burn_status,
                             levels = c("prefire", "burned_post", "unburned_post"))

ss_int_pred <- predict(ssT_global,
                       type = 'psi',
                       newdata = ss_int,
                       interval = "confidence")

ss_int_pred$tsf <- c(-22:0, 0:22, 0:22)
ss_int_pred$burn_status <- ss_int$burn_status

ss_int_pred <- split(ss_int_pred, ss_int_pred$burn_status)

#plot
plot(1~1, type = "n", xlab = "Seasons since fire", ylab = "Detection",
     xlim = c(-22, 23), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(ss_int_pred$prefire$tsf, rev(ss_int_pred$prefire$tsf)),
  y = c(ss_int_pred$prefire$lower, rev(ss_int_pred$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = ss_int_pred$prefire$tsf, y = ss_int_pred$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(ss_int_pred$burned_post$tsf, rev(ss_int_pred$burned_post$tsf)),
  y = c(ss_int_pred$burned_post$lower, rev(ss_int_pred$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = ss_int_pred$burned_post$tsf, y = ss_int_pred$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(ss_int_pred$unburned_post$tsf, rev(ss_int_pred$unburned_post$tsf)),
  y = c(ss_int_pred$unburned_post$lower, rev(ss_int_pred$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = ss_int_pred$unburned_post$tsf, y = ss_int_pred$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topleft", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

s9 <- recordPlot() # occupancy
s10 <- recordPlot() # detection

#### ndvi ####
ndvi_seq <- seq(0, 0.30, 0.01)

ss_ndvi_real <- data.frame(
  burn_status = factor("burned_post",
                       levels = c("prefire", "burned_post", "unburned_post")
  ),
  tsf = -3,
  ndvi = ndvi_seq,
  onroad = 0
)

# scaling
ss_ndvi_scaled <- ss_ndvi_real
ss_ndvi_scaled$ndvi <- (
  ss_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_ss <- predict(
  object = ssT_global,
  type = "psi",
  newdata = ss_ndvi_scaled
)

# add on covariate data
ndvi_ss <- data.frame(
  ndvi_ss,
  ss_ndvi_real
)

# plot
ggplot(ndvi_ss, aes(x = ndvi, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Vegetation Greenness", y = "Occupancy", title = "Striped Skunk") + 
  ylim(0,1) +
  theme_classic(18)

#### burn_status*ndvi ####
ndvi_seq <- seq(0, 0.35, 0.01)

ss_ndvi_real <- data.frame(
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
ss_ndvi_scaled <- ss_ndvi_real
ss_ndvi_scaled$ndvi <- (
  ss_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_ss <- predict(
  object = ssT_global,
  type = "psi",
  newdata = ss_ndvi_scaled
)

# add on covariate data
ndvi_ss <- data.frame(
  ndvi_ss,
  ss_ndvi_real
)

ndvi_ss <- split(ndvi_ss, ndvi_ss$burn_status)

#plot
plot(1~1, type = "n", xlab = "Vegetation Biomass", ylab = "Occupancy",
     xlim = c(0, 0.3), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(ndvi_ss$prefire$ndvi, rev(ndvi_ss$prefire$ndvi)),
  y = c(ndvi_ss$prefire$lower, rev(ndvi_ss$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = ndvi_ss$prefire$ndvi, y = ndvi_ss$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(ndvi_ss$burned_post$ndvi, rev(ndvi_ss$burned_post$ndvi)),
  y = c(ndvi_ss$burned_post$lower, rev(ndvi_ss$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = ndvi_ss$burned_post$ndvi, y = ndvi_ss$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(ndvi_ss$unburned_post$ndvi, rev(ndvi_ss$unburned_post$ndvi)),
  y = c(ndvi_ss$unburned_post$lower, rev(ndvi_ss$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = ndvi_ss$unburned_post$ndvi, y = ndvi_ss$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topleft", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

s11 <- recordPlot() # occupancy
s12 <- recordPlot() # detection

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

estimate_h <- round(ssS_het@estimates$Est, 3)
estimate_r <- round(ssS_rdnbr@estimates$Est, 3)
estimate_d <- round(ssS_dp@estimates$Est, 3)

SE_h <- round(ssS_het@estimates$SE, 3)
SE_r <- round(ssS_rdnbr@estimates$SE, 3)
SE_d <- round(ssS_dp@estimates$SE, 3)

p_h <- round(ssS_het@estimates$p, 3)
p_r <- round(ssS_rdnbr@estimates$p, 3)
p_d <- round(ssS_dp@estimates$p, 3)

# fire heterogeneity 
ssHet_tab <- bind_cols(param_h, estimate_h, SE_h, p_h)
ssHet_tab <- ssHet_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
ssHet_tab <- ssHet_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

ss_h <- ssHet_tab %>% 
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
  tab_spanner(label = md("Fire Heterogeneity Model (AIC = 1047.8)"),
              columns = everything())

ss_h %>% 
  gtsave(filename = "ss_h.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# fire severity
ssR_tab <- bind_cols(param_r, estimate_r, SE_r, p_r)
ssR_tab <- ssR_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
ssR_tab <- ssR_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

ss_r <- ssR_tab %>% 
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
  tab_spanner(label = md("Fire Severity Model (AIC = 1046.2)"),
              columns = everything())

ss_r %>%  
  gtsave(filename = "ss_r.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# distance to fire perimeter
ssD_tab <- bind_cols(param_d, estimate_d, SE_d, p_d)
ssD_tab <- ssD_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
ssD_tab <- ssD_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

ss_d <- ssD_tab %>% 
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
  tab_spanner(label = md("Distance to Fire Perimeter Model (AIC = 1056.8)"),
              columns = everything())

ss_d %>% 
  gtsave(filename = "ss_d.html",
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

estimate <- ssT_global@estimates$Est
estimate <- round(estimate, 3)

SE <- ssT_global@estimates$SE
SE <- round(SE, 3)

p <- ssT_global@estimates$p
p <- round(p, 3)

ssGlobal_tab <- bind_cols(parameter, estimate, SE, p)
ssGlobal_tab <- ssGlobal_tab %>% 
  rename("Parameter" = ...1, 
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
ssGlobal_tab <- ssGlobal_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

ss_g <- ssGlobal_tab %>% 
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

ss_g %>% 
  gtsave(filename = "ss_g.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")




