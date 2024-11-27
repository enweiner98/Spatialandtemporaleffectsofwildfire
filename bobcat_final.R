#### BOBCAT FINAL ANALYSIS #####
# single-species, multi-factor models
# approach: temporal & spatial global models

#### LIBRARIES ####
library(tidyverse)
library(caret)
library(stats)
library(autoOcc)
library(beepr)
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

#### BOBCAT DATAFRAMES ####
bob_all <- all_occu$bobcat  # all bobcat detections
bob_post <- post_occu$bobcat  # post-fire bobcat detections

## summary statistics
# temporal: total detections
bob_all %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 1259 total detections

# detections at burned v. unburned sites 
bob_bs <- bob_all
bob_bs$burn_status <- var_cov$burned

bs_count <- bob_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 323 detections 
# unburned: 551 detections 
# prefire: 385 detections

# temporal: total sites 
tot_sites <- bob_all %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 67 sites 

# spatial: total detections
bob_post %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 369 total detections

# detections at burned v. unburned sites 
bob_bs <- bob_post
bob_bs$burn_status <- post_varCov$burned

bs_count <- bob_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 218 detections 
# unburned: 151 detections 

# spatial: total sites 
tot_sites <- bob_post %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 46 sites 

## creating dataframes for autoOcc
# post-fire only
bobcat_S <- format_y(
  x = bob_post,
  site_column = "Site",
  time_column = "Season",
  history_column = "Week"
)

# includes pre-fire and post-fire data
bobcat_T <- format_y(        
  x = bob_all,
  site_column = "Site",
  time_column = "Season",
  history_columns = "Week" 
)

#### SPATIAL ANALYSIS ####
#### null model ####
bobS_null <- auto_occ(
  formula = ~1~1,
  y = bobcat_S
)

# overall expected occupancy: 0.533
(
  intercept_preds_psi <- predict(
    bobS_null,
    type = "psi"
  )
)

# overall expected detection: 0.409
(
  intercept_preds_rho <- predict(
    bobS_null, 
    type = "rho"
  )
)

#### p1: fire heterogeneity ####
bobS_het <- auto_occ(
  formula = ~onroad + burn_status * rdnbr_het
  ~burn_status * rdnbr_het,
  y = bobcat_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(bobS_het)

#### p2: fire severity ####
bobS_rdnbr <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad
  ~burn_status + rdnbr + rdnbrQuad,
  y = bobcat_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(bobS_rdnbr)

#### p3: distance to fire perimeter ####
bobS_dp <- auto_occ(
  formula = ~onroad + burn_status * dist_perim
  ~burn_status * dist_perim,
  y = bobcat_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(bobS_dp)

#### 2 factor: severity ####
bobS_sev_dp <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim
  ~burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim,
  y = bobcat_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(bobS_sev_dp)

#### 2 factor: heterogeneity ####
bobS_het_dp <- auto_occ(
  formula = ~onroad + burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim
  ~burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim,
  y = bobcat_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(bobS_het_dp)

#### SPATIAL MODEL SELECTION ####
compare_models(list(bobS_null, bobS_het, bobS_rdnbr, bobS_dp, bobS_sev_dp, bobS_het_dp),
               digits = 2) # best fit model: m2 - rdnbr_het

#### SPATIAL MODEL PREDICTIONS ####

#### fire heterogeneity ####
bob_rdnbrHet_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  
  rdnbr_het = rdnbrHet,
  dist_perim = 0,
  onroad = 0
)

# scaling
bob_rdnbrHet_scaled <- bob_rdnbrHet_real
bob_rdnbrHet_scaled$rdnbr_het <- (
  bob_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_bob <- predict(
  object = bobS_het_dp,
  type = "rho",
  newdata = bob_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_bob <- data.frame(
  rdnbrHet_bob,
  bob_rdnbrHet_real
)

# plot 
ggplot(rdnbrHet_bob, aes(rdnbr_het, estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Fire Heterogeneity", y = "Detection") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*fire heterogeneity ####
bob_rdnbrHet_real <- data.frame(
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
bob_rdnbrHet_scaled <- bob_rdnbrHet_real
bob_rdnbrHet_scaled$rdnbr_het <- (
  bob_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_bob <- predict(
  object = bobS_het_dp,
  type = "psi",
  newdata = bob_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_bob <- data.frame(
  rdnbrHet_bob,
  bob_rdnbrHet_real
)

rdnbrHet_bob <- split(rdnbrHet_bob, rdnbrHet_bob$burn_status)

#plot
plot(1~1, type = "n", xlab = "Fire Heterogeneity", ylab = "Occupancy",
     xlim = c(20, 1750), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(rdnbrHet_bob$burned$rdnbr_het, rev(rdnbrHet_bob$burned$rdnbr_het)),
  y = c(rdnbrHet_bob$burned$lower, rev(rdnbrHet_bob$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = rdnbrHet_bob$burned$rdnbr_het, y = rdnbrHet_bob$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(rdnbrHet_bob$unburned$rdnbr_het, rev(rdnbrHet_bob$unburned$rdnbr_het)),
  y = c(rdnbrHet_bob$unburned$lower, rev(rdnbrHet_bob$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = rdnbrHet_bob$unburned$rdnbr_het, y = rdnbrHet_bob$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topleft", c("Unburned", "Burned"), lwd = 2, lty = c(2, 1), col = c("darkgreen", "darkorange"), bty = "n", cex=1.25)

b1 <- recordPlot()
b2 <- recordPlot()

b1 # occupancy
b2 # detection

#### fire severity ####
bob_rdnbr_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("unburned", "burned")
  ),
  rdnbr = rdnbr_seq,
  rdnbrQuad = 0,
  onroad = 0
)

# scaling
bob_rdnbr_scaled <- bob_rdnbr_real
bob_rdnbr_scaled$rdnbr <- (
  bob_rdnbr_scaled$rdnbr - mean(post_varCov$rdnbr_abs)
) / sd(post_varCov$rdnbr_abs)
view(bob_rdnbr_scaled)

# the model prediction
rdnbr_bob <- predict(
  object = bobS_rdnbr,
  type = "rho",
  newdata = bob_rdnbr_scaled
)

# add on covariate data
rdnbr_bob <- data.frame(
  rdnbr_bob,
  bob_rdnbr_real
)

b4 <- ggplot(rdnbr_bob, aes(x = rdnbr, y = estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Fire Severity (RdNBR)", y = "Detection") + 
  ylim(0,1) +
  theme_classic(18)

b3 # occupancy
b4 # detection

#### distance to fire perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

bob_dp_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  dist_perim = dp_seq,
  rdnbr_het = 0,
  onroad = 0
)

# scaling
bob_dp_scaled <- bob_dp_real
bob_dp_scaled$dist_perim <- (
  bob_dp_scaled$dist_perim - mean(post_varCov$dist_perim)
) / sd(post_varCov$dist_perim)

# the model prediction
dp_bob <- predict(
  object = bobS_het_dp,
  type = "psi",
  newdata = bob_dp_real
)

# add on covariate data
dp_bob <- data.frame(
  dp_bob,
  bob_dp_real
)

ggplot(dp_bob, aes(x = dist_perim, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Distance to Fire Perimeter (km)", y = "Occupancy", title = "Bobcat") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*distance to perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

bob_dp_real <- data.frame(
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
dp_bob <- predict(
  object = bobS_het_dp,
  type = "psi",
  newdata = bob_dp_real
)

dp_bob$dist_perim <- c(-6:0, 0:9)
dp_bob$burn_status <- bob_dp_real$burn_status

dp_bob <- split(dp_bob, bob_dp_real$burn_status)

#plot
plot(1~1, type = "n", xlab = "Distance to Fire Perimeter (km)", ylab = "Occupancy",
     xlim = c(-5, 5), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(dp_bob$burned$dist_perim, rev(dp_bob$burned$dist_perim)),
  y = c(dp_bob$burned$lower, rev(dp_bob$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = dp_bob$burned$dist_perim, y = dp_bob$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(dp_bob$unburned$dist_perim, rev(dp_bob$unburned$dist_perim)),
  y = c(dp_bob$unburned$lower, rev(dp_bob$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = dp_bob$unburned$dist_perim, y = dp_bob$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("bottomright", c("Burned", "Unburned", "Fire Perimeter"), lwd = 2, lty = c(1, 2, 1), col = c("darkorange", "darkgreen", "black"), bty = "n", cex=1.25)

b5 <- recordPlot()
b6 <- recordPlot()

b5 # occupancy
b6 # detection

#### TEMPORAL ANALYSIS ####

#### null model ####
bobT_null <- auto_occ(
  formula = ~1~1,
  y = bobcat_T
)

# overall expected occupancy: 0.535
(
  intercept_preds_psi <- predict(
    bobT_null,
    type = "psi"
  )
)

# overall expected detection: 0.397
(
  intercept_preds_rho <- predict(
    bobT_null, 
    type = "rho"
  )
)

#### global model ####
bobT_global <- auto_occ(
  formula = ~onroad + burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi
  ~burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi,
  y = bobcat_T,
  det_covs = covFrame_T,
  occ_covs = covFrame_T
)
summary(bobT_global)

#### TEMPORAL MODEL SELECTION ####
compare_models(list(bobT_null, bobT_global), digits = 2)

#### TEMPORAL MODEL PREDICTIONS ####

#### burn status ####
burned_seq <- c("prefire", "burned_post", "unburned_post")
season_levels <- levels(cov_frame$season$V1)

bob_burn_real <- data.frame(
  burn_status = factor(burned_seq, 
                       levels = burned_seq),
  tsf = -3,
  ndvi = 0, 
  ndvi_het = 0,
  onroad = 0
)

# the model prediction
burn_bob <- predict(
  object = bobT_global,
  type = "rho",
  newdata = bob_burn_real
)

# add on the covariate data
burn_bob <- data.frame(
  burn_bob,
  bob_burn_real
)

# plot
b8 <- ggplot(burn_bob, aes(x = factor(burn_status,
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

b7 # occupancy
b8 # detection

#### burn_status*tsf ####
bob_int <- data.frame(
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
bob_int$burn_status <- factor(bob_int$burn_status,
                             levels = c("prefire", "burned_post", "unburned_post"))

bob_int_pred <- predict(bobT_global,
                       type = 'rho',
                       newdata = bob_int,
                       interval = "confidence")

bob_int_pred$tsf <- c(-22:0, 0:22, 0:22)
bob_int_pred$burn_status <- bob_int$burn_status

bob_int_pred <- split(bob_int_pred, bob_int_pred$burn_status)

#plot
plot(1~1, type = "n", xlab = "Seasons since fire", ylab = "Detection",
     xlim = c(-22, 23), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(bob_int_pred$prefire$tsf, rev(bob_int_pred$prefire$tsf)),
  y = c(bob_int_pred$prefire$lower, rev(bob_int_pred$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = bob_int_pred$prefire$tsf, y = bob_int_pred$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(bob_int_pred$burned_post$tsf, rev(bob_int_pred$burned_post$tsf)),
  y = c(bob_int_pred$burned_post$lower, rev(bob_int_pred$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = bob_int_pred$burned_post$tsf, y = bob_int_pred$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(bob_int_pred$unburned_post$tsf, rev(bob_int_pred$unburned_post$tsf)),
  y = c(bob_int_pred$unburned_post$lower, rev(bob_int_pred$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = bob_int_pred$unburned_post$tsf, y = bob_int_pred$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topright", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

b9 <- recordPlot()
b10 <- recordPlot()

b9 # occupancy
b10 # detection

#### ndvi ####
ndvi_seq <- seq(0, 0.30, 0.01)

bob_ndvi_real <- data.frame(
  burn_status = factor("burned_post",
                       levels = c("prefire", "burned_post", "unburned_post")
  ),
  tsf = -3,
  ndvi = ndvi_seq,
  ndvi_het = 0,
  onroad = 0
)

# scaling
bob_ndvi_scaled <- bob_ndvi_real
bob_ndvi_scaled$ndvi <- (
  bob_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_bob <- predict(
  object = bobT_global,
  type = "psi",
  newdata = bob_ndvi_scaled
)

# add on covariate data
ndvi_bob <- data.frame(
  ndvi_bob,
  bob_ndvi_real
)

# plot
ggplot(ndvi_bob, aes(x = ndvi, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "NDVI", y = "Detection") + 
  ylim(0,1) +
  theme_classic(18)

#### burn_status*ndvi ####
ndvi_seq <- seq(0, 0.35, 0.01)

bob_ndvi_real <- data.frame(
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
bob_ndvi_scaled <- bob_ndvi_real
bob_ndvi_scaled$ndvi <- (
  bob_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_bob <- predict(
  object = bobT_global,
  type = "rho",
  newdata = bob_ndvi_scaled
)

# add on covariate data
ndvi_bob <- data.frame(
  ndvi_bob,
  bob_ndvi_real
)

ndvi_bob <- split(ndvi_bob, ndvi_bob$burn_status)

#plot
plot(1~1, type = "n", xlab = "NDVI", ylab = "Detection",
     xlim = c(0, 0.3), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(ndvi_bob$prefire$ndvi, rev(ndvi_bob$prefire$ndvi)),
  y = c(ndvi_bob$prefire$lower, rev(ndvi_bob$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = ndvi_bob$prefire$ndvi, y = ndvi_bob$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(ndvi_bob$burned_post$ndvi, rev(ndvi_bob$burned_post$ndvi)),
  y = c(ndvi_bob$burned_post$lower, rev(ndvi_bob$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = ndvi_bob$burned_post$ndvi, y = ndvi_bob$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(ndvi_bob$unburned_post$ndvi, rev(ndvi_bob$unburned_post$ndvi)),
  y = c(ndvi_bob$unburned_post$lower, rev(ndvi_bob$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = ndvi_bob$unburned_post$ndvi, y = ndvi_bob$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topright", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

b11 <- recordPlot()
b12 <- recordPlot()

b11 # occupancy
b12 # detection

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

estimate_h <- round(bobS_het@estimates$Est, 3)
estimate_r <- round(bobS_rdnbr@estimates$Est, 3)
estimate_d <- round(bobS_dp@estimates$Est, 3)

SE_h <- round(bobS_het@estimates$SE, 3)
SE_r <- round(bobS_rdnbr@estimates$SE, 3)
SE_d <- round(bobS_dp@estimates$SE, 3)

p_h <- round(bobS_het@estimates$p, 3)
p_r <- round(bobS_rdnbr@estimates$p, 3)
p_d <- round(bobS_dp@estimates$p, 3)

# fire heterogeneity 
bobHet_tab <- bind_cols(param_h, estimate_h, SE_h, p_h)
bobHet_tab <- bobHet_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
bobHet_tab <- bobHet_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

bh <- bobHet_tab %>% 
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
  tab_spanner(label = md("Fire Heterogeneity Model (AIC = 1552.1)"),
              columns = everything())

bh %>% 
  gtsave(filename = "bh.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# fire severity
bobR_tab <- bind_cols(param_r, estimate_r, SE_r, p_r)
bobR_tab <- bobR_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
bobR_tab <- bobR_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

br <- bobR_tab %>% 
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
  tab_spanner(label = md("Fire Severity Model (AIC = 1560.4)"),
              columns = everything())

br %>%  
  gtsave(filename = "br.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# distance to fire perimeter
bobD_tab <- bind_cols(param_d, estimate_d, SE_d, p_d)
bobD_tab <- bobD_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
bobD_tab <- bobD_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

bd <- bobD_tab %>% 
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
  tab_spanner(label = md("Distance to Fire Perimeter Model (AIC = 1554.8)"),
              columns = everything())

bd %>% 
  gtsave(filename = "bd.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

#### temporal ####
parameter <- bobT_global@estimates$parameter
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

estimate <- bobT_global@estimates$Est
estimate <- round(estimate, 3)

SE <- bobT_global@estimates$SE
SE <- round(SE, 3)

p <- bobT_global@estimates$p
p <- round(p, 3)

bobGlobal_tab <- bind_cols(parameter, estimate, SE, p)
bobGlobal_tab <- bobGlobal_tab %>% 
  rename("Parameter" = ...1, 
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
bobGlobal_tab <- bobGlobal_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

bg <- bobGlobal_tab %>% 
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

bg %>% 
  gtsave(filename = "bg.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")
