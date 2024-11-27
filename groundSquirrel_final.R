#### CA GROUND SQUIRREL FINAL ANALYSIS #####
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

#### CA GROUND SQUIRREL DATAFRAMES ####
cags_all <- all_occu$CA_ground_squirrel  # all CA ground squirrel detections
cags_post <- post_occu$CA_ground_squirrel  # post-fire CA ground squirrel detections

#### summary statistics ####
# temporal: total detections
cags_all %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 161 total detections

# detections at burned v. unburned sites 
cags_bs <- cags_all
cags_bs$burn_status <- var_cov$burned

bs_count <- cags_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 20 detections 
# unburned: 86 detections 
# prefire: 55 detections

# temporal: total sites 
tot_sites <- cags_all %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 26 sites 

# spatial: total detections
cags_post %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 61 total detections

# detections at burned v. unburned sites 
cags_bs <- cags_post
cags_bs$burn_status <- post_varCov$burned

bs_count <- cags_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 51 detections 
# unburned: 10 detections 

# spatial: total sites 
tot_sites <- cags_post %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 13 sites

## creating dataframes for autoOcc
# post-fire only
cags_S <- format_y(
  x = cags_post,
  site_column = "Site",
  time_column = "Season",
  history_column = "Week"
)

# includes pre-fire and post-fire data
cags_T <- format_y(        
  x = cags_all,
  site_column = "Site",
  time_column = "Season",
  history_columns = "Week" 
)

#### SPATIAL ANALYSIS ####
#### null model ####
cagsS_null <- auto_occ(
  formula = ~1~1,
  y = cags_S
)

# overall expected occupancy: 0.091
(
  intercept_preds_psi <- predict(
    cagsS_null,
    type = "psi"
  )
)

# overall expected detection: 0.345
(
  intercept_preds_rho <- predict(
    cagsS_null, 
    type = "rho"
  )
)

#### p1: fire heterogeneity ####
cagsS_het <- auto_occ(
  formula = ~onroad + burn_status * rdnbr_het
  ~burn_status * rdnbr_het,
  y = cags_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(cagsS_het)

#### p2: fire severity ####
cagsS_rdnbr <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad
  ~burn_status + rdnbr + rdnbrQuad,
  y = cags_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(cagsS_rdnbr)

#### p3: distance to fire perimeter ####
cagsS_dp <- auto_occ(
  formula = ~onroad + burn_status * dist_perim
  ~burn_status * dist_perim,
  y = cags_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(cagsS_dp)

#### global model ####
cagsS_global <- auto_occ(
  formula = ~onroad + burn_status + rdnbr_het + rdnbr + rdnbrQuad + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim
  ~burn_status + rdnbr_het + rdnbr + rdnbrQuad + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim,
  y = cags_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(cagsS_global)

#### 2 factor: severity ####
cagsS_sev_dp <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim
  ~burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim,
  y = cags_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(cagsS_sev_dp)

#### 2 factor: heterogeneity ####
cagsS_het_dp <- auto_occ(
  formula = ~onroad + burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim
  ~burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim,
  y = cags_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(cagsS_het_dp)

#### SPATIAL MODEL SELECTION ####
compare_models(list(cagsS_null, cagsS_het, cagsS_rdnbr, cagsS_dp, cagsS_sev_dp, cagsS_het_dp),
               digits = 2) # best fit model: rdnbr

#### SPATIAL MODEL PREDICTIONS ####

#### fire heterogeneity ####
cags_rdnbrHet_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  
  rdnbr_het = rdnbrHet,
  onroad = 0
)

# scaling
cags_rdnbrHet_scaled <- cags_rdnbrHet_real
cags_rdnbrHet_scaled$rdnbr_het <- (
  cags_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_cags <- predict(
  object = cagsS_het,
  type = "rho",
  newdata = cags_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_cags <- data.frame(
  rdnbrHet_cags,
  cags_rdnbrHet_real
)

# plot 
ggplot(rdnbrHet_cags, aes(rdnbr_het, estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Fire Heterogeneity", y = "Detection", title = "California Ground Squirrel") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*fire heterogeneity ####
cags_rdnbrHet_real <- data.frame(
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
cags_rdnbrHet_scaled <- cags_rdnbrHet_real
cags_rdnbrHet_scaled$rdnbr_het <- (
  cags_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_cags <- predict(
  object = cagsS_het,
  type = "rho",
  newdata = cags_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_cags <- data.frame(
  rdnbrHet_cags,
  cags_rdnbrHet_real
)

rdnbrHet_cags <- split(rdnbrHet_cags, rdnbrHet_cags$burn_status)

#plot
plot(1~1, type = "n", xlab = "Fire Heterogeneity", ylab = "Detection",
     xlim = c(20, 1750), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(rdnbrHet_cags$burned$rdnbr_het, rev(rdnbrHet_cags$burned$rdnbr_het)),
  y = c(rdnbrHet_cags$burned$lower, rev(rdnbrHet_cags$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = rdnbrHet_cags$burned$rdnbr_het, y = rdnbrHet_cags$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(rdnbrHet_cags$unburned$rdnbr_het, rev(rdnbrHet_cags$unburned$rdnbr_het)),
  y = c(rdnbrHet_cags$unburned$lower, rev(rdnbrHet_cags$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = rdnbrHet_cags$unburned$rdnbr_het, y = rdnbrHet_cags$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topright", c("Unburned", "Burned"), lwd = 2, lty = c(2, 1), col = c("darkgreen", "darkorange"), bty = "n", cex=1.25)

ca1 <- recordPlot() # occupancy
ca2 <- recordPlot() # detection

#### fire severity ####
cags_rdnbr_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("unburned", "burned")
  ),
  rdnbr = rdnbr_seq,
  rdnbrQuad = 0,
  dist_perim = 0,
  onroad = 0
)

# scaling
cags_rdnbr_scaled <- cags_rdnbr_real
cags_rdnbr_scaled$rdnbr <- (
  cags_rdnbr_scaled$rdnbr - mean(post_varCov$rdnbr_abs)
) / sd(post_varCov$rdnbr_abs)

# the model prediction
rdnbr_cags <- predict(
  object = cagsS_sev_dp,
  type = "psi",
  newdata = cags_rdnbr_scaled
)

# add on covariate data
rdnbr_cags <- data.frame(
  rdnbr_cags,
  cags_rdnbr_real
)

ca3 <- ggplot(rdnbr_cags, aes(x = rdnbr, y = estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Fire Severity (RdNBR)", y = "Occupancy") + 
  ylim(0,1) +
  theme_classic(18)

ca3 # occupancy
ca4 # detection

#### distance to fire perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

cags_dp_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  dist_perim = dp_seq,
  rdnbr = 0, 
  rdnbrQuad = 0,
  onroad = 0
)

# scaling
cags_dp_scaled <- cags_dp_real
cags_dp_scaled$dist_perim <- (
  cags_dp_scaled$dist_perim - mean(post_varCov$dist_perim)
) / sd(post_varCov$dist_perim)

# the model prediction
dp_cags <- predict(
  object = cagsS_sev_dp,
  type = "rho",
  newdata = cags_dp_real
)

# add on covariate data
dp_cags <- data.frame(
  dp_cags,
  cags_dp_real
)

ggplot(dp_cags, aes(x = dist_perim, y = estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Distance to Fire Perimeter (km)", y = "Occupancy") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*distance to perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

cags_dp_real <- data.frame(
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
dp_cags <- predict(
  object = cagsS_sev_dp,
  type = "rho",
  newdata = cags_dp_real
)

dp_cags$dist_perim <- c(-6:0, 0:9)
dp_cags$burn_status <- cags_dp_real$burn_status

dp_cags <- split(dp_cags, cags_dp_real$burn_status)

#plot
plot(1~1, type = "n", xlab = "Distance to Fire Perimeter (km)", ylab = "Detection",
     xlim = c(-5, 5), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(dp_cags$burned$dist_perim, rev(dp_cags$burned$dist_perim)),
  y = c(dp_cags$burned$lower, rev(dp_cags$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = dp_cags$burned$dist_perim, y = dp_cags$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(dp_cags$unburned$dist_perim, rev(dp_cags$unburned$dist_perim)),
  y = c(dp_cags$unburned$lower, rev(dp_cags$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = dp_cags$unburned$dist_perim, y = dp_cags$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topright", c("Burned", "Unburned", "Fire Perimeter"), lwd = 2, lty = c(1, 2, 1), col = c("darkorange", "darkgreen", "black"), bty = "n", cex=1.25)

ca5 <- recordPlot()
ca6 <- recordPlot()

#### TEMPORAL ANALYSIS ####

#### null model ####
cagsT_null <- auto_occ(
  formula = ~1~1,
  y = cags_T
)

# overall expected occupancy: 0.070
(
  intercept_preds_psi <- predict(
    cagsT_null,
    type = "psi"
  )
)

# overall expected detection: 0.343
(
  intercept_preds_rho <- predict(
    cagsT_null, 
    type = "rho"
  )
)

#### global model ####
cagsT_global <- auto_occ(
  formula = ~onroad + burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi
  ~burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi,
  y = cags_T,
  det_covs = covFrame_T,
  occ_covs = covFrame_T
)
summary(cagsT_global)

#### TEMPORAL MODEL SELECTION ####
compare_models(list(cagsT_null, cagsT_global),
               digits = 2) # global better than null

#### TEMPORAL MODEL PREDICTIONS ####

#### burn status ####
burned_seq <- c("prefire", "burned_post", "unburned_post")
season_levels <- levels(cov_frame$season$V1)

cags_burn_real <- data.frame(
  burn_status = factor(burned_seq, 
                       levels = burned_seq),
  tsf = -3,
  ndvi = 0, 
  ndvi_het = 0,
  onroad = 0
)

# the model prediction
burn_cags <- predict(
  object = cagsT_global,
  type = "rho",
  newdata = cags_burn_real
)

# add on the covariate data
burn_cags <- data.frame(
  burn_cags,
  cags_burn_real
)

# plot
ggplot(burn_cags, aes(x = factor(burn_status,
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

ca7 # occupancy
ca8 # detection

#### burn_status*tsf ####
cags_int <- data.frame(
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
cags_int$burn_status <- factor(cags_int$burn_status,
                              levels = c("prefire", "burned_post", "unburned_post"))

cags_int_pred <- predict(cagsT_global,
                        type = 'psi',
                        newdata = cags_int,
                        interval = "confidence")

cags_int_pred$tsf <- c(-22:0, 0:22, 0:22)
cags_int_pred$burn_status <- cags_int$burn_status

cags_int_pred <- split(cags_int_pred, cags_int_pred$burn_status)

#plot
plot(1~1, type = "n", xlab = "Seasons Since Fire", ylab = "Occupancy",
     xlim = c(-22, 23), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(cags_int_pred$prefire$tsf, rev(cags_int_pred$prefire$tsf)),
  y = c(cags_int_pred$prefire$lower, rev(cags_int_pred$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = cags_int_pred$prefire$tsf, y = cags_int_pred$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(cags_int_pred$burned_post$tsf, rev(cags_int_pred$burned_post$tsf)),
  y = c(cags_int_pred$burned_post$lower, rev(cags_int_pred$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = cags_int_pred$burned_post$tsf, y = cags_int_pred$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(cags_int_pred$unburned_post$tsf, rev(cags_int_pred$unburned_post$tsf)),
  y = c(cags_int_pred$unburned_post$lower, rev(cags_int_pred$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = cags_int_pred$unburned_post$tsf, y = cags_int_pred$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topright", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

ca9 <- recordPlot() # occupancy
ca10 <- recordPlot() # detection

#### ndvi ####
ndvi_seq <- seq(0, 0.30, 0.01)

cags_ndvi_real <- data.frame(
  burn_status = factor("burned_post",
                       levels = c("prefire", "burned_post", "unburned_post")
  ),
  tsf = -3,
  ndvi = ndvi_seq,
  onroad = 0
)

# scaling
cags_ndvi_scaled <- cags_ndvi_real
cags_ndvi_scaled$ndvi <- (
  cags_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_cags <- predict(
  object = cagsT_global,
  type = "psi",
  newdata = cags_ndvi_scaled
)

# add on covariate data
ndvi_cags <- data.frame(
  ndvi_cags,
  cags_ndvi_real
)

# plot
ggplot(ndvi_cags, aes(x = ndvi, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Vegetation Greenness", y = "Detection", title = "California Ground Squirrel") + 
  ylim(0,1) +
  theme_classic(18)

#### burn_status*ndvi ####
ndvi_seq <- seq(0, 0.35, 0.01)

cags_ndvi_real <- data.frame(
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
cags_ndvi_scaled <- cags_ndvi_real
cags_ndvi_scaled$ndvi <- (
  cags_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_cags <- predict(
  object = cagsT_global,
  type = "rho",
  newdata = cags_ndvi_scaled
)

# add on covariate data
ndvi_cags <- data.frame(
  ndvi_cags,
  cags_ndvi_real
)

ndvi_cags <- split(ndvi_cags, ndvi_cags$burn_status)

#plot
plot(1~1, type = "n", xlab = "NDVI", ylab = "Detection",
     xlim = c(0, 0.3), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(ndvi_cags$prefire$ndvi, rev(ndvi_cags$prefire$ndvi)),
  y = c(ndvi_cags$prefire$lower, rev(ndvi_cags$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = ndvi_cags$prefire$ndvi, y = ndvi_cags$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(ndvi_cags$burned_post$ndvi, rev(ndvi_cags$burned_post$ndvi)),
  y = c(ndvi_cags$burned_post$lower, rev(ndvi_cags$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = ndvi_cags$burned_post$ndvi, y = ndvi_cags$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(ndvi_cags$unburned_post$ndvi, rev(ndvi_cags$unburned_post$ndvi)),
  y = c(ndvi_cags$unburned_post$lower, rev(ndvi_cags$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = ndvi_cags$unburned_post$ndvi, y = ndvi_cags$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topleft", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

ca11 <- recordPlot() # occupancy
ca12 <- recordPlot() # detection

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

estimate_h <- round(cagsS_het@estimates$Est, 3)
estimate_r <- round(cagsS_rdnbr@estimates$Est, 3)
estimate_d <- round(cagsS_dp@estimates$Est, 3)

SE_h <- round(cagsS_het@estimates$SE, 3)
SE_r <- round(cagsS_rdnbr@estimates$SE, 3)
SE_d <- round(cagsS_dp@estimates$SE, 3)

p_h <- round(cagsS_het@estimates$p, 3)
p_r <- round(cagsS_rdnbr@estimates$p, 3)
p_d <- round(cagsS_dp@estimates$p, 3)

# fire heterogeneity 
cagsHet_tab <- bind_cols(param_h, estimate_h, SE_h, p_h)
cagsHet_tab <- cagsHet_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
cagsHet_tab <- cagsHet_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

cags_h <- cagsHet_tab %>% 
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
  tab_spanner(label = md("Fire Heterogeneity Model (AIC = 382.6)"),
              columns = everything())

cags_h %>% 
  gtsave(filename = "cags_h.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# fire severity
cagsR_tab <- bind_cols(param_r, estimate_r, SE_r, p_r)
cagsR_tab <- cagsR_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
cagsR_tab <- cagsR_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

cags_r <- cagsR_tab %>% 
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
  tab_spanner(label = md("Fire Severity Model (AIC = 379.4)"),
              columns = everything())

cags_r %>%  
  gtsave(filename = "cags_r.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# distance to fire perimeter
cagsD_tab <- bind_cols(param_d, estimate_d, SE_d, p_d)
cagsD_tab <- cagsD_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
cagsD_tab <- cagsD_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

cags_d <- cagsD_tab %>% 
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
  tab_spanner(label = md("Distance to Fire Perimeter Model (AIC = 382.4)"),
              columns = everything())

cags_d %>% 
  gtsave(filename = "cags_d.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

#### temporal ####
parameter <- c(
  paste(greeks("psi"), "- Intercept"), paste(greeks("psi"), "- Burn Status (Burned)"), 
  paste(greeks("psi"), "- Burn Status (Unburned)"), paste(greeks("psi"), "- Time Since Fire"),
  paste(greeks("psi"), "- NDVI"), paste(greeks("psi"), "- Burned x Time Since Fire"), 
  paste(greeks("psi"), "- Unburned x Time Since Fire"), paste(greeks("psi"), "- Burned x NDVI"),
  paste(greeks("psi"), "- Unburned x NDVI"), paste(greeks("psi"), "-", greeks("theta")), 
  paste(greeks("rho"), "- Intercept"), paste(greeks("rho"), "- Microsite Attractant"), paste(greeks("rho"), "- Burn Status (Burned)"), 
  paste(greeks("rho"), "- Burn Status (Unburned)"), paste(greeks("rho"), "- Time Since Fire"),
  paste(greeks("rho"), "- NDVI"), paste(greeks("rho"), "- Burned x Time Since Fire"), 
  paste(greeks("rho"), "- Unburned x Time Since Fire"), paste(greeks("rho"), "- Burned x NDVI"),
  paste(greeks("rho"), "- Unburned x NDVI")
)

estimate <- cagsT_global@estimates$Est
estimate <- round(estimate, 3)

SE <- cagsT_global@estimates$SE
SE <- round(SE, 3)

p <- cagsT_global@estimates$p
p <- round(p, 3)

cagsGlobal_tab <- bind_cols(parameter, estimate, SE, p)
cagsGlobal_tab <- cagsGlobal_tab %>% 
  rename("Parameter" = ...1, 
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
cagsGlobal_tab <- cagsGlobal_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

cags_g <- cagsGlobal_tab %>% 
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

cags_g %>% 
  gtsave(filename = "cags_g.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")



