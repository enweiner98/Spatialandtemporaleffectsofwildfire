#### RABBIT FINAL ANALYSIS #####
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

#### RABBIT DATAFRAMES ####
rab_all <- all_occu$rabbit_sp  # all rabbit detections
rab_post <- post_occu$rabbit_sp  # post-fire rabbit detections

#### summary statistics ####
# temporal: total detections
rab_all %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 361 total detections

# detections at burned v. unburned sites 
rab_bs <- rab_all
rab_bs$burn_status <- var_cov$burned

bs_count <- rab_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 12 detections 
# unburned: 186 detections 
# prefire: 163 detections

# temporal: total sites 
tot_sites <- rab_all %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 28 sites 

# spatial: total detections
rab_post %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 154 total detections

# detections at burned v. unburned sites 
rab_bs <- rab_post
rab_bs$burn_status <- post_varCov$burned

bs_count <- rab_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 109 detections 
# unburned: 45 detections 

# spatial: total sites 
tot_sites <- rab_post %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 17 sites 

## creating dataframes for autoOcc
# post-fire only
rabbit_S <- format_y(
  x = rab_post,
  site_column = "Site",
  time_column = "Season",
  history_column = "Week"
)

# includes pre-fire and post-fire data
rabbit_T <- format_y(        
  x = rab_all,
  site_column = "Site",
  time_column = "Season",
  history_columns = "Week" 
)

#### SPATIAL ANALYSIS ####
#### null model ####
rabS_null <- auto_occ(
  formula = ~1~1,
  y = rabbit_S
)

# overall expected occupancy: 0.162
(
  intercept_preds_psi <- predict(
    rabS_null,
    type = "psi"
  )
)

# overall expected detection: 0.516
(
  intercept_preds_rho <- predict(
    rabS_null, 
    type = "rho"
  )
)

#### p1: fire heterogeneity ####
rabS_het <- auto_occ(
  formula = ~onroad + burn_status * rdnbr_het
  ~burn_status * rdnbr_het,
  y = rabbit_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(rabS_het)

#### p2: fire severity ####
rabS_rdnbr <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad
  ~burn_status + rdnbr + rdnbrQuad,
  y = rabbit_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(rabS_rdnbr)

#### p3: distance to fire perimeter ####
rabS_dp <- auto_occ(
  formula = ~onroad + burn_status * dist_perim
  ~burn_status * dist_perim,
  y = rabbit_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(rabS_dp)

#### 2 factor: severity ####
rabS_sev_dp <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim
  ~burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim,
  y = rabbit_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(rabS_sev_dp)

#### 2 factor: heterogeneity ####
rabS_het_dp <- auto_occ(
  formula = ~onroad + burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim
  ~burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim,
  y = rabbit_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(rabS_het_dp)

#### SPATIAL MODEL SELECTION ####
compare_models(list(rabS_null, rabS_het, rabS_rdnbr, rabS_dp, rabS_sev_dp, rabS_het_dp),
               digits = 2) # best fit model: rdnbr

#### SPATIAL MODEL PREDICTIONS ####

#### fire heterogeneity ####
rab_rdnbrHet_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  
  rdnbr_het = rdnbrHet,
  dist_perim = 0,
  onroad = 0
)

# scaling
rab_rdnbrHet_scaled <- rab_rdnbrHet_real
rab_rdnbrHet_scaled$rdnbr_het <- (
  rab_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_rab <- predict(
  object = rabS_het_dp,
  type = "rho",
  newdata = rab_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_rab <- data.frame(
  rdnbrHet_rab,
  rab_rdnbrHet_real
)

# plot 
ggplot(rdnbrHet_rab, aes(rdnbr_het, estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Fire Heterogeneity", y = "Detection") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*fire heterogeneity ####
rab_rdnbrHet_real <- data.frame(
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
rab_rdnbrHet_scaled <- rab_rdnbrHet_real
rab_rdnbrHet_scaled$rdnbr_het <- (
  rab_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_rab <- predict(
  object = rabS_het,
  type = "rho",
  newdata = rab_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_rab <- data.frame(
  rdnbrHet_rab,
  rab_rdnbrHet_real
)

rdnbrHet_rab <- split(rdnbrHet_rab, rdnbrHet_rab$burn_status)

#plot
plot(1~1, type = "n", xlab = "Fire Heterogeneity", ylab = "Detection",
     xlim = c(20, 1750), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(rdnbrHet_rab$burned$rdnbr_het, rev(rdnbrHet_rab$burned$rdnbr_het)),
  y = c(rdnbrHet_rab$burned$lower, rev(rdnbrHet_rab$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = rdnbrHet_rab$burned$rdnbr_het, y = rdnbrHet_rab$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(rdnbrHet_rab$unburned$rdnbr_het, rev(rdnbrHet_rab$unburned$rdnbr_het)),
  y = c(rdnbrHet_rab$unburned$lower, rev(rdnbrHet_rab$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = rdnbrHet_rab$unburned$rdnbr_het, y = rdnbrHet_rab$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topleft", c("Unburned", "Burned"), lwd = 2, lty = c(2, 1), col = c("darkgreen", "darkorange"), bty = "n", cex = 1.25)

rab1 <- recordPlot() # occupancy
rab2 <- recordPlot() # detection

#### fire severity ####
rab_rdnbr_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("unburned", "burned")
  ),
  rdnbr = rdnbr_seq,
  rdnbrQuad = 0,
  dist_perim = 0,
  onroad = 0
)

# scaling
rab_rdnbr_scaled <- rab_rdnbr_real
rab_rdnbr_scaled$rdnbr <- (
  rab_rdnbr_scaled$rdnbr - mean(post_varCov$rdnbr_abs)
) / sd(post_varCov$rdnbr_abs)

# the model prediction
rdnbr_rab <- predict(
  object = rabS_rdnbr,
  type = "rho",
  newdata = rab_rdnbr_scaled
)

# add on covariate data
rdnbr_rab <- data.frame(
  rdnbr_rab,
  rab_rdnbr_real
)

rab4 <- ggplot(rdnbr_rab, aes(x = rdnbr, y = estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Fire Severity (RdNBR)", y = "Detection") + 
  ylim(0,1) +
  theme_classic(18)

rab3 # occupancy
rab4 # detection

#### distance to fire perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

rab_dp_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  dist_perim = dp_seq,
  onroad = 0
)

# scaling
rab_dp_scaled <- rab_dp_real
rab_dp_scaled$dist_perim <- (
  rab_dp_scaled$dist_perim - mean(post_varCov$dist_perim)
) / sd(post_varCov$dist_perim)

# the model prediction
dp_rab <- predict(
  object = rabS_dp,
  type = "rho",
  newdata = rab_dp_real
)

# add on covariate data
dp_rab <- data.frame(
  dp_rab,
  rab_dp_real
)

ggplot(dp_rab, aes(x = dist_perim, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Distance to Fire Perimeter (km)", y = "Detection", title = "Rabbit") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*distance to perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

rab_dp_real <- data.frame(
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
dp_rab <- predict(
  object = rabS_het_dp,
  type = "psi",
  newdata = rab_dp_real
)

dp_rab$dist_perim <- c(-6:0, 0:9)
dp_rab$burn_status <- rab_dp_real$burn_status

dp_rab <- split(dp_rab, rab_dp_real$burn_status)

#plot
plot(1~1, type = "n", xlab = "Distance to Fire Perimeter (km)", ylab = "Occupancy",
     xlim = c(-5, 5), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(dp_rab$burned$dist_perim, rev(dp_rab$burned$dist_perim)),
  y = c(dp_rab$burned$lower, rev(dp_rab$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = dp_rab$burned$dist_perim, y = dp_rab$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(dp_rab$unburned$dist_perim, rev(dp_rab$unburned$dist_perim)),
  y = c(dp_rab$unburned$lower, rev(dp_rab$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = dp_rab$unburned$dist_perim, y = dp_rab$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topright", c("Burned", "Unburned", "Fire Perimeter"), lwd = 2, lty = c(1, 2, 1), col = c("darkorange", "darkgreen", "black"), bty = "n", cex=1.25)

rab5 <- recordPlot() # occupancy
rab6 <- recordPlot() # detection

#### TEMPORAL ANALYSIS ####

#### null model ####
rabT_null <- auto_occ(
  formula = ~1~1,
  y = rabbit_T
)

# overall expected occupancy: 0.121
(
  intercept_preds_psi <- predict(
    rabT_null,
    type = "psi"
  )
)

# overall expected detection: 0.507
(
  intercept_preds_rho <- predict(
    rabT_null, 
    type = "rho"
  )
)

#### global model ####
rabT_global <- auto_occ(
  formula = ~onroad + burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi
  ~burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi,
  y = rabbit_T,
  det_covs = covFrame_T,
  occ_covs = covFrame_T
)
summary(rabT_global)

#### TEMPORAL MODEL SELECTION ####
compare_models(list(rabT_null, rabT_global), digits = 2)

#### TEMPORAL MODEL PREDICTIONS ####

#### burn status ####
burned_seq <- c("prefire", "burned_post", "unburned_post")
season_levels <- levels(cov_frame$season$V1)

rab_burn_real <- data.frame(
  burn_status = factor(burned_seq, 
                       levels = burned_seq),
  tsf = -3,
  ndvi = 0, 
  ndvi_het = 0,
  onroad = 0
)

# the model prediction
burn_rab <- predict(
  object = rabT_global,
  type = "rho",
  newdata = rab_burn_real
)

# add on the covariate data
burn_rab <- data.frame(
  burn_rab,
  rab_burn_real
)

# plot
rab8 <- ggplot(burn_rab, aes(x = factor(burn_status,
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

rab7 # occupancy
rab8 # detection

#### burn_status*tsf ####
rab_int <- data.frame(
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
rab_int$burn_status <- factor(rab_int$burn_status,
                               levels = c("prefire", "burned_post", "unburned_post"))

rab_int_pred <- predict(rabT_global,
                         type = 'rho',
                         newdata = rab_int,
                         interval = "confidence")

rab_int_pred$tsf <- c(-22:0, 0:22, 0:22)
rab_int_pred$burn_status <- rab_int$burn_status

rab_int_pred <- split(rab_int_pred, rab_int_pred$burn_status)

#plot
plot(1~1, type = "n", xlab = "Seasons Since Fire", ylab = "Detection",
     xlim = c(-22, 23), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(rab_int_pred$prefire$tsf, rev(rab_int_pred$prefire$tsf)),
  y = c(rab_int_pred$prefire$lower, rev(rab_int_pred$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = rab_int_pred$prefire$tsf, y = rab_int_pred$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(rab_int_pred$burned_post$tsf, rev(rab_int_pred$burned_post$tsf)),
  y = c(rab_int_pred$burned_post$lower, rev(rab_int_pred$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = rab_int_pred$burned_post$tsf, y = rab_int_pred$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(rab_int_pred$unburned_post$tsf, rev(rab_int_pred$unburned_post$tsf)),
  y = c(rab_int_pred$unburned_post$lower, rev(rab_int_pred$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = rab_int_pred$unburned_post$tsf, y = rab_int_pred$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topleft", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

rab9 <- recordPlot() # occupancy
rab10 <- recordPlot() # detection

#### ndvi ####
ndvi_seq <- seq(0, 0.30, 0.01)

rab_ndvi_real <- data.frame(
  burn_status = factor("burned_post",
                       levels = c("prefire", "burned_post", "unburned_post")
  ),
  tsf = -3,
  ndvi = ndvi_seq,
  ndvi_het = 0,
  onroad = 0
)

# scaling
rab_ndvi_scaled <- rab_ndvi_real
rab_ndvi_scaled$ndvi <- (
  rab_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_rab <- predict(
  object = rabT_global,
  type = "psi",
  newdata = rab_ndvi_scaled
)

# add on covariate data
ndvi_rab <- data.frame(
  ndvi_rab,
  rab_ndvi_real
)

# plot
ggplot(ndvi_rab, aes(x = ndvi, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Vegetation Greenness", y = "Occupancy", title = "Rabbit") + 
  ylim(0,1) +
  theme_classic(18)

#### burn_status*ndvi ####
ndvi_seq <- seq(0, 0.35, 0.01)

rab_ndvi_real <- data.frame(
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
rab_ndvi_scaled <- rab_ndvi_real
rab_ndvi_scaled$ndvi <- (
  rab_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_rab <- predict(
  object = rabT_global,
  type = "rho",
  newdata = rab_ndvi_scaled
)

# add on covariate data
ndvi_rab <- data.frame(
  ndvi_rab,
  rab_ndvi_real
)

ndvi_rab <- split(ndvi_rab, ndvi_rab$burn_status)

#plot
plot(1~1, type = "n", xlab = "Vegetation Biomass", ylab = "Detection",
     xlim = c(0, 0.3), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(ndvi_rab$prefire$ndvi, rev(ndvi_rab$prefire$ndvi)),
  y = c(ndvi_rab$prefire$lower, rev(ndvi_rab$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = ndvi_rab$prefire$ndvi, y = ndvi_rab$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(ndvi_rab$burned_post$ndvi, rev(ndvi_rab$burned_post$ndvi)),
  y = c(ndvi_rab$burned_post$lower, rev(ndvi_rab$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = ndvi_rab$burned_post$ndvi, y = ndvi_rab$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(ndvi_rab$unburned_post$ndvi, rev(ndvi_rab$unburned_post$ndvi)),
  y = c(ndvi_rab$unburned_post$lower, rev(ndvi_rab$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = ndvi_rab$unburned_post$ndvi, y = ndvi_rab$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topright", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

rab11 <- recordPlot() # occupancy
rab12 <- recordPlot() # detection

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

estimate_h <- round(rabS_het@estimates$Est, 3)
estimate_r <- round(rabS_rdnbr@estimates$Est, 3)
estimate_d <- round(rabS_dp@estimates$Est, 3)

SE_h <- round(rabS_het@estimates$SE, 3)
SE_r <- round(rabS_rdnbr@estimates$SE, 3)
SE_d <- round(rabS_dp@estimates$SE, 3)

p_h <- round(rabS_het@estimates$p, 3)
p_r <- round(rabS_rdnbr@estimates$p, 3)
p_d <- round(rabS_dp@estimates$p, 3)

# fire heterogeneity 
rabHet_tab <- bind_cols(param_h, estimate_h, SE_h, p_h)
rabHet_tab <- rabHet_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
rabHet_tab <- rabHet_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

rab_h <- rabHet_tab %>% 
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
  tab_spanner(label = md("Fire Heterogeneity Model (AIC = 601.8)"),
              columns = everything())

rab_h %>% 
  gtsave(filename = "rab_h.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# fire severity
rabR_tab <- bind_cols(param_r, estimate_r, SE_r, p_r)
rabR_tab <- rabR_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
rabR_tab <- rabR_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

rab_r <- rabR_tab %>% 
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
  tab_spanner(label = md("Fire Severity Model (AIC = 599.6)"),
              columns = everything())

rab_r %>%  
  gtsave(filename = "rab_r.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# distance to fire perimeter
rabD_tab <- bind_cols(param_d, estimate_d, SE_d, p_d)
rabD_tab <- rabD_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
rabD_tab <- rabD_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

rab_d <- rabD_tab %>% 
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
  tab_spanner(label = md("Distance to Fire Perimeter Model (AIC = 607.2)"),
              columns = everything())

rab_d %>% 
  gtsave(filename = "rab_d.html",
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

estimate <- rabT_global@estimates$Est
estimate <- round(estimate, 3)

SE <- rabT_global@estimates$SE
SE <- round(SE, 3)

p <- rabT_global@estimates$p
p <- round(p, 3)

rabGlobal_tab <- bind_cols(parameter, estimate, SE, p)
rabGlobal_tab <- rabGlobal_tab %>% 
  rename("Parameter" = ...1, 
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
rabGlobal_tab <- rabGlobal_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

rab_g <- rabGlobal_tab %>% 
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

rab_g %>% 
  gtsave(filename = "rab_g.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")





